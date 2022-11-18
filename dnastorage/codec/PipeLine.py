'''
PipeLine Class file which builds the encoding/decoding pipeline based on a set of codecs. 
'''
import sys
from dnastorage.codec.codecfile import *
from dnastorage.codec.base import *
from dnastorage.exceptions import *
from dnastorage.codec.base_conversion import *
from dnastorage.strand_representation import *
from dnastorage.primer.primer_util import *
from dnastorage.util.mpi_utils import *
from dnastorage.codec_types import *
from io import *
import random
import math
logger = logging.getLogger('dna.codec.pipeline')
logger.addHandler(logging.NullHandler())

def get_previous_probes(l,index):
    out_list = []
    if index==0: return out_list
    for x in l[index-1::-1]:
        if isinstance(x,Probe):
            out_list.append(x)
        else:
            break
    return out_list[::-1]

#builds cascades of components 
def cascade_build(component_list):
    previous=None
    for c in component_list[::-1]:
        c.set_object(previous)
        c.reverse=False #this flag is a quick work around to make sure the BaseCodec class works ok with the pipeline-based processing direction
        previous=c
    return component_list[0]

class PipeLine(EncodePacketizedFile,DecodePacketizedFile):
    def __init__(self,components,packetsize_bytes,
                 basestrand_bytes, DNA_upper_bound, final_decode_iterations,dna_consolidator=None,cw_consolidator=None,packetizedfile=None,
                 barcode=tuple(),constant_index_bytes=None):
        
        EncodePacketizedFile.__init__(self,None)
        DecodePacketizedFile.__init__(self,None)
 
        self._packetsize_bytes=packetsize_bytes
        self._packetizedFile=packetizedfile
        if packetizedfile !=None:
            self._packetizedFile.packetSize=self._packetsize_bytes
        self._outer_codecs=[]
        self._inner_codecs=[]
        self._cw_to_DNA=[] #These processes are those that turn CW information into DNA
        self._DNA_to_DNA=[] #These processes are those that post-process dna strand information, e.g. primers, cut sites, etc.
        self._basestrand_bytes=basestrand_bytes #base number of bytes in a strand before encoding occurs
        self._DNA_upper_bound=DNA_upper_bound #upper bound DNA length
        self._decode_strands=[] #strands that will be decoded
        self._final_decode_iterations=final_decode_iterations #how many complete final decoding processes should be done
        self._index_bytes = constant_index_bytes #allow constant index_bytes to be used, useful for controlling for DNA strand size
        
        #support the ability for 2 consolidators, dna consolidators may actually need a cw_consolidator to filter repeat indexes
        self._dna_consolidator = dna_consolidator 
        self._cw_consolidator = cw_consolidator 

        assert dna_consolidator!=None or cw_consolidator!=None #should have at least 1
        
        self._final_decode_run=False #flag that tracks the final decode run_hedges
        self._filtered_strands=[] #strands that have been filtered out due to physical pipeline or incoherence during the decode process
        self.mpi=None #mpi not used on default
        self._barcode=barcode
        for index,component in enumerate(components):
            if index>0: prev_component=components[index-1]
            if isinstance(component,BaseOuterCodec):
                if not index==0 and not isinstance(prev_component,BaseOuterCodec):
                    raise PipeLineConstructionError("Error with Outer Codec connecting")
                self._outer_codecs.append(component)            
            elif isinstance(component,CWtoCW):
                if not index==0 and (not isinstance(prev_component,BaseOuterCodec) and not isinstance(prev_component,CWtoCW) and not isinstance(prev_component,Probe)):
                    raise PipeLineConstructionError("Error with inner codec constructions")
                self._inner_codecs+=get_previous_probes(components,index) #add probes
                self._inner_codecs.append(component)
            elif isinstance(component,CWtoDNA):
                if not index==0 and (not isinstance(prev_component,CWtoCW) and  not isinstance(prev_component,CWtoDNA) and not isinstance(prev_component,Probe) and not isinstance(prev_component,BaseOuterCodec)):
                    raise PipeLineConstructionError("Error with DNA to CW constructions")
                self._cw_to_DNA+=get_previous_probes(components,index) #add probes
                self._cw_to_DNA.append(component)

            elif isinstance(component,DNAtoDNA):
                if not index==0 and (not isinstance(prev_component,CWtoDNA) and not isinstance(prev_component,DNAtoDNA) and not isinstance(prev_component,Probe)):
                    raise PipeLineConstructionError("Error with DNA to DNA construction")
                self._DNA_to_DNA+=get_previous_probes(components,index) #add probes 
                self._DNA_to_DNA.append(component)
            elif isinstance(component,Probe):
                #probe should go into next cascade to avoid placement in locations like the outer codec
                if component is components[-1]:
                    #PipeLineConstructionError("Probe placed at end of pipeline")
                    self._DNA_to_DNA.append(component) #allow probe to go on end of pipeline
                continue
            else:
                raise PipeLineConstructionError("Error with type in pipeline construction")
 
        self._outer_cascade = cascade_build(self._outer_codecs)
        #check the outer codecs, make sure that there is a proper conclusion for indexing
        if self._outer_cascade.remainder(int(math.ceil(self._packetsize_bytes/self._basestrand_bytes))) > 1:
            final_divisor = self._outer_cascade.remainder(int(math.ceil(self._packetsize_bytes/self._basestrand_bytes)))
            #print("final divisor {}".format(final_divisor))
            strand_indexer=BaseOuterCodec(final_divisor)
            self._outer_codecs[-1].set_object(strand_indexer) #this should make sure every strand is given an individual index if

        if len(self._inner_codecs)==0: self._inner_codecs.append(BaseCodec())
        self._inner_cascade = cascade_build(self._inner_codecs)
        if len(self._cw_to_DNA)==0: raise PipeLineConstructionError("Need transformations to go from codewords to DNA")
        self._cw_to_DNA_cascade = cascade_build(self._cw_to_DNA)
        if len(self._DNA_to_DNA)==0: self._DNA_to_DNA.append(BaseCodec())
        self._dna_to_dna_cascade = cascade_build(self._DNA_to_DNA)

        
    def _encode_pipeline(self,packet):
        #Break the packet down into the basic strands with the packet and data/parity index
        packet_strands=[]
        strand_bytes=[]
        for b in packet[1]:
            strand_bytes.append(b)
            if (len(strand_bytes)%self._basestrand_bytes)==0:
                packet_strands.append(BaseDNA(codewords=strand_bytes,index_ints=(packet[0],)))
                strand_bytes=[]
      
        if len(strand_bytes)>0:
            while len(strand_bytes)<self._basestrand_bytes: strand_bytes.append(0)
            packet_strands.append(BaseDNA(codewords=strand_bytes,index_ints=(packet[0],)))

        #encode the packet using our set of cascade encoders
        total_packet_strands=self._outer_cascade.encode(packet_strands)
       
        #Grab indexing information
        self._index_bit_set= (math.ceil(math.log(math.ceil(self._packetizedFile.numberOfPackets),2)),)+self._outer_cascade.get_index_bits() #need to take into account the packet index
        self._index_bit_set = (8,)*len(self._barcode)+self._index_bit_set
        total_bits= sum(self._index_bit_set) 
        index_bytes = total_bits//8
        if not total_bits%8==0: index_bytes+=1
        if self._index_bytes!=None: assert self._index_bytes>=index_bytes
        else: self._index_bytes=index_bytes
        final_DNA_strands=[]
        index_pad_array=[0]*(self._index_bytes-index_bytes) #pad array just in case someone uses a larger index bytes than inferred
        for packet_strand in total_packet_strands:
            packet_strand.index_ints= self._barcode+packet_strand.index_ints
            packet_strand.index_bytes=self._index_bytes
            packet_strand.codewords = pack_bits_to_bytes(packet_strand.index_ints,self._index_bit_set)+index_pad_array+packet_strand.codewords #turn index integers to bytes
            #pass through the next steps   
            self._inner_cascade.encode(packet_strand)
            self._cw_to_DNA_cascade.encode(packet_strand)
            self._dna_to_dna_cascade.encode(packet_strand)
           
            if len(packet_strand.dna_strand)>self._DNA_upper_bound:
                raise DNAStrandPayloadWrongSize("Pipeline Encoding, DNA strand is length {}".format(len(packet_strand.dna_strand)))
            final_DNA_strands.append(packet_strand)

        return final_DNA_strands

    def _strand_valid_index(self,strand): #test if strand has valid index
        return not (None in strand.codewords[0:self._index_bytes] or len(strand.codewords)<self._index_bytes)
           
    def _inner_pipeline(self,strand):
        self._cw_to_DNA_cascade.decode(strand)
        self._inner_cascade.decode(strand)
      
    def final_decode(self):
        #performs the final decode on the streamed in strands 
        self._final_decode_run=True
        if isinstance(self._dna_consolidator,DNAConsolidate):
            self._dna_consolidator.mpi = self.mpi #hand off mpi to consolidator
            self._decode_strands = self._dna_consolidator.decode(self._decode_strands)
        after_inner=[]
        for strand in self._decode_strands:
            self._inner_pipeline(strand)
            if not self._strand_valid_index(strand):
                self.filter_strand(strand)
                continue
            strand.index_ints = unpack_bytes_to_indexes(strand.codewords[0:self._index_bytes],self._index_bit_set)
            if len(self._barcode)>0 and tuple(strand.index_ints[:len(self._barcode)])!=self._barcode:
                self.filter_strand(strand) #strand does not match barcode
                continue
            strand.index_ints=strand.index_ints[len(self._barcode):]
            after_inner.append(strand)
        self._decode_strands=after_inner
        #IF MPI: should get strands back to master
        if self.mpi:
            self._decode_strands=object_gather(self._decode_strands,self.mpi)

        packets={}
        #let rank 0 cut up packets
        if not self.mpi or self.mpi and self.mpi.Get_rank()==0:
            if isinstance(self._cw_consolidator,CWConsolidate):
                self._decode_strands=self._cw_consolidator.decode(self._decode_strands)
            #Now we need to collect packets together and get the packet indexes before unrolling the outer encoding
            for s in self._decode_strands:
                s.codewords=s.codewords[self._index_bytes:] #strip off indexing before outer decoding
                if self._outer_cascade.is_zero(s.index_ints) or not self._outer_cascade.valid(s.index_ints) or len(s.codewords)<self._basestrand_bytes:
                    continue
                packets[s.index_ints[0]]=packets.get(s.index_ints[0],[])+[s]

        #parallleize packet computation
        if self.mpi:
            logger.info("Scatter packets")
            packets=object_scatter([_ for  _ in packets.items()],self.mpi)
            packets={x:y for x,y in packets}
        
        #at this point we have our main packets, so lets run through outer decoding
        total_out_datas=[]
        for p in packets:
            output_packet=[]
            for i in range(0,self._final_decode_iterations):
                output_packet=self._outer_cascade.decode(packets[p])
            #data should be ordered correctly inherently coming out of decoding,filters Nones that may come from dropouts
            out_data = [ c if c!=None else 0 for x in output_packet for c in x.codewords]
            total_out_datas.append((p,out_data))

        if self.mpi:
            total_out_datas = object_gather(total_out_datas,self.mpi)
            if self.mpi.Get_rank()!=0:
                logger.info("Rank {} leaving pipeline".format(self.mpi.Get_rank()))
                return #mpi support only up to the outer code

        for p,out_data in total_out_datas:
            self.writeToFile(p,out_data) #write out packet
        self.write()

    def decode(self,strand):
        if self._final_decode_run is True:
            self._final_decode_run=False
            self._decode_strands=[]
            self._filtered_strands=[]
        strand.before_decode=strand.dna_strand
        strand.is_reversed=False
        #perform the stream portion of the encoding pipeline, which processes the physical cascade. This should include stuff like removing physical portions
        self._dna_to_dna_cascade.decode(strand)
        if strand.dna_strand == None:
            strand.dna_strand=strand.before_decode
            try:
                is_rna = strand.is_RNA
            except: #if it is not RNA it should throw exception
                strand.is_reversed=True
                #try reverse_complement
                strand.dna_strand=reverse_complement(strand.dna_strand)
                self._dna_to_dna_cascade.decode(strand)
                if strand.dna_strand ==None:
                    self.filter_strand(strand)
                    return
            else:#no sense in doing reverse complements for RNA 
                assert is_rna
                self.filter_strand(strand)
                return

        strand.index_bytes = self._index_bytes
        strand.index_bit_set =  self._index_bit_set
        self._decode_strands.append(strand)
   

    def filter_strand(self,strand):
        strand.dna_strand=strand.before_decode #single point to revert filtered strand for future processing by pipelines
        strand.is_reversed=False
        self._filtered_strands.append(strand)
        
    def get_filtered(self): #get strands that are not confirmed to be a part of this pipeline
        return self._filtered_strands
    
    def set_read_pf(self,read_pf):
        assert(isinstance(read_pf,ReadPacketizedFilestream))
        read_pf.packetSize = self._packetsize_bytes
        self._packetizedFile=read_pf
        
    def set_write_pf(self,write_pf):
        assert(isinstance(write_pf,WritePacketizedFilestream))
        write_pf.packetSize=self._packetsize_bytes
        self._packetizedFile=write_pf

    def encode(self):
        block = self._encode()
        try:
            block = (block[0],[ ord(_) for _ in block[1] ])
        except Exception as e:
            # This is pretty horrible, kids don't try this at home!
            block = (block[0],[ _ for _ in block[1] ])
        return self._encode_pipeline(block) # get entire block


    def encode_header_data(self):
        #encodes indexing information into an array of bytes
        data=[]
        data+=convertIntToBytes(self._index_bytes,1)
        data+=convertIntToBytes(len(self._index_bit_set),1)
        for i in self._index_bit_set:
            data+=convertIntToBytes(i,1)

        #get data that codecs want to be preserved across instantiations
        data+=self._inner_cascade.encode_header()
        data+=self._cw_to_DNA_cascade.encode_header()
        data+=self._dna_to_dna_cascade.encode_header()
        data+=self._outer_cascade.encode_header()

        return data

    def decode_header_data(self,buff):
        pos=0
        self._index_bytes=convertBytesToInt(buff[pos:pos+1])
        pos+=1
        len_index_bit_set = convertBytesToInt(buff[pos:pos+1])
        pos+=1
        index_bit_set=[]
        for i in range(len_index_bit_set):
            index_bit_set+=[convertBytesToInt(buff[pos:pos+1])]
            pos+=1

        self._index_bit_set=tuple(index_bit_set)
        buf = buff[pos::]
        #decode header data for all the cascades
        buf=self._inner_cascade.decode_header(buf)
        buf=self._cw_to_DNA_cascade.decode_header(buf)
        buf=self._dna_to_dna_cascade.decode_header(buf)
        buf=self._outer_cascade.decode_header(buf)
        return buf

    @property
    def mpi(self):
        return self._mpi
    @mpi.setter
    def mpi(self,comm):
        if not ("mpi4py" in sys.modules or "mpi4py.MPI" in sys.modules):
            raise SystemError("mpi4py has not been loaded")
        self._mpi = comm #note, this should be a communicator of processes that will work together during decoding
        
if __name__=="__main__":
    from dnastorage.codec.consolidation import *
    from dnastorage.codec.block import *
    from dnastorage.codec.strand import *
    from dnastorage.codec.phys import *
    from dnastorage.codec.commafreecodec import *
    
    
    #perform simple encode-decode on some data
    mem_buffer=BytesIO()
    pf = ReadPacketizedFilestream(mem_buffer)

    
    #define encoding architecture
    strand_length_bytes=10
    dna_length=500
    decode_iterations=1
    #packet_in_strands = 185
    #packet_size_bytes = strand_length_bytes*packet_in_strands
    packet_size_bytes=1000*4
    rs1_outer = 33
    rs2_outer= 10

    rs1Outer = ReedSolomonOuterPipeline(rs1_outer,10)
    rs2Outer = ReedSolomonOuterPipeline(rs2_outer,10)
    rsInner = ReedSolomonInnerCodecPipeline(2)
    commafree = CommaFreeCodecPipeline(0)
    append = AppendSequencePipeline("A"*20)
    prepend = PrependSequencePipeline("T"*20)
    consolidator = SimpleMajorityVote()

    pipeline = PipeLine((rs1Outer,rs2Outer,rsInner,commafree,append,prepend),consolidator,packet_size_bytes,strand_length_bytes,
                        dna_length,1)

    pipeline.set_read_pf(pf)

    #do some dummy writes

    for i in range(1000):
        buff = bytearray(convertIntToBytes(i,4))
        tell = mem_buffer.tell()
        mem_buffer.seek(0,2)
        # convert string
        mem_buffer.write(buff)
        mem_buffer.seek(tell,0)

    strands=[]
    for b in pipeline:
        for s in b:
            strands.append(s)
    mem_buffer_after_decode=BytesIO()
    #reverse the strands 
    write_pf = WritePacketizedFilestream(mem_buffer_after_decode,1000*4,minKey=0,zeroFillMissing=True,packetSize=0)
    pipeline.set_write_pf(write_pf)

    for s in strands:
        pipeline.decode(s)

    pipeline.final_decode()

    mem_buffer.seek(0)
    mem_buffer_after_decode.seek(0)
    while True:
        s1= mem_buffer.read(4)
        s2= mem_buffer_after_decode.read(4)
        if len(s1)==0 or len(s2)==0:
            break
        else:
            int1 = convertBytesToInt([x for x in s1])
            int2 = convertBytesToInt([x for x in s2])
            #print("Int1 {} Int2 {}".format(int1,int2))
            assert(int1==int2)

    
    #test index encodings
    print(pipeline._index_bit_set)
    print(pipeline._index_bytes)
    buf=pipeline.encode_header_data()
    pipeline._index_bytes=0
    pipeline._index_bit_set=()
    print(pipeline.decode_header_data(buf))
    print(pipeline._index_bit_set)
    print(pipeline._index_bytes)
            
    print("Done")
        
    
    
