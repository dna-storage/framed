from dnastorage.system.formats import *
from dnastorage.system.header_class import *
from dnastorage.strand_representation import *
import io
import sys
import copy

import logging
logger = logging.getLogger('dna.storage.system.dnafile')
logger.addHandler(logging.NullHandler())


DATA_BARCODE=0xCC


class DNAFilePipeline:
    def __init__(self):
        return
    
    @classmethod
    def open(self, filename, op, format_name="", write_incomplete_file=False, header_version="0.1",fsmd_abbrev='FSMD',fsmd_header_filename=None,input_strands=None,encoder_params={},header_params={}):
        # check if we are reading or writing
        if op=="rf":
            return  SegmentedReadDNAFile(write_incomplete_file=write_incomplete_file,\
                                        fsmd_abbrev=fsmd_abbrev,\
                                        strands=strands)
        if op=="r":            
            # 1. filename is the input file of strands
            # 2. primer_info optionally tells us the primers we're looking,
            #    if not specified, just assume first first 20 bases are
            #    are the primer.
            # 3. format is optional for reading and often ignored. Instead we look
            #    at the header to deduce what to do.
            #
            
            # Ugly: but we have to find the header to know what to do!
            if filename is not None and os.path.exists(filename):
                fd = open(filename,"r")
                logger.debug("open {} for reading.".format(filename))
                strands = get_strands(fd)
            elif input_strands is not None:
                strands=input_strands

          
            assert strands is not None
            header = Header(header_version,header_params)
            if fsmd_header_filename!=None:
                with open(fsmd_header_filename,"rb") as serialized_header_pipeline_data:
                    header.set_pipeline_data(serialized_header_pipeline_data.read())
                    #TODO: this should be changed to something like getting an encoder but for the header
            h = header.decode_file_header(copy.deepcopy(strands)) 

            if h is None: return None #couldnt decode header, don't really have much else to go off of

            logger.debug("decoded header: {}".format(h))
            if h['formatid'] == file_system_formatid_by_abbrev("Segmented"):
                logger.debug("SegmentedReadDNAFile({},{},{})".format(filename,primer5,primer3))
                return SegmentedReadDNAFilePipeline(input=filename,\
                                            write_incomplete_file=write_incomplete_file,\
                                            fsmd_abbrev=fsmd_abbrev)
            else:
                return ReadDNAFilePipeline(input=filename,
                                           fsmd_abbrev=fsmd_abbrev,input_strands=input_strands,encoder_params=encoder_params,
                                           header_version=header_version,header_params=header_params,fsmd_header_filename=fsmd_header_filename)
            
        elif "w" in op and "s" in op:
            return SegmentedWriteDNAFilePipeline(output=filename,
                                                 format_name=format_name,fsmd_abbrev=fsmd_abbrev,fsmd_header_filename=fsmd_header_filename)       
        elif op=="w":
            return WriteDNAFilePipeline(output=filename,
                                        format_name=format_name,fsmd_abbrev=fsmd_abbrev,encoder_params=encoder_params,
                                        header_version=header_version,header_params=header_params,fsmd_header_filename=fsmd_header_filename)
        else:
            return None

    def flush(self):
        return

    def close(self):
        return

    def readable(self):
        assert False
    def writable(self):
        assert False

def get_strands(in_fd):
    strands = []
    while True:
        s = in_fd.readline()        
        if len(s) == 0:
            break
        s = s.strip()
        if s.startswith('%'):
            continue
        strands.append(BaseDNA(dna_strand=s))
    return strands

class ReadDNAFilePipeline(DNAFilePipeline):
    # ReadDNAFile reads a set of strands from a file.  It finds the header,
    # determines compatibility and encoding type, and then decodes the file.
    #
    def __init__(self,**kwargs):
        DNAFilePipeline.__init__(self)
        self._enc_opts=kwargs["encoder_params"]
        self._header_params=kwargs["header_params"]
        self._header_version=kwargs["header_version"]
        if 'input' in kwargs and kwargs["input"]!=None:
            self.input_filename = kwargs['input']
            self.in_fd = open(self.input_filename,"r")
            strands = get_strands(self.in_fd)
        elif 'in_fd' in kwargs:
            self.in_fd = kwargs['in_fd']
            self.input_filename = ""
        elif "input_strands" in kwargs:
            strands=kwargs["input_strands"]
        else:
            assert 0
            
        if not ('fsmd_abbrev' in kwargs):
            self.fsmd_abbrev = 'FSMD'
        else:
            self.fsmd_abbrev = kwargs['fsmd_abbrev']

        header = Header(self._header_version,self._header_params)

        if 'fsmd_header_filename' in kwargs and kwargs["fsmd_header_filename"]!=None: #read into the header pipeline data that described its encoding process
            with open(kwargs["fsmd_header_filename"],"rb") as serialized_header_pipeline_data:
                header.set_pipeline_data(serialized_header_pipeline_data.read())
                
                
        h = header.decode_file_header(strands)
        self.strands = header.pick_nonheader_strands()

        
        self.formatid = h['formatid']
        self.header = h 
        self.size = h['size']
        
        # set up mem_buffer 
        self.mem_buffer = BytesIO()
        
        if self.formatid == 0x1000:
            # let sub-classes handle initialization
            return

        constructor_function = file_system_decoder(self.formatid)
            
        self.mem_buffer = BytesIO()
        self.pf = WritePacketizedFilestream(self.mem_buffer,self.size,0)

        self.pipe = constructor_function(self.pf,**self._enc_opts,barcode=(DATA_BARCODE,))
        self.pipe.decode_header_data(self.header["other_data"])
        
        for s in self.strands:
            self.pipe.decode(s)
        self.pipe.final_decode()

        self.mem_buffer.seek(0,0) # set read point at beginning of buffer
        return

    def read(self, n=1):        
        return self.mem_buffer.read(n)
    def readline(self, n=-1):
        return self.mem_buffer.readline(n)

    def readable(self):
        return True
    def writable(self):
        return False
        
class WriteDNAFilePipeline(DNAFilePipeline):
    # WriteDNAFile writes a set of strands.
    def __init__(self,**kwargs):
        DNAFilePipeline.__init__(self)
        self._enc_opts=kwargs["encoder_params"]
        self._header_params=kwargs["header_params"]
        self._header_version=kwargs["header_version"]
        self.out_fd=None
        self._header_fd=None
        self.output_filename="none.dna"
        if "fsmd_header_filename" in kwargs and kwargs["fsmd_header_filename"]!=None:
            self._header_fd=open(kwargs["fsmd_header_filename"],"wb+")
        
        if 'formatid' in kwargs:            
            enc_func = file_system_encoder(kwargs['formatid'])
            self.formatid = kwargs['formatid']
        elif 'format_name' in kwargs:
            enc_func = file_system_encoder_by_abbrev(kwargs['format_name'])
            self.formatid = file_system_formatid_by_abbrev(kwargs['format_name'])

        if not ('fsmd_abbrev' in kwargs):
            self.fsmd_abbrev = 'FSMD'
        else:
            self.fsmd_abbrev = kwargs['fsmd_abbrev']

        self.mem_buffer = BytesIO()
        self.pf= ReadPacketizedFilestream(self.mem_buffer)
        self.pipe = enc_func(self.pf,**self._enc_opts,barcode=(DATA_BARCODE,))
        if 'output' in kwargs and kwargs["output"] !=None :
            self.output_filename = kwargs['output']
            self.out_fd = open(self.output_filename,"w")
        elif 'out_fd' in kwargs:
            self.out_fd = kwargs['out_fd']
            self.output_filename = ""
            
        self.size = 0
        self.strands = []
      
        return

    def _encode_buffer(self):
        for block in self.pipe:
            if type(block) == list:
                for s in block:
                    self.strands.append(s)
            else:
                self.strands.append(block)

    def read(self, size):
        assert False and "Do not read at the same time as writing"
        
    def writable(self):
        return True
    def readable(self):
        return False
    
    def write(self, buff):
        self.size += len(buff)
        tell = self.mem_buffer.tell()
        self.mem_buffer.seek(0,2)

        try:
            # convert string
            self.mem_buffer.write(buff)
            self.mem_buffer.seek(tell,0)
        except Exception as e:
            buff = bytearray([x for x in buff])
            self.mem_buffer.write(buff)
            self.mem_buffer.seek(tell,0)
        return
    
    def flush(self):
        self._encode_buffer()
        return

    def close(self):
        self.flush()
        header = Header(self._header_version,self._header_params)

        
        header_dict={}
        header_dict["filename"]=self.output_filename
        header_dict["formatid"]=self.formatid
        header_dict["size"]=self.size
        hdr_strands= header.encode_file_header(header_dict,self.pipe.encode_header_data())

        for i,h in enumerate(hdr_strands):
            self.strands.insert(i,h)

        #write out meta data for header to digital data for recovery later
        if self._header_fd!=None:
            self._header_fd.write(bytearray(header.get_header_pipeline_data()))
            self._header_fd.close()
            
        if self.out_fd is None: return
        comment = header.encode_file_header_comments(header_dict)
        comment+="% Primer 5 : {}\n".format(self._enc_opts["primer5"])
        comment+="% Primer 3 : {}\n".format(self._enc_opts["primer3"])
        self.out_fd.write(comment)
        for s in self.strands:
            self.out_fd.write("{}\n".format(s.dna_strand))

        if self.out_fd != sys.stdout and self.out_fd != sys.stderr:
            self.out_fd.close()
        return
    
class SegmentedWriteDNAFilePipeline(WriteDNAFilePipeline):
    # SegmentedWriteDNAFile writes a set of strands.
    def __init__(self,**kwargs):     
        WriteDNAFilePipeline.__init__(self,**kwargs)
        self.segments = []
        self.beginIndex = 0
        return

    def _record_segment(self):
        self.segments += [[ self.formatid, self.size, self.primer5, self.primer3, self.beginIndex, self.flanking_primer5, self.flanking_primer3 ]]

    def new_segment(self, format_name, primer5, primer3, flanking_primer5="", flanking_primer3=""):
        self._encode_buffer()  # write everything in the buffer to the file
        self._record_segment() # remember new segment
        self.primer5 = primer5
        self.primer3 = primer3
        self.flanking_primer5 = flanking_primer5
        self.flanking_primer3 = flanking_primer3
        # get current index + 1
        self.beginIndex = self.enc.index
        #print "beginIndex={}".format(self.enc.index)
        enc_func = file_system_encoder_by_abbrev(format_name)
        self.formatid = file_system_formatid_by_abbrev(format_name)
        # we consumed the prior buffer, so just make a new one to avoid
        # cursor positioning problems (JMT: not sure if this is the best way)
        self.mem_buffer = BytesIO()
        self.pf = ReadPacketizedFilestream(self.mem_buffer)
        self.enc = enc_func(self.pf,flanking_primer5+primer5,flanking_primer3+primer3,bIndex=self.beginIndex)
        self.size=0

    def encode_segments_header(self,segments):
        oprimer5 = segments[0][2]
        oprimer3 = segments[0][3]
        assert len(segments) <= 256
        #if len(segments) == 0:
        #    return []
        hdr = [ len(segments) ]
        logger.debug("encode segments header : {}".format(hdr))
        for s in segments:
            hdr += convertIntToBytes(s[0],2)
            hdr += encode_size_and_value( s[1] )
            hdr += encode_size_and_value( s[4] )
            hdr += encode_primer_diff(oprimer5,s[2])
            hdr += encode_primer_diff(oprimer3,s[3])

        return hdr

    def encode_segments_header_comments(self,segments):
        comment = "% segment descriptions\n"        
        tab = "%    "
        for i,s in enumerate(segments):
            comment += tab + "{}. ".format(i) + file_system_format_description(s[0]) + "\n"
            comment += tab + "    size = {}".format(s[1]) + "\n"
            comment += tab + "    5' = {}".format(s[2]) + "\n"
            comment += tab + "    3' = {}".format(s[3]) + "\n"
            comment += tab + "    beginIndex = {}".format(s[4]) + "\n"            
            comment += "%\n"
        return comment

    def close(self):
        logger.debug("WriteSegmentedDNAFile.close")
        
        self.flush()
        self._record_segment() # record last segment
        
        hdr_other = self.encode_segments_header(self.segments)

        formatid = file_system_formatid_by_abbrev("Segmented")

        #print "formatid=",formatid
        
        size = sum([x[1] for x in self.segments])
        primer5 = self.segments[0][2]
        primer3 = self.segments[0][3]
        flanking5 = self.segments[0][-2]
        flanking3 = self.segments[0][-1]
        
        hdr = encode_file_header(self.output_filename,formatid,size,hdr_other,flanking5+primer5,flanking3+primer3,fsmd_abbrev=self.fsmd_abbrev)

        for i,h in enumerate(hdr):
            self.strands.insert(i,h)

        comment = encode_file_header_comments(self.output_filename,formatid,\
                                              size,hdr_other,primer5,primer3)
        self.out_fd.write(comment)
        comment = self.encode_segments_header_comments(self.segments)
        self.out_fd.write(comment)
        for ss in self.strands:
            if type(ss) is list:
                for s in ss:
                    self.out_fd.write("{}\n".format(s))
            else:
                self.out_fd.write("{}\n".format(ss))
                    
        if self.out_fd != sys.stdout and self.out_fd != sys.stderr:
            self.out_fd.close()
        return

class SegmentedReadDNAFilePipeline(ReadDNAFilePipeline):

    def decode_segments_header(self,other_data):
        val = other_data
        numSeg = val[0]
        if numSeg==0:
            return
        pos = 1
        allSegs = []
        for i in range(numSeg):
            # get format of this segment
            seg = [convertBytesToInt(val[pos:pos+2])]
            pos+=2

            # get size in bytes
            v,p = decode_size_and_value(val,pos)
            pos += p
            seg += [v]

            # get begin index
            v,p = decode_size_and_value(val,pos)
            pos += p
            seg += [v]

            primer5,p = decode_primer_diff(val[pos:], self.primer5)

            pos += p
            primer3,p = decode_primer_diff(val[pos:], self.primer3)

            pos += p
            seg.append(primer5)
            seg.append(primer3)
            allSegs.append(seg)

        return allSegs



    # ReadDNAFile reads a set of strands from a file.  It finds the header,
    # determines compatibility and encoding type, and then decodes the file.
    #    
    def __init__(self,**kwargs):     
        ReadDNAFilePipeline.__init__(self,**kwargs)

        logger.debug("sizeof other_data = {}".format(len(self.header['other_data'])))
        
        if len(self.header['other_data'])==0:
            return

        # restore cursor to end of buffer for writing
        self.mem_buffer.seek(0,2)

        #print self.header['other_data']
        
        segs = self.decode_segments_header(self.header['other_data'])
        self.segments = segs

        #print "segments=",segs

        for s in segs:
            logger.debug("formatid={} size={} bindex={} primer5={} primer3={}".format(s[0],s[1],s[2],s[3],s[4]))                                             
            formatid = s[0]
            size = s[1]
            bindex = s[2]
            primer5 = s[3]
            primer3 = s[4]
            if 'use_single_primer' in kwargs and kwargs['use_single_primer']==True:
                # strands from sequencing should all have the same primer
                primer5 = self.primer5
                primer3 = self.primer3

            dec_func = file_system_decoder(formatid)
            #self.mem_buffer = BytesIO()
            self.pf = WritePacketizedFilestream(self.mem_buffer,size,\
                                                file_system_format_packetsize(formatid),\
                                                minKey=bindex)

            #print primer5, primer3, bindex
            self.dec = dec_func(self.pf,primer5,primer3,bindex)

            for s in self.strands:
                if s.find(primer5)!=-1:
                    #print ("dnafile.py",self.dec.decode_from_phys_to_strand(s))
                    self.dec.decode(s)

            self.dec.write()
            write_anyway = ('write_incomplete_file' in kwargs) and \
                kwargs['write_incomplete_file']==True
            #if self.dec.complete or write_anyway:
            #    self.dec.write()
            if not write_anyway:
                assert self.dec.complete
            
        #print [x for x in self.mem_buffer.getvalue()]
        self.mem_buffer.seek(0,0) # set read point at beginning of buffer        
        return

    def read(self, n=1):        
        return self.mem_buffer.read(n)
    def readline(self, n=-1):
        return self.mem_buffer.readline(n)

    def readable(self):
        return True
    def writable(self):
        return False

    
if __name__ == "__main__":
    import os
    import tempfile
    encoder_params={
        "primer3":"A"*10,
        "primer5":"T"*10,
        "blockSizeInBytes":150,
        "strandSizeInBytes":15,
        "innerECC":2,
        "outerECCStrands":20,
        "dna_length":300,
    }

    temp_header_pipeline_data=tempfile.NamedTemporaryFile(mode="wb+",delete=False)

    wf = DNAFilePipeline.open("out.dna","w",format_name='Pipe-RS+CFC8',fsmd_abbrev='FSMD-Pipe',encoder_params=encoder_params)
    
    wf._header_fd=temp_header_pipeline_data
    
    for i in range(1000):
        wf.write( bytearray(convertIntToBytes(i,4)) )

    wf.close()


    rf = DNAFilePipeline.open("out.dna","r",format_name='Pipe-RS+CFC8',fsmd_abbrev='FSMD-Pipe',encoder_params=encoder_params,fsmd_header_filename=temp_header_pipeline_data.name)

    print ("Should print out 0 to 1000: ")
    while True:
        s = rf.read(4)
        if len(s)==0:
            break
        n = convertBytesToInt([x for x in s])
        print (n)

    print ("Done.")
    os.remove("out.dna")
    sys.exit(0)
    
        

