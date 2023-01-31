from dnastorage.system.formats import *
from dnastorage.system.header_class import *
from dnastorage.strand_representation import *
from dnastorage.util.stats import stats
from dnastorage.util.strandinterface import *
from dnastorage.util.mpi_utils import *
import io
import os
import sys
import copy
import traceback
import logging
logger = logging.getLogger('dna.storage.system.dnafile')
logger.addHandler(logging.NullHandler())

DATA_BARCODE=0xCC

def check_mpi(comm): #check that required modules are loaded for mpi communication
    if comm:
        if not ("mpi4py" in sys.modules or "mpi4py.MPI" in sys.modules):
            raise SystemError("mpi4py is not loaded, needed for pipeline mpi support")

class DNAFilePipeline:
    def __init__(self):
        return
    @classmethod
    def open(self,op, format_name="", write_incomplete_file=False, header_version="0.1",fsmd_header_filename=None,
             payload_header_filename = None,encoder_params={},header_params={},do_write = True, file_barcode=tuple(),
             mpi=None,strand_interface=None,dna_file_name = None,store_header=True):
        # check if we are reading or writing
        if op=="r":
            assert strand_interface!=None
            strands=[]
            check_mpi(mpi)
            strands = strand_interface.strands
            if mpi: #communicate out strands to different processes
                strands = object_scatter(strands,mpi)
            assert (strands is not None) and len(strands)>0
            header = Header(header_version,header_params,barcode_suffix=file_barcode,mpi=mpi)
            #initialize header pipeline
            try:
                with open(fsmd_header_filename,"rb") as serialized_header_pipeline_data:
                    header.set_pipeline_data(serialized_header_pipeline_data.read())
            except Exception as e:
                logger.fatal("Could not initialize header pipeline: {}".format(e))
                exit(1)
            #attempt to get header data
            try:
                h = header.decode_file_header(strands) 
                if mpi: #broadcast the actual decoded header
                    logger.info("Rank {} Broadcasting header".format(mpi.rank))
                    h = mpi.bcast(h,root=0)
                    header.set_header_dict(h)
                if h!=None: logger.info("Able to decode header from DNA")
                if h==None:
                    if not mpi or mpi.rank==0: stats.inc("dead_header",1) #only increment this in one rank, all other rank counts of this are redundant
                    raise ValueError("Header failed to decode, trying file")  
            except Exception as e:
                try:
                    traceback.print_exc()
                    logger.fatal("{}".format(e))
                    #now try to decode with bytes
                    with open(payload_header_filename,"rb") as serialized_payload_pipeline_data:
                        logger.info("using binary file for payload header instead of DNA")
                        h = header.header_from_bytes(serialized_payload_pipeline_data.read())
                    if h==None: raise ValueError("Header is None when attempting to decode from bytes")
                except Exception as e:
                    logging.warning("Could not recover payload header: {}".format(e))
                    return None

            logger.debug("decoded header: {}".format(h))

            if h['main_pipeline_formatid'] == file_system_formatid_by_abbrev("Segmented"):
                assert 0 and "Segmented files not supported at this moment with pipelines"
            else:
                return ReadDNAFilePipeline(encoder_params=encoder_params,
                                           header=header,file_barcode=file_barcode,mpi=mpi)
        elif op=="w":
            return WriteDNAFilePipeline(output=dna_file_name,
                                        format_name=format_name,encoder_params=encoder_params,
                                        header_version=header_version,header_params=header_params,fsmd_header_filename=fsmd_header_filename,
                                        payload_header_filename=payload_header_filename,file_barcode=file_barcode,do_write=do_write,store_header=store_header)
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

class ReadDNAFilePipeline(DNAFilePipeline):
    # ReadDNAFile reads a set of strands from a file.  It finds the header,
    # determines compatibility and encoding type, and then decodes the file.
    def __init__(self,**kwargs):
        DNAFilePipeline.__init__(self)
        self._enc_opts=kwargs["encoder_params"]
        mpi=kwargs.get("mpi",None) #allow mpi-parallel decoding
        self._file_barcode=kwargs.get("file_barcode",tuple())
        header_class = kwargs["header"] #header was already decoded just grab it
        self.header = header_class.header_dict() # store the header dictionary3
        self.strands = header_class.pick_nonheader_strands()
        self.formatid = self.header['main_pipeline_formatid']
 
        self.size = self.header['size']
        # set up mem_buffer 
        self.mem_buffer = BytesIO()
    
        if self.formatid == 0x1000:
            # let sub-classes handle initialization
            return
        constructor_function = file_system_decoder(self.formatid)  
        self.mem_buffer = BytesIO()
        self.pf = WritePacketizedFilestream(self.mem_buffer,self.size,0)
        self.pipe = constructor_function(self.pf,**self._enc_opts,barcode=(DATA_BARCODE,)+self._file_barcode)
        self.pipe.decode_header_data(self.header["other_data"])
        self.pipe.mpi = mpi #attach the communicator to the pipeline
        for s in self.strands:
            self.pipe.decode(s)
        self.pipe.final_decode()
        if mpi: logger.info("Rank {} leaving dnafile".format(mpi.rank))
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
    def get_returned_strands(self):
        #return strands that were kicked out of the pipeline(s)
        try:
            return self.pipe.get_filtered()
        except Exception as e:
            logger.fatal("Issue with getting pipeline filter strands: {}".format(e))
            exit(1)
    
class WriteDNAFilePipeline(DNAFilePipeline):
    # WriteDNAFile writes a set of strands.
    def __init__(self,**kwargs):
        DNAFilePipeline.__init__(self)
        self._enc_opts=kwargs["encoder_params"]
        self._header_params=kwargs["header_params"]
        self._header_version=kwargs["header_version"]
        self._file_barcode = kwargs["file_barcode"]
        self._do_write = kwargs.get("do_write",True)
        self.out_fd=None
        self._header_fd=None
        self._payload_header_fd = None
        self.output_filename="none.dna"
        self._store_header=kwargs.get("store_header",True)
        
        if "fsmd_header_filename" in kwargs and kwargs["fsmd_header_filename"]!=None:
            self._header_fd=open(kwargs["fsmd_header_filename"],"wb+")

        if "payload_header_filename" in kwargs and kwargs["payload_header_filename"]!=None:
            self._payload_header_fd = open(kwargs["payload_header_filename"],"wb+")
            
        if 'formatid' in kwargs:            
            enc_func = file_system_encoder(kwargs['formatid'])
            self.formatid = kwargs['formatid']
        elif 'format_name' in kwargs:
            enc_func = file_system_encoder_by_abbrev(kwargs['format_name'])
            self.formatid = file_system_formatid_by_abbrev(kwargs['format_name'])

        self.mem_buffer = BytesIO()
        self.pf= ReadPacketizedFilestream(self.mem_buffer)

        self.pipe = enc_func(self.pf,**self._enc_opts,barcode=(DATA_BARCODE,)+self._file_barcode)
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
        header = Header(self._header_version,self._header_params,barcode_suffix=self._file_barcode)

        header_dict={}
        header_dict["filename"]=os.path.basename(self.output_filename)
        header_dict["main_pipeline_formatid"]=self.formatid
        header_dict["size"]=self.size
        logger.info("Right before encoding the header block")
        hdr_strands= header.encode_file_header(header_dict,self.pipe.encode_header_data())
        
        if self._payload_header_fd!=None:
            self._payload_header_fd.write(header.encoded_header_bytes) #write down the bytes for the header, useful for debugging purposes 
            self._payload_header_fd.close()
            
        if self._store_header: #allow headers to be skipped from being stored
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
        
        if self._do_write:
            for s in self.strands:
                self.out_fd.write("{}:{}\n".format(s.index_ints,s.dna_strand))
            
        if self.out_fd != sys.stdout and self.out_fd != sys.stderr:
            self.out_fd.close()
        return
    
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

    wf = DNAFilePipeline.open("out.dna","w",format_name='Pipe-RS+CFC8',encoder_params=encoder_params,file_barcode=(1,2,3))
    
    wf._header_fd=temp_header_pipeline_data
    
    for i in range(1000):
        wf.write( bytearray(convertIntToBytes(i,4)) )

    wf.close()


    rf = DNAFilePipeline.open("out.dna","r",format_name='Pipe-RS+CFC8',encoder_params=encoder_params,fsmd_header_filename=temp_header_pipeline_data.name,file_barcode=(1,2,3))

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
    
        

