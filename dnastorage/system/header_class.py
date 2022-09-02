from dnastorage.codec.base_conversion import convertIntToBytes,convertBytesToInt
import editdistance as ed
from io import BytesIO
from dnastorage.util.packetizedfile import *
import math
import struct
from dnastorage.system.formats import *
import hashlib
from importlib_metadata import version
import sys

import logging
logger = logging.getLogger("dnastorage.system.header_class")
logger.addHandler(logging.NullHandler())

def version_0_1():
    #returns a formatting dictionary for the given version
    version_dict={}
    version_dict["major"]=(0,1,int) #(value,number of bytes)
    version_dict["minor"]=(1,1,int) 
    version_dict["size"]=(None,8,int) #None values are determined from user input
    version_dict["filename"]=(None,"v",str) #assume filenames no larger than 50 bytes
    version_dict["main_pipeline_formatid"]=(None,2,int)
    version_dict["header_barcode"]=(0xEE,1,int)
    version_dict["pipeline_barcode_ID"]=(None,"v",bytes) 
    version_dict["decoding_format"]=('FSMD-Pipe',"v",str)
    version_dict["dnastorage_module_version"]=(version("dnastorage"),"v",str) #version of core decoding code, note pre-defined values only matter during decoding
    return version_dict



def version_0_2():
    #returns a formatting dictionary for the given version
    version_dict={}
    version_dict["major"]=(0,1,int) #(value,number of bytes)
    version_dict["minor"]=(2,2,int) 
    version_dict["size"]=(None,8,int) #None values are determined from user input
    version_dict["filename"]=(None,"v",str) #assume filenames no larger than 50 bytes
    version_dict["main_pipeline_formatid"]=(None,2,int)
    version_dict["header_barcode"]=(0xEE,1,int)
    version_dict["pipeline_barcode_ID"]=(None,"v",bytes)
    version_dict["decoding_format"]=('CustomSDC',"v",str)
    version_dict["dnastorage_module_version"]=(version("dnastorage"),"v",str)
    return version_dict


def version_0_3():
    #returns a formatting dictionary for the given version
    version_dict={}
    version_dict["major"]=(0,1,int) #(value,number of bytes)
    version_dict["minor"]=(3,2,int) 
    version_dict["size"]=(None,8,int) #None values are determined from user input
    version_dict["filename"]=(None,"v",str) #assume filenames no larger than 50 bytes
    version_dict["main_pipeline_formatid"]=(None,2,int)
    version_dict["header_barcode"]=(0xEE,1,int)
    version_dict["pipeline_barcode_ID"]=(None,"v",bytes)
    version_dict["decoding_format"]=('CustomPipe-RS+Codeword+Hedges',"v",str)
    version_dict["dnastorage_module_version"]=(version("dnastorage"),"v",str)
    return version_dict


def version_0_4():
    #returns a formatting dictionary for the given version
    version_dict={}
    version_dict["major"]=(0,1,int) #(value,number of bytes)
    version_dict["minor"]=(4,2,int) 
    version_dict["size"]=(None,8,int) #None values are determined from user input
    version_dict["filename"]=(None,"v",str) #assume filenames no larger than 50 bytes
    version_dict["main_pipeline_formatid"]=(None,2,int)
    version_dict["header_barcode"]=(0xEE,1,int)
    version_dict["pipeline_barcode_ID"]=(None,"v",bytes)
    version_dict["decoding_format"]=('CustomPipe-RS+CFC8',"v",str)
    version_dict["dnastorage_module_version"]=(version("dnastorage"),"v",str)
    return version_dict

def version_0_5():
    #returns a formatting dictionary for the given version
    version_dict={}
    version_dict["major"]=(0,1,int) #(value,number of bytes)
    version_dict["minor"]=(5,2,int) 
    version_dict["size"]=(None,8,int) #None values are determined from user input
    version_dict["filename"]=(None,"v",str) #assume filenames no larger than 50 bytes
    version_dict["main_pipeline_formatid"]=(None,2,int)
    version_dict["header_barcode"]=(0xEE,1,int)
    version_dict["pipeline_barcode_ID"]=(None,"v",bytes)
    version_dict["decoding_format"]=('BasicHedges',"v",str)
    version_dict["dnastorage_module_version"]=(version("dnastorage"),"v",str)
    return version_dict



global versional_formats
version_formats={"0.1": version_0_1,
                 "0.2": version_0_2,
                 "0.3": version_0_3,
                 "0.4": version_0_4,
                 "0.5": version_0_5
}


class Header(object): 
    def __init__(self, version_number,encoder_params,barcode_suffix=tuple(),mpi=None):
        assert version_number in version_formats
        self._version=version_number
        self._format_dict= version_formats[self._version]()
        enc_func = file_system_encoder_by_abbrev(self._format_dict["decoding_format"][0])
        barcode = (self._format_dict["header_barcode"][0],)+barcode_suffix
        self._pipeline = enc_func(None,**encoder_params,barcode=barcode)
        self._format_dict["pipeline_barcode_ID"]=(barcode_suffix,)+self._format_dict["pipeline_barcode_ID"][1::]
        self.encoded_header_bytes = None 
        self._mpi = mpi
        
    def get_header_pipeline_data(self):
        buff=self._pipeline.encode_header_data()
        buff+=self._encoded_header_hash
        buff+=convertIntToBytes(self._header_length,4)
        return buff
    
    def set_pipeline_data(self,buff):
        buff=self._pipeline.decode_header_data(buff)
        self._encoded_header_hash=buff[:hashlib.md5().digest_size]
        pos = hashlib.md5().digest_size
        self._header_length=convertBytesToInt(buff[pos:pos+4])
        
    def encode_file_header(self,encode_dict,other_data,**kwargs):
        data=[]
        for field in self._format_dict:
            t = self._format_dict[field][2]
            if self._format_dict[field][0]==None:
                assert field in encode_dict
                value = encode_dict[field]
            else:
                value=self._format_dict[field][0]
            assert value is not None
            if t is int:
                data+=convertIntToBytes(value,self._format_dict[field][1])
            elif t is str:
                if self._format_dict[field][1]=="v":
                    data += convertIntToBytes(len(value)+1,2)
                if self._format_dict[field][2] is str:
                    data += [ ord(_) for _ in value ] + [0]
            elif t is bytes:
                data += convertIntToBytes(len(value),2)
                data += bytearray(value)
                    
        data = bytearray(data)+bytearray(convertIntToBytes(len(other_data),2))+bytearray(other_data)
        
        dhash = hashlib.md5()
        dhash.update(data)

        self._encoded_header_hash = dhash.digest()
        self._header_length = len(data)

        b  = BytesIO(data)
        pf = ReadPacketizedFilestream(b)
        self._pipeline.set_read_pf(pf)
        self.encoded_header_bytes = b.getvalue() #bytes that were encoded as the header for a pipeline
        strands = []
        
        for block in self._pipeline:
            for s in block:
                strands.append(s)
        return strands

    def encode_file_header_comments(self,encode_dict):
        comment=""
        for field in self._format_dict:
            if self._format_dict[field][0] != None:
                value = self._format_dict[field][0]
            else:
                value = encode_dict[field]
            comment += "% {} : {}\n".format(field,value)
        return comment


    def header_from_bytes(self,data): #returns the header dictionary from bytes
        logger.info("Header length {} expected length {}".format(len(data),self._header_length))
        if self._encoded_header_hash!=None: #quick check to determine if header survived correctly
            dhash = hashlib.md5()
            dhash.update(bytearray(data[:self._header_length]))
            if dhash.digest() != self._encoded_header_hash:
                return None
        pos = 0
        #now read out bytes
        return_header={}
        for field in self._format_dict:
            is_variable=False
            t = self._format_dict[field][2]
            length = self._format_dict[field][1]
            if length=="v":
                is_variable=True
                size_bytes = data[pos:pos+2]
                length = convertBytesToInt(size_bytes)
                pos+=2 #positions take two bytes
            value_bytes = data[pos:pos+length]
            pos+=length
            if t is int:
                return_header[field]=convertBytesToInt(value_bytes)
            elif t is bytes and is_variable:
                return_header[field]=tuple([int(x) for x in value_bytes])
            elif t is str and is_variable:
                return_header[field]="".join([chr(x) for x in value_bytes])
        #get other data
        size_other_data = convertBytesToInt(data[pos:pos+2])
        pos += 2
        return_header['other_data'] = [ x for x in data[pos:pos+size_other_data] ]
        self._decode_header=return_header
        return return_header
    
    def decode_file_header(self,strands):
        b = BytesIO()
        pf = WritePacketizedFilestream(b,self._header_length,0)
        self._pipeline.set_write_pf(pf)
        self._non_header_strands=[]
        self._pipeline.mpi=self._mpi
        for s in strands:
            self._pipeline.decode(s)
        #should be able to finish decoding here
        self._pipeline.final_decode()
        try:
            data = [ ord(x) for x in b.getvalue() ]
        except Exception as e:
            data = [ x for x in b.getvalue() ]
        if self._mpi!=None and self._mpi.Get_rank()!=0: return dict() 
        return self.header_from_bytes(data)
        
    def pick_nonheader_strands(self):
        return self._pipeline.get_filtered()

    def set_header_dict(self,h):
        self._decode_header=h

    def header_dict(self):
        if self._decode_header !=None:
            return self._decode_header
        else:
            return None


if __name__=="__main__":
   

    enc_dict ={"size":400,"filename":"testing.filename",
               "main_pipeline_formatid":0xbeef}

    encoding_params={"primer5":"T"*10,
                     "primer3":"A"*10
                     }

    header =Header("0.1",encoding_params,(0,1,2))
    
    strands = header.encode_file_header(enc_dict,[0,1,2])
    
    header=None

    header2 =Header("0.1",encoding_params,(1,2,3))
    for s in strands:
        s.codewords=[]
        s.index_ints=()
        s.index_bytes=0
        
    decoded_header = header2.decode_file_header(strands)

    print(decoded_header)
  
