from dnastorage.codec.base_conversion import convertIntToBytes,convertBytesToInt
import editdistance as ed
from io import BytesIO
from dnastorage.util.packetizedfile import *
import math
import struct
from dnastorage.system.formats import *
import hashlib
def version_0_1():
    #returns a formatting dictionary for the given version
    version_dict={}
    version_dict["major"]=(0,1,int) #(value,number of bytes)
    version_dict["minor"]=(1,1,int) 
    version_dict["size"]=(None,8,int) #None values are determined from user input
    version_dict["filename"]=(None,"v",str) #assume filenames no larger than 50 bytes
    version_dict["formatid"]=(None,2,int)
    #version_dict["magic"]=('ATCGATGC',"v",str)
    version_dict["magic"]=('ATCGATGCGACGAGGA',"v",str)
    version_dict["decoding_format"]=('FSMD-Pipe',"v",str)
    return version_dict



global versional_formats
version_formats={"0.1": version_0_1}


class Header(object): 
    def __init__(self, version_number,primer5="",primer3=""):
        assert version_number in version_formats
        self._version=version_number
        self._format_dict= version_formats[self._version]()
        enc_func = file_system_encoder_by_abbrev(self._format_dict["decoding_format"][0])
        self._pipeline = enc_func(None,magic=self._format_dict["magic"][0],primer3=primer3,primer5=primer5)

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
                data += [ ord(_) for _ in value ] + [0]
        data = bytearray(data)+bytearray(convertIntToBytes(len(other_data),2))+bytearray(other_data)
        
        dhash = hashlib.md5()
        dhash.update(data)
        self._encoded_header_hash = dhash.digest()
        self._header_length = len(data)
        
        
        pf = ReadPacketizedFilestream(BytesIO(data))
        self._pipeline.set_read_pf(pf)
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
            
    
    def decode_file_header(self,strands):
        b = BytesIO()
        fid = file_system_formatid_by_abbrev(self._format_dict["decoding_format"][0])
        packetsize = file_system_format_packetsize(fid)
        pf = WritePacketizedFilestream(b,packetsize,packetsize)
        self._pipeline.set_write_pf(pf)
        
        self._non_header_strands=[]
        for s in strands:
            returned_strand = self._pipeline.decode(s)
            if returned_strand is not None:
                self._non_header_strands.append(returned_strand)
                continue
        #should be able to finish decoding here
        self._pipeline.final_decode()

        try:
            data = [ ord(x) for x in b.getvalue() ]
        except Exception as e:
            data = [ x for x in b.getvalue() ]


        if self._encoded_header_hash!=None: #quick check to determine if header survived correctly
            dhash = hashlib.md5()
            dhash.update(bytearray(data[:self._header_length]))
            if dhash.digest() != self._encoded_header_hash:
                return None
        
        pos = 0
        #now read out bytes
        return_header={}
        for field in self._format_dict:
            t = self._format_dict[field][2]
            length = self._format_dict[field][1]
            if length is "v":
                size_bytes = data[pos:pos+2]
                length = convertBytesToInt(size_bytes)
                pos+=2 #positions take two bytes
            value_bytes = data[pos:pos+length]
            pos+=length
            if t is int:
                return_header[field]=convertBytesToInt(value_bytes)
            elif t is str:
                return_header[field]="".join([chr(x) for x in value_bytes])

        #get other data
        size_other_data = convertBytesToInt(data[pos:pos+2])
        pos += 2
        return_header['other_data'] = [ x for x in data[pos:pos+size_other_data] ]
        
        return return_header
        
    def pick_nonheader_strands(self):
        #gets non-header strands, simply those not picked for decoding
        assert self._non_header_strands
        return self._non_header_strands


if __name__=="__main__":
    header =Header("0.1",primer5 = "T"*10, primer3 = "A"*10)

    enc_dict ={"size":400,"filename":"testing.filename",
               "formatid":0xbeef}

    strands = header.encode_file_header(enc_dict,[0,1,2])
    
    header=None

    header2 =Header("0.1",primer5 = "T"*10, primer3 = "A"*10)
    for s in strands:
        s.codewords=[]
        s.index_ints=()
        s.index_bytes=0
        
    decoded_header = header2.decode_file_header(strands)

    print(decoded_header)
  
