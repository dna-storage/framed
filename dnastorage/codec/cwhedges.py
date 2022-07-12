from dnastorage.codec.base_conversion import *
from dnastorage.codec.base import *
from dnastorage.codec_types import *
from math import log2,floor
import bitarray

from dnastorage.codec.hedges import hedges_state

import dnastorage.codec.codewordhedges as codewordhedges


class CodewordHedgesPipeline(BaseCodec,CWtoDNA):
    def __init__(self,codebook,guess_limit=100000,CodecObj=None,Policy=None):
        self._bits_per_cw = math.floor(log2(len(codebook)))
        assert self._bits_per_cw <=32 and "Codewords are limited to representing 32 bits at most, see C++ implementation if it needs to be changed"
        self._codebook = {}
        #take the first 2**(self._bits_per_cw) codewords from the codebook
        for key,DNA in sorted(codebook.items(),key= lambda x: x[0]):
            if key >= 2**(self._bits_per_cw): break
            self._codebook[key]=DNA
        self._hedges_state = hedges_state(rate=self._bits_per_cw,seq_bytes=0,pad_bits=0,prev_bits=0)
        CWtoDNA.__init__(self)
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        #initialize codebooks
        codewordhedges.codebook_init(self._codebook,"codewords")
        self._guess_limit=guess_limit
        
    def _encode(self,strand):
        self._hedges_state.set_message_bytes(len(strand.codewords))
        adjustment_bits = (self._bits_per_cw - 8*len(strand.codewords)%self._bits_per_cw)#need adjustment bits for padding to have the right size for codewords
        if adjustment_bits == self._bits_per_cw: adjustment_bits=0
        total_bits = adjustment_bits+8*len(strand.codewords)
        #we're just going to encode here, should be fast enough w/o c++ to just do codeword lookups
        b_array = bitarray.bitarray(endian = "little") #little endian is what the c code uses
        b_array.frombytes(bytearray(strand.codewords))
        b_array = b_array + bitarray.bitarray('0',endian="little")*adjustment_bits
        assert len(b_array) == total_bits and "Bit arrray length mismatch"
        #go through the bit array and get codewords
        codeword_list=[]
        for i in range(0,len(b_array),self._bits_per_cw):
            start_bit = i
            end_bit = i+self._bits_per_cw
            bits = b_array[start_bit:end_bit]
            key = int.from_bytes(bits.tobytes(),"little") #need to use little endian to be consistent with the rest of the C++ library
            codeword_list.append(self._codebook[key])
        strand.dna_strand = "".join(codeword_list)
        return strand

    def _decode(self,strand):
        strand.codewords = codewordhedges.decode(strand.dna_strand,self._hedges_state,self._guess_limit)
        return strand

    def __del__(self):
        codewordhedges.codebook_destroy("codewords") #destroy the codebook that is stored as a c-pointer

    #store some pertinent information like bit lengths of data seen to be able to reinstantiate the decoder in a correct state    
    def _encode_header(self):
        data = []
        data+=convertIntToBytes(self._hedges_state.message_bytes,4)
        return data
    
    def _decode_header(self,buff):
        pos=0
        message_bytes = convertBytesToInt(buff[pos:pos+4])
        pos+=4
        self._hedges_state.set_message_bytes(message_bytes)
        return buff[pos:]



if __name__ == "__main__":
    import random
    from dnastorage.strand_representation import *
    #test case for codeword hedges
    test_bytes = [random.randint(0,255) for _ in range(0,100)]
    '''
    test_codebook = {
        0: "AGAGAACT",
        1: "TCAGCTTT",
        2: "TCATTTTT",
        3: "ATATTAAA",
        4: "TGAAAAAA",
        5: "GAAAAAAA",
        6: "CTATATAA",
        7: "CTATAGAA"
    }'''

    from commafreecodec import cfc_all

    test_codebook={}
    for i in range(0,len(cfc_all)): #using cfc codebook to test
        test_codebook[i]=cfc_all[i]

    print(len(test_codebook))
    print(test_codebook)
    cwhedge  = CodewordHedgesPipeline(test_codebook)

    test_DNA = BaseDNA(codewords=test_bytes)

    print("DNA Bytes {}".format(test_DNA.codewords))
    
    cwhedge.encode(test_DNA)
    
    print("DNA after encoding {}".format(test_DNA.dna_strand))
    
    cwhedge.decode(test_DNA)

    print("Bytes after decoding {}".format(test_DNA.codewords))

    assert(test_bytes == list(test_DNA.codewords) and "Error Bytes don't match")
