from dnastorage.codec.base_conversion import *
from dnastorage.codec.base import *
from dnastorage.codec_types import *
from dnastorage.codec.codebooks import *
from math import log2,floor
import bitarray
import bitarray.util as bit_util

from dnastorage.codec.hedges import hedges_state

import dnastorage.codec.codewordhedges as codewordhedges

import logging
logger = logging.getLogger("dnastorage.codec.cwhedges")
logger.addHandler(logging.NullHandler())



class CodewordHedgesPipeline(BaseCodec,CWtoDNA):
    def __init__(self,codebook,syncbook=None,sync_period=0, parity_period=0, pad_bits=0, parity_history = 0, custom_reward=-0.31, guess_limit=100000,CodecObj=None,Policy=None):
        self._bits_per_cw = math.floor(log2(len(codebook)))
        logger.info("CWHEDGES BITS PER CW {}".format(self._bits_per_cw))
        logger.info("CWHEDGES Codebook Size {}".format(len(codebook)))
        assert self._bits_per_cw <=32 and "Codewords are limited to representing 32 bits at most, see C++ implementation if it needs to be changed"
        self._codebook = {}
        self._syncbook = syncbook
        #take the first 2**(self._bits_per_cw) codewords from the codebook
        for key,DNA in sorted(codebook.items(),key= lambda x: x[0]):
            if key >= 2**(self._bits_per_cw): break
            self._codebook[key]=DNA
        self._hedges_state = hedges_state(rate=self._bits_per_cw,seq_bytes=0,pad_bits=pad_bits,prev_bits=0,sync_period=sync_period,
                                          parity_period=parity_period,parity_history=parity_history,custom_reward=custom_reward)
        CWtoDNA.__init__(self)
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        #initialize codebooks
        codewordhedges.codebook_init(self._codebook,"codewords")
        if syncbook!=None: codewordhedges.syncbook_init(self._syncbook)
        self._guess_limit=guess_limit
        if self._hedges_state.parity_period>0 and self._hedges_state.parity_period==1: assert 0 and "parity_period should be set greater than 1 if !=0"


    def _set_parity(self,b_array): #NOTE: you can extend this to actually do different parities as well
        if self._hedges_state.parity_period==0: return b_array #no parity to add
        index=0
        while index!=len(b_array):
            if (index)%self._hedges_state.parity_period==0 and index!=0: #place parity bits at regular intervals
                start_bit = index-self._hedges_state.parity_history
                if index<self._hedges_state.parity_history or self._hedges_state.parity_history==0:
                    start_bit = 0
                parity_bit = bitarray.bitarray(str(bit_util.parity(b_array[start_bit:index])),endian="little")
                b_array = b_array[0:index]+parity_bit+b_array[index::]
            index+=1
        return b_array

 
    def _encode(self,strand):
        #we're just going to encode here, should be fast enough w/o c++ to just do codeword lookups
        b_array = bitarray.bitarray(endian = "little") #little endian is what the c code uses
        b_array.frombytes(bytearray(strand.codewords))
        b_array = b_array + bitarray.bitarray('0',endian = "little")*self._hedges_state.pad_bits
        b_array = self._set_parity(b_array)
        #set hedges state
        self._hedges_state.set_message_bytes(len(strand.codewords))
        
        #need adjustment bits for padding to have the right size for codewords
        adjustment_bits = (self._bits_per_cw - len(b_array)%self._bits_per_cw)
        if adjustment_bits == self._bits_per_cw: adjustment_bits=0
        b_array = b_array + bitarray.bitarray('0',endian="little")*adjustment_bits

        #go through the bit array and get codewords
        codeword_list=[]
        for i in range(0,len(b_array),self._bits_per_cw):
            start_bit = i
            end_bit = i+self._bits_per_cw
            bits = b_array[start_bit:end_bit]
            key = int.from_bytes(bits.tobytes(),"little") #need to use little endian to be consistent with the rest of the C++ library
            codeword_list.append(self._codebook[key])
            
        final_cw_list=[]
        if self._syncbook!=None and self._hedges_state.cw_sync_period!=0: #bake in the synchronization points for decoding simulation purposes
            sync_counter=0
            for i in range(0,len(codeword_list)):
                final_cw_list.append(codeword_list[i])
                if (i+1)%self._hedges_state.cw_sync_period==0:
                    final_cw_list.append(self._syncbook[sync_counter%len(self._syncbook)])
                    sync_counter+=1
        else:
            final_cw_list=codeword_list
        strand.dna_strand = "".join(final_cw_list)
        return strand

    def _decode(self,strand):
        decode_result= codewordhedges.decode(strand.dna_strand,self._hedges_state,self._guess_limit)
        strand.codewords = decode_result["return_bytes"]
        print("Score {}".format(decode_result["score"]))
        return strand

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
    from commafreecodec import cfc_all

    #Test out basic codeword, no parity or syncbooks
    test_codebook=CFC_DUMMY()
    
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
    
    #Test out the decoder with synbooks
    test_DNA = BaseDNA(codewords=test_bytes)
    syncbook=TEST_SYNC_BOOK()
    cwhedge = CodewordHedgesPipeline(test_codebook,syncbook=syncbook,sync_period=5)
    cwhedge.encode(test_DNA)
    print("DNA Bytes {}".format(test_DNA.codewords))
    print("DNA after encoding with synchronization points: {}".format(test_DNA.dna_strand))
    cwhedge.decode(test_DNA)
    print("Bytes after decoding with synchronization points: {}".format(test_DNA.codewords))
    assert(test_bytes == list(test_DNA.codewords) and "Error Bytes don't match for synchronization points")


    #Test out the decoder with parity
    print("Test bytes {}".format(test_bytes))
    test_DNA = BaseDNA(codewords=test_bytes)
    cwhedge = CodewordHedgesPipeline(test_codebook,parity_period=3,parity_history=0)
    cwhedge.encode(test_DNA)
    print("DNA Bytes {}".format(test_DNA.codewords))
    print("DNA after encoding with parity: {}".format(test_DNA.dna_strand))
    cwhedge.decode(test_DNA)
    print("Bytes after decoding with parity: {}".format(test_DNA.codewords))
    assert(test_bytes == list(test_DNA.codewords) and "Error Bytes don't match for parity")


    #Test out the decoder with parity and history
    print("Test bytes {}".format(test_bytes))
    test_DNA = BaseDNA(codewords=test_bytes)
    cwhedge = CodewordHedgesPipeline(test_codebook,parity_period=3,parity_history=8)
    cwhedge.encode(test_DNA)
    print("DNA Bytes {}".format(test_DNA.codewords))
    print("DNA after encoding with parity and history: {}".format(test_DNA.dna_strand))
    cwhedge.decode(test_DNA)
    print("Bytes after decoding with parity and history: {}".format(test_DNA.codewords))
    assert(test_bytes == list(test_DNA.codewords) and "Error Bytes don't match for parity and history")
