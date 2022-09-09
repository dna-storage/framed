from dnastorage.codec.base import *
from dnastorage.exceptions import *
from dnastorage.codec.reedsolomon.rs import ReedSolomon,get_reed_solomon,ReedSolomonError
from random import randint
from dnastorage.util.stats import stats
from dnastorage.codec_types import *
from dnastorage.strand_representation import *


import logging
logger = logging.getLogger('dna.storage.codec.strand')
logger.addHandler(logging.NullHandler())

class RandomizeCodec(BaseCodec):
    def __init__(self,numRandBytes,CodecObj=None,Policy=None):
        super(RandomizeCodec,self).__init__(CodecObj=CodecObj,Policy=Policy)
        self._numRandBytes = numRandBytes

    def _encode(self, array):
        for i in range(self._numRandBytes):
            array.append( randint(0,255) )

    def _decode(self, array):
        return array[:self._numRandBytes]


class ReedSolomonInnerCodec(BaseCodec):
    """
    ReedSolomonInnerCodec takes a key,value pair as input to the _encode function and
    produces a Reed-Solomon encoded message as a byte array. Both key and value are
    protected. If the value is shorter than expected, it is padded to padWidth with
    random data.

    This is an "Inner" Codec because it only can correct errors within a strand.

    This class is hard coded to use GF(256).
    """
    def __init__(self,numberECCBytes,c_exp=8,CodecObj=None,Policy=None):
        super(ReedSolomonInnerCodec,self).__init__(CodecObj=CodecObj,Policy=Policy)

        self.rs = get_reed_solomon(c_exp=c_exp)
        self._numberECCBytes = numberECCBytes


    """ Accepts list/array of bytes. Return the Reed-Solomon encoded byte array.
    """
    def _encode(self,array):
        try:
            assert len(array) <= self.rs.field_charac
        except:
            #print "Failed RS check"
            return array
        # usually not needed, but makes the code a little more tolerant for use with
        # a variety of codecs that may pass in strings, bytearrays, or lists:
        message = [x for x in array]

        try:
            # encoded the message using the RS library
            mesecc = self.rs.rs_encode_msg(message, self._numberECCBytes)
        except ReedSolomonError as e:
            raise DNACodingError("Error while encoding Reed-Solomon Inner Codec.")
        except ZeroDivisionError as e:
            pass
        
        #print "RSInner:",len(mesecc)
        return mesecc

    """
    This function expects a list of unsigned integers in the GF. For now, erasures are
    denoted with -1.
    """
    def _decode(self,array):
        #print array
        message = [x for x in array] 
        # find the -1s in the list
        erasures = [ i for i in range(len(message)) if (message[i]==-1 or message[i]==None) ]
        # correct the message
        if self._numberECCBytes == 0:
            return array
        
        try:
            corrected_message, corrected_ecc = self.rs.rs_correct_msg(message,self._numberECCBytes, erase_pos=erasures)
            value = corrected_message
        except ReedSolomonError as e:
            if self._Policy.allow(e):
                # leave erasures, may be helpful for outer decoder
                value = message[0:-(self._numberECCBytes)]
                pass # nothing to do
            else:
                #print (str(e))
                raise DNACodingError("RSInnerCodec failed to correct message.")
            # just proceed without further error correction
            pass
        except ZeroDivisionError as e:
            #stats.inc("RSInnerCodec.ZeroDivision")
            if self._Policy.allow(e):
                # leave erasures, may be helpful for outer decoder
                value = message[0:(self._numberECCBytes)]
                pass # nothing to do
            else:
                #print (str(e))
                raise DNACodingError("RSInnerCodec failed to correct message.")

        return value



#Wrapper class so that the reedsolomon inner codec can be used with the pipeline infrastructure
class ReedSolomonInnerCodecPipeline(ReedSolomonInnerCodec,CWtoCW):
    def __init__(self,numberECCBytes,c_exp=8,CodecObj=None,Policy=None):
        ReedSolomonInnerCodec.__init__(self,numberECCBytes,c_exp=c_exp,CodecObj=CodecObj,Policy=Policy)
        CWtoCW.__init__(self)
    def _encode(self,strand):
        strand.codewords=ReedSolomonInnerCodec._encode(self,strand.codewords)
        return strand
    def _decode(self,strand):
        strand.codewords=ReedSolomonInnerCodec._decode(self,strand.codewords)
        return strand

class CRC8(BaseCodec,CWtoCW):
    def __init__(self,CodecObj=None,Policy=None):
        CWtoCW.__init__(self)
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)

        self._memo_array = [0]*256
        self._polynomial = 0x131
        self._checksum = 0 
        #initialize memor array
        for i in range(0,256):
            crc=i
            for j in range(0,8)[::-1]:
                if(0x80 & crc): crc = ((crc<<1)&0xff) ^ (self._polynomial&0xff)
                else:  crc = (crc<<1)&0xff
            self._memo_array[i]=crc

    def _crc(self,byte_array):
        crc = 0
        for b in byte_array:
            crc=crc^b
            crc = self._memo_array[crc]
        return crc
            
    def _encode(self,strand):
        crc = self._crc(strand.codewords)
        strand.codewords.append(crc)
        return strand
    
    def _decode(self,strand):
        if None in strand.codewords:
            strand.codewords=[None]*(len(strand.codewords)-1) #remove strand from consideration
            return strand
        crc = self._crc(strand.codewords)
        if crc!=self._checksum:
            strand.codewords=[None]*(len(strand.codewords)-1) #remove strand from consideration
        else:
            strand.codewords=strand.codewords[0:len(strand.codewords)-1]
        return strand

class CRC8_Index(CRC8,CWtoCW): #Instead of CRCing the whole strand, we just CRC the index to make sure we don't miss-place indices, useful for matching strands correctly from sequencing experiments
    def __init__(self,CodecObj=None,Policy=None):
        CWtoCW.__init__(self)
        CRC8.__init__(self)

    def _encode(self,strand): #insert CRC close to index as possible
        crc = self._crc(strand.codewords[0:strand.index_bytes])
        strand.codewords = strand.codewords[0:strand.index_bytes] + [crc] + strand.codewords[strand.index_bytes::]
        return strand

    def _decode(self,strand):
        if None in strand.codewords[0:strand.index_bytes+1]:
            strand.codewords=[None]*(len(strand.codewords)-1) #remove strand from consideration
            return strand
        crc = self._crc(strand.codewords[0:strand.index_bytes+1])
        #logger.info("CRC {} index+crc {}".format(crc,strand.codewords[0:strand.index_bytes+1]))
        #logger.info("Entire strand {}".format(strand.codewords))
        if crc!=self._checksum:
            strand.codewords=[None]*(len(strand.codewords)-1) #remove strand from consideration
        else:
            strand.codewords=strand.codewords[0:strand.index_bytes] +strand.codewords[strand.index_bytes+1::] #remove embedded crc
        return strand


    
if __name__=="__main__":
    import copy
    codewords=[5,6,7,8,9,214]
    original_codewords = copy.copy(codewords)
    s = BaseDNA(codewords=codewords)
    crc = CRC8()
    crc.encode(s)
    print(s.codewords)
    crc.decode(s)
    print(s.codewords)

    assert(s.codewords==original_codewords)
    
    print("Testing Index CRC")
    crc = CRC8_Index()

    s.index_bytes = 2
    crc.encode(s)
    print(s.codewords)
    crc.decode(s)
    print(s.codewords)

    assert(s.codewords==original_codewords)
    
    
