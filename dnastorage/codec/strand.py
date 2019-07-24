from dnastorage.codec.base import *
from dnastorage.exceptions import *
from dnastorage.codec.reedsolomon.rs import ReedSolomon,get_reed_solomon,ReedSolomonError
from random import randint

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

        assert len(array) <= self.rs.field_charac

        # usually not needed, but makes the code a little more tolerant for use with
        # a variety of codecs that may pass in strings, bytearrays, or lists:
        message = [x for x in array]

        try:
            # encoded the message using the RS library
            mesecc = self.rs.rs_encode_msg(message, self._numberECCBytes)
        except ReedSolomonError, e:
            raise DNACodingError("Error while encoding Reed-Solomon Inner Codec.")
        except ZeroDivisionError, e:
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
        erasures = [ i for i in range(len(message)) if message[i]==-1 ]
        # correct the message
        try:
            corrected_message, corrected_ecc = self.rs.rs_correct_msg(message,self._numberECCBytes, erase_pos=erasures)
            value = corrected_message
        except ReedSolomonError, e:
            stats.inc("RSInnerCodec.ReedSolomonError")
            if self._Policy.allow(e):
                # leave erasures, may be helpful for outer decoder
                value = message[0:(self._numberECCBytes)]
                pass # nothing to do
            else:
                print str(e)
                raise DNACodingError("RSInnerCodec failed to correct message.")
            # just proceed without further error correction
            pass
        except ZeroDivisionError, e:
            stats.inc("RSInnerCodec.ZeroDivision")
            if self._Policy.allow(e):
                # leave erasures, may be helpful for outer decoder
                value = message[0:(self._numberECCBytes)]
                pass # nothing to do
            else:
                print str(e)
                raise DNACodingError("RSInnerCodec failed to correct message.")

        return value