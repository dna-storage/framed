from dnastorage.codec.base import *
from math import log, ceil
from dnastorage.codec import base_conversion
from dnastorage.codec.reedsolomon.rs import ReedSolomon,get_reed_solomon,ReedSolomonError
from collections import Counter
from dnastorage.strand_representation import *
from dnastorage.codec_types import *
import logging
logger = logging.getLogger('dna.storage.codec.block')
logger.addHandler(logging.NullHandler())


class ReedSolomonOuterPipeline(BaseOuterCodec):
    def __init__(self,packet_divisor,parity_packets,c_exp=8,OuterCodecObj=None,Policy=None):
        BaseOuterCodec.__init__(self,packet_divisor,OuterCodecObj=OuterCodecObj,Policy=Policy)
        self._rs = get_reed_solomon(c_exp=c_exp)
        self._parity_packets=parity_packets
        assert(packet_divisor+parity_packets<=self._rs.field_charac)
    def _encode(self,packets):
        parity_packets=[]
        #initialize parity_packets
        for i in range(0,self._parity_packets):
            new= []
            parity_packets.append(new)
        #encode a set of packets
        for i in range(0,len(packets[0])):
            strand_set=[]
            for j in range(0,len(packets)):
                strand_set.append(packets[j][i])
            #use strand set to calculate a message
            assert(len(strand_set)>0)
            strand_length = len(strand_set[0].codewords)
            for i in range(0,self._parity_packets): parity_packets[i].append(BaseDNA(codewords=[])) #make a new strand in each parity packet
            #go column by column, calculating a RS message and splitting the ECC across strands in each packet
            for k in range(0,strand_length):
                message=[]
                for j in range(0,len(strand_set)):
                    message.append(strand_set[j].codewords[k]) #get bytes
                assert(len(message)>0)
                mesecc = self._rs.rs_encode_msg(message, self._parity_packets)
                ecc = mesecc[len(message):]
                for ei, e in enumerate(ecc):
                    parity_packets[ei][-1].codewords.append(e)
            for p in parity_packets:
                assert len(p[-1].codewords)==strand_length
        return packets+parity_packets
    
    def _decode(self,packets):
        #decode a set of packets 
        #Need to get packets into the right order to unroll the encoding, basically will get messages the same way as encoding
        assert len(packets)>0
        for i in range(0,len(packets[0])):
            strands=[]
            for j in range(0,self._total_sub_packets):
                strands.append(packets[j][i])
            #now get messages
            assert len(strands)>0
            for byte_index in range(0,len(strands[0].codewords)):
                message=[]
                for strand_index in range(0,len(strands)):
                    message.append(strands[strand_index].codewords[byte_index])
                assert len(message)==self._total_sub_packets
                try:
                    #perform correction
                    erasures = [i for i in range(0,len(message)) if message[i]==None]
                    corrected_message, corrected_ecc = \
                        self._rs.rs_correct_msg(message, \
                                                self._parity_packets, \
                                                erase_pos=erasures)
                except ReedSolomonError as e:
                    #just gonna assume we can't fix this, hopefully someone else will
                    #stats.inc("RSOuterCodec::ReedSolomonError")
                    wr_e = DNAReedSolomonOuterCodeError(msg=\
                             "RSOuterCodec found error at sub-packet index".format(strands[0].index_ints[:self._level]))
                    corrected_message = message
                except Exception as e:
                    if self._Policy.allow(e):
                        corrected_message = message                        
                        pass
                    else:
                        raise e

                #write the corrected message back into strand classes
                for s,m in zip(strands,corrected_message):
                    s.codewords[byte_index] = m
        return packets
    
