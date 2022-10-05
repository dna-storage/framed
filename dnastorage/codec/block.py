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

block_count = 0

def reportBlockStatus(blocks,bindex,interIndexSize,intraIndexSize):
    global block_count
    for j,b in enumerate(blocks):
        block_count+=1
        contents = [ '_' for i in range(255) ]
        for s in b[1]:
            intra = base_conversion.convertBytesToInt(s[interIndexSize:interIndexSize+intraIndexSize])
            contents[intra] = '*'
        
        #print b[0]
        #stats["reportBlockStatus::blocks({}<{}>)".format(b[0],block_count)] = b[1]
        stats["reportBlockStatus({})::block_size({}<{}>)".format(bindex,b[0],block_count)] = len(b[1])
        stats["reportBlockStatus({})::block_profile({}<{}>)".format(bindex,b[0],block_count)] = "".join(contents)


def partitionStrandsIntoBlocks(strands, interIndexSize=2):
    logger.debug("partitionStrandsIntoBlocks:strands={}".format(len(strands)))
    blocksD = {}
    for s in strands:
        idx = base_conversion.convertBytesToInt(s[0:interIndexSize])
        blocksD[idx] = blocksD.get(idx,[]) + [s]
    blocks = blocksD.items()
    #blocks.sort()
    return blocks



def doMajorityVote(strands, indexBytes=3):
    key_value = {}
    strand_array = []
    for s in strands:
        #print "get key,value for strand {}".format(s)
        key = base_conversion.convertBytesToInt(s[0:indexBytes])
        #print value
        if key is not None:
            key_value[key] = key_value.get(key,[]) + [s]
            
    for key in key_value:
        data=[]
        if len(key_value[key])==1:
            # no reason to vote on a single strand
            strand_array.append(key_value[key][0])
            continue
        #print "processing key values"
        mx = max([len(data_strand) for data_strand in key_value[key]])
        for x in range(0,mx):
            #get a list of values that belong to the same location
            # if there's a smaller strand, due to error, this assert
            # will fire before attempting to access an illegal index
            same_position_values = []
            for data_strand in key_value[key]:
                if x < len(data_strand):
                    same_position_values.append(data_strand[x])
            #print same_position_values
            cnt=Counter(same_position_values)
            most_common_value=cnt.most_common(1)[0][0]
            #print most_common_value
            data.append(most_common_value)
            #put the key and data into an array
        strand_array.append(data)
        
    return strand_array


class NormalizeBlock(BaseCodec):
    def __init__(self, blockSizeInBytes):
        self.blockSize = blockSizeInBytes
        def _encode(self, packet):
            if len(packet[0]) < self.blockSize:
                block = packet[0] + [0] * (len(packet[0])-self.blockSize)
            else:
                block = packet[0]
            return packet[0],block

        def _decode(self, packet):
            return packet

class BlockToStrand(BaseCodec):
    """ This class encodes an (index,block) pair into a strand that contains 
        it's own index follwed by the data.  After this conversion, the index
        simpy appears as part of the strand.
        The _encode returns a list of strands.
    
        The _decode function expects to receive a set of strands that belong to
        the same block. Then they are assembled into a cohesive block and returned.     
    """

    def __init__(self, strandSizeInBytes, blockSizeInBytes, intraIndexSize=1, interIndexSize=2, nSyms=256, CodecObj=None, Policy=None, filterZeroes=False):
        super(BlockToStrand,self).__init__(CodecObj=CodecObj,Policy=Policy)
        self._strandSizeInBytes = strandSizeInBytes
        self._blockSize = blockSizeInBytes
        self._interIndex = interIndexSize
        self._intraIndex = intraIndexSize
        self._nSyms = nSyms
        self._removeZeroes = filterZeroes
        logger.info("filterZeroes = {}".format(filterZeroes))
        
    def allzeroes(self, l):
        for z in l:
            if z != 0:
                return False
        return True

    def _filter_zeroes(self, block):
        if self._removeZeroes == False:
            return [ False for _ in range(0,len(block),self._strandSizeInBytes) ]

        zeroes = []
        for i in range(0,len(block),self._strandSizeInBytes): 
            strand = block[i:i+self._strandSizeInBytes]
            if self.allzeroes(strand):
                zeroes.append(True)
            else:
                zeroes.append(False)
                    
        fil = [ False for _ in range(len(zeroes)) ]
        assert len(fil) == len(zeroes)
        for i,z in enumerate(zeroes):
            if i > 0 and i < len(zeroes)-1:
                if zeroes[i-1]==True and zeroes[i+1]==True and zeroes[i]==True:
                    fil[i] = True

        #print fil
        return fil
                    
    def _encode(self, packet):
        strands = []
        bindex = base_conversion.convertIntToBytes(packet[0],self._interIndex)
        block = packet[1]
        #print block

        fil = self._filter_zeroes(block)
            
        #assert len(block) <= self._blockSize
        
        assert len(bindex) <= self._interIndex
        
        for i in range(0,len(block),self._strandSizeInBytes):
            #if allzeroes(block[i:i+self._strandSizeInBytes]):
            #    continue
            if fil[int(i/self._strandSizeInBytes)]==True: # don't emit zeroes
                stats.inc("BlockToStrand::filterAllZeroStrand")
                continue
            
            sindex = base_conversion.convertIntToBytes(i/self._strandSizeInBytes,self._intraIndex)
            assert len(sindex) <= self._intraIndex
            s = bindex + sindex  + block[i:i+self._strandSizeInBytes]            
            strands.append(s)
            
        return strands

    def _decode(self, packet):        
        """ accepts a list of strands and coverts them into a contiguous block """

        d = {} # track all strands we find

        stats.inc("BlockToStrand::_decode::numberOfBlocks")
        stats.append("BlockToStrand::_decode::indices",packet[0])
        
        prior_bindex = packet[0]
        strands = packet[1]
        max_index = 0
        for s in strands:
            bindex = base_conversion.convertBytesToInt(s[0:self._interIndex])
            sindex = base_conversion.convertBytesToInt(s[self._interIndex:self._interIndex+self._intraIndex])

            if bindex != prior_bindex:
                e = DNABlockBadIndex("Bad index {}!={}".format(bindex,prior_bindex))
                #print "Bad index {}!={}".format(bindex,prior_bindex)
                if self._Policy.allow(e):
                    # ignore this strand
                    if sindex in d:
                        # already found a strand with this sindex, don't use this one
                        continue
                    else:
                        pass
                else:
                    raise e

            if len(s[self._interIndex+self._intraIndex:]) != self._strandSizeInBytes:            
                # handle exception condition here!
                err = DNABlockPayloadWrongSize("Strand payload size should be {} bytes but was {}.".format(self._strandSizeInBytes,len(s[self._interIndex+self._intraIndex:])))
                if self._Policy.allow(err):
                    if len(s[self._interIndex+self._intraIndex:]) > self._strandSizeInBytes:
                        d[sindex] = s[self._interIndex+self._intraIndex:self._interIndex+self._intraIndex+self._strandSizeInBytes]
                    else:
                        d[sindex] = s[self._interIndex+self._intraIndex:]+[-1]*(self._strandSizeInBytes-len(s[self._interIndex+self._intraIndex:]))
                else:
                    raise err
            else:
                d[sindex] = s[self._interIndex+self._intraIndex:]
                
        data = []
        max_index = int(self._blockSize / self._strandSizeInBytes)
        contents = [ '_' for _ in range(max_index) ]
        for i in range(max_index):
            if not i in d:
                if self._removeZeroes:
                    if not (i-1 in d) and not ((i+1) in d):
                        # guess that it's all zeros
                        stats.inc("BlockToStrand::decode::guessAllZeroStrand")
                        data += [ 0 for _ in range(self._strandSizeInBytes) ]
                        contents[i] = '0'
                    elif (i-1) in d and self.allzeroes(d[i-1]) and not (i+1 in d):
                        # guess that it's all zeros
                        stats.inc("BlockToStrand::decode::guessAllZeroStrand")
                        data += [ 0 for _ in range(self._strandSizeInBytes) ]
                        contents[i] = '0'
                    elif (i+1) in d and self.allzeroes(d[i+1]) and not (i-1 in d):
                        # guess that it's all zeros
                        stats.inc("BlockToStrand::decode::guessAllZeroStrand")
                        data += [ 0 for _ in range(self._strandSizeInBytes) ]
                        contents[i] = '0'
                    else:
                        err = DNABlockMissingIndex("Block missing index = {}".format(i))
                        stats.inc("BlockToStrand::decode::guessMissingStrand")
                        if self._Policy.allow(err):
                            #print "Block missing index = {}".format(i)
                            data += d.get(i,[-1 for _ in range(self._strandSizeInBytes)])
                            contents[i] = 'M'
                        else:
                            raise err
                else:
                    err = DNABlockMissingIndex("Block missing index = {}".format(i))
                    if self._Policy.allow(err):
                        #print "Block missing index = {}".format(i)
                        data += d.get(i,[-1 for _ in range(self._strandSizeInBytes)])
                    else:
                        raise err
            else:
                contents[i] = '-'
                data += d[i]
        stats.unique("BlockToStrand::_decode::block_profile({})".format(prior_bindex),"".join(contents))
                
        # smaller is allowed due to end of file, larger is a potential problem
        if len(data) > self._blockSize:
            err = DNABlockTooLargeError("Block too large, is {} should be {}".format(len(data),self._blockSize))
            if self._Policy.allow(err):                
                data = data[:self._blockSize]
            else:
                raise err

        #print data
        #stats.append("BlockToStrand::_decode::data[{}]".format(prior_bindex),data)
        
        return prior_bindex,data
        

class DoNothingOuterCodec(BaseCodec):
    """This is a simple do nothing Outer Codec to bypass this step when testing new designs"""
    def __init__(self,packetSize,payloadSize,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._packetSize=packetSize
        self._payloadSize=payloadSize
    def _encode(self,packet):
        """simply check length of packet to make sure it is a multiple of payload size"""
        index = packet[0]
        data = packet[1][:]
        
        assert len(packet[1])<=self._packetSize

        # hopefully the last packet if it's not a full packetSize
        # so, pad it out with zeros to make an even number of strands if
        # isn't already
        
        if len(data) < self._packetSize:
            stats.inc("DoNothingOuterCodec.padPacket")
            # normalize to multiple of payloadSize
            rem = len(data) % self._payloadSize
            data += [0]*(self._payloadSize-rem)

        #print "Data length {}".format(len(data))
        assert len(data)==self._packetSize
        return packet[0],data



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
    
class ReedSolomonOuterCodec(BaseCodec):
    """ReedSolomonOuterCodec takes a block of data as input and produces a new
    block of data that's error encoded. 

    The encoding input is a linear array of length payloadSize * nCol bytes,
    where nRow is ideally set to the payload per strand. However, it's
    not a strict requirement.

    Bcol must be less than rs.field_charac, which is either 2**8 or 2**16.

    """
    def __init__(self,packetSize,errorSymbols,payloadSize,c_exp=8,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)

        self._rs = get_reed_solomon(c_exp=c_exp)

        self._packetSize = packetSize
        self._errorSymbols = int(errorSymbols)
        self._payloadSize = payloadSize
        
        self._lengthMessage = packetSize / payloadSize + errorSymbols

        assert(self._lengthMessage <= self._rs.field_charac) # requirement of RS
        
    """
    packet is a tuple with the key at position 0 and value at position 1. This should return
    the Reed-Solomon encoded byte array.
    """
    def _encode(self,packet):

        index = packet[0]
        data = packet[1][:]
        
        assert packet[0] < 256**3
        assert len(packet[1])<=self._packetSize

        # hopefully the last packet if it's not a full packetSize
        # so, pad it out with zeros to make an even number of strands if
        # isn't already
        
        if len(data) < self._packetSize:
            stats.inc("RSOuterCodec.padPacket")
            # normalize to multiple of payloadSize
            rem = len(data) % self._payloadSize
            data += [0]*(self._payloadSize-rem)

        #print data
        assert len(data) % self._payloadSize == 0
        rows = len(data) / self._payloadSize
        # construct message

        ecc = []        
        for i in range(self._payloadSize):
            message = data[i:len(data):self._payloadSize]
            mesecc = self._rs.rs_encode_msg(message, self._errorSymbols)
            ecc += mesecc[len(message):]
            #print "encode message=",mesecc
            #print message, ecc
            #print  len(message), len(ecc), self._errorSymbols
            #assert len(ecc) == self._errorSymbols

        tranpose = []
        for i in range( int(self._errorSymbols) ):
            syms = ecc[i:len(ecc):self._errorSymbols]
            data += syms
                    
        # separate out key and value
        # convert mesecc into a string
        return packet[0],data

    """
    This function expects a list of unsigned integers in GF(256). For now, erasures are
    denoted with -1.
    """
    def _decode(self,packet):
        #print packet[1]
        data = [x for x in packet[1]]
        #print data
        rows = len(data) / self._payloadSize
        for i in range(self._payloadSize):
            message = data[i:len(data):self._payloadSize]            
            erasures = [ j for j in range(len(message)) if message[j]==-1 ]
            try:
                # correct the message
                corrected_message, corrected_ecc = \
                    self._rs.rs_correct_msg(message, \
                                            self._errorSymbols, \
                                            erase_pos=erasures)
                #debug help
                #print corrected_message
                #print corrected_ecc
                data[i:len(data):self._payloadSize] = corrected_message + corrected_ecc
                stats.inc("RSOuterCodec::correct")
            except ReedSolomonError as e:
                #print "couldn't correct block {}".format(message)
                stats.inc("RSOuterCodec::ReedSolomonError")
                # wrap exception into a library specific one when checking policy:
                wr_e = DNAReedSolomonOuterCodeError(msg=\
                             "RSOuterCodec found error at index={}".format(packet[0]))
                if self._Policy.allow(wr_e):
                    # fix -1
                    corrected_message = [ max(0,_) for _ in message ]
                    data[i:len(data):self._payloadSize] = corrected_message
                    pass # nothing to do
                else:
                    # policy doesn't allow handling, raise exception
                    raise wr_e
                # just proceed without further error correction
                pass
            except Exception as e:
                if self._Policy.allow(e):
                    pass
                else:
                    raise e

        # discard outer correction codes
        data = data[0:len(data)-self._errorSymbols*self._payloadSize]        
        return packet[0],data


if __name__ == "__main__":
    import sys
    from random import randint
    from dnastorage.exceptions import NoTolerance,AllowAll

    b2s = BlockToStrand(20,80,Policy=AllowAll())

    x = [ randint(0,255) for x in range(4*20) ]
    print (x)
    
    r = b2s.encode( [2323,x] )
    print (r)

    k = randint(0,len(r)-2)

    r[3] = r[3]+[32,235, 22]
    
    index,y = b2s.decode ( r[:k]+r[k+1:] )
    print (index, y)
    
    print (sum([ (a-b)**2 for a,b in zip(x,y) ]))
    
