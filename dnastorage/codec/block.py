from dnastorage.codec.base import *
from math import log, ceil
from dnastorage.codec import base_conversion

def partitionStrandsIntoBlocks(strands, interIndexSize=2):
    blocksD = {}
    for s in strands:
        idx = base_conversion.convertBytesToInt(s[0:interIndexSize])
        blocksD[idx] = blocksD.get(idx,[]) + [s]
    blocks = blocksD.items()
    blocks.sort(cmp=lambda x,y: cmp(x[0],y[0]))
    return blocks


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

    def __init__(self, strandSizeInBytes, blockSizeInBytes, intraIndexSize=1, interIndexSize=2, nSyms=256, CodecObj=None, Policy=None):
        super(BlockToStrand,self).__init__(CodecObj=CodecObj,Policy=Policy)
        self._strandSizeInBytes = strandSizeInBytes
        self._blockSize = blockSizeInBytes
        self._interIndex = interIndexSize
        self._intraIndex = intraIndexSize
        self._nSyms = nSyms

    def allzeroes(self, l):
        for z in l:
            if z != 0:
                return False
        return True
        
    def _encode(self, packet):
        strands = []
        bindex = base_conversion.convertIntToBytes(packet[0],self._interIndex)
        block = packet[1]

        
        #assert len(block) <= self._blockSize
        
        assert len(bindex) <= self._interIndex
        
        for i in range(0,len(block),self._strandSizeInBytes):
            #if allzeroes(block[i:i+self._strandSizeInBytes]):
            #    continue

            sindex = base_conversion.convertIntToBytes(i/self._strandSizeInBytes,self._intraIndex)
            assert len(sindex) <= self._intraIndex
            s = bindex + sindex  + block[i:i+self._strandSizeInBytes]            
            strands.append(s)
            
        return strands

    def _decode(self, packet):        
        """ accepts a list of strands and coverts them into a contiguous block """
        d = {} # track all strands we find
        
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
                    if d.has_key(sindex):
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
        max_index = self._blockSize / self._strandSizeInBytes
        for i in range(max_index):
            if not d.has_key(i):
                err = DNABlockMissingIndex("Block missing index = {}".format(i))
                if self._Policy.allow(err):
                    #print "Block missing index = {}".format(i)
                    data += d.get(i,[-1 for _ in range(self._strandSizeInBytes)])
                else:
                    raise err
            else:
                data += d[i]
                
        # smaller is allowed due to end of file, larger is a potential problem
        if len(data) > self._blockSize:
            err = DNABlockTooLargeError("Block too large, is {} should be {}".format(len(data),self._blockSize))
            if self._Policy.allow(err):                
                data = data[:self._blockSize]
            else:
                raise err
            
        return prior_bindex,data
        


if __name__ == "__main__":
    import sys
    from random import randint
    from dnastorage.exceptions import NoTolerance,AllowAll

    b2s = BlockToStrand(20,80,Policy=AllowAll())

    x = [ randint(0,255) for x in range(4*20) ]
    print x
    
    r = b2s.encode( [2323,x] )
    print r

    k = randint(0,len(r)-2)

    r[3] = r[3]+[32,235, 22]
    
    index,y = b2s.decode ( r[:k]+r[k+1:] )
    print index, y
    
    print sum([ (a-b)**2 for a,b in zip(x,y) ])
    
