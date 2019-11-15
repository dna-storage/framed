#!/usr/bin/python
from dnastorage.primer.primer_util import *
from dnastorage.primer import nextera
from dnastorage.codec import dense
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import base_conversion
from dnastorage.codec import fountain
from dnastorage.codec.codecfile import *
from dnastorage.codec.base import *
from dnastorage.util.packetizedfile import *
from copy import *

import logging
logger = logging.getLogger('dna.storage.arch.strand')
logger.addHandler(logging.NullHandler())


class Checksum(BaseCodec):
    def __init__(self,nbytes=1,CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        self._nbytes = nbytes

    def _encode(self,packet):
        key = packet[0]
        value = [x for x in bytearray(packet[1])]
        chk = convertIntToBytes(sum(value) % 256, self._nbytes)
        value = value + chk
        # double check to make sure decoder will work
        assert convertBytesToInt(value[-self._nbytes:]) == sum(value[:-self._nbytes])% 256
        value = str(bytearray(value))
        return (key,value)

    def _decode(self,s):
        key = s[0]
        value = [x for x in bytearray(s[1])]
        if convertBytesToInt(value[-self._nbytes:]) == sum(value[:-self._nbytes])% 256:
            return (key,str(bytearray(value[:-self._nbytes])))
        #dont break on a mismatching checksum
        else:
            return None,None

class CheckLength(BaseCodec):
    def __init__(self,nbytes=1):
        BaseCodec.__init__(self,None)
        self._nbytes = nbytes

    def _encode(self,packet):
        key = packet[0]
        value = [x for x in bytearray(packet[1])]
        assert len(value)==self._nbytes
        return packet

    def _decode(self,s):
        if s[1] == None: # check for end
            return s
        value = [x for x in bytearray(s[1])]
        if len(value) > self._nbytes:
            #JMT: Should we just truncate? or someting smarter?
            value = value[0:self._nbytes]
            return (s[0],str(bytearray(value)))
        if len(value) < self._nbytes:
            # JMT: fix me
            # should we just take the first ones? or something smarter?
            return (s[0],str(bytearray(value[0:self._nbytes])))
        return s


class StrandPrimers(BaseCodec):
    def __init__(self,bprimer,eprimer,CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        self.begin_primer = bprimer
        self.end_primer = eprimer
        self.rc_end_primer = reverse_complement(eprimer)
    def _encode(self,packet):
        assert isinstance(packet,str)
        s = self.begin_primer + packet  + self.rc_end_primer
        return s
    def _decode(self,strand):
        assert isinstance(strand,str)
        if not (self.begin_primer in strand[0:len(self.begin_primer)]):
            #assert False and "Begin primer not found"
            logger.debug("Didn't find 5' primer: {} in {}".format(self.begin_primer,strand[0:len(self.begin_primer)]))
        if not (self.rc_end_primer in strand[-len(self.end_primer):]):
            #assert False and "End primer not found"
            logger.debug("Didn't find 3' primer: {} in {}".format(self.rc_end_primer,strand[-len(self.rc_end_primer):]))
        return strand[len(self.begin_primer):-len(self.end_primer)]

class EncodeNaiveStrand(EncodePacketizedFile):
    def __init__(self,packetizedFile,CodecObj):
        assert CodecObj!=None
        EncodePacketizedFile.__init__(self,packetizedFile,CodecObj)
        self.index = 0

    def __iter__(self):
        EncodePacketizedFile.__iter__(self)
        self.index = 0
        return self

    def _encode(self):
        #next packet from file
        value = self._packetizedFile.next()
        # get index
        key = self.index
        self.index+=1
        #return key,value
        return (key,value)

class DecodeNaiveStrand(DecodePacketizedFile):
    def __init__(self,packetizedFile,CodecObj):
        assert CodecObj!=None
        DecodePacketizedFile.__init__(self,packetizedFile,CodecObj)

    def _decode(self, key, value):
        # pass key,value to file writer
        DecodePacketizedFile.writeToFile(self,key,value)

class EncodeGoldmanStrand(EncodePacketizedFile):
    def __init__(self,packetizedFile,overlap,CodecObj):
        EncodePacketizedFile.__init__(self,packetizedFile,CodecObj)
        self.index = 0
        self._overlap = overlap

    def __iter__(self):
        self.index = 0
        return self

    def next(self):
        if self.index >= self._packetizedFile.numberOfPackets:
            raise StopIteration()
        packet = self.encode()
        return packet

    def _encode(self):
        mod = self._packetizedFile.numberOfPackets
        s = ""
        for i in range (self._overlap):
            s += self._packetizedFile[(self.index+i)%mod]
        key = self.index
        self.index+=1
        return (key,s)

class DecodeGoldmanStrand(DecodePacketizedFile):
    def __init__(self,packetizedFile,overlap,CodecObj):
        DecodePacketizedFile.__init__(self,packetizedFile,CodecObj)
        self._overlap = overlap


    def _decode(self, key, value):
        self._data[key] = value
        packetSize = self._packetizedFile.packetSize
        mod = self._packetizedFile.numberOfPackets
        for i in range (self._overlap):
            if not self._packetizedFile.has_key((key+i)%mod):
                self._packetizedFile[(key+i)%mod] = value[i*packetSize:(i+1)*packetSize]

    def write(self):
        assert self._packetizedFile.complete
        DecodePacketizedFile.write(self)
    #dummy write for Goldman encoding
    def dummy_write(self):
        return DecodePacketizedFile.dummy_write(self)

class EncodeXORStrand(EncodePacketizedFile):
    def __init__(self,packetizedFile,CodecObj):
        EncodePacketizedFile.__init__(self,packetizedFile,CodecObj)
        self.index = 0

    def __iter__(self):
        self.index = 0
        return self

        # 1 -> 1
        # 2 -> 3
        # 3 -> 5

    def next(self):
        if self.index >= 2 * self._packetizedFile.numberOfPackets - 1:
            raise StopIteration()
        packet = self.encode()
        return packet

    def _encode(self):
        if self.index % 2 == 0:
            s = self._packetizedFile[self.index/2]
            l = [x for x in bytearray(s)]
        else:
            s1 = self._packetizedFile[self.index/2]
            s2 = self._packetizedFile[self.index/2+1]
            l = [x^y for x,y in zip(bytearray(s1),bytearray(s2))]
        key = self.index
        self.index+=1
        return (key,str(bytearray(l)))

class DecodeXORStrand(DecodePacketizedFile):
    def __init__(self,packetizedFile,CodecObj):
        DecodePacketizedFile.__init__(self,packetizedFile,CodecObj)
        self._num_XOR_used=0
    def _decode(self, key, value):
        if key % 2 == 0:
            #self.decode(key,value)
            self._data[key] = value
            self._packetizedFile[key/2] = value
        else:
            self._data[key] = value

    def assemble(self,m1,m2):
        s1 = self._data[m1]
        s2 = self._data[m2]
        l = [x^y for x,y in zip(bytearray(s1),bytearray(s2))]
        return str(bytearray(l))

    def write(self):
        changes=True
        if not self._packetizedFile.complete:
            while changes:
                changes=False
                missing = self._packetizedFile.getMissingKeys()
                missing = [ 2*m for m in missing ]
                for m in missing:
                    assert self._data.has_key(m)==False
                    # strategy one: check m-1 and m-2
                    if self._data.has_key(m-1) and self._data.has_key(m-2):
                        changes=True
                        value = self.assemble(m-1,m-2)
                        self._data[m] = value
                        self._packetizedFile[m/2] = value
                        # strategy two: check m+1 and m+2
                    elif self._data.has_key(m+1) and self._data.has_key(m+2):
                        changes=True
                        value = self.assemble(m+1,m+2)
                        self._data[m] = value
                        self._packetizedFile[m/2] = value
                    else:
                        assert False and "Could not recover file."
            assert self._packetizedFile.complete
        DecodePacketizedFile.write(self)


    #dummy writes are needed for all architectures to be evaluated by FI tools
    def dummy_write(self):
        changes=True
        if not self._packetizedFile.complete:
            while changes:
                changes=False
                missing = self._packetizedFile.getMissingKeys()
                missing = [ 2*m for m in missing ]
                for m in missing:
                    assert self._data.has_key(m)==False
                    # strategy one: check m-1 and m-2
                    if self._data.has_key(m-1) and self._data.has_key(m-2):
                        changes=True
                        value = self.assemble(m-1,m-2)
                        self._data[m] = value
                        self._packetizedFile[m/2] = value
                        self._num_XOR_used+=1
                        # strategy two: check m+1 and m+2
                    elif self._data.has_key(m+1) and self._data.has_key(m+2):
                        changes=True
                        value = self.assemble(m+1,m+2)
                        self._data[m] = value
                        self._packetizedFile[m/2] = value
                        self._num_XOR_used+=1
                #print self._num_XOR_used
        return(DecodePacketizedFile.dummy_write(self))


class EncodeFountainStrand(EncodePacketizedFile):
    def __init__(self,packetizedFile,factor,CodecObj):
        EncodePacketizedFile.__init__(self,packetizedFile,CodecObj)
        self.index = 0
        self.fountain = fountain.UnlimitedFountain(self._packetizedFile.numberOfPackets)
        self.maxStrands = int(factor*self._packetizedFile.numberOfPackets)

    def __iter__(self):
        self.index = 0
        return self

    def getTable(self):
        return self.fountain.getTable()

    def next(self):
        if self.fountain.numStrands >= self._packetizedFile.numberOfPackets:
            can = self.fountain.canDecode()
            if can and self.fountain.numStrands >= self.maxStrands:
                raise StopIteration()
        packet = self.encode()
        return packet

    def _encode(self):
        assert False

    def encode(self):
        fails = 1
        while True:
                #self.index += 1
            enc1,strands = self.tryEncode()

            enc = self._Codec.encode(enc1)

            if fails % 10000 == 0:
                print ("10000 fails...")
                print (enc)

            if hasSingleRun(enc):
                fails+=1
                continue

            if not nextera_strand_comparison(enc,5):
                #print "Discarded strand..."
                self.index += 1
                fails += 1
                continue

            self.fountain.insert(self.index,strands)
            self.index+=1
            return enc

    def tryEncode(self):
        strands = []
        if self.fountain.numStrands >= self._packetizedFile.numberOfPackets:
            if not self.fountain.canDecode():
                m = self.fountain.getMissing()
                if len(m) > 0:
                    strands = m[0:1]
            else:
                strands = self.fountain.next()
        else:
            strands = self.fountain.next()

        l = bytearray(self._packetizedFile[strands[0]])
        if len(strands) > 1:
            for s in strands[1:]:
                l2 = self._packetizedFile[ s ]
                xor = [x^y for x,y in zip(bytearray(l),bytearray(l2))]
                l = bytearray(xor)

        return (self.index,str(bytearray(l))),strands

class DecodeFountainStrand(DecodePacketizedFile):
    def __init__(self,packetizedFile,table,CodecObj):
        DecodePacketizedFile.__init__(self,packetizedFile,CodecObj)
        self._lookup = { x[0] : copy(x[1]) for x in table}

    def _updateAllData(self):
        visit = [ i for i in range(self._packetizedFile.numberOfPackets) ]
        while len(visit)>0:
            nextVisit = []
            for i in visit[:]:
                if self._packetizedFile.has_key(i):
                    val = self._packetizedFile[i]
                    nextVisit += self._updateData(i,val)
            visit = nextVisit

    def _updateData(self,index,value):
        visited = []
        for item in copy(self._lookup.items()):
            if self._data.has_key(item[0]):
                for i in copy(item[1]):
                    if i == index:
                        val = bytearray(self._data[item[0]])
                        val2 = bytearray(value)
                        xor = [x^y for x,y in zip(bytearray(val),bytearray(val2))]
                        self._data[item[0]] = bytearray(xor)
                        item[1].remove(i)
                if len(item[1])==1:
                    self._packetizedFile[item[1][0]] = self._data[item[0]]
                    self._lookup.pop(item[0])
                    visited.append(item[1][0])
        return visited

    def _decode(self, key, value):
        self._data[key] = value
        assert self._lookup.has_key(key)
        l = self._lookup[key]
        if len(l)==1:
            self._packetizedFile[l[0]] = value
            self._updateData(l[0],value)

    def write(self):
        if not self._packetizedFile.complete:
            self._updateAllData()
            assert self._packetizedFile.complete
        DecodePacketizedFile.write(self)


if __name__ == "__main__":
    import os
    import sys
    from random import randint


    illini = illinois.IllinoisCodec('AGGTCGGACAACGCCTTAAG',150,Checksum())

    d = dense.DenseCodec(Checksum())
    #h = huffman.HuffmanRotateCodec(huffman.HuffmanCodec(23,Checksum()))

    #enc = EncodeFountainStrand(sys.argv[1],22,1.5,h)

    b = binary.BinaryRotateCodec(binary.BinaryCodec(Checksum()))
    h = huffman.HuffmanRotateCodec(huffman.HuffmanCodec(21,Checksum()))

    #enc = EncodeGoldmanStrand(sys.argv[1],5,4,b)
    #dec = DecodeGoldmanStrand("out.d",enc.fileSize,5,4,b)

    #enc = EncodeXORStrand(sys.argv[1],22,h)
    #dec = DecodeXORStrand("out.d",enc.fileSize,22,h)

    enc = EncodeXORStrand(sys.argv[1],149,illini)
    dec = DecodeXORStrand("out.d",enc.fileSize,149,illini)

    i = 0
    strands = []
    for s in enc:
        strands.append(s)
        print ("{}. ({}) - {}".format(i,len(s),s))
        #if i % 5 != 0 or i == 0 or i > 500:
        #    dec.decodeStrand(s)
        i += 1

    #dec = DecodeFountainStrand("out.d",enc.fileSize,22,enc.getTable(),h)
    for s in strands:
        dec.decode(s)
    dec.write()

    #enc = EncodeXORHuffmanStrand(sys.argv[1],136)
    #enc1 = EncodeXORHuffmanStrand("out.d",136)
    #i = 0
    #for x,y in zip(enc,enc1):
    #    if x != y:
    #        print i,x
    #        print i,y
    #    i += 1
