#!/usr/bin/python
from copy import copy
import random
from dnastorage.codec.base_conversion import *
from dnastorage.primer.primer_util import *
from dnastorage.codec.base import *

ibases = ['A', 'T', 'C']

ivalues = { 'A' : 0 ,
           'C' : 2 ,
           'G' : 3 ,  
           'T' : 1   }

def rangeWithExclusion(length, primer):
    excluded = []
    plist = [b for b in primer]
    for i in range(len(primer)):
        b = copy(ibases)
        if plist[i] != 'G':
            b.remove(plist[i])
        excluded.append(b)    
    size = 1
    for i in range(length):
        if i < len(primer):
            size *= len(excluded[i])
        else:
            size *= len(ibases)
    return size


def encodeWithExclusionHelper(dec,s,excluded):
    if len(s) < len(excluded): 
        b = len(excluded[len(s)])
        m = dec % b
        q = dec / b
        s += excluded[len(s)][m]
    else:
        m = dec % 3
        q = dec / 3
        s += ibases[m]

    if q > 0:
        return encodeWithExclusionHelper(q,s,excluded)
    else:
        return s

def encodeWithExclusion(dec,length,primer):
    excluded = []
    plist = [b for b in primer]
    for i in range(len(primer)):
        b = copy(ibases)
        if plist[i] != 'G':
            b.remove(plist[i])
        excluded.append(b)
    #print excluded
    s = encodeWithExclusionHelper(dec,"",excluded)
    if len(s) < len(primer):
        for i in range(len(s),min(length,len(primer))):
            s += excluded[i][0]
    s = s.ljust(length,ibases[0])
    return s


def decodeWithExclusion(s,primer):
    val = 0L
    power = 1L
    excluded = []
    plist = [b for b in primer]
    for i in range(len(primer)):
        b = copy(ibases)
        if plist[i]!='G':
            b.remove(plist[i])
        excluded.append( { b[i] : i for i in range(len(b)) }  )
    #print excluded
    base = 3
    for i in range(len(s)):
        if i < len(excluded):
            base = len(excluded[i])
            val += excluded[i][s[i]] * power
        else:
            base = 3
            val += ivalues[s[i]] * power
        power *= base
    return val


class IllinoisCodec(BaseCodec):
    def __init__(self,primer,numberBytes,CodecObj=None,keyWidth=20):
        BaseCodec.__init__(self,CodecObj)
        self._keyWidth=keyWidth
        self._primer = primer
        self._numberBytes = numberBytes

    def _encode(self,packet):
        key = encodeWithExclusion(packet[0],self._keyWidth,self._primer)
        array = bytearray(packet[1])
        #print [x for x in array]
        assert len(array) == self._numberBytes
        i = 0
        enc = []
        while i < len(array): 
            val = long(array[i])
            if i+1 < len(array):
                val += long(array[i+1]) << 8
            if i+2 < len(array):
                val += long(array[i+2]) << 16
            if i+3 < len(array):
                val += long(array[i+3]) << 24
            #print val,encodeWithExclusion(val,26,self._primer)
            enc.append(encodeWithExclusion(val,26,self._primer))
            i += 4                
        return key+"".join(enc)

    def _decode(self,s):
        key = decodeWithExclusion(s[0:self._keyWidth],self._primer)        
        i = 0
        rest = s[self._keyWidth:]
        assert len(rest)%26==0
        value = []
        while i < len(rest):
            val = decodeWithExclusion(rest[i:i+26],self._primer)             
            #print rest[i:i+26],val
            for _ in range(4):
                value.append(val&255)
                val = long(val >> 8)
            i += 26            
        assert len(value) >= self._numberBytes
        #print [x for x in bytearray(value[0:self._numberBytes])]
        return key,bytearray(value[0:self._numberBytes])


if __name__ == "__main__":
    import math
    import sys
    primer = "AGGTCGGACAACGCCTTAAG"
    N = rangeWithExclusion(26,primer)
    print N, math.log(N,2)
    sys.exit(0)
    for i in range(1,2**20):
        r = random.randint(0,N)
        s = encodeWithExclusion(r,24,primer)
        if correlation_distance(primer,s) > 3:
            show_correlation(primer,s)
            assert False
        print "{} = {}  ({})".format(r,s,decodeWithExclusion(s,primer))
        assert r == decodeWithExclusion(s,primer)
