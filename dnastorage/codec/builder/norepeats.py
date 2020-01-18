#!/usr/bin/python
from dnastorage.primer.design import *
from dnastorage.codec.base_conversion import *
import random
import sys
from math import log, trunc
from dnastorage.codec.builder.codecbuilder import TableCodecBuilder


def _makeNoRepeatStrandsHelper(oList,allowedRepeats=0):
    nList = []
    for o in oList:
        for b in ['A', 'C', 'G', 'T']:
            copy = [ _ for _ in o ]
            copy.append(b)
            if repetition(copy)<=allowedRepeats:
                nList.append("".join(copy))
    return nList

def makeNoRepeatStrands(length=5,allowedRepeats=0):
    good = []
    for i in range(4**min(5,length)):
        #print "in here!"
        dd = convertBase(4,i,min(5,length))
        if repetition(dd) <= allowedRepeats:
            good.append(dd)
    #cw = 4**trunc(log(len(good),4))
    #print size,len(good), cw, log(cw,4)/size
    if length <= 5:
        return good
    else:
        for _ in range(6,length+1):
            good = _makeNoRepeatStrandsHelper(good,allowedRepeats)
        return good

def repetition(s,verbose=False):
    r = sum([ int(a==b) for a,b in zip([_ for _ in s[:-1]],[_ for _ in s[1:]]) ])
    if verbose:
        zipped = zip([_ for _ in s[:-1]],[_ for _ in s[1:]])
        idx = [ int(a==b) for a,b in zipped ]
        ss =  " "
        for i,x in enumerate(idx):
            if x==1:
                ss = ss+'^'
            else:
                ss = ss+' '
        print (s)
        print (ss)
    return r


def build_norepeats_codec(filename,class_name):
    list = makeNoRepeatStrands(10)
    final = []
    for l in list:
        if not (l[0]=='T' and l[-1]=='T') and not(l[0]=='G' and l[-1]=='G'):
            final.append(l)
    random.shuffle(final)
    final = final[0:2**16]
    assert len(final) == 2**16
    t = TableCodecBuilder(filename,class_name,final)
    t.build()
    return

if __name__ == "__main__":
    for i in range(1,11):
        l = makeNoRepeatStrands(i)
        #print i,l[0:5],len(l)
        cw = 4**trunc(log(len(l),4))
        print ("{:3d} {:8d} {:8d} {:3.2f}".format(i,len(l), cw, log(cw,4)/i))
        #print l

    list = makeNoRepeatStrands(10)
    hist = {}

    for i,l in enumerate(list):
        s = l[0] + l[-1]
        if hist.has_key(s):
            hist[s] = hist[s] + 1
        else:
            hist[s] = 1

    print (hist)
    total = len(list)-hist['AA']-hist['TT']-hist['GG']-hist['CC']
    print (total)

    d = {}
    for i,l in enumerate(list):
        d[l] = i

    assert( len(d.keys()) == len(list) )
    #build_norepeats_codec("norepeatscodec.py","NoRepeatCodec")
