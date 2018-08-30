#!/usr/bin/python
import string
import time
import argparse

from dnastorage.primer.primer_util import *
from dnastorage.primer.design import *
from dnastorage.codec.base_conversion import convertQuarnary
from dnastorage.codec.commafreecodec import cfc
from random import randint,shuffle
import editdistance as ed

"""
If s can be rotated to create a value equal to s, then it has a self cycle
"""
def self_cycle_check(s):
    t = s+s
    # create all rotations of s
    l = [ t[i:i+6] for i in range(1,6) ]
  
    # if s is present in l, there's a cycle
    if s in l:
       # print s
        #print l
        return True
    return False

"""
Make sure that b does not appear in a cycle of a
"""
def same_cycle(a,b):
    t = a+a
    l = [ t[i:i+6] for i in range(1,6) ]
    return b in l

def print_characteristics(all):
    base = { 'A' : [],
             'C' : [],
             'G' : [],
             'T' : []  } 
    gc = {}
    for s in all:
        base[s[0]].append(s)
        c = countGC(s)
        if not gc.has_key(c):
            gc[c] = [ ]
        #print "{} {}".format(d,c)
        gc[c].append(s)

    for key,item in gc.items():
        print "{} - {}".format(key,len(item))

    for key,item in base.items():
        print "{}:{} {}".format(key,len(item),item[0:4])


def comma_free_check(all):
    # pick random starting point
    while True:
        i = randint(0,len(all))
        if self_cycle_check(all[i]):
            continue        
        start = [ all[i] ]
        break

    #randomize list order
    shuffle(all)

    for i,s in enumerate(all):
        rs = reverse_complement(s)
        
        if (s in start) or (rs in start):
            continue

        if self_cycle_check(s):
            continue

        found = True
        for a in start:
            if same_cycle(s,a):
                found = False
                break
            #if hamming_distance(a,s) == 1:
            #    found = False
            #    break
            t = a + a
            if t.find(s) >= 0 or t.find(rs) >= 0:
                found = False
                break

        if not found:
            continue

        for k,a in enumerate(start):
            for l in range(k+1,len(start)):
                b = start[l]
                if a == b:
                    continue
                t = a + b + a
                if t.find(s) >= 0: #  or t.find(rs) >= 0:
                    found = False
                    break

                t = a + s + a
                if t.find(b) >= 0:
                    found = False
                    break

                t = b + s + b
                if t.find(a) >= 0:
                    found = False
                    break

            if not found:
                break
            
        if not found:
            if i%100 == 0:
                print "{}% tested".format(float(i)/len(all)*100.0)
            continue
        else:
            print "{} {} - {}".format(s,len(start), start[-5:])
            start.append(s)

    return start


def countGC(s):
    l = [ _ for _ in s if _ == 'G' or _ == 'C' ]
    return len(l)


# design a comma free codec
def comma_free_codec_design():

    # can't match these sequences
    index1 = [ 'CAAGCAGAAGACGGCATACGAGAT',  'GTCTCGTGGGCTCGG' ]
    index2 = [ 'AATGATACGGCGACCACCGAGATCTACAC', 'TCGTCGGCAGCGTC' ]
    i7 = [ 'TCGCCTTA', 'CTAGTACG', 'TTCTGCCT', 'GCTCAGGA', 'AGGAGTCC', 'CATGCCTA', 'GTAGAGAG', 'CCTCTCTG', 'AGCGTAGC', 'CAGCCTCG', 'TGCCTCTT', 'TCCTCTAC' ]
    i5 = [ 'TAGATCGC', 'CTCTCTAT', 'TATCCTCT', 'AGAGTAGA', 'GTAAGGAG', 'ACTGCATA', 'AAGGAGTA', 'CTAAGCCT', 'GCGTAAGA' ]
    
    # make all sequences
    index1ex = [ index1[0]+x+index1[1] for x in i7 ]
    index2ex = [ index2[0]+x+index2[1] for x in i5 ]

    # put them in here in 8 nit chunks
    allbytes = []

    for a in index1ex:
        allbytes += [a[i:i+6] for i in range(len(a)-5)]
        
    for a in index2ex:
        allbytes += [a[i:i+6] for i in range(len(a)-5)]

    allbytes += [ reverse_complement(a) for a in allbytes ]

    all = []

    # rule out certain harmful sequences
    for i in range(2**12):
        d = convertQuarnary(i,6)
        if hasRepeat(d):
            continue
        if hasShortDimerRun(d):
            continue
        #if d[0] == 'G':
        #    continue
        #if d[-1] == 'C':
        #    continue
        if not ( (d[0] == 'C' and d[-1]=='A') or ((d[0] == 'C' and d[-1]=='T'))):
            continue
        #if d[0] == d[-1]:
        #    continue
        if d in allbytes:
            continue
        if reverse_complement(d)==d:
            continue
        if self_cycle_check(d):
            continue

        all.append(d)

    
    new_cfc = comma_free_check(all)
    # sort list to prefer balanced GC content
    new_cfc.sort(cmp=lambda x,y: (countGC(x)-4)**2 < (countGC(y)-4)**2 )
    # take first 256, this could be smarter
    return new_cfc[0:256]

# design codes
cfc = comma_free_codec_design()
# print GC content and other stats
print_characteristics(cfc)

# dump list that can be copied into commafreecodec.py
print "cfc = {}".format(cfc)

# exit program, comment out for ED analysis
sys.exit(0)

# analyze edit distance of selected codes
ED = {}
for a in range(len(cfc)):
    for b in range(a,len(cfc)):
        if cfc[a] != cfc[b]:
            d = ed.eval(cfc[a],cfc[b])
            ED.get(d,0)
            if not ED.has_key(d):
                ED[d] = 0
            ED[d] += 1
            #if d == 1:
            #    print "{} - {}".format(cfc[a],cfc[b])

print ED


