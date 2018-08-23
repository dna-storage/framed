#!/usr/bin/python

import string
import time
import argparse

from dnastorage.primer.primer_util import *
from dnastorage.primer.design import *
from dnastorage.codec.base_conversion import convertQuarnary
from dnastorage.codec.commafreecodec import create_cfc_inv,cfc,cfc_inv
from random import randint,shuffle
import editdistance as ed


# perform substition on strand
def subst(ss):
    r = randint(0,20*8-1)    
    copy = [ _ for _ in ss ]
    syms = ['A', 'C', 'G', 'T']
    copy[r] = syms[ randint(0,3) ]
    ss = "".join(copy)
    return ss

# perform insertion on strand
def insert(ss):
    r = randint(0,20*8-1)    
    copy = [ _ for _ in ss ]
    syms = ['A', 'C', 'G', 'T']
    nt = syms[ randint(0,3) ]
    print "Insert {} at {}".format(nt,r)
    copy.insert(r, nt)
    ss = "".join(copy)
    return ss

# perform deletion on strand
def delete(ss):
    r = randint(0,20*8-1)    
    copy = [ _ for _ in ss ]
    print "Delete at {}".format(r)
    copy = copy[:r] + copy[r+1:]
    ss = "".join(copy)
    return ss

# naive decode
def decode(cfc_inv, s):
    if cfc_inv.has_key(s):
        return cfc_inv[s]
    return None

# strand decode
def cfc_decode(s,l):
    split = [ s[i:i+8] for i in range(0,l*8) ]
    for i in range(len(split)):
        split[i] = decode(cfc_inv,split[i])

    dec = [ None for _ in range(l) ]
    i = 0
    prior = True
    n = 0
    for n in range(l):
        if i >= l*8:
            break
        if split[i] != None:
            dec[n] = split[i]
            i += 8
            prior = True
        else:
            # maybe have an insertion, deletion, or substitution
            if split[i+8] != None:
                # guess substitution, leave entry as None
                i += 8
                prior = False
            elif split[i+8] == None: # somehow we are off track
                if prior: 
                    # this is the first one, so look forward from i for next non-None
                    k = i
                else:
                    k = i-4

                while k < i+4:
                    if split[k] != None:
                        i = k;
                        dec[n] = split[i]
                        prior = True
                        break
                    k += 1
                # whether found or not, move to next possible symbol
                if dec[n] == None:
                    prior = False
                i += 8                        
    return dec
                                

# get comma-free encoding and produce inverse map
create_cfc_inv()

# make a random strand
l = [ randint(0,255) for x in range(20) ]
# turn it into CFC
s = [ cfc[x] for x in l ]
ss = "".join(s)


# perform some random modifications
s = delete(ss)
s = insert(s)
#s = subst(s)
print s
print ss

# split up strand into sliding 8 nt chunks
split = [ s[i:i+8] for i in range(0,20*8) ]

# decode each piece
for i in range(len(split)):
    split[i] = decode(cfc_inv,split[i])

# print it out
print split

# just take the non-None values
split = [ x for x in split if x != None ]

# compare to original
print l
print split

# compare to smarter decode algorithm
r = cfc_decode(s,20)
print "{} {}".format(r, len(r))
