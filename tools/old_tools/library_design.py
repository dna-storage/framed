#!/usr/bin/python
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import random
import string
import subprocess
import re
import os
import argparse

from dnastorage.primer.nextera import *
from nupack.mfe import *
from nupack.complexes import *
from dnastorage.primer.primer_util import *

bases = ['A', 'C', 'G', 'T']

nextBases = { 'A' : [ 'C', 'G', 'T' ],
              'C' : [ 'A', 'G', 'T' ],
              'G' : [ 'A', 'C', 'T' ],
              'T' : [ 'A', 'C', 'G' ]
}

count = 0

def genPrimers(seq,n,l): 
    global count

    if len(seq)>=4 and seq[-1] == seq[-2] and seq[-2] == seq[-3] and seq[-3] == seq[-4]:
        B = nextBases[ seq[-1] ]
    else:
        B = bases[:]
        
    for b in B:
        s = seq + b
        if len(s) == l:
            if hasSingleRun(seq) or hasDimerRun(seq):
                continue
            r = Seq.reverse_complement(Seq(seq))
            if check(seq) and check(r):
                count = count + 1
                print seq
        else:
            genPrimers(s,n+1,l)


def design_rules_met(s,L,nupack,nextera_binding):
    #print s
    if repetitionScore(s) < 0.99:
        #print "reptition score", repetitionScore(s)
        return False
    if hasSingleRun(s) or hasDimerRun(s): # or hasDimer(seq,5):
        #print "runs"
        return False

    # ending of primer should not have too much GC content
    if not checkGC(s[-5:],(0,60)):
        #print s[-5:],"GC short content",GC(s[-5:])
        return False

    for l in L:     
        if correlation_distance(s,l)>4:
            #print "correlation"
            return False
        if hamming_distance(s,l) < 10:
            #print "hamming"
            return False    
        if nupack:
            # repeat s to search for homo-dimers
            c = checkComplexes([s,l])
            if len(c) > 0:
                print "*****Heterodimer: {} vs {} -- {}".format(l,s, c[0]['pattern'])
                return False

            c = checkComplexes([l,reverse_complement(s)])
            if len(c) > 0:
                print "*****False binding: {} vs {} -- {}".format(l,s, c[0]['pattern'])
                return False

            c = checkComplexes([s,reverse_complement(l)])
            if len(c) > 0:
                print "*****False binding: {} vs {} -- {}".format(l,s, c[0]['pattern'])
                return False

        #if reverse_correlation_distance(s,l)>3:
        #    show_correlation(s,reverse(l))
        #    return False

    if not checkTm(s,(50.0,60.0)):
        #print "TM"
        return False

    if not checkGC(s,(40.0,60.0)):
        #print "GC full"
        return False

    if not nextera_strand_comparison(s,3):
        #print "nextera"
        return False

    if nupack and nextera_binding:
        nextera_primers = get_nextera_primers()
        for l in nextera_primers:
            # repeat s to search for homo-dimers
            c = checkComplexes([s,l])
            if len(c) > 0:
                print "*****Nextera Heterodimer: {} vs {} -- {}".format(l,s, c[0]['pattern'])
                return False

            c = checkComplexes([l,reverse_complement(s)])
            if len(c) > 0:
                print "*****Nextera False binding: {} vs {} -- {}".format(l,s, c[0]['pattern'])
                return False

            c = checkComplexes([s,reverse_complement(l)])
            if len(c) > 0:
                print "*****Nextera Complement False binding: {} vs {} -- {}".format(l,s, c[0]['pattern'])
                return False


    if not check_old_strands(s):
        #print "check old"
        return False

    if nupack:
        # repeat s to search for homo-dimers
        c = checkComplexes([s,s])
        if len(c) > 0:
            print "*****Homodimer: {} {}".format(s, c[0]['pattern'])
            return False

    return True

def montecarlo(args):
    N = args.simulations
    use_nupack = args.use_nupack
    distance = args.distance
    count = 0
    timeout = args.timeout

    L = [] # list of primers
    D = {} # dictionary for tracking which primers we have

    primers = []
    if args.primers != None:
        f = open(args.primers,"r")
        tmp = f.readlines()
        f.close()
        kk = 0
        for l in tmp:
            l = l.strip()
            if len(l) == 0:
                continue
            if design_rules_met(l,L,use_nupack,args.check_nextera_binding):
                L.append(l)
                if l[-1] != 'G':
                    print "Input primer ({}) does not end in G!".format(l)
                D[l] = 1
            else:
                #print "Removing {} from list.".format(l)
                kk = kk + 1
        print "Removed {} primers from list.".format(kk)

    i = 0 # number of attempts
    while count < N and i < timeout:
        if args.fast:
            s = fast_generator(size=20,chars="ATGC")
        else:
            s = id_generator(size=20,chars="ATGC")
        if not D.has_key(s):
            i = i+1
            if s[-1] != 'G':
                continue
            found = True
            if len(L) and s[0] == L[0][-1]:
                found = False
            else:
                found = design_rules_met(s,L,use_nupack,args.check_nextera_binding)

            if found:
                L.append(s)
                D[s] = 1
                count = count+1

    return (L,count,i)


def check(seq):
    s = seq
    gc = GC(s)
    if gc >= 40 and gc <= 60:            
        t = mt.Tm_NN(s)
        if t >= 50 and t <= 60:
            return True        
    return False


complement = { 'A' : 'T',
               'C' : 'G',
               'G' : 'C',
               'T' : 'A' }

def checkComplexes(seqs):
    prefix = create_mfe_input(seqs,2)    
    complexes(complexes_args(prefix))
    #c = read_mfe_output(prefix+".ocx-mfe")
    #print prefix+".ocx-mfe"
    #problems = [cc for cc in c if cc['deltaG'] < threshold]
    c = read_ocx_output(prefix+".ocx")
    return c

def printComplex(l,s,d):
    print "*****Nextera Heterodimer: {} vs {} -- {}".format(l,s, c[0]['pattern'])


parser = argparse.ArgumentParser(description="Analyze hybridization potential of given primers w.r.t. to the given library.")

parser.add_argument('--use-nupack',dest="use_nupack",action="store_true",help="Perform analysis using nupack.")

parser.add_argument('--primers',dest="primers",action="store",default=None, help="Previously selected primers.")
parser.add_argument('--library',dest="library",action="store",default=None, help="comma separated string of files containing the DNA library")


args = parser.parse_args()

if args.primers != None and args.library != None:
    print "Analyzing primers and strands."
    Primers = read_primers(args.primers)
    Library = []

    library_files = args.library.split(',')
    for f in library_files:
        # read_primers just grabs the strands without concern for what they really are
        # I use it here for simplicity
        Library += read_primers(f)

    for l in Library:
        print l
        five_prime = l[0:20]
        three_prime = reverse_complement(l[-20:])         

        pSet = [x for x in Primers if not x in [five_prime,three_prime]]
        #pSet = Primers
        #for p in Primers:
        #    print "Next primer: ",p
        #    if p == five_prime or p == three_prime:
        #        continue
            
        

        strands = [l]+pSet
        c = checkComplexes(strands)
        
        deltaGbase = c['complexes'][1]['deltaG']
        bad = [ ll for i,ll in c['complexes'].items() if ll['order'] == 220 and ll['deltaG']-deltaGbase < -12 ]
        if len(bad)>0:        
            print ">>Strand conflicts with one or more primers."
            for b in bad:
                print "   {}-{}: {}".format(b['deltaG'],deltaGbase,b['sequences'][1])


        strands = [reverse_complement(l)]+pSet
        c = checkComplexes(strands)
        deltaGbase = c['complexes'][1]['deltaG']
        bad = [ ll for i,ll in c['complexes'].items() if ll['order'] == 220 and ll['deltaG']-deltaGbase < -12 ]
        if len(bad)>0:        
            print ">>Anti-sense strand conflicts with one or more primers."
            for b in bad:
                print "   {}-{}: {}".format(b['deltaG'],deltaGBase,b['sequences'][1])