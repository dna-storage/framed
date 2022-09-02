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

from nextera import *
from nupack.mfe import *
from nupack.complexes import *
from dnastorage.primer.primer_util import *

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def mutate_sequence(seq, d):
    used = {}
    for i in range(0,d):
        while True:
            r = random.randint(0,len(seq)-2)
            if not used.has_key(r):
                break
        used[r] = 1
        c = random.choice('ACGT')
        seq[r] = c

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


def design_rules_met(s,L,nupack):
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
        return False

    if not check_old_strands(s):
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
            if design_rules_met(l,L,use_nupack):
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
        s = id_generator(size=20,chars="ATGC")
        if not D.has_key(s):
            i = i+1
            if s[-1] != 'G':
                continue
            found = True
            if len(L) and s[0] == L[0][-1]:
                found = False
            else:
                found = design_rules_met(s,L,use_nupack)

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
    prefix = create_mfe_input(seqs)    
    complexes(complexes_args(prefix))
    c = read_mfe_output(prefix+".ocx-mfe")
    problems = [cc for cc in c if cc['deltaG'] < -10.0]
    return problems

def checkFold(seq):
    f = open(seq+".fasta","w")
    f.write(seq)
    f.close()
    return checkFoldFile(seq)    

def checkFoldFile(seq):
    fnull = open(os.devnull, "w")
    subprocess.call(['mfold','SEQ='+seq+'.fasta','NA=DNA'],stdout=fnull, stderr=fnull)
    f = open(seq+".out","r")
    s = f.read()
    f.close()
    rx = re.compile("dG =[ \t]*(-?\d+\.?\d+)")
    g = rx.search(s)
    if g!=None:
        d = float(g.groups()[0])
        if d < -3.0:
            return False
    return True

def runAll():
    global count
    genPrimers("",0,20)
    print "Found: ",count, " out of ", 4**20
    print "{}%".format(float(count)/(4**20)*100)


parser = argparse.ArgumentParser(description="Select a method for computing number of primers.")
parser.add_argument('--mc',help="Monte Carlo method", action='store_true')
parser.add_argument('--s',type=int,dest="simulations",action="store",default=100, help="Number of MonteCarlo simulations")
parser.add_argument('--distance',type=int,dest="distance",action="store",default=10, help="Hamming distance between primers")
parser.add_argument('--repetition',type=float,dest="repetition",action="store",default=0.80, help="Fraction of length-5 rolling window with repeating nucleotides")
parser.add_argument('--use-nupack',dest="use_nupack",action="store_true",help="Perform analysis using nupack.")


parser.add_argument('--primers',dest="primers",action="store",default=None, help="Previously selected primers.")

parser.add_argument('--timeout',type=int,dest="timeout",action="store",default=10000, help="If no primers are produced after timeout tries, give up.")

parser.add_argument('--o',dest="o",action="store",default=None, help="Output file.")


args = parser.parse_args()

if args.mc:
    (L,count,i) = montecarlo(args)

    if args.o != None:
        f = open(args.o,"w")
        for l in L:
            #print l,"  Tm=",mt.Tm_NN(l)
            f.write(l+"\n")
        f.close()
    else:
        for l in L:
            print l

    print "Simulation: ",count, " out of ", i
    print "Estimated primer population: ",4**20 * float(count)/i    
    if args.use_nupack:
        goodPairing = []
        used = {}
        for i in range(0,len(L)):
            for j in range(i+1,len(L)):
                if used.has_key(j):
                    continue
                pair = [ L[i], L[j] ]
                problems = checkComplexes(pair)
                if len(problems) == 0:
                    used[i] = 1
                    used[j] = 1
                    goodPairing.append(pair)
                    break
            if not used.has_key(i):
                print "Did not find pair for ",L[i],"."
                
        if len(goodPairing)>0:
            f = open(args.o+".pairs","w")
            i = 1
            for p in goodPairing:
                #f.write("% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n")
                f.write("% {}. Pair \n".format(i))
                i+=1
                f.write("{}\n{}\n".format(p[0],p[1]))
            f.close()
                
