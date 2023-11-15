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

def checkComplexes(seqs,deltaG):
    prefix = create_mfe_input(seqs)    
    complexes(complexes_args(prefix))
    c = read_mfe_output(prefix+".ocx-mfe")
    problems = [cc for cc in c if cc['deltaG'] < deltaG]
    return problems

def checkFold(seq):
    f = open(seq+".fasta","w")
    f.write(seq);
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


def verifyPairing(s,l,deltaG):
    # repeat s to search for homo-dimers
    c = checkComplexes([s,l],deltaG)
    if len(c) > 0:
        print "*****Heterodimer: {} vs {} -- {}".format(l,s, c[0]['pattern'])
        return False

    c = checkComplexes([l,reverse_complement(s)],deltaG)
    if len(c) > 0:
        print "*****False binding: {} vs {} -- {}".format(l,s, c[0]['pattern'])
        return False

    c = checkComplexes([s,reverse_complement(l)],deltaG)
    if len(c) > 0:
        print "*****False binding: {} vs {} -- {}".format(l,s, c[0]['pattern'])
        return False

def createPairing(primers,goodPairs,deltaG,verify):
    for g in goodPairs:
        used[g[0]] = 1
        used[g[1]] = 1
    used = {}
    dropped = []
    for i in range(0,len(primers)):
        if used.has_key(i):
            continue
        for j in range(i+1,len(primers)):
            if used.has_key(j):
                continue
            if verify and verifyPairing(primers[i],primers[j],deltaG):
                used[i] = 1
                used[j] = 1
                goodPairs.append(pair)
                break
        if not used.has_key(i):
            dropped.append(i)
    return goodPairs,dropped


def maxPairing(primers,goodPairs,deltaG,verify):
    for g in goodPairs:
        used[g[0]] = 1
        used[g[1]] = 1
    used = {}
    dropped = []
    for i in range(0,len(primers)):
        if used.has_key(i):
            continue
        for j in range(i+1,len(primers)):
            if used.has_key(j):
                continue
            if verify and verifyPairing(primers[i],primers[j],deltaG):
                used[i] = 1
                used[j] = 1
                goodPairs.append(pair)
                break
        if not used.has_key(i):
            dropped.append(i)
    return goodPairs,dropped


parser = argparse.ArgumentParser(description="Pair selected primers.")

parser.add_argument('primers',nargs=1,help="Previously selected primers.")

parser.add_argument('--o',dest="o",required=True,action="store",default=None, help="Output file.")
parser.add_argument('--deltaG',dest="deltaG",type=int,action="store",default=-10.0, help="deltaG")
parser.add_argument('--verify',dest="verify",type=int,action="store_true",default=False, help="Verify pairings are appropriate.")
parser.add_argument('--max',dest="max",type=int,action="store_true",default=False, help="Look for largest pairwise grouping.")

args = parser.parse_args()
primers = read_primers(args.primers[0])

if args.max:
    # find max possible pairings

else:
    goodPairing,dropped = createPairing(primers,[],args.deltaG,args.verify)
    if len(goodPairing)>0:
        f = open(args.o+".pairs","w")
        i = 1
        for p in goodPairing:
            #f.write("% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n")
            f.write("% {}. Pair \n".format(i))
            i+=1
            f.write("{}\n{}\n".format(p[0],p[1]))
        f.close()
                
