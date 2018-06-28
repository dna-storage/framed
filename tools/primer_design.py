#!/usr/bin/python
import string
import time
import argparse

from dnastorage.primer.primer_util import *
from dnastorage.primer.design import *

def montecarlo(args):
    use_nupack = args.use_nupack
    distance = args.distance
    count = 0
    length = args.primer_length
    L = [] # list of primers

    design_rules = build_standard_design_rules(L,use_nupack)

    if args.fast:
        pg = LikelyPrimerGenerator(chars="AGCT",length=length)
    else:
        pg = UniquePrimerGenerator(chars="AGCT",length=length)

    if args.primers != None:
        f = open(args.primers,"r")
        tmp = f.readlines()
        f.close()
        kk = 0
        for l in tmp:
            l = l.strip()
            if len(l) == 0:
                continue
            if not args.skip:
                if design_rules.check(l):
                    L.append(l)
                    pg.append(l) # do not generate this strand as a possible one
                    if l[-1] != 'G':
                        print "Input primer ({}) does not end in G!".format(l)
                else:
                    #print "Removing {} from list.".format(l)
                    kk = kk + 1
            else:
                L.append(l)
                pg.append(l)
        print "Removed {} primers from list.".format(kk)

    t = time.time();
    i = 0
    while count < args.n and (time.time() - t <= args.timeout):
        s = pg.get()
        i += 1
        if design_rules.check(s):
            L.append(s)
            count = count+1

    print design_rules
    return (L,count,i)


parser = argparse.ArgumentParser(description="Select a method for computing number of primers.")

parser.add_argument('--distance',type=int,dest="distance",action="store",default=10, help="Hamming distance between primers")

parser.add_argument('--primer-length',dest="primer_length",type=int,action="store",default=20, help="Primer length.")

parser.add_argument('--use-nupack',dest="use_nupack",action="store_true",help="Perform analysis using nupack.")

parser.add_argument('--primers',dest="primers",action="store",default=None, help="Previously selected primers.")
parser.add_argument('--skip-check', dest="skip", help="Skip checking of old primers", action="store_true", default=False)


parser.add_argument('--mc',help="Monte Carlo method", action='store_true')

parser.add_argument('--n',type=int,dest="n",action="store",default=100, help="Number of new primers sought")

parser.add_argument('--fast',dest="fast",action="store_true",help="Use a faster search that doesn't generate primers randomly.")

parser.add_argument('--timeout',type=int,dest="timeout",action="store",default=60, help="If no primers are produced after timeout tries, give up.")

parser.add_argument('--o',dest="o",action="store",default=None, help="Output file.")

args = parser.parse_args()

if args.mc:
    (L,count,i) = montecarlo(args)

    if args.o != None:
        f = open(args.o,"w")
        for l in L:
            # create verbose command line option?
            #print l,"  Tm=",mt.Tm_NN(l)
            f.write(l+"\n")
        f.close()
    else:
        for l in L:
            print l

    print "Simulation: ",count, " out of ", i
    print "Estimated primer population: ",4**20 * float(count)/i                    
