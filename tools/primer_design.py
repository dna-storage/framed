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

    design_rules = build_standard_design_rules(L,use_nupack, distance=distance)

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

    t = time.time()
    i = 0
    while (count < args.n):
        if(args.timeout != None): 
          if(time.time() - t >= args.timeout):
             break 
        s = pg.get()
        i += 1
        if design_rules.check(s):
            L.append(s)
            count = count+1

    print design_rules
    return (L,count,i)

def ranged(args):
    use_nupack = args.use_nupack
    distance = args.distance
    length = args.primer_length
    L = [] # list of primers
   
    design_rules = build_standard_design_rules(L,use_nupack)
    
    pg=LinearPrimerGenerator(args.ranged)

    if(args.primers != None): 
        L.append(read_resumed_runs(args.primers))
        

    t = time.time()
    i=args.ranged[0]
    while (i<(args.ranged[0]+args.ranged[1]) and i<4**20): 
       if(args.timeout != None): 
          if(time.time() - t >= args.timeout):
             break 
       s=pg.get(i)
       i+=1
       found=True
       if len(L) and s[0] == L[0][-1]:
            found = False
       else:
            found = design_rules.check(s)
       if found:
            L.append(s)

    f = open("runInfo.txt", "a")
    f.write(str(design_rules))
    f.close()
    return L, i-args.ranged[0]

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

parser.add_argument('--ranged',type=int, nargs=2,dest="ranged",action="store",help="Checks certain range of primers. First arg is starting point and second is range from that point")


parser.add_argument('--o',dest="o",action="store",default=None, help="Output file.")

parser.add_argument('--notime',action="store_true", help="Get rid of time limitations for run")

args = parser.parse_args()

if args.notime: 
   args.timeout = None

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

elif args.ranged != None:
     (L,count)=ranged(args)

     if args.o != None:
        f = open(args.o,"w")
        if(count>=args.ranged[1]):
            f.write("finished Wooo\n")
        else:
            f.write(str(args.ranged[0]+count) +" "+ str(args.ranged[1]-count)+"\n")

        for l in L:
            #print l,"  Tm=",mt.Tm_NN(l)
            f.write(l+"\n")
        f.close()
        os.system("mv "+str(args.o)+" ./Primers")
     else:
        for l in L:
            print l
                   
