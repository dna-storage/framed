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
from primer_util import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rank primers")
    parser.add_argument('--primers',dest="primers",action="store",default=None, help="Selected primers for library.")
    parser.add_argument('--distance',type=int,dest="distance",action="store",default=4, help="Minimum disallowable correlation distance")
    parser.add_argument('--o',dest="o",action="store",default=None, help="Output file.")
    args = parser.parse_args()

    primers = read_primers(args.primers)

    l = []
    for p in primers:
        l.append( [p,repetitionScore(p)] )
        
    l.sort(cmp=lambda x,y: cmp(y[1],x[1]))

    i = 1
    for ll in l:
        if ll[1] > .99:
            print "{}".format(ll[0])
        i = i+1

    if args.o != None:
        ofile = open(args.o, "w")
        ofile.write("Rank,Primer\n")
        for ll in l:
            ofile.write("{:5.2},{}\n".format(ll[1],ll[0]))
        ofile.close()
