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

#illumina_primers = [ 
#    'CAAGCAGAAGACGGCATACGAGAT',
#    'GTCTCGTGGGCTCGG',
#    'AATGATACGGCGACCACCGAGATCTACAC'
#    'TCGTCGGCAGCGTC'
#]
illumina_primers = [ 
    'GTCTCGTGGGCTCGG',
    'TCGTCGGCAGCGTC'
]

def nextera_chunks(size):
    l = {}
    for p in illumina_primers:
        for i in range(0,len(p)-size+1):
            l[ p[i:i+size] ] = 1
    return l.keys()

def get_nextera_primers():
    return illumina_primers[:]

def nextera_comparison(args):
    primers = read_primers(args.primers)
    errors = False
    distance = args.distance    
    #i5 = read_primers(args.i5)
    #i7 = read_primers(args.i7)    
    for p in primers:        
        print '-'*80 
        print "{}:".format(p)
        for i in illumina_primers:
            d = correlation_distance(p,i)
            if d > distance:
                errors = True
                print "Warning: {}-way correlated to {}".format(d,i)
                show_correlation(p,i)
                
    if not errors:
        print "All primers Ok."

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check against Nextera sequences")
    parser.add_argument('--primers',dest="primers",action="store",default=None, help="Selected primers for library.")
    #parser.add_argument('--i5',dest="i5",action="store",default=None, help="Name of the file with Illumina Nextera i5 sequences we compare against.")
    #parser.add_argument('--i7',dest="i7",action="store",default=None, help="Name of the file with Illumina Nextera i7 sequences we compare against.")
    parser.add_argument('--distance',type=int,dest="distance",action="store",default=4, help="Minimum disallowable correlation distance")

    args = parser.parse_args()
    nextera_comparison(args)

