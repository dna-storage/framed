"""
Small tool to convert maps to fastq files
"""

import os
import pickle
from Bio.SeqIO.QualityIO import FastqGeneralIterator

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="convert map to fastq")
    parser.add_argument('--original_fastq',default=None,required=True,help="<path> to original fastq file")
    parser.add_argument('--map_path', default=None,required=True,help="<path> to mapping file")
    parser.add_argument('--out_name',default=None,required=True,help="<path> to output")
    args = parser.parse_args()

    out_fasq=open(args.out_name,"w+")

    #need to load in maps
    m = pickle.load(open(args.map_path,"rb"))
    tmp_m = {}
    for sub_m in m:
        tmp_m={**tmp_m,m[sub_m]}
    m = tmp_m
    
    record_set = set(m.keys())
    m=None

    for title, seq, qual in FastqGeneralIterator(open(args.original_fastq,"r")):
        if title in record_set:
            out_fasq.write("@{}\n{}\n+\n{}\n".format(title,seq,qual))
    out_fasq.close()
    
