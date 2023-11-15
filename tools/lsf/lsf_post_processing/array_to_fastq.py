"""
Small tool to arrays of strands to fastq. Useful if using dumped strands in IGV.
"""

import os
import pickle
from Bio.SeqIO.QualityIO import FastqGeneralIterator

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="convert map to fastq")
    parser.add_argument('--original_fastq',default=None,required=True,help="<path> to original fastq file")
    parser.add_argument('--array_path', default=None,required=True,help="<path> to array file")
    parser.add_argument('--out_name',default=None,required=True,help="<path> to output")
    args = parser.parse_args()

    out_fasq=open(args.out_name,"w+")

    #need to load in maps
    array = pickle.load(open(args.map_path,"rb"))
    

    for s in array:
        qual = "("*len(s)
        title="@"
        if title in record_set:
            out_fasq.write("@{}\n{}\n+\n{}\n".format(title,seq,qual))
    out_fasq.close()
    
