"""
Prints out a random subsample of record IDs that were mapped for a given sequencing run.
"""

import os
import pickle
import random

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="convert map to fastq")
    parser.add_argument('--map_path', default=None,required=True,help="<path> to mapping file")
    parser.add_argument('--number_maps',default=10,type=int,required=False,help="number of mapped sequencing IDs to output")
    parser.add_argument('--seed',default=0,type=int,required=False,help="seed to use for randomizing")
    args = parser.parse_args()
    
    #need to load in maps
    m = pickle.load(open(args.map_path,"rb"))
    tmp_m = {}
    for sub_m in m:
        tmp_m={**tmp_m,**m[sub_m]}
    m = tmp_m
    record_set = list(m.keys())
    rand = random.Random(args.seed)
    output_set = rand.sample(record_set,args.number_maps)
    for i in output_set: print("{}".format(i))
