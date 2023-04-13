"""
Small tool to add strand length distributions of indexed and total sequenced strands to the sequencing.stats file.
Presents raw sizes, e.g. no stripping of primers
Useful for verifying how much the mapped distribution covers the sequencing file if needed.
"""

import os
import json
import shutil
import pickle
from Bio import SeqIO
import numpy as np

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="sequencing data analysis path for framed")
    parser.add_argument('--framed_seq_path',default=None,required=True,help="<path> to sequencing directory of framed output")
    args = parser.parse_args()

    assert os.path.exists(args.framed_seq_path) and os.path.isdir(args.framed_seq_path)

    sequencing_params_path = os.path.join(args.framed_seq_path,"sequencing_params.json")
    sequencing_map_path = os.path.join(args.framed_seq_path,"sequencing.map.pickle")
    sequencing_stats_path = os.path.join(args.framed_seq_path,"sequencing.pickle")
    
    assert os.path.exists(sequencing_stats_path) and os.path.exists(sequencing_params_path) and os.path.exists(sequencing_params_path)


    #load in indexed strands fastq IDs
    m = pickle.load(open(sequencing_map_path,"rb"))
    tmp_m = {}
    for sub_m in m:
        tmp_m={**tmp_m,**m[sub_m]}
    m = tmp_m    
    record_set = set(m.keys())
    m=None
    
    original_sequencing_path=json.load(open(sequencing_params_path,"r"))["sequencing_data_path"]
    assert os.path.exists(original_sequencing_path)
    print("Getting records from {}".format(original_sequencing_path))
    

    index_strand_lengths = []
    full_sequencing_file_lengths = []
    for record in SeqIO.parse(original_sequencing_path,"fastq"):
        if record.id in record_set:
            index_strand_lengths.append(len(record.seq))
        full_sequencing_file_lengths.append(len(record.seq))

    #add the distributions to the sequencing.pickle dictionary

    
    backup_path = os.path.join(args.framed_seq_path,"sequencing.stats.backup")
    
    if not os.path.exists(backup_path): shutil.copyfile(sequencing_stats_path,backup_path) #don't overwrite backup path until it is gone
    
    sequencing_stats = pickle.load(open(sequencing_stats_path,'rb'))

    sequencing_stats["indexed:dist"] = np.histogram(index_strand_lengths,bins="auto")
    sequencing_stats["total_seq:dist"] = np.histogram(full_sequencing_file_lengths,bins="auto")
    
    out_sequencing_stats_path = os.path.join(args.framed_seq_path,"sequencing.pickle")
    pickle.dump(sequencing_stats,open(out_sequencing_stats_path,'wb+'))
    
