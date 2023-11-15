"""
Small tool to add sequence computation time to pickled results of FrameD runs. Looks at log files
produced by each process, deduces the number of strands it is responsible for and determines how long it takes to
do the primary payload analysis on the strands.
"""

import os
import json
import shutil
import pickle
import numpy as np
import re
import datetime

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="sequencing data analysis path for framed")
    parser.add_argument('--framed_path',default=None,required=True,help="<path> to sequencing directory of framed output")
    args = parser.parse_args()


    #print("\n\nAnalyzing path {}\n\n".format(args.framed_path))
    assert os.path.exists(args.framed_path) and os.path.isdir(args.framed_path)

    stats_path = os.path.join(args.framed_path,"fi.pickle")
    strand_latency = 0
    total_strands=0
    for f in os.listdir(args.framed_path):
        file_search = re.search("^fi_info_[0-9]+.log$",f)
        if file_search is None: continue
        #print(f)
        num_strands=0
        start_datetime=None
        end_datetime=None
        scatter_flag=False
        with open(os.path.join(args.framed_path,f),'r') as log_file:
            for line in log_file:
                rank_strands = re.search("Rank ([0-9]+) has ([0-9]+) objects after complete scatter",line)
                inner_code_start_line = re.search("Header length",line)
                inner_code_stop_line = re.search("Rank [0-9]+ communicating [0-9]+ objects back to rank [0-9]+",line)  
                if "Scatter packets" in line: 
                    scatter_flag=True
                    continue  
                if rank_strands is not None and not scatter_flag:
                    num_strands = int(rank_strands.group(2))
                    assert num_strands>0
                    continue
                if inner_code_start_line:
                    start_time=re.search("([0-9]+):([0-9]+):([0-9]+),([0-9]+)",line).group(0).replace(",",".")+"000"
                    start_date = re.search("[0-9]+-[0-9]+-[0-9]+",line).group(0)
                    start_datetime=datetime.datetime.fromisoformat(start_date +" " + start_time)
                    continue
                if inner_code_stop_line:
                    if scatter_flag: 
                        scatter_flag=False
                        continue
                    end_time=re.search("([0-9]+):([0-9]+):([0-9]+),([0-9]+)",line).group(0).replace(",",".")+"000"
                    end_date = re.search("[0-9]+-[0-9]+-[0-9]+",line).group(0)
                    end_datetime=datetime.datetime.fromisoformat(end_date +" " + end_time)
                    continue
        if num_strands>0 and start_datetime is not None and end_datetime is not None:
            strand_latency+=(end_datetime-start_datetime).total_seconds()
            total_strands+=num_strands
    print("For Path {} latency is {}".format(args.framed_path,strand_latency/total_strands))


    sequencing_stats = pickle.load(open(stats_path,'rb'))
    sequencing_stats["average_strand_compute_time::optimistic"]=strand_latency/total_strands 
    sequencing_stats["average_strand_compute_time::pessimistic"]=strand_latency/sequencing_stats["payload::hedges::total_strands"]
    out_sequencing_stats_path = os.path.join(args.framed_path,"fi.pickle.runtime")
    pickle.dump(sequencing_stats,open(out_sequencing_stats_path,'wb+'))
    