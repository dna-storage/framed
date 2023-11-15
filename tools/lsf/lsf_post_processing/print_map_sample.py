import os
import pickle
import json

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="sequencing data analysis path for framed")
    parser.add_argument('--framed_seq_path',default=None,required=True,help="<path> to sequencing directory of framed output")
    parser.add_argument('--total_per_map',default=10,help="Number of strands to dump per map")
    args = parser.parse_args()

    assert os.path.exists(args.framed_seq_path) and os.path.isdir(args.framed_seq_path)

    sequencing_params_path = os.path.join(args.framed_seq_path,"sequencing_params.json")
    sequencing_map_path = os.path.join(args.framed_seq_path,"sequencing.map.pickle")
    
    assert os.path.exists(sequencing_params_path) and os.path.exists(sequencing_params_path)


    #load in indexed strands fastq IDs
    m = dict(pickle.load(open(sequencing_map_path,"rb")))
    m_sub_map_keys = m.keys()
    sub_map_counter={_:0 for _ in m_sub_map_keys if "payload" in _}
    
    original_sequencing_path=json.load(open(sequencing_params_path,"r"))["sequencing_data_path"]
    assert os.path.exists(original_sequencing_path)
    print("Getting records from {}".format(original_sequencing_path))
    
    for record in SeqIO.parse(original_sequencing_path,"fastq"):
        for k in m:
            if record.id in m[k] and sub_map_counter[k]<args.total_per_map:
                sub_map_counter[k]+=1
                print("Map Key {} \n {}".format(k,record.seq))


    
    

    
    
