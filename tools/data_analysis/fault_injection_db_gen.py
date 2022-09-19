'''
This script goes through a path where fault injection data has been dumped. Generates a dataframe and dumps it into a dataframe that can be 
accessed to generate data.
'''

import os
import pandas as pd
import numpy as np
import json
import pickle

class NpEncoder(json.JSONEncoder): #this class should help with encoding numpy data types that will arise in the different ranges
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def find_data_paths(path,out_list):
    pickle_path = False
    for p in os.listdir(path):
        total_path = os.path.join(path,p)
        if os.path.isdir(total_path):
            find_data_paths(total_path,out_list)
        if "fi.pickle" in total_path:
            pickle_path=True
    if pickle_path: out_list.append(path)


def load_json_dict(path,prefix=None):
    if not os.path.exists(path):
        return dict()
    if prefix:
        old_dict = json.load(open(path,'r'))
        new_dict = {"{}::{}".format(prefix,key):old_dict[key] for key in old_dict}
        return new_dict
    else:
        return json.load(open(path,'r'))

def load_pickle_dict(path):
    assert os.path.exists(path) and "Path does not exist when loading pickle file {}".format(path)
    data_dict =  pickle.load(open(path,'rb'))
    #expand arrays to make pandas processing easier
    additional_entries={}
    name_deletions=[]
    for item in data_dict:
        if isinstance(data_dict[item],np.ndarray):
            base_name = item
            for index, i in enumerate(data_dict[item]):
                name = "{}::{}".format(base_name,index)
                additional_entries[name]=i
            name_deletions.append(item)
    #for _ in name_deletions: del data_dict[_]
    data_dict.update(additional_entries)
    return data_dict

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Bring fault injection data together into a dataframe")                                                                                   
    parser.add_argument('--path',required=True,help="Path to generated data")

    args = parser.parse_args()

    out_list =[]
    assert(os.path.exists(args.path))

    find_data_paths(args.path,out_list)

    print("Directories with data {}".format(out_list))

    final_dicts=[]
    for p in out_list:
        try:
            all_dicts=[]
            all_dicts.append(load_json_dict(os.path.join(p,"encoder_params.json"),prefix="encoder"))
            all_dicts.append(load_json_dict(os.path.join(p,"distribution_params.json")))
            all_dicts.append(load_json_dict(os.path.join(p,"dna_process.json")))
            all_dicts.append(load_json_dict(os.path.join(p,"fault_params.json")))
            all_dicts.append(load_pickle_dict(os.path.join(p,"fi.pickle")))
            all_dicts.append(load_json_dict(os.path.join(p,"header_params.json"),prefix="header"))
            all_dicts.append(load_json_dict(os.path.join(p,"sim_params.json")))
            complete_dict = {}
            for x in all_dicts: complete_dict.update(x)
            final_dicts.append(complete_dict)
        except:
            continue
        
    out_df = pd.DataFrame(final_dicts,dtype=object).fillna('')
    out_df.to_csv(os.path.join(args.path,"dataframe.csv"),index=False)
    out_df.to_pickle(os.path.join(args.path,"dataframe.pickle"))
    
