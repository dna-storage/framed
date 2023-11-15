'''
This script goes through a path where fault injection data has been dumped. Generates a dataframe and dumps it into a dataframe that can be 
accessed to generate data.
'''

import os
import pandas as pd
import numpy as np
import json
import pickle
import re

class NpEncoder(json.JSONEncoder): #this class should help with encoding numpy data types that will arise in the different ranges
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def find_data_paths(path,out_list,match_file):
    pickle_path = False
    for p in os.listdir(path):
        total_path = os.path.join(path,p)
        if os.path.isdir(total_path):
            find_data_paths(total_path,out_list,match_file)
        if match_file in total_path:
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
    data_dict.update(additional_entries)
    return data_dict


def load_fi_json(pickle_name,path,all_dicts):
    all_dicts.append(load_json_dict(os.path.join(path,"encoder_params.json"),prefix="encoder"))
    all_dicts.append(load_json_dict(os.path.join(path,"distribution_params.json")))
    all_dicts.append(load_json_dict(os.path.join(path,"dna_process.json")))
    all_dicts.append(load_json_dict(os.path.join(path,"fault_params.json")))
    all_dicts.append(load_json_dict(os.path.join(path,"header_params.json"),prefix="header"))
    all_dicts.append(load_json_dict(os.path.join(path,"sim_params.json")))
    all_dicts.append(load_pickle_dict(os.path.join(path,pickle_name)))

def load_seq_json(path,all_dicts):
    decoder_paths=[]
    header_paths=[]
    other_paths=[]
    for p in os.listdir(path):
        find_decoder_path=re.search("(decoder\_[0-9]+)(\.json)",p)
        find_header_path=re.search("(header\_[0-9]+)(\.json)",p)
        find_other_path=re.search("(other_params\_[0-9]+)(\.json)",p)
        if find_decoder_path:
            decoder_name = find_decoder_path.group(1)
            decoder_paths.append((p,decoder_name))
        if find_header_path:
            header_name=find_header_path.group(1)
            header_paths.append((p,header_name))
        if find_other_path:
            other_name=find_other_path.group(1)
            other_paths.append((p,other_name))
    for (dp,decoder_name) in decoder_paths: all_dicts.append(load_json_dict(os.path.join(path,dp),prefix="{}".format(decoder_name)))
    for (hp,header_name) in header_paths: all_dicts.append(load_json_dict(os.path.join(path,hp),prefix="{}".format(header_name)))
    for (op,other_name) in decoder_paths: all_dicts.append(load_json_dict(os.path.join(path,op),prefix="{}".format(other_name)))
    all_dicts.append(load_json_dict(os.path.join(path,"sequencing_params.json")))
    all_dicts.append(load_pickle_dict(os.path.join(path,"sequencing.pickle")))

if __name__=="__main__":
    import argparse
    import functools
    parser = argparse.ArgumentParser(description="Bring fault injection data together into a dataframe")                                                                                   
    parser.add_argument('--path',required=True,help="Path to generated data")
    parser.add_argument('--name',required=True,help="name of dataframe")
    parser.add_argument('--type',default="fi",choices=["fi","sequencing"],help="fi: data is from fault injection runs, sequencing: data is from sequencing runs")
    parser.add_argument('--ext',default="",help="extension to add to pickle files to look for")
    args = parser.parse_args()
    out_list =[]
    assert(os.path.exists(args.path))

    load_function=None
    if args.type=="fi":
        match_file = "fi.pickle"+args.ext
        load_function=functools.partial(load_fi_json,match_file)
    elif args.type=="sequencing":
        match_file="sequencing.pickle"+args.ext
        load_function=load_seq_json
    else:
        assert 0
    assert load_function and match_file

    find_data_paths(args.path,out_list,match_file)

#   print("Directories with data {}".format(out_list))
    print("Total directories with data {}".format(len(out_list)))
    final_dicts=[]
    for p in out_list:
        try:
            all_dicts=[]
            load_function(p,all_dicts)
            complete_dict = {}
            for x in all_dicts: complete_dict.update(x)
            final_dicts.append(complete_dict)
        except Exception as e:
            print(e)
            print("Error on path {}, data will not be loaded from here".format(p))
            continue
        
    out_df = pd.DataFrame(final_dicts,dtype=object).fillna('')
    out_df.to_csv(os.path.join(args.path,"{}.csv".format(args.name)),index=False)
    out_df.to_pickle(os.path.join(args.path,"{}.pickle".format(args.name)))
    
