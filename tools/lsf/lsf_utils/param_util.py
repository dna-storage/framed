"""
Parsing help for lsf json parameter scripts.
"""
import json
import numpy as np
import copy
import os
import hashlib
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

def regex_path_search(root_path,regexp,is_file):
    #search for paths from the root that meet the regexp
    reg_exp_prog = re.compile(regexp)
    paths=[]
    out_paths=[]
    for (root,dirs,files) in os.walk(root_path,topdown=True): 
        if is_file:
            for f in files:
                path=os.path.join(root,f)
                match = reg_exp_prog.search(path)
                if match is None: continue
                out_paths.append(path)
        else:
            for d in dirs:
                path=os.path.join(root,d)
                match = reg_exp_prog.search(path)
                if match is None: continue
                out_paths.append(path)
    return out_paths
    
def param_aggregation(l,out_list=None,previous=[]):
    if len(l)==0:
        #take previous names and values and create the parameter strings and directory strings
        param_string=""
        directory_string=""
        param_dict={}
        for p in previous:
            name,value = p[0],p[1]
            param_string+=" --{} {} ".format(name,value)
            if isinstance(value,str) and os.path.exists(os.path.expandvars(value)): 
                value = os.path.expandvars(value) #This allows for env variables in jsons
                md_hash=hashlib.md5()
                md_hash.update(value.encode("utf-8"))
                directory_string+="{}={}___".format(name,md_hash.hexdigest()[0:8])
            else: 
                directory_string+="{}={}___".format(name,value)
            param_dict[name]=value
        out_list.append((param_string,directory_string,param_dict))
        return
    
    if out_list is None:
        return_list=[]
    else:
        return_list=out_list
    entry_name= l[0][0]
    params = l[0][1]
    if type(params) is list:
        value_type=l[0][1][0]
        if value_type == "value_list":
            param_values=params[1:] 
        elif value_type=="dir_name":
            param_values = regex_path_search(os.path.expandvars(params[1]),params[2],False)
        elif value_type=="file_name":
            param_values = regex_path_search(os.path.expandvars(params[1]),params[2],True)
        elif value_type=="range":
            param_values = np.arange(*(params[1:]))
        else:
            raise ValueError("value type {} not supported or specified".format(value_type))
    else:
        param_values=[l[0][1]]
    for v in param_values:
        copy_previous=copy.copy(previous)
        current_value=v
        current_name=entry_name
        copy_previous.append((current_name,current_value))
        param_aggregation(l[1:],out_list=return_list,previous=copy_previous)
    if out_list is None:
        return return_list
