"""
Parsing help for lsf json parameter scripts.
"""
import json
import numpy as np
import copy
import os
import hashlib

class NpEncoder(json.JSONEncoder): #this class should help with encoding numpy data types that will arise in the different ranges
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


def param_aggregation(l,out_list=None,previous=[]):
    if len(l)==0:
        #take previous names and values and create the parameter strings and directory strings
        param_string=""
        directory_string=""
        param_dict={}
        for p in previous:
            param_string+=" --{} {} ".format(p[0],p[1])
            if isinstance(p[1],str) and os.path.exists(p[1]): #TODO: make sure paths are handled sensibly
                md_hash=hashlib.md5()
                md_hash.update(p[1].encode("utf-8"))
                directory_string+="{}={}___".format(p[0],md_hash.hexdigest()[0:8])
            else: 
                directory_string+="{}={}___".format(p[0],p[1])
            param_dict[p[0]]=p[1]
        out_list.append((param_string,directory_string,param_dict))
        return
    
    if out_list is None:
        return_list=[]
    else:
        return_list=out_list
    entry_name= l[0][0]
    if type(l[0][1]) is list:
        if l[0][1][0] == "value_list":
            param_values=l[0][1][1:] #not a range
        else:
            param_values = np.arange(*(l[0][1]))
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
