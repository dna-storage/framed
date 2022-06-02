import os
from lsf_submit import *
import json
import copy
import hashlib
import numpy as np


class NpEncoder(json.JSONEncoder): #this class should help with encoding numpy data types that will arise in the different ranges
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


def recursive_param_aggregation(l,out_list=None,previous=[]):
    if len(l)==0:
        #take previous names and values and create the parameter strings and directory strings
        param_string=""
        directory_string=""
        param_dict={}
        for p in previous:
            param_string+=" --{} {} ".format(p[0],p[1])
            directory_string+="{}={}:".format(p[0],p[1])
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
        recursive_param_aggregation(l[1:],out_list=return_list,previous=copy_previous)
    if out_list is None:
        return return_list


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Script to generate fault injection jobs")
    parser.add_argument('--params',type=str,required=True,action="store",help="Path to json file with parameters to perform fault injection with")
    parser.add_argument('--memory',default=8,type=int,action="store",help="Memory for each job")
    parser.add_argument('--cores',default=4,type=int,action="store",help="Cores for each job")
    parser.add_argument('--time',default=10,type=int,action="store",help="Time allowed for each job")
    parser.add_argument('--queue',default="tuck",type=str,action="store",help="Queue to use for jobs")
    parser.add_argument('--dump_dir',required=True,help="path to store results")
    args = parser.parse_args()

    job = LSFJob()
    job.queue=args.queue
    job.time=args.time
    job.cores=args.cores
    job.memory=args.memory
    assert os.path.exists(args.params)

    with open(args.params,"r") as json_fd:
        params_dict= json.load(json_fd)

    assert params_dict

    #need to at least have some basic params in the params file
    try:
        encoder_params = params_dict["encoder_params"]
        distribution_params=params_dict["distribution_params"]
        fault_model_params=params_dict["fault_params"]
        header_params = params_dict["header_params"]
        
        encoding_architecture = params_dict["arch"]
        fault_model = params_dict["fault_model"]
        file_path = params_dict["file"]
        distribution=params_dict["strand_distribution"]
       
    except:
        raise ValueError("Parameter file does not have right parameters")


    fault_injection_path = "/tuck_data/kvolkel/dnastorage/tools/fault_injection.py"
    
    
    base_file = os.path.basename(file_path).split(".")[0]
    #need to make up parts of directories where things are going to exist
    top_experiment_portion="{}:{}:{}".format(base_file,fault_model,distribution)

    run_path = os.path.join(args.dump_dir,top_experiment_portion)
    run_path = os.path.join(run_path,encoding_architecture)

    fault_model_param_list=recursive_param_aggregation([(_[0],_[1]) for _ in fault_model_params.items()])
    distribution_param_list=recursive_param_aggregation([(_[0],_[1]) for _ in distribution_params.items()])
    encoder_param_list=recursive_param_aggregation([(_[0],_[1]) for _ in encoder_params.items()])
    header_instance=recursive_param_aggregation([(_[0],_[1]) for _ in header_params.items()])[0]

    
    for distribution_instance in distribution_param_list:
        dist_run_path=os.path.join(run_path,distribution_instance[1])
        for fault_model_instance in fault_model_param_list:
            fault_run_path=os.path.join(dist_run_path,fault_model_instance[1])
            for encoder_instance in encoder_param_list:
                #Now build the lsf job and set it off
                md_hash = hashlib.md5()
                md_hash.update(json.dumps(encoder_instance[2],cls=NpEncoder).encode())
                final_run_path = os.path.join(fault_run_path,"encoder:{}".format(md_hash.hexdigest()))
                os.makedirs(final_run_path,exist_ok=True)                
                #round out the param string
                complete_param_string= distribution_instance[0]+fault_model_instance[0]+"--enc_params "+os.path.join(final_run_path,"encoder_params.json")
                complete_param_string+=" --out_dir {} ".format(final_run_path) + " --cores {} ".format(args.cores)
                complete_param_string+=" --header_params {} ".format(os.path.join(final_run_path,"header_params.json"))
            
                #round out the param string with stuff from the params dictionary
                for param in params_dict:
                    if not isinstance(params_dict[param],dict):
                        complete_param_string+=" --{} {} ".format(param,params_dict[param]) 
                #dump params 
                json.dump(distribution_instance[2],open(os.path.join(final_run_path,"dist_params.json"),"w+"),cls=NpEncoder)
                json.dump(fault_model_instance[2],open(os.path.join(final_run_path,"fault_params.json"),"w+"),cls=NpEncoder)
                json.dump(encoder_instance[2],open(os.path.join(final_run_path,"encoder_params.json"),"w+"),cls=NpEncoder)
                json.dump(header_instance[2],open(os.path.join(final_run_path,"header_params.json"),"w+"),cls=NpEncoder)
                
                print(final_run_path)
                
                command = "python "+ fault_injection_path + complete_param_string
                job.command = command
                job.run_path = final_run_path

                job.submit() #finally submit the experiment
