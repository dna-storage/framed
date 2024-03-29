"""
User level script to generate fault injection jobs through LSF.
"""

import os
from lsf_utils.submit import *
from lsf_utils.param_util import *
import hashlib

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Script to generate fault injection jobs")
    parser.add_argument('--params',type=str,required=True,action="store",help="Path to json file with parameters to perform fault injection with")
    parser.add_argument('--memory',default=None,action="store",help="Memory for each job")
    parser.add_argument('--cores',default=4,type=int,action="store",help="Cores for each job, this is the active number of working processes during fi execution")
    parser.add_argument('--time',default=10,type=int,action="store",help="Time allowed for each job,in hours")
    parser.add_argument('--queue',default="tuck",type=str,action="store",help="Queue to use for jobs")
    parser.add_argument('--dump_dir',required=True,help="path to store results")
    parser.add_argument('--core_depth',default=1,type=int,action="store",help="the number of additional cores to request per core. Example use case is if you want multiple mpi cores per mpi task")
    parser.add_argument('--job_name',default="dnastorage_fi",action="store",help="name for jobs that will be spawned")
    parser.add_argument('--avoid_hosts',default=None,nargs='+',help="hosts to avoid when running jobs")
    parser.add_argument('--experiment_prefix',default="",type=str,action="store",help="custom prefix to combine with top level directory name")
    parser.add_argument('--no',action='store_true', help="don't run bsub command, just test everything.")
    parser.add_argument('--conda_env_path',action="store",default=None,help="conda env to load")
    parser.add_argument('--modules',default=list(), nargs='+',help="modules to load")
    parser.add_argument('--submission', default = "LSF",choices = ["LSF","shell"])
    args = parser.parse_args()

    #Job parameter setup 
    if args.submission=="LSF":
        job = LSFJob()
        job.queue=args.queue
        job.time=args.time
        job.memory=args.memory
        job.one_host=False
        job.job_name = args.job_name
        if args.avoid_hosts!=None:
            job.avoid_hosts=args.avoid_hosts
    elif args.submission=="shell":
        #shell execution
        job = TcshJob()
    else:
        print("Did not specify correct execution backend")
        exit(1)
    job.cores=args.cores*args.core_depth+1 #plus one just to make sure we have enough cores
    job.load_modules=args.modules
    job.using_conda_env=args.conda_env_path
    job.using_ncsu_mpi = True #want to use ncsu's MPI enviroment

    #loading parameters for jobs
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
        fi_env_params = params_dict["fi_env_params"]
        encoding_architecture = params_dict["arch"]
        file_path = os.path.expandvars(params_dict["file"])
    except Exception as e:
        raise ValueError("Parameter file does not have right parameters: {}".format(e))

    #KV: Fixed hard-coded fault injection path, note you need to source the dnastorage environment in the base project directory
    if "DNASTORAGE_TOOLS" not in os.environ:
        print("Env Variable DNASTORAGE_TOOLS not found, make sure to source dnastorage.env in base project directory")

    fault_injection_path = os.environ["DNASTORAGE_TOOLS"] 
    fault_injection_path = os.path.join(fault_injection_path,"fault_injection.py")

    assert os.path.exists(fault_injection_path)

    base_file = os.path.basename(file_path).split(".")[0]
    #need to make up parts of directories where things are going to exist
    top_experiment_portion="{}___{}___{}___{}".format(args.experiment_prefix,base_file,fi_env_params["fault_model"],fi_env_params["distribution"])

    run_path = os.path.join(args.dump_dir,top_experiment_portion)
    run_path = os.path.join(run_path,encoding_architecture)

    fault_model_param_list=param_aggregation([(_[0],_[1]) for _ in fault_model_params.items()])
    distribution_param_list=param_aggregation([(_[0],_[1]) for _ in distribution_params.items()])
    encoder_param_list=param_aggregation([(_[0],_[1]) for _ in encoder_params.items()])
    header_instance=param_aggregation([(_[0],_[1]) for _ in header_params.items()])[0]
    env_instance=param_aggregation([(_[0],_[1]) for _ in fi_env_params.items()])[0]

    if "dna_processing" in params_dict:
        dna_proc_dict = params_dict["dna_processing"]
    else:
        dna_proc_dict = None
    
    for distribution_instance in distribution_param_list:
        dist_run_path=os.path.join(run_path,distribution_instance[1])
        for fault_model_instance in fault_model_param_list:
            fault_run_path=os.path.join(dist_run_path,fault_model_instance[1])
            for encoder_instance in encoder_param_list:
                #Now build the lsf job and set it off
                md_hash = hashlib.md5()
                md_hash.update(json.dumps(encoder_instance[2],cls=NpEncoder).encode())
                final_run_path = os.path.join(fault_run_path,"encoder___{}".format(md_hash.hexdigest()[0:8]))
                os.makedirs(final_run_path,exist_ok=True)                
                #Create the paramater string
                complete_param_string= " --enc_params \"{}\" ".format(os.path.join(final_run_path,"encoder_params.json"))
                complete_param_string+=" --header_params \"{}\" ".format(os.path.join(final_run_path,"header_params.json"))
                complete_param_string+=" --fi_env_params \"{}\" ".format(os.path.join(final_run_path,"fi_env_params.json"))
                complete_param_string+=" --distribution_params \"{}\" ".format(os.path.join(final_run_path,"distribution_params.json"))
                complete_param_string+=" --fault_params \"{}\" ".format(os.path.join(final_run_path,"fault_params.json"))
                complete_param_string+=" --out_dir \"{}\" ".format(final_run_path) + " --cores {} ".format(args.cores)
                if dna_proc_dict is not None: complete_param_string+=" --dna_process \"{}\"".format(os.path.join(final_run_path,"dna_process.json"))
                
                sim_param_dict = {} #parameters related to the overall simulater
                #round out the param string with stuff from the params dictionary
                for param in params_dict:
                    if not isinstance(params_dict[param],dict) and not isinstance(params_dict[param],list):
                        complete_param_string+=" --{} {} ".format(param,params_dict[param])
                        sim_param_dict[param]=params_dict[param]
                    elif isinstance(params_dict[param],list):
                        #expand tuples (needed for bacodes)
                        param_string = "".join([" {} ".format(x) for x in params_dict[param]])
                        complete_param_string +=" --{} {} ".format(param,param_string)
                        sim_param_dict[param]=params_dict[param]
                        
                #dump params 
                json.dump(sim_param_dict,open(os.path.join(final_run_path,"sim_params.json"),"w+"),cls=NpEncoder)
                json.dump(distribution_instance[2],open(os.path.join(final_run_path,"distribution_params.json"),"w+"),cls=NpEncoder)
                json.dump(fault_model_instance[2],open(os.path.join(final_run_path,"fault_params.json"),"w+"),cls=NpEncoder)
                json.dump(encoder_instance[2],open(os.path.join(final_run_path,"encoder_params.json"),"w+"),cls=NpEncoder)
                json.dump(header_instance[2],open(os.path.join(final_run_path,"header_params.json"),"w+"),cls=NpEncoder)
                json.dump(env_instance[2],open(os.path.join(final_run_path,"fi_env_params.json"),"w+"),cls=NpEncoder)
                if dna_proc_dict is not None: json.dump(dna_proc_dict,open(os.path.join(final_run_path,"dna_process.json"),"w+"),cls=NpEncoder)

                print(final_run_path)
                
                command = "python "+ fault_injection_path + complete_param_string
                job.command = command
                job.run_path = final_run_path

                job.submit(args.no) #finally submit the experiment
