"""
User level script to generate sequencing jobs through LSF.

"""

import os
import errno
from lsf_utils.lsf_submit import *
from lsf_utils.param_util import *
import hashlib
import itertools

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Script to generate fault injection jobs")
    parser.add_argument('--sequencing_path',required=True,action="store",help="Top path to sequencing data to analyze, can be a directory or a single file")
    parser.add_argument('--params',type=str,required=True,action="store",help="Path to json file with parameters to perform sequencing analysis with")
    parser.add_argument('--memory',default=8,type=int,action="store",help="Memory for each job")
    parser.add_argument('--cores',default=4,type=int,action="store",help="Cores for each job, this is the active number of working processes during fi execution")
    parser.add_argument('--time',default=10,type=int,action="store",help="Time allowed for each job,in hours")
    parser.add_argument('--queue',default="tuck",type=str,action="store",help="Queue to use for jobs")
    parser.add_argument('--dump_dir',required=True,help="path to store results")
    parser.add_argument('--no',action='store_true', help="don't run bsub command, just test everything.")
    args = parser.parse_args()

    job = LSFJob(modules=['PrgEnv-intel'])
    job.queue=args.queue
    job.time=args.time
    job.cores=args.cores 
    job.memory=args.memory
    job.one_host=False
    job.using_ncsu_mpi = True #want to use ncsu's MPI enviroment

    job.load_module = ['conda']
    assert os.path.exists(args.params)

    try:
        with open(args.params,"r") as json_fd: params_dict= json.load(json_fd)
    except Exception as e:
        print("Issue opening param file: {}".format(e))
        exit(1)
        
    #need to at least have some basic params in the params file
    try:
        #Get mandatory parameters  
        decoders = params_dict["file_list"]
    except Exception as e:
        raise ValueError("Parameter file does not have right parameters: {}".format(e))
    
    #KV: Fixed hard-coded fault injection path, note you need to source the dnastorage environment in the base project directory
    if "DNASTORAGE_TOOLS" not in os.environ:
        print("Env Variable DNASTORAGE_TOOLS not found, make sure to source dnastorage.env in base project directory")
    
    
    sequencing_tool_path = os.environ["DNASTORAGE_TOOLS"] 
    sequencing_tool_path = os.path.join(sequencing_tool_path,"sequencing_analysis.py")

    assert os.path.exists(sequencing_tool_path)

    #need to make up parts of directories where things are going to exist
    top_experiment_portion="Sequencing___{}".format(os.path.basename(args.sequencing_path.rstrip("/")))
    run_path = os.path.join(args.dump_dir,top_experiment_portion)

    #Need to go through each file in the file list and create parameter sets for each one
    decoder_param_lists=[]
    for d in decoders:
        assert isinstance(d,dict)
        decoder_param_lists.append(param_aggregation([(_[0],_[1]) for _ in d["encoder_params"].items()]))
        #TODO: may also need to make it so that header parameters can be easily varried
        

    decoder_prod = itertools.product(*decoder_param_lists)
    decoder_combinations=[_ for _ in decoder_prod]

    #get sequencing files
    sequencing_paths = []
    assert os.path.exists(args.sequencing_path)
    if os.path.isdir(args.sequencing_path):
        for p in os.listdir(args.sequencing_path):
            sequencing_paths.append(os.path.join(args.sequencing_path,p))
    elif os.path.isfile(args.sequencing_path):
        sequencing_paths.append(args.sequencing_path) 
    else:
        raise ValueError("Issue with sequencing data paths")

    #set up commands for sequencing analysis
    for path in sequencing_paths:
        sequencing_run_path = os.path.join(run_path,os.path.basename(path))
        for decoder_combination in decoder_combinations:
            complete_param_string=" --sequencing_data_path {} ".format(path)
            #build a hash to differentiate sequencing analysis runs
            md_hash = hashlib.md5()
            decoder_dicts = [_[2] for _ in decoder_combination]
            assert len(decoder_dicts)==len(decoders)
            for decoder in decoder_dicts: md_hash.update(json.dumps(decoder,cls=NpEncoder).encode())
            final_run_path = os.path.join(sequencing_run_path,"decoder___{}".format(md_hash.hexdigest()))
            os.makedirs(final_run_path,exist_ok=True)

            complete_param_string+=" --dna_file_params {} ".format(os.path.join(final_run_path,"file_params.json"))
            complete_param_string+=" --out_dir {} ".format(final_run_path)

            analysis_dict={}
            analysis_dict["file_list"]=[]
            for decoder_index, decoder in enumerate(decoder_dicts): #dump some dictionaries for documentation purposes
                json.dump(decoder,open(os.path.join(final_run_path,"decoder_{}.json".format(decoder_index)),"w+"),cls=NpEncoder)
                file_dict=decoders[decoder_index]
                file_dict["encoder_params"]=decoder
                analysis_dict["file_list"].append(file_dict)
                
            json.dump(analysis_dict,open(os.path.join(final_run_path,"file_params.json"),"w+"),cls=NpEncoder)
            json.dump({"sequencing_data_path":path},open(os.path.join(final_run_path,"sequencing_params.json"),"w+"))

            try:
                os.symlink(path,os.path.join(final_run_path,"sequencing_data_link")) 
            except OSError as e: #force symlink
                if e.errno==errno.EEXIST:
                    os.remove(os.path.join(final_run_path,"sequencing_data_link"))
                    os.symlink(path,os.path.join(final_run_path,"sequencing_data_link"))

            print(final_run_path)
            command = "python "+ sequencing_tool_path + complete_param_string
            job.command = command
            job.run_path = final_run_path
            job.submit(args.no) #finally submit the experiment

