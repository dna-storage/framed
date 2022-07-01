from dnastorage.fi.fi_env import *
import dnastorage.fi.dna_processes as dna_process
from dnastorage.system.pipeline_dnafile import *
from dnastorage.system.formats import *
import generate
import logging
import numpy as np
import sys
import os
import time
import copy
import json
from joblib import Parallel, delayed
from multiprocessing import Queue
from logging.handlers import QueueHandler, QueueListener


logging_queue=Queue() #queue to enable parallel logging

def setup_main_logger():
    global logging_queue
    #set up handlers for multiprocessing logging
    handler=logging.StreamHandler()
    ql = QueueListener(logging_queue,handler)
    logger = logging.getLogger()
    logger.addHandler(handler)
    return ql

def setup_process_logger():
    global logging_queue
    qh=QueueHandler(logging_queue)
    logger=logging.getLogger()
    logger.addHandler(qh)



def _monte_kernel(monte_start,monte_end,args): #function that will run per process
    #We cant use the data_keeper object here, should be private per process, need data structures to propagate results back up to parent process
    results=[] #each element will be an object that encapsulates the entire statistics    
    stats.clear()
    file_to_fault_inject= open(args.file, "rb")
    file_to_fault_inject.seek(0,2)
    file_to_inject_size=file_to_fault_inject.tell()
    file_to_fault_inject.seek(0,0)

    header_data_path = os.path.join(args.out_dir,"header_pipeline{}.header".format(monte_end))
    
    if os.path.exists(args.enc_params):
        with open(args.enc_params,'rb') as param_file:
                encoding_params = json.load(param_file) #load in params for an encoding architecture through a json file
                
    if encoding_params is None:
        encoding_params={}

        
    if os.path.exists(args.header_params):
        with open(args.header_params,'rb') as param_file:
                header_params = json.load(param_file) #load in params for an encoding architecture through a json file


    if args.dna_process is not None and os.path.exists(args.dna_process):
        with open(args.dna_process,'rb') as param_file:
            strand_proc_dict = json.load(param_file)
                

    with open(args.fault_params,'rb') as param_file:
        fault_args = json.load(param_file)
    with open(args.distribution_params,'rb') as param_file:
        dist_args = json.load(param_file)

    assert fault_args !=None and dist_args !=None

    if monte_end!=args.num_sims: #save last header data incase we want to keep it
        write_dna = DNAFilePipeline.open(None,"w",format_name=args.arch,header_params=header_params,header_version=args.header_version,
                                         encoder_params=encoding_params,fsmd_header_filename=header_data_path,file_barcode=tuple(args.file_barcode))
    else:
        base_file_path  = os.path.basename(args.file)
        write_dna = DNAFilePipeline.open("{}.dna".format(base_file_path),"w",format_name=args.arch,header_params=header_params,header_version=args.header_version,
                                         encoder_params=encoding_params,fsmd_header_filename=header_data_path,file_barcode=tuple(args.file_barcode))
        
    write_dna.write(file_to_fault_inject.read())
    write_dna.close()

    
    stats["total_encoded_strands"]=len(write_dna.strands)
    stats["header_strand_length"]=len(write_dna.strands[0].dna_strand) #header strands should be first
    stats["payload_strand_length"]=len(write_dna.strands[-1].dna_strand)#strands with real payload should be at the end
    stats["dead_header"]=0
    stats["file_size_bytes"]=file_to_inject_size

    #do some basic processes on the DNA strands
    if strand_proc_dict is not None:
        #do some strand processing
        for s in write_dna.strands:
            for process in strand_proc_dict:
                func = getattr(dna_process,process)
                s.dna_strand = func(s.dna_strand,*strand_proc_dict[process])
        
    #Encode file we are fault injecting on
    fault_environment =  Fi_Env(args.strand_distribution, args.fault_model,write_dna.strands,dist_args,fault_args)
    
    for sim_number in range(monte_start,monte_end):
        logging.info("Monte Carlo Sim: {}".format(sim_number))
        fault_environment.Run()
        stats.inc("total_file_data_bytes",stats["file_size_bytes"])
        stats.inc("total_strands_analyzed",len(fault_environment.get_strands()))

        read_dna = DNAFilePipeline.open(None,"r",input_strands=fault_environment.get_strands(),header_params=header_params,
                                        header_version=args.header_version,format_name=args.arch,encoder_params=encoding_params,fsmd_header_filename=header_data_path,
                                        file_barcode=tuple(args.file_barcode))
        '''
        else:
            print("reading from file")
            read_dna = DNAFilePipeline.open("{}.dna".format(base_file_path),"r",header_params=header_params,
                                            header_version=args.header_version,format_name=args.arch,encoder_params=encoding_params,fsmd_header_filename=header_data_path,
                                            file_barcode=tuple(args.file_barcode))
        ''' 
        if read_dna is None:
            #dead header
            stats.inc("dead_header",1)
            stats.inc("error",1)
            stats.inc("total_mismatch_bytes",stats["file_size_bytes"])
            continue
        #calculate missing bytes
        file_to_fault_inject.seek(0,0)
        total_mismatch_data=0
        length_fi_data=0
        index=0
        while True:
            fi_data=read_dna.read(1)
            if len(fi_data)==0: break
            original_data=file_to_fault_inject.read(1)
            if len(original_data)==0: break
            length_fi_data+=1
            index+=1
            if fi_data==original_data:
                continue
            else:
                total_mismatch_data+=1
        stats.inc("total_mismatch_bytes",total_mismatch_data)
        stats.inc("file_size_difference_bytes",abs(file_to_inject_size-length_fi_data))
        if stats["total_mismatch_bytes"]>0:
            stats.inc("error",1)
        else:
            stats.inc("error",0)
        logging.info("Finished decoding erroneous file")
    results.append(copy.deepcopy(stats))
    
    if monte_end!=args.num_sims: #save last header data incase we want to keep it
        os.remove(header_data_path)
    file_to_fault_inject.close()    
    return results


def _monte_parallel_wrapper(monte_start,monte_end,args): #wrapper for parallelization
    setup_process_logger()
    return _monte_kernel(monte_start,monte_end,args)

#function for running monte carlo simulations for a fixed rate fault model
def run_monte(args):
    #run many simulations to perform statistical analysis
    stats_file_path =os.path.join(args.out_dir,"fi.stats")
    stats_pickle_path = os.path.join(args.out_dir,"fi.pickle")
    stats_fd = open(stats_file_path,'w+')
    pickle_fd=open(stats_pickle_path,'wb+')
    dist=[]
    monteIters=args.num_sims//args.cores
    last_monte_iter = args.num_sims%args.cores
    processIters=0
    tasks=[] #list of arguments for each task
    parallel=Parallel(n_jobs=args.cores,backend="multiprocessing")
    for i in range(args.cores):
        if i==args.cores-1:
            tasks.append((processIters,processIters+monteIters+last_monte_iter,args))
        else:
            tasks.append((processIters,processIters+monteIters,args))
        processIters+=monteIters #set up next range of montecarlo simulations
    results=parallel(delayed(_monte_parallel_wrapper)(*t) for t in tasks )
    total_results=[]
    for job_res in results:#unpack results from jobs
        for res in job_res:
            total_results.append(res)

    stats.clear()
    stats.set_fd(stats_fd)
    stats.set_pickle_fd(pickle_fd)
    #need to aggregate data across all results
    for s in total_results:
        stats.aggregate(s,["total_encoded_strands","header_strand_length","payload_strand_length",
                           "file_size_bytes"]) #just want to copy information about total strands/strand_length
    stats["total_runs"]=args.num_sims
    stats.persist()
    stats_fd.close()
    pickle_fd.close()



'''
Main file for injecting faults into an input DNA file, 3 available options are available

miss_strand --- missing strand fault model, e.g. strand is not available to the decoder
strand_fault --- faults within strand fault model, e.g. insert deletions, insertions, substitutions
fixed_rate --- fixed fault rate applied to each nucleotide in each strand.
Combo fault mode will not be supported immediately for fault injection due to the current run-time constraints
miss_strand will eventually be supported

'''



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Inject faults into input file and perform analysis")
    parser.add_argument('--fault_model',required=True,choices=fault_injection_modes(),help="fault modes to use")
    parser.add_argument('--strand_distribution',required=True,choices=distribution_functions(),help="distribution function to use")    
    parser.add_argument('--simulation_runs',dest="num_sims",action="store",type=int,default=1000,help="Number of simulations to run")
    
    parser.add_argument('--arch',required=True,choices=file_system_formats(),help="Encoding/decoding architecture")
    parser.add_argument('--file',required=True,help = "File to fault inject on")
    parser.add_argument('--cores', type=int, action="store",default=1,help="Number of threads to run monte carlo simulations")

    parser.add_argument('--enc_params',type=str,required=True,action="store",help="Path to json file with parameters used to set error correction in encoding")
    parser.add_argument('--header_params',type=str,required=True,action="store",help="Path to json file with parameters used to set error correction in encoding for the header")
    parser.add_argument('--header_version',type=str,required=False,action="store",default="0.1",help="Header version")
    parser.add_argument('--distribution_params',required=True,action="store",help="Parameters used to specify the distribution")
    parser.add_argument('--fault_params',required=True,action="store",help="Parameters used to specify the fault model")
    parser.add_argument('--dna_process',required=False,default=None,help="set of processing steps to do on dna strands before fault injection")
    parser.add_argument('--out_dir',type=str,required=True,action="store",help="Directory where data will be dumped")
    parser.add_argument('--file_barcode',required=False,default=tuple(),nargs="+",type=int,help="Barcode for the file")
    

    args = parser.parse_args()

    assert os.path.isdir(args.out_dir)

    logging.basicConfig(filename=os.path.join(args.out_dir,"fi_info.log"),level=logging.INFO,
                        format='%(asctime)s.%(msecs)03d %(process)s %(levelname)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    generate.seed()
    queue_listener = setup_main_logger()
    queue_listener.start()
    run_monte(args)
    queue_listener.stop()
