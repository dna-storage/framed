"""
Main script for analysis on real sequencing data through the use of encode/decode pipelines.
"""
from dnastorage.system.pipeline_dnafile import *
from dnastorage.system.formats import *
from dnastorage.util.strandinterface import *
from dnastorage.util.mpi_logger import *
from dnastorage.util.mpi_utils import *
from dnastorage.codec.base_conversion import convertIntToBytes,convertBytesToInt

import logging
import numpy as np
import sys
import os
import time
import copy
import json
import numpy as np
from mpi4py import MPI

logger = logging.getLogger()                                                                                                                                     

def is_master(comm): return comm.rank==0

if __name__=="__main__":
    import argparse

    world_comm=MPI.COMM_WORLD
 
    parser = argparse.ArgumentParser(description="Analyze sequencing data for given encoders")
    parser.add_argument('--dna_file_params',type=str,required=True,action="store",help="Path to json file that describes the encoder parameters for files in the sequencing data")
    parser.add_argument('--out_dir',type=str,required=True,action="store",help="Directory where data will be dumped")
    parser.add_argument('--sequencing_data_path',required=True,action="store",help="Path to sequencing data")
    
    args = parser.parse_args()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')                                                                                           
    mpi_handler = MPIFileHandler(os.path.join(args.out_dir,"sequencing_info_{}.log".format(world_comm.rank)))
    mpi_handler.setFormatter(formatter)                                                                                                                                             
    logger.addHandler(mpi_handler)                                                                                                                                                  
    logger.setLevel(logging.DEBUG)
    assert os.path.isdir(args.out_dir)                                                                                                                                             

    strand_interface = BaseStrandInterface()
    if is_master(world_comm): #only rank 0 loading in data
        #try to figure out what format the sequencing data is in through some basic checks
        try:
            if os.path.isfile(args.sequencing_data_path):
                if ".fastq" in args.sequencing_data_path:
                    strand_interface = BaseStrandInterface.open("fastq",args.sequencing_data_path)
                elif ".fast5" in args.sequencing_data_path:
                    strand_interface = BaseStrandInterface.open("fast5",args.sequencing_data_path)
            elif os.path.isdir(args.sequencing_data_path):
                strand_interface = BaseStrandInterface.open("fast5",args.sequencing_data_path)
        except Exception as e:
            logger.fatal("Could not find sequencing data: {}".format(e))
            exit(1)
        stats["total_sequencing_strands"]=len(strand_interface.strands)
        logger.info("Total Reads Read in {}".format(len(strand_interface.strands)))
    #At this point, we have sequencing data loaded into strand interface, now load dna file params to launch
    try:
        with open(args.dna_file_params,'rb') as param_file:
            dna_file_params = json.load(param_file)
        if "file_list" not in dna_file_params: raise ValueError("file_list not in sequencing parameter file")
    except Exception as e:
        logger.fatal("Could not open and read encoding parameter_file: {}".format(e))
        exit(1)

    file_list = dna_file_params["file_list"]
    
    for file_index,file_params in enumerate(dna_file_params["file_list"]): #Inner loop, allows for multiple encoders to analyze the sequencing data in one run, could be useful for multi-filed data sets
        #load up required values, throw exceptions as necessary
        try:
            arch = file_params["arch"]
            header_version = file_params["header_version"]
            encoder_params = file_params["encoder_params"]
            header_params = file_params["header_params"]
            payload_header_file = file_params.get("payload_header_path",None) #if this does not exist, need to be able to decode it in DNAFileInterface
            header_header_file = file_params["header_header_path"] #this needs to exist
            barcode = tuple(file_params.get("barcode",tuple()))
            if not isinstance(encoder_params,dict): raise ValueError("encoder_params not a dictionary")
            if not isinstance(header_params,dict): raise ValueError("header_params not a dictionary")
            if "barcode" in file_params: #if barcoded, make sure to suffix the pipeline titles just so that data is kept separate
                barcode_bytes = file_params["barcode"]
                barcode_integer = convertBytesToInt(barcode_bytes)
                file_suffix =str(barcode_integer)
                encoder_params["title"] = "{}_{}".format(encoder_params.get("title",""),file_suffix)
                header_params["title"] = "{}_{}".format(header_params.get("title",""),file_suffix)
        except Exception as e:
            logger.fatal("Fatal issue in a file_params dictionary : {}".format(e))
            exit(1)
        #Now, construct the DNA file and allow the decoder to run
        logger.info("Analyzing File {} ".format(file_index))
        read_dna = DNAFilePipeline.open("r",format_name=arch,header_params=header_params,
                                        header_version=header_version,encoder_params=encoder_params,
                                        fsmd_header_filename=header_header_file,file_barcode=barcode,
                                        payload_header_filename = payload_header_file,mpi=world_comm,strand_interface=strand_interface)

        #recoordinate strands
        gather_strands=object_gather(read_dna.get_returned_strands(),world_comm)
        strand_interface.strands=gather_strands
        logger.info("Recoordinated strands")
        
    #merge together stats that were collecting during decoding for all processes
    logger.info("Gathering stats")
    all_stats = world_comm.gather(stats,root=0)
    logger.info("Stats have been gathered")
    if is_master(world_comm):
        stats.clear()
        for s in all_stats: stats.aggregate(s)
        stats["strands_not_indexed"]  = len(strand_interface.strands)
        filtered_strand_lengths =[len(s.dna_strand) for s in strand_interface.strands]
        stats["unindexed_length:average"]=np.average(filtered_strand_lengths)
        stats["uindexed_length:stdev"]=np.std(filtered_strand_lengths)
        #write out statistics
        stats_file_path =os.path.join(args.out_dir,"sequencing.stats")
        stats_pickle_path = os.path.join(args.out_dir,"sequencing.pickle")
        stats_fd = open(stats_file_path,'w+')
        pickle_fd=open(stats_pickle_path,'wb+')
        stats.set_fd(stats_fd)
        stats.set_pickle_fd(pickle_fd)
        stats.persist()
        stats_fd.close()
        pickle_fd.close()
    
    
