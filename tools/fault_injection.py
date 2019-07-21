from dnastorage.fi import fault_injector
from dnastorage.codec import dense
from dnastorage.codec import commafreecodec
from dnastorage.codec import norepeatscodec
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec.rscodec import *
from dnastorage.arch.strand import *
from dnastorage.arch.builder import *
from dnastorage.primer.primer_util import *
from dnastorage.util.neg_binomial_gen import *
from dnastorage.handle_strands.strand_handlers import *
from dnastorage.handle_strands.cluster import *
from dnastorage.util.data_fi import *
import numpy as np
import sys
import os
import time
from joblib import Parallel, delayed

#!/usr/bin/env python

def build_decode_xx_architecture(arch, pf, primer5, primer3, fountain_table=None):
    if arch == "UW+MSv1":
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 20
        dec = DecodeXORStrand(pf,p)
        return dec

    elif arch == "Illinois":
        illini = illinois.IllinoisCodec(primer5,150,Checksum())
        p = StrandPrimers(primer5, primer3, illini)
        pf.packetSize = 149
        enc = DecodeNaiveStrand(pf,p)
        return enc

    elif arch == 'Goldman':
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 5
        dec = DecodeGoldmanStrand(pf,4,p)
        return dec

    elif arch == 'Binary':
        b = binary.BinaryRotateCodec(binary.BinaryCodec(Checksum()))
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 17
        enc = DecodeNaiveStrand(pf,p)
        return enc

    elif arch == 'Fountain':
        #assert False and "Not fully implemented"
        h = huffman.RotateCodec(huffman.HuffmanCodec(23,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 22
        enc = DecodeFountainStrand(pf,fountain_table,p)
        return enc

    elif arch == 'NCState':
        assert False and "Not fully implemented"
        return None

    elif arch == 'NRDense':
        #b = norepeatscodec.NoRepeatCodec(Checksum(2))
        b = norepeatscodec.NoRepeatCodec(Checksum(2,CheckLength(28)))
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 28
        dec = DecodeNaiveStrand(pf,p)
        return dec

    elif arch == 'RS+CFC8':
        b = commafreecodec.CommaFreeCodec(13,None,2)
        p = StrandPrimers(primer5, primer3, b)
        dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        return dec
    elif arch == 'RS+ROT':
        h = huffman.RotateCodec(huffman.HuffmanCodec(11,None,11,99))
        p = StrandPrimers(primer5, primer3, h)
        enc = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        return enc


def read_header(dec_file):
    f = dec_file
    header = {}
    while True:
        l = f.readline()
        if l.startswith('%') and 'bytes encoded' in l:
            words = l[1:].split(' ')
            try:
                size = int(words[0])
                header['size'] = size
            except:
                continue
        elif l.startswith('%') and "- primer 5' end" in l:
            words = l[1:].split(' ')
            header['primer5'] = words[0]
        elif l.startswith('%') and "- primer 3' end" in l:
            words = l[1:].split(' ')
            header['primer3'] = words[0]
        elif not l.startswith('%'):
            break

    return header

#function that will add arguments to args by reading the header if needed, also returns a clean packetizedFile, clean Decoder object and a set of clean strands
#required for determining correctness when compared to faulty runs
def clean_run(args):
    if args.input_file != sys.stdin:
        h = read_header(args.input_file)

    if args.filesize == 0:
        if h.has_key('size'):
            args.filesize = h['size']
        else:
            print "Please provide the file's size."
            sys.exit(-1)

    if args.primer3 == None:
        if h.has_key('primer3'):
            args.primer3 = h['primer3']
        else:
            print "Please provide primer3."
            sys.exit(-1)

    if args.primer5 == None:
        if h.has_key('primer5'):
            args.primer5 = h['primer5']
        else:
            print "Please provide primer3."
            sys.exit(-1)

    table = []
    if args.arch == 'Fountain':
        if args.input_file == sys.stdin:
            tfile = open("default.table","r")
        else:
            tfile = open(args.input_file.name+".table","r")
        if tfile == None:
            print "Please fountain table."
            sys.exit(-1)
        while True:
            l = tfile.readline()
            if l == "":
                break
            x = [int(y) for y in l.split(',')]
            table.append( [x[0],x[1:]] )
        tfile.close()


    clean_packetizedFile = WritePacketizedFilestream(args.o,args.filesize,20)

    clean_Decoder = build_decode_architecture(args.arch, clean_packetizedFile, args.primer5, args.primer3,table)

    ii = 0
    #initial_strands will hold the initial strands of the input file before fault injection
    initial_strands=[]
    while not clean_Decoder.complete:
        s = args.input_file.readline()
        ii += 1
        if len(s) == 0:
            break
        s = s.strip()
        if s.startswith('%'):
            continue

        initial_strands.append(s)
        clean_Decoder.decode(s)

    clean_file=clean_Decoder.dummy_write()

    return clean_file,initial_strands,clean_Decoder


#use the negative binomial distribution to generate a new pool of strands
#will only return 1 copy of each strand if negative binomial distribution is disabled
def distribute_reads(strands,neg_bin_randomizer,tuple_format=True):
    read_count=[]
    new_pool=[]
    dist=[]
    #Tuple format accelerates simulation of fault models of nucleotide strand faults and missing strand faults
    if tuple_format==False:
        for ind,s in enumerate(strands):
            if neg_bin_randomizer is not None:
                read_count.append(neg_bin_randomizer.gen())
            else:
                read_count.append(1)
            dist.append((ind,read_count[ind]))
        for index, s in enumerate(strands):
            new_pool+=[s for i in range(0,read_count[index])]
        return new_pool,dist
    else:
        if neg_bin_randomizer is not None:
            #make a simple array of tuples in format (strand,count for that strand)
            new_pool=[(s,neg_bin_randomizer.gen()) for s in strands]
        else:
            new_pool=[(s,1) for s in strands]
        pool_size=0
        for strand in new_pool:
            pool_size+=strand[1]
        return new_pool,pool_size

#builds the correct class for handling a set of dirty strands
def build_strand_handler(strand_handler,Codec,cluster_algorithm):
    if strand_handler == "data_vote_simple":
        return data_vote_simple(Codec)
    elif strand_handler == "cluster_BMA":
        return clustering_handler(Codec,"BMA",cluster_algorithm)
    elif strand_handler=="cluster_BMA_ED":
        return clustering_handler(Codec,"BMA_ED",cluster_algorithm)
    elif strand_handler=="cluster_ED":
        return clustering_handler(Codec,"ED",cluster_algorithm)

#builds the class used by strand handlers that contains the interfaces
#to the clustering algorithms
def build_cluster_algorithm(cluster_algorithm):
    if cluster_algorithm=="starcode_MP":
        return starcode("MP")


#function to wrap the process of running many strand_fault and missing_strand fault model simulations
def run_monte_MS(args,desired_faulty_count,clean_strands,clean_file,data_keeper,strand_handler,fault_model,desired_faults_per_strand=None):
    #run many simulations to perform statistical analysis
    table=[]

    if desired_faults_per_strand is not None:
        fault_model.set_faults_per_strand(desired_faults_per_strand)


    for sim_number in range(0,args.num_sims):

        #build a "dirty" decoder and packetizedFile for each run
        dirty_packetizedFile= WritePacketizedFilestream(args.o,args.filesize,20)
        dirty_Decoder = build_decode_architecture(args.arch, dirty_packetizedFile, args.primer5, args.primer3,table)


        #get a new pool of strands based on the negative binomial distribution
        multiple_strand_pool,pool_size=distribute_reads(clean_strands,read_randomizer)


        # make sure we do not select more strands than pool size
        if pool_size<desired_faulty_count:
            fault_model.set_faulty_missing(pool_size,pool_size)
        else:
            fault_model.set_faulty_missing(desired_faulty_count,desired_faulty_count)
        #set the fault_model's library to the new strand pool
        fault_model.set_library(multiple_strand_pool)


        strands_after_faults=fault_model.Run_tuple()
        #print len(strands_after_faults)
        #call the strand handler, will return either key,value pairs or strands
        processed_strands=strand_handler.process_tuple(strands_after_faults)

        #Hand the processed strands off to decoding
        for proc in processed_strands:
            if type(proc) is tuple:
                dirty_Decoder.decode(None,bypass=True,input_key=proc[0],input_value=proc[1])
            else:
                dirty_Decoder.decode(proc)
        #perform a dummy write
        bad_file=dirty_Decoder.dummy_write()
        percent_correct=data_keeper.calculate_correctness(bad_file,clean_file)
        data_keeper.insert_correctness_result(percent_correct)
        
    #grab an example of the distribution that was used
    data_keeper.set_distribution(multiple_strand_pool)
    #keep track of correctness results
    lower,middle,upper=data_keeper.calculate_midpoint()
    prob=data_keeper.get_probability_value()
    if desired_faults_per_strand is None:
        data_point=(desired_faulty_count,lower,middle,upper)
        data_keeper.insert_probability_point((desired_faulty_count,prob))
    else:
        data_point=(desired_faults_per_strand,desired_faulty_count,lower,middle,upper)
        data_keeper.insert_probability_point((desired_faults_per_strand,desired_faulty_count,prob))
    data_keeper.insert_data_point(data_point)
    data_keeper.clear_correctness_results()
    data_keeper.clear_probability_results()





def _monte_rate_kernel(args,clean_strands,clean_file,data_keeper,strand_handler,fault_model,error_rate,monte_start,monte_end): #function that will run per process
    #We cant use the data_keeper object here, should be private per process, need data structures to propagate results back up to parent process
    results=[] #each element will be a tuple ---> ((float) percent_correct, (bool) File_Correct)
    table=[]
    for sim_number in range(monte_start,monte_end+1):
        print "Monte Carlo Sim: {}".format(sim_number)
        #build a "dirty" decoder and packetizedFile for each run
        dirty_packetizedFile= WritePacketizedFilestream(args.o,args.filesize,20)
        dirty_Decoder = build_decode_architecture(args.arch, dirty_packetizedFile, args.primer5, args.primer3,table)
        #get a new pool of strands based on the negative binomial distribution
        multiple_strand_pool,dist=distribute_reads(clean_strands,read_randomizer,False)

        #set the fault_model's library to the new strand pool
        fault_model.set_library(multiple_strand_pool)
        strands_after_faults=fault_model.Run()

        #call the strand handler, will return either key,value pairs or strands
        processed_strands=strand_handler.process(strands_after_faults)

        #Hand the processed strands off to decoding
        for proc in processed_strands:
            if type(proc) is tuple:
                dirty_Decoder.decode(None,bypass=True,input_key=proc[0],input_value=proc[1])
            else:
                dirty_Decoder.decode(proc)
        print "Finished decoding erroneous file"
        #perform a dummy write
        bad_file=dirty_Decoder.dummy_write()
        results.append(data_keeper.calculate_correctness(bad_file,clean_file))
    return results

def _monte_rate_parallel_wrapper(args): #wrapper for parallelization
    return _monte_rate_kernel(*args)

#function for running monte carlo simulations for a fixed rate fault model
def run_monte_rate(args,clean_strands,clean_file,data_keeper,strand_handler,fault_model,error_rate):
    #run many simulations to perform statistical analysis
    dist=[]
    fault_model.set_fault_rate(error_rate)
    monteIters=args.num_sims/args.cores
    processIters=0
    tasks=[] #list of arguments for each task
    parallel=Parallel(n_jobs=args.cores)
    for i in range(args.cores):
        tasks.append((args,clean_strands,clean_file,data_keeper,strand_handler,fault_model,error_rate,processIters,processIters+(monteIters-1)))
        processIters+=monteIters #set up next range of montecarlo simulations
    
    print "Running Monte Carlo sim for error rate {}".format(error_rate)
    
    results=parallel(delayed(_monte_rate_parallel_wrapper)(t) for t in tasks )
    _r=[]
    for job_res in results:#unpack results from jobs
        for res in job_res:
            _r.append(res)
    #add results to the data keeper after the parallel jobs are finished
    for r in _r:
        data_keeper.insert_correctness_result(r[0])
        if r[1]: data_keeper.inc_num_correct()
    
    #keep track of correctness results
    lower,middle,upper=data_keeper.calculate_midpoint()
    data_point=(error_rate,lower,middle,upper)
    data_keeper.insert_data_point(data_point)
    data_keeper.clear_correctness_results()
    prob=data_keeper.get_probability_value()
    data_keeper.insert_probability_point((error_rate,prob))
    data_keeper.clear_probability_results()



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
    parser.add_argument('--fault_rate_file', dest="fault_file",action="store", default=None, help='file to have fault rate per nucleotide')
    parser.add_argument('--fault_model',required=True,choices=['miss_strand','strand_fault','combo','fixed_rate'])
    parser.add_argument('--run',action="store_true",help="Should strand_fault faults be in a run")
    parser.add_argument('--faulty_count',dest="faulty",action="store",default='10-10',help="range to sweep the number of faulty/missing strands")
    parser.add_argument('--fail_count',dest="fails",action="store",default='2-2',help="range to sweep the number of nucleotide errors ")
    parser.add_argument('--faulty_step',dest="faulty_step",action="store",type=int,default=1,help="Step size through strand range")
    parser.add_argument('--fail_step',dest="fails_step",type=int,action="store",default=1,help="Step size through nucleotide range")
    parser.add_argument('--fault_rate',dest="rate",action="store",default='0.001-0.001',help="range to sweep the fault rate")
    parser.add_argument('--rate_step',dest="rate_step",type=float,action="store",default=0.001,help="Step size through fault rate range")
    parser.add_argument('--filesize',type=int,dest="filesize",action="store",default=0, help="Size of file to decode.")
    parser.add_argument('--primer5',dest="primer5",action="store",default=None, help="Beginning primer.")
    parser.add_argument('--primer3',dest="primer3",action="store",default=None, help="Ending primer.")
    parser.add_argument('--primer1_length',dest="p1",action="store",type=int,default=20,help="Length of primer1 (total, so if hierarchy set to 40 instead of 20)")
    parser.add_argument('--primer2_length',dest="p2",action="store",type=int,default=20,help="Length of primer2, for all our experiments this will be set to 20")
    parser.add_argument('--seq_mean',dest="mean",action="store",type=float,default=0,help="Mean of the sequencing distribution,set to 0 to deactivate distribution")
    parser.add_argument('--seq_var',dest="var",action="store",type=float,default=4.3,help="Variance of the sequencing distribution, set to 0 to deactivate distribution")
    parser.add_argument('--simulation_runs',dest="num_sims",action="store",type=int,default=10000,help="Number of simulations to run")
    parser.add_argument('--o',nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="Output file.")
    parser.add_argument('--arch',required=True,choices=available_file_architectures())
    parser.add_argument('input_file', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Clean strands file')
    parser.add_argument('--strand_handler',required=True,choices=['data_vote_simple','cluster_BMA','cluster_BMA_ED','cluster_ED'])
    parser.add_argument('--cluster_algorithm',required=False,choices=['starcode_MP'],default=None)
    parser.add_argument('--fileID',action="store",default='1',help="ID for the input file")
    parser.add_argument('--cores', type=int, action="store",default='1',help="Number of threads to run monte carlo simulations")


    args = parser.parse_args()

    #Parse the ranges for the nucleotide errors count and number of strands to be selected
    nuc_lower,nuc_upper=args.fails.split('-')
    nuc_lower=int(nuc_lower)
    nuc_upper=int(nuc_upper)+1
    strand_lower,strand_upper=args.faulty.split('-')
    strand_lower=int(strand_lower)
    strand_upper=int(strand_upper)+1

    #parse the range for the fault rate range
    rate_lower,rate_upper=args.rate.split('-')
    rate_lower=float(rate_lower)
    rate_upper=float(rate_upper)+args.rate_step

    model_name= args.fault_model

    #set up fault injection arguments, note this is just to initialize some variables in the fault model object
    #things like injector.fails/faulty/rate will be changed durin simulation sweeps
    injector_args=fault_injector.fault_injector_arguments()
    injector_args.o=args.o
    injector_args.input_file=None
    injector_args.fault_file=args.fault_file
    injector_args.fault_model=args.fault_model
    injector_args.run=args.run
    injector_args.faulty=0
    injector_args.fails=0
    injector_args.p1=args.p1
    injector_args.p2=args.p2

    # check if the selected fault model class exists
    if not hasattr(fault_injector, model_name):
        print "Couldn't run fault model '%s'. Class '%s' doesn't exist in fault_injector.py" % \
        (args.fault_model, model_name)
        sys.exit(1)
    #clean_Decoder will buffer the clean_packetizedFile we will use to compare with the faulty file
    clean_file,clean_strands,clean_Decoder=clean_run(args)
    # instantiate a fault model object
    fault_model = eval('fault_injector.'+model_name+'(injector_args)')
    #instantiate the negative binomial random number generator, used for distributing the number of reads
    if(args.mean !=0 and args.var !=0):
        read_randomizer=neg_bin(args.mean,args.var)
    else:
        read_randomizer=None


    
    
    #instantiate the clustering algorithm
    clustering_algorithm=build_cluster_algorithm(args.cluster_algorithm)
    #instantiate the strand handler
    strand_handler=build_strand_handler(args.strand_handler,clean_Decoder._Codec, clustering_algorithm)   

    
    
    #determine which model to run
    if args.fault_model == "strand_fault" and args.fault_file is None:
        #initialize the data helper object
        data_keeper=data_helper(args.fileID,args.arch,args.strand_handler,args.mean,args.var,args.fault_model,args.num_sims,strand_range=args.faulty,nuc_range=args.fails,rate_range=args.rate)
        for n in range(nuc_lower,nuc_upper,args.fails_step):
            for s in range(strand_lower,strand_upper,args.faulty_step):
                run_monte_MS(args,s,clean_strands,clean_file,data_keeper,strand_handler,fault_model,n)

    elif args.fault_model == "fixed_rate":
        #initialize the data helper object
        data_keeper=data_helper(args.fileID,args.arch,args.strand_handler,args.mean,args.var,args.fault_model,args.num_sims,strand_range=args.faulty,nuc_range=args.fails,rate_range=args.rate)
        #run a set of simulations for each failure rate
        for error_rate in np.arange(rate_lower,rate_upper,args.rate_step):
             run_monte_rate(args,clean_strands,clean_file,data_keeper,strand_handler,fault_model,error_rate)

    data_keeper.plot_pmf(read_randomizer) #make a plot of the probability mass function
    data_keeper.dump()
