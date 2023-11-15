#Classes for each of the fault models
import os
import random
import csv
import pickle
import re
import time
from Bio import SeqIO
import scipy
import numpy as np
import math

from dnastorage.util.stats import stats
import dnastorage.util.generate as generate
from dnastorage.strand_representation import *
from dnastorage.fi.fault_strand_representation import *

import logging
logger = logging.getLogger("dnastorage.fi.fault_injector")
logger.addHandler(logging.NullHandler())


#substitution dictionary used for substitution errors
sub_dict={'A':['G','C','T'],'G':['C','A','T'], 'T':['G','C','A'], 'C':['G','T','A']}
#nucleotide list for insertion errors
nuc_list=['A','C','T','G']

#Base class with some common functions to all of the fault models
class BaseFI:
    def __init__(self):
        pass
    @classmethod
    def open(self,fault_injector,**kwargs):
        if fault_injector=="fixed_rate":
            return fixed_rate(**kwargs)
        elif fault_injector=="position_fixed_rate":
            return position_fixed_rate(**kwargs)
        elif fault_injector=="strand_fault_compressed":
            return strand_fault_compressed(**kwargs)
        elif fault_injector=="sequencing_experiment":
            return sequencing_experiment(**kwargs)
        elif fault_injector=="sequencing_experiment_downsample":
            return sequencing_experiment_downsample(**kwargs)
        elif fault_injector=="sequencing_experiment_coverage":
            return sequencing_experiment_coverage(**kwargs)
        elif fault_injector=="pattern_fixed_rate":
            return pattern_fixed_rate(**kwargs)
        elif fault_injector=="DNArSim":
            return DNArSim(**kwargs)
        else:
            raise ValueError()
    #These two setting functions allow easier altertion of the input library to fault injection, and parameters around fault injection
    def set_library(self,input_strands):
        self._input_library=input_strands
    def Run(self):
        raise NotImplementedError()    
    def read_csv(self,file_name):
        _file=open(file_name,'r')
        csv_parsed=csv.reader(_file,delimiter=',')
        return csv_parsed
 
#use sequencing data to perform fault injection
class sequencing_experiment(BaseFI):
    def __init__(self,**args):
        BaseFI.__init__(self)
        if "sequencing_data_path" in args and "mapping_path" in args:
            assert os.path.exists(args["sequencing_data_path"]) and os.path.exists(args["mapping_path"])
            self._sequencing_data_path = args["sequencing_data_path"]
            self._map  = pickle.load(open(args["mapping_path"],"rb")) #map: sequencing_record-->index_ints
        elif "experiment_path" in args:
            #need to derive sequencing_data_path and mapping_path from an experiment_path
            assert os.path.exists(args["experiment_path"]) and os.path.isdir(args["experiment_path"])
            self._sequencing_data_path = os.path.join(args["experiment_path"],"sequencing_data_link")
            map_path = os.path.join(args["experiment_path"],"sequencing.map.pickle")
            assert os.path.exists(map_path) and os.path.exists(self._sequencing_data_path)
            self._map=pickle.load(open(map_path,"rb"))
        else:
            raise ValueError("Invalid arguments for sequencing_experiment fault injection")
        tmp_map = {}
        for sub_map in self._map: #push things into one map that should be indexed by index_ints of the DNAStrand
            tmp_map={**tmp_map,**self._map[sub_map]}
        self._map=tmp_map      
    def sequence_run_injection(self):
        out_list=[]
        input_library_map={}
        for strand in self._input_library:
            input_library_map[strand.index_ints] = strand
        #assume fastq file for now
        record_counter=0
        for record in SeqIO.parse(self._sequencing_data_path,"fastq"):
            record_counter+=1
            seq_strand,subs = re.subn("U","T",str(record.seq))
            seq_strand,n_subs = re.subn("N","A",seq_strand)
            seq_strand=seq_strand.rstrip('\x00')
            index_ints = self._map.get(record.id,None)
            if index_ints is None or index_ints not in input_library_map: continue #dead strand or is unwanted by the given studied file
            lib_strand=input_library_map[index_ints]
            out_list.append(FaultDNA(lib_strand,seq_strand))
            out_list[-1].record_id=record.id
        stats["sequencing_experiment::total_seq_file_strands"] = record_counter
        assert len(out_list)>0
        return out_list
    def Run(self):
        return self.sequence_run_injection()
    
class sequencing_experiment_coverage(sequencing_experiment):
    def __init__(self,**args):
        sequencing_experiment.__init__(self,**args) 
        self._coverage = args["coverage"]
        self._rng = np.random.default_rng()
    def sequence_run_injection(self):
        out_list=[]
        input_library_map={}
        for strand in self._input_library:
            input_library_map[strand.index_ints] = strand
        strands_to_sample = self._coverage*len(self._input_library)
        records=[]
        for record in SeqIO.parse(self._sequencing_data_path,"fastq"):
            records.append(record)
        records_to_use = self._rng.choice(np.array(records,dtype=object),int(strands_to_sample))
        for record in records_to_use:
            seq_strand,subs = re.subn("U","T",str(record.seq))
            seq_strand,n_subs = re.subn("N","A",seq_strand)
            seq_strand=seq_strand.rstrip('\x00')
            index_ints = self._map.get(record.id,None)
            if index_ints is None or index_ints not in input_library_map: continue #dead strand or is unwanted by the given studied file
            lib_strand=input_library_map[index_ints]
            out_list.append(FaultDNA(lib_strand,seq_strand))
        return out_list
    def Run(self):
        return self.sequence_run_injection()

#builds on sequencing experiment FI to consider a subset of the strands that are mapped
class sequencing_experiment_downsample(sequencing_experiment):
    def __init__(self,**args):
        sequencing_experiment.__init__(self,**args)
        self._down_sample = args.get("down_sample",None) #for percentages, down_sample is the percent you want removed, for integer inputs, downsample is what you want to keep
        if type(self._down_sample) is str: 
            if not os.path.exists(self._down_sample):
                raise ValueError("down_sample should point to a path indicating what down_sample the experiment should be subject to")
            experiment_downsample_file = open(self._down_sample,'r')
            for line in experiment_downsample_file.readlines():
                experiment_name = line.split()[0]
                if experiment_name in self._sequencing_data_path:
                    self._down_sample=float(line.split()[1])
                    break
            if not type(self._down_sample) is float: raise ValueError("down sample file did not include entry that matched input sequencing data path") 
        else: #Use a raw value 
            self._down_sample=float(self._down_sample)
        assert type(self._down_sample) is float  and self._down_sample>=0 #should have a downsample value at this point
        self._rng = np.random.default_rng(seed=0)
    def Run(self):
        strands = sequencing_experiment.Run(self)
        if self._down_sample is None: return strands #no downsampling
        self._rng.shuffle(strands)
        if self._down_sample<=1.0:
            keep_percent = 1.0-self._down_sample
            strands_to_keep = int(math.ceil(keep_percent*len(strands))) #round up the number of strands to keep
        elif self._down_sample>=1.0:
            strands_to_keep=int(math.ceil(self._down_sample)) #just keep this number of strands
        return strands[0:strands_to_keep]

    
#this class applies a fixed error rate to each nucleotide
#Strands are a list, not a tuple representation
class fixed_rate(BaseFI):
    def __init__(self,**args):
        BaseFI.__init__(self)
        self.fault_rate=args.get("fault_rate",0)
        
    def Run(self):
        self._injection={}
        self._injection=self.injection_sites()
        return self.inject_faults(self._injection)

    #go through each nucleotide in each strand and apply a flat fault rate
    def injection_sites(self):
        injection_sites={}
        #go through the strands and pick faults
        for strand_index,strand in enumerate(self._input_library):
            injection_sites[strand_index]=[]
            for nuc_index,nuc in enumerate(strand.dna_strand):
                inject_fault=generate.rand()
                if inject_fault<=self.fault_rate:
                    #need to inject a fault at this nucleotide, fault types are equally probable
                    injection_sites[strand_index].append((nuc_index,str(generate.rand_in_range(0,2))))
        return injection_sites
                    
    #inject errors in the list of strands, each 
    def inject_faults(self,inject_sites):
        out_list=self._input_library[:]
        for strand_indexes in inject_sites:
            library_strand=self._input_library[strand_indexes].dna_strand
            new_strand = self._input_library[strand_indexes].dna_strand
            for (fault_indexes,error) in sorted(inject_sites[strand_indexes],reverse=True,key=lambda x: x[0]):
                #substitution error
                if error == '0':
                    #chose a random nucleotide that is different from the current one
                    sub_nucleotide=random.choice(sub_dict[new_strand[fault_indexes]])
                    #add on some extra information to the injection sites that indicates the nucleotide used for substitution
                    new_strand=new_strand[0:fault_indexes]+sub_nucleotide+new_strand[fault_indexes+1:len(new_strand)]
                #deletion error
                elif error == '1':
                    #add on some extra information to the injection sites, append the nucleotide that was removed from the original strand 
                    new_strand=new_strand[0:fault_indexes]+new_strand[fault_indexes+1:len(new_strand)]
                #insertion error
                elif error == '2':
                    insert_nucleotide=random.choice(nuc_list)
                    new_strand=new_strand[0:fault_indexes]+insert_nucleotide+new_strand[fault_indexes:len(new_strand)]
                else:
                    raise ValueError()
            assert len(new_strand)>0
            out_list[strand_indexes]=FaultDNA(self._input_library[strand_indexes],new_strand)
        return out_list


class position_fixed_rate(fixed_rate): #does fixed rate error rates, except on a per-base basis based on experimental data
    def __init__(self,**args):
        fixed_rate.__init__(self,**args)
        path = args.get("error_rate_path",None)
        try:
            self._rate_data = pickle.load(open(path,"rb"))
            logger.info(self._rate_data)
        except Exception as e:
            logger.warning(e)
            exit()
    #go through each nucleotide in each strand and apply a flat fault rate
    def injection_sites(self):
        injection_sites={}
        #go through the strands and pick faults
        for strand_index,strand in enumerate(self._input_library):
            injection_sites[strand_index]=[]
            for nuc_index,nuc in enumerate(strand.dna_strand):
                inject_fault=generate.rand()
                if inject_fault<=self._rate_data[nuc_index]:
                    #need to inject a fault at this nucleotide, fault types are equally probable
                    injection_sites[strand_index].append((nuc_index,str(generate.rand_in_range(0,2))))
        return injection_sites


class pattern_fixed_rate(fixed_rate): #allows the injection of patterns rather than single errors
    def __init__(self,**args):
        fixed_rate.__init__(self,**args)
        rate_path = args.get("error_rate_path",None) #path to rates for each base indicating probability of pattern occuring
        try:
            rate_dict=pickle.load(open(rate_path,"rb"))
            self._rate_data = rate_dict["error_dist"]
            pattern_data = rate_dict["pattern_dist"]
        except Exception as e:
            logger.warning(e)
            exit(1)
        items = sorted(pattern_data.items(),key=lambda x: x[1],reverse=True)
        indexed_probs = [(i,_[1]) for i,_ in enumerate(items)] #distribution will generate indices to patterns
        self._pattern_map = [_[0] for _ in items] #map indices back to patterns
        self._pattern_gen = scipy.stats.rv_discrete(name="pattern",values=(list(zip(*indexed_probs))[0],list(zip(*indexed_probs))[1]))
            
    #go through each nucleotide in each strand and apply a flat fault rate
    def injection_sites(self):
        injection_sites={}
        #go through the strands and pick faults
        for strand_index,strand in enumerate(self._input_library):
            injection_sites[strand_index]=[]
            skip_to_index=-1
            for nuc_index,nuc in enumerate(strand.dna_strand):
                if nuc_index<skip_to_index: continue #allow skipping to get burst errors
                inject_fault=generate.rand()
                if inject_fault<=self._rate_data[nuc_index]:
                    pattern_index = self._pattern_gen.rvs(size=1)
                    pattern = self._pattern_map[int(pattern_index)]
                    inject_index=nuc_index
                    for e in pattern:
                        if e=="I":
                            injection_sites[strand_index].append((inject_index,"2"))
                        elif e=="D":
                            if inject_index>=len(strand.dna_strand): break
                            injection_sites[strand_index].append((inject_index,"1"))
                            inject_index+=1
                        elif e=="R":
                            if inject_index>=len(strand.dna_strand): break
                            injection_sites[strand_index].append((inject_index,"0"))
                            inject_index+=1
                    skip_to_index=inject_index+1 #if we replaced or deleted in this burst error, skip to next position that can take an error
        return injection_sites



class strand_fault_compressed(BaseFI):
    def __init__(self,**args):
        BaseFI.__init__(self)
        self.faulty=args["faulty"]
        self.error_run=args["error_run"]
        self.fails=args["fails"]
        
    def Run(self):
        self._injection_sites={}
        self._strands_after_errors=[]
        self._injection_sites=self.injection_sites_tuple(self._input_library)
        self._strands_after_errors=self.inject_errors_tuple(self._injection_sites,self._input_library)
        return self._strands_after_errors

    #This function converts an index generated in the injection_sites_tuple funciton to an index into the actual library
    def convert_index(self,index_to_convert,library,last_index,sum_to_last_index):
        index=last_index
        _sum=sum_to_last_index
        while index_to_convert>_sum:
            index+=1
            _sum+=library[index][1]
        assert index>=0 and index<len(library)
        return index,_sum
    
    #get strands to inject faults at and nucleotides within the strand
    def injection_sites_tuple(self,input_library):
        converted_index=0
        sum_to_last_index=input_library[0][1]
        #list of strand indexes to chose from
        sum_values=0
        for s in input_library:
            sum_values+=s[1]
        #print sum_values
        strand_indexes=range(sum_values)
        fault_list={} 
        
        strand_locations_before_conversion=sorted(random.sample(strand_indexes,self.faulty))
        strand_locations=[]
        #convert the locations down to the original library to see what strand we are actually injecting 
        for samples in strand_locations_before_conversion:
            converted_index,sum_to_last_index=self.convert_index(samples,input_library,converted_index,sum_to_last_index)
            strand_locations.append(converted_index)
        #print "Strand locations length {}".format(len(strand_locations))
        for index,strand_index in enumerate(strand_locations):
            strand_ID=str(strand_index)+'-'+str(strand_locations_before_conversion[index])
            #print "Strand Id {}".format(strand_ID)
            fault_list[strand_ID]={}
            if self.error_run is True:
                start_point=random.randint(0,len(input_library[strand_index][0])-self.fails)
                nucleotide_indexes=range(start_point,start_point+self.fails)
            else:
                nucleotide_indexes=random.sample(range(0,len(input_library[strand_index][0])),self.fails)
            for nuc_ind in nucleotide_indexes:
                fault_type=generate.rand_in_range(0,2)
                fault_list[strand_ID][nuc_ind]=str(fault_type)
        #print "fault_list length {}".format(len(fault_list))
        return fault_list
       #inject the errors into the selected strands and nucleotides
    def inject_errors_tuple(self,inject_sites,input_library):
        out_list=input_library[:]
        #print "Length of injection sites {}".format(len(inject_sites))
        for index in inject_sites:
            strand_indexes=int(index.split('-')[0])
            error_strand=out_list[strand_indexes][0]
            for fault_indexes in sorted(inject_sites[index],reverse=True):
                #substitution error
                if inject_sites[index][fault_indexes] == '0':
                    #chose a random nucleotide that is different from the current one
                    sub_nucleotide=random.choice(sub_dict[out_list[strand_indexes][0][fault_indexes]])
                    #add on some extra information to the injection sites that indicates the nucleotide used for substitution
                    inject_sites[index][fault_indexes]='0-'+sub_nucleotide
                    #calculate the error strand
                    error_strand=error_strand[0:fault_indexes]+sub_nucleotide+error_strand[fault_indexes+1:len(error_strand)]
                #deletion error
                elif inject_sites[index][fault_indexes] == '1':
                    #add on some extra information to the injection sites, append the nucleotide that was removed from the original strand
                    inject_sites[index][fault_indexes]='1-'+out_list[strand_indexes][0][fault_indexes]
                    error_strand=error_strand[0:fault_indexes]+error_strand[fault_indexes+1:len(error_strand)]
                #insertion error
                elif inject_sites[index][fault_indexes] == '2':
                    insert_nucleotide=random.choice(nuc_list)
                    inject_sites[index][fault_indexes]='2-'+insert_nucleotide
                    error_strand=error_strand[0:fault_indexes]+insert_nucleotide+error_strand[fault_indexes:len(error_strand)]
            #subtract a clean strand from the list location indicated by strand_indexes
            out_list[strand_indexes]=(out_list[strand_indexes][0],out_list[strand_indexes][1]-1)
            #append the error strand to the output list
            out_list.append((error_strand,1))
        return out_list        




from julia.api import Julia
jl = Julia(compiled_modules=False)    
from julia import Main    

#Python interface to DNArSim that is implemented in Julia, the injection module directly calls the DNArSim fault injector to generate nanopore-based error profiles
class DNArSim(BaseFI):
    def __init__(self,**args):
        assert os.environ['DNArSimPath']
        self.DNArSimPath=os.environ['DNArSimPath']
        BaseFI.__init__(self)
        self.kmer_length = args.get("kmer_length",6)
        if "probability_path" not in args:
            raise ValueError("Path to probability path for DNArSim does not exist, please specify")
        #set up julia environment, this should really need to be done once each entire batch run
        Main.eval("""using DelimitedFiles""")
        Main.include(os.path.join(self.DNArSimPath,"functions.jl"))
        Main.include(os.path.join(self.DNArSimPath,"channel.jl"))
        Main.include(os.path.join(self.DNArSimPath,"interface.jl"))
        Main.load_parameters(self.kmer_length,args["probability_path"])
        Main.include(os.path.join(self.DNArSimPath,"loadProb.jl"))
    def julia_run_injection(self):
        inject_set=[x.dna_strand for x in self._input_library]
        out_list=Main.channel(self.kmer_length,inject_set)
        assert len(inject_set)==len(out_list)
        out_list=[FaultDNA(x,y) for x,y in zip(self._input_library,out_list)]
        Main.GC.gc()
        return out_list
    def Run(self):
        return self.julia_run_injection()
