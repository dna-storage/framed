#Classes for each of the fault models
import os
import random
import csv
import pickle
import re

import dnastorage.util.generate as generate
import time
from dnastorage.strand_representation import *
from dnastorage.fi.fault_strand_representation import *
from Bio import SeqIO
from scipy import stats

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
        elif fault_injector=="missing_strands":
            return miss_strand(**kwargs)
        elif fault_injector=="strand_fault_compressed":
            return strand_fault_compressed(**kwargs)
        elif fault_injector=="strand_fault":
            return strand_fault(**kwargs)
        elif fault_injector=="distribution_rate":
            return distribution_rate(**kwargs)
        elif fault_injector=="sequencing_experiment":
            return sequencing_experiment(**kwargs)
        elif fault_injector=="pattern_fixed_rate":
            return pattern_fixed_rate(**kwargs)
        elif fault_injector=="combo":
            return combo(**kwargs)
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
        for record in SeqIO.parse(self._sequencing_data_path,"fastq"):
            seq_strand,subs = re.subn("U","T",str(record.seq))
            index_ints = self._map.get(record.id,None)
            if index_ints is None: continue #dead strands
            lib_strand=input_library_map[index_ints]
            out_list.append(FaultDNA(lib_strand,seq_strand))
        assert len(out_list)>0
        return out_list
    def Run(self):
        return self.sequence_run_injection()
  
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
        strand_rate_path = args.get("error_rate_path",None) #path to rates for each base indicating probability of pattern occuring
        pattern_dist_path = args.get("pattern_dist_path",None) #distribution of error patterns, given an error ocurred
        try:
            self._rate_data = pickle.load(open(strand_rate_path,"rb"))
            pattern_data = pickle.load(open(pattern_dist_path,"rb"))
        except Exception as e:
            logger.warning(e)
            exit(1)
        items = sorted(pattern_data.items(),key=lambda x: x[1],reverse=True)
        indexed_probs = [(i,_[1]) for i,_ in enumerate(items)] #distribution will generate indices to patterns
        self._pattern_map = [_[0] for _ in items] #map indices back to patterns
        self._pattern_gen = stats.rv_discrete(name="pattern",values=(list(zip(*indexed_probs))[0],list(zip(*indexed_probs))[1]))
            
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

   
#class that contains functionality for missing strands fault model
class miss_strand(BaseFI):
    def __init__(self,**args):
        BaseFI.__init__(self)
        self._missing=args["missing"]
        
    def Run(self):
        removal_sites=[]
        self._faulty_strands=[]
        removal_sites=self.remove_sites(self._input_library)
        self._faulty_strands=self.remove_strands(self._removal_sites,self._input_library)
        return self._faulty_strands
      
    #gets a collection of indexes to remove from the input file, all indexes are relative to the original file e.g removal_index=1 means we will remove the second strand (indexing starting at 0)
    def remove_sites(self,input_library):
        return random.sample(range(len(input_library)),self.missing)

 #remove the strands indicated by the sites list from the input library
    def remove_strands(self,sites,input_library):
        out_list=input_library[:]
        #remove in reversye order so that we do not have to change the index of some sites after deleting
        for strand_index in sorted(sites,reverse=True):
            del out_list[strand_index]
        return out_list
            



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

    
    

class distribution_rate(BaseFI):
    def __init__(self,**args):
        BaseFI.__init__(self)
        self.fault_file=args["fault_file"]
        self._csv_data=self.read_csv(self.fault_file)

    def Run(self,write_out=False):
        self._injection={}
        self._error_strands=[]
        self._error_strands=self.inject_distribution(self._input_library,self._csv_data)
        return self._error_strands
        
    #These functions are used for testing code
    def get_fault_spread(self):
        return self._fault_spread
    def get_del_spread(self):
        return self._del_spread
    def get_ins_spread(self):
        return self._ins_spread
    def get_sub_spread(self):
        return self._sub_spread
    def get_fault_rate(self):
        return self._fault_rate
    def get_del_rate(self):
        return self._del_rate
    def get_ins_rate(self):
        return self._ins_rate
    def get_sub_rate(self):
        return self._sub_rate

    #function to insert errors based on the distribution data in the csv files 
    def inject_distribution(self,input_strands,prob_data):
        #variables for collecting data for testing 
        self._fault_spread={}
        self._del_spread={}
        self._ins_spread={}
        self._sub_spread={}
        self._fault_rate=[]
        self._del_rate=[]
        self._ins_rate=[]
        self._sub_rate=[]
        time0=time.time()
        out_list=input_strands[:]
        #inject faults throughout the input strands, no subset is chosen, unlike the other fault model, maybe add that in later if wanted?
        overall_error=[]
        del_giv_error=[]
        inser_giv_error=[]
        sub_giv_error=[]
        #go through the prob_data data structure and grab relevant data and scale to make it easier for random number generate
        for row_of_data in prob_data:
            #each row is a row from the input spreadsheet
            if row_of_data[0] == "Overall Error":
                overall_error=[int(float(prob)*10000) for prob in row_of_data[1:len(row_of_data[1:])+1]]
                self._fault_rate=overall_error
                
            elif row_of_data[0] == "Del/Error":
                del_giv_error=[int(float(prob)*1000) for prob in row_of_data[1:len(row_of_data[1:])+1]]
                self._del_rate=del_giv_error
            elif row_of_data[0] == "Ins/Error":
                inser_giv_error=[int(float(prob)*1000) for prob in row_of_data[1:len(row_of_data[1:])+1]]
                self._ins_rate = inser_giv_error
            elif row_of_data[0] == "Sub/Error":
                sub_giv_error=[int(float(prob)*1000) for prob in row_of_data[1:len(row_of_data[1:])+1]]
                self._sub_rate = sub_giv_error
                
        assert(len(overall_error)>0 and len(del_giv_error)>0 and len(inser_giv_error)>0 and len(sub_giv_error)>0)
        #go through each strand and nucleotide and inject errors
        for strand_index, strand in enumerate(out_list):
            for nuc_index, nucleotide in enumerate(strand[0:len(strand)]):
                #generate random numbers and apply a appropriate error if error has occured
                inject_error=random.randint(1,10000)
                #should inject error 
                if inject_error <= overall_error[nuc_index]:
                    #collect the amount of errors injected 
                    if nuc_index not in self._fault_spread:
                        self._fault_spread[nuc_index]=0
                    else:
                         self._fault_spread[nuc_index]=self._fault_spread[nuc_index]+1
                    del_boundary_lower=1
                    del_boundary_upper=del_giv_error[nuc_index]
                    inser_boundary_lower=del_boundary_upper+1
                    inser_boundary_upper=del_boundary_upper+inser_giv_error[nuc_index]
                    sub_boundary_lower=inser_boundary_upper+1
                    sub_boundary_upper=inser_boundary_upper+sub_giv_error[nuc_index]
                    
                    #random value to select amongst the error type
                    error_type=random.randint(1,sub_boundary_upper)
                    
                    #inject deletion error
                    if error_type >= del_boundary_lower and error_type <= del_boundary_upper:
                        #collect the amount of deletion errors injected 
                        if nuc_index not in self._del_spread:
                            self._del_spread[nuc_index]=0
                        else:
                            self._del_spread[nuc_index]=self._del_spread[nuc_index]+1
                        out_list[strand_index]=out_list[strand_index][0:nuc_index]+out_list[strand_index][nuc_index+1:len(out_list[strand_index])]
                        
                    #inject insertion error
                    elif error_type >= inser_boundary_lower and error_type <= inser_boundary_upper:
                        #collect the amount of insertion errors injected 
                        if nuc_index not in self._ins_spread:
                            self._ins_spread[nuc_index]=0
                        else:
                            self._ins_spread[nuc_index]=self._ins_spread[nuc_index]+1           
                        insert_nucleotide=random.choice(nuc_list)
                        out_list[strand_index]=out_list[strand_index][0:nuc_index]+insert_nucleotide+out_list[strand_index][nuc_index:len(out_list[strand_index])]
                    #substitution error
                    elif error_type >= sub_boundary_lower and error_type <= sub_boundary_upper:
                        #collect the amount of substitution errors injected 
                        if nuc_index not in self._sub_spread:
                            self._sub_spread[nuc_index]=0
                        else:
                            self._sub_spread[nuc_index]=self._sub_spread[nuc_index]+1
                        sub_nucleotide=random.choice(sub_dict[out_list[strand_index][nuc_index]])
                        out_list[strand_index]= out_list[strand_index][0:nuc_index]+sub_nucleotide+out_list[strand_index][nuc_index+1:len(out_list[strand_index])]
                    else:
                        assert(0)
        return out_list

    

class strand_fault(BaseFI):
    def __init__(self,**args):
        BaseFI.__init__(self)
        self.faulty=args["faulty"]
        self.fails=args["fails"]
        self.error_run=args["error_run"]
        
    def Run(self,write_out=False):
        self._injection={}
        self._error_strands=[]
        #if there is no csv file, chose random spots and errors
        self._injection=self.injection_sites(self._input_library)
        self._error_strands=self.inject_errors(self._injection,self._input_library)
        return self._error_strands
       
    #get strands to inject faults at and nucleotides within the strand
    def injection_sites(self,input_library):
        #list of strand indexes to chose from 
        strand_indexes=range(len(input_library))
        fault_list={} 

        print ("injection sites:",len(input_library))
        
        strand_locations=random.sample(range(len(input_library)),self.faulty)
        for strand_index in strand_locations:
            fault_list[strand_index]={}
            if self.run is True:
                start_point=random.randint(0,len(input_library[strand_index])-self.fails)
                nucleotide_indexes=range(start_point,start_point+self.fails)
            else:
                nucleotide_indexes=random.sample(range(0,len(input_library[strand_index])),self.fails)
            for nuc_ind in nucleotide_indexes:
                fault_type=generate.rand_in_range(0,2)
                fault_list[strand_index][nuc_ind]=str(fault_type)
                
        return fault_list
 
    #inject the errors into the selected strands and nucleotides
    def inject_errors(self,inject_sites,input_library):
        #time0=time.time()
        out_list=input_library[:]
        #time1=time.time()
        #print "copy time {}".format(time1-time0)
        for strand_indexes in inject_sites:
            for fault_indexes in sorted(inject_sites[strand_indexes],reverse=True):
                #substitution error
                library_strand=self._input_library[strand_indexes].dna_strand
                new_strand=""
                if inject_sites[strand_indexes][fault_indexes] == '0':
                    #chose a random nucleotide that is different from the current one
                    sub_nucleotide=random.choice(sub_dict[out_list[strand_indexes][fault_indexes]])
                    #add on some extra information to the injection sites that indicates the nucleotide used for substitution
                    inject_sites[strand_indexes][fault_indexes]='0-'+sub_nucleotide
                    new_strand=library_strand[0:fault_indexes]+sub_nucleotide+library_strand[fault_indexes+1:len(library_strand)]

                #deletion error
                elif inject_sites[strand_indexes][fault_indexes] == '1':
                     #add on some extra information to the injection sites, append the nucleotide that was removed from the original strand 
                    inject_sites[strand_indexes][fault_indexes]='1-'+out_list[strand_indexes][fault_indexes]
                    new_strand=library_strand[0:fault_indexes]+library_strand[fault_indexes+1:len(library_strand)]
                #insertion error
                elif inject_sites[strand_indexes][fault_indexes] == '2':
                    insert_nucleotide=random.choice(nuc_list)
                    inject_sites[strand_indexes][fault_indexes]='2-'+insert_nucleotide
                    new_strand=library_strand[0:fault_indexes]+insert_nucleotide+library_strand[fault_indexes:len(library_strand)]

                else:
                    raise ValueError()
                out_list[strand_indexes]=FaultDNA(self._input_library[strand_indexes],new_strand)
        return out_list


#This class implements a combination of missing strands and error strands 
class combo(BaseFI):
    def __init__(self,**args):
        BaseFI.__init__(self)
        self.strand_faults=strand_fault(**args)
        self.missing_strands=miss_strand(**args)
        
    def Run(self,write_out=False):
        strands_after_missing=[]
        missing_sites=[]
        missing_strands_with_errors=[]
        error_sites={}
        self.missing_strands.remove_sites(missing_sites,self._input_library)
        strands_after_missing=self.missing_strands.remove_strands(missing_sites,self._input_library)
        missing_strands_with_errors=self.strand_faults.inject_distribution(strands_after_missing,self._csv_data)
        return missing_strands_with_errors


