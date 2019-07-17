
#Classes for each of the fault models
import os
import random
import csv
import generate
import time
#structure to bundle together the parser's results, will allow the classes in fault_injector.py be more flexible
class fault_injector_arguments:
    o = None
    input_file = None
    fault_file = None
    fault_model = None
    fault_rate=None
    run = False
    missing = 0
    faulty = 0
    fails = 0
    p1 = 0
    p2 = 0
    

#substitution dictionary used for substitution errors
sub_dict={'A':['G','C','T'],'G':['C','A','T'], 'T':['G','C','A'], 'C':['G','T','A']}
#nucleotide list for insertion errors
nuc_list=['A','C','T','G']


#Base class with some common functions to all of the fault models
class BaseModel:
    def __init__(self,args):
        self._args=args
    #read input file and store it into an array
    def read_file(self):
        if self._args.input_file is not None:
            _file=open(self._args.input_file,'r')
            self._input_library=[]
            for line in _file:
                line=line.strip('\n')
                self._input_library.append(line)
            _file.close()

    #These two setting functions allow easier altertion of the input library to fault injection, and parameters around fault injection
    def set_library(self,input_strands):
        self._input_library=input_strands
        
    def set_faulty_missing(self,desired_faulty_count,desired_missing_count):
        self._args.faulty=desired_faulty_count
        self._args.missing=desired_missing_count

    #NOTE:deprecated function, random.sample is used for generating indexes now
    def random_array_element(self,input_set):
        random_index=random.randint(0,len(input_set)-1)
        return random_index,input_set[random_index]

    #write out the new DNA file with the faults
    def write_out(self,faulty_strands):
        out_file=open(self._args.o,'w+')
        for strand in faulty_strands:
            out_file.write("{}\n".format(strand))
        out_file.close()
    def read_csv(self,file_name):
        _file=open(file_name,'r')
        csv_parsed=csv.reader(_file,delimiter=',')
        return csv_parsed
    def set_faults_per_strand(self,desired_faults_per_strand):
        self._args.fails=desired_faults_per_strand
    def set_fault_rate(self,desired_fault_rate):
        self._args.fault_rate=desired_fault_rate

#this class applies a fixed error rate to each nucleotide
#Strands are a list, not a tuple representation
class fixed_rate(BaseModel):
    def __init__(self,args):
        BaseModel.__init__(self,args)
        self.read_file()
    def Run(self,write_out=False):
        self._injection={}
        self._error_strands=[]
        self._injection=self.injection_sites(self._input_library)
        self._error_strands=self.inject_faults(self._input_library,self._injection)
        if write_out is False:
            return self._error_strands
        else:
            self.write_out(self._error_strands)
            return self._error_strands

    #go through each nucleotide in each strand and apply a flat fault rate
    def injection_sites(self,input_library):
        generate.seed()
        injection_sites={}
        #go through the strands and pick faults
        for strand_index,strand in enumerate(input_library):
            injection_sites[strand_index]={}
            for nuc_index,nuc in enumerate(strand):
                inject_fault=generate.rand()
                if inject_fault<=self._args.fault_rate:
                    #need to inject a fault at this nucleotide, fault types are equally probable
                    #FIX if want to have variable error rates per type
                    injection_sites[strand_index][nuc_index]=str(generate.rand_in_range(0,2))
        return injection_sites
                    
    #inject errors in the list of strands, each 
    def inject_faults(self,input_library,inject_sites):
        out_list=input_library[:]
        #print inject_sites
        for strand_indexes in inject_sites:
            for fault_indexes in sorted(inject_sites[strand_indexes],reverse=True):
                #substitution error
                if inject_sites[strand_indexes][fault_indexes] == '0':
                    #chose a random nucleotide that is different from the current one
                    sub_nucleotide=random.choice(sub_dict[out_list[strand_indexes][fault_indexes]])
                    #add on some extra information to the injection sites that indicates the nucleotide used for substitution
                    inject_sites[strand_indexes][fault_indexes]='0-'+sub_nucleotide
                    out_list[strand_indexes]= out_list[strand_indexes][0:fault_indexes]+sub_nucleotide+out_list[strand_indexes][fault_indexes+1:len(out_list[strand_indexes])]
                #deletion error
                elif inject_sites[strand_indexes][fault_indexes] == '1':
                     #add on some extra information to the injection sites, append the nucleotide that was removed from the original strand 
                    inject_sites[strand_indexes][fault_indexes]='1-'+out_list[strand_indexes][fault_indexes]
                    out_list[strand_indexes]=out_list[strand_indexes][0:fault_indexes]+out_list[strand_indexes][fault_indexes+1:len(out_list[strand_indexes])]
                #insertion error
                elif inject_sites[strand_indexes][fault_indexes] == '2':
                    insert_nucleotide=random.choice(nuc_list)
                    inject_sites[strand_indexes][fault_indexes]='2-'+insert_nucleotide
                    out_list[strand_indexes]=out_list[strand_indexes][0:fault_indexes]+insert_nucleotide+out_list[strand_indexes][fault_indexes:len(out_list[strand_indexes])]
        return out_list

#class that contains functionality for missing strands fault model
class miss_strand(BaseModel):
    def __init__(self,args):
        BaseModel.__init__(self,args)
        self.read_file()
        
    def Run(self,write_out=False):
        removal_sites=[]
        self._faulty_strands=[]
        removal_sites=self.remove_sites(self._input_library)
        self._faulty_strands=self.remove_strands(self._removal_sites,self._input_library)
        if write_out is False:
            return self._faulty_strands
        else:
            self.write_out(self._faulty_strands)    
    #gets a collection of indexes to remove from the input file, all indexes are relative to the original file e.g removal_index=1 means we will remove the second strand (indexing starting at 0)
    def remove_sites(self,input_library):
        return random.sample(range(len(input_library)),self._args.missing)

 #remove the strands indicated by the sites list from the input library
    def remove_strands(self,sites,input_library):
        out_list=input_library[:]
        #remove in reversye order so that we do not have to change the index of some sites after deleting
        for strand_index in sorted(sites,reverse=True):
            del out_list[strand_index]
        return out_list
            


        
#class that contains fucntionality for injecting faults into strands, only injects faults into the data fields, not the primers
class strand_fault(BaseModel):
    def __init__(self,args):
        BaseModel.__init__(self,args)
        self.read_file()
        #read csv file if it is available
        if self._args.fault_file is not None:
            self._csv_data=self.read_csv(self._args.fault_file)
    def Run(self,write_out=False):
        self._injection={}
        self._error_strands=[]
        #if there is no csv file, chose random spots and errors
        if self._args.fault_file is None:
            self._injection=self.injection_sites(self._input_library)
            self._error_strands=self.inject_errors(self._injection,self._input_library)
        #if there is a csv file, inject errors based off the data
        elif self._args.fault_file is not None:
            self._error_strands=self.inject_distribution(self._input_library,self._csv_data)
        if write_out is False:
            return self._error_strands
        else:
            self.write_out(self._error_strands)


    #this function complies with a tuple format for the library
    def Run_tuple(self):
        self._injection_sites={}
        self._strands_after_errors=[]
        self._injection_sites=self.injection_sites_tuple(self._input_library)
        self._strands_after_errors=self.inject_errors_tuple(self._injection_sites,self._input_library)
        return self._strands_after_errors
        


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
        strand_index_start=self._args.p1
        strand_locations_before_conversion=sorted(random.sample(strand_indexes,self._args.faulty))
        strand_locations=[]
        #need to seed the C generator
        generate.seed()
        #convert the locations down to the original library to see what strand we are actually injecting 
        for samples in strand_locations_before_conversion:
            converted_index,sum_to_last_index=self.convert_index(samples,input_library,converted_index,sum_to_last_index)
            strand_locations.append(converted_index)
        #print "Strand locations length {}".format(len(strand_locations))
        for index,strand_index in enumerate(strand_locations):
            strand_ID=str(strand_index)+'-'+str(strand_locations_before_conversion[index])
            #print "Strand Id {}".format(strand_ID)
            fault_list[strand_ID]={}
            if self._args.run is True:
                start_point=random.randint(strand_index_start,len(input_library[strand_index][0])-self._args.p2-self._args.fails)
                nucleotide_indexes=range(start_point,start_point+self._args.fails)
            else:
                nucleotide_indexes=random.sample(range(strand_index_start,len(input_library[strand_index][0])-self._args.p2),self._args.fails)
            for nuc_ind in nucleotide_indexes:
                fault_type=generate.rand_in_range(0,2)
                fault_list[strand_ID][nuc_ind]=str(fault_type)
        #print "fault_list length {}".format(len(fault_list))
        return fault_list

    #This function converts an index generated in the injection_sites_tuple funciton to an index into the actual library
    def convert_index(self,index_to_convert,library,last_index,sum_to_last_index):
        index=last_index
        _sum=sum_to_last_index
        while index_to_convert>_sum:
            index+=1
            _sum+=library[index][1]
        assert index>=0 and index<len(library)
        return index,_sum
        
    
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
        strand_index_start=self._args.p1
        #go through the prob_data data structure and grab relevant data and scale to make it easier for random number generate
        for row_of_data in prob_data:
            #each row is a row from the input spreadsheet
            if row_of_data[0] == "Overall Error":
                overall_error=[int(float(prob)*10000) for prob in row_of_data[1+strand_index_start:len(row_of_data[1:])-self._args.p2+1]]
                self._fault_rate=overall_error
                
            elif row_of_data[0] == "Del/Error":
                del_giv_error=[int(float(prob)*1000) for prob in row_of_data[1+strand_index_start:len(row_of_data[1:])-self._args.p2+1]]
                self._del_rate=del_giv_error
            elif row_of_data[0] == "Ins/Error":
                inser_giv_error=[int(float(prob)*1000) for prob in row_of_data[1+strand_index_start:len(row_of_data[1:])-self._args.p2+1]]
                self._ins_rate = inser_giv_error
            elif row_of_data[0] == "Sub/Error":
                sub_giv_error=[int(float(prob)*1000) for prob in row_of_data[1+strand_index_start:len(row_of_data[1:])-self._args.p2+1]]
                self._sub_rate = sub_giv_error
                
        assert(len(overall_error)>0 and len(del_giv_error)>0 and len(inser_giv_error)>0 and len(sub_giv_error)>0)
        assert(len(overall_error)==(len(out_list[0])-self._args.p1-self._args.p2))
        #go through each strand and nucleotide and inject errors
        for strand_index, strand in enumerate(out_list):
            for nuc_index, nucleotide in enumerate(strand[strand_index_start:len(strand)-self._args.p2]):
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




        
    #get strands to inject faults at and nucleotides within the strand
    def injection_sites(self,input_library):
        #list of strand indexes to chose from 
        strand_indexes=range(len(input_library))
        fault_list={} 
        strand_index_start=self._args.p1
        strand_locations=random.sample(range(len(input_library)),self._args.faulty)
        #need to seed the C generator
        generate.seed()
        for strand_index in strand_locations:
            fault_list[strand_index]={}
            if self._args.run is True:
                start_point=random.randint(strand_index_start,len(input_library[strand_index])-self._args.p2-self._args.fails)
                nucleotide_indexes=range(start_point,start_point+self._args.fails)
            else:
                nucleotide_indexes=random.sample(range(strand_index_start,len(input_library[strand_index])-self._args.p2),self._args.fails)
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
                if inject_sites[strand_indexes][fault_indexes] == '0':
                    #chose a random nucleotide that is different from the current one
                    sub_nucleotide=random.choice(sub_dict[out_list[strand_indexes][fault_indexes]])
                    #add on some extra information to the injection sites that indicates the nucleotide used for substitution
                    inject_sites[strand_indexes][fault_indexes]='0-'+sub_nucleotide
                    out_list[strand_indexes]= out_list[strand_indexes][0:fault_indexes]+sub_nucleotide+out_list[strand_indexes][fault_indexes+1:len(out_list[strand_indexes])]
                #deletion error
                elif inject_sites[strand_indexes][fault_indexes] == '1':
                     #add on some extra information to the injection sites, append the nucleotide that was removed from the original strand 
                    inject_sites[strand_indexes][fault_indexes]='1-'+out_list[strand_indexes][fault_indexes]
                    out_list[strand_indexes]=out_list[strand_indexes][0:fault_indexes]+out_list[strand_indexes][fault_indexes+1:len(out_list[strand_indexes])]
                #insertion error
                elif inject_sites[strand_indexes][fault_indexes] == '2':
                    insert_nucleotide=random.choice(nuc_list)
                    inject_sites[strand_indexes][fault_indexes]='2-'+insert_nucleotide
                    out_list[strand_indexes]=out_list[strand_indexes][0:fault_indexes]+insert_nucleotide+out_list[strand_indexes][fault_indexes:len(out_list[strand_indexes])]            
        return out_list


#This class implements a combination of missing strands and error strands 
class combo(BaseModel):
    def __init__(self,args):
        BaseModel.__init__(self,args)
        self.strand_faults=strand_fault(args)
        self.missing_strands=miss_strand(args)
        self.read_file()
        #read csv file if it is available
        if self._args.fault_file is not None:
            self._csv_data=self.read_csv(self._args.fault_file)
    def Run(self,write_out=False):
        strands_after_missing=[]
        missing_sites=[]
        missing_strands_with_errors=[]
        error_sites={}
        self.missing_strands.remove_sites(missing_sites,self._input_library)
        strands_after_missing=self.missing_strands.remove_strands(missing_sites,self._input_library)
        if self._args.fault_file is None:
            error_sites=self.strand_faults.injection_sites(strands_after_missing)
            missing_strands_with_errors=self.strand_faults.inject_errors(error_sites,strands_after_missing)
        elif self._args.fault_file is not None:
            missing_strands_with_errors=self.strand_faults.inject_distribution(strands_after_missing,self._csv_data)
        if write_out is False:
            return missing_strands_with_errors
        else:
            self.write_out(missing_strands_with_errors)


                    
#test code for the fault injections 
if __name__ == "__main__":

    #test the missing strands fault injections
    arguments=fault_injector_arguments()
    arguments.input_file="test_dna.txt"
    arguments.missing=10
    arguments.p1=20
    arguments.p2=20
    #instantiate class that will remove strands 
    missing_strands=miss_strand(arguments)
    missing_strands.read_file()
    clean_strands=missing_strands._input_library
    removal_locations=missing_strands.remove_sites(clean_strands)
    final_strands=missing_strands.remove_strands(removal_locations,clean_strands)
    for site in removal_locations:
        if clean_strands[site] in final_strands:
            print "strand {} still in final strands".format(site)
            assert False
    
