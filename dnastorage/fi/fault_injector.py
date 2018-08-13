#Classes for each of the fault models
import os
import random
import csv

#substitution dictionary used for substitution errors
sub_dict={'A':['G','C','T'],'G':['C','A','T'], 'T':['G','C','A'], 'C':['G','T','A']}
#nucleotide list for insertion errors
nuc_list=['A','C','T','G']


#Base class with some common functions to all of the fault models
class BaseModel:
    def __init__(self,args):
        self._args=args
    #read input file
    def read_file(self):
        _file=open(self._args.input_file[0],'r')
        self._input_library=[]
        
        for line in _file:
            line=line.strip('\n')
            self._input_library.append(line)
        _file.close()
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

#class that contains functionality for missing strands fault model
class miss_strand(BaseModel):
    def __init__(self,args):
        BaseModel.__init__(self,args)
        self.read_file()
        
    def Run(self):
        self._removal_sites=[]
        self._faulty_strands=[]
        self.remove_sites(self._removal_sites,self._input_library)
        self._faulty_strands=self.remove_strands(self._removal_sites,self._input_library)
        self.write_out(self._faulty_strands)

        
    #gets a collection of indexes to remove from the input file, all indexes are relative to the original file e.g removal_index=1 means we will remove the second strand (indexing starting at 0)
    def remove_sites(self,removal_sites,input_library):
        strand_indexes=range(len(input_library))
        for fault_ID in range(self._args.missing):
            index,removal_strand_index=self.random_array_element(strand_indexes)
            del strand_indexes[index]
            removal_sites.append(removal_strand_index)

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
    def Run(self):
        self._injection={}
        self._error_strands=[]
        #if there is no csv file, chose random spots and errors
        if self._args.fault_file is None:
            self._injection=self.injection_sites(self._input_library)
            self._error_strands=self.inject_errors(self._injection,self._input_library)
        #if there is a csv file, inject errors based off the data
        elif self._args.fault_file is not None:
            self._error_strands=self.inject_distribution(self._input_library,self._csv_data)
        self.write_out(self._error_strands)



    #function to insert errors based on the distribution data in the csv files 
    def inject_distribution(self,input_strands,prob_data):
        out_list=input_strands[:]
        #inject faults throughout the input strands, no subset is chosen, unlike the other fault model, maybe add that in later if wanted?
        overall_error=[]
        del_giv_error=[]
        inser_giv_error=[]
        sub_giv_error=[]
        strand_index_start=self._args.p1
        
        #go through the prob_data data structure and grab relevant data and condition to make it easier for random number generate
        for row_of_data in prob_data:
            #each row is a row from the input spreadsheet
            if row_of_data[0] == "Overall Error":
                overall_error=[int(float(prob)*10000) for prob in row_of_data[1+strand_index_start:len(row_of_data[1:])-self._args.p2+1]]
            elif row_of_data[0] == "Del/Error":
                del_giv_error=[int(float(prob)*1000) for prob in row_of_data[1+strand_index_start:len(row_of_data[1:])-self._args.p2+1]]
            elif row_of_data[0] == "Ins/Error":
                inser_giv_error=[int(float(prob)*1000) for prob in row_of_data[1+strand_index_start:len(row_of_data[1:])-self._args.p2+1]]
            elif row_of_data[0] == "Sub/Error":
                sub_giv_error=[int(float(prob)*1000) for prob in row_of_data[1+strand_index_start:len(row_of_data[1:])-self._args.p2+1]]
        assert(len(overall_error)>0 and len(del_giv_error)>0 and len(inser_giv_error)>0 and len(sub_giv_error)>0)
        assert(len(overall_error)==(len(out_list[0])-self._args.p1-self._args.p2))
        #go through each strand and nucleotide and inject errors
        for strand_index, strand in enumerate(out_list):
            for nuc_index, nucleotide in enumerate(strand[strand_index_start:len(strand)-self._args.p2]):
                #generate random numbers and apply a appropriate error if error has occured
                inject_error=random.randint(1,10000)
                #should inject error 
                if inject_error <= overall_error[nuc_index]:
                    print "error"
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
                        print "deletion"
                        out_list[strand_index]=out_list[strand_index][0:nuc_index]+out_list[strand_index][nuc_index+1:len(out_list[strand_index])]
                    #inject insertion error
                    elif error_type >= inser_boundary_lower and error_type <= inser_boundary_upper:
                        print "insertion" 
                        insert_nucleotide=random.choice(nuc_list)
                        out_list[strand_index]=out_list[strand_index][0:nuc_index]+insert_nucleotide+out_list[strand_index][nuc_index:len(out_list[strand_index])]
                    #substitution error
                    elif error_type >= sub_boundary_lower and error_type <= sub_boundary_upper:
                        print "sub"
                        sub_nucleotide=random.choice(sub_dict[out_list[strand_index][nuc_index]])
                        out_list[strand_index]= out_list[strand_index][0:nuc_index]+sub_nucleotide+out_list[strand_index][nuc_index+1:len(out_list[strand_index])]
                    else:
                        print del_boundary_upper
                        print inser_boundary_lower
                        print inser_boundary_upper
                        print sub_boundary_lower
                        print sub_boundary_upper
                        print error_type
                        #should not come here
                        assert(0)
        return out_list




        
    #get strands to inject faults at and nucleotides within the strand
    def injection_sites(self,input_library):
        #list of strand indexes to chose from 
        strand_indexes=range(len(input_library))
        fault_list={} 
        strand_index_start=self._args.p1
        for strand_num in range(self._args.faulty):
            #chose a strand
            index,chosen_strand_index=self.random_array_element(strand_indexes)
            #delete chosen strands from list so that it is not chosen again
            del strand_indexes[index]
            strand=input_library[chosen_strand_index]
            
            fault_list[chosen_strand_index]={}
            #create list of nucleotide indexes so that we can avoid choosing the same index
            nucleotide_indexes=range(strand_index_start,len(strand)-self._args.p2)
            if self._args.run is True:
                start_point=random.randint(strand_index_start,len(strand)-self._args.p2-self._args.fails)
                nucleotide_indexes=range(start_point,start_point+self._args.fails)
            #now find some nucleotides within the strand to inject faults
            for nucleotide_fault in range(self._args.fails):
                if self._args.run is not True:
                    index,nucleotide_index=self.random_array_element(nucleotide_indexes)
                    del nucleotide_indexes[index]
                elif self._args.run is True:
                    nucleotide_index=nucleotide_indexes[nucleotide_fault]
                fault_type=random.randint(0,2)
                fault_list[chosen_strand_index][nucleotide_index]=fault_type
        print fault_list
        return fault_list

     

        
    #inject the errors into the selected strands and nucleotides
    def inject_errors(self,inject_sites,input_library):
        out_list=input_library[:]
        for strand_indexes in inject_sites:
            for fault_indexes in sorted(inject_sites[strand_indexes],reverse=True):
                #substitution error
                if inject_sites[strand_indexes][fault_indexes] == 0:
                    #chose a random nucleotide that is different from the current one
                    sub_nucleotide=random.choice(sub_dict[out_list[strand_indexes][fault_indexes]])
                    out_list[strand_indexes]= out_list[strand_indexes][0:fault_indexes]+sub_nucleotide+out_list[strand_indexes][fault_indexes+1:len(out_list[strand_indexes])]
                #deletion error
                elif inject_sites[strand_indexes][fault_indexes] == 1:
                    out_list[strand_indexes]=out_list[strand_indexes][0:fault_indexes]+out_list[strand_indexes][fault_indexes+1:len(out_list[strand_indexes])]
                #insertion error
                elif inject_sites[strand_indexes][fault_indexes] == 2:
                    insert_nucleotide=random.choice(nuc_list)
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
    def Run(self):
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
        self.write_out(missing_strands_with_errors)


                    
