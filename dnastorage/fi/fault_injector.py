#Classes for each of the fault models
import os
import random

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


#class that contains functionality for missing strands fault model
class miss_strand(BaseModel):
    def __init__(self,args):
        BaseModel.__init__(self,args)
        self.read_file()
        
    def Run(self):
        self._removal_sites=[]
        self._faulty_strands=[]
        self.remove_sites(self._removal_sites,self._input_library)
        self._faulty_strands=self.remove_strands(self._removal_sites)
        self.write_out(self._faulty_strands)

        
    #gets a collection of indexes to remove from the input file, all indexes are relative to the original file e.g removal_index=1 means we will remove the second strand (indexing starting at 0)
    def remove_sites(self,removal_sites,input_library):
        strand_indexes=range(len(input_library))
        for fault_ID in range(self._args.missing):
            index,removal_strand_index=self.random_array_element(strand_indexes)
            del strand_indexes[index]
            removal_sites.append(removal_strand_index)

 #remove the strands indicated by the sites list from the input library
    def remove_strands(self,sites):
        out_list=self._input_library[:]
        #remove in reverse order so that we do not have to change the index of some sites after deleting
        for strand_index in sorted(sites,reverse=True):
            del out_list[strand_index]
        return out_list
            


        
#class that contains fucntionality for injecting faults into strands 
class strand_fault(BaseModel):
    def __init__(self,args):
        BaseModel.__init__(self,args)
        self.read_file()
    def Run(self):
        self._injection={}
        self._error_strands=[]
        self._injection=self.injection_sites()
        self._error_strands=self.inject_errors(self._injection)
        self.write_out(self._error_strands)
    #get strands to inject faults at and nucleotides within the strand
    def injection_sites(self):
        #list of strand indexes to chose from 
        strand_indexes=range(len(self._input_library))
        fault_list={} 
        strand_index_start=self._args.p1
        for strand_num in range(self._args.faulty):
            #chose a strand
            index,chosen_strand_index=self.random_array_element(strand_indexes)
            #delete chosen strands from list so that it is not chosen again
            del strand_indexes[index]
            strand=self._input_library[chosen_strand_index]
            
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
    def inject_errors(self,inject_sites):
        out_list=self._input_library
        for strand_indexes in inject_sites:
            for fault_indexes in sorted(inject_sites[strand_indexes],reverse=True):
                #substitution error
                if inject_sites[strand_indexes][fault_indexes] == 0:
                    #chose a random nucleotide that is different from the current one
                    sub_nucleotide=random.choice(sub_dict[out_list[strand_indexes][fault_indexes]])
                    out_list[strand_indexes]=out_list[strand_indexes][0:fault_indexes]+sub_nucleotide+out_list[strand_indexes][fault_indexes+1:len(out_list[strand_indexes])]
                #deletion error
                elif inject_sites[strand_indexes][fault_indexes] == 1:
                    out_list[strand_indexes]=out_list[strand_indexes][0:fault_indexes]+out_list[strand_indexes][fault_indexes+1:len(out_list[strand_indexes])]
                #insertion error
                elif inject_sites[strand_indexes][fault_indexes] == 2:
                    insert_nucleotide=random.choice(nuc_list)
                    out_list[strand_indexes]=out_list[strand_indexes][0:fault_indexes]+insert_nucleotide+out_list[strand_indexes][fault_indexes:len(out_list[strand_indexes])]
                    
        return out_list
                



                    
