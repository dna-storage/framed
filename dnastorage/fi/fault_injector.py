#Classes for each of the fault models
import os
import random

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


class miss_strand(BaseModel):
    def __init__(self,args):
        BaseModel.__init__(self,args)
        self.read_file()
    def Run(self):
        self._removal_sites=[]
        self.remove_sites(self._removal_sites)
        self.write_out()
        
    def remove_sites(self,removal_sites):
        index_range=len(self._input_library)
        for fault_ID in range(self._args.missing):
            removal_index=random.randint(0,index_range-1)
            if removal_index not in removal_sites:
                print removal_index
                print "{} {}".format(fault_ID,len(self._input_library[removal_index]))
                removal_sites.append(removal_index)
    
    def write_out(self):
        out_file=open(self._args.o,'w+')
        for index, strand in enumerate(self._input_library):
            if index not in self._removal_sites:
                out_file.write("{}\n".format(strand))
        out_file.close()
