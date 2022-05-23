'''
Fi_Env class ecapsulates all of the pieces needed for fault injection 
'''
from readdist import *
from fault_injector import *
from dnastorage.strand_representation import *
import copy




def fault_injection_modes():
    return ["fixed_rate","missing_strands","strand_fault_compressed","strand_fault","distribution_rate","combo"]

def distribution_functions():
    return ["negative_binomial","poisson"]


class Fi_Env(object):
    def __init__(self,distribution,fault_injection_mode,clean_strands,**kwargs):
        self._read_distributor=ReadDistribution.open(distribution,**kwargs)
        self._fault_mode=BaseFI.open(fault_injection_mode,**kwargs)
        self._clean_strands=clean_strands #clean strand library
        self._og_strands=clean_strands
        self._fault_strands=[]

    def _distribute_reads(self):
        read_count=[]
        new_pool=[]
        dist=[]
        #Tuple format accelerates simulation of fault models of nucleotide strand faults and missing strand faults
        if not isinstance(self._fault_mode,strand_fault_compressed):
            for ind,s in enumerate(self._og_strands):
                read_count=self._read_distributor.gen()
                dist.append((ind,read_count))
                new_pool+=[copy.copy(s) for i in range(0,read_count)]
            return new_pool,dist
        else:
            #make a simple array of tuples in format (strand,count for that strand)
            new_pool=[(copy.copy(s),self._read_distributor.gen()) for s in self._og_strands]
            pool_size=0
            for strand in new_pool:
                pool_size+=strand[1]
            return new_pool,pool_size

    def Run(self):
        #run an instance of fault injection and simulation
        pool,size = self._distribute_reads()
        self._fault_mode.set_library(pool)
        self._fault_mode.Run()
        self._fault_strands=pool #pool should be modified with FaultDNA types
        #make sure everything goes in as a fault strand
        for s_index,s in enumerate(self._fault_strands):
            if not isinstance(s,FaultDNA):
                self._fault_strands[s_index]=FaultDNA(self._fault_strands[s_index],self._fault_strands[s_index].dna_strand) #simple copy with no errors inejcted

    def __next__(self):
        if self._iter_count<self._num_strands:
            res=self._fault_strands[self._iter_count]
            self._iter_count+=1
            return res
        else:
            raise StopIteration
        
    def __iter__(self):
        self._iter_count=0
        self._num_strands=len(self._fault_strands)
        return self 

    def get_strands(self):
        return self._fault_strands
    
if __name__=="__main__":
    import itertools
    import random
    import dnastorage.strand_representation
    #make some fake strands and inject faults just to test data flows of infrastructures

    dna_list=["A","G","C","T"]
    strands=[]
    for i in range(2000):
        strand=[]
        for j in range(200):
            strand.append(dna_list[random.randint(0,3)])            
        strands.append(BaseDNA(dna_strand="".join(strand)))
        print("Base Strand {} {}".format(i,strands[-1].dna_strand))

    #isntantiate fault injection environment
    fault_environment = Fi_Env("poisson","fixed_rate",strands,fault_rate=0.1,mean=10)
    fault_environment.Run()
    for i in fault_environment.get_strands():
        print("Error Strand {}".format(i.dna_strand))
    print(len(fault_environment.get_strands()))

    
