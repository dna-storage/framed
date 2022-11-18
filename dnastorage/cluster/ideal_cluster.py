import numpy as np
import Levenshtein as ld
from dnastorage.cluster.basecluster import *

'''
Simple ideal cluster that uses fault injector information
'''
class IdealCluster(BaseCluster):
    def __init__(self):
        BaseCluster.__init__(self)        
    def Run(self,strands):
        if not self.is_mpi_master: return []
        locator_dict={}
        for s in strands:
            if not hasattr(s,"encoded_index_ints"):
                raise ValueError("Strand does not have fault injection information in IdealCluster")
            index = getattr(s,"encoded_index_ints")
            locator_dict[index]=locator_dict.get(index,[])+[s]
        return [ y for x,y in locator_dict.items()]
