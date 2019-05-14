import starcode_bindings


'''
filename: cluster.py

Description: This file contains classes that provide the interface to the required functionalities
of a clustering algorithm. Namely, a 'cluster' method must exist, and must calculate the set of clusters based on the input library 'strands' array. A 'get_cluster' method must also exist that returns a cluster, and upon no remaining clusters to analyze returns the None object. It does not matter how those clustering algorithms are implemented (e.g. may be in C or in Python), as long as these interfaces are defined.
'''


class starcode:
    def __init__(self,cluster_algorithm):
        self._cluster_algorithm=cluster_algorithm
    def cluster(self,strands):
        #invoke the c extension of the top level starcode function to start starcode clustering process
        if self._cluster_algorithm=='MP':
            starcode_bindings.starcode_top(10,4,0,1,strands,len(strands))            
    def get_cluster(self):
        #invoke the c extension of the getter function for starcode
        return starcode_bindings.starcode_get_cluster();

