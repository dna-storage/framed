#This file holds classes for handling multiple reads of the same strands
from collections import Counter
from cluster import *
from cluster_analyzers import *
import faulthandler; faulthandler.enable(); 
"""
data_vote_simple is a policy that collects data with the same index, and
uses a majority voting policy in order to infer what the correct data is.
The resulting key,value pairs will be passed to the upper level of the decoder
e.g. for the RS decoder, the level that applies ECC codes to inner and outer
portions of data.
"""

class data_vote_simple:
    def __init__(self,Codec):
        self._Codec=Codec
    def process(self,strands):
        key_value={}
        key_val_array=[]
        #print "processing strands"
        #seperare the strands up to their indexes
        for s in strands:
            #print "get key,value for strand {}".format(s)
            key,value=self._Codec.decode(s)
            #print value
            if key not in key_value and key is not None and value is not None:
                key_value[key]=[]
                key_value[key].append(value)
            elif key is not None and value is not None:
                key_value[key].append(value)
            #print "stored key,value"
        #go through each index and figure out the majority data
        for key in key_value:
            data=[]
            #print "processing key values"
            mx = max([len(data_strand) for data_strand in key_value[key]])
            for x in range(0,mx):
                #get a list of values that belong to the same location
                # if there's a smaller strand, due to error, this assert
                # will fire before attempting to access an illegal index
                #if x >= min([len(data_strand) for data_strand in key_value[key]]):
                #    print x,"[",",".join([str(len(data_strand)) for data_strand in key_value[key]]),"]"
                #assert x < min([len(data_strand) for data_strand in key_value[key]])
                same_position_values = []
                for data_strand in key_value[key]:
                    if x < len(data_strand):
                        same_position_values.append(data_strand[x])
                #print same_position_values
                cnt=Counter(same_position_values)
                most_common_value=cnt.most_common(1)[0][0]
                #print most_common_value
                data.append(most_common_value)
            #put the key and data into an array
            key_val_array.append((key,data))
        return key_val_array
    def process_tuple(self,strands):
        key_value={}
        key_val_array=[]
        #seperare the strands up to their indexes
        for s in strands:
            if(s[1]>0):
                key,value=self._Codec.decode(s[0])
                if key is not None and value is not None:
                    if key not in key_value:
                        key_value[key]=[]
                        for i in range(0,s[1]):
                            key_value[key].append(value)
                    else:
                        for i in range(0,s[1]):
                            key_value[key].append(value)

        #go through each index and figure out the majority data
        for key in key_value:
            data=[]
            for x in range(0,len(key_value[key][0])):
                #get a list of values that belong to the same location
                same_position_values=[data_strand[x] for data_strand in key_value[key]]
                #print same_position_values
                cnt=Counter(same_position_values)
                most_common_value=cnt.most_common(1)[0][0]
                #print most_common_value
                data.append(most_common_value)
            #put the key and data into an array
            key_val_array.append((key,data))
            #print key
            #print data
        return key_val_array



'''
The clusterig_handler is a version of strand handler that leverages a clustering algorithm
specified by the cluster_algorithm __init__ argument. The cluster_algorithm must be a class
that supports 2 interfaces: 1. To start the clustering process and finish, and 2. a method
that iteratively returns cluster by cluster. For each cluster the clustering_handler
will utilize the cluster analysis algorithm that is specified by the cluster_analysis_algorithm
argument to the __init__ function for the clustering_handler class. The cluster_analysis_algorithm
should return a single strand that is supposed to act as the representative of that cluster.
'''
class clustering_handler:
    def __init__(self,Codec,cluster_analysis_algorithm,cluster_algorithm):
        self._cluster_algorithm=cluster_algorithm #class which implements the type of clustering
        self._cluster_analysis_algorithm=cluster_analysis_algorithm #how to interpret the cluster
        self._data_vote=data_vote_simple(Codec) #simple data vote will be used to clean up results from clustering and post clustering analysis
    def process(self,strands):
        corrected_strand_array=[]
        #have the cluster algorithm cluster the strands 
        
        self._cluster_algorithm.cluster(strands)
        strand_cluster=self._cluster_algorithm.get_cluster()
        while strand_cluster is not None:
            #print strand_cluster 
            if self._cluster_analysis_algorithm == "BMA":
                s=cluster_analyze_majority(strand_cluster,200)
            elif self._cluster_analysis_algorithm=="BMA_ED":
                s=cluster_analyze_majority(strand_cluster,200)
            elif self._cluster_analysis_algorithm=="ED":
                s=cluster_analyze_ed(strand_cluster,200,2)
            corrected_strand_array.append(s)
            strand_cluster=self._cluster_algorithm.get_cluster()
        #print "analyzing corrected strands"
        key_value_data=self._data_vote.process(corrected_strand_array)
        #print key_value_data
        return key_value_data
 

