#This file holds classes for handling multiple reads of the same strands
from collections import Counter
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
        #seperare the strands up to their indexes
        for s in strands:
            key,value=self._Codec.decode(s)
            #print value
            if key not in key_value:
                key_value[key]=[]
                key_value[key].append(value)
            else:
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
