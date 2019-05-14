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
            if key not in key_value and key is not None and value is not None:
                key_value[key]=[]
                key_value[key].append(value)
            elif key is not None and value is not None:
                key_value[key].append(value)
        #go through each index and figure out the majority data
        for key in key_value:
            data=[]
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
            #print key
            #print data
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
