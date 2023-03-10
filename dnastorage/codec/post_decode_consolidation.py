from dnastorage.codec.base import *
from dnastorage.codec_types import *
from dnastorage.strand_representation import *
import dnastorage.codec.base_conversion as base_conversion
from collections import Counter
import copy
'''
Consolidation is the process of clustering, e.g. taking multiple copies of the same strand and providing a representative
'''

class SimpleMajorityVote(BaseCodec,CWConsolidate):
    def __init__(self):
        BaseCodec.__init__(self)
        CWConsolidate.__init__(self)
    def _encode(self):
        raise NotImplemented()
    
    def _decode(self,strands):
        key_value = {}
        strand_array = []
        for s in strands:
            #print "get key,value for strand {}".format(s)
            key =tuple(s.index_ints)
            #print value
            if hasattr(s,"alignment_weight"): #allow alignment weights to be 
                key_value[key]=key_value.get(key,[]) + [s]*s.alignment_weight
            else: 
                key_value[key]=key_value.get(key,[]) + [s]
        for key in key_value:
            data=[]
            if len(key_value[key])==1:
                # no reason to vote on a single strand
                strand_array.append(key_value[key][0])
                continue
            #print "processing key values"
            mx = max([len(s.codewords) for s in key_value[key]])
            for x in range(0,mx):
                #get a list of values that belong to the same location
                # if there's a smaller strand, due to error, this assert
                # will fire before attempting to access an illegal index
                same_position_values = []
                for s in key_value[key]:
                    if x < len(s.codewords):
                        same_position_values.append(s.codewords[x])
                #print same_position_values
                cnt=Counter(same_position_values)
                most_common_value=cnt.most_common(1)[0][0]
                #print most_common_value
                data.append(most_common_value)
                #put the key and data into an array
            #try to find the representive strand to copy over class contents
            min_distance=0xFFFFFFFFFFFFFFFF
            min_strand=key_value[key][0]
            for s in key_value[key]:
                distance=0
                for i,v in enumerate(data):
                    if i>=len(s.codewords):
                        distance+=len(data)-len(s.codewords)
                        break
                    if not v==s.codewords[i]:
                        distance+=1
                if distance<min_distance:
                    min_strand=s
            new_strand=copy.copy(min_strand)
            strand_array.append(new_strand)
   
        return strand_array
