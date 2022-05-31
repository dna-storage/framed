'''
Class that implements probes that can be used to check fault injection results throughout the pipeline.
'''
from dnastorage.codec.base import *
from dnastorage.fi.fault_strand_representation import *
from dnastorage.codec_types import *
from dnastorage.util.stats import dnastats
import copy
import numpy as np



class CodewordErrorRateProbe(BaseCodec,Probe):
    #This probe should be placed after a single strand codec so that analysis can be performed 
    probe_id=0
    def __init__(self,probe_name="",CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        if probe_name=="":
            self._name = "codeword_probe_{}".format(CodewordErrorRateProbe.probe_id)
        else:
            self._name = probe_name
        self._initial_codeword_attr = "{}::init_codewords".format(self._name)
        self._total_error_rate_key = "{}::total_errors".format(self._name)
        self._strands_seen_key = "{}::total_strands".format(self._name)
        CodewordErrorRateProbe.probe_id+=1

    def _encode(self,s):
        setattr(s,self._initial_codeword_attr,copy.copy(s.codewords)) #take a snapshot of the codewords under an attribute specific to this isntantiation
        return s

    def _decode(self,s):
        if not isinstance(s,FaultDNA): #only analyze fault strands
            return
        if not isinstance(s,BaseDNA): #only analyze strands
            return
        if not hasattr(s,self._initial_codeword_attr): return s #might have a leak over from some other encoding/decoding pipeline, cannot rely on this being true always, best compromise is to return
        base_codewords = getattr(s,self._initial_codeword_attr)
        fault_codewords = s.codewords
        #now analyze error rates b/w the base and fault
        for i in range(len(base_codewords)):
            if i>=len(fault_codewords) or base_codewords[i]!=fault_codewords[i]:
                stats.inc(self._total_error_rate_key,dflt=np.zeros((len(base_codewords),)),coords=i)
            
        stats.inc(self._strands_seen_key)

        return s
        

