'''
Class that implements probes that can be used to check fault injection results throughout the pipeline.
'''
from dnastorage.codec.base import *
from dnastorage.fi.fault_strand_representation import *
from dnastorage.codec_types import *
from dnastorage.util.stats import dnastats
from dnastorage.codec.base_conversion import unpack_bytes_to_indexes
from dnastorage.primer.primer_util import *
import copy
import numpy as np
import Levenshtein as ld
import sys, traceback


import logging
logger = logging.getLogger("dnastorage.fi.probes")
logger.addHandler(logging.NullHandler())

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
        self._correct_key="{}::correct_strands".format(self._name)
        self._incorrect_key="{}::incorrect_strands".format(self._name)
        self._incorrect_not_none="{}::incorrect_strands_not_none".format(self._name)
        self._first_byte_error="{}::first_error".format(self._name)
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
        first=True
        for i in range(len(base_codewords)):
            if i>=len(fault_codewords) or base_codewords[i]!=fault_codewords[i]:
                stats.inc(self._total_error_rate_key,dflt=np.zeros((len(base_codewords),)),coords=i)
                if first:
                    first=False
                    stats.inc(self._first_byte_error,dflt=np.zeros((len(base_codewords),)),coords=i)
            if i<len(fault_codewords) and base_codewords[i]!=fault_codewords[i] and fault_codewords[i]!=None:
                stats.inc(self._incorrect_not_none,dflt=np.zeros((len(base_codewords),)),coords=i)
        if(fault_codewords==base_codewords):
            stats.inc(self._correct_key)
        else:
            stats.inc(self._incorrect_key)
        stats.inc(self._strands_seen_key)
        
        return s


class DNAErrorProbe(BaseCodec,Probe):
    #This probe should be placed after a single strand codec so that analysis can be performed 
    probe_id=0
    def __init__(self,probe_name="",CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        if probe_name=="":
            self._name = "dna_probe_{}".format(DNAErrorProbe.probe_id)
        else:
            self._name = probe_name
        self._initial_dna_attr = "{}::init_dna".format(self._name)
        self._total_error_rate_key = "{}::total_hamming_errors".format(self._name)
        self._total_ed_rate_key = "{}::total_edit_errors".format(self._name)
        self._strands_seen_key = "{}::total_strands".format(self._name)
        self._correct_key="{}::correct_strands".format(self._name)
        self._incorrect_key="{}::incorrect_strands".format(self._name)
        DNAErrorProbe.probe_id+=1

    def _encode(self,s):
        setattr(s,self._initial_dna_attr,copy.copy(s.dna_strand)) #take a snapshot of the dna under an attribute specific to this isntantiation
        return s

    def _decode(self,s):
        if not hasattr(s,self._initial_dna_attr) or s.dna_strand is None:
            return s 
        base_dna = getattr(s,self._initial_dna_attr)
        fault_dna = s.dna_strand
        #now analyze error rates b/w the base and fault
        for i in range(len(base_dna)):
            if i>=len(fault_dna) or base_dna[i]!=fault_dna[i]:
                stats.inc(self._total_error_rate_key,dflt=np.zeros((len(base_dna),)),coords=i)
        if(fault_dna==base_dna):
            stats.inc(self._correct_key)
        else:
            stats.inc(self._incorrect_key)
        stats.inc(self._strands_seen_key)

        #do edit distance analysis
        editops= ld.editops(base_dna,fault_dna)
        for operation,base_index,fault_index in editops:
            if base_index<len(base_dna):
                stats.inc(self._total_ed_rate_key,dflt=np.zeros((len(base_dna),)),coords=base_index)
        return s


class FilteredDNACounter(BaseCodec,Probe):
    #This probe should be placed after DNAtoDNA probes to count filtered strands 
    probe_id=0
    def __init__(self,probe_name="",CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        if probe_name=="":
            self._name = "dna_filtered_probe_{}".format(FilteredDNACounter.probe_id)
        else:
            self._name = probe_name
        self._strands_filtered_key = "{}::filtered_strand".format(self._name)
        #stats[self._strands_filtered_key]=0
        FilteredDNACounter.probe_id+=1

    def _encode(self,s):
        return s

    def _decode(self,s):
        if s.dna_strand is None and s.is_reversed:
            stats.inc(self._strands_filtered_key)
        elif s.dna_strand is None and hasattr(s,"is_nanopore") and s.is_nanopore:
            stats.inc(self._strands_filtered_key)
        return s

class IndexDistribution(BaseCodec,Probe):
    #This probe should be placed after a single strand codec so that analysis can be performed 
    probe_id=0
    def __init__(self,probe_name="",CodecObj=None,prefix_to_match=tuple()):
        BaseCodec.__init__(self,CodecObj)

        if probe_name=="":
            self._name = "index_dist_probe_{}".format(IndexDistribution.probe_id)
        else:
            self._name = probe_name
        #stats collected during encoding, gives baseline
        self._index_dist_probe_key_encode = "{}::index_dist_encode".format(self._name)
        self._total_indexes_encode = "{}::total_indexes_encode".format(self._name)
        #stats collected during decoding, measures reality
        self._index_dist_probe_key_decode = "{}::index_dist_decode".format(self._name)
        self._total_indexes_decode = "{}::total_indexes_decode".format(self._name)
        self._total_indexes_lost = "{}::total_indexes_lost".format(self._name)
        self._seq_map = "{}::seq_map".format(self._name)
        stats[self._seq_map]={}
        self._prefix_to_match = prefix_to_match
        self._correct_key="{}::correct_indexes".format(self._name)
        self._incorrect_key="{}::incorrect_indexes".format(self._name)
        self._initial_index_ints_attr = "{}::initial_index_attribute".format(self._name)
        #register sequencing mappings
        stats.register_file(self._seq_map,"sequencing.map")
        IndexDistribution.probe_id+=1

    def _encode(self,s):
        stats.inc(self._total_indexes_encode)
        stats.inc(self._index_dist_probe_key_encode,dflt=dict(),coords=s.index_ints)
        setattr(s,self._initial_index_ints_attr,copy.copy(s.index_ints)) #take a snapshot of the indices
        return s

    def _decode(self,s):
        #pipeline doesn't make index_ints until AFTER the whole inner CW to CW cascade runs, do our own translation for now
        try:
            index_ints = tuple(unpack_bytes_to_indexes(s.codewords[0:s.index_bytes],s.index_bit_set))
        except Exception as e:
            stats.inc(self._total_indexes_lost)
            return s
        if self._prefix_to_match != index_ints[:len(self._prefix_to_match)]:
            stats.inc(self._total_indexes_lost)
        else:
            stats.inc(self._total_indexes_decode)
            stats.inc(self._index_dist_probe_key_decode,dflt=dict(),coords=index_ints)
            #if the strand we are seeing is from sequencing, record its index
            if hasattr(s,"record_id"): 
                stats[self._seq_map][s.record_id]=index_ints
            if hasattr(s,self._initial_index_ints_attr):
                if index_ints == tuple(getattr(s,self._initial_index_ints_attr)):
                    stats.inc(self._correct_key)
                else:
                    stats.inc(self._incorrect_key)
        return s

class StrandCheckProbe(BaseCodec,Probe):
    #This probe should be placed after a single strand codec so that analysis can be performed 
    probe_id=0
    def __init__(self,strands,CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        self._strands = strands
    def _encode(self,s):
        for check_strand in self._strands:
            forward_key = "{}:forward_min".format(check_strand)
            reverse_key = "{}:reverse_min".format(check_strand)
            stats[forward_key]=(0xffffffffffffffff,None,None,None)
            stats[reverse_key]=(0xffffffffffffffff,None,None,None)
            for i in range(0,len(s.dna_strand)-len(check_strand)):
                #print(s.dna_strand)
                forward_distance = hamming_distance(check_strand,s.dna_strand[i:i+len(check_strand)])
                reverse_distance = hamming_distance(reverse_complement(check_strand),s.dna_strand[i:i+len(check_strand)])
                if forward_distance == 0 or reverse_distance ==0: #im assuming any match is due to a region that should be there...
                    continue
                new_tuple_reverse = (reverse_distance,s.index_ints,i,s.dna_strand[i:i+len(check_strand)]) #get more information about the closest match
                stats[reverse_key]=min(new_tuple_reverse,stats[reverse_key],key=lambda x: x[0])
                new_tuple_forward = (forward_distance,s.index_ints,i,s.dna_strand[i:i+len(check_strand)])
                stats[forward_key]=min(new_tuple_forward,stats[forward_key], key=lambda x: x[0])
        return s

    def _decode(self,s):
        return s
