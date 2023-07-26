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


def calculate_forward_burst(editops,index,edit_type):
    current_pos = index
    next_pos = current_pos+1
    total_burst=0
    while next_pos<len(editops) and editops[next_pos][1]==(editops[current_pos][1]+1) and editops[next_pos][0]==edit_type:
        total_burst+=1
        current_pos=next_pos
        next_pos+=1
    return total_burst

def calculate_reverse_burst(editops,index,edit_type):
    total_burst=0
    current_pos=index
    next_pos=current_pos-1
    while next_pos>0 and editops[next_pos][1]==(editops[current_pos][1]-1) and editops[next_pos][0]==edit_type:
        total_burst+=1
        current_pos=next_pos
        next_pos-=1
    return total_burst


class DNAErrorProbe(BaseCodec,Probe):
    #This probe should be placed after a single strand codec so that analysis can be performed 
    probe_id=0
    count=0
    
    @property
    def dna_attr(self):
        return self._initial_dna_attr
    @dna_attr.setter
    def dna_attr(self,attr:str):
        self._initial_dna_attr=attr
    def __init__(self,probe_name="",calculate_hamming=False,CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        if probe_name=="":
            self.name = "dna_probe_{}".format(DNAErrorProbe.probe_id)
        else:
            self.name = probe_name
        self.dna_attr = "{}::init_dna".format(self.name)
        #keys for edit distance
        self._total_ed_rate_key = "{}::total_edit_errors".format(self.name) #number of edits at a position of all types
        self._inser_ed_rate_key = "{}::insertion_edit_errors".format(self.name) #number of insertions at a position
        self._del_ed_rate_key = "{}::deletion_edit_errors".format(self.name) #number of deletions at a position
        self._sub_ed_rate_key = "{}::substitution_edit_errors".format(self.name)#number of substitution at a position
        self._DNA_strands_seen_key = "{}::ed_total_strands".format(self.name)#total strands that pass through counter
        self._DNA_correct_key="{}::ed_correct_strands".format(self.name)#number of strands with 0 base errors
        self._DNA_incorrect_key="{}::ed_incorrect_strands".format(self.name)#number of strands with at least 1 base error
        self._DNA_del_burst_len_key = "{}::del_burst_length".format(self.name)#length of error for deletions (considers forward and backward from a location)
        self._DNA_del_burst_rate_key = "{}::del_burst_rate".format(self.name)#given a deletion error, the number of dels that are bursts
        self._DNA_strand_length_key="{}::strand_length_hist".format(self.name) #histogram for strand lengths
        self._DNA_strand_error_key="{}::strand_error_hist".format(self.name) #histogram for strand error rates
        self._DNA_max_kmer_error_key = "{}::kmer_max_error_hist".format(self.name) #histogram for max kmer error window
        self._DNA_error_pattern_key = "{}::error_pattern_hist".format(self.name) #histogram, but it will be a dictionary given we are mapping 
        self._DNA_pattern_start_rate_key = "{}::pattern_rate_key".format(self.name) #tracks rate of error patterns in the strand, instead of individual errors
        self._DNA_pattern_location_key="{}::pattern_locations".format(self.name) #tracks the location that error patterns occur throughout the strand
        self._DNA_hamming_distance_key = "{}::hamming_distance_key".format(self.name) #tracks hamming distance for strands that are the same length
        self._DNA_hamming_distance_strand_count_key = "{}::total_hamming_distance_strands".format(self.name)
        
        stats.register_hist(self._DNA_strand_length_key)
        stats.register_hist(self._DNA_strand_error_key)
        stats.register_hist(self._DNA_max_kmer_error_key)
        if not self._DNA_error_pattern_key in stats: stats[self._DNA_error_pattern_key]={}
        if not self._DNA_pattern_location_key in stats: stats[self._DNA_pattern_location_key]={}
        DNAErrorProbe.probe_id+=1
        self._calculate_hamming=calculate_hamming
    def _encode(self,s):
        setattr(s,self.dna_attr,copy.copy(s.dna_strand)) #take a snapshot of the dna under an attribute specific to this isntantiation
        return s
    def _decode(self,s):
        if not hasattr(s,self.dna_attr) or s.dna_strand is None:
            return s
        stats.append(self._DNA_strand_length_key,len(s.dna_strand))
        base_dna = getattr(s,self.dna_attr)
        fault_dna = s.dna_strand
        if(fault_dna==base_dna):
            stats.inc(self._DNA_correct_key)
        else:
            stats.inc(self._DNA_incorrect_key)
        if len(fault_dna)==len(s.dna_strand) and self._calculate_hamming:
            stats.inc(self._DNA_hamming_distance_strand_count_key,dflt=0)
            for index,(i,j) in enumerate(zip(fault_dna,base_dna)):
                if i!=j: stats.inc(self._DNA_hamming_distance_key,dflt=np.zeros((len(fault_dna),)),coords=index)
            
        stats.inc(self._DNA_strands_seen_key)
        #do edit distance analysis
        editops= ld.editops(base_dna,fault_dna)
        edit_strand_vis,applied_edits,max_kmer_edits,pattern_dist,pattern_starts= calculate_edit_list(editops,len(base_dna),kmer_length=50,pattern=True)
        for edit_op_index,(operation,base_index,fault_index) in enumerate(editops):
            if base_index<len(base_dna):
                stats.inc(self._total_ed_rate_key,dflt=np.zeros((len(base_dna),)),coords=base_index)
                if operation=="delete":
                    stats.inc(self._del_ed_rate_key,dflt=np.zeros((len(base_dna),)),coords=base_index)
                    forward_burst = calculate_forward_burst(editops,edit_op_index,"delete")
                    reverse_burst = calculate_reverse_burst(editops,edit_op_index,"delete")
                    total_burst=forward_burst+reverse_burst
                    stats.inc(self._DNA_del_burst_len_key,i=total_burst,dflt=np.zeros((len(base_dna),)),coords=base_index)
                    if total_burst>0: stats.inc(self._DNA_del_burst_rate_key,dflt=np.zeros((len(base_dna),)),coords=base_index)
                elif operation=="insert":
                    stats.inc(self._sub_ed_rate_key,dflt=np.zeros((len(base_dna),)),coords=base_index)
                elif operation=="replace":
                    stats.inc(self._inser_ed_rate_key,dflt=np.zeros((len(base_dna),)),coords=base_index)
            else: break
        stats.append(self._DNA_strand_error_key,applied_edits) #collect histogram of errors per strand
        stats.append(self._DNA_max_kmer_error_key,max_kmer_edits) #collect histogram of max errors in kmer window
        for pattern in pattern_dist: stats[self._DNA_error_pattern_key][pattern] = stats[self._DNA_error_pattern_key].get(pattern,0)+pattern_dist[pattern]
        for start,pattern in pattern_starts: #sorted(pattern_starts,key=lambda x: x[1]):
            stats.inc(self._DNA_pattern_start_rate_key,dflt=np.zeros((len(base_dna),)),coords=start)
            #if pattern not in stats[self._DNA_pattern_location_key]: stats[self._DNA_pattern_location_key][pattern]=np.zeros((len(base_dna),))
            #stats[self._DNA_pattern_location_key][pattern][start]+=1 #increment location that pattern shows up
        return s

class CodewordErrorRateProbe(BaseCodec,Probe):
    #This probe should be placed after a single strand codec so that analysis can be performed 
    probe_id=0
    @property
    def cw_attr(self):
        return self._initial_codeword_attr
    @cw_attr.setter
    def cw_attr(self,attr: str):
        self._initial_codeword_attr=attr

    @property
    def dna_attr(self):
        return self._initial_dna_attr
    @dna_attr.setter
    def dna_attr(self,attr:str):
        self._initial_dna_attr=attr
        #create new dna error probes for the given dna attribute
        self._correct_decode_strand_dna_probe =DNAErrorProbe("{}::correct".format(self.name))
        self._correct_decode_strand_dna_probe.dna_attr=attr
        self._incorrect_decode_strand_dna_probe=DNAErrorProbe("{}::incorrect".format(self.name))
        self._incorrect_decode_strand_dna_probe.dna_attr=attr
    def __init__(self,probe_name="",dna_hook=None,CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        if probe_name=="":
            self.name = "codeword_probe_{}".format(CodewordErrorRateProbe.probe_id)
        else:
            self.name = probe_name
        #instantiate some DNAErrorProbes that we can activate in decode/not decoded scenarios
        self._incorrect_decode_strand_dna_probe=None
        self._correct_decode_strand_dna_probe=None
        if dna_hook:
            self._correct_decode_strand_dna_probe =DNAErrorProbe("{}::correct".format(self.name))
            self._correct_decode_strand_dna_probe.dna_attr=dna_hook
            self._incorrect_decode_strand_dna_probe=DNAErrorProbe("{}::incorrect".format(self.name))
            self._incorrect_decode_strand_dna_probe.dna_attr=dna_hook
        self.cw_attr = "{}::init_codewords".format(self.name)
        self._total_error_rate_key = "{}::total_errors".format(self.name)
        self._strands_seen_key = "{}::total_strands".format(self.name)
        self._correct_key="{}::correct_strands".format(self.name)
        self._incorrect_key="{}::incorrect_strands".format(self.name)
        self._incorrect_not_none="{}::incorrect_strands_not_none".format(self.name)
        self._first_byte_error="{}::first_error".format(self.name)
        self._decoded_strands_key="{}::decoded_strands".format(self.name)#actual sequences that decoded
        self._failed_strands_key="{}::failed_decode_strands".format(self.name) #sequences that failed
        stats.register_file(self._decoded_strands_key,"decoded_strands.txt")
        stats.register_file(self._failed_strands_key,"failed_strands.txt")
        CodewordErrorRateProbe.probe_id+=1
    def _encode(self,s):
        setattr(s,self.cw_attr,copy.copy(s.codewords)) #take a snapshot of the codewords under an attribute specific to this isntantiation
        return s
    def _decode(self,s):
        if not hasattr(s,self.cw_attr): return s #might have a leak over from some other encoding/decoding pipeline, cannot rely on this being true always, best compromise is to return
        base_codewords = getattr(s,self.cw_attr)
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
        if(fault_codewords==base_codewords): #strand decoded bytes correctly
            stats.inc(self._correct_key)
            #stats.append(self._decoded_strands_key,s.dna_strand)
            if self._correct_decode_strand_dna_probe:
                self._correct_decode_strand_dna_probe._decode(s) #call error rate analysis
        else:
            #stats.append(self._failed_strands_key,s.dna_strand)
            if self._incorrect_decode_strand_dna_probe:#strand decoded bytes incorrectly
                self._incorrect_decode_strand_dna_probe._decode(s) #call error rate analysis
            stats.inc(self._incorrect_key)
        stats.inc(self._strands_seen_key)
        return s

class FilteredDNACounter(BaseCodec,Probe):
    #This probe should be placed after DNAtoDNA probes to count filtered strands 
    probe_id=0
    def __init__(self,probe_name="",CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        if probe_name=="":
            self.name = "dna_filtered_probe_{}".format(FilteredDNACounter.probe_id)
        else:
            self.name = probe_name
        self._strands_filtered_key = "{}::filtered_strand".format(self.name)
        FilteredDNACounter.probe_id+=1
    def _encode(self,s):
        return s
    def _decode(self,s):
        if s.dna_strand is None and s.is_reversed:
            stats.inc(self._strands_filtered_key)
        elif s.dna_strand is None and hasattr(s,"is_RNA") and s.is_RNA:
            stats.inc(self._strands_filtered_key)
        return s

class IndexDistribution(BaseCodec,Probe):
    #This probe should be placed after a single strand codec so that analysis can be performed 
    probe_id=0
    def __init__(self,probe_name="",CodecObj=None,prefix_to_match=tuple()):
        BaseCodec.__init__(self,CodecObj)
        if probe_name=="":
            self.name = "index_dist_probe_{}".format(IndexDistribution.probe_id)
        else:
            self.name = probe_name
        #stats collected during encoding, gives baseline
        self._index_dist_probe_key_encode = "{}::index_dist_encode".format(self.name)
        self._total_indexes_encode = "{}::total_indexes_encode".format(self.name)
        #stats collected during decoding, measures reality
        self._index_dist_probe_key_decode = "{}::index_dist_decode".format(self.name)
        self._total_indexes_decode = "{}::total_indexes_decode".format(self.name)
        self._total_indexes_lost = "{}::total_indexes_lost".format(self.name)
        self._total_indexes_lost_barcode_mismatch = "{}::total_indexes_lost_barcode_mismatch".format(self.name)
        self._seq_map = "{}::seq_map".format(self.name)
        stats[self._seq_map]={}
        self._prefix_to_match = prefix_to_match
        self._correct_key="{}::correct_indexes".format(self.name)
        self._incorrect_key="{}::incorrect_indexes".format(self.name)
        self._initial_index_ints_attr = "{}::initial_index_attribute".format(self.name)
        #register sequencing mappings
        self._ideal_index_dist_attr = stats.get_next_name("{}::ideal_decode_dist".format(self.name))
        stats.register_file(self._seq_map,"sequencing.map")
        IndexDistribution.probe_id+=1

    def _encode(self,s):
        stats.inc(self._total_indexes_encode)
        stats.inc(self._index_dist_probe_key_encode,dflt=dict(),coords=s.index_ints)
        setattr(s,self._initial_index_ints_attr,copy.copy(s.index_ints)) #take a snapshot of the indices
        return s

    def _decode(self,s):
        #if hasattr(s,self._initial_index_ints_attr):
        #    stats.inc(self._ideal_index_dist_attr,dflt=dict(),coords=getattr(s,self._initial_index_ints_attr))
        #pipeline doesn't make index_ints until AFTER the whole inner CW to CW cascade runs, do our own translation for now
        try:
            index_ints = tuple(unpack_bytes_to_indexes(s.codewords[0:s.index_bytes],s.index_bit_set))
        except Exception as e:
            stats.inc(self._total_indexes_lost)
            return s
        if self._prefix_to_match != index_ints[:len(self._prefix_to_match)]:
            stats.inc(self._total_indexes_lost_barcode_mismatch)
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
    def __init__(self,strand_path,CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        self._strands = []
        with open(strand_path,'r') as strands:
            for s in strands:
                self._strands.append(s.strip())
    def _encode(self,s):
        for check_strand in self._strands:
            forward_key = "{}:forward_min".format(check_strand)
            reverse_key = "{}:reverse_min".format(check_strand)
            stats[forward_key]=(0xffffffffffffffff,None,None,None)
            stats[reverse_key]=(0xffffffffffffffff,None,None,None)
            for i in range(0,len(s.dna_strand)-len(check_strand)):
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

class HookProbe(BaseCodec,Probe):
    """
    HookProbe does not do any real analysis, instead it is used to provide a hook to an attribute that can be used by another probe
    Useful for analysis that relies on snapshotting information that does not exist yet for another probe.
    Example: An analysis wants to analyze only DNA strands that pass through encoding correctly.
    strand_attr: strand attribute to hook into during encoding
    probe_name: name of the snapshot attribute
    """
    probe_id=0
    def __init__(self,strand_attr,probe_name="",CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        if probe_name=="":
            self.name = "hook_probe_{}".format(HookProbe.probe_id)
        else:
            self.name = probe_name
        self.name = "{}::{}".format(self.name,strand_attr)
        self._strand_attr=strand_attr
        HookProbe.probe_id+=1
    def _encode(self,s):
        assert hasattr(s,self._strand_attr)
        setattr(s,self.name,copy.copy(getattr(s,self._strand_attr))) #take a snapshot of an attribute in the dnastrand
        return s
    def _decode(self,s):
        return s
