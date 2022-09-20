#!/usr/bin/python
import os
import logging
from Bio import SeqIO
from dnastorage.strand_representation import *
import h5py
import re

logger = logging.getLogger("dnastorage.util.strandinterface")
logger.addHandler(logging.NullHandler())

"""
Provide utilities for loading strands from different formats, array, fastq, fast5, etc.
"""
class BaseStrandInterface:
    def __init__(self):
        self._strands = []
    @property
    def strands(self):
        return self._strands
    @strands.setter
    def strands(self,s):
        self._strands = s
    @classmethod
    def open(self,format_type,path=""):
        if format_type=="fast5":
            return Fast5Interface(path)
        elif format_type=="array":
            return ArrayInterface()
        elif format_type=="fastq":
            return FastqInterface(path)
        elif format_type=="DNA":
            return DNAFileInterface(path)
        else:
            raise ValueError("No matching format to extract information")
        
class ArrayInterface(BaseStrandInterface):
    def __init__(self):
        BaseStrandInterface.__init__(self)

class FastqInterface(BaseStrandInterface):
    def __init__(self,path):
        BaseStrandInterface.__init__(self)
        assert os.path.exists(path)
        tmp_strands=[]
        for record in SeqIO.parse(path,"fastq"):
            tmp_strands.append(BaseDNA(dna_strand=str(record.seq)))
            (tmp_strands[-1].dna_strand,cnt) = re.subn("U","T",tmp_strands[-1].dna_strand) #replace all U with T, should work even if strand is DNA
            tmp_strands[-1].record_id=record.id
            if cnt>0: tmp_strands[-1].is_RNA=True
        self.strands=tmp_strands

class DNAFileInterface(BaseStrandInterface):
    def __init__(self,path):
        BaseStrandInterface.__init__(self)
        assert os.path.exists(path)
        tmp_strands=[]
        with open(path,'r') as in_fd:
            while True:
                s = in_fd.readline()        
                if len(s) == 0:
                    break
                s = s.strip()
                if s.startswith('%'):
                    continue
                tmp_strands.append(BaseDNA(dna_strand=s))
        self.strands = tmp_strands

class Fast5Interface(BaseStrandInterface):
    def __init__(self,path):
        BaseStrandInterface.__init__(self)
        assert os.path.exists(path)

        if os.path.isfile(path):
            self.strands=self._extract_fast5_file(path)
        elif os.path.isdir(path):
            #Assumed that if a direct file is not given, then we are using Oxford Nanopore's fast5 directory layout
            for p in os.listdir(path):
                if p!="fast5_pass" and p!="fast5_fail": raise ValueError("Nanopore fast5 directory format not followed")
                fast5_dir = os.path.join(path,p)
                for f5 in os.listdir(fast5_dir):
                    if not ".fast5" in f5: raise ValueError("Non fast5 file in fast5 directory")
                    self.strands = self.strands + self._extract_fast5_file(os.path.join(fast5_dir,f5))
        else:
            raise ValueError("Fast5 path is not a file or directory")
        
    def _extract_fast5_file(self,path):
        tmp_strands =[]
        assert os.path.exists(path) and os.path.isfile(path)
        #extract strands from a given file, along with other interesting meta data
        hd_fd=h5py.File(path,'r',libver="latest")
        root_group = hd_fd["/"]
        for key in root_group.keys():
            if not isinstance(root_group[key],h5py.Group): continue
            fastq_read = root_group[key]["Analyses/Basecall_1D_000/BaseCalled_template/Fastq"][()].decode("utf-8")
            parsed_fastq=fastq_read.split("\n")
            quality_scores = parsed_fastq[-2]
            DNA = parsed_fastq[1].replace("U","T")
            tmp_strands.append(BaseDNA(dna_strand=DNA))
            tmp_strands[-1].record_id = root_group[key]["Raw"].attrs["read_id"].decode("utf-8")
            tmp_strands[-1].quality_scores=quality_scores
            tmp_strands[-1].is_RNA=True
            tmp_strands[-1].channel_id = int(root_group[key]["channel_id"].attrs["channel_number"])
        return tmp_strands

if __name__=="__main__":
    #test out the interface on fast5 files
    path = "/tuck_data/kvolkel/paul_fastq/211013_cgaa_ivt_qch.fastq"
    strand_interface = BaseStrandInterface.open("fastq",path)
    assert(isinstance(strand_interface,FastqInterface))
    print("Total Number of strands Read {}".format(len(strand_interface.strands)))
    for i in range(10):
        strand = strand_interface.strands[i]
        print("===========Printing Strand {}=========".format(i)) 
        print(strand.dna_strand)
        print(strand.record_id)
      
    
