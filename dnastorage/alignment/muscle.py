from dnastorage.alignment.basealignment import *
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import os
import tempfile
import random
import subprocess

"""
Performs MSA using muscle
num_strands: maximum number of strands to use in multi sequence alignment
"""
class MuscleAlign(BaseAlignment):
    def __init__(self,num_strands):
        BaseAlignment.__init__(self)
        #add some parameters for muscle alignment
        self._muscle_exe = "muscle"
        self._num_strands = num_strands
    def Run(self, cluster):
        in_fd,in_path = tempfile.mkstemp()
        os.close(in_fd)
        cluster_to_align = cluster
        if len(cluster)>self._num_strands:
            cluster_to_align = random.choices(cluster,k=self._num_strands)
        with open(in_path,"w+") as tmpfile:#os.fdopen(in_fd,'w') as tmpfile:
            #need to write the clsuter to a temporary file, best we can do right now given the interfaces to muscle
            for i,s in enumerate(cluster_to_align): #each element in cluster should at least be BaseDNA
                tmpfile.write(">S{}\n{}\n".format(i,s.dna_strand))
        out_fd,out_path = tempfile.mkstemp()
        os.close(out_fd)
        #launch muscle application
        stdout_arg = subprocess.PIPE
        stderr_arg = subprocess.PIPE
        stdin = None
        child_process = subprocess.Popen(
            [self._muscle_exe, "-align",in_path,"-output",out_path],
            stdin=subprocess.PIPE,
            stdout=stdout_arg,
            stderr=stderr_arg,
            universal_newlines=True,
            cwd=None,
            env=None
        )
        stdout_str, stderr_str = child_process.communicate(stdin)
        try:
            msa =  AlignIO.read(out_path,"fasta")
        except Exception as e:
            print(e)
            print("Length of Cluster: {}".format(len(cluster)))
            exit(1)
        aligned_cluster = []
        for aligned_strand,strand in zip(msa,cluster_to_align):
            strand.dna_strand=str(aligned_strand.seq)
            aligned_cluster.append(strand) #write into the dna strands their aligned versions
        os.unlink(in_path)
        os.unlink(out_path)
        return aligned_cluster



if __name__=="__main__":
    import copy
    import random
    from dnastorage.strand_representation import *
    from dnastorage.codec.DNAConsolidatemodels import *
    import Levenshtein as ld
    base_sub = {"A":["G","C","T"],"G":["A","C","T"],"C":["A","G","T"],"T":["G","A","C"]}
    bases = ["A","G","C","T"]
    #generate a couple random strands 
    strand1 =  "".join(random.choices(["A","C","G","T"],k=200))
    print("Original Strand {}".format(strand1))
    total_strands=[]
    #create some copies
    for i in range(0,4):
        s = BaseDNA(dna_strand=strand1)
        locations = random.choices(range(len(s.dna_strand)),k=20)
        for i in sorted(locations,reverse=True):
            if i >=len(s.dna_strand): continue
            base = s.dna_strand[i]
            x=random.choices([0,1,2],k=1)[0]
            if x==0:
                #sub
                s.dna_strand = s.dna_strand[0:i]+random.choices(base_sub[base],k=1)[0]+s.dna_strand[i+1::]
            elif x==1:
                #del
                s.dna_strand = s.dna_strand[0:i]+s.dna_strand[i+1::]
            elif x==2:
                #ins
                s.dna_strand = s.dna_strand[0:i]+random.choices(bases,k=1)[0]+s.dna_strand[i::]
        total_strands.append(s)
    #shuffle strands around
    random.shuffle(total_strands)
    align = MuscleAlign(len(total_strands))
    alignment = align.Run(total_strands)
    print("Alignment Done")
    #print out the clusters
    for s in alignment:
        print(s.dna_strand)
    print("Consensus vs original")
    x = BasicDNAClusterModel(None,None)
    s = x._consensus_from_alignment(alignment)
    print(s)
    print(strand1)
    print(ld.distance(s,strand1))
