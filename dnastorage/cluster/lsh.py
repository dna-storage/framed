import numpy as np
import Levenshtein as ld
from dnastorage.cluster.basecluster import *
from dnastorage.util.mpi_utils import *

'''
Clusters strands based on locality sensitive hashing based on algorithm from ***Low cost DNA data storage using
photolithographic synthesis and advanced information reconstruction and error correction***.

m: number of min hashes to compute
k: kmer length to calculate shingles for
k_lsh: number of min hashes to use per band
ell_lsh: total number of bands to try out
l: sub-strand length to use for calculating kmers
'''
class LocalitySensitiveHashCluster(BaseCluster):
    base_table = {"A":0,"G":1,"C":2,"T":3}
    def __init__(self,m,k,k_lsh,ell_lsh,l):
        BaseCluster.__init__(self)
        self._m = m
        self._k = k
        self._k_lsh = k_lsh
        self._ell_lsh = ell_lsh
        self._length = l
        self._rng = np.random.RandomState(69) #should seed rng to agreed upon state
        
        
        
    def Run(self,strands):
        if self.mpi:
            strands = object_scatter(strands,self.mpi)
        assert len(strands)>0
        #need to process kmer min hashes
        strands=self.calculate_min_hashes(strands) #min_hashes is an strandxm matrix of min hashes
        if self.mpi:
            strands = object_gather(strands,self.mpi)
        #let rank 0 process bands for each strand
        if not self.is_mpi_master: return []
        clusters={} #need to make sure clusters do not have any repeats
        for pairs in self.calculate_pairs(strands):
            for pair in pairs:
                p0 = min(pair,key=id)
                p1 = max(pair,key=id)
                if p0  in clusters:
                    clusters[p0].add(p1)
                if p1 in clusters:
                    clusters[p1].add(p0)
                if p0 not in clusters and p1 not in clusters:
                    clusters[p0]=set([p0])
        #merge clusters
        dead_set=set([])
        return_list=[]
        for c in clusters:
            if c in dead_set: continue
            merged_cluster=set([])
            self._hcluster(clusters,c,dead_set,merged_cluster)
            if len(merged_cluster)>3: return_list.append(list(merged_cluster))
        return return_list


    def _hcluster(self,clusters,c,dead_set,merged_cluster):
        merged_cluster.add(c)
        if c in dead_set or c not in clusters: return
        dead_set.add(c)
        for strand in clusters[c]:
            if strand is c: continue
            self._hcluster(clusters,strand,dead_set,merged_cluster)

    def calculate_pairs(self,strand_and_sigs):
        for i in range(self._ell_lsh):
            min_sig_indexes = self._rng.permutation(self._m)[:self._k_lsh]
            d = {}
            for strand, min_hashes in strand_and_sigs:
                lsh_hash = 0
                for i,index in enumerate(min_sig_indexes):
                    lsh_hash += min_hashes[index]*(4**self._k)**i
                d[lsh_hash] = d.get(lsh_hash,[]) + [strand]
            #make candidate pairs
            pairs = set()
            for group in d.values():
                s_group = sorted(group,key=id)
                g0 = s_group[0]
                pairs=pairs.union(set([(g0,_) for _ in s_group[1::] if ld.distance(g0.dna_strand,_.dna_strand)<0.338*len(g0.dna_strand)]))
            yield pairs
                
    def calculate_min_hashes(self,dna_strands):
        tables = [self._rng.permutation(4**self._k) for i in range(self._m)]
        min_hashes = []
        for s in dna_strands:
            kmers = self.calculate_kmers(s.dna_strand[:self._length])
            kmers = [self.convert_dna_to_index(_) for _ in kmers] #convert kmers to indexes
            min_hashes.append((s,[min([table[_] for _ in kmers]) for table in tables]))
        return min_hashes

    def convert_dna_to_index(self,strand):
        index=0
        for i,base in enumerate(strand):
            index+=LocalitySensitiveHashCluster.base_table[base]*4**i
        return index

    def calculate_kmers(self,strand):
        kmers=[]
        for i in range(0,len(strand)-self._k+1):
            kmers.append(strand[i:i+self._k])
        return kmers


if __name__=="__main__":
    import copy
    import random
    from dnastorage.strand_representation import *
    base_sub = {"A":["G","C","T"],"G":["A","C","T"],"C":["A","G","T"],"T":["G","A","C"]}
    bases = ["A","G","C","T"]
    #generate a couple random strands 
    strand1 =  "".join(random.choices(["A","C","G","T"],k=200))
    strand2 =  "".join(random.choices(["A","C","G","T"],k=200))
    total_strands=[]
    original_strands=[strand1,strand2]
    #create some copies
    for index,d in enumerate(original_strands):
        for i in range(0,150):
            s = BaseDNA(dna_strand=d)
            s.index=index
            s.strand_index = index*10+i
            locations = random.choices(range(len(s.dna_strand)),k=30)
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
    #instantiate cluster
    lsh_m_sigs = 50 #number of signatures to make per strand
    lsh_kmer = 5 #kmer length to build a signature out of
    lsh_sim = 0.5 #simularity of strands
    lsh_sig_samples = 4 #samples to make on signatures
    lsh_sample_length = 500  #length of strand to consider when hashing
    cluster = LocalitySensitiveHashCluster(lsh_m_sigs,lsh_kmer,lsh_sig_samples,int(1/(lsh_sim**lsh_sig_samples)),lsh_sample_length)
    clusters = cluster.Run(total_strands)
    print("Clustering Done")
    #print out the clusters
    for c in clusters:
        print("Printing Cluster length {}".format(len(c)))
        for s in c:
            print("{} {}".format(s.index,s.strand_index))
