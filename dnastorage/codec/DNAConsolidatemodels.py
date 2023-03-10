import copy
from dnastorage.codec_types import *
from dnastorage.codec.base import *
from dnastorage.util.mpi_utils import *
from dnastorage.cluster.basecluster import *
from dnastorage.alignment.basealignment import *
from dnastorage.util.stats import *

class BasicDNAClusterModel(DNAConsolidate,BaseCodec):
    id_num=0
    def __init__(self,cluster,align,name=None):
        if cluster!=None:assert isinstance(cluster,BaseCluster)
        if align!=None: assert isinstance(align,BaseAlignment)
        self._cluster=cluster
        self._align=align
        BaseCodec.__init__(self)
        DNAConsolidate.__init__(self)
        basename = ""
        if name is None:
            basename="DNA_cluster_{}".format(BasicDNAClusterModel.id_num)
        else:
            basename=name
        self._cluster_purity_key="{}::purity".format(basename)
        self._cluster_size_key="{}::cluster_size".format(basename)
        stats.register_hist(self._cluster_purity_key)
        stats.register_hist(self._cluster_size_key)    
        BasicDNAClusterModel.id_num+=1
    def _encode(self):
        raise NotImplemented()

    def _decode(self,strands):
        self._align.mpi = self.mpi
        self._cluster.mpi=self.mpi        
        return_strands=[]
        if self.mpi: #reset all strands to rank 0 in the communicator, this is inefficient, but if single thread cluster algorithms are being developed, this will be best for those
            strands=object_gather(strands,self.mpi)
        if not self.mpi or (self.mpi and not self._cluster.single_thread) or (self.mpi and self._cluster.single_thread and self.mpi.rank==0):
            clusters=self._cluster.Run(strands)  
        #distribute clusters among ranks, cluster set must be in rank 0
        if self.mpi:
            clusters=object_scatter(clusters,self.mpi,100) 
        for c in clusters:
            consensus_ID = self._consensus_id(c)
            consensus_strand = self._consensus_from_alignment(self._align.Run(c))
            new_strand = copy.copy(consensus_ID)
            new_strand.alignment_weight = len(c)
            new_strand.dna_strand = consensus_strand
            return_strands.append(new_strand)
        #strands should already be distributed based on building conesensus strands
        return return_strands

    def _consensus_id(self,cluster):
        #finds a consensus ID for the purposes of copying fault injection information
        id_counter={}
        max_ID=0
        max_strand=None
        for s in cluster:
            if not hasattr(s,"encoded_index_ints"): break
            id_counter[s.encoded_index_ints] = id_counter.get(s.encoded_index_ints,0)+1
            if id_counter[s.encoded_index_ints]>max_ID:
                max_strand=s
                max_ID=id_counter[s.encoded_index_ints]
        stats.append(self._cluster_purity_key,max_ID/len(cluster))
        stats.append(self._cluster_size_key,len(cluster))
        assert max_strand!=None
        return max_strand

    def _consensus_from_alignment(self,alignment):
        #derive a consensus strand from the aligned strands
        #assume strands are the same length
        final_strand = []
        for i in range(len(alignment[0].dna_strand)):
            majority_scores = {"A":0,"C":0,"T":0,"G":0,"-":0}
            for s in alignment:
                majority_scores[s.dna_strand[i]]+=1
            character = max(majority_scores,key=majority_scores.get)
            if not character=="-": #majority - means that there some strand might have an insertion
                final_strand.append(max(majority_scores,key=majority_scores.get))
        return "".join(final_strand)
