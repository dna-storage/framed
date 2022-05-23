from dnastorage.strand_representation import BaseDNA


#FaultDNA allows for a BaseDNA strand to also be tracked with information such as its encoded codewords and and faulty strand counter-part
class FaultDNA(BaseDNA):
    def __init__(self, basedna: BaseDNA,fault_strand):
        BaseDNA.__init__(self,dna_strand=fault_strand)
        #copy over apriori information from the base strand while leaving a clean slate for decoding evaluation
        self.encoded_index_ints=basedna.index_ints
        default_base_dna = BaseDNA()
        #copy over attributes embedded by probes
        for attribute in basedna.__dict__:
            if attribute not in default_base_dna.__dict__:
                setattr(self,attribute,basedna.__dict__[attribute])
                
    @property
    def encoded_index_ints(self):
        return self._encoded_index_ints

    @encoded_index_ints.setter
    def encoded_index_ints(self,ints):
        self._encoded_index_ints = ints
        
