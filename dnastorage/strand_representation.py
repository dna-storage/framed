'''
filename: strand_representation.py

This file implements classes that represent DNA strands that will be encoded and decoded. The purpose of enveloping DNA strands in this manner is to allow for more flexibility in how expressive the representation is. For example, the class may natively represent a DNA strand in any manner that it wants, but there will be interfaces that convert the representations appropriate for the encoding/decoding pipeline. For example, depending on whether fault injection is being run, a fault injection DNA strand can express more information about the DNA's origin.

'''


#BaseDNA class holds a basic representation of DNA, assumed to be string
class BaseDNA(object):
    def __init__(self,dna_strand = "", codewords=[],index_ints=(),before_inner=[],index_bytes=1):
        self._dna_strand=dna_strand
        self._codewords=codewords
        self._index_ints=index_ints
        self._index_bytes=index_bytes
    @property
    def dna_strand(self):
        return self._dna_strand

    @dna_strand.setter
    def dna_strand(self,DNA):
        self._dna_strand=DNA

    @property
    def codewords(self):
        return self._codewords

    @codewords.setter
    def codewords(self,c):
        self._codewords=c

    @property
    def index_ints(self):
        return self._index_ints

    @index_ints.setter
    def index_ints(self,index_ints):
        self._index_ints=index_ints
        
    @property
    def index_bytes(self):
        return self._index_bytes

    @index_bytes.setter
    def index_bytes(self,index_bytes):
        self._index_bytes=index_bytes
