#!/usr/bin/python

from dnastorage.codec.base import *
from dnastorage.codec import base_conversion

import unittest
class HuffmanTableTests(unittest.TestCase):
    def test_dense_16(self):
        syms = [ x for x in range (16) ]    
        ht = HuffmanTable(2, ['0','1'], syms)
        enc,dec = ht.get_tables()
        for e in dec.keys():
            assert len(e) == 4

    def test_dense_9(self): 
        syms = [ x for x in range (9) ]    
        ht = HuffmanTable(3, ['0','1', '2'], syms)
        enc,dec = ht.get_tables()
        for e in dec.keys():
            assert len(e) == 2


class HuffmanTable:
    class Node:
        def __init__(self, nbase, symbol, weight=None, childlist=None):
            self._nbase = nbase
            self.symbol = symbol
            self.weight = weight
            self.enc = None
            self._childlist = childlist
            if weight == None and len(childlist)>0:
                self.weight = sum([ x.weight for x in childlist]) 

        def __cmp__(self, other):
            return cmp(self.weight, other.weight)
        
        def symbol(self):
            return self.symbol

        def __str__(self):
            return "{}:({},{})".format(self.enc,self.symbol,self.weight)

        def assign_enc(self, bsyms, enc):
            self.enc = enc
            if self._childlist != None:
                self._enc = enc
                for s,n in zip(bsyms,self._childlist):
                    n.assign_enc(bsyms,enc+s)
            #else:
            #    print str(self)
                

    def _build_tree(self):
        while len(self._queue) > 1:
            nodes = []
            i = 0
            while len(self._queue) > 0 and i < self._nbase:
                nodes.append( self._queue.pop(0) )
                i += 1
            new = HuffmanTable.Node(self._nbase,None,None,nodes)
            self._nodes.append(new)
            self._queue.append(new)
            # Note, this is ineffecient.  We should replace this dumb sort
            # with something more efficient like a min-heap
            self._queue.sort()

        assert len(self._queue) == 1
        r = self._queue[0]
        r.assign_enc(self._base_syms,"")

    def get_tables(self):
        table = []
        for n in self._nodes:
            if n.symbol != None:
                table.append( [n.enc, n.symbol] )

        table.sort(cmp=lambda x,y: cmp(x[1],y[1]) )

        enc = { x[1] : x[0] for x in table }
        dec = { x[0] : x[1] for x in table } 

        return (enc, dec)

    def average_length(self):
        L = 0.0
        for n in self._nodes:
            if n.symbol != None:
                L += n.weight * len(n.enc)
        return L

    def histogram(self):
        h = {}
        #H = {}
        for n in self._nodes:
            if n.symbol != None:
                l = len(n.enc)
                h[l] = h.get(l,0) + 1
                #H[l] = H.get(l,[]) + [ n.enc ]
                #H[l].sort()
        return h

    def __init__(self, nbase, base_syms, symbols, weights=None):
        """ Build a huffman encoder/decoder table under the following assumptions:    """
        """     nbase: number of symbols in codeword (4 for DNA, 2 for binary)        """
        """     base_syms: symbols used in codeword (['A','C','G','T'] or ['0', '1']) """
        """     symbols: the symbols we'll replace with huffman codewords             """
        """     weights: frequency of the symbols (same order as symbols)             """
        self._nbase = nbase
        self._base_syms = base_syms
        self._symbols = symbols
        self._weights = weights
        self._nodes = []
        if weights == None:
            self._weights = [ 1.0 / len(symbols) for _ in range(len(symbols)) ]
        else:
            W = sum(weights)
            self._weights = [ float(x) / W for x in weights ]
        assert nbase == len(base_syms)
        assert len(self._symbols) == len(self._weights)
        self._queue = []
        for s,w in zip(self._symbols,self._weights):
            n = HuffmanTable.Node(nbase,s,w)
            self._queue.append( n )
            self._nodes.append( n )
        self._queue.sort()
        self._build_tree()


if __name__ == "__main__":
    import random
    R = random.Random()
    l = [ HuffmanTable.Node(2,R.randint(0,256),R.random()) for _ in range(256) ] 
    l.sort()
    
    #for ll in l:
        #print str(ll)

    syms = [ x for x in range (256) ]
    weights = [ R.random() for _ in range(256) ]

    #ht = HuffmanTable(3, ['0','1', '2'], syms, weights)

    syms = [ 23, 14, 16, 99, 5 ]
    weights = None #[ .55, .20, .15, .10 ]

    ht2 = HuffmanTable(2, ['0', '1'], syms, weights)
    
    t = ht2.get_tables()
    print t
    print ht2.average_length()

    # The following code produces encoder table and decoder dict equivalent to 
    # the ones defined manually in huffman.py

    syms = [ x for x in range(256) ]
    w = [ 0.1 for x in range(256) ]
    for i in range(ord('A'),ord('z'),1):
        w[i] = 0.2
    w[0] = 0.2

    ht3 = HuffmanTable(3, ['0','1', '2'], syms, w)
    print ht3.get_tables()
    print ht3.average_length()
    print ht3.histogram()[0]
    print ht3.histogram()[1]

    
