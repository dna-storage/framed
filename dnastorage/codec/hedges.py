from bitarray import bitarray
from copy import copy
from dnastorage.codec.base_conversion import *
from dnastorage.codec.base import *
import heapq
from dnastorage.codec_types import*
from math import log10, log2, ceil, sqrt
import dnastorage.codec.fasthedges as fasthedges
from random import randint
from dnastorage.primer.primer_util import reverse_complement


try:
    from random import randbytes
except:
    def randbytes(b):
        return bytes([ randint(0,255) for _ in range(b) ])

import sys

from inspect import currentframe

from dnastorage.codec.base_conversion import *
from Bio.SeqUtils import GC

class BitTree(object):
    ''' BitTree lets us store only 1 bit in each node of the search tree. '''
    ''' We do this purely for memory efficiency. Otherwise, the search    '''
    ''' will have so many copies of the message that it will add up to    '''
    ''' many GBs on long searches.                                        '''
    def __init__(self, parent, bit):
        self.parent = parent
        self.bit = bit

    def get_prev_n(self, nbits):
        assert nbits > 0
        if self.parent != None:
            return self.parent.get_n(nbits)
        else:
            return bitarray('')

    def get_n(self, nbits):
        assert nbits > 0
        b = bitarray('')
        t = self
        while nbits > 0 and t != None:
            b = t.bit + b
            t = t.parent
            nbits -= 1
        return b
            
    def get_all(self):
        b = bitarray('')
        t = self
        while t != None:
            b = t.bit + b
            t = t.parent
        return b

class StrandTree(object):
    ''' StrandTree lets us store only 1 nt in each node of the search tree. '''    
    def __init__(self, parent, nt):
        self.parent = parent
        self.nt = nt

    def get_prev_n(self, n):
        assert n > 0
        if self.parent != None:
            return self.parent.get_n(n)
        else:
            return ''

    def get_n(self, nbits):
        assert nbits > 0
        b = ''
        t = self
        while nbits > 0 and t != None:
            b = t.nt + b
            t = t.parent
            nbits -= 1
        return b
            
    def get_all(self):
        b = ''
        t = self
        while t != None:
            b = t.nt + b
            t = t.parent
        return b



class _Heap(object):
   ''' Heap of objects with arbitrary lamda function for heapification. '''
   def __init__(self, initial=None, key=lambda x:x):
       self.key = key
       self.index = 0
       if initial:
           self._data = [(key(item), i, item) for i, item in enumerate(initial)]
           self.index = len(self._data)
           heapq.heapify(self._data)
       else:
           self._data = []

   def push(self, item):
       heapq.heappush(self._data, (self.key(item), self.index, item))
       self.index += 1

   def pop(self):
       return heapq.heappop(self._data)[2]

   def top(self):
       if len(self._data) == 0:
           return None
       return self._data[0][2]

   def empty(self):
       return len(self._data)==0

class GreedySearch:
    def __init__(self, root=None, max_guesses=100000):
        self.guesses = 0
        self.max_guesses = max_guesses
        self.root = root
        self.winner = _Heap(key=lambda x: x.score)
        self.longest = None
                                
    def search(self):
        node = self.root
        h = _Heap(key=lambda x: x.score)
        h.push(node)

        while not h.empty() and self.guesses < self.max_guesses:
            node = h.pop()
            if self.longest == None or node.index > self.longest.index:
                self.longest = node
                
            self.guesses+=1
            if node.decodeIsWinner():
                self.winner.push(node)
                break
            else:
                if node.decodeIncomplete() and not node.giveUp():
                    guesses = node.makeGuesses()                
                    for g in guesses:
                        h.push(g)
                if node != self.longest:
                    del node
                
        node = self.winner.top()
        return node
        

def ranhash(u):
    ''' Logic borrwed and adapted from Press et al., Numerical Recipes '''
    v = u * 3935559000370003845 + 2691343689449507681
    v ^= v >> 21
    v ^= v << 37
    v ^= v >> 4;
    v *= 4768777513237032717
    v ^= v << 20
    v ^= v >> 41
    v ^= v << 5
    return  v

def digest(bits, n_bits, index, index_bits, salt, n_salt, mod, verbose=False):
    b_mask = (1 << n_bits) - 1
    index_mask = (1 << index_bits) - 1
    t =  ((((index&index_mask) << n_bits) | (bits & b_mask)) << n_salt) | salt    
    t = ranhash(t) % mod 
    #t = ranhash(((((bits & b_mask) << index_bits) | index) << n_salt) | salt) % mod 
    #if verbose:
    #    print("{}: result={} bits={} nbits={} salt={} n_salt={} mod={}"\
    #          .format(index,t,bits,n_bits,salt,n_salt,mod))
    return t

DNA = { 0: 'A',
        1: 'C',
        2: 'G',
        3: 'T' }

def hasRepeat(seq):
    if seq.find("AA") >= 0:
        return True
    if seq.find("TT") >= 0:
        return True
    if seq.find("GG") >= 0:
        return True
    if seq.find("CC") >= 0:
        return True
    return False

def hasSingleRun(seq):
    if seq.find("AAA") >= 0:
        return True
    if seq.find("TTT") >= 0:
        return True
    if seq.find("GGG") >= 0:
        return True
    if seq.find("CCC") >= 0:
        return True
    return False

def hasLongRun(seq):
    if seq.find("AAAA") >= 0:
        return True
    if seq.find("TTTT") >= 0:
        return True
    if seq.find("GGGG") >= 0:
        return True
    if seq.find("CCCC") >= 0:
        return True
    return False

class Node:
    insertion_penalty = 1
    deletion_penalty = 1
    substitution_penalty = 1
    match_reward = 1
    max_penalty = max(insertion_penalty,deletion_penalty,substitution_penalty)
    good_threshold = 40

    @classmethod
    def init(cls, rate):
        Node.match_reward = cls._get_reward(rate)      

    @classmethod
    def _get_reward(cls, rate):
        if rate == 1.0:
            return -0.025
        elif rate == 0.75:
            return -0.035
        elif rate == 0.600:
            return -0.082
        elif rate == 0.500:
            return -0.127
        elif rate == 1/3:
            return -0.229
        elif rate == 0.25:
            return -0.265
        elif rate == 1/6:
            return -0.324
        elif rate == 0.125:
            return -.410
        elif rate == 0.0625:
            return -.440
        else:
            assert False
        
    def __init__(self, hedge, observed, offset, index, seqnum_bt, message_bt, corrected_bt, score, fails):
        self.hedge = hedge 
        self.observed = observed
        self.offset = offset
        self.index = index
        # These are trees that hold only one bit
        self.seqnum_bt = seqnum_bt
        self.message_bt = message_bt
        self.corrected_bt = corrected_bt
        # self.seqnum_vector = seqnum_vector
        # self.message_vector = message_vector
        # self.corrected = copy(corrected)
        self.message_vector = None
        self.seqnum_vector = None
        self.corrected = None        
        self.score = score
        self.fails = fails

    def get_message_vector(self):
        if self.message_vector == None:
            self.message_vector = self.message_bt.get_all()
        return self.message_vector
    
    def get_seqnum_vector(self):
        if self.seqnum_vector == None:
            self.seqnum_vector = self.seqnum_bt.get_all()
        return self.seqnum_vector

    def get_corrected(self):
        if self.corrected == None:
            self.corrected = self.corrected_bt.get_all()
        return self.corrected

    def setup(self):
        self.get_message_vector()
        self.get_seqnum_vector()
        self.get_corrected()
        
    def cleanup(self):
        del self.corrected
        del self.seqnum_vector
        del self.message_vector
        self.corrected = None
        self.seqnum_vector = None
        self.message_vector = None

    def atObserved(self, at):
        if at >= len(self.observed):
            return 'A'
        else:
            return self.observed[at]        
    
    def giveUp(self):
        return self.offset >= len(self.observed)
        
    def decodeIncomplete(self):
        return self.index < self.hedge.index_range  

    def checkPad(self):
        if self.hedge.pad_bits == 0:
            return True
        bits,i = self.hedge._get_nbits(self.index)        
        pad = self.hedge.pad_bits + self.hedge.pad_message
        # make message vector for checking
        message_vector = self.get_message_vector()
        pad0 = self.hedge._get_prev_bits(len(message_vector),message_vector,pad)
        #print (self.message_vector)
        return pad0 == 0
    
    def decodeIsWinner(self):
        winner = not self.decodeIncomplete() and self.checkPad()
        return winner

    def keepGoing(self):
        # TODO check pad bits
        return self.score > Node.good_threshold

    def guessBit(self,val):              
        if val==0:
            b = bitarray('0')
        else:
            b = bitarray('1')

        if self.hedge._inSeqNumber(self.index):
            s = self.hedge.get_symbol(self.index, self.seqnum_vector+b, self.message_vector, self.corrected)
        else:
            s = self.hedge.get_symbol(self.index, self.seqnum_vector, self.message_vector+b, self.corrected)

        return s

    def guess2Bit(self,val):              
        if val==0:
            b = bitarray('00')
        elif val==1:
            b = bitarray('01')
        elif val==2:
            b = bitarray('10')            
        else:
            b = bitarray('11')

        if self.hedge._inSeqNumber(self.index):
            s = self.hedge.get_symbol(self.index, self.seqnum_vector+b, self.message_vector, self.corrected)
        else:
            s = self.hedge.get_symbol(self.index, self.seqnum_vector, self.message_vector+b, self.corrected)

        return s
    

    def _addMatch(self, s, b):
        if self.hedge._inSeqNumber(self.index):
            match = Node(self.hedge,self.observed,self.offset+1,self.index+1,\
                         BitTree(self.seqnum_bt,b), self.message_bt, StrandTree(self.corrected_bt,s),\
                         self.score+Node.match_reward, self.fails)
        else:
            match = Node(self.hedge,self.observed,self.offset+1,self.index+1,\
                         self.seqnum_bt, BitTree(self.message_bt,b), StrandTree(self.corrected_bt,s),\
                         self.score+Node.match_reward, self.fails)
        return [match]               

    def _addSubst(self, s, b):
        if self.hedge._inSeqNumber(self.index):
            sub = Node(self.hedge,self.observed,self.offset+1,self.index+1,\
                         BitTree(self.seqnum_bt,b), self.message_bt, StrandTree(self.corrected_bt,s),\
                         self.score+Node.substitution_penalty, self.fails+1)
        else:
            sub = Node(self.hedge,self.observed,self.offset+1,self.index+1,\
                         self.seqnum_bt, BitTree(self.message_bt,b), StrandTree(self.corrected_bt,s),\
                         self.score+Node.substitution_penalty, self.fails+1)
        return [sub]

    def _addDel(self, s, b):
        if self.hedge._inSeqNumber(self.index):
            de = Node(self.hedge,self.observed,self.offset,self.index+1,\
                      BitTree(self.seqnum_bt,b), self.message_bt, StrandTree(self.corrected_bt,s),\
                      self.score + Node.deletion_penalty, self.fails+1)
        else:
            de = Node(self.hedge,self.observed,self.offset,self.index+1,\
                      self.seqnum_bt, BitTree(self.message_bt,b), StrandTree(self.corrected_bt,s),\
                      self.score + Node.deletion_penalty, self.fails+1)
        return [de]

    def _addIns(self):
        insert = Node(self.hedge,self.observed,self.offset+1,self.index,\
                      self.seqnum_bt, self.message_bt, self.corrected_bt,\
                      self.score + Node.insertion_penalty, self.fails+1)                        
        return [insert]

    def _addIns2(self, s, b):
        if self.atObserved(self.offset+1) == s:
            if self.hedge._inSeqNumber(self.index):
                insert = Node(self.hedge,self.observed,self.offset+2,self.index+1,\
                              BitTree(self.seqnum_bt,b), self.message_bt, StrandTree(self.corrected_bt,s),\
                              self.score+Node.insertion_penalty+Node.match_reward, self.fails+1)
            else:
                insert = Node(self.hedge,self.observed,self.offset+2,self.index+1,\
                              self.seqnum_bt, BitTree(self.message_bt,b), StrandTree(self.corrected_bt,s),\
                              self.score+Node.insertion_penalty+Node.match_reward, self.fails+1)
        else:
            if self.hedge._inSeqNumber(self.index):
                insert = Node(self.hedge,self.observed,self.offset+2,self.index+1,\
                              BitTree(self.seqnum_bt,b), self.message_bt, StrandTree(self.corrected_bt,s),\
                              self.score+Node.substitution_penalty+Node.insertion_penalty, self.fails+1)
            else:
                insert = Node(self.hedge,self.observed,self.offset+2,self.index+1,\
                              self.seqnum_bt, BitTree(self.message_bt,b), StrandTree(self.corrected_bt,s),\
                              self.score+Node.substitution_penalty+Node.insertion_penalty, self.fails+1)
            
        return [insert]
        

    def makeGuesses(self):
        self.setup()
        
        bits,i = self.hedge._get_nbits(self.index)
        if bits == 0:
            assert self.hedge.rate < 0.5
            g = self.make0bitGuesses()
        elif bits == 1:
            g = self.make1bitGuesses()
        else:
            assert bits == 2
            g = self.make2bitGuesses()

        self.cleanup()
        return g
            
    def make0bitGuesses(self):
        # Generate all possible guesses from this index assuming we know
        # it must be a 0.
        
        # matches get the best score, but it's possible a "match"
        # is erroneous, so also guess substitution, deletion, and
        # insertion. If match is ruled out, the others will be
        # searched. If match is correct, the others will remain low
        # priority in the queue and likely will never be searched.
        b0 = bitarray('') # do not add a bit to the message
        s_0 = self.guessBit(0)

        children = []
        
        if s_0 == self.atObserved(self.offset):
            # Actual match
            children += self._addMatch(s_0,b0)
        else:
            children += self._addSubst(s_0,b0)

        children += self._addDel(s_0,b0)
        children += self._addIns2(s_0, b0)
        
        # if self.atObserved(self.offset+1) == s_0:
        # else:
        #     children += self._addIns()
                
        return children
        
    def make1bitGuesses(self):
        # Generate all possible guesses from this index. We do it this way
        # so that we never need to visit this node again once it's popped
        # from the heap. But, this comes at the expense of using memory on nodes
        # we will likely never explore.
        
        # matches get the best score, but it's possible a "match"
        # is erroneous, so also guess substitution, deletion, and
        # insertion. If match is ruled out, the others will be
        # searched. If match is correct, the others will remain low
        # priority in the queue and likely will never be searched.
        b0 = bitarray('0')
        b1 = bitarray('1')
        s_0 = self.guessBit(0)
        s_1 = self.guessBit(1)

        children = []
        
        if s_0 == self.atObserved(self.offset):
            # Actual match
            children += self._addMatch(s_0,b0)
        else:
            children += self._addSubst(s_0,b0)

        children += self._addDel(s_0,b0)

        if s_1 == self.atObserved(self.offset):
            children += self._addMatch(s_1,b1)
        else:
            children += self._addSubst(s_1,b1)

        children += self._addDel(s_1,b1)

        children += self._addIns2(s_0, b0)
        children += self._addIns2(s_1, b1)
        
        # For insertions, take into account the next position. If they match, it increases
        # the odds that the insertion actually occured here.
        # if self.offset + 1 < len(self.observed) and (self.atObserved(self.offset+1) == s_0 or \
        #                                              self.atObserved(self.offset+1) == s_1):
        #     if self.atObserved(self.offset+1) == s_0:                
        #         children += self._addIns2(s_0, b0)
        #     else:
        #         children += self._addIns2(s_1, b1)
        # else:
        #     children += self._addIns()
                
        return children

    def make2bitGuesses(self):
        # Generate all possible guesses from this index. We do it this way
        # so that we never need to visit this node again once it's popped
        # from the heap. But, this comes at the expense of using memory on nodes
        # we will likely never explore.
        
        # matches get the best score, but it's possible a "match"
        # is erroneous, so also guess substitution, deletion, and
        # insertion. If match is ruled out, the others will be
        # searched. If match is correct, the others will remain low
        # priority in the queue and likely will never be searched.
        b00 = bitarray('00')
        b01 = bitarray('01')
        b10 = bitarray('10')
        b11 = bitarray('11')
        s_00 = self.guess2Bit(0)
        s_01 = self.guess2Bit(1)
        s_10 = self.guess2Bit(2)
        s_11 = self.guess2Bit(3)

        children = []
        
        children += self._addIns()
        
        if s_00 == self.atObserved(self.offset):
            # Actual match
            children += self._addMatch(s_00,b00)
        else:
            children += self._addSubst(s_00,b00)
            children += self._addDel(s_00,b00)
                        
        if s_01 == self.atObserved(self.offset):
            children += self._addMatch(s_01,b01)
        else:
            children += self._addSubst(s_01,b01)
            children += self._addDel(s_01,b01)

        if s_10 == self.atObserved(self.offset):
            children += self._addMatch(s_10,b10)
        else:
            children += self._addSubst(s_10,b10)
            children += self._addDel(s_10,b10)

        if s_11 == self.atObserved(self.offset):
            children += self._addMatch(s_11,b11)
        else:
            children += self._addSubst(s_11,b11)
            children += self._addDel(s_11,b11)
            
        return children
        
class HEDGE:
    '''
    HEDGE is a class that performs encoding or decoding over a fixed set of parameters:
     1. rate = codec rate [1.0, .75, 0.5, 1/3, 0.25, 1/6, 1/8, 1/16]         
     2. index_bits = number of bits devoted to indexing the strand.          
     3. message_bits = length of payload in bits.                            
     4. pad_bits = extra 0 bits at end.
     5. prev_bits = # of previous message bits used in hash calculation. 
     6. salt_bits = # of previous index bits used in hash calculation.
    '''
    rates = [1.0, .75, 0.5, 1/3, 0.25, 1/6, 1/8, 1/16]
    
    def __init__(self, rate, pad_bits, prev_bits, seqnum_bits=8, message_bits=80):
        self.rate = rate
        assert self.rate in HEDGE.rates       
        self.seqnum_bits = seqnum_bits
        self.message_bits = message_bits
        self.pad_bits = pad_bits
        self.prev_bits = prev_bits
        self.salt_bits = seqnum_bits
        self._bits_hash = {}
        self._set_pad()


    def set_bit_sizes(self,seqnum_bits,message_bits): #interface to allow for easing setting of seqnum_bits
        self.seqnum_bits = seqnum_bits
        if seqnum_bits < 32:
            self.salt_bits = seqnum_bits
        else:
            # unlikely to have more than 4 bytes of index, but even if we do, just use the lower
            # 32 bits. We can revisit this later.
            self.salt_bits = 32 
        self.message_bits = message_bits
        self._set_pad()

    def _set_pad(self):
        self.pad_seqnum = self._check_if_padding_needed(self.rate,self.seqnum_bits)
        self.pad_message = self._check_if_padding_needed(self.rate,self.seqnum_bits+\
                                                         self.pad_seqnum+self.message_bits+self.pad_bits)       
        self.adj_seqnum_bits = self.seqnum_bits + self.pad_seqnum
        self.adj_message_bits = self.message_bits + self.pad_bits + self.pad_message        
        self.index_bits = int(ceil(log2((self.adj_seqnum_bits + self.adj_message_bits) / 4 / rate)))
        self.index_range = self._get_index_range(self.rate, self.adj_seqnum_bits + self.adj_message_bits)

        if self.adj_seqnum_bits<self.salt_bits:
            self.salt_bits=self.adj_seqnum_bits #KV: upper bound off the salt bits

        assert self.salt_bits + self.index_bits + self.prev_bits <= 64
            
    def _get_pattern(self, rate):
        if rate == 1.0:
            return [2]
        elif rate == 0.75:
            return [2,1]
        elif rate == 0.5:
            return [1]
        elif rate == 1/3:
            return [1,1,0]
        elif rate == 0.25:
            return [1,0]
        elif rate == 1/6:
            return [1,0,0]
        elif rate == 1/8:
            return [1,0,0,0]
        elif rate == 1/16:
            return [1,0,0,0,0,0,0,0]
        else:
            assert False and "No such rate."

    def _insert_zeroes(self, vector):
        if self.rate >= 0.5:
            return vector
        pattern = self._get_pattern(self.rate)
        new_vector = bitarray('')
        i = 0
        j = 0
        while i < len(vector):
            if pattern[j] == 1:
                new_vector += vector[i:i+1]
                i += 1
            elif pattern[j] == 2:
                new_vector += vector[i:i+2]
                i += 2
            elif pattern[j] == 0:
                new_vector += bitarray('0')
            else:
                assert False and "Illegal case detected."
                
            j = (j+1) % len(pattern)
        return new_vector
    
    def _get_index_range(self, rate, bits):
        pat = self._get_pattern(rate)
        i = 0
        j = 0
        while bits > 0:
            bits -= pat[j % len(pat)]
            i += 1
            j += 1
        if j % len(pat) != 0:
            i += len(pat) - j % len(pat)
        assert bits == 0
        return i
            
    def _check_if_padding_needed(self, rate, bits):
        pat = self._get_pattern(rate)
        s = sum(pat)
        if bits % s == 0:
            return 0
        else:
            return s - bits%s
        #i = 0
        #while bits > 0:
        #    bits -= pat[i]
        #    i = (i+1) % len(pat)
        #return - bits
            
    def _get_prev_bits(self, i, bmessage, prev_bits):
        if prev_bits == 0:
            return 0
        # indexing here is inclusive
        assert i >= prev_bits-1
        # i need the value to be in the low bits, so if prev_bits < 8
        # add a buffer of 0s to make sure value stays in the low bits
        buf = 0    
        if prev_bits % 8 > 0:
            buf = 8 - prev_bits % 8
        prev = bitarray("0"*(buf))+bmessage[i-prev_bits:i]
        return convertBytesToInt(prev.tobytes())

    def _analyze_strand_constraints(self, strand, n_bits):
        ''' 
            This function assumes that the strand is partially filled up to the current
            index that's being encoded/decoded. So, it only needs to consider
            the end of the DNA string to determine the coding constraints.
        '''
        DNA = [ 'A', 'C', 'G', 'T' ]
        if len(strand) < 4:
            return DNA

        if n_bits <= 1:
            '''
               For the case of n_bits = 1, we can remove up to 2 bases without loss
               of information.
            '''
            if hasSingleRun(strand[-3:]) or hasRepeat(strand[-2:]):
                DNA.remove(strand[-1])
                return DNA
            elif GC(strand) > 75:
                DNA.remove('C')
                DNA.remove('G')
            elif GC(strand) < 25:
                DNA.remove('A')
                DNA.remove('T')
            return DNA
        else:
            '''
               For the case of n_bits = 2, we cannot any bases.
            '''            
            assert n_bits==2
            return DNA

    def _get_nbits_helper(self, index):
        if index in self._bits_hash:
            return self._bits_hash[index]
        pattern = self._get_pattern(self.rate)
        t = 0
        b = 0
        while t < index:
            b += pattern[t%len(pattern)]
            t += 1
        self._bits_hash[index] = (pattern[index%len(pattern)],b)
        return (pattern[index%len(pattern)],b)
                
    def _get_nbits(self, index):
        if self.rate < 0.5:
            return self._get_nbits_helper(index)
        if self.rate == 0.5:
            # for rates less than or equal to 0.5, we insert 0s into the message for simplicity
            # but the decoder knows they have to be 0s
            return 1,index
        elif self.rate == 0.75:
            pattern = [2,1]
            return pattern[index%len(pattern)], \
                sum(pattern)*index//len(pattern) + sum(pattern[0:(index%len(pattern))])     
        elif self.rate == 1.0:
            return 2,index*2 


    def _inSeqNumber(self, index):
        bits,i = self._get_nbits(index)                
        return i < self.adj_seqnum_bits
        
    def _get_salt_for_index(self, index, seqnum_vector):
        if self._inSeqNumber(index):
            return 0
        assert self.salt_bits <= len(seqnum_vector)
        pbits =  self._get_prev_bits(len(seqnum_vector),seqnum_vector,self.salt_bits)
        #print (pbits)
        #print (seqnum_vector)
        return pbits
        
    def _get_prev_message_bits_for_index(self, index, message_vector):
        if self._inSeqNumber(index):
            return 0
        bits,i = self._get_nbits(index)
        i = i - self.adj_seqnum_bits
        prev_bits = min(i, self.prev_bits)
        prev = self._get_prev_bits(i, message_vector, prev_bits)
        return prev

        
    def get_symbol(self, index, seqno_vector, message_vector, strand, verbose=False):
        ''' Determine next symbol of the strand. This function does not produce 
            side-effects so that it can be used in both encoding and decoding. '''
        bits,i = self._get_nbits(index)
        if self._inSeqNumber(index):
            sym = self._get_prev_bits(i+bits,seqno_vector,bits)
            prev = self._get_prev_bits(i,seqno_vector,min(i,self.prev_bits))
            salt = 0            
        else:            
            sym = self._get_prev_bits(i-self.adj_seqnum_bits+bits,message_vector,bits)
            prev = self._get_prev_message_bits_for_index(index, message_vector)
            salt = self._get_salt_for_index(index, seqno_vector)
            
        DNA = self._analyze_strand_constraints(strand,bits)
        res = digest(prev, self.prev_bits, index, self.index_bits, salt, self.salt_bits, len(DNA), verbose=verbose)
        res = (res+sym) % len(DNA)
        return DNA[res]    

    def encode(self, seqnum_vector, message_vector):
        """ seqnum_vector:  bitarray that matches index_bits in length. """
        """ message_vector: bytes array truncated to match self.message_bits """
        
        message_vector_adj = message_vector + bitarray("0"*(self.pad_bits+self.pad_message))
        seqnum_vector_adj = seqnum_vector + bitarray("0"*self.pad_seqnum)

        assert  self.adj_message_bits == len(message_vector_adj)
        assert  self.adj_seqnum_bits == len(seqnum_vector_adj)

        strand = ""
        
        for index in range(0,self.index_range):
            s = self.get_symbol(index,seqnum_vector_adj, message_vector_adj, strand)
            strand += s
            
        return strand

    def decode(self, observed, best_guess=False, max_guesses=100000):
        strand_pos = 0
        index = 0        
        Node.init(self.rate)
        node = Node(self, observed, strand_pos, index, BitTree(None,bitarray("")), BitTree(None,bitarray("")), StrandTree(None,""), self.index_range, 0)
        gs = GreedySearch(node,max_guesses=max_guesses)
        node = gs.search()
        if node == None and best_guess==True:
            node = gs.longest
        
        if node != None and node.decodeIsWinner():
            padding = self.pad_message + self.pad_bits
            if padding == 0:
                return node.get_seqnum_vector(), node.get_message_vector(), node.get_corrected()                
            else:                
                return node.get_seqnum_vector(), node.get_message_vector()[0:-padding], node.get_corrected()
        elif node != None:
            return node.get_seqnum_vector(), node.get_message_vector(), node.get_corrected()                            
        else:
            return None,None,None


class PyHedgesPipeline(BaseCodec,CWtoDNA):
    def __init__(self,rate,pad_bits,prev_bits,guess_limit=100000,CodecObj=None,Policy=None):
        self._hedges = HEDGE(rate,pad_bits,prev_bits)
        self._guess_limit=guess_limit
        CWtoDNA.__init__(self)
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)

    def _encode(self,strand):
        data_bits = 8*len(strand.codewords[strand.index_bytes:])
        index_bytes = strand.codewords[:strand.index_bytes]
        index_bits = 8*len(index_bytes)
        #convert codewords to bits
        data_bit_array = bitarray(buffer=bytearray(strand.codewords[strand.index_bytes:]))
        index_bit_array=bitarray(buffer=bytearray(index_bytes))
        
        self._hedges.set_bit_sizes(index_bits,data_bits)
        strand.dna_strand = self._hedges.encode(index_bit_array,data_bit_array)
        return strand
        
    def _decode(self,strand):
        index,data,corr = self._hedges.decode(strand.dna_strand,best_guess=True,max_guesses=self._guess_limit)
        strand.codewords =[_ for _ in index.tobytes()+data.tobytes()]
        return strand

    #store some pertinent information like bit lengths of data seen to be able to reinstantiate the decoder in a correct state    
    def _encode_header(self):
        data = []
        data+=convertInttoBytes(self._hedges.seqnum_bits,4)
        data+=convertIntoBytes(self._hedges.message_bits,4)
        return data
    def _decode_header(self,buff):
        pos=0
        seqnum_bits = convertBytestoInt(buff[pos:pos+4])
        pos+=4
        message_bits = convertBytestoInt(buff[pos:pos+4])
        pos+=4
        self._hedges.set_bit_sizes(seqnum_bits,message_bits)
        return buff[pos:]




class hedges_state:
    def __init__(self, rate=1.0/4, seq_bytes=2, message_bytes=14, pad_bits=8, prev_bits = 8,sync_period=0,parity_period=0,parity_history=0):
        self.rate = float(rate)
        self.seq_bytes = seq_bytes
        self.message_bytes = message_bytes
        self.pad_bits = pad_bits #this is user-defined pad bits, not those determined under the hood to make patterns work out
        self.prev_bits = prev_bits
        self.salt_bits = seq_bytes * 8
        self.parity_period=parity_period #period in which to consider built in parity
        self.cw_sync_period=sync_period #sync period is the period in which we put syncpoints
        if self.salt_bits > 32:
            self.salt_bits = 32
        self.parity_history = parity_history
        if self.parity_history>64:
            self.parity_history=64
        
    def set_message_bytes(self,message_bytes):
        self.message_bytes = message_bytes

    def set_seqnum_bytes(self,seq_bytes):
        self.seq_bytes = seq_bytes
        self.salt_bits = seq_bytes * 8
        if self.salt_bits > 32:
            self.salt_bits = 32
    
class FastHedgesPipeline(BaseCodec,CWtoDNA):
    def __init__(self,rate,pad_bits=8,prev_bits=8,guess_limit=100000,CodecObj=None,Policy=None,try_reverse = False):
        self._hedges_state = hedges_state(rate=rate,pad_bits=pad_bits,prev_bits=prev_bits)
        self._guess_limit=guess_limit
        self._try_reverse=try_reverse #option to try reverse complement, should do this if DNA not guarenteed to be in right position
        CWtoDNA.__init__(self)
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)

    def _encode(self,strand):
        self._hedges_state.set_message_bytes(len(strand.codewords)-strand.index_bytes)
        self._hedges_state.set_seqnum_bytes(strand.index_bytes)
        strand.dna_strand = fasthedges.encode(bytes(strand.codewords),self._hedges_state)
        return strand
        
    def _decode(self,strand):
        reverse = False
        if self._try_reverse: #if we need to check reverse, try a small number of guesses
            reverse_codewords = fasthedges.decode(reverse_complement(strand.dna_strand), self._hedges_state, 1000)
            forward_codewords = fasthedges.decode(strand.dna_strand, self._hedges_state, 1000)
            reverse_none = sum([1 if _==None else 0 for _ in reverse_codewords])
            forward_none = sum([1 if _==None else 0 for _ in forward_codewords])
            if reverse_none<forward_none: #utilize reverse complement
                reverse=True
        #print (self._hedges_state.seq_bytes, self._hedges_state.message_bytes)
        if not reverse: strand.codewords = fasthedges.decode(strand.dna_strand, self._hedges_state, self._guess_limit)
        else: strand.codewords = fasthedges.decode(reverse_complement(strand.dna_strand), self._hedges_state, self._guess_limit)
        return strand
    
    #store some pertinent information like bit lengths of data seen to be able to reinstantiate the decoder in a correct state    
    def _encode_header(self):
        data = []
        data+=convertIntToBytes(self._hedges_state.seq_bytes,4)
        data+=convertIntToBytes(self._hedges_state.message_bytes,4)
        return data
    
    def _decode_header(self,buff):
        pos=0
        seqnum_bytes = convertBytesToInt(buff[pos:pos+4])
        pos+=4
        message_bytes = convertBytesToInt(buff[pos:pos+4])
        pos+=4
        self._hedges_state.set_message_bytes(message_bytes)
        self._hedges_state.set_seqnum_bytes(seqnum_bytes)
        return buff[pos:]

    
def inject(strand,rate):
    strand_list = [ _ for _ in strand ]
    output = []
    strand_len = len(strand)

    faults = 0
    i = 0
    while i < len(strand):
        s = strand_list[i]
        DNA = ['A', 'C', 'G', 'T']
        r = randint(0,1000)
        #print (r/1000)
        if r/1000 < rate:
            c = randint(0,2)            
            if c==0: # substitute
                DNA.remove(s)
                output.append(DNA[randint(0,len(DNA)-1)])
                i += 1
                faults += 1
                pass
            elif c==1: # insert
                output.append(DNA[randint(0,len(DNA)-1)])
                #output.append(s)
                #i += 1
                faults += 1                
                pass
            elif c==2:
                # delete
                i += 1
                faults += 1                
                pass
        else:
            output.append(s)
            i += 1

    #print ("Inserted {} faults in {} length strand.".format(faults,strand_len))
    return "".join(output)



def diff(A,B,leading=""):
    diff = [ ' ' if a==b else "|" for a,b in zip(A,B) ]
    print ("{}{}".format(leading,A))
    print ("{}{}".format(leading,"".join(diff)))
    print ("{}{}".format(leading,B))


def run_hedges(rate,fi_rate=0.01,length=10):
    #h = HEDGE(rate, 8, 16, 8, seqnum_bits= 8, message_bits=length*8)
    h = FastHedgesPipeline(rate)
    rb = buffer=randbytes(length)
    #print([_ for _ in rb])
    s = BaseDNA(codewords=[_ for _ in rb])
    s.index_bytes = 4
    h.encode(s)
    s.dna_strand = inject(s.dna_strand,fi_rate)
    header = h.encode_header()
    h = FastHedgesPipeline(rate)
    h.decode_header(header)
    h.decode(s)
    #
    if s.codewords==[_ for _ in rb]:
        return True
    else:
        print ([_ for _ in rb])
        print([_ for _ in s.codewords])        
        print (s.dna_strand)
        return False

    
if __name__ == "__main__":
    import sys
    import dnastorage.strand_representation    
    
    if len(sys.argv) >= 4:
        rate = sys.argv[1]
        try:
            rate = float(sys.argv[1])
        except:
            if rate == "3/4":
                rate = 3/4
            elif rate == "1/3":
                rate = 1/3
            elif rate == "1/4":
                rate = 1/4
            elif rate == "1/6":
                rate = 1/6
            elif rate == "1/8":
                rate = 1/8
            elif rate == "1/16":
                rate = 1/16
            else:
                print ("assuming rate = 1/2")
                rate = 1/2                
        fi = float(sys.argv[2])                
        length = int(sys.argv[3])
    else:
        rate = 0.5
        fi = 0.01
        length = 10

    matched = 0
    trials = 1000
    total = 0
    while trials > 0:
        if run_hedges(rate,fi_rate=fi,length=length):
            matched += 1
        trials -= 1
        total += 1
        #print (matched, total)

    print ("rate={} fi={} {} bytes {}%".format(rate,fi, length, matched/total * 100)) 
