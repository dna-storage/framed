#!/usr/bin/python
from dnastorage.codec.base_conversion import *
from dnastorage.primer.primer_util import *
from random import *
from copy import deepcopy


def create_simple_xor_table(size):
    table = []
    for i in range (0,size-1):
        table.append( [i] )
        table.append( [i,i+1] )

    table.append( [size-1] )
    table.append( [0,size-1] )
    return table


def create_fountain_distribution_table(size,n,limit):
    table = []
    selected = {}

    done = False
    first = True
    while True:
        r = uniform(0,1)
        if r < 0.50:
            rr = randint(0,size-1)
            if not selected.has_key(rr):
                table.append([rr])
                selected[rr] = 1
        elif r > 0.10:
            # make a large packet
            indices = {randint(0,size-1):1 for j in range(0,min(2*n,size))}
            k = indices.keys()
            k.sort()
            if len(k)>0:
                table.append(k)            
        else:
            indices = {randint(0,size-1):1 for j in range(0,max(2,randint(1,n)))}
            k = indices.keys()
            k.sort()
            if len(k)>0:
                table.append(k)

        if len(table) >= size and first:
            first = False
            missing = fountain_can_decode(table,size)
            if len(missing)>0:
                for m in missing:
                    table.append([m])
                #table.append(missing)
                assert len(fountain_can_decode(table,size))==0
        elif len(table) > 1.08*size:
            if len(fountain_can_decode(table,size))==0:
                break
                
    #table.append( [x for x in range(0,size)] )

    while len(table) < limit:
        indices = {randint(0,size-1):1 for j in range(0,randint(1,n))}
        k = indices.keys()
        k.sort()
        if len(k)>0:
            table.append(k)

    return table

def fountain_coverage(table,size):
    coverage = {x:[0,0] for x in range(0,size)}
    for entry in table:
        if len(entry)==1:
            for e in entry:
                coverage[e][0] += 1
        else:
            for e in entry:
                coverage[e][1] += 1

    print coverage

def fountain_can_decode(table,size):
    decoded = {}
    waiting = {}

    table_copy = deepcopy(table)

    for t in table_copy[:]:
        if len(t)==1:
            decoded[t[0]] = 1
            table_copy.remove(t)

    changed = True
    extra = []
    while changed:
        changed = False
        for copy in table_copy[:]:
            for tt in copy[:]:
                if decoded.has_key(tt):
                    copy.remove(tt)
            if len(copy)==1 and not decoded.has_key(copy[0]):
                decoded[copy[0]] = 1
                changed = True
                table_copy.remove(copy)
            elif len(copy)==0:
                table_copy.remove(copy)

        if changed == False and len(table_copy)>0:
            if len(table_copy[0]) > 1:
                m = table_copy[0][0]
                extra.append(m)
                decoded[m] = 1
                changed = True

    missing = [x for x in range(0,size) if not decoded.has_key(x)] + extra
    missing = {x:1 for x in missing}.keys()
    #print "Table copy: ", table_copy
    return missing
           
def evaluate_reliability(t, size, trials, erasures):
    success = {x:[0,0] for x in range(0,erasures)}
    for i in range(0,trials):
        for d in range (1,erasures):
            D = []
            table = t[:]        
            for j in range(0,d):
                r = randint(0,max(len(table)-1,1))
                if r < len(table):
                    #print "delete:",table[r]
                    D.append(table[r])
                    del table[r]
            dd = fountain_can_decode(table,size)
            if len(dd)==0:
                success[d][0]+=1
                success[d][1]+=1
            else:
                if d==1:
                    print "Deleted: ",D
                    print "Missing: ",dd
                success[d][1]+=1
    return success

class UnlimitedFountain:
    def __init__(self,numPackets,d=5):
        self._numPackets = numPackets
        self._d = d
        self.index = 0
        self._selected = {}
        self._table = []
        self._awaiting = []
        self._decoded = {}
        self._missing = { x:1 for x in range(0,numPackets) }

    def _updateDecodedAndMissing(self,strands):    
        if len(strands)==1:
            s = strands[0]
            self._decoded[s] = 1
            if self._missing.has_key(s):
                self._missing.pop(s)
            changed = True
            while changed:
                changed = False
                for a in self._awaiting[:]:
                    for i in a[:]:
                        if self._decoded.has_key(i):
                            a.remove(i)
                    if len(a) == 1:
                        if not self._decoded.has_key(a[0]):
                            self._decoded[a[0]] = 1
                            changed = True
                            self._missing.pop(a[0])
                            self._awaiting.remove(a)
                        else:
                            self._awaiting.remove(a)
                    elif len(a) == 0:
                        self._awaiting.remove(a)
        else:
            self._awaiting.append(strands[:])
        
        

    def _formNextStrand(self):
        strand = []
        r = uniform(0,1)
        if r < 0.50:
            while True:
                rr = randint(0,self.numPackets-1)
                if not self._selected.has_key(rr):
                    strand.append(rr)
                    return strand
        elif r > 0.10:
            # make a large packet
            indices = {randint(0,self.numPackets-1):1 
                       for j in range(0,min(4*self._d,self.numPackets))}
            k = indices.keys()
            k.sort()
            strand = k
        else:
            indices = {randint(0,size-1):1 for j in range(0,max(2,randint(1,n)))}
            k = indices.keys()
            k.sort()
            strand = k
        return strand

    def getTable(self):
        return self._table

    def insert(self,index,strand):
        if len(strand)==1:
            if self._selected.has_key(strand[0]):
                return
        self._table.append([index,strand])
        self._updateDecodedAndMissing(strand)

    def getMissing(self):
        return self._missing.keys()
        return fountain_can_decode(self._table,self.numPackets)

    def canDecode(self):
        return len(self._missing.keys())==0

    # Number of packets file is divided into
    @property
    def numPackets(self):
        return self._numPackets

    @property
    def numStrands(self):
        return len(self._table)

    def __iter__(self):
        self.index = 0
        self._table = []
        return self

    def next(self):
        strand = self._formNextStrand()
        return strand


if __name__ == "__main__":
    import sys
    
    #F = UnlimitedFountain(10,5)
    #for f in F:
    #    F.insert(f)
    #    print f
    #    if F.canDecode() and F.numStrands >= 15:
    #        break
    #    if F.numStrands > 10 and not F.canDecode():
    #        missing = F.getMissing()
    #        for m in missing:
    #            print [m]
    #            F.insert([m])
    #
    #print F.getTable()

    #sys.exit(0)

    average = 0
    a = []

    #for i in range(0,100):        
    #    t = create_fountain_distribution_table(10,5,15)
    #    average += len(t)
    #    a.append(len(t))

    #print a
    #print min(a), max(a)

    

    print "*"*100
    t = create_fountain_distribution_table(100,4,200)
    #print len(t),t
    success = evaluate_reliability(t,100,100,100)
    #print success
    for x,y in success.items():
        if y[1] != 0:
            print "{}\t{}".format(x,y[0]/float(y[1])) 

    print sum([y[0] for (x,y) in success.items()])/float(sum([y[1] for x,y in success.items()]))

    print "*"*100
    t = create_simple_xor_table(100)
    #print len(t),t
    success = evaluate_reliability(t,100,100,100)
    for x,y in success.items():
        if y[1] != 0:
            print "{}\t{}".format(x,y[0]/float(y[1])) 
    print sum([y[0] for (x,y) in success.items()])/float(sum([y[1] for x,y in success.items()]))
    

    #print "Making big table..."
    #t = create_fountain_distribution_table(4**7,20,4**7+4**6)
    #print len(t),len(t)/float(4**7)
