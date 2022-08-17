from dnastorage.codec.base import *
from dnastorage.primer.primer_util import *
from dnastorage.exceptions import *
from dnastorage.codec_types import *
from dnastorage.strand_representation import *
from Bio import pairwise2

import logging
logger = logging.getLogger('dna.storage.system.dnafile')
logger.addHandler(logging.NullHandler())


class CombineCodewords(BaseCodec):
    def __init__(self,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)

    def _encode(self, codeword_list):
        return "".join(codeword_list)

    def _decode(self, s):
        assert ("not used for decoding")




class NormalizeStrandLength(BaseCodec):
    def __init__(self,length,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self.length = length
        
    def _encode(self, phys_s):
        if len(phys_s) > self.length:
            e = DNAStrandPayloadWrongSize("NormalizeStrandLength: Strand is too long ({})".format(len(phys_s)))
            if self._Policy.allow(e):
                pass
            else:
                raise e
        elif len(phys_s) < self.length:
            add = self.length - len(phys_s)
            phys_s = phys_s + "".join([ random.choice('AGCT') for _ in range(add) ])
            
        return phys_s

    def _decode(self, s):
        assert ("not used for decoding")

        
        
class InsertMidSequence(BaseCodec):
    def __init__(self,seq,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._seq = seq

    def _encode(self,strand):
        if strand.find(self._seq) != -1:
            err = DNAStrandPoorlyFormed("Found sequence already present while prepending {}"\
                                        .format(self._seq))
            if self._Policy.allow(err):
                pass
            else:
                raise err

        middle = int(len(strand)/2)
        return strand[0:middle]+self._seq+strand[middle:]

    def _decode(self,strand):
        index = strand.find(self._seq)
        if index != -1:
            return strand[0:index]+strand[index+len(self._seq):]
        else:
            err = DNAStrandMissingSequence("{} should have had {} inside it.".format(strand,self._seq))
            # there could be errors in the cut preventing us from seeing it
            # with an exact match, so now we look for an inexact match
            if self._Policy.allow(err):
                middle = len(strand)//2
                slen = len(self._seq)
                res = []
                for m in range(middle-slen,middle):
                    sli = strand[m:m+slen]
                    res.append( ed.eval(sli,self._seq) )
                mn = min(res)
                if mn < slen//3:
                    place = middle - slen + res.index(mn)
                    return strand[0:place]+strand[place+slen:]
                else:
                    # just leave the strand along, and hopefully codewords can
                    # still be extracted properly
                    return strand
            else:
                raise err



def find_ed(strand,seq,distance):
    results=[]
    for strand_index, base in enumerate(strand):
        check = strand[strand_index:len(seq)+strand_index]
        results.append(ed.eval(check,seq))
    if len(results)==0: return None
    m = min(results)
    if m <=distance: return results.index(m)
    else: return None
    
            
            

class PrependSequence(BaseCodec):
    def __init__(self,seq,CodecObj=None,Policy=None,isPrimer=False,handler="ed",search_range=50):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._seq = seq[:]
        self.is_primer = isPrimer
        self._handler=handler
        self._search_range=search_range
    def _encode(self,strand):
        if strand.find(self._seq) != -1:
            err = DNAStrandPoorlyFormed("Found sequence already present while prepending {}"\
                                        .format(self._seq))
            if self._Policy.allow(err):
                pass
            else:
                raise err
        return self._seq + strand

    def _decode(self,strand):
        slen = len(self._seq)
        index = strand.find(self._seq)
        if index != -1: # expected at beginning
            return strand[index+len(self._seq):]
        elif self._handler=="ed":
            err = DNAStrandMissingSequence("{} should have had {} at beginning.".format(strand,self._seq))
            # there could be errors in the cut preventing us from seeing it
            # with an exact match, so now we look for an inexact match
            if self._Policy.allow(err):
                idx = find_ed(strand[0:self._search_range+slen],self._seq,5)
                if idx!=None:
                    return strand[idx+slen:]
                else:
                    if self.is_primer:
                        raise DNAMissingPrimer("Missing primer {}".format(self._seq))
                    # just leave the strand along, and hopefully codewords can
                    # still be extracted properly
                    return strand 
            else:
                raise err
        elif self._handler=="align":
            align = pairwise2.align.localms(strand[0:self._search_range+slen],self._seq,1,-1,-1,-1,one_alignment_only=True)
            if len(align)==0: return strand
            align=align[0]
            score= align[2]
            score = sum([1 if _==_2 else 0 for _,_2 in zip(align[1],align[0])])
            if score<(len(self._seq) - (len(self._seq)*0.3)): #KV made change for score to count total matches after realigning 
                #logger.info("Strand ID {} \n Searched Sequence: {} \n Desired Sequence: {}".format(id(strand),strand[0:self._search_range+slen],self._seq))
                return strand #finding alignment unsuccessful
            else:
                return strand[align[4]:]



class AppendSequence(BaseCodec):
    def __init__(self,seq,CodecObj=None,Policy=None,isPrimer=False,handler="ed",search_range=50):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._seq = seq
        self.is_primer = isPrimer
        self._handler=handler
        self._search_range=search_range
    def _encode(self,strand):
        if strand.find(self._seq) != -1:
            err = DNAStrandPoorlyFormed("Found sequence already present while appending {}"\
                                        .format(self._seq))
            if self._Policy.allow(err):
                pass
            else:
                raise err
        
        return strand + self._seq

    def _decode(self,strand):
        index = strand.find(self._seq)
        slen = len(self._seq)
        if index != -1: # expected at end
            return strand[:index]
        elif self._handler=="ed":
            err = DNAStrandMissingSequence("{} should have had {} at end.".format(strand,self._seq))
            # there could be errors in the cut preventing us from seeing it
            # with an exact match, so now we look for an inexact match
            if self._Policy.allow(err):
                idx = find_ed(strand[(len(strand)-self._search_range):len(strand)],self._seq,5)
                if idx!=None:
                    idx+=(len(strand)-self._search_range)
                    return strand[:idx+len(strand)-self._search_range]
                else:
                    if self.is_primer:
                        raise DNAMissingPrimer("Missing primer".format(self._seq))
                    # just leave the strand along, and hopefully codewords can
                    # still be extracted properly
                    return strand
            else:
                raise err
        elif self._handler=="align":
            align = pairwise2.align.localms(strand[len(strand)-len(self._seq)-self._search_range:len(strand)],self._seq,1,-1,-1,-1,one_alignment_only=True)
            if len(align)==0:return strand
            align=align[0]
            score= align[2]
            score = sum([1 if _==_2 else 0 for _,_2 in zip(align[1],align[0])])
            if score<(len(self._seq) - (len(self._seq)*0.3)):
                #logger.info("Strand ID {} \n Searched Sequence: {} \n Desired Sequence: {}".format(id(strand),strand[len(strand)-len(self._seq)-self._search_range:len(strand)],self._seq))
                return strand
            else:
                return strand[0:align[3]+len(strand)-len(self._seq)-self._search_range]

class PrependSequencePipeline(PrependSequence,DNAtoDNA):
    def __init__(self,seq,CodecObj=None,Policy=None,isPrimer=False,ignore=False,handler="ed",search_range=50):
        PrependSequence.__init__(self,seq,CodecObj=CodecObj,Policy=Policy,isPrimer=isPrimer,handler=handler,search_range=search_range)
        DNAtoDNA.__init__(self)
        self._ignore=ignore
    def _encode(self,strand):
        #this is a wrapper around the basic PrependSequence _encode so that the strand interface
        strand.dna_strand=PrependSequence._encode(self,strand.dna_strand)
        return strand
    def _decode(self,strand):
        if strand.dna_strand is None or self._ignore or self._seq=="": return strand
        strand_before =strand.dna_strand
        strand.dna_strand=PrependSequence._decode(self,strand.dna_strand)
        if strand_before==strand.dna_strand and self._seq != "":
            strand.dna_strand = None
        return strand


class ReversePipeline(BaseCodec,DNAtoDNA):
    def __init__(self,CodecObj=None,Policy=None):
        DNAtoDNA.__init__(self)
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
    def _encode(self,strand):
        #this is a wrapper around the basic PrependSequence _encode so that the strand interface
        strand.dna_strand=strand.dna_strand[::-1]
        return strand
    def _decode(self,strand):
        if strand.dna_strand is None: return strand
        strand.dna_strand = strand.dna_strand[::-1] 
        return strand


    
class AppendSequencePipeline(AppendSequence,DNAtoDNA):
    def __init__(self,seq,CodecObj=None,Policy=None,isPrimer=False,ignore=False,handler="ed",search_range=50):
        AppendSequence.__init__(self,seq,CodecObj=CodecObj,Policy=Policy,isPrimer=isPrimer,handler=handler,search_range=search_range)
        DNAtoDNA.__init__(self)
        self._ignore=ignore
    def _encode(self,strand):
        #this is a wrapper around the basic AppendSequence _encode so that the strand interface
        strand.dna_strand=AppendSequence._encode(self,strand.dna_strand)
        return strand
    def _decode(self,strand):
        if strand.dna_strand is None or self._ignore or self._seq=="": return strand
        strand_before =strand.dna_strand
        strand.dna_strand=AppendSequence._decode(self,strand.dna_strand)
        if strand_before==strand.dna_strand:
            strand.dna_strand = None
        return strand
 

if __name__ == "__main__":
    import random
    
    #cut = InsertMidSequence('AGATATAGGG',Policy=NoTolerance())
    #pre = PrependSequence('AAAAAAAAAAAG',CodecObj=cut,Policy=NoTolerance())
    #app = AppendSequence('CAAAAAAAAAAA',CodecObj=pre,Policy=NoTolerance())

    cut = InsertMidSequence('AGATATAGGG')
    pre = PrependSequence('TAAAGGAAAAAG',CodecObj=cut)
    app = AppendSequence( 'CAAAATATAAAA',CodecObj=pre)
    app = app
    
    match = 0
    for _ in range(10000):
        strand = [ random.choice('AGCT') for _ in range(100) ]
        strand = "".join(strand)

        original = strand

        strand = app.encode(strand)

        copy = [ random.choice('AGCT') for _ in range(10) ] + [ _ for _ in strand ] + [ random.choice('AGCT') for _ in range(10) ]
        copy[ random.randint(0,9) ] = random.choice('AGCT')
        copy[ random.randint(140,150) ] = random.choice('AGCT')
        copy[ random.randint(75,78) ] = random.choice('AGCT')
        copy = "".join(copy)

        
        new =  app.decode(copy)
        if new == original:
            match += 1
        else:
            print ("--",len(original)-len(new))
            print (original)
            print (new)
            
    print (match / 10000.0 * 100)
