from dnastorage.codec.base import *
from dnastorage.codec import base_conversion
from dnastorage.codec_types import *
from dnastorage.strand_representation import *
from random import randint
import editdistance as ed

from dnastorage.util.stats import stats


#In this order, no codewords have a hamming distance
#smaller than 2, and only two consecutive codewords have a hamming
#distance of 2. This helps make sure that consecutive indices vary by
#at least two positions

cfc_all = ['CGTGAGCA', 'ACACGATG', 'ATAGCTGC', 'ATAGAGAT', 'CAGTCATC',
           'TCATGCAT', 'TCGACGCA', 'CTATGTCG', 'ACGAGTCA', 'CACATGCG',
           'ATGTGTGC', 'TGTGACAG', 'TCATAGCA', 'TCAGATAG', 'ACGATCTC',
           'CGTATATA', 'CAGACTGA', 'CTCAGCTA', 'CTACTACA', 'ATGAGAGT',
           'ACTGTGCT', 'CATATGAC', 'GTCAGAGT', 'GTACGCAC', 'CTAGACTC',
           'TACGTGTC', 'CATAGAGT', 'TCGAGTCG', 'ATCATGCG', 'CGAGTGTG',
           'GCTCTCTA', 'TCGTATAC', 'TGTCTGCA', 'GTGTATCG', 'TCACATAC',
           'ATCAGACA', 'CACAGCGA', 'CTACTGAC', 'ATCATGTC', 'GCGTGATA',
           'ACGAGTGC', 'TCGAGATG', 'CGTAGTCA', 'ATACATCT', 'GTGACGCG',
           'CAGTCACG', 'CTCGTGCG', 'GCAGAGTG', 'CTCGCGCA', 'TCACTGTC',
           'GCTCTGAG', 'ATCAGTCG', 'ATCGAGCA', 'TATCATCT', 'CTCACTGA',
           'CTGTATAT', 'ATCGACAG', 'ATCAGCGA', 'ACTGTATG', 'GTAGTGTG',
           'CGTGTGCT', 'CGCGACAG', 'GCAGATAT', 'CTCGTACA', 'GCTCGCAT',
           'CTGTGCAT', 'AGTAGCGC', 'CTACTGCG', 'CTCGCACT', 'TCTCGCTA',
           'CTATGCAC', 'GTCAGCTG', 'CTATCTGA', 'GCAGATCG', 'GAGACTCT',
           'ACTGCACT', 'GCGATGCT', 'TGAGTAGA', 'CTAGATAT', 'GTACGCTA',
           'TCATACAC', 'GTAGTATA', 'TACGTACT', 'CGTAGACG', 'GCAGACAG',
           'TGCGAGAT', 'ACATAGAT', 'GTCTACAG', 'GTCTCTGA', 'GTGAGTCA',
           'CATGACAT', 'CGACTCTG', 'CGATCACG', 'CAGCGTAC', 'GTCTCTCT',
           'ATGTATGA', 'ACGTCGTC', 'TCGAGTAC', 'TCGATCGC', 'GTCTAGCT',
           'GCATCTCT', 'ATAGCGCA', 'CATATGCT', 'CACAGTAT', 'AGTAGCTG',
           'TACGACAT', 'ACTAGAGC', 'GACACGCT', 'CTATGTGC', 'AGACTGCG',
           'TCGTGCAC', 'GTGAGTGC', 'CACACGTG', 'GAGAGATG', 'ATCGTGCT',
           'CACAGTCG', 'ACACACTC', 'GTGATCTC', 'GTCAGTAT', 'TATCACGT',
           'ACTGTAGC', 'CTCGACTA', 'CAGACGCT', 'CAGAGATC', 'TACTGATA',
           'CATGCTGC', 'ATGTGATA', 'CGACGATG', 'GCTCTGCT', 'TGCGAGCA',
           'ATCACAGA', 'GCTCGCTG', 'GCTCTGTC', 'GCATCTGA', 'AGACTACA',
           'CGATCATC', 'TAGTAGAG', 'GTCATGAT', 'CACATGTC', 'ATGTGACT',
           'CTGTAGTC', 'GTGTATGC', 'CATGCACG', 'GCAGCTCA', 'GCAGCACG',
           'GTCGATCT', 'CTGTGTAC', 'GAGACATC', 'CGTGATGC', 'CAGCTCTA',
           'ATGAGCTG', 'TAGTGACA', 'GCATGTGC', 'GCGTAGCT', 'TACGAGCT',
           'ACGACTGA', 'GACGACGC', 'CTACTCTC', 'CACACAGA', 'CACATGAT',
           'ACGTAGCA', 'ATAGCGTC', 'ACTGCTCA', 'TATCTGCT', 'CTGTACAC',
           'CGTAGAGA', 'CAGTGATA', 'CTGTACTA', 'ACATGACT', 'TACGTGCA',
           'TGTCTGAT', 'ACACGAGC', 'GCAGATGC', 'TAGTAGCA', 'TCACATCT',
           'ACGTACAG', 'TGTATAGA', 'TCAGTACA', 'GTCGCGTC', 'CGACGTCT',
           'CACTACGC', 'TAGTGTAC', 'GACGACAG', 'TATAGCAC', 'CTATGCTA',
           'CTGTGCGC', 'ACGTATGC', 'GACACTGA', 'CGATAGAT', 'CGACACTC',
           'CTAGATGC', 'CAGCGCAT', 'TCGTGTCA', 'CACACATC', 'TCATATCG',
           'ATCAGATG', 'CACATATA', 'TCTCGATG', 'GTCTGCTA', 'AGTGAGCG',
           'CGTGAGAT', 'CGCGATCT', 'CGAGTATA', 'GTAGTGCA', 'CTACGAGA',
           'CAGAGCTA', 'TCAGCGTC', 'GTGTAGCA', 'TCGACGTG', 'CATGTGTC',
           'TCTAGAGT', 'CAGATCTG', 'CGTATCTC', 'TGTATCTA', 'TCAGTGCG',
           'ACGATATA', 'TGAGAGCG', 'GCATGCTA', 'ATCGATAT', 'GAGAGCAC',
           'CACTACAG', 'AGACTAGC', 'CAGTGTCG', 'AGTAGTCT', 'GCAGACTC',
           'TACTCTGA', 'TATGTAGA', 'CTGTAGAG', 'CGATACGC', 'GCGAGAGT',
           'TCGCGCTG', 'TCACACGC', 'TAGTGATG', 'CACTAGCG', 'CATATCTA',
           'TCATAGTC', 'CTATCTCT', 'GTGAGACG', 'TGACGCTG', 'TCTCGAGC',
           'TCAGCACT', 'ATAGTACA', 'ACACAGCA', 'CTACTAGC', 'CTAGACAG',
           'CAGTGCAC', 'CTAGTGTC', 'GCGTGACT', 'TGTGAGCT', 'GCAGCGCT',
           'CAGTGACT', 'GCAGAGCA', 'GCGCGCTA', 'TCGACATC', 'TATAGCTA',
           'GTAGTACT']

cfc_extra = ['TACGATGA', 'CACTGCTG', 'CGTAGATC', 'TCACTCTA']
cfc_special_header = 'ATCGATGC'

#cfc = ['CTACACGA', 'CATGCAGT', 'CTATGCGA', 'CACGAGCT', 'CTGCGTCT',
#'CGACGTGA', 'CATGTAGA', 'CGTGCACT', 'CAGCATCA', 'CGCACTGA',
#'CGCGCTGA', 'CAGATCGT', 'CTCTATAT', 'CATCTGCA', 'CACTACGT',
#'CTAGATGA', 'CTGAGTAT', 'CGCACGTA', 'CATGATCA', 'CGCATAGA',
#'CGACGACA', 'CACATAGA', 'CATAGCAT', 'CTGATACT', 'CTAGTACA',
#'CTAGTATA', 'CTCGTATA', 'CTAGACGA', 'CTATCTCT', 'CGATCTCA',
#'CGTCAGTA', 'CGCGAGAT', 'CGTAGTAT', 'CGTGCTGT', 'CACGTCGA',
#'CGACGTCT', 'CAGCGCAT', 'CACACGTA', 'CGTGCTCT', 'CTGCTCAT',
#'CTGCTAGT', 'CTGATATA', 'CGCACTCT', 'CATAGCGT', 'CTATACGT',
#'CACAGAGA', 'CTACAGTA', 'CAGCATAT', 'CAGCAGTA', 'CATCGAGT',
#'CTCGATAT', 'CAGCATCT', 'CTGAGACT', 'CACAGTCT', 'CATGTCGA',
#'CACGAGAT', 'CACATACA', 'CGCATGTA', 'CGCTAGAT', 'CTGTGCGT',
#'CTAGTGCA', 'CTAGAGAT', 'CTACATGA', 'CACTATAT', 'CAGCATGA',
#'CTGTACAT', 'CACTGTCA', 'CACTGCGA', 'CTCTATGA', 'CTGCTCGA',
#'CGATGTGT', 'CATCGAGA', 'CTGTCTCA', 'CTGACTCA', 'CTATCACT',
#'CGAGCATA', 'CAGATGCA', 'CAGACTCT', 'CATGAGCT', 'CTGATGTA',
#'CGTGTGCA', 'CAGACAGT', 'CGCGTCGA', 'CTACACAT', 'CTCTGCGA',
#'CTAGCTCA', 'CTCGCTGT', 'CTGTATAT', 'CGTAGACT', 'CAGAGCGT',
#'CATGCATA', 'CGCGTAGT', 'CAGTGTCT', 'CAGATCAT', 'CAGCGAGT',
#'CTGCGTAT', 'CGATCTGT', 'CTATCTCA', 'CACATGCT', 'CTACATCA',
#'CACAGTAT', 'CGAGCACT', 'CATATACA', 'CGACGAGA', 'CAGAGTCA',
#'CTCGAGCT', 'CAGCAGCT', 'CGCGACAT', 'CGCATGCA', 'CGATGCAT',
#'CGCACTCA', 'CTGATCAT', 'CTATCGCT', 'CACGCTGA', 'CACTGTGT',
#'CTAGAGTA', 'CGAGCGTA', 'CGCAGCGA', 'CTATAGTA', 'CGACTACT',
#'CGCGCTCT', 'CTGTAGCT', 'CGCTAGCA', 'CTAGATGT', 'CGACGTAT',
#'CTATCAGT', 'CAGCGTCT', 'CACTACGA', 'CACGCTGT', 'CACACTGA',
#'CAGAGTGT', 'CGCGATGA', 'CAGTGCGA', 'CTATAGCA', 'CGATAGTA',
#'CAGCTACT', 'CTCGTACA', 'CGCTGCGA', 'CAGTGACT', 'CTCTACGA',
#'CTAGTAGT', 'CGCGTACA', 'CTATCTGT', 'CGCGCTCA', 'CATATAGT',
#'CGTGTACT', 'CACGATCA', 'CACGCAGT', 'CATGTCGT', 'CATCGTCT',
#'CGATGAGT', 'CTCTGTCA', 'CGCGTCGT', 'CGAGCGCT', 'CTAGTGCT',
#'CTACAGCT', 'CGCGCTGT', 'CTATGCAT', 'CGTCATCA', 'CGATGCGT',
#'CGACGAGT', 'CTGTATCT', 'CGCTACGT', 'CTGAGATA', 'CTCGCAGT',
#'CTGTATCA', 'CGATGATA', 'CTGCGTGT', 'CGCATACA', 'CGCATGCT',
#'CAGCAGCA', 'CTAGCTGA', 'CTGTACGT', 'CTAGCTGT', 'CTGATCGT',
#'CATCTACT', 'CTCGAGAT', 'CTAGCGTA', 'CGCTACGA', 'CGTGCAGT',
#'CACATACT', 'CATGTGCA', 'CGTGCGCT', 'CAGAGCAT', 'CGATAGAT',
#'CTCGATCA', 'CATAGTCA', 'CATATGAT', 'CGTAGCAT', 'CATCGTGA',
#'CGCACGCT', 'CTGCTACT', 'CGTGCTCA', 'CACAGTCA', 'CGCTATAT',
#'CACACTCT', 'CTGTATGA', 'CGATACGA', 'CAGATACA', 'CACGATAT',
#'CTATGTCT', 'CGCGTGAT', 'CATATGCT', 'CTGATAGT', 'CGTGCTGA',
#'CGATCTCT', 'CTCTGCAT', 'CTATAGAT', 'CTAGCATA', 'CTGAGTGT',
#'CATGAGAT', 'CATGATGA', 'CTATACGA', 'CTAGTAGA', 'CTACGCAT',
#'CATGCGCT', 'CTATCTGA', 'CTATCATA', 'CTGTAGAT', 'CGATAGCT',
#'CTGTATGT', 'CACGTGAT', 'CTAGCGCT', 'CAGACTGT', 'CAGCGTCA',
#'CTAGTGAT', 'CTGATGCA', 'CGCGTATA', 'CTACATAT', 'CTGCGCAT',
#'CATATGTA', 'CATAGTCT', 'CGATGTAT', 'CATATGCA', 'CGCTGCAT',
#'CTCGAGTA', 'CGTAGAGA', 'CATGAGTA', 'CGCGCAGT', 'CGCGAGCT',
#'CTCGCTCA', 'CGTGCATA', 'CTGAGCGT', 'CTGCGTCA', 'CTGCTGCA',
#'CGCGTACT', 'CAGCGTAT', 'CGATGCGA', 'CTAGACGT', 'CGATCGTA',
#'CTAGCTCT', 'CTGAGTGA', 'CAGCTGCA', 'CTGTGAGT', 'CTGCTCGT',
#'CTCGTACT']

cfc = cfc_all[0:256]

cfc_inv = {}
def create_cfc_inv():
    global cfc_inv
    for i,c in enumerate(cfc):
        cfc_inv[c] = i
    cfc_inv[None] = -1000



    
class CommaFreeCodec(BaseCodec):
    def __init__(self,numberBytes,CodecObj=None,keyBytes=3):
        BaseCodec.__init__(self,CodecObj)
        self._keyWidthInBytes=keyBytes
        self._numberBytes = numberBytes
        global create_cfc_inv
        create_cfc_inv()

    def _encode(self,packet):
        key = base_conversion.convertIntToBytes(packet[0],self._keyWidthInBytes)
        value = packet[1]
        strand = key + value
        enc_strand = [ cfc[s] for s in strand ] #lookup byte s
        return "".join(enc_strand)

    def _decode_cfc(self, s):
        global cfc_inv
        if s in cfc_inv:
            return cfc_inv[s]
        return None

    def _decode_helper(self,s):
        l = self._numberBytes
        split = [ s[i:i+8] for i in range(0,l*8) ]
        for i in range(len(split)):
            split[i] = self._decode_cfc(split[i])
        dec = [ None for _ in range(l) ]
        i = 0
        prior = True
        n = 0
        for n in range(l):
            if i >= l*8:
                break
            if split[i] != None:
                dec[n] = split[i]
                i += 8
                prior = True
            elif split[i] == None:
                if i+8>=l*8:
                    break
                # maybe have an insertion, deletion, or substitution
                if split[i+8] != None:
                    # guess substitution, leave entry as None
                    i += 8
                    prior = False
                elif split[i+8] == None: # somehow we are off track
                    if prior: 
                        # this is the first one, so look forward from i for next non-None
                        k = i
                    else:
                        k = i-4

                    while k < i+4:
                        if split[k] != None:
                            i = k;
                            dec[n] = split[i]
                            prior = True
                            break
                        k += 1
                    # whether found or not, move to next possible symbol
                    if dec[n] == None:
                        prior = False
                    i += 8
        return dec


    def _decode(self,s):
        dec = self._decode_helper(s)
        #have a None object, generate random byte values
        key_array=[]
        #need to handle case of not being able to come up with a complete key
        if None in dec[0:self._keyWidthInBytes]:
            for i in dec[0:self._keyWidthInBytes]:
                if i is None:
                    key_array.append(randint(0,255))
                else:
                    key_array.append(i)
        else:
            key_array=dec[0:self._keyWidthInBytes]    
        key = base_conversion.convertBytesToInt(key_array)
        #return key,"".join([ chr(x) for x in dec[self._keyWidthInBytes:]])
        
        #Return an array of values rather than joining them as characters, possibly FIX_ME?
        return key,dec[self._keyWidthInBytes:]

#wrapper so that the old comma free codec plays nice with pipeline based processing
class CommaFreeCodecPipeline(CommaFreeCodec,CWtoDNA): 
    def __init__(self,numberBytes=15,CodecObj=None,Policy=None):
        CommaFreeCodec.__init__(self,numberBytes,CodecObj=CodecObj,keyBytes=1)
        CWtoDNA.__init__(self)
    def _encode(self,strand):
        self._numberBytes=len(strand.codewords)        
        strand.dna_strand=CommaFreeCodec._encode(self,(strand.codewords[0],strand.codewords[1:]))        
        return strand
    def _decode(self,strand):
        a, b=CommaFreeCodec._decode(self,strand.dna_strand)
        strand.codewords = [a]+b
        return strand

    def _encode_header(self):
        return convertIntToBytes(self._numberBytes,1)

    def _decode_header(self,buf):
        self._numberBytes=convertBytesToInt(buf[0:1])
        return buf[1::]
    

class CommaFreeCodewords(BaseCodec):
    def __init__(self,numberSymbols,CodecObj=None,Policy=None):
        BaseCodec.__init__(self,CodecObj,Policy=Policy)
        self._numberBytes = numberSymbols
        global create_cfc_inv
        create_cfc_inv()

    def _encode(self,strand):
        enc_strand = [ cfc[s] for s in strand ]
        return enc_strand

    def _decode_cfc(self, s):
        global cfc_inv
        if s in cfc_inv:
            return cfc_inv[s]
        return None

    def _decode(self,s):
        dec = self._decode_helper(s)
        return dec

    def exact_vote(self, s):
        exact = []        
        for i in range(len(s)):
            exact.append(self._decode_cfc(s[i:i+8]))
        return exact

    def inexact_vote(self, s):
        #print s
        global cfc
        D = {}
        for x in cfc:
            d = ed.eval(s,x)
            D[d] = D.get(d,[]) + [x]
        if 0 in D:
            assert len(D[0])==1
            return [100,D[0][0],0]        
        elif 1 in D:
            if len(D[1])==1:
                return [100,D[1][0],8]
            elif len(D[1])>1:
                p = D[1][ randint(0,len(D[1])-1) ]
                return [1.0/len(D[1])*100,p,8]
        elif 2 in D:
            if len(D[2])==1:
                return [100,D[2][0],8]
            elif len(D[2])>1:
                p = D[2][ randint(0,len(D[2])-1) ]
                return [1.0/len(D[2])*100,p,8]
        else:
                return [0.0,None,8]


    def get_next_exact(self, i, exact):
        while i < len(exact):
            if exact[i] != None:
                return i
            i += 1
        return i

    def _decode_helper(self, s):
        numSyms = self._numberBytes
        exact = self.exact_vote(s)

        #print len(s)
        
        new_strand = []
        i = 0
        while i < len(s) and len(new_strand) < numSyms:
            #print i
            if not(exact[i] is None):
                if i%8 != 0:
                    stats.inc("CFC8::misAlignedCodeword")
                new_strand.append(exact[i])
                #print i,"->",exact[i]
                i += 8
            else:
                #print "missing codeword:",exact[i:i+8]
                err = DNABadCodeword("Missing expected CFC8 symbol")
                if self._Policy.allow(err):
                    j = self.get_next_exact(i,exact)
                    skipped = int(round((j-i)/8.0))
                    i_tmp = i
                    for k in range(skipped):
                        #new_strand.append(-1)
                        v = self.inexact_vote(s[i_tmp:i_tmp+8])
                        if v[0]==100:
                            new_strand.append(cfc_inv[v[1]])
                            i_tmp+=v[2]
                        else:
                            i_tmp+=8
                            new_strand.append(-1)
                    i = j
                else:
                    raise err
        
        if len(new_strand) != numSyms:
            err = DNAStrandPayloadWrongSize("Payload wrong size in CFC8")
            stats.inc("CFC8:payloadWrongSize")
            if self._Policy.allow(err):            
                if len(new_strand) < numSyms:
                    while len(new_strand) < numSyms:
                        new_strand.append(-1)
                else:
                    new_strand = new_strand[0:numSyms]
            else:
                raise err


        cnt=0
        found=False
        for _ in new_strand:
            stats.inc("CFC8::TotalCodewords")
            if _ == -1:
                found=True
                stats.inc("CFC8::UnidentifiedCodeword")
            else:
                stats.inc("CFC8::IdentifiedCodeword")
        if found:
            #print "found error"
            stats.inc("CFC8::strandsWithErrors")
            stats.inc("CFC8::strandsTotal")
        else:
            stats.inc("CFC8::strandsTotal")
            
        return new_strand
                
            # if j<3:
            #     i += (j-i)
            #     continue
            # v = inexact_vote(s[i:i+10])
            # if v[0] == 100:
            #     new_strand.append(v[1])
            #     if j < len(s):
            #         i += (j-i)
            #     else:
            #         i += v[0]
            # else:
            #     i += 1

#wrapper so that the old comma free codec plays nice with pipeline based processing
class CommaFreeCodewordsPipeline(CommaFreeCodewords,CWtoDNA):
    def __init__(self,numberBytes=15,CodecObj=None,Policy=None):
        CommaFreeCodewords.__init__(self,numberBytes,CodecObj=CodecObj)
        CWtoDNA.__init__(self)
    
    def _encode(self,strand):
        self._numberBytes=len(strand.codewords)
        strand.dna_strand="".join(CommaFreeCodewords._encode(self,strand.codewords))
        print (strand.dna_strand)        
        return strand
    def _decode(self,strand):
        strand.codewords=CommaFreeCodewords._decode(self,strand.dna_strand)
        print (strand.codewords)
        
        return strand

    def _encode_header(self):
        return convertIntToBytes(self._numberBytes,1)

    def _decode_header(self,buf):
        self._numberBytes=convertBytesToInt(buf[0:1])
        return buf[1::]
    


    
if __name__ == "__main__":
    from dnastorage.primer.primer_util import *
    from dnastorage.primer.design import *
    from random import randint,shuffle
    import sys
    
    def countGC(s):
        l = [ _ for _ in s if _ == 'G' or _ == 'C' ]
        return len(l)


    """
    If s can be rotated to create a value equal to s, then it has a self cycle
    """
    def self_cycle_check(s,l=8):
        t = s+s
        # create all rotations of s
        all = [ t[i:i+l] for i in range(1,l) ]
        # if s is present in l, there's a cycle
        if s in all:
            return True
        return False

    """
    Make sure that b does not appear in a cycle of a
    """
    def same_cycle(a,b,l=8):
        t = a+a
        ll = [ t[i:i+l] for i in range(1,l) ]
        return b in ll

    def comma_free_check(all,l):
        # pick random starting point
        while True:
            i = randint(0,len(all))
            if self_cycle_check(all[i],l):
                continue        
            start = [ all[i] ]
            break

        #randomize list order
        shuffle(all)

        for i,s in enumerate(all):
            if (s in start):
                continue

            if self_cycle_check(s,l):
                continue

            found = True
            for a in start:
                if same_cycle(s,a,l):
                    found = False
                    break

            if not found:
                continue

            for k,a in enumerate(start):
                for j in range(k+1,len(start)):
                    b = start[j]
                    if a == b:
                        continue
                    t = a + b + a
                    if t.find(s) >= 0: #  or t.find(rs) >= 0:
                        found = False
                        break

                    t = a + s + a
                    if t.find(b) >= 0:
                        found = False
                        break

                    t = b + s + b
                    if t.find(a) >= 0:
                        found = False
                        break

                if not found:
                    break

            if not found:
                if i%100 == 0:
                    print ("{}% tested".format(float(i)/len(all)*100.0))
                continue
            else:
                print ("{} {} - {}".format(s,len(start), start[-5:]))
                start.append(s)

        return start

    
    def print_characteristics(all):
        base = { 'A' : [],
                 'C' : [],
                 'G' : [],
                 'T' : []  } 
        gc = {}
        for s in all:
            base[s[0]].append(s)
            c = countGC(s)
            if not (c in gc):
                gc[c] = [ ]
            #print "{} {}".format(d,c)
            gc[c].append(s)

        for key,item in gc.items():
            print ("{} - {}".format(key,len(item)))

        for key,item in base.items():
            print ("{}:{} {}".format(key,len(item),item[0:4]))

        D = {}
        similar=[]
        for j,a in enumerate(all):
            for k,b in enumerate(all):
                if k<=j:
                    continue
                d = hamming_distance(a,b)
                D[d] = D.get(d,0)+1
                if d == 1:
                    similar.append( [j,k] )
        print (D)
        print (similar)


        
    create_cfc_inv()

    print_characteristics(cfc)
    #comma_free_check(cfc,8)
    print (len(cfc))
    
    arr = [ randint(0,255) for _ in range(20) ]


    r = 0
    for i in range(10000):
        a = cfc[randint(0,255)]
        b = cfc[randint(0,255)]
        if a[-1:] == b[0:1]:
            r += 1
    print ("repeat rate = {}%".format(r/10000.00 * 100))
        
    sys.exit(0)
    
    import editdistance as ed
        
    r = 0
    for i in range(10000):
        a = cfc[randint(0,255)]
        a_copy = a
        a = [ _ for _ in a ]
        rr = randint(0,7)
        sym = a[rr]
        while a[rr] == sym:
            a[ rr ] = random.choice('ACGT')
        D = {}
        for x in cfc:
            d = ed.eval("".join(a),x)
            D[d] = D.get(d,[]) + [x]
        if 0 in D:
            pass # didn't find it
        elif 1 in D:
            if len(D[1])==1 and a_copy in D[1]:
                r += 1
            elif len(D[1])>1 and a_copy in D[1]:
                p = D[1][ randint(0,len(D[1])-1) ]
                if p == a_copy:
                    r += 1
    print ("substitution rate = {}%".format(r/10000.00 * 100))

    r = 0
    for i in range(10000):
        a = cfc[randint(0,255)]
        a_copy = a
        a = [ _ for _ in a ]
        rr = randint(0,7)
        a = a[0:rr] + [random.choice('ACGT')] + a[rr:]
        #print a,a_copy,ed.eval("".join(a),a_copy)
        D = {}
        for x in cfc:
            d = ed.eval("".join(a),x)
            D[d] = D.get(d,[]) + [x]
        if 0 in D:
            pass # didn't find it
        elif 1 in D:
            if len(D[1])==1 and a_copy in D[1]:
                r += 1
            elif len(D[1])>1 and a_copy in D[1]:
                p = D[1][ randint(0,len(D[1])-1) ]
                if p == a_copy:
                    r += 1
            
    print ("insertion rate = {}%".format(r/10000.00 * 100))

    r = 0
    for i in range(10000):
        a = cfc[randint(0,255)]
        a_copy = a
        a = [ _ for _ in a ]
        rr = randint(0,7)
        a = a + [random.choice('ACGT')]
        #print a,a_copy,ed.eval("".join(a),a_copy)
        D = {}
        for x in cfc:
            d = ed.eval("".join(a),x)
            D[d] = D.get(d,[]) + [x]
        if 0 in D:
            pass # didn't find it
        elif 1 in D:
            if len(D[1])==1 and a_copy in D[1]:
                r += 1
            elif len(D[1])>1 and a_copy in D[1]:
                p = D[1][ randint(0,len(D[1])-1) ]
                if p == a_copy:
                    r += 1
            
    print ("insertion at end rate = {}%".format(r/10000.00 * 100))


    r = 0
    for i in range(10000):
        a = cfc[randint(0,255)]
        a_copy = a
        a = [ _ for _ in a ]
        rr = randint(0,7)
        a = a[0:rr]  + a[rr+1:] 
        #print a,a_copy,ed.eval("".join(a),a_copy)
        D = {}
        for x in cfc:
            d = ed.eval("".join(a),x)
            D[d] = D.get(d,[]) + [x]
        if 0 in D:
            pass # didn't find it
        elif 1 in D:
            if len(D[1])==1 and a_copy in D[1]:
                r += 1
            elif len(D[1])>1 and a_copy in D[1]:
                p = D[1][ randint(0,len(D[1])-1) ]
                if p == a_copy:
                    r += 1

        elif 2 in D:
            if len(D[2])==1 and a_copy in D[2]:
                r += 1
            elif len(D[2])>1 and a_copy in D[2]:
                p = D[2][ randint(0,len(D[2])-1) ]
                if p == a_copy:
                    r += 1
        elif 3 in D:
            if len(D[3])==1 and a_copy in D[3]:
                r += 1
            elif len(D[3])>1 and a_copy in D[3]:
                p = D[3][ randint(0,len(D[3])-1) ]
                if p == a_copy:
                    r += 1


                    
    print ("deletion rate = {}%".format(r/10000.00 * 100))


    r = 0
    for i in range(10000):
        a = cfc[randint(0,255)]
        a_copy = a
        a = [ _ for _ in a ]
        rr = randint(0,7)
        a = a[0:-1]
        #print a,a_copy,ed.eval("".join(a),a_copy)
        D = {}
        for x in cfc:
            d = ed.eval("".join(a),x)
            D[d] = D.get(d,[]) + [x]
        if 0 in D:
            pass # didn't find it
        elif 1 in D:
            if len(D[1])==1 and a_copy in D[1]:
                r += 1
            elif len(D[1])>1 and a_copy in D[1]:
                p = D[1][ randint(0,len(D[1])-1) ]
                if p == a_copy:
                    r += 1

        elif 2 in D:
            if len(D[2])==1 and a_copy in D[2]:
                r += 1
            elif len(D[2])>1 and a_copy in D[2]:
                p = D[2][ randint(0,len(D[2])-1) ]
                if p == a_copy:
                    r += 1
        elif 3 in D:
            if len(D[3])==1 and a_copy in D[3]:
                r += 1
            elif len(D[3])>1 and a_copy in D[3]:
                p = D[3][ randint(0,len(D[3])-1) ]
                if p == a_copy:
                    r += 1
                    
    print ("deletion at end rate = {}%".format(r/10000.00 * 100))

    
    cf = CommaFreeCodewords(20)

    enc_arr = cf.encode(arr)

    print (enc_arr)
