#!/usr/bin/python
from dnastorage.codec.base import *
from dnastorage.codec import base_conversion
import random

dense_enc_table = [
'AAAA','AAAC','AAAT','AAAG','AACA','AACC','AACT','AACG',
'AATA','AATC','AATT','AATG','AAGA','AAGC','AAGT','AAGG',
'ACAA','ACAC','ACAT','ACAG','ACCA','ACCC','ACCT','ACCG',
'ACTA','ACTC','ACTT','ACTG','ACGA','ACGC','ACGT','ACGG',
'ATAA','ATAC','ATAT','ATAG','ATCA','ATCC','ATCT','ATCG',
'ATTA','ATTC','ATTT','ATTG','ATGA','ATGC','ATGT','ATGG',
'AGAA','AGAC','AGAT','AGAG','AGCA','AGCC','AGCT','AGCG',
'AGTA','AGTC','AGTT','AGTG','AGGA','AGGC','AGGT','AGGG',
'CAAA','CAAC','CAAT','CAAG','CACA','CACC','CACT','CACG',
'CATA','CATC','CATT','CATG','CAGA','CAGC','CAGT','CAGG',
'CCAA','CCAC','CCAT','CCAG','CCCA','CCCC','CCCT','CCCG',
'CCTA','CCTC','CCTT','CCTG','CCGA','CCGC','CCGT','CCGG',
'CTAA','CTAC','CTAT','CTAG','CTCA','CTCC','CTCT','CTCG',
'CTTA','CTTC','CTTT','CTTG','CTGA','CTGC','CTGT','CTGG',
'CGAA','CGAC','CGAT','CGAG','CGCA','CGCC','CGCT','CGCG',
'CGTA','CGTC','CGTT','CGTG','CGGA','CGGC','CGGT','CGGG',
'TAAA','TAAC','TAAT','TAAG','TACA','TACC','TACT','TACG',
'TATA','TATC','TATT','TATG','TAGA','TAGC','TAGT','TAGG',
'TCAA','TCAC','TCAT','TCAG','TCCA','TCCC','TCCT','TCCG',
'TCTA','TCTC','TCTT','TCTG','TCGA','TCGC','TCGT','TCGG',
'TTAA','TTAC','TTAT','TTAG','TTCA','TTCC','TTCT','TTCG',
'TTTA','TTTC','TTTT','TTTG','TTGA','TTGC','TTGT','TTGG',
'TGAA','TGAC','TGAT','TGAG','TGCA','TGCC','TGCT','TGCG',
'TGTA','TGTC','TGTT','TGTG','TGGA','TGGC','TGGT','TGGG',
'GAAA','GAAC','GAAT','GAAG','GACA','GACC','GACT','GACG',
'GATA','GATC','GATT','GATG','GAGA','GAGC','GAGT','GAGG',
'GCAA','GCAC','GCAT','GCAG','GCCA','GCCC','GCCT','GCCG',
'GCTA','GCTC','GCTT','GCTG','GCGA','GCGC','GCGT','GCGG',
'GTAA','GTAC','GTAT','GTAG','GTCA','GTCC','GTCT','GTCG',
'GTTA','GTTC','GTTT','GTTG','GTGA','GTGC','GTGT','GTGG',
'GGAA','GGAC','GGAT','GGAG','GGCA','GGCC','GGCT','GGCG',
'GGTA','GGTC','GGTT','GGTG','GGGA','GGGC','GGGT','GGGG'
]

dense_dec_table = {
'AAAA' : 0,
'AAAC' : 1,
'AAAT' : 2,
'AAAG' : 3,
'AACA' : 4,
'AACC' : 5,
'AACT' : 6,
'AACG' : 7,
'AATA' : 8,
'AATC' : 9,
'AATT' : 10,
'AATG' : 11,
'AAGA' : 12,
'AAGC' : 13,
'AAGT' : 14,
'AAGG' : 15,
'ACAA' : 16,
'ACAC' : 17,
'ACAT' : 18,
'ACAG' : 19,
'ACCA' : 20,
'ACCC' : 21,
'ACCT' : 22,
'ACCG' : 23,
'ACTA' : 24,
'ACTC' : 25,
'ACTT' : 26,
'ACTG' : 27,
'ACGA' : 28,
'ACGC' : 29,
'ACGT' : 30,
'ACGG' : 31,
'ATAA' : 32,
'ATAC' : 33,
'ATAT' : 34,
'ATAG' : 35,
'ATCA' : 36,
'ATCC' : 37,
'ATCT' : 38,
'ATCG' : 39,
'ATTA' : 40,
'ATTC' : 41,
'ATTT' : 42,
'ATTG' : 43,
'ATGA' : 44,
'ATGC' : 45,
'ATGT' : 46,
'ATGG' : 47,
'AGAA' : 48,
'AGAC' : 49,
'AGAT' : 50,
'AGAG' : 51,
'AGCA' : 52,
'AGCC' : 53,
'AGCT' : 54,
'AGCG' : 55,
'AGTA' : 56,
'AGTC' : 57,
'AGTT' : 58,
'AGTG' : 59,
'AGGA' : 60,
'AGGC' : 61,
'AGGT' : 62,
'AGGG' : 63,
'CAAA' : 64,
'CAAC' : 65,
'CAAT' : 66,
'CAAG' : 67,
'CACA' : 68,
'CACC' : 69,
'CACT' : 70,
'CACG' : 71,
'CATA' : 72,
'CATC' : 73,
'CATT' : 74,
'CATG' : 75,
'CAGA' : 76,
'CAGC' : 77,
'CAGT' : 78,
'CAGG' : 79,
'CCAA' : 80,
'CCAC' : 81,
'CCAT' : 82,
'CCAG' : 83,
'CCCA' : 84,
'CCCC' : 85,
'CCCT' : 86,
'CCCG' : 87,
'CCTA' : 88,
'CCTC' : 89,
'CCTT' : 90,
'CCTG' : 91,
'CCGA' : 92,
'CCGC' : 93,
'CCGT' : 94,
'CCGG' : 95,
'CTAA' : 96,
'CTAC' : 97,
'CTAT' : 98,
'CTAG' : 99,
'CTCA' : 100,
'CTCC' : 101,
'CTCT' : 102,
'CTCG' : 103,
'CTTA' : 104,
'CTTC' : 105,
'CTTT' : 106,
'CTTG' : 107,
'CTGA' : 108,
'CTGC' : 109,
'CTGT' : 110,
'CTGG' : 111,
'CGAA' : 112,
'CGAC' : 113,
'CGAT' : 114,
'CGAG' : 115,
'CGCA' : 116,
'CGCC' : 117,
'CGCT' : 118,
'CGCG' : 119,
'CGTA' : 120,
'CGTC' : 121,
'CGTT' : 122,
'CGTG' : 123,
'CGGA' : 124,
'CGGC' : 125,
'CGGT' : 126,
'CGGG' : 127,
'TAAA' : 128,
'TAAC' : 129,
'TAAT' : 130,
'TAAG' : 131,
'TACA' : 132,
'TACC' : 133,
'TACT' : 134,
'TACG' : 135,
'TATA' : 136,
'TATC' : 137,
'TATT' : 138,
'TATG' : 139,
'TAGA' : 140,
'TAGC' : 141,
'TAGT' : 142,
'TAGG' : 143,
'TCAA' : 144,
'TCAC' : 145,
'TCAT' : 146,
'TCAG' : 147,
'TCCA' : 148,
'TCCC' : 149,
'TCCT' : 150,
'TCCG' : 151,
'TCTA' : 152,
'TCTC' : 153,
'TCTT' : 154,
'TCTG' : 155,
'TCGA' : 156,
'TCGC' : 157,
'TCGT' : 158,
'TCGG' : 159,
'TTAA' : 160,
'TTAC' : 161,
'TTAT' : 162,
'TTAG' : 163,
'TTCA' : 164,
'TTCC' : 165,
'TTCT' : 166,
'TTCG' : 167,
'TTTA' : 168,
'TTTC' : 169,
'TTTT' : 170,
'TTTG' : 171,
'TTGA' : 172,
'TTGC' : 173,
'TTGT' : 174,
'TTGG' : 175,
'TGAA' : 176,
'TGAC' : 177,
'TGAT' : 178,
'TGAG' : 179,
'TGCA' : 180,
'TGCC' : 181,
'TGCT' : 182,
'TGCG' : 183,
'TGTA' : 184,
'TGTC' : 185,
'TGTT' : 186,
'TGTG' : 187,
'TGGA' : 188,
'TGGC' : 189,
'TGGT' : 190,
'TGGG' : 191,
'GAAA' : 192,
'GAAC' : 193,
'GAAT' : 194,
'GAAG' : 195,
'GACA' : 196,
'GACC' : 197,
'GACT' : 198,
'GACG' : 199,
'GATA' : 200,
'GATC' : 201,
'GATT' : 202,
'GATG' : 203,
'GAGA' : 204,
'GAGC' : 205,
'GAGT' : 206,
'GAGG' : 207,
'GCAA' : 208,
'GCAC' : 209,
'GCAT' : 210,
'GCAG' : 211,
'GCCA' : 212,
'GCCC' : 213,
'GCCT' : 214,
'GCCG' : 215,
'GCTA' : 216,
'GCTC' : 217,
'GCTT' : 218,
'GCTG' : 219,
'GCGA' : 220,
'GCGC' : 221,
'GCGT' : 222,
'GCGG' : 223,
'GTAA' : 224,
'GTAC' : 225,
'GTAT' : 226,
'GTAG' : 227,
'GTCA' : 228,
'GTCC' : 229,
'GTCT' : 230,
'GTCG' : 231,
'GTTA' : 232,
'GTTC' : 233,
'GTTT' : 234,
'GTTG' : 235,
'GTGA' : 236,
'GTGC' : 237,
'GTGT' : 238,
'GTGG' : 239,
'GGAA' : 240,
'GGAC' : 241,
'GGAT' : 242,
'GGAG' : 243,
'GGCA' : 244,
'GGCC' : 245,
'GGCT' : 246,
'GGCG' : 247,
'GGTA' : 248,
'GGTC' : 249,
'GGTT' : 250,
'GGTG' : 251,
'GGGA' : 252,
'GGGC' : 253,
'GGGT' : 254,
'GGGG' : 255
}

def dense_encode_byte(byte):
    return dense_enc_table[byte]

def dense_decode_byte(byte):
    assert dense_dec_table.has_key(byte)==True
    return dense_dec_table[byte]

def dense_encode(packet):
    array = bytearray(packet)
    strand = []
    for b in array:
        strand.append( dense_encode_byte(b) )
    return "".join(strand)

def dense_decode(strand):
    i = 0
    array = []
    while i < len(strand):
        s = strand[i:i+4]
        byte = dense_decode_byte(s)
        array.append(chr(byte))
        i+=4
    packet = bytearray(array)
    return packet


class RandomizeDenseCodec(BaseCodec):
    def __init__(self,CodecObj=None):
        BaseCodec.__init__(self,CodecObj)
        self._table = [ "AT" , "ACGT", "AGCATACG" ]
        self._normalized = False

    def _countRepetitions(self, ss):
        reps = 0
        last = 'A'
        for i,s in enumerate(ss):
            if i>0 and last == s:
                reps = reps + 1
            last = s
        return reps

    def _createMask(self, s, r):
        ss = [_ for _ in s]
        rr = [_ for _ in r]
        mask = [self._xor(a,b) for a,b in zip(ss,rr)]
        return mask

    def _check(self, s):
        matrix = []
        mins = []
        for i,r in enumerate(self._table):
            mi = self._createMask(s,r)
            reps = self._countRepetitions(mi)
            mrow = []
            for j,rr in enumerate(self._table):
                mj = self._createMask(mi,rr)
                nreps = reps - self._countRepetitions(mj)
                mrow.append(nreps)
            matrix.append(mrow)
            mins.append( [mrow.index(min(mrow)),min(mrow)])
        #print self._countRepetitions(s),matrix
        #print mins
        first = True
        smallest = -1
        for i,val in enumerate(mins):
            if i == val[0]:
                if first:
                    smallest = i
                    first = False
                else:
                    if val[1] < mins[smallest][1]:
                        smallest = i
        #print smallest
        #print ""


    def _randomize(self, s):
        lreps = []
        self._normalize(len(s))
        ss = [_ for _ in s]
        best = ss
        best_i = -1
        reps = self._countRepetitions(ss)
        min_reps = 0
        lreps.append( min_reps )
        for i,r in enumerate(self._table):
            rr = [_ for _ in r]
            n = [self._xor(a,b) for a,b in zip(ss,rr)]
            nreps = self._countRepetitions(n) - reps
            lreps.append(nreps)
            if nreps < min_reps:
                best = n
                best_i = i
                min_reps = nreps
        #print "Reps: ",lreps
        self._check(s)
        return "".join(best)

    def _normalize(self, le):
        if self._normalized == True:
            return
        for i,s in enumerate(self._table):
            l = len(s)
            if le % l == 0:
                s = s * (le / l)
                #print s
                self._table[i] = s
            else:
                r = le % l
                s = s * (le / l)
                s += s[0:r]
                #print s
                self._table[i] = s
            assert len(self._table[i]) == le

        self._table.append(''.join(random.choice("ACGT") for _ in range(le)))
        self._table.append(''.join(random.choice("ACGT") for _ in range(le)))

        self._normalized = True

    def _xor(self, a, b):
        if a == 'A':
            return b
        elif a == 'C':
            xor = { 'A' : 'C',
                    'C' : 'A',
                    'G' : 'T',
                    'T' : 'G'}
            return xor[b]
        elif a == 'G':
            xor = { 'A' : 'G',
                    'C' : 'T',
                    'G' : 'A',
                    'T' : 'C'}
            return xor[b]
        elif a == 'T':
            xor = { 'A' : 'T',
                    'C' : 'G',
                    'G' : 'C',
                    'T' : 'A'}
            return xor[b]
        else:
            return '?'

    def _encode(self,s):
        assert isinstance(s,str) or "Didn't get expected type."
        return self._randomize(s)

    def _decode(self,s):
        assert isinstance(s,str) or "Didn't get expected type."
        return self._unrandomize(s)

class DenseTableCodec(TableCodec):
    def __init__(self,CodecObj=None,keyWidth=16):
        TableCodec.__init__(self,CodecObj,keyWidth,keyWidth/4,4,1)

    def _enctab(self, val):
        #print "_enctab: ",val,dense_encode_byte(val)
        return dense_encode_byte(val)

    def _dectab(self, s):
        #print "_dectab: ",s,dense_decode_byte(s)
        return dense_decode_byte(s)

class DenseCodec(BaseCodec):
    def __init__(self,CodecObj=None,keyWidth=20):
        BaseCodec.__init__(self,CodecObj)
        self._keyWidth=keyWidth

    def _encode(self,packet):
        key = base_conversion.convertBase(4,packet[0],self._keyWidth)
        value = dense_encode(packet[1])
        return key+value

    def _decode(self,s):
        key = base_conversion.convertFromBase(4,s[0:self._keyWidth])
        value = dense_decode(s[self._keyWidth:])
        return key,value


if __name__ == "__main__":
    from base_conversion import *
    enc = []
    dec = {}
    for i in range(256):
        s = convertBase(4,i,4)
        s = [x for x in s]
        s.reverse()
        s = "".join(s)
        enc.append("'"+s+"'")
        dec[s] = i

    print "dense_enc_table = ["
    for i in range(len(enc)/8):
        if i==len(enc)/8-1:
            print ",".join(enc[i*8:i*8+8])
        else:
            print ",".join(enc[i*8:i*8+8])+","
    print "]"

    print "dense_dec_table = {"
    for i in range(len(enc)):
        s =  "{} : {}".format(enc[i],i)
        if i!=len(enc)-1:
            s += ","
        print s
    print "}"

    s = "12345678"
    print dense_encode(s)
    assert s == dense_decode(dense_encode(s))
