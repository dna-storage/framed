from dnastorage.codec.base import *
from dnastorage.codec import base_conversion

cfc = ['CTACACGA', 'CATGCAGT', 'CTATGCGA', 'CACGAGCT', 'CTGCGTCT', 'CGACGTGA', 'CATGTAGA', 'CGTGCACT', 'CAGCATCA', 'CGCACTGA', 'CGCGCTGA', 'CAGATCGT', 'CTCTATAT', 'CATCTGCA', 'CACTACGT', 'CTAGATGA', 'CTGAGTAT', 'CGCACGTA', 'CATGATCA', 'CGCATAGA', 'CGACGACA', 'CACATAGA', 'CATAGCAT', 'CTGATACT', 'CTAGTACA', 'CTAGTATA', 'CTCGTATA', 'CTAGACGA', 'CTATCTCT', 'CGATCTCA', 'CGTCAGTA', 'CGCGAGAT', 'CGTAGTAT', 'CGTGCTGT', 'CACGTCGA', 'CGACGTCT', 'CAGCGCAT', 'CACACGTA', 'CGTGCTCT', 'CTGCTCAT', 'CTGCTAGT', 'CTGATATA', 'CGCACTCT', 'CATAGCGT', 'CTATACGT', 'CACAGAGA', 'CTACAGTA', 'CAGCATAT', 'CAGCAGTA', 'CATCGAGT', 'CTCGATAT', 'CAGCATCT', 'CTGAGACT', 'CACAGTCT', 'CATGTCGA', 'CACGAGAT', 'CACATACA', 'CGCATGTA', 'CGCTAGAT', 'CTGTGCGT', 'CTAGTGCA', 'CTAGAGAT', 'CTACATGA', 'CACTATAT', 'CAGCATGA', 'CTGTACAT', 'CACTGTCA', 'CACTGCGA', 'CTCTATGA', 'CTGCTCGA', 'CGATGTGT', 'CATCGAGA', 'CTGTCTCA', 'CTGACTCA', 'CTATCACT', 'CGAGCATA', 'CAGATGCA', 'CAGACTCT', 'CATGAGCT', 'CTGATGTA', 'CGTGTGCA', 'CAGACAGT', 'CGCGTCGA', 'CTACACAT', 'CTCTGCGA', 'CTAGCTCA', 'CTCGCTGT', 'CTGTATAT', 'CGTAGACT', 'CAGAGCGT', 'CATGCATA', 'CGCGTAGT', 'CAGTGTCT', 'CAGATCAT', 'CAGCGAGT', 'CTGCGTAT', 'CGATCTGT', 'CTATCTCA', 'CACATGCT', 'CTACATCA', 'CACAGTAT', 'CGAGCACT', 'CATATACA', 'CGACGAGA', 'CAGAGTCA', 'CTCGAGCT', 'CAGCAGCT', 'CGCGACAT', 'CGCATGCA', 'CGATGCAT', 'CGCACTCA', 'CTGATCAT', 'CTATCGCT', 'CACGCTGA', 'CACTGTGT', 'CTAGAGTA', 'CGAGCGTA', 'CGCAGCGA', 'CTATAGTA', 'CGACTACT', 'CGCGCTCT', 'CTGTAGCT', 'CGCTAGCA', 'CTAGATGT', 'CGACGTAT', 'CTATCAGT', 'CAGCGTCT', 'CACTACGA', 'CACGCTGT', 'CACACTGA', 'CAGAGTGT', 'CGCGATGA', 'CAGTGCGA', 'CTATAGCA', 'CGATAGTA', 'CAGCTACT', 'CTCGTACA', 'CGCTGCGA', 'CAGTGACT', 'CTCTACGA', 'CTAGTAGT', 'CGCGTACA', 'CTATCTGT', 'CGCGCTCA', 'CATATAGT', 'CGTGTACT', 'CACGATCA', 'CACGCAGT', 'CATGTCGT', 'CATCGTCT', 'CGATGAGT', 'CTCTGTCA', 'CGCGTCGT', 'CGAGCGCT', 'CTAGTGCT', 'CTACAGCT', 'CGCGCTGT', 'CTATGCAT', 'CGTCATCA', 'CGATGCGT', 'CGACGAGT', 'CTGTATCT', 'CGCTACGT', 'CTGAGATA', 'CTCGCAGT', 'CTGTATCA', 'CGATGATA', 'CTGCGTGT', 'CGCATACA', 'CGCATGCT', 'CAGCAGCA', 'CTAGCTGA', 'CTGTACGT', 'CTAGCTGT', 'CTGATCGT', 'CATCTACT', 'CTCGAGAT', 'CTAGCGTA', 'CGCTACGA', 'CGTGCAGT', 'CACATACT', 'CATGTGCA', 'CGTGCGCT', 'CAGAGCAT', 'CGATAGAT', 'CTCGATCA', 'CATAGTCA', 'CATATGAT', 'CGTAGCAT', 'CATCGTGA', 'CGCACGCT', 'CTGCTACT', 'CGTGCTCA', 'CACAGTCA', 'CGCTATAT', 'CACACTCT', 'CTGTATGA', 'CGATACGA', 'CAGATACA', 'CACGATAT', 'CTATGTCT', 'CGCGTGAT', 'CATATGCT', 'CTGATAGT', 'CGTGCTGA', 'CGATCTCT', 'CTCTGCAT', 'CTATAGAT', 'CTAGCATA', 'CTGAGTGT', 'CATGAGAT', 'CATGATGA', 'CTATACGA', 'CTAGTAGA', 'CTACGCAT', 'CATGCGCT', 'CTATCTGA', 'CTATCATA', 'CTGTAGAT', 'CGATAGCT', 'CTGTATGT', 'CACGTGAT', 'CTAGCGCT', 'CAGACTGT', 'CAGCGTCA', 'CTAGTGAT', 'CTGATGCA', 'CGCGTATA', 'CTACATAT', 'CTGCGCAT', 'CATATGTA', 'CATAGTCT', 'CGATGTAT', 'CATATGCA', 'CGCTGCAT', 'CTCGAGTA', 'CGTAGAGA', 'CATGAGTA', 'CGCGCAGT', 'CGCGAGCT', 'CTCGCTCA', 'CGTGCATA', 'CTGAGCGT', 'CTGCGTCA', 'CTGCTGCA', 'CGCGTACT', 'CAGCGTAT', 'CGATGCGA', 'CTAGACGT', 'CGATCGTA', 'CTAGCTCT', 'CTGAGTGA', 'CAGCTGCA', 'CTGTGAGT', 'CTGCTCGT', 'CTCGTACT']

cfc_inv = {}
def create_cfc_inv():
    global cfc_inv
    for i,c in enumerate(cfc):
        cfc_inv[c] = i

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
        enc_strand = [ cfc[s] for s in strand ]
        return "".join(enc_strand)

    def _decode_cfc(self, s):
        global cfc_inv
        if cfc_inv.has_key(s):
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
            else:
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
        key = base_conversion.convertBytesToInt(dec[0:self._keyWidthInBytes])        
        return key,"".join([ chr(x) for x in dec[self._keyWidthInBytes:]])
