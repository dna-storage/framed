#!/usr/bin/python
from base_conversion import *
from primer_util import *
import random
from nextera import *
from brds_ecc import addECC

def create_codes(base,length):
    codes = [convertBase(base,x,length) for x in range(0,base**length)]
    return codes
    
def code_filter(codes):
    if len(codes)==0:
        return codes
    chunks = { x:1 for x in nextera_chunks(len(codes[0])) }
    new_codes = []
    for c in codes:
        if chunks.has_key(c):
            continue
        if hasRepeat(c):
            continue
        if not (c[-1] == 'G' or c[-1] == 'C'):
            continue
        if (c[0] == 'G' or c[0] == 'C'):
            continue
        new_codes.append(c)
    return new_codes

def digital_sum(code,alphabetScore):
    return sum(alphabetScore[s] for s in code)

def compute_brds_table(symbols,codes,sums,alphabetScore):
    sums.sort()
    print sums
    if len(sums) * len(symbols) > len(codes):
        print "Warning: succcessful table unlikely!"
    min_sums = min(sums)
    max_sums = max(sums)
    brds_table = {}
    for s in symbols:
        brds_table[s] = {}
        for m in sums:
            brds_table[s][m] = ""
    
    partition = {}
    for c in codes:
        s = digital_sum(c,alphabetScore)
        if not partition.has_key(s):
            partition[s] = []
        partition[s].append(c)
    
    #print "Partition"
    #print partition 

    h = {}

    for s in symbols:
        prev = ""
        for m in sums:
            if prev!="" and brds_table[s][m] == "":
                ds = digital_sum(prev,alphabetScore)
                if m+ds >= min_sums and m+ds <= max_sums:
                    brds_table[s][m] = prev

            if brds_table[s][m] == "":
                found = False
                for p in sums:
                    val = p - m
                    if partition.has_key(val) and len(partition[val])>0:
                        brds_table[s][m] = partition[val].pop(random.randint(0,len(partition[val])-1))
                        if not h.has_key(val):
                            h[val] = 0
                        h[val] = h[val]+1
                        found = True
                        break
                if found:
                    prev = brds_table[s][m]
            
    print brds_table
    print h
    for s in symbols:
        for m in sums:
            if brds_table[s][m] == "":
                return None

    return brds_table                        


if __name__ == "__main__":
    print digital_sum("AGCT",{'A':-1,'C':1,'G':1,'T':-1})

    score = {'A':-1,'C':1,'G':1,'T':-1}
    even = [-4,-2,2, 4]
    odd = [-2,-1,1,2]

    symbols = [x for x in range(0,256)]
    #codes = [convertQuarnary(x,7) for x in range(0,4**7)]

    #for k in range(3,9):
    #    print "len({})={}".format(k,len(code_filter(create_codes(4,k))))

    codes = create_codes(4,5)
    #codes = code_filter(create_codes(4,7))
    #ecc_codes = code_filter([addECC(x,4) for x in create_codes(4,4)])

    #print len(ecc_codes[0])

    #print len(codes), len(ecc_codes)

    table = compute_brds_table(symbols,codes,odd,score)
    print table
