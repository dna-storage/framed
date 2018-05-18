#!/usr/bin/python
import os
import subprocess
import sys

def complexes_args(filename):
    #return ['complexes','-T','52','-material','dna','-ordered','-pairs',filename]
    return ['complexes','-material','dna','-ordered','-pairs','-mfe',filename]

def complexes(args):    
    devnull = open(os.devnull,"w")
    result = subprocess.call(args,stdout=devnull,stderr=devnull,shell=False)
    devnull.close()
    return result

def read_ocx_output(filename):
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    i = 0
    complexes = []
    l = lines

    output = {}
    strands = {}
    
    # look for strands ids
    found = False
    while i < len(l):
        if "% Number of strands: " in l[i]:
            strand_count = int(l[i][21:])
            i+=2
            for ii in range(i,i+strand_count):
                ll = l[ii][1:].split()
                strands[ int(ll[0].strip()) ] = ll[1].strip()
            found = True
            i += strand_count
            output['strands'] = strands
            break
        i += 1
    
    if "% T = " in l[i]:
        Temp = float(l[i].split()[3])
        output['Temp'] = Temp
        i += 1        

    output['complexes'] = {}
    while i < len(l):
        #print l[i]
        ll = l[i].split()
        ll2 = [int(x) for x in ll[:-1]]
        ll2.append(float(ll[-1]))        
        indexes = [ x-1 for x in range(2,len(ll2)-1) if ll2[x]==1 ]
        indexes = []
        for x in range(2,len(ll2)-1):
            if ll2[x]==2:
                indexes.append(x-1)
                indexes.append(x-1)
            elif ll2[x]==1:
                indexes.append(x-1)
        
        d = {}
        d['structure'] = ll2[0]
        #print indexes
        d['order'] = sum([ len(strands[j]) for j in indexes ])
        d['indexes'] = indexes
        d['sequences'] = [ strands[j] for j in indexes ]
        d['deltaG'] = ll2[-1]        
        output[ 'complexes' ][ ll2[0] ] = d
        i += 1
    
    #print output
    return output


if __name__ == "__main__":
    from mfe import *
    #filename = "TTGCGGAATTTAATCCCGGG2"
    #r = complexes(complexes_args(filename))
    
    prefix = create_mfe_input(['TTGCGGAATTTAATCCCGGG',
                               'CCCGGGATTAAATTCCGCAA'])
    
    complexes(complexes_args(prefix))
    c = read_mfe_output(prefix+".ocx-mfe")
    print [cc for cc in c if cc['deltaG'] < 10.0]
