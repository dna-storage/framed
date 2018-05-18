#!/usr/bin/python

from nupack.mfe import *
from nupack.complexes import *
from dnastorage.primer.primer_util import *
import string
import argparse

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def checkComplexes(seqs):
    prefix = create_mfe_input(seqs,2)    
    complexes(complexes_args(prefix))
    c = read_mfe_output(prefix+".ocx-mfe")
    problems = [cc for cc in c if cc['deltaG'] < -10.0]
    return problems

def chanceOfBonding(N):
    count = 0
    for i in range (0,N):
        p1 = id_generator(size=20,chars="AGCT")
        p2 = id_generator(size=20,chars="AGCT")
        #p2 = reverse_complement(p2)
        c = checkComplexes([p1,reverse_complement(p2)])
        if len(c) > 0:
            count += 1
        else:
            c = checkComplexes([p2,reverse_complement(p1)])
            if len(c) > 0:
                count += 1
            #print c
    return 100.0*float(count)/10.0

def probabilityOfBonding(primers,N,ed):
    data = []
    plen = len(primers[0])
    for i in range(1,plen-3):
        count = 0
        trials = 0
        for p in primers:
            for j in range(0,N):
                if ed:
                    m = mutate_sequence_ed(p,i)
                else:
                    m = mutate_sequence(p,i)
                assert not m in p
                c = checkComplexes([reverse_complement(p),m])
                if len(c) > 0:
                    #print c
                    count += 1                    
                else:
                    c = checkComplexes([reverse_complement(m),p])
                    if len(c) > 0:
                        #print c
                        count += 1                                        
                trials+=1

        data.append( [i,count,trials] )
        print data
    return data


if __name__ == "__main__":
    import sys

    parser = argparse.ArgumentParser(description="Analyze primers.")
    parser.add_argument('--o',dest="o",action="store",default=None, help="Output file.")    
    parser.add_argument('--base-primer',dest="base_primer",action="store",default=None, help="Base primer to use for comparison.")

    parser.add_argument('--primer-file',dest="primer_file",action="store",default=None, help="List of primers.")


    args = parser.parse_args()
    
    primers = []
    if args.primer_file != None:
        primers = read_primers(args.primer_file)        
    else:
        print "No primers specified."
        sys.exit(0)

    details = []
    if args.base_primer != None:
        H = { x:0 for x in range(0,20) }
        E = { x:0 for x in range(0,20) }
        bH = { x:0 for x in range(0,20) }
        bE = { x:0 for x in range(0,20) }
        for p1 in primers:
            e = ed.eval(p1,args.base_primer);
            h = hamming_distance(p1,args.base_primer);
            H[h] += 1
            E[e] += 1
            c = checkComplexes([reverse_complement(p1),args.base_primer])
            if len(c) > 0:
                bH[h] += 1
                bE[e] += 1                
                details.append([p1,h,e,1,GC(p1),getTm(p1)])
            else:
                details.append([p1,h,e,0,GC(p1),getTm(p1)])
            print " {} {} - {}".format(h,e,p1);
        
        for d in details:            
            print "\t".join([str(x) for x in d])            

        for x in range(0,20):
            if H[x] > 0:
                print "{}\t{}".format(x,float(bH[x])/H[x]*100)

        for x in range(0,20):
            if E[x] > 0:
                print "{}\t{}".format(x,float(bE[x])/E[x]*100)
            
    else:
        all = []
        for p1 in primers:
            l = []
            for p2 in primers:
                e = ed.eval(p1,p2);
                h = hamming_distance(p1,p2);
                l.append([e,h])
            all.append(l)
            print l

