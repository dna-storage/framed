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
        #print data
    return data


if __name__ == "__main__":
    import sys

    parser = argparse.ArgumentParser(description="Create primers with high likelihood of false binding.")
    parser.add_argument('--N',type=int,dest="N",action="store",default=20, help="Trials per hamming distance value.")    
    parser.add_argument('--o',dest="o",action="store",default=None, help="Output file.")    
    parser.add_argument('--primer-file',dest="primer_file",action="store",default=None, help="List of primers.")
    parser.add_argument('--primer',dest="primer",action="store",default=None, help="Single primer to analyze.")

    parser.add_argument('--sweep',dest="sweep",action="store_true",default=True,help="Analyze likelihood of bonding w.r.t. to primers.")   
    parser.add_argument('--ed',dest="ed",action="store_true",default=False, help="Use edit distance rather than hamming distance during sweep.")

    parser.add_argument('--mc',dest="mc",action="store_true",default=False,help="Analyze likelihood of bonding between two random primers using Monte Carlo method.")    
    parser.add_argument('--csv',dest="csv",action="store_true",help="Dump to screen of file in CSV format")

    args = parser.parse_args()
    
    if args.mc:
        pr = chanceOfBonding(args.N)
        print "Chance of bonding = {}\%".format(pr)
    elif args.sweep:    
        primers = []
        if args.primer_file != None:
            primers = read_primers(args.primer_file)        
            if args.primer != None:
                print "Ignoring primer ",args.primer, " and using input file only"
        elif args.primer != None:
            primers.append(args.primer)
        
        data = probabilityOfBonding(primers,args.N,args.ed)
        if args.csv:
            s = [ "{},{}\n".format(d[0],d[1]/float(d[2])) for d in data ]
            s = "".join(s)
        else:
            s = [ "{}\t{}\n".format(d[0],d[1]/float(d[2])) for d in data ]
            s = "".join(s)
        if args.o:
            ofile = open(args.o,"w")
            ofile.write(s)
            ofile.close()
        else:
            print s
