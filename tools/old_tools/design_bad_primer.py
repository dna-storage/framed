#!/usr/bin/python

from nupack.mfe import *
from nupack.complexes import *
from dnastorage.primer.primer_util import *
import string
import argparse
import editdistance as ed

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def checkComplexes(seqs):
    prefix = create_mfe_input(seqs)    
    complexes(complexes_args(prefix))
    c = read_mfe_output(prefix+".ocx-mfe")
    problems = [cc for cc in c if cc['deltaG'] < -10.0]
    return problems

def chanceOfBonding():
    count = 0
    for i in range (0,10):
        p1 = id_generator(size=20,chars="AGCT")
        p2 = id_generator(size=20,chars="AGCT")
        p2 = reverse_complement(p2)
        c = checkComplexes([p1,p2])
        if len(c) > 0:
            count += 1
            #print c
    print "Bonding: {}".format(100.0*float(count)/10.0)


def int_list(arg):
    l = arg.split(",")
    i = [ int(ll) for ll in l ]
    return i

if __name__ == "__main__":
    import sys

    parser = argparse.ArgumentParser(description="Create primers with high likelihood of false binding.")
    parser.add_argument('--N',type=int,dest="N",action="store",default=100, help=".")    
    parser.add_argument('--o',dest="o",action="store",default=None, help="Output file.")    
    parser.add_argument('--primer-file',dest="primer_file",action="store",default=None, help="List of primers.")
    parser.add_argument('--primer',dest="primer",action="store",default=None, help="Single primer to use as starting point.")
    parser.add_argument('--range',type=int_list,dest="range",action="store",default=[2,20], help="Range of Hamming distance to bad primers")
    parser.add_argument('--step',type=int,dest="step",action="store",default=1, help="Step between of Hamming distances")
    parser.add_argument('--use-nupack',dest="use_nupack",action="store_true",default=True,help="Perform analysis using nupack.")

    args = parser.parse_args()

    #chanceOfBonding()
    data = []
    primers = []
    if args.primer_file != None:
        primers = read_primers(args.primer_file)        
        if args.primer != None:
            print "Ignoring primer ",args.primer, " and using input file only"
    elif args.primer != None:
        primers.append(args.primer)
        
    bad_primers = []
    perStep = args.N/((args.range[1]-args.range[0])/args.step)
    
    for i in range(args.range[0],args.range[1],args.step):
        for p in primers:
            for j in range(0,perStep):
                while True:
                    m = mutate_sequence(p,i)
                    
                    if ed.eval(p,m) != i:
                        continue
                        
                    if m in bad_primers:
                        continue

                    assert not m in p
                    #m = create_anti_sequence(p)
                    #if args.use_nupack:
                    #    c1 = checkComplexes([reverse_complement(p),m])
                    #    c2 = checkComplexes([reverse_complement(m),p])                        
                    if checkTm(m,(51.5,52.5)):
                        bad_primers.append(m)
                        break
                    #    break
                    #else:
                    #    bad_primers.append(m)
                    #    break

    if args.o != None:
        out = open(args.o,"w")
        for p in bad_primers:
            out.write(p+"\n")
        out.close()
    #else:
    print "Input: ",args.primer,getTm(args.primer)
    for p in bad_primers:
        #print p,"\t",getTm(p),"\t",hamming_distance(p,args.primer)
        print hamming_distance(p,args.primer),",", p,",", getTm(p)
