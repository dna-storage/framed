#!/usr/bin/python

from dnastorage.codec import dense
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.arch.strand import *
from dnastorage.primer.primer_util import *
import sys
import os

#read through the library and get the strands
def read_library(dna_library):
    return_library=[]
    with open(dna_library, 'rb') as lib:
        strand=lib.readline()
        strand=strand.strip('\n')
        #make sure that the strands are 200 in length
        assert len(strand) == 200
        return_library.append(strand)
        while strand:
            strand=lib.readline()
            if strand:
                strand=strand.strip('\n')
                assert len(strand)==200
                return_library.append(strand)
    lib.close()
    return return_library

#strip off the beginning and end primers of DNA strands
def strip_primers(dna_strands):
    return_strands=[]
    for strand in dna_strands:
        return_strands.append(strand[20:-20])
    return return_strands

#returns the smallest hamming distance between the input sub strand and sub strands within the full strand under consideration 
def hamming(sub_strand,full_strand,sub_strand_start,sub_strand_home,full_ID):
    min=-1
    for full_strand_start in range(140):
        #continue if the sub strand in the full strand is the current input sub strand
        if sub_strand_home == full_ID and sub_strand_start==full_strand_start:
            assert sub_strand==full_strand[full_strand_start:full_strand_start+20]
            continue
        distance=sum(c1!=c2 for c1,c2 in zip(sub_strand,full_strand[full_strand_start:full_strand_start+20]))

        if min==-1 or min>distance:
            min=distance
            
    return min
    



#this code searches a library of strands in order to find strand sequences that do not repeat for attacking data


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Search for 20 bp strings in library to find unique attack regions")

    parser.add_argument('--o',dest="o",action="store",default="attack.out", help="Output file.")     
    parser.add_argument('input_library', nargs=1, help='library to be searched')
    parser.add_argument('--lower_range',dest='low',type=int, action="store", default=0, help="Lower Bound of strands to search")
    parser.add_argument('--upper_range',type=int,dest='high', action="store", default=6000, help="Upper Bound of strands to search")
    parser.add_argument('--max_secret_size',dest='max_secret',action='store',default=70,type=int,help='Max size of the secret region')
    
    args = parser.parse_args()
        
    if len(args.input_library) < 1:
        print "No input library specified."
        sys.exit(-1)
        
    dna_library=[]
    dna_library=read_library(args.input_library[0])
  

    #search through the dna library and look for 2 data regions while keeping the resultant product greater than 100
    # and less than 150 bp. Make sure that the target data strands are 20 bp. Make sure that no other 
    dna_library=strip_primers(dna_library)
    assert len(dna_library[0]) == 160

    
    #dictionary to keep count for each data strand pairs reoccurance, we will select the one that re-occurs less
    strand_count={}
    strand_hamming_A={}
    strand_hamming_B={}
    data_len=20
    strand_ID=0
    #look through all the strands and find sub strand pairs that are not in other strands, keep a count
    for strand_ID,  strand in enumerate(dna_library[args.low:args.high]):
        #print strand_ID
        secret_length_start=args.max_secret
        #A start and B start is where the two sub strings start for a particular strand 
        A_start_beginning=0
        B_start_beginning=A_start_beginning+secret_length_start+20

        while B_start_beginning+20<=160:
            secret_length_start=args.max_secret
            B_start=B_start_beginning
            while secret_length_start>=20:
                secret_length_search=secret_length_start
                A_start=A_start_beginning
                while secret_length_search>=20:
                    A_substrand=strand[A_start:A_start+20]
                    B_substrand=strand[B_start:B_start+20]
                    #create a unique key for each substrand choice
                    strand_key=str(strand_ID+args.low)+':'+str(A_start)+':'+str(B_start)
                    #make sure we do not double search a key
                    if strand_key not in strand_count:
                        strand_count[strand_key]=0
                        #look through all strands
                        hamming_A_dist={}
                        hamming_B_dist={}
                        for search_ID, search_strand in enumerate(dna_library):
                            reverse_strand=search_strand[::-1]
                            #calculate hamming distance
                            min_A=hamming(A_substrand,search_strand,A_start,strand_ID+args.low,search_ID)
                            min_B=hamming(B_substrand,search_strand,B_start,strand_ID+args.low,search_ID)
                            min_A_rev=hamming(A_substrand,reverse_strand,None,None,None)
                            min_B_rev=hamming(B_substrand,reverse_strand,None,None,None)
                            #build up a distribution for hamming distances
                            if min_A not in hamming_A_dist:
                                hamming_A_dist[min_A]=1
                            else:
                                hamming_A_dist[min_A]=hamming_A_dist[min_A]+1
                            if min_B not in hamming_B_dist:
                                hamming_A_dist[min_B]=1
                            else:
                                hamming_B_dist[min_B]=hamming_B_dist[min_B]+1
                            if min_A_rev not in hamming_A_dist:
                                hamming_A_dist[min_A_rev]=1
                            else:
                                hamming_A_dist[min_A_rev]=hamming_A_dist[min_A_rev]+1
                            if min_B_rev not in hamming_B_dist:
                                hamming_B_dist[min_B_rev]=1
                            else:
                                hamming_B_dist[min_B_rev]=hamming_B_dist[min_B_rev]+1
                                
                                
                            #update the dictionary to track the smallest hamming distance for both substrands of a candidate
                            if strand_key not in strand_hamming_A or strand_hamming_A[strand_key]>min_A :
                                strand_hamming_A[strand_key]=min_A
                            if strand_key not in strand_hamming_B or strand_hamming_B[strand_key]>min_B:
                                strand_hamming_B[strand_key]=min_B

                            if strand_hamming_A[strand_key]>min_A_rev:
                                strand_hamming_A[strand_key]=min_A_rev
                            if strand_hamming_B[strand_key]>min_B_rev:
                                strand_hamming_B[strand_key]=min_B_rev

                            
                            
                            #dont count the considered sub-strand as a repitition
                            if search_ID == strand_ID:
                                strand_count[strand_key]=strand_count[strand_key]+(search_strand.count(A_substrand)-1)+(search_strand.count(B_substrand)-1)+reverse_strand.count(A_substrand)+reverse_strand.count(B_substrand)
                            elif search_ID != strand_ID:
                                strand_count[strand_key]=strand_count[strand_key]+search_strand.count(A_substrand)+search_strand.count(B_substrand)+reverse_strand.count(A_substrand)+reverse_strand.count(B_substrand)

                        if strand_count[strand_key] == 0:
                            print '\n'
                            print strand_key
                            print strand_hamming_A[strand_key]
                            print strand_hamming_B[strand_key]
                            print '\n'
                            print 'A hamming Distribution\n'
                            for key in hamming_A_dist:
                                print '{} ====== {}\n'.format(key,hamming_A_dist[key])
                            print 'B hamming Distribution\n'
                            for key in hamming_B_dist:
                                print '{} ====== {}\n'.format(key,hamming_B_dist[key])
                        #print strand_key
                        #print strand_count[strand_key]
                    A_start=A_start+1
                    secret_length_search=secret_length_search-1
                B_start=B_start-1
                secret_length_start=secret_length_start-1
            A_start_beginning=A_start_beginning+1
            B_start_beginning=B_start_beginning+1
            

            

    #write out to file the sorted dictionary
    out_results=open(args.o,'w')
    sorted_strands=sorted(strand_count, key=lambda x: strand_count[x])
    for strand_names in sorted_strands:
        #split up key to print A and B substrings
        split_key=strand_names.split(':')
        A=dna_library[int(split_key[0])][int(split_key[1]):int(split_key[1])+20]
        B=dna_library[int(split_key[0])][int(split_key[2]):int(split_key[2])+20]
        out_results.write("{} ---- {} ---- A: {}   B: {}   A_HAM: {}    B_HAM: {} \n".format(strand_names,strand_count[strand_names],A,B,strand_hamming_A[strand_names],strand_hamming_B[strand_names]))
 
