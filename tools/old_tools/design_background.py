#This script aims to create a 160 bp background strand in the aims that it is not similar to any of the strands in the library
from dnastorage.codec import dense
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.arch.strand import *
from dnastorage.primer.primer_util import *
import sys
import os
import math

complement_table={'A':'T',
                  'G':'C',
                  'C':'G',
                  'T':'A'}

background_strand=''

def compute_complement(strand):
    return_strand=[]
    for nucleotide in strand:
        return_strand.append(complement_table[nucleotide])
    #join the array and return it
    return ''.join(return_strand)




#scans the strands and notes what sequences are already used
def scan_strands(library, length):
    for strand in library:
        start=0
        while start+length<=len(strand):
            sub_strand=strand[start:start+length]
            nucleotide_set[sub_strand]='used'
            start=start+1

    


#read through the library and get the strands
def read_library(dna_library):
    return_library=[]
    with open(dna_library, 'rb') as lib:
        strand=lib.readline()
        strand=strand.strip('\n')
        return_library.append(strand)
        while strand:
            strand=lib.readline()
            if strand:
                strand=strand.strip('\n')
                #need to append the 4 strands of concern
                return_library.append(strand)
                reverse=strand[::-1]
                return_library.append(reverse)
                complement=compute_complement(strand)
                return_library.append(complement)
                reverse_complement=compute_complement(reverse)
                return_library.append(reverse_complement)
    lib.close()
    return return_library



#no_repeat_table={'A':['G','C','T'],
  #               'G':['A','C','T'],
 #                'C':['A','G','T'],
   #              'T':['A','G','C']
    #             }



no_repeat_table={'A':['G','C','T'],
                 'G':['A','C','T'],
                 'C':['A','G','T'],
                 'T':['A','G','C']
                 }

nucs=['A','G','C','T']

#recursively add on to end of strand
def recursive_strand(strand,max_len,nuc_set):
    if(len(strand)+1==max_len):
        for n in no_repeat_table[strand[len(strand)-1]]:
            #save all of the found strands
            nuc_set[strand+n]='clear'
    else:
        for n in no_repeat_table[strand[len(strand)-1]]:
            recursive_strand(strand+n,max_len,nuc_set)
        

    
def create_set(max_run,nuc_set):
    return_set={}
    for n in nucs:
        #send the beginning 
        recursive_strand(n,max_run,nuc_set)



#extract the good strands from the dictionary
def find_clear(input_set):
    output_set=[]
    for key in input_set:
        if input_set[key]=='clear':
            output_set.append(key)
    return output_set



#recursively build all combinations of the good strands
def search_strands(set_of_strands,strand,max_run):
    global background_strand
    if len(set_of_strands)==1:
        _strand=strand+set_of_strands[0]
        if'AA' not in _strand and 'GG' not in _strand and 'CC' not in _strand and 'TT' not in _strand and search_candidates(_strand,max_run):
            print ' Full strand '
            background_strand=strand+set_of_strands[0]
            return True
        else:
            return False
        
    else:
        for n,  s in enumerate(set_of_strands):
            _strand=strand+s
            _set_of_strands=set_of_strands[:]
            del _set_of_strands[n]
            if 'AA' not in _strand and 'GG' not in _strand and 'CC' not in _strand and 'TT' not in _strand:
                if(search_candidates(_strand,max_run)):
                    if(search_strands(_set_of_strands,_strand,max_run)):
                        return True
    return False
            
    

#function to find candidate background strands 
def find_background(good_strands,max_run):
    output_strand_list=[]
    num_subs=160/max_run
    start=0
    while start+num_subs<=len(good_strands):
        strand_set=good_strands[start:start+num_subs]
        if(search_strands(strand_set,'',max_run)):
            break
        start=start+1
   


#search for candidates that give us what we need
def search_candidates(candidate,max_len):
    strand=candidate
    #print strand
    #search the strand to make sure its ok
    start=0
    while start+max_len<=len(strand):
        if nucleotide_set[strand[start:start+max_len]] != 'clear':
            return False
        else:
            start=start+1
            
    return True
        

def pad_background(max_run):
    global background_strand
    pad_length=160-len(background_strand)
    pad_set={}
    create_set(pad_length,pad_set)
    for pad_seq in pad_set:
        strand=background_strand+pad_seq
        if 'AA' in strand or 'GG' in strand or 'CC' in strand or 'TT' in strand:
            continue
        else:
            if search_candidates(strand,max_run):
                background_strand=strand
                print 'padded'
                return True
            else:
                continue
    return False



        
if __name__ == "__main__":
    import argparse
    #dictionary to keep track of what non repeating strands are used in the library, will be used to meet
    #the requirement that our background strand does not repeat this
    nucleotide_set={}


    parser = argparse.ArgumentParser(description="Design a 160 bp background strand")
    parser.add_argument('--o',dest="o",action="store",default="background.out", help="Output file.")     
    parser.add_argument('input_library', nargs=1, help='library to be searched')
    parser.add_argument('--max_run',dest='max_run',type=int, action="store", default=0, help="Max # of nucleotides to match")
    
    
    args = parser.parse_args()
        
    if len(args.input_library) < 1:
        print "No input library specified."
        sys.exit(-1)
        
    dna_library=[]
    #read the DNA library with all of the strands 
    dna_library=read_library(args.input_library[0])
    create_set(args.max_run,nucleotide_set)


    #scan through every 5 lenth strand in the library, fill in the sub strand table as it goes
    scan_strands(dna_library,args.max_run)


    #find the strings that do not repeat    
    good_strands=find_clear(nucleotide_set)


    #try to find background strand
    find_background(good_strands,args.max_run)


    #pad out background strand if it is less than 160 bp
    if background_strand == '':
        background_strand='No available background strand'
    elif len(background_strand)<160:
        if not pad_background(args.max_run):
            background_strand='No available background strand'
        

    out_results=open(args.o,'w')
    out_results.write("Background Strand ---- {}\n".format(background_strand))
 
