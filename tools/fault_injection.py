from dnastorage.fi import fault_injector
from dnastorage.codec import dense
from dnastorage.codec import commafreecodec
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec.rscodec import *
from dnastorage.arch.strand import *
from dnastorage.primer.primer_util import *
from dnastorage.util.neg_binomial_gen import *
from dnastorage.handle_strands.strand_handlers import *
import sys
import os

#!/usr/bin/python

def build_decode_architecture(arch, pf, primer5, primer3, fountain_table=None):
    if arch == "UW+MSv1":
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 20
        dec = DecodeXORStrand(pf,p)
        return dec

    elif arch == "Illinois":
        illini = illinois.IllinoisCodec(primer5,150,Checksum())
        p = StrandPrimers(primer5, primer3, illini)
        pf.packetSize = 149
        enc = DecodeNaiveStrand(pf,p)
        return enc

    elif arch == 'Goldman':
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 5
        dec = DecodeGoldmanStrand(pf,4,p)
        return dec

    elif arch == 'Binary':
        b = binary.BinaryRotateCodec(binary.BinaryCodec(Checksum()))
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 17
        enc = DecodeNaiveStrand(pf,p)
        return enc

    elif arch == 'Fountain':
        #assert False and "Not fully implemented"
        h = huffman.RotateCodec(huffman.HuffmanCodec(23,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 22
        enc = DecodeFountainStrand(pf,fountain_table,p)    
        return enc

    elif arch == 'NCState':
        assert False and "Not fully implemented"
        return None

    elif arch == 'RS+CFC8':
        b = commafreecodec.CommaFreeCodec(13,None,2)
        p = StrandPrimers(primer5, primer3, b)
        dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        return dec


def read_header(dec_file):
    f = dec_file
    header = {}
    while True:
        l = f.readline()
        if l.startswith('%') and 'bytes encoded' in l:
            words = l[1:].split(' ')
            try:
                size = int(words[0])
                header['size'] = size
            except:
                continue
        elif l.startswith('%') and "- primer 5' end" in l:
            words = l[1:].split(' ')
            header['primer5'] = words[0]
        elif l.startswith('%') and "- primer 3' end" in l:
            words = l[1:].split(' ')
            header['primer3'] = words[0]
        elif not l.startswith('%'):
            break

    return header

#function that will add arguments to args by reading the header if needed, also returns a clean packetizedFile, clean Decoder object and a set of clean strands
#required for determining correctness when compared to faulty runs   
def clean_run(args): 
    if args.input_file != sys.stdin:
        h = read_header(args.input_file)

    if args.filesize == 0:
        if h.has_key('size'):
            args.filesize = h['size']
        else:
            print "Please provide the file's size."
            sys.exit(-1)

    if args.primer3 == None: 
        if h.has_key('primer3'):
            args.primer3 = h['primer3']
        else:
            print "Please provide primer3."
            sys.exit(-1)

    if args.primer5 == None:
        if h.has_key('primer5'):
            args.primer5 = h['primer5']
        else:
            print "Please provide primer3."
            sys.exit(-1)                

    table = []        
    if args.arch == 'Fountain':
        if args.input_file == sys.stdin:
            tfile = open("default.table","r")
        else:
            tfile = open(args.input_file.name+".table","r")
        if tfile == None:
            print "Please fountain table."
            sys.exit(-1)                                            
        while True:
            l = tfile.readline()
            if l == "":
                break
            x = [int(y) for y in l.split(',')]
            table.append( [x[0],x[1:]] )
        tfile.close()
        #print table

    clean_packetizedFile = WritePacketizedFilestream(args.o,args.filesize,20)
    
    clean_Decoder = build_decode_architecture(args.arch, clean_packetizedFile, args.primer5, args.primer3,table)

    ii = 0
    #initial_strands will hold the initial strands of the input file before fault injection 
    initial_strands=[]
    while not clean_Decoder.complete:
        s = args.input_file.readline()
        ii += 1
        if len(s) == 0:
            break
        s = s.strip()
        if s.startswith('%'):
            continue
        
        initial_strands.append(s)
        clean_Decoder.decode(s)  

    return clean_Decoder,clean_packetizedFile,initial_strands


#use the negative binomial distribution to generate a new pool of strands
def distribute_reads(strands,neg_bin_ranomizer):
    read_count=[]
    new_pool=[]
    for s in strands:
        read_count.append(neg_bin_ranomizer.gen())
    for index, s in enumerate(strands):
        new_pool+=[s for i in range(0,read_count[index])]
    return new_pool


def build_strand_handler(strand_handler,Codec):
    if strand_handler == "data_vote_simple":
        return data_vote_simple(Codec)






'''
Main file for injecting faults into an input DNA file, 3 available options are available

miss_strand --- missing strand fault model, e.g. strand is not available to the decoder
strand_fault --- faults within strand fault model, e.g. insert deletions, insertions, substitutions
combo --- combine both missing strands and within strand fault model

'''



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Inject faults into input file and perform analysis")
    parser.add_argument('--fault_rate_file', dest="fault_file",action="store", default=None, help='file to have fault rate per nucleotide')
    parser.add_argument('--fault_model',required=True,choices=['miss_strand','strand_fault','combo'])
    parser.add_argument('--run',action="store_true",help="Should strand_fault faults be in a run")
    parser.add_argument('--missing_count',dest="missing",action="store",type=int,default=10,help="Number of missing strands")
    parser.add_argument('--faulty_count',dest="faulty",action="store",type=int,default=10,help="Number of faulty strands")    
    parser.add_argument('--fail_count',dest="fails",action="store",type=int,default=2,help="Number of faults within faulty strands")
    parser.add_argument('--filesize',type=int,dest="filesize",action="store",default=0, help="Size of file to decode.")
    parser.add_argument('--primer5',dest="primer5",action="store",default=None, help="Beginning primer.")
    parser.add_argument('--primer3',dest="primer3",action="store",default=None, help="Ending primer.")
    parser.add_argument('--primer1_length',dest="p1",action="store",type=int,default=20,help="Length of primer1 (total, so if hierarchy set to 40 instead of 20)")    
    parser.add_argument('--primer2_length',dest="p2",action="store",type=int,default=20,help="Length of primer2, for all our experiments this will be set to 20")    
    parser.add_argument('--seq_mean',dest="mean",action="store",type=float,default=0,help="Mean of the sequencing distribution")    
    parser.add_argument('--seq_var',dest="var",action="store",type=float,default=4.3,help="Variance of the sequencing distribution") 
    parser.add_argument('--simulation_runs',dest="num_sims",action="store",type=int,default=10000,help="Number of simulations to run")    
    parser.add_argument('--o',nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="Output file.")    
    parser.add_argument('--arch',required=True,choices=['UW+MSv1','Illinois','Binary','Goldman','Fountain','RS+CFC8'])
    parser.add_argument('input_file', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Clean strands file') 
    parser.add_argument('--strand_handler',required=True,choices=['simple','data_vote_simple','nuc_vote_simple'])

    
    args = parser.parse_args()

    model_name= args.fault_model

    #set up fault injection arguments
    injector_args=fault_injector.fault_injector_arguments()
    injector_args.o=args.o
    injector_args.input_file=None
    injector_args.fault_file=args.fault_file
    injector_args.fault_model=args.fault_model
    injector_args.run=args.run
    injector_args.missing=args.missing
    injector_args.faulty=args.faulty
    injector_args.fails=args.fails
    injector_args.p1=args.p1
    injector_args.p2=args.p2

    # check if the selected fault model class exists
    if not hasattr(fault_injector, model_name):
        print "Couldn't run fault model '%s'. Class '%s' doesn't exist in fault_injector.py" % \
        (args.fault_model, model_name)
        sys.exit(1)
    
    #clean_Decoder will buffer the clean_packetizedFile we will use to compare with the faulty file
    clean_Decoder,clean_packetizedFile,clean_strands=clean_run(args)

    # instantiate a fault model object
    fault_model = eval('fault_injector.'+model_name+'(injector_args)')
    #instantiate the negative binomial random number generator, used for distributing the number of reads
    read_randomizer=neg_bin(args.mean,args.var)


    #instantiate the strand handler
    strand_handler=build_strand_handler(args.strand_handler,clean_Decoder._Codec)
    #FIX ME: need a method for collecitng data

    #data=data_helper()

    #will need to add 2 more loops to sweep across numbers of faults and numbers of strands

    desired_missing_count=args.missing
    desired_faulty_count=args.faulty
    
    #run many simulations to get an average 
    for sim_number in range(0,args.num_sims):
        #build a "dirty" decoder and packetizedFile for each run
        dirty_packetizedFile= WritePacketizedFilestream(args.o,args.filesize,20)
        table=[]
        dirty_Decoder = build_decode_architecture(args.arch, dirty_packetizedFile, args.primer5, args.primer3,table)


        #get a new pool of strands based on the negative binomial distribution
        multiple_strand_pool=distribute_reads(clean_strands,read_randomizer)

    


        # set the number of strands accordingly for the fault model class
        if len(multiple_strand_pool)<desired_faulty_count or len(multiple_strand_pool)<desired_missing_count:
            fault_model.set_faulty_missing(len(multiple_strand_pool),len(multiple_strand_pool))
        else:
            fault_model.set_faulty_missing(desired_faulty_count,desired_missing_count)
            

        
        #set the fault_model's library equal to the new multiple strand pool
        fault_model.set_library(multiple_strand_pool)
        #inject the faults
        strands_after_faults=fault_model.Run()
        #call the strand handler, will return either key,value pairs or strands
        processed_strands=strand_handler.process(strands_after_faults)
        #Hand the processed strands off to decoding
        for proc in processed_strands:
            if type(proc) is tuple:
               dirty_Decoder.decode(None,bypass=True,input_key=proc[0],input_value=proc[1])
            else:
                dirty_Decoder.decode(proc)
        #if there is a final decoding try, try it
        if hasattr(dirty_Decoder,'attempt_final_decoding'):
            dirty_Decoder.attempt_final_decoding()
        
