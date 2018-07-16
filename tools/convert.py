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


def build_encode_architecture(arch, filename, primer5, primer3):    
    if arch == "UW+MSv1":
        h = huffman.HuffmanRotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        enc = EncodeXORStrand(filename,20,p)
        return enc

    elif arch == "Illinois":
        illini = illinois.IllinoisCodec(primer5,150,Checksum())
        p = StrandPrimers(primer5, primer3, illini)
        enc = EncodeNaiveStrand(filename,149,p)
        return enc

    elif arch == 'Goldman':
        h = huffman.HuffmanRotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        enc = EncodeGoldmanStrand(filename,5,4,p)
        return enc

    elif arch == 'Fountain':
        #assert False and "Not fully implemented"
        h = huffman.HuffmanRotateCodec(huffman.HuffmanCodec(23,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        enc = EncodeFountainStrand(filename,22,1.5,p)    
        return enc

    elif arch == 'Binary':
        b = binary.BinaryRotateCodec(binary.BinaryCodec(Checksum()))
        p = StrandPrimers(primer5, primer3, b)
        enc = EncodeNaiveStrand(filename,17,p)
        return enc

    elif arch == 'NCState':
        assert False and "Not fully implemented"
        return None

def build_decode_architecture(arch, filename, filesize, primer5, primer3, fountain_table=None):
    if arch == "UW+MSv1":
        h = huffman.HuffmanRotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        dec = DecodeXORStrand(filename,filesize,20,p)
        return dec

    elif arch == "Illinois":
        illini = illinois.IllinoisCodec(primer5,150,Checksum())
        p = StrandPrimers(primer5, primer3, illini)
        enc = DecodeNaiveStrand(filename,filesize,149,p)
        return enc

    elif arch == 'Goldman':
        h = huffman.HuffmanRotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        dec = DecodeGoldmanStrand(filename,filesize,5,4,p)
        return dec

    elif arch == 'Binary':
        b = binary.BinaryRotateCodec(binary.BinaryCodec(Checksum()))
        p = StrandPrimers(primer5, primer3, b)
        enc = DecodeNaiveStrand(filename,filesize,17,p)
        return enc

    elif arch == 'Fountain':
        #assert False and "Not fully implemented"
        h = huffman.HuffmanRotateCodec(huffman.HuffmanCodec(23,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        enc = DecodeFountainStrand(filename,filesize,22,fountain_table,p)    
        return enc

    elif arch == 'NCState':
        assert False and "Not fully implemented"
        return None

def read_header(dec_file):
    f = open(dec_file,"r")
    header = {}
    while True:
        l = f.readline()
        if l.startswith('%') and 'file size' in l:
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
    f.close()
    return header


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create primers with high likelihood of false binding.")

    parser.add_argument('--o',dest="o",action="store",default="a.out", help="Output file.")    
    parser.add_argument('--encode',dest="encode",action="store_true",default=False,help="Encode input file into DNA.")    
    parser.add_argument('--decode',dest="decode",action="store_true",default=False,help="Decode input file from DNA into original binary form.")    

    parser.add_argument('--arch',required=True,choices=['UW+MSv1','Illinois','Binary','Goldman','Fountain','NCState'])
    parser.add_argument('--filesize',type=int,dest="filesize",action="store",default=0, help="Size of file to decode.")

    parser.add_argument('--primer5',dest="primer5",action="store",default=None, help="Beginning primer.")
    parser.add_argument('--primer3',dest="primer3",action="store",default=None, help="Ending primer.")
    parser.add_argument('input_file', nargs=1, help='file to be converted')
    args = parser.parse_args()
    
    # Debbuging: did I get the args?
    #print args.o
    #print args.input_file
    #print args.primer5

    if args.encode and args.decode:
        print "Can't encode and decode at the same time."
        print "Assuming encode and continuing."
    
    if len(args.input_file) < 1:
        print "No input file specified."
        sys.exit(-1)

    if args.encode:
        enc = build_encode_architecture(args.arch, args.input_file[0], args.primer5, args.primer3)
        ofile = open(args.o,"w")
        ofile.write("%{}\n".format(args.input_file[0]))
        ofile.write("%{} - file size\n".format(enc.fileSize))
        ofile.write("%{} - primer 5' end\n".format(args.primer5))
        ofile.write("%{} - primer 3' end\n".format(args.primer3))        
        for e in enc:
            ofile.write("{}\n".format(e))
        ofile.close()

        if args.arch == 'Fountain':
            t = enc.getTable()
            tablefile = open(args.o+".table","w")
            for x,y in t:
                entries = [ str(e) for e in y ]
                entries = ",".join(entries)
                tablefile.wriemate("{},{}\n".format(x,entries))
            tablefile.close()
                            
    elif args.decode:
        h = read_header(args.input_file[0]) 
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
            tfile = open(args.input_file[0]+".table","r")
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

        Decoder = build_decode_architecture(args.arch, args.o, args.filesize, args.primer5, args.primer3,table)
        strands = read_primers(args.input_file[0])
        for s in strands:
            Decoder.decode(s)
        Decoder.write()
            
