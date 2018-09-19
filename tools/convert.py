#!/usr/bin/python

from dnastorage.codec import dense
from dnastorage.codec import commafreecodec
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec.rscodec import *
from dnastorage.arch.strand import *
from dnastorage.primer.primer_util import *
import sys
import os

def build_encode_architecture(arch, pf, primer5, primer3):    
    if arch == "UW+MSv1":
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        ## build packetized file encode
        pf.packetSize = 20
        enc = EncodeXORStrand(pf,p)
        return enc

    elif arch == "Illinois":
        illini = illinois.IllinoisCodec(primer5,150,Checksum())
        p = StrandPrimers(primer5, primer3, illini)
        pf.packetSize = 149
        enc = EncodeNaiveStrand(pf,p)
        return enc

    elif arch == 'Goldman':
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 5
        enc = EncodeGoldmanStrand(pf,4,p)
        return enc

    elif arch == 'Fountain':
        #assert False and "Not fully implemented"
        h = huffman.RotateCodec(huffman.HuffmanCodec(23,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 22
        enc = EncodeFountainStrand(pf,1.5,p)    
        return enc

    elif arch == 'Binary':
        b = binary.BinaryRotateCodec(binary.BinaryCodec(Checksum()))
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 17
        enc = EncodeNaiveStrand(pf,p)
        return enc

    elif arch == 'RS+CFC8':
    
        b = commafreecodec.CommaFreeCodec(13,None,2)
        p = StrandPrimers(primer5, primer3, b)
        enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        return enc

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
        #print pf.packetSize
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


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create primers with high likelihood of false binding.")

    parser.add_argument('--o',nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="Output file.")    
    parser.add_argument('--encode',dest="encode",action="store_true",default=False,help="Encode input file into DNA.")    
    parser.add_argument('--decode',dest="decode",action="store_true",default=False,help="Decode input file from DNA into original binary form.")    

    parser.add_argument('--arch',required=True,choices=['UW+MSv1','Illinois','Binary','Goldman','Fountain','RS+CFC8'])
    parser.add_argument('--filesize',type=int,dest="filesize",action="store",default=0, help="Size of file to decode.")

    parser.add_argument('--primer5',dest="primer5",action="store",default="", help="Beginning primer.")
    parser.add_argument('--primer3',dest="primer3",action="store",default="", help="Ending primer.")
    parser.add_argument('input_file', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='file to be converted')
    args = parser.parse_args()
    
    # Debbuging: did I get the args?
    #print args.o
    #print args.input_file
    #print args.primer5

    if args.encode and args.decode:
        print "Can't encode and decode at the same time."
        print "Assuming encode and continuing."
    
    if args.encode:
        packetizedFile = ReadPacketizedFilestream(args.input_file)
        enc = build_encode_architecture(args.arch, packetizedFile, args.primer5, args.primer3)
        strands = []
        #try:
        for e in enc:
            strands.append(e)
        #except:
        #    pass

        ofile = args.o
        ofile.write("%{}\n".format(args.input_file.name))
        if args.input_file == sys.stdin:
            ofile.write("%{} - bytes encoded\n".format(enc.bytes_encoded))
        else:
            ofile.write("%{} - bytes encoded\n".format(packetizedFile.size))
        ofile.write("%{} - primer 5' end\n".format(args.primer5))
        ofile.write("%{} - primer 3' end\n".format(args.primer3))        
        for s in strands:
            ofile.write("{}\n".format(s))
        ofile.close()

        if args.arch == 'Fountain':
            t = enc.getTable()
            if args.o == sys.stdout:
                tablefile = open("default.table","w")
            else:
                tablefile = open(args.o+".table","w")
            for x,y in t:
                entries = [ str(e) for e in y ]
                entries = ",".join(entries)
                tablefile.wriemate("{},{}\n".format(x,entries))
            tablefile.close()
                            
    elif args.decode:
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

        packetizedFile = WritePacketizedFilestream(args.o,args.filesize,20)

        Decoder = build_decode_architecture(args.arch, packetizedFile, args.primer5, args.primer3,table)

        ii = 0
        while not Decoder.complete:
            s = args.input_file.readline()
            ii += 1
            if len(s) == 0:
                break
            s = s.strip()
            if s.startswith('%'):
                continue
            Decoder.decode(s)  

        if args.o == sys.stdout:
            print "".join([ '-' for _ in range(0,80) ])

        Decoder.write()

        if args.o == sys.stdout:
           print ""
        


