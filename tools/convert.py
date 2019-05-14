#!/usr/bin/python

from dnastorage.codec import dense
from dnastorage.codec import norepeatscodec
from dnastorage.codec import commafreecodec
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec.rscodec import *
from dnastorage.arch.strand import *
from dnastorage.arch.builder import *
from dnastorage.primer.primer_util import *

import sys
import os

def read_header(dec_file):
    f = dec_file
    header = {}
    while True:
        # ugly hack to restore seek position if we fail to match a header line
        tell = f.tell()
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
            f.seek(tell) # restore position to beginning of previous read
            break

    return header


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Encode and decode files.")

    parser.add_argument('--o',nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="Output file.")
    parser.add_argument('--encode',dest="encode",action="store_true",default=False,help="Encode input file into DNA.")
    parser.add_argument('--decode',dest="decode",action="store_true",default=False,help="Decode input file from DNA into original binary form.")

    parser.add_argument('--arch',required=True,choices=available_file_architectures())
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
