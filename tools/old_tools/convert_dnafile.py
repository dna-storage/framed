#!/usr/bin/env python
from dnastorage.codec import norepeatscodec
from dnastorage.codec import dense
from dnastorage.codec import commafreecodec
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec.deprecated.rscodec import *
from dnastorage.arch.strand import *
from dnastorage.system.formats import file_system_formats
from dnastorage.system.dnafile import *
#from dnastorage.arch.builder import *
#from dnastorage.primer.primer_util import *

import sys
import os
import tempfile

import logging
logger = logging.getLogger('dna')
logger.setLevel(logging.DEBUG)
_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
_ch = logging.FileHandler("dna.debug.log",mode='w')
_ch.setFormatter(_formatter)
logger.addHandler(_ch)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Encode and decode DNA files.")

    parser.add_argument('--o', dest="o", default="out.dna", help="Output file.")

    parser.add_argument('--encode',dest="encode",action="store_true",default=False,help="Encode input file into DNA.")

    parser.add_argument('--decode',dest="decode",action="store_true",default=False,help="Decode input file from DNA into original binary form.")

    parser.add_argument('--format',required=True,choices=file_system_formats())

    parser.add_argument('--primer5',dest="primer5",action="store",default="AAAA", help="Beginning primer.")
    parser.add_argument('--primer3',dest="primer3",action="store",default="AAAA", help="Ending primer.")
    parser.add_argument('input_file', nargs='?', default="", help='file to be converted')

    args = parser.parse_args()

    if args.encode and args.decode:
        print ("Can't encode and decode at the same time.")
        sys.exit(1)
        
    elif args.encode:

        wf = DNAFile.open(args.o, "w", \
                          args.primer5,\
                          args.primer3,\
                          format_name=args.format)

        with open(args.input_file,"rb") as input_file:
            while True:
                b = input_file.read(1000)
                if len(b)==0:
                    break
                wf.write(b)

        wf.close()
        
        
    elif args.decode:

        rf = DNAFile.open(args.input_file,"r", \
                          args.primer5,\
                          args.primer3)

        with open(args.o,"wb") as of:
            while True:
                b = rf.read(1000)
                if len(b)==0:
                    break
                of.write(b)

    else:
        print ("Neither encode nor decode was selected.")
        print ("No action completed.")
        sys.exit(1)
