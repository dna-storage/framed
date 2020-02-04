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
import filecmp

import logging
logger = logging.getLogger('dna')
logger.setLevel(logging.DEBUG)
_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
_ch = logging.FileHandler("dna.debug.log",mode='w')
_ch.setFormatter(_formatter)
logger.addHandler(_ch)

if __name__ == "__main__":


    import argparse

    parser = argparse.ArgumentParser(description="Test all formats for encoding and decoding DNA files.")


    parser.add_argument('--primer5',dest="primer5",action="store",default="AAAA", help="Beginning primer.")
    parser.add_argument('--primer3',dest="primer3",action="store",default="AAAA", help="Ending primer.")
    parser.add_argument('input_file', nargs='?', default="", help='file to be converted')

    args = parser.parse_args()


    fails = 0
    succeeds = 0
    
    formats = file_system_formats()
    for f in formats:
    
        tmp_fd,tmpname = tempfile.mkstemp()
        os.close(tmp_fd)
        # don't need this

        orig_tmp_fd,orig_tmp = tempfile.mkstemp()

        try:
            wf = DNAFile.open(tmpname, "w", \
                              args.primer5,\
                              args.primer3,\
                              format_name=f)

            with open(args.input_file,"rb") as input_file:
                while True:
                    b = input_file.read(1000)
                    if len(b)==0:
                        break
                    wf.write(b)
            wf.close()
        
            rf = DNAFile.open(tmpname,"r", \
                              primer5=args.primer5,\
                              primer3=args.primer3)

            while True:
                b = rf.read(1000)
                if len(b)==0:
                    break
                os.write(orig_tmp_fd,b)
            os.close(orig_tmp_fd)            
                                
        except Exception as e:
            print (e)
            print ("Exception when using {}.".format(f))
            
        if filecmp.cmp(args.input_file, orig_tmp,shallow=False)==False:
            print ("{} failed.".format(f))
            fails += 1
        else:
            print ("{} succeeded.".format(f))
            succeeds += 1

        # Remove tmp files
        try:
            os.remove(tmpname)
        except:
            print ("Failed to remove {}".format(tmpname)) 
            os.remove(tmpname)
        try:
            os.remove(orig_tmp)
        except:
            print ("Failed to remove {}".format(orig_tmp)) 
            os.remove(orig_tmp)

            
    print ("{} failures, {} successes.".format(fails,succeeds))
