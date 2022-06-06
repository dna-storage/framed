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
from dnastorage.system.formats import file_system_formats, file_system_encoder_by_abbrev
from dnastorage.system.formats import file_system_formatid_by_abbrev, file_system_abbrev
from dnastorage.system.dnafile import DNAFile 
from dnastorage.system.pipeline_dnafile import DNAFilePipeline
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


def get_dnafile_class(formatid):
    if formatid >= 0x2000 or formatid < 0x100:
        return DNAFile
    else:
        #print ("Use pipeline for {}".format(file_system_abbrev(formatid)))
        return DNAFilePipeline


if __name__ == "__main__":


    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser(description="Test all formats for encoding and decoding DNA files.")


    parser.add_argument('--primer5',dest="primer5",action="store",default="AAAA", help="Beginning primer.")
    parser.add_argument('--primer3',dest="primer3",action="store",default="AAAA", help="Ending primer.")
    parser.add_argument('input_file', nargs='?', default="", help='file to be converted')

    args = parser.parse_args()


    stats = {}
    
    fails = 0
    succeeds = 0
    
    formats = file_system_formats()
    #print (formats)
    
    for f in formats:        
        if file_system_encoder_by_abbrev(f) == None:
            continue

        id = file_system_formatid_by_abbrev(f)
        if id < 0x20:
            continue

        DNAFilecls = get_dnafile_class(id)
        
        stats[f] = {}
        stats[f]['encoded'] = "Fail"
        stats[f]['decoded'] = "Fail"
        stats[f]['matched'] = "Fail"
        stats[f]['dnafile'] = DNAFilecls.__name__        
        
                
        tmp_fd,tmpname = tempfile.mkstemp()
        os.close(tmp_fd)
        # don't need this

        orig_tmp_fd,orig_tmp = tempfile.mkstemp()

        header_fd,header_data_path = tempfile.mkstemp()
        #print (header_data_path)
        os.close(header_fd);
        
        encoder_params = {
            "primer3":args.primer3,
            "primer5":args.primer5
        }

        
        try:
            if DNAFilecls == DNAFile:
                wf = DNAFilecls.open(tmpname, "w", \
                                  args.primer5,\
                                  args.primer3,\
                                  format_name=f)
            else:
                wf = DNAFilecls.open(tmpname, "w", \
                                     encoder_params=encoder_params,\
                                     format_name=f,\
                                     fsmd_abbrev='FSMD-Pipe',\
                                     fsmd_header_filename=header_data_path)                

            with open(args.input_file,"rb") as input_file:
                while True:
                    b = input_file.read(1000)
                    if len(b)==0:
                        break
                    wf.write(b)
            wf.close()

            
            
            stats[f]["encoded"] = "Success"
            #print ("Encode {} succeeded.".format(f))

            if DNAFilecls == DNAFile:
                rf = DNAFilecls.open(tmpname,"r", \
                                  primer5=args.primer5,\
                                primer3=args.primer3)
            else:
                rf = DNAFilecls.open(tmpname,"r", \
                                     encoder_params=encoder_params,\
                                     fsmd_abbrev='FSMD-Pipe',\
                                     fsmd_header_filename=header_data_path)
                if rf == None:
                    raise DNAStorageError(msg="Couldn't make a {}.".format(DNAFilecls.__name__))
                
            while True:
                b = rf.read(1000)
                if len(b)==0:
                     break
                # else:
                #     #print (b)
                #     #print ("got some data!", len(b))
                os.write(orig_tmp_fd,b)
                
            os.close(orig_tmp_fd)

            #print (os.stat(orig_tmp))
            
            stats[f]["decoded"] = "Success"
                                
        except Exception as e:
            #print (type(e),e)
            print ("Exception when using {}.".format(f))

        try:
            if filecmp.cmp(args.input_file, orig_tmp,shallow=False)==False:
                # print ("{} failed.".format(f))
                # print (orig_tmp)
                fails += 1
                stats[f]["matched"] = "Fail"
                #break
            else:
                # print ("{} succeeded.".format(f)) #
                succeeds += 1
                stats[f]["matched"] = "Success"                
                                
        except:
            stats[f]["matched"] = 'Exception'
            print ("Exception {} failed.".format(f))
            fails += 1
                
        # Remove tmp files
        try:
            os.remove(tmpname)
            os.remove(orig_tmp)
            os.remove(header_data_path)
        except:
            pass
            #print ("Failed to remove {}".format(tmpname)) 


    df = pd.DataFrame.from_dict(stats)
    print (df)
    print ("{} failures, {} successes.".format(fails,succeeds))
    
