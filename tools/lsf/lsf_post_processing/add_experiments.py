"""
Tool to add experiments to sequencing data.
"""


import os
import pickle
import random
import re
import distutils.dir_util
import shutil

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="add experiments to sequencing data")
    parser.add_argument('-p','--sequencing_experiment_top_path', default=None,required=True,help="target directory for Framed experiments")
    parser.add_argument('-e','--regex',default="",type=str,required=True,help="regular expression to look for, creates a copy of an already existing experiment if the substring matches a directory in the sequencing experiment data path and also matches a fastq file name input to this tool which represents a new experiment.")
    parser.add_argument('-t','--target_directory',default="",type=str,required=True,help="path to target directory containing fastqs")
    args = parser.parse_args()

    assert os.path.exists(args.target_directory) and os.path.exists(args.sequencing_experiment_top_path)

    target_directory_abs = os.path.abspath(args.target_directory)


    experiments=[]
    for s in os.listdir(args.sequencing_experiment_top_path):
        if os.path.isdir(os.path.join(args.sequencing_experiment_top_path,s)):
            experiments.append(os.path.join(args.sequencing_experiment_top_path,s))


    for f in os.listdir(target_directory_abs):
        if ".fastq" in f:
            fastq_name=f
            print("Making new directory for {}".format(fastq_name))
            fastq_path=os.path.join(target_directory_abs,f)
            substring_match = re.search(args.regex,fastq_name)
            if substring_match is None: continue
            substring=substring_match[0]
            for s in experiments:
                if os.path.basename(s)==fastq_name: continue #copy from yourself is weird if redoing
                if substring in os.path.basename(s):
                    #we are going to create a new directory that mirrors the experiment and change the top name to the fastq file name
                    new_directory = os.path.join(args.sequencing_experiment_top_path,fastq_name)
                    print("Making new directory {}".format(new_directory))
                    os.makedirs(new_directory,exist_ok=True)
                    root,files,dirs=next(os.walk(new_directory,topdown=False))
                    for fd in dirs: #remove some old sym links 
                        if os.path.islink(os.path.join(root,fd)):# or os.path.isfile(os.path.join(root,fd)):
                            print("Remove {}".format(os.path.join(root,fd)))
                            os.remove(os.path.join(root,fd))

                    print("Copy from source {} to destination {}".format(s,new_directory))
                    shutil.copytree(s,new_directory,symlinks=True,dirs_exist_ok=True)
                    root,files,dirs = next(os.walk(new_directory,topdown=False))
                    assert os.path.islink(os.path.join(root,"sequencing_data_link"))
                    print("Removing sym link {}".format(os.path.join(root,"sequencing_data_link")))
                    os.remove(os.path.join(root,"sequencing_data_link"))
                    os.symlink(fastq_path,os.path.join(root,"sequencing_data_link"))
                    break
