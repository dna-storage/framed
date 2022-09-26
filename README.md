# dnastorage
![testing](https://github.ncsu.edu/dna-based-storage/dnastorage/actions/workflows/makefile.yml/badge.svg)

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Quick Start](#quick-start)
- [License](#hardware-requirements)
- [Issues](https://github.ncsu.com/dna-based-storage/dnastorage/issues)

# Overview

Core encoding, decoding, and file manipulation support for modeling DNA-based information storage systems.

# Documentation

As documentation for the software becomes available, it will be placed under the docs folder.

# System Requirements

## Hardware Requirements

The dnastorage module requires only a standard computer with enough RAM and compute power to support the needed operations. However, encoding or decoding large files may perform poorly, necessitating a more capable system. The infrastructure to encode, decode, and analyze sequencing data will perform best if using an HPC system with an MPI implementation installed.

**Given the nature of evaluating dnastorage systems, we not only leverage MPI parallelism, but we also rely on batch level parallelism. Our infrastructure currently assumes an IBM LSF managed platform to launch many independent jobs. If your system is different, e.g. Cobalt, you will need to write in your support for launching such jobs in the [lsf](tools/lsf) directory.**

## Software Requirements
### OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:

+ macOS: Catalina 10.15.3
+ Linux: Ubuntu 18.04.3
+ Linux: CentOS 7

Note that most OSes will support our software by using Docker. You will also need:
- git (optional but recommended)
- Python
- C++ compiler
- pip, python package installer
- conda, package management
- MPI implementation, currently tested with Intel MPI Library for Linux OS, Version 2017 Update 1 Build 20161016


### Python Dependences

Our code has been tested on python versions 3.6 to 3.10. A working environment with all dependencies that can be installed via conda can be found at [.yaml file](dnastorage.yaml).

# Installation Guide

If you already have python 3 installed on your system, the simplest thing to do is download or checkout the code from GitHub.  Then, run the following commands:

    git clone [this repository] dnastorage
    cd dnastorage
    
I recommend making a conda environment with the required dependencies:

    conda env create -f dnastorage.yaml
    conda activate <path to conda environment>
    
Make sure to source environment variables for the project:
    source dnastorage.env

To install dnastorage package for local development:

    python setup.py develop

Check to make sure it's (mostly) working properly. These tests need some further development to thoroughly validate the infrastructure:

    pytest new_tests

Note: activating the environment and loading environment variables should be done every new user session.

# Quick Start 

## Running Fault Injection Jobs

More detail on how to write new encoding/decoding pipelines and how to configure them will be more extensively covered in the Wiki documentation for this project. The following is a quick start to running a fault injection job. 

First open the example configuration in a text editor:
    
    emacs tools/lsf/examples/fault_injection/basic_hedges.json
    
Make a change to the following line in the example json file, this "file" option should be changed to a path to a file on your system you want to encode/decode:

    "file":"path/to/file"
    
Once this is changed, perform the following commands to run fault injection experiments:
    
    cd $DNASTORAGE_HOME
    cd tools/lsf/
    python generate_fi_jobs.py --params examples/fault_injection/basic_hedges.json --memory 8 --cores 100 --core_depth 1 --dump_dir `pwd` --queue standard --time 4
    
This command will generate a number of jobs, 1 for each encoding and fault injection enviroment combination. The command specifies a run time of 4 hours, using 100 cores for fault injection parallelization, and 8 GB of memory for each node. The results will be placed in a directory with a name related to the file name used in the previous step. From there you will find a directory tree where each leaf in the tree is a directory containing results for each unique run
 

## Compiling Fault Injection Results

After the fault injection experiments are run, there will be a number of statistics and json files in each leaf directory. In these leaves, you will also find a human readable "fi.stats" file containing statistics gathered during fault injection. To make analysis more interesting by comparing the results of different run configurations, you can compile all the data together for a given file that was fault injected with the following commands.

    cd $DNASTORAGE_HOME
    cd tools/data_analysis
    python fault_injection_db_gen.py --path <top result path> --name fault_injection_dataframe
 
The final command generates a Pandas dataframe, and stores it at the top directory path of the collection of results, e.g. <top result path>. This data frame can be loaded and analyzed any way in which you prefer.


# License

This software is released under the LGPLv3 license.

# Acknowledgment

This work was supported by the National Science Foundation (Grants CNS-1650148, CNS-1901324, ECCS 2027655) and a North Carolina State University Research and Innovation Seed Funding Award (Grant 1402-2018-2509).




