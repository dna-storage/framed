[![DOI](https://zenodo.org/badge/611912849.svg)](https://zenodo.org/badge/latestdoi/611912849)

# SDC Project

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [FrameD Analysis](#running-framed-analysis)
- [License](#license)

# Overview

Core encoding, decoding, and file manipulation support for modeling DNA-based information storage systems as published at https://academic.oup.com/bioinformatics/article/39/10/btad572/7274858?utm_source=advanceaccess&utm_campaign=bioinformatics&utm_medium=email.

# Documentation

In depth documentation can be found in the [wiki](https://github.com/dna-storage/framed/wiki).

# System Requirements

## Hardware Requirements

The dnastorage module requires only a standard computer with enough RAM and compute power to support the needed operations. However, encoding or decoding large files may perform poorly, necessitating a more capable system. The infrastructure to encode, decode, and analyze sequencing data will perform best if using an HPC system, or system with considerble number of cores at its disposal, with an MPI implementation installed.

**Given the nature of evaluating dnastorage systems, we not only leverage MPI parallelism, but we also rely on batch level parallelism. Our infrastructure currently assumes an IBM LSF managed platform to launch many independent jobs. If your system is different, e.g. Cobalt, Slurm, etc you will need to write in your support for launching such jobs by extending the TcshJob class in [submit.py](tools/lsf/submit.py) directory.**

**For simplicity, we also support background jobs run in a sub-shell as a quick way of testing a small set of encoding/decoding excersizes.**

## Software Requirements
### OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:

+ Linux: CentOS 7

Note that most OSes will support our software by using Docker. **You will also need the following before completing the installation steps**:
- git 
- C shell
- C++ compiler
- pip, python package installer
- conda, package management
- MPI implementation
	- tested with Intel MPI Library for Linux OS, Version 2017 Update 1 Build 20161016
 	- tested with open MPI Library for Linux OS, Version 4.1.5
- Julia, currently tested with 1.6.2 (2021-07-14)


### Python Dependences

Our code has been tested on python versions 3.6 to 3.10. A working environment with all dependencies that can be installed via conda can be found at [.yaml file](dnastorage.yml). This installation is handled automatically, as demostrated in the following guide on installation.

# Installation Guide

## Installing natively 

Given that the previous dependencies have been installed and loaded, the first step is to get the code from GitHub and initialize the enviroment with additional dependencies with the following commands.
   
    tcsh 
    git clone [this repository] framed
    cd framed
    source dnastorage.env
    make init

To install dnastorage package for local development:
    
    make develop

**Note: activating the environment and loading environment variables should be done every new user session.**


## Installing via Docker Image

You can build our analysis environemnt that is based on FrameD by running the command 

 	docker build -t framed-image .

No further installation steps should be required, and you should have the FrameD analysis suite installed within the ``framed-image`` image. 

# Running FrameD Analysis 

# Non-HPC Example Fault Injection

We recognize that HPC systems may not be immediately available, so we include 2 small examples that can be run to excersize the pipelines and error models evaluated in the FrameD paper on commodity hardware. We do this by spawning subshell background processes that run job scripts to emulate the HPC-based job submission process and output. These examples can be found in [examples/small](examples/small). The following command will generate two shell jobs, with a total of 12 cores being used.

    cd $DNASTORAGE_HOME
    tcsh examples/small/run_small.csh

Within the core limitations of the user's hardware, this approach can be used to perform larger simulations if desired and only requires an implementation of MPI to be installed on the machine.

## Compiling Fault Injection Results

After the fault injection experiments are run, there will be a number of statistics and json files in each leaf directory. In these leaves, you will also find a human readable "fi.stats" file containing statistics gathered during fault injection. To make analysis more interesting by comparing the results of different run configurations, you can compile all the data together for a given file that was fault injected with the following commands.

    cd $DNASTORAGE_HOME
    cd tools/data_analysis
    python db_gen.py --path <top result path> --name fault_injection_dataframe
 
The final command generates a Pandas dataframe, and stores it at the top directory path of the collection of results, e.g. <top result path>. This data frame can be loaded and analyzed any way in which you prefer. An example of anlayzing data frames can be found in the notebook [framed.ipynb](notebooks/framed.ipynb) which produces the figures of the paper using raw data collected using `db_gen.py`.



# License

This software is released under the LGPLv3 license.

# Acknowledgment

This work was supported by the National Science Foundation (Grants CNS-1650148, CNS-1901324, ECCS 2027655) and a North Carolina State University Research and Innovation Seed Funding Award (Grant 1402-2018-2509).






