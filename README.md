[![DOI](https://zenodo.org/badge/611912849.svg)](https://zenodo.org/badge/latestdoi/611912849)

# FrameD

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [FrameD Analysis](#running-framed-analysis)
- [License](#license)

# Overview

Core encoding, decoding, and file manipulation support for modeling DNA-based information storage systems.

# Documentation

In depth documentation can be found in the [wiki](https://github.com/dna-storage/framed/wiki).

# System Requirements

## Hardware Requirements

The dnastorage module requires only a standard computer with enough RAM and compute power to support the needed operations. However, encoding or decoding large files may perform poorly, necessitating a more capable system. The infrastructure to encode, decode, and analyze sequencing data will perform best if using an HPC system with an MPI implementation installed.

**Given the nature of evaluating dnastorage systems, we not only leverage MPI parallelism, but we also rely on batch level parallelism. Our infrastructure currently assumes an IBM LSF managed platform to launch many independent jobs. If your system is different, e.g. Cobalt, you will need to write in your support for launching such jobs by extending the TcshJob class in [submit.py](tools/lsf/submit.py) directory.**

## Software Requirements
### OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:

+ macOS: Catalina 10.15.3
+ Linux: Ubuntu 18.04.3
+ Linux: CentOS 7

Note that most OSes will support our software by using Docker. **You will also need the following before completing the installation steps**:
- git 
- C shell
- C++ compiler
- pip, python package installer
- conda, package management
- MPI implementation, currently tested with Intel MPI Library for Linux OS, Version 2017 Update 1 Build 20161016
- Julia, currently tested with 1.6.2 (2021-07-14)


### Python Dependences

Our code has been tested on python versions 3.6 to 3.10. A working environment with all dependencies that can be installed via conda can be found at [.yaml file](dnastorage.yml). This installation is handled automatically, as demostrated in the following guide on installation.

# Installation Guide

Given that the previous dependencies have been installed and loaded, the first step is to get the code from GitHub and initialize the enviroment with additional dependencies with the following commands.
   
    tcsh 
    git clone [this repository] framed
    cd framed
    source dnastorage.env
    make init

To install dnastorage package for local development:
    
    make develop

**Note: activating the environment and loading environment variables should be done every new user session.**

# Running FrameD Analysis 

## Complete Data Generation Configuration

This repository comes with ocnfiguration files that were used to generate the data in the paper. These can be found under [examples/](examples). There are 6 total `.json` files used for the paper, 3 under [examples/iid](examples/iid) for simulating the i.i.d error model, and 3 under [examples/DNArSim](examples/DNArSim) for simulating a nanopore channel using the DNArSim simulator. Each directory contains a C-shell script which calls the job generation tools using the jsons as input, running these will launch all 240 simulations. 

**NOTE: This will only be tractable on HPC systems, and at the moment these commands will only work on HPC systems with an LSF scheduler. Extending to other schedulers, like SLURM, should be relatively easy by following the LSFJob class in [sumbit.py](tools/lsf/lsf_utils/submit.py).**

## Non-HPC Examples

We recognize that HPC systems may not be immediately available, so we include 2 small examples that can be run to excersize the pipelines and error models evaluated in the FrameD paper on commodity hardware. We do this by spawning subshell background processes that run job scripts to emulate the HPC-based job submission process and output. These examples can be found in [examples/small](examples/small). The following command will generate two shell jobs, with a total of 12 cores being used.

    cd $DNASTORAGE_HOME
    tcsh examples/small/run_small.csh

Within the core limitations of the user's hardware, this approach can be used to perform larger simulations if desired and only requires an implementation of MPI to be installed on the machine.

## Compiling Fault Injection Results

After the fault injection experiments are run, there will be a number of statistics and json files in each leaf directory. In these leaves, you will also find a human readable "fi.stats" file containing statistics gathered during fault injection. To make analysis more interesting by comparing the results of different run configurations, you can compile all the data together for a given file that was fault injected with the following commands.

    cd $DNASTORAGE_HOME
    cd tools/data_analysis
    python db_gen.py --path <top result path> --name fault_injection_dataframe
 
The final command generates a Pandas dataframe, and stores it at the top directory path of the collection of results, e.g. <top result path>. This data frame can be loaded and analyzed any way in which you prefer. An example of anlayzing data frames can be found in the notebook [bioinformatics.pynb](notebooks/bioinformatics.pynb) which produces the figures of the paper using raw data collected using `db_gen.py`.

## FrameD Notebooks

Included in this repo is a directory with a Jupyter notebook that was used to generate figures for the manuscript. This notebook is at the path `notebooks/framed.ipynb`. The raw data from simulation outputs that are used to generate figures are at the paths `notebooks/framed_DNArSim.pickle` and `notebooks/framed_iid.pickle`. The notebook can be run by starting in the top directory for the project and running:

	conda activate framed_conda
	cd notebooks
	jupyter-lab framed.ipynb
	
From there run the notebook to run the full analysis.


# License

This software is released under the LGPLv3 license.

# Acknowledgment

This work was supported by the National Science Foundation (Grants CNS-1650148, CNS-1901324, ECCS 2027655) and a North Carolina State University Research and Innovation Seed Funding Award (Grant 1402-2018-2509).




