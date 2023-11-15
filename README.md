[![DOI](https://zenodo.org/badge/611912849.svg)](https://zenodo.org/badge/latestdoi/611912849)

# SDC Project

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Running Sequencing Analysis](#running-sequencing-analysis)
- [License](#license)

# Overview

Core encoding, decoding, and file manipulation support for modeling DNA-based information storage systems.

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

We will include a pre-built docker image at https://hub.docker.com/repository/docker/kvolkel/framed/ when it becomes available.

In the meantime, with this branch, you can try to build an image yourself and run it with the following 2 commands

	docker build -t sdc_submission .
 	docker run -it sdc_submission

When a pre-built image becomes available, you will be able to run the following command to run the container 

	docker run -it kvolkel/framed:sdc_submission

This should set up an appropriate environment in which the code for this project can at least be run, although it may not be tractible unless on an HPC system.

# Running Sequencing Analysis 
TODO: Finish this section
## Get Raw Data

## Commands For Performing Sequencing Analysis


**NOTE: This will only be tractable on HPC systems, and at the moment these commands will only work on HPC systems with an LSF scheduler. Extending to other schedulers, like SLURM, should be relatively easy by following the LSFJob class in [sumbit.py](tools/lsf/lsf_utils/submit.py).**

## Non-HPC Examples

We recognize that HPC systems may not be immediately available, so we include some small examples that analyze small portions of the reads that we have mapped from our sequencing runs. 

TODO: Finish this section


# License

This software is released under the LGPLv3 license.

# Acknowledgment

This work was supported by the National Science Foundation (Grants CNS-1650148, CNS-1901324, ECCS 2027655) and a North Carolina State University Research and Innovation Seed Funding Award (Grant 1402-2018-2509).




