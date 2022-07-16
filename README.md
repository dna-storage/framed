# dnastorage
![testing](https://github.ncsu.edu/dna-based-storage/dnastorage/actions/workflows/makefile.yml/badge.svg)

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [License](#license)
- [Issues](https://github.ncsu.com/dna-based-storage/dnastorage/issues)

# Overview

Core encoding, decoding, and file manipulation support for modeling DNA-based information storage systems.

# Documentation

As documentation for the software becomes available, it will be placed under the docs folder.

# System Requirements

## Hardware Requirements
The dnastorage module requires only a standard computer with enough RAM and compute power to support the needed operations. However, encoding or decoding large files may perform poorly, necessitating a more capable system.

## Software Requirements
### OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:

+ macOS: Catalina 10.15.3
+ Linux: Ubuntu 18.04.3

Note that most OSes will support our software by using Docker. You will also need:
- git (optional but recommended)
- Python
- C++ compiler
- pip, python package installer


### Python Dependences

Our code has been tested on python versions 3.6 to 3.8. It has the following dependences:

```
nose
sphinx
biopython
editdistance
statistics
matplotlib
numpy
scipy
Pillow
python-Levenshtein
joblib
editdistance
bitarray
pytest
```

# Installation Guide

If you already have python 3 installed on your system, the simplest thing to do is download or checkout the code from GitHub.  Then, run the following commands:

    git clone [this repository] dnastorage
    cd dnastorage
    
I recommend making a virtual environment first, but this is optional:

    python -m venv venv
    source venv/bin/activate
    
Then, use pip to install the requirements and packages:

    pip install -r requirements.txt

To install dnastorage package for local development:

    python setup.py develop

Check to make sure it's (mostly) working properly. These tests need some further development to thoroughly validate the infrastructure:

    pytest new_tests

   
# License

This software is released under the LGPLv3 license.

# Acknowledgment

This work was supported by the National Science Foundation (Grants CNS-1650148, CNS-1901324, ECCS 2027655) and a North Carolina State University Research and Innovation Seed Funding Award (Grant 1402-2018-2509).




