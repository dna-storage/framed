# dnastorage
![testing](https://github.ncsu.edu/dna-based-storage/dnastorage/actions/workflows/makefile.yml/badge.svg)


Python modules and C++ tools to support modeling and simulation of a DNA-based information storage system.

## Getting Started

### Dependencies

To use this software, you need:
- git (optional but recommended)
- Python
- C++ compiler
- pip, python package installer

Linux computers have this software by default. Macs and Windows may not have them. (There are a variety of ways to get Python, more details coming soon.)

### Installation Steps for Python Tools

1. Clone or download the repository into a working directory

      1. cd /full/path/to/some/directory
      2. git clone https://github.ncsu.edu/jtuck/dnastorage.git

2. Run make to install package dependencies and build nupack software

      - make init

3. Set two environment variables. You will need to do this each time you want to use this software. Or, add these commands to your shell profile to have them run automatically each time you open terminal.

      - export PATH="/full/path/to/some/directory/dnastorage/other_software/nupack3.0.6/bin:$PATH"
      - export PYTHONPATH="/full/path/to/some/directory/dnastorage/"

  
4. Run tests to make sure it works okay. (Note, still adding more meaningful tests.  This may pass even though something isn't installed properly.)

      - make test

### Installation steps for C++ DNA Storage System Simulator

1. Run the makefile in system_sim directory to get simulator dependencies.

    1. cd system_sim (resulting path should be "/full/path/to/some/directory/dnastorage/system_sim")
    2. make init

2. Run the makefile in system_sim directory to compile the simulator
   2. make

3. To delete the executable and object files run the makefile in the system_sim directory
   1. make clean


The compiled executable will be in the directory "/full/path/to/some/directory/dnastorage/dnastorage/system_sim/bin" with the name "system_sim".
The compiled object files will be in the directory "/full/path/to/some/directory/dnastorage/dnastorage/system_sim/build".
The source code for the system simulator can be found in the directory "/full/path/to/some/directory/dnastorage/dnastorage/system_sim/src".
