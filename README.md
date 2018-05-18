# dnastorage

Python module to support modeling and simulation of a DNA-based information storage system.

## Getting Started

### Dependencies

To use this software, you need:
- Python
- C++ compiler
- pip, python package installer

Linux computers have this software by default. Macs and Windows may not have them. (There are a variety of ways to get Python, more details coming soon.)

### 

1. Clone or download the repository into a working directory

  1. cd /full/path/to/some/directory
  2. git clone https://github.ncsu.edu/jtuck/dnastorage.git

2. Run make to build nupack  

  - make init

3. Set two environment variables. You will need to do this each time you want to use this software. Or, add these commands to your shell profile to have them run automatically each time you open terminal.

  - export PATH="/full/path/to/some/directory/dnastorage/other_software/nupack3.0.6/bin:$PATH"
  - export PYTHONPATH="/full/path/to/some/directory/dnastorage/"

  
4. Run tests to make sure it works okay. (Note, still adding more meaningful tests.  This may pass even though something isn't installed properly.)

   - make test