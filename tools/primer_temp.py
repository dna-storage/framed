#!/usr/bin/python
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import random
import string
import subprocess
import re
import os
import argparse

from dnastorage.primer.nextera import *
from nupack.mfe import *
from nupack.complexes import *
from dnastorage.primer.primer_util import *

import sys

primers = read_primers(sys.argv[1])

for p in primers:
    t = mt.Tm_NN(p)
    print "{}\t{}\t{}".format(p,t,hamming_distance(p,"CGTGGCAATATGACTACGGA"))

print "P3 Tm=",mt.Tm_NN("CAGGTACGCAGTTAGCACTC")
print "P4 Tm=",mt.Tm_NN("CGTGGCAATATGACTACGGA")
