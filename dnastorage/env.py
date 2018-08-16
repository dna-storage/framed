#!/usr/bin/python
import os
import csv

def get_home_directory():
    if os.environ.has_key('DNASTORAGEHOME'):
        return os.environ['DNASTORAGEHOME']
    else:
        path =  os.getcwd()
        # assume we are somewhere in the path
        i = path.find('dna-storage-sim')
        if i >= 0:
            return path[:(i+len('dna-storage-sim'))]
        else:
            return path

def get_data_directory():
    d = get_home_directory()
    return os.sep.join([d,"oligos"])

def get_oligo_csv():
    oligo = get_data_directory() + os.sep + "oligos.csv"
    f = open(oligo,"r")
    reader = csv.reader(f)
    oligos = {}
    for row in reader:
        oligos[row[0]] = row[1]
    f.close()
    return oligos

if __name__ == "__main__":
     print get_home_directory()
     print get_data_directory()
