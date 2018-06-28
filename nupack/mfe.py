#!/usr/bin/python
import os
import tempfile

def create_mfe_input(sequences,num):
    if len(sequences) == 0:
        return None
    (fd,filename) = tempfile.mkstemp(suffix=".in",prefix=sequences[0][0:5],text=True)
    os.close(fd)
    prefix = filename[:filename.find(".")]
    f = open(filename,"w")
    f.write( "{}\n".format(len(sequences)) )
    for s in sequences:
        f.write( "{}\n".format(s) )
    f.write( "{}\n".format(num) )
    f.close()
    return prefix

def read_mfe_output(filename):
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    i = 0
    complexes = []
    l = lines
    while i < len(l):
        if "% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %" in l[i]:
            i+=1
            # mfe file will not have a name, but an ocx-mfe file will have one
            if l[i][0] == '%':
                name = l[i+1][2:].strip()
                i+=1
            else:
                name = "unknown"
            try:
                order = int(l[i].strip())
                deltaG = float(l[i+1].strip())
                pattern = l[i+2].strip()
                i += 3
            except:
                order = 0
                deltaG = -100
                pattern = ""
            pairs = []
            while not("% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %" in l[i]):
                p = [ int(x) for x in l[i].split() ]
                pairs.append(p)
                i+=1
                
            # capture information about complexes found
            complexes.append( { "name": name, 
                                "order": order,
                                "deltaG": deltaG,
                                "pattern":pattern,
                                "pairs":pairs } )                
        i += 1
    return complexes
    

if __name__ == "__main__":
    import sys
    c = read_output(sys.argv[1])
    
    for cc in c:
        if "((((((((((((((((((((+))))))))))))))))))))" in cc["pattern"]:
            continue
        if cc["deltaG"] < -10.0:
            print "{} ({}) -- {}".format(cc["name"],cc["deltaG"],cc["pattern"])

    print c
