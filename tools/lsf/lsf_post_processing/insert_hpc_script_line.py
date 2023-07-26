#Script to allow new lines to be inserted into a set of files, useful for adding more parameters to an hpc script

import os

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Script to generate reruns")
    parser.add_argument('--input',type=str,required=True,action="store",help="Path to file with paths to rerun,paths should be absolute")
    parser.add_argument('--line_number',type=int,required=True,action="store",help="line number to put injected string")
    parser.add_argument('--s',type=str,required=True,action="store",help="line to inject")
    parser.add_argument('--t',type=str,required=True,choices=["insert","sub"])
    args = parser.parse_args()
    with open(args.input,'r') as f:
        for line in f:
            p = line.strip()
            if os.path.isfile(p):
                p=os.path.dirname(p)
            p = os.path.join(p,"hpc_lsfjob.csh")
            print("Injecting in {}".format(p))
            assert os.path.exists(p)
            lines=[]
            with open(p,"r") as f: lines = f.readlines()
            if args.t=="insert": lines.insert(args.line_number,"{}\n".format(args.s))
            elif args.t=="sub": lines[args.line_number]="{}\n".format(args.s)
            #print(lines)
            with open(p,"w+") as f:
                for l in lines: f.write(l)
            
