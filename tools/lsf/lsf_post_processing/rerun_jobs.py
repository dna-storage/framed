import os

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Script to generate reruns")
    parser.add_argument('--input',type=str,required=True,action="store",help="Path to file with paths to rerun,paths should be absolute")
    args = parser.parse_args()

    with open(args.input,'r') as f:
        for line in f:
            p = line.strip()
            if os.path.isfile(p):
                p=os.path.dirname(p)
            print("Running in {}".format(p))

            os.chdir(p)
            os.system('bsub < hpc_lsfjob.csh')
