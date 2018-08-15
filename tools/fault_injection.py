from dnastorage.fi import fault_injector





'''
Main file for injecting faults into an input DNA file, 3 available options are available

miss_strand --- missing strand fault model, e.g. strand is not available to the decoder
strand_fault --- faults within strand fault model, e.g. insert deletions, insertions, substitutions
combo --- combine both missing strands and within strand fault model

'''



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Inject faults into input file")

    parser.add_argument('--o',dest="o",action="store",default="fi.dna", help="Output file.")    
    parser.add_argument('input_file', nargs=1, help='file to have injections inserted')
    parser.add_argument('--fault_rate_file', dest="fault_file",action="store", default=None, help='file to have fault rate per nucleotide')

    parser.add_argument('--fault_model',required=True,choices=['miss_strand','strand_fault','combo'])

    parser.add_argument('--run',action="store_true",help="Should strand_fault faults be in a run")
  
    
    parser.add_argument('--missing_count',dest="missing",action="store",type=int,default=10,help="Number of missing strands")
    
    parser.add_argument('--faulty_count',dest="faulty",action="store",type=int,default=10,help="Number of faulty strands")    
    parser.add_argument('--fail_count',dest="fails",action="store",type=int,default=2,help="Number of faults within faulty strands")    

    parser.add_argument('--primer1_length',dest="p1",action="store",type=int,default=20,help="Length of primer1 (total, so if hierarchy set to 40 instead of 20)")    
    parser.add_argument('--primer2_length',dest="p2",action="store",type=int,default=20,help="Length of primer2, for all our experiments this will be set to 20")    

    args = parser.parse_args()

    model_name= args.fault_model

    injector_args=fault_injector.fault_injector_arguments()
    injector_args.o=args.o
    injector_args.input_file=args.input_file[0]
    injector_args.fault_file=args.fault_file
    injector_args.fault_model=args.fault_model
    injector_args.run=args.run
    injector_args.missing=args.missing
    injector_args.faulty=args.faulty
    injector_args.fails=args.fails
    injector_args.p1=args.p1
    injector_args.p2=args.p2

    
    # check if the selected fault model class exists
    if not hasattr(fault_injector, model_name):
        print "Couldn't run fault model '%s'. Class '%s' doesn't exist in fault_injector.py" % \
        (args.fault_model, model_name)
        sys.exit(1)

    # instantiate a fault model object
    fault_model = eval('fault_injector.'+model_name+'(injector_args)')

    #run the fault model
    fault_model.Run()
    
    
