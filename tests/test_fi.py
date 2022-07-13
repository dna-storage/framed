
# -*- coding: utf-8 -*-

from .context import dnastorage

import unittest
import statistics
from dnastorage.fi import fault_injector

def check_strand_errors(injection_sites,clean_strands,final_strands):
    #check to make sure that the injected errors are correct
    for strand_index in injection_sites:
        sorted_nuc=sorted(injection_sites[strand_index],reverse=True)
        for nuc_index in sorted_nuc:
            fault=injection_sites[strand_index][nuc_index]
            fault_type_nuc=fault.split('-')
            fault_type=fault_type_nuc[0]
            fault_nuc=fault_type_nuc[1]

            clean_selected_strand=clean_strands[strand_index]
            error_selected_strand=final_strands[strand_index]
            #adjust the nucleotide index for erroroneos strands, need to because of insertions/deletions
            adjusted_nucleotide_index=nuc_index
            for n in sorted_nuc:
                if n < nuc_index:
                    if injection_sites[strand_index][n].split('-')[0] == '1':
                        #deletion, so subtract from true nucleotide index
                        adjusted_nucleotide_index=adjusted_nucleotide_index-1
                    elif injection_sites[strand_index][n].split('-')[0] == '2':
                        #insertion, so add to the from nucleotide index
                        adjusted_nucleotide_index=adjusted_nucleotide_index+1

            #check substitution and insertion
            if fault_type == '0' or fault_type == '2':
                assert error_selected_strand[adjusted_nucleotide_index] == fault_nuc
                #check deletion, check the removed nucleotide with the clean strand's nucleotide
            elif fault_type == '1':
                print ("clean {} fault{}".format(clean_selected_strand[nuc_index],fault_nuc))
                assert clean_selected_strand[nuc_index] == fault_nuc 


#Calculates percent difference for the generated rates and the spread sheet retes
#overall rate is used for computing the ins/del/sub given an error rate
#multiplier is used to match the input data

def calc_percent_difference(raw_value,num_strands,real_rate,multiplier, overall_rate=None):
    #analyze overall results
    if overall_rate is not None:
        measured_rate=int(float((float(raw_value)/float(overall_rate*(num_strands)))*multiplier))
    else:
        measured_rate=int(float((float(raw_value)/float((num_strands)))*multiplier))
    difference = measured_rate-real_rate
    percent_difference=abs(100*(float(difference)/float(real_rate)))
    return percent_difference

    
def allUnique(x):
    seen = set()
    return not any(i in seen or seen.add(i) for i in x)

class FiTestSuite(unittest.TestCase):
    """Fault_Injection test cases."""

    #Test strand removal
    def test_strand_removal(self):
        #test the missing strands fault injections
        arguments=fault_injector.fault_injector_arguments()
        arguments.input_file="tests/test_dna.txt"
        arguments.missing=10
        arguments.p1=20
        arguments.p2=20

        #instantiate class that will remove strands 
        missing_strands=fault_injector.miss_strand(arguments)
        missing_strands.read_file()
        clean_strands=missing_strands._input_library
        removal_locations=missing_strands.remove_sites(clean_strands)
        final_strands=missing_strands.remove_strands(removal_locations,clean_strands)

        assert len(removal_locations) ==  arguments.missing

        #make sure selections are unique
        strand_indexes=[]
        for strand_index in removal_locations:
            strand_indexes.append(strand_index)
        assert allUnique(strand_indexes) 

        
        for site in removal_locations:
            if clean_strands[site] in final_strands:
                print ("strand {} still in final strands".format(site))
                assert False

    #Test errors randomly put throughout random strands
    def test_within_dna_random_spots(self):
        arguments=fault_injector.fault_injector_arguments()
        arguments.input_file="tests/test_dna.txt"
        arguments.faulty=10
        arguments.fails=6
        arguments.p1=20
        arguments.p2=20

        model=fault_injector.strand_fault(arguments)
        model.read_file()
        clean_strands=model._input_library
        injection_sites=model.injection_sites(clean_strands)
        final_strands=model.inject_errors(injection_sites,clean_strands)

        print (injection_sites)

        
        #make sure selections for strands and nucleotides within strands are unique, and the number of faults generated is correct
        strand_indexes=[]
        for strand_index in injection_sites:
            strand_indexes.append(strand_index)
            nuc_indexes=[]
            for nuc_index in injection_sites[strand_index]:
                nuc_indexes.append(nuc_index)
            assert allUnique(nuc_indexes)
            assert len(nuc_indexes) == arguments.fails
        assert allUnique(strand_indexes)
        assert len(strand_indexes) == arguments.faulty

        check_strand_errors(injection_sites,clean_strands,final_strands)
    #Test runs of errors
    def test_within_dna_runs(self):
        arguments=fault_injector.fault_injector_arguments()
        arguments.input_file="tests/test_dna.txt"
        arguments.faulty=10
        arguments.fails=6
        arguments.p1=20
        arguments.p2=20
        arguments.run=True
        model=fault_injector.strand_fault(arguments)
        model.read_file()
        
        clean_strands=model._input_library
        injection_sites=model.injection_sites(clean_strands)
        final_strands=model.inject_errors(injection_sites,clean_strands)
        print (injection_sites)
        
        #check some properties of the injection sites  
        strand_indexes=[]
        for strand_index in injection_sites:
            strand_indexes.append(strand_index)
            nuc_indexes=[]
            for nuc_index in injection_sites[strand_index]:
                nuc_indexes.append(nuc_index)
            assert allUnique(nuc_indexes)
            #make sure that the nuc indexes are consecutive
            assert sorted(nuc_indexes) == range(min(nuc_indexes),max(nuc_indexes)+1) 
            assert len(nuc_indexes) == arguments.fails
        assert allUnique(strand_indexes)
        assert len(strand_indexes) == arguments.faulty

        #check the errors inserted into strands
        check_strand_errors(injection_sites,clean_strands,final_strands)

        
    def test_error_distribution(self):
        arguments=fault_injector.fault_injector_arguments()
        arguments.input_file="tests/test_dna.txt"
        arguments.p1=20
        arguments.p2=20
        arguments.fault_file="tests/test_rate.csv"
        model=fault_injector.strand_fault(arguments)
        model.read_file()
        csv_data=model.read_csv(arguments.fault_file)
        clean_strands=model._input_library
        temp=clean_strands[:]

        #generate a large strand pool to test the probabilistic nature of this code
        #if testing is too long, change number of strands from 100000 to something less
        while len(clean_strands)<100000:
            for strand in temp:
                clean_strands.append(strand)
            
        final_strands=model.inject_distribution(clean_strands,csv_data)

        fault_spread=model.get_fault_spread()
        del_spread = model.get_del_spread()
        ins_spread = model.get_ins_spread()
        sub_spread = model.get_sub_spread()



        real_fault_rate=model.get_fault_rate()
        real_del_rate = model.get_del_rate()
        real_ins_rate = model.get_ins_rate()
        real_sub_rate = model.get_sub_rate()
        
        #make sure the spread is only over the inner data 
        assert len(real_fault_rate) ==  (len(clean_strands[0])-arguments.p1-arguments.p2)

        print ("Error rate results")
        percent_difference_array=[]
        del_difference_array=[]
        ins_difference_array=[]
        sub_difference_array=[]

        
        for index, nuc in enumerate(sorted(fault_spread)):
            #avoid key errors if a certain nucleotide never reached an error 
            if nuc not in del_spread:
                del_spread[nuc]=0
            if nuc not in ins_spread:
                ins_spread[nuc]=0
            if nuc not in sub_spread:
                sub_spread=0


            overall_rate=(float(fault_spread[nuc])/float(len(final_strands)))
            percent_difference_array.append(calc_percent_difference(fault_spread[nuc],len(final_strands),real_fault_rate[index],10000))
            del_difference_array.append(calc_percent_difference(del_spread[nuc],len(final_strands),real_del_rate[index],1000,overall_rate))
            ins_difference_array.append(calc_percent_difference(ins_spread[nuc],len(final_strands),real_ins_rate[index],1000,overall_rate))
            sub_difference_array.append(calc_percent_difference(sub_spread[nuc],len(final_strands),real_sub_rate[index],1000,overall_rate))
            
            print ("index: {} Overall: {}  Del: {} Ins: {} Sub: {}".format(nuc,percent_difference_array[index],del_difference_array[index], ins_difference_array[index], sub_difference_array[index]))

        print ("")
        print ("")
        print ("Average percent difference Overall: {}".format(statistics.mean(percent_difference_array)))
        print ("Average percent difference Del: {}".format(statistics.mean(del_difference_array)))
        print ("Average percent difference Ins: {}".format(statistics.mean(ins_difference_array)))
        print ("Average percent difference Sub: {}".format(statistics.mean(sub_difference_array)))

        #make sure that the percent difference is less that 20 percent for each component 
        #if assertion is thrown, increase the number of strands tested to see if the percent difference decreases
        assert statistics.mean(percent_difference_array)<20 
        assert statistics.mean(del_difference_array)<20 
        assert statistics.mean(ins_difference_array)<20 
        assert statistics.mean(sub_difference_array)<20 
       
            
            

        
if __name__ == '__main__':
    unittest.main()

