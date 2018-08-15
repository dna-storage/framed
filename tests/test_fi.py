
# -*- coding: utf-8 -*-

from .context import dnastorage

import unittest

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
                assert clean_selected_strand[nuc_index] == fault_nuc




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
                print "strand {} still in final strands".format(site)
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

        print injection_sites

        
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

        print injection_sites



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
        

if __name__ == '__main__':
    unittest.main()

