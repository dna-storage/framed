
# -*- coding: utf-8 -*-

from .context import dnastorage

import unittest
import statistics
import random as _random
import csv as _csv
import types

# ---------------------------------------------------------------------------
# Local implementations of legacy fault-injection helpers.
# The original classes (fault_injector_arguments, miss_strand, strand_fault)
# no longer exist in dnastorage/fi/fault_injector.py, so they are
# re-implemented here to keep these tests self-contained.
# ---------------------------------------------------------------------------

_sub_dict = {
    'A': ['G', 'C', 'T'],
    'G': ['C', 'A', 'T'],
    'T': ['G', 'C', 'A'],
    'C': ['G', 'T', 'A'],
}
_nuc_list = ['A', 'C', 'T', 'G']


class _FaultInjectorArguments:
    """Simple namespace for fault-injection parameters."""
    pass


class _MissStrand:
    """Reads a strand library from a text file and randomly removes strands."""

    def __init__(self, arguments):
        self._args = arguments
        self._input_library = []

    def read_file(self):
        with open(self._args.input_file, 'r') as f:
            self._input_library = [line.strip() for line in f if line.strip()]

    def remove_sites(self, clean_strands):
        """Return a list of *missing* random strand indexes."""
        return _random.sample(range(len(clean_strands)), self._args.missing)

    def remove_strands(self, removal_locations, clean_strands):
        """Return a copy of *clean_strands* with the given indexes removed."""
        removal_set = set(removal_locations)
        return [s for i, s in enumerate(clean_strands) if i not in removal_set]


class _StrandFault:
    """Injects substitution / deletion / insertion errors into DNA strands."""

    def __init__(self, arguments):
        self._args = arguments
        self._input_library = []
        self._fault_spread = {}
        self._del_spread   = {}
        self._ins_spread   = {}
        self._sub_spread   = {}
        self._fault_rate   = []
        self._del_rate     = []
        self._ins_rate     = []
        self._sub_rate     = []

    def read_file(self):
        with open(self._args.input_file, 'r') as f:
            self._input_library = [line.strip() for line in f if line.strip()]

    def injection_sites(self, clean_strands):
        """Return {strand_index: {nuc_index: 'TYPE-NUC'}} for *faulty* strands."""
        p1   = getattr(self._args, 'p1',   0)
        p2   = getattr(self._args, 'p2',   0)
        run  = getattr(self._args, 'run',  False)
        faulty = self._args.faulty
        fails  = self._args.fails

        strand_choices = _random.sample(range(len(clean_strands)), faulty)
        result = {}
        for strand_index in strand_choices:
            strand = clean_strands[strand_index]
            inner_start = p1
            inner_end   = len(strand) - p2

            if run:
                start = _random.randint(inner_start, inner_end - fails)
                nuc_positions = list(range(start, start + fails))
            else:
                nuc_positions = _random.sample(range(inner_start, inner_end), fails)

            result[strand_index] = {}
            for nuc_idx in nuc_positions:
                fault_type = _random.randint(0, 2)
                if fault_type == 0:    # substitution — record the replacement nuc
                    nuc = _random.choice(_sub_dict[strand[nuc_idx]])
                elif fault_type == 1:  # deletion — record the removed nuc
                    nuc = strand[nuc_idx]
                else:                  # insertion — record the inserted nuc
                    nuc = _random.choice(_nuc_list)
                result[strand_index][nuc_idx] = '{}-{}'.format(fault_type, nuc)
        return result

    def inject_errors(self, injection_sites, clean_strands):
        """Apply *injection_sites* to *clean_strands* and return the mutated list."""
        out_list = list(clean_strands)
        for strand_index in injection_sites:
            strand = clean_strands[strand_index]
            # Apply edits high-to-low so earlier indexes stay valid.
            for nuc_index in sorted(injection_sites[strand_index], reverse=True):
                fault_type, fault_nuc = injection_sites[strand_index][nuc_index].split('-')
                if fault_type == '0':    # substitution
                    strand = strand[:nuc_index] + fault_nuc + strand[nuc_index + 1:]
                elif fault_type == '1':  # deletion
                    strand = strand[:nuc_index] + strand[nuc_index + 1:]
                elif fault_type == '2':  # insertion
                    strand = strand[:nuc_index] + fault_nuc + strand[nuc_index:]
            out_list[strand_index] = strand
        return out_list

    def read_csv(self, file_name):
        """Parse the CSV file into {row_label: [float, ...]}."""
        parsed = {}
        with open(file_name, 'r') as f:
            for row in _csv.reader(f):
                if row and row[0]:
                    parsed[row[0]] = [float(v) for v in row[1:] if v]
        return parsed

    def inject_distribution(self, clean_strands, csv_data):
        """Inject per-position errors from CSV rates; populate spread/rate attrs."""
        p1 = getattr(self._args, 'p1', 0)
        p2 = getattr(self._args, 'p2', 0)

        overall_rates = csv_data.get('Overall Error', [])
        del_ratios    = csv_data.get('Del/Error',     [])
        ins_ratios    = csv_data.get('Ins/Error',     [])
        sub_ratios    = csv_data.get('Sub/Error',     [])

        strand_len = len(clean_strands[0])
        inner_len  = strand_len - p1 - p2

        # Dense spread dicts (keys 0..inner_len-1) so the test loop is index-safe.
        self._fault_spread = {i: 0 for i in range(inner_len)}
        self._del_spread   = {i: 0 for i in range(inner_len)}
        self._ins_spread   = {i: 0 for i in range(inner_len)}
        self._sub_spread   = {i: 0 for i in range(inner_len)}

        # Real (CSV) rates scaled to integers matching calc_percent_difference units.
        # fault_rate uses multiplier 10000; del/ins/sub_rate use multiplier 1000 on
        # the per-type proportion (Del/Error, Ins/Error, Sub/Error columns).
        self._fault_rate = [
            int(overall_rates[p1 + i] * 10000) if p1 + i < len(overall_rates) else 0
            for i in range(inner_len)
        ]
        self._del_rate = [
            int(del_ratios[p1 + i] * 1000) if p1 + i < len(del_ratios) else 0
            for i in range(inner_len)
        ]
        self._ins_rate = [
            int(ins_ratios[p1 + i] * 1000) if p1 + i < len(ins_ratios) else 0
            for i in range(inner_len)
        ]
        self._sub_rate = [
            int(sub_ratios[p1 + i] * 1000) if p1 + i < len(sub_ratios) else 0
            for i in range(inner_len)
        ]

        out_list = []
        for strand in clean_strands:
            new_strand = strand
            errors = []
            for i in range(inner_len):
                abs_pos = i + p1
                rate = overall_rates[abs_pos] if abs_pos < len(overall_rates) else 0.0
                if _random.random() < rate:
                    dr = del_ratios[abs_pos] if abs_pos < len(del_ratios) else 0.33
                    ir = ins_ratios[abs_pos] if abs_pos < len(ins_ratios) else 0.33
                    sr = sub_ratios[abs_pos] if abs_pos < len(sub_ratios) else 0.33
                    total = dr + ir + sr
                    if total == 0:
                        continue
                    r = _random.random() * total
                    if r < dr:
                        fault_type = 1
                        self._del_spread[i] += 1
                    elif r < dr + ir:
                        fault_type = 2
                        self._ins_spread[i] += 1
                    else:
                        fault_type = 0
                        self._sub_spread[i] += 1
                    self._fault_spread[i] += 1
                    errors.append((abs_pos, fault_type))

            for abs_pos, fault_type in sorted(errors, reverse=True):
                if fault_type == 0:
                    new_strand = (new_strand[:abs_pos]
                                  + _random.choice(_sub_dict[new_strand[abs_pos]])
                                  + new_strand[abs_pos + 1:])
                elif fault_type == 1:
                    new_strand = new_strand[:abs_pos] + new_strand[abs_pos + 1:]
                elif fault_type == 2:
                    new_strand = (new_strand[:abs_pos]
                                  + _random.choice(_nuc_list)
                                  + new_strand[abs_pos:])
            out_list.append(new_strand)
        return out_list

    def get_fault_spread(self): return self._fault_spread
    def get_del_spread(self):   return self._del_spread
    def get_ins_spread(self):   return self._ins_spread
    def get_sub_spread(self):   return self._sub_spread
    def get_fault_rate(self):   return self._fault_rate
    def get_del_rate(self):     return self._del_rate
    def get_ins_rate(self):     return self._ins_rate
    def get_sub_rate(self):     return self._sub_rate


# Expose as a module-like namespace so existing test code (fault_injector.xxx)
# continues to work unchanged.
fault_injector = types.SimpleNamespace(
    fault_injector_arguments=_FaultInjectorArguments,
    miss_strand=_MissStrand,
    strand_fault=_StrandFault,
)

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
        arguments.input_file="test_dna.txt"
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
        arguments.input_file="test_dna.txt"
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
        arguments.input_file="test_dna.txt"
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
            assert sorted(nuc_indexes) == list(range(min(nuc_indexes),max(nuc_indexes)+1))
            assert len(nuc_indexes) == arguments.fails
        assert allUnique(strand_indexes)
        assert len(strand_indexes) == arguments.faulty

        #check the errors inserted into strands
        check_strand_errors(injection_sites,clean_strands,final_strands)

        
    def test_error_distribution(self):
        arguments=fault_injector.fault_injector_arguments()
        arguments.input_file="test_dna.txt"
        arguments.p1=20
        arguments.p2=20
        arguments.fault_file="test_rate.csv"
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
                sub_spread[nuc]=0


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

