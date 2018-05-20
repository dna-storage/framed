from dnastorage.primer import primer_util
import nupack
from dnastorage.primer import nextera
import time

class Rule:
    def __init__(self, r, name=""):
        self.r = r
        self.name = name
        self.passed = 0
        self.total = 0
        self.total_time = 0

    def run_rule(self, strand):
        return self.r(strand)

    def run(self, strand):
        self.total += 1 
        
        t1 = time.time()
        res = self.run_rule(strand)
        t2 = time.time()
        self.total_time += (t2-t1)

        if res == True:
            self.passed += 1
            return True
        else:
            return False

    def __str__(self):
        return "{:50}{:>5} / {:<7} \t Time = {:<.2e} s".format(self.name,self.passed,self.total, self.total_time/self.total) 

class LibraryRule(Rule):
    def __init__(self, r, name="", L=[]):
        Rule.__init__(self,r,name)
        self.Library = L

    def run_rule(self, strand):
        for l in self.Library:
            if self.r(strand, l) == False:
                return False
        return True
     

def build_last_must_be_g_rule():
    return Rule(lambda s: s[-1] == 'G', "Last base must be G")

def build_repetition_rule(fraction):
    return Rule(lambda s: not (primer_util.repetitionScore(s) < fraction),"Limit repetition to {}".format(fraction))
   

def build_singlerun_rule():
    return Rule(lambda s: not primer_util.hasSingleRun(s), "Has no run")

def build_dimerrun_rule():
    return Rule(lambda s: not primer_util.hasDimerRun(s), "Has no dimer run")

def build_GC_rule(f_lo,f_hi):
    return Rule(lambda s: primer_util.checkGC(s,(f_lo,f_hi)), "GC content between {} and {}".format(f_lo,f_hi))

def build_GC_at_end_rule():
    return Rule(lambda s: primer_util.checkGC(s[-5:],(0,60)), "GC at end does not exceed 60%")

def build_Tm_rule(T_lo, T_hi):
    return Rule(lambda s: primer_util.checkTm(s,(T_lo,T_hi)), "Tm between {} and {}".format(T_lo,T_hi))

def build_check_old_strands_rule():
    return Rule(primer_util.check_old_strands, "Check compatibility with s1/s2/s3")

def build_nupack_homodimer_rule():
    return Rule(primer_util.nupack_check_homodimer, "Avoid no homodimer")

def build_nupack_nextera_bindings_rule():
    return Rule(primer_util.nupack_check_nextera_bindings, "Avoid Nextera bindings")

def build_nextera_comparison_rule():
    return Rule(lambda s: primer_util.nextera_strand_comparison(s,3), "Avoid similarity withe Nextera primers")

def build_correlation_distance_library_rule(L):
    return LibraryRule(lambda s,l: not (primer_util.correlation_distance(s,l) > 4), "Correlation distance <= 4", L)

def build_hamming_distance_library_rule(distance,L):
    return LibraryRule(lambda s,l: not (primer_util.hamming_distance(s,l) < distance), "Hamming distance > {}".format(distance), L)

def build_nupack_nonspecific_bindings_library_rule(L):
    return LibraryRule(primer_util.nupack_check_complexes, "Avoid non-specific binding with libary", L)


class DesignRules:
    def __init__(self, name=""):
        self.rules = []
        self.name = name

    def getRules(self):
        return self.rules[:]

    def addRule(self, r):
        self.rules.append(r)

    def check(self, strand):
        for r in self.rules:
            if r.run(strand) == False:
                return False        
        return True

    def __str__(self):
        return "{} Results \n\t{}".format(self.name,"\n\t".join(str(self.rules[i]) for i in range(0,len(self.rules))))

def build_standard_design_rules(Library, with_nupack=True):
    dr = DesignRules("Standard Design Rules")
    dr.addRule(build_last_must_be_g_rule())
    dr.addRule(build_singlerun_rule())
    dr.addRule(build_dimerrun_rule())
    dr.addRule(build_GC_rule(40,60))
    dr.addRule(build_repetition_rule(0.99))
    dr.addRule(build_GC_at_end_rule())
    dr.addRule(build_Tm_rule(50,60))
    dr.addRule(build_hamming_distance_library_rule(10,Library))
    dr.addRule(build_check_old_strands_rule())
    dr.addRule(build_nextera_comparison_rule())
    dr.addRule(build_correlation_distance_library_rule(Library))
    if with_nupack:
        dr.addRule(build_nupack_homodimer_rule())
        dr.addRule(build_nupack_nonspecific_bindings_library_rule(Library))
    return dr

if __name__ == "__main__":
    import string
    import random
    import nupack

    def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    s = id_generator(size=20, chars="ACGT")
    print s
    print primer_util.repetitionScore(s)

    r = build_repetition_rule(.99)
    r.run(s)
    print str(r)

    dr = build_standard_design_rules([],True)
    dr.addRule(r)
    dr.check(s)
    print dr
