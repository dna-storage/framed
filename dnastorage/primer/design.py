from dnastorage.primer import primer_util
import nupack
from dnastorage.primer import nextera
import time
import random

class RandomPrimerGenerator:
    def __init__(self, chars="AGCT", length=20):
        self.chars = chars
        self.len = length

    def get(self):
        return ''.join(random.choice(self.chars) for _ in range(self.len))

class UniquePrimerGenerator(RandomPrimerGenerator):
    def __init__(self, chars="AGCT", length=20, library=[]):
        RandomPrimerGenerator.__init__(self,chars,length)
        self.avoid = library[:]

    def append(self, l):
        self.avoid.append(l)

    def _get_helper(self):
        return RandomPrimerGenerator.get(self)

    def get(self):
        while True:
            s = self._get_helper()
            if not (s in self.avoid):
                self.avoid.append(s)
                return s

class SpecificPrimerGenerator:
    def __init__(self,value, chars="AGCT", length=20):
        self.chars = chars
        self.len = length

    def get(self,v):
        l=''
        value=v
        for i in range(19,-1,-1):
            if(value>=3*4**i):
              l += "T"
              value-=3*4**i
            elif(value>=2*4**i):
              l += "C"
              value-=2*4**i
            elif(value>=1*4**i):
              l += "G"
              value-=1*4**i
            else:
              l += "A"
        #print value
        # print l
        return l

class LinearPrimerGenerator(SpecificPrimerGenerator):
    def __init__(self, chars="AGCT", length=20, library=[]):
        SpecificPrimerGenerator.__init__(self,chars,length)
        self.avoid = library[:]

    def append(self, l):
        self.avoid.append(l)

    def _get_helper(self,i):
        return SpecificPrimerGenerator.get(self,i)

    def get(self,i):
        while True:
            s = self._get_helper(i)
            if not (s in self.avoid):
                self.avoid.append(s)
                return s

class LikelyPrimerGenerator(UniquePrimerGenerator):
    def __init__(self, chars="AGCT", length=20, library=[], last='G', repl_factor=5):
        UniquePrimerGenerator.__init__(self,chars,length,library)
        self.last = last
        self.repl_factor = repl_factor

    def _get_helper(self):
        l = [random.choice(self.chars) for _ in range(self.len)]
        # force replication
        for j in range(self.len/self.repl_factor+1):
            i = random.randint(-2,3)
            if i+j*self.repl_factor+1 < self.len and i+j*self.repl_factor > 0:
                l[i+j*5+1] = l[i+j*5]

        l[-1] = self.last
        return "".join(l)


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
        return "{:50}{:>5} / {:<7} \t Time = {:<.2e} s".format(self.name,self.passed,self.total, self.total_time/(max(self.total,1))) 

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

def build_norepeat_rule():
    return Rule(lambda s: not primer_util.hasRepeat(s), "Has no repeat")

def build_longrun_rule():
    return Rule(lambda s: not primer_util.hasLongRun(s), "Has no run of 4 or longer")

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
    return LibraryRule(lambda s,l: not (primer_util.hamming_distance(s,l) <= distance), "Hamming distance > {}".format(distance), L)

def build_nupack_nonspecific_bindings_library_rule(L, Tm=50):
    return LibraryRule(lambda s,l: primer_util.nupack_check_complexes(s,l,Tm), "Avoid non-specific binding with libary (Tm={})".format(Tm), L)

def build_nupack_heterodimer_bindings_library_rule(L, Tm=50):
    return LibraryRule(lambda s,l: primer_util.nupack_check_heterodimer_complexes(s,l,Tm), "Avoid dimer-dimer binding with libary (Tm={})".format(Tm), L)

def build_nupack_allow_nonspecific_bindings_library_rule(L, Tm=50):
    return LibraryRule(lambda s,l: not primer_util.nupack_check_complexes(s, l, Tm), "Allow non-specific binding with library (Tm={})".format(Tm), L)

def build_nupack_allow_heterodimer_bindings_library_rule(L, Tm=50):
    return LibraryRule(lambda s,l: not primer_util.nupack_check_heterodimer_complexes(s, l, Tm), "Allow dimer-dimer binding with library (Tm={})".format(Tm), L)

class DesignRules:
    def __init__(self, name=""):
        self.rules = []
        self.name = name

    def get_rules(self):
        return self.rules[:]

    def add_rule(self, r):
        self.rules.append(r)

    def check(self, strand):
        for r in self.rules:
            if r.run(strand) == False:
                return False        
        return True

    def __str__(self):
        return "{} Results \n\t{}".format(self.name,"\n\t".join(str(self.rules[i]) for i in range(0,len(self.rules))))

def build_standard_design_rules(Library, with_nupack=True, distance=10):

    dr = DesignRules("Standard Design Rules")
    dr.add_rule(build_last_must_be_g_rule())
    dr.add_rule(build_singlerun_rule())
    dr.add_rule(build_dimerrun_rule())
    dr.add_rule(build_GC_rule(40,60))
    dr.add_rule(build_repetition_rule(0.99))
    dr.add_rule(build_GC_at_end_rule())
    dr.add_rule(build_Tm_rule(50,60))
    dr.add_rule(build_hamming_distance_library_rule(distance,Library))
    dr.add_rule(build_check_old_strands_rule())
    dr.add_rule(build_nextera_comparison_rule())
    dr.add_rule(build_correlation_distance_library_rule(Library))
    if with_nupack:
        dr.add_rule(build_nupack_homodimer_rule())
        dr.add_rule(build_nupack_nonspecific_bindings_library_rule(Library))
        dr.add_rule(build_nupack_heterodimer_bindings_library_rule(Library))
    return dr



###
### Deprecated function -- DO NOT USE 
### Will be deleted soon 
###
def design_rules_met(s,L,nupack,nextera_binding):
    #print s
    if primer_util.repetitionScore(s) < 0.99:
        #print "reptition score", repetitionScore(s)
        return False
    if primer_util.hasSingleRun(s) or primer_util.hasDimerRun(s): # or hasDimer(seq,5):
        #print "runs"
        return False

    # ending of primer should not have too much GC content
    if not primer_util.checkGC(s[-5:],(0,60)):
        #print s[-5:],"GC short content",GC(s[-5:])
        return False

    for l in L:     
        if primer_util.correlation_distance(s,l)>4:
            #print "correlation"
            return False
        if primer_util.hamming_distance(s,l) < 10:
            #print "hamming"
            return False    
        if nupack:
            # repeat s to search for homo-dimers
            c = primer_util.checkComplexes([s,l])
            if len(c) > 0:
                print ("*****Heterodimer: {} vs {} -- {}".format(l,s, c[0]['pattern']))
                return False

            c = primer_util.checkComplexes([l,primer_util.reverse_complement(s)])
            if len(c) > 0:
                print ("*****False binding: {} vs {} -- {}".format(l,s, c[0]['pattern']))
                return False

            c = primer_util.checkComplexes([s,primer_util.reverse_complement(l)])
            if len(c) > 0:
                print ("*****False binding: {} vs {} -- {}".format(l,s, c[0]['pattern']))
                return False

        #if reverse_correlation_distance(s,l)>3:
        #    show_correlation(s,reverse(l))
        #    return False

    if not primer_util.checkTm(s,(50.0,60.0)):
        #print "TM"
        return False

    if not primer_util.checkGC(s,(40.0,60.0)):
        #print "GC full"
        return False

    if not primer_util.nextera_strand_comparison(s,3):
        #print "nextera"
        return False

    if nupack and nextera_binding:
        nextera_primers = nextera.get_nextera_primers()
        for l in nextera_primers:
            # repeat s to search for homo-dimers
            c = primer_util.checkComplexes([s,l])
            if len(c) > 0:
                print ("*****Nextera Heterodimer: {} vs {} -- {}".format(l,s, c[0]['pattern']))
                return False

            c = primer_util.checkComplexes([l,primer_util.reverse_complement(s)])
            if len(c) > 0:
                print ("*****Nextera False binding: {} vs {} -- {}".format(l,s, c[0]['pattern']))
                return False

            c = primer_util.checkComplexes([s,primer_util.reverse_complement(l)])
            if len(c) > 0:
                print ("*****Nextera Complement False binding: {} vs {} -- {}".format(l,s, c[0]['pattern']))
                return False


    if not primer_util.check_old_strands(s):
        #print "check old"
        return False

    if nupack:
        # repeat s to search for homo-dimers
        c = primer_util.checkComplexes([s,s])
        if len(c) > 0:
            print ("*****Homodimer: {} {}".format(s, c[0]['pattern']))
            return False

    return True




if __name__ == "__main__":
    import string
    import random
    import nupack

    def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    s = id_generator(size=20, chars="ACGT")
    print (s)
    print (primer_util.repetitionScore(s))

    r = build_repetition_rule(.99)
    r.run(s)
    print (str(r))

    dr = build_standard_design_rules([],True)
    dr.add_rule(r)
    dr.check(s)
    print (dr)
