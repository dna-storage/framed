import logging
from random import choice
import string
import pickle
import numpy as np
import os

logger = logging.getLogger('dna.util.stats')
logger.addHandler(logging.NullHandler())


def random_string(l=5):
    return ''.join([choice(string.ascii_letters + string.digits) for n in range(l)])

class dnastats(object):
    """ collect stats regarding the dnastorage system simulation """
    def __init__(self, msg="", fd=None,pickle_fd=None):
        self.fd = fd
        self.pickle_fd=pickle_fd
        self.msg = msg
        self.all_stats = {} #dictionary for holding all stats
        self.formats = {} #formatting for statistics
        self.file_register={} #registration dictionary for separating outputs into different files
        self.hist_register={} #registration dictionary for condensing data into histograms, we use numpy histograms
        self.experiment_counter=0 #counter that is used to separate results 
        self._name_to_counter={}
    def set_fd(self, fd):
        self.fd = fd
    def set_pickle_fd(self,fd):
        self.pickle_fd=fd
    def register_file(self,stats_name,filename):
        self.file_register[stats_name]=filename #note, this just registerse a statistic with a file, should check stat exists before persisting it
    def register_hist(self,stats_name):
        self.hist_register[stats_name]=None
    def get_next_name(self,name):
        if name not in self._name_to_counter:
            self._name_to_counter[name]=self._experiment_counter+1
        else:
            self._name_to_counter[name]+=1
        return "{}::{}".format(name,self._name_to_counter[name]-1)

    # Use inc to track and dump a counter. example:
    #    stats.inc("block::errors")
    def inc(self, name, i=1, dflt=0,coords=None):
        if coords is None:self.all_stats[name] = self.all_stats.get(name,dflt) + i
        else:
            if name not in self.all_stats:
                self.all_stats[name]=dflt
            if type(self.all_stats[name]) is dict:
                self.all_stats[name][coords] = self.all_stats[name].get(coords,0)+i
            else:
                self.all_stats[name][coords] = self.all_stats[name][coords]+i
    #use clear to reset stat counters
    def clear(self):
        self.all_stats={}
        self.formats={}

    @property
    def experiment_counter(self):
        return self._experiment_counter
    @experiment_counter.setter
    def experiment_counter(self,x):
        self._name_to_counter={}
        self._experiment_counter=x
    
    # Use append to record a list of data
    def append(self, name, s, dflt=[]):
        self.all_stats[name] = self.all_stats.get(name,dflt) + [s]

    def unique(self, name, val):
        while name in self.all_stats:
            name += "::"+random_string(5)
        self.all_stats[name] = val            
        
    def __getitem__(self, name):
        return self.all_stats[name]

    def __setitem__(self, name, value):
        self.all_stats[name] = value

    def __contains__(self,name):
        return name in self.all_stats
        
    def format(self, name, f):
        self.formats[name] = f

    def persist(self):
        main_stats = {}
        for key in self.all_stats:
            if key in self.hist_register: self.all_stats[key]=np.histogram(self.all_stats[key],bins="auto")
            if key not in self.file_register: main_stats[key]=self.all_stats[key]
        if not (self.fd is None):
            if len(self.msg) > 0:
                self.fd.write(self.msg+"\n")
            items = main_stats.items()
            items=sorted(items,key=lambda x: x[0])
            for k,v in items:
                fmt = "{}="+"{}\n".format(self.formats.get(k,"{}"))
                self.fd.write(fmt.format(k,v))
        else:
            items = main_stats.items()
            items=sorted(items,key=lambda x: x[0])
            for k,v in items:
                fmt = "{}="+"{}".format(self.formats.get(k,"{}"))
        if not (self.pickle_fd is None): #dumping pickle of stats 
            pickle.dump(main_stats,self.pickle_fd)
        self._write_file_register(os.path.dirname(self.fd.name))

    def _write_file_register(self,output_dir):
        """
        Write out the file register to separate files. Stats that dump to the same external file
        outside of the regular stats file will be collected together and then dumped
        """
        file_dictionaries={}
        for stat_key in self.file_register:
            if stat_key not in self.all_stats: continue #stat was never filled in
            file_key = self.file_register[stat_key]
            if file_key not in file_dictionaries:
                file_dictionaries[file_key]={}
            file_dictionaries[file_key][stat_key]=self.all_stats[stat_key]
        for f in file_dictionaries:
            dump_dict = file_dictionaries[f]
            fd = open(os.path.join(output_dir,f),"w+")
            pickle_fd = open(os.path.join(output_dir,"{}.pickle".format(f)),"wb+")
            pickle.dump(dump_dict,pickle_fd)
            items = dump_dict.items()
            items=sorted(items,key=lambda x: x[0])
            for k,v in items:
                fmt = "{}="+"{}\n".format(self.formats.get(k,"{}"))
                fd.write(fmt.format(k,v))
            pickle_fd.close()
            
    def aggregate(self,other_stats,copy_list=[]): #aggregate other_stats into this stats
        assert isinstance(other_stats,dnastats)
        self.file_register={**self.file_register,**other_stats.file_register}
        self.hist_register={**self.hist_register,**other_stats.hist_register}
        """
        copy_list is understood as a list of keys that should just be copied between two stat.
        In order to allow the most flexibility these keys are understood as sub-strings. That is,
        if anything in copy_list matches a substring of a key in other.all_stats then that is called a match
        in which values are only copied and nothing is accumulated.
        """
        skip_list=[]
        for c in copy_list:
            for d in other_stats.all_stats:
                if c in d: skip_list.append(d)
        for x in other_stats.all_stats:
            if x not in self.all_stats:
                self.all_stats[x] = other_stats.all_stats[x]
            elif x not in skip_list:
                if type(other_stats.all_stats[x]) is dict:
                    self.all_stats[x]=self._aggregate_dict(self.all_stats[x],other_stats.all_stats[x])
                else:
                    self.all_stats[x]=self.all_stats[x]+other_stats.all_stats[x]

    def _aggregate_dict(self,d1,d2):
        merged_dict = d1
        for key,value in d2.items():
            if key in d1:
                merged_dict[key]+=d2[key]
            else:
                merged_dict[key]=d2[key]
        return merged_dict
    
            
    def __del__(self):
        return
        #self.persist()
    

global stats
stats = dnastats(msg="global stats collection for dnastorage project")

if __name__ == "__main__":
    import sys
    from dnastorage.util.stats import dnastats,stats
    logger.setLevel(logging.DEBUG)
    _formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    _ch = logging.FileHandler("dna.debug.log",mode='w')
    _ch.setFormatter(_formatter)
    logger.addHandler(_ch)
#    stats.set_fd(sys.stdout)
    stats.inc("unique_strands")
    stats.inc("unique_strands")
    stats.inc("unique_strands")
    
