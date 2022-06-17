import logging
from random import choice
import string
import pickle
import numpy as np


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
        self.all_stats = {}
        self.formats = {}

    def set_fd(self, fd):
        self.fd = fd
    def set_pickle_fd(self,fd):
        self.pickle_fd=fd

    # Use inc to track and dump a counter. example:
    #    stats.inc("block::errors")
    # If called 5 times, this will produce a log entry like this:
    #    block::errors=5
    #
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

    def format(self, name, f):
        self.formats[name] = f

    def persist(self):
        if not (self.fd is None):
            if len(self.msg) > 0:
                self.fd.write(self.msg+"\n")
                logger.info(self.msg.strip("\n"))

            items = self.all_stats.items()
            items=sorted(items,key=lambda x: x[0])
                
            for k,v in items:
                fmt = "{}="+"{}\n".format(self.formats.get(k,"{}"))
                self.fd.write(fmt.format(k,v))
                logger.info(fmt.format(k,v))
        else:
            if len(self.msg) > 0:
                logger.info(self.msg)
            items = self.all_stats.items()
            items=sorted(items,key=lambda x: x[0])
            for k,v in items:
                fmt = "{}="+"{}".format(self.formats.get(k,"{}"))
                logger.info(fmt.format(k,v))
                
        if not (self.pickle_fd is None): #dumping pickle of stats 
            pickle.dump(self.all_stats,self.pickle_fd)


    def aggregate(self,other_stats,copy_list=[]):
        assert isinstance(other_stats,dnastats)
        #aggregate other_stats into this stats
        for x in other_stats.all_stats:
            if x not in self.all_stats:
                self.all_stats[x] = other_stats.all_stats[x]
            elif x not in copy_list and not type(other_stats.all_stats[x]) is dict:
                self.all_stats[x]=self.all_stats[x]+other_stats.all_stats[x]    
            
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
    
