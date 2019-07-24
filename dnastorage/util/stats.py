import logging

logger = logging.getLogger('dna.util.stats')
logger.addHandler(logging.NullHandler())

class dnastats:
    """ collect stats regarding the dnastorage system simulation """
    def __init__(self, msg="", fd=None):
        self.fd = fd
        self.msg = msg
        self.all_stats = {}
        self.formats = {}

    def set_fd(self, fd):
        self.fd = fd
    
    def inc(self, name, i=1, dflt=0):
        self.all_stats[name] = self.all_stats.get(name,dflt) + i
        
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
                logger.info(self.msg)
            for k,v in self.all_stats.items():
                fmt = "{},"+"{}\n".format(self.formats.get(k,"{}"))
                self.fd.write(fmt.format(k,v))
                logger.info(fmt.format(k,v))
                
    def __del__(self):
        self.persist()

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
    
