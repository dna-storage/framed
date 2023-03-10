#base class for defining a cluster
class BaseCluster:
    def __init__(self):
        self.mpi=None
        self.single_thread=False #default is to run multthreaded if possible (mpi)
    def Run(self,strands):
        raise NotImplemented()
    @property
    def single_thread(self):
        return self._single_thread
    @single_thread.setter
    def single_thread(self,x):
        self._single_thread=x
    @property
    def mpi(self):
        return self._mpi
    @mpi.setter
    def mpi(self,x):
        self._mpi=x
    @property
    def is_mpi_master(self):
        return self.mpi==None or self.mpi.rank==0
