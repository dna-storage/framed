#base class for defining a cluster
class BaseAlignment:
    def __init__(self):
        self.mpi=None
    def Run(self,cluster):
        raise NotImplemented()
    @property
    def mpi(self):
        return self._mpi
    @mpi.setter
    def mpi(self,x):
        self._mpi=x
