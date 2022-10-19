"""
Logging handlers/streams adapted from https://gist.github.com/sixy6e/ed35ea88ba0627e0f7dfdf115a3bf4d1.

Changes to Note made by Kevin Volkel:
1) MPI shared file pointers was converted to private file pointers, and log files are separated by rank. 
   Reason to write individual log files is that there is a concurrency bug regarding the Atomicity operation of
   shared file pointers in the implementation of NCSU HPC's MPI implementation Intel MPI 2017, as per
   the documented issues here - https://groups.google.com/g/mpi4py/c/SaNzc8bdj6U.
"""


import logging
from mpi4py import MPI
import structlog
import os

class MPIIOStream(object):

    """
    A very basic MPI stream handler for synchronised I/O.
    """

    def __init__(self, filename, comm, mode):
        self._file = MPI.File.Open(MPI.COMM_SELF,filename, mode)
        print(self._file.Get_info())
    def write(self, msg):
        # if for some reason we don't have a unicode string...
        try:
            msg = msg.encode()
        except AttributeError:
            pass
        self._file.Write_shared(msg)

    def sync(self):
        """
        Synchronise the processes
        """
        self._file.Sync()

    def close(self):
        self.sync()
        self._file.Close()


class MPIFileHandler(logging.StreamHandler):

    """
    A basic MPI file handler for writing log files.
    Internally opens a synchronised MPI I/O stream via MPIIOStream.
    Ideas and some code from:
    * https://groups.google.com/forum/#!topic/mpi4py/SaNzc8bdj6U
    * https://gist.github.com/JohnCEarls/8172807
    * https://stackoverflow.com/questions/45680050/cannot-write-to-shared-mpi-file-with-mpi4py
    """

    def __init__(self, filename,
                 mode=MPI.MODE_WRONLY|MPI.MODE_CREATE, comm=MPI.COMM_WORLD):
        self.filename = filename
        if os.path.exists(self.filename) and os.path.isfile(self.filename):
            os.remove(self.filename)
        self.mode = mode
        self.comm = comm

        super(MPIFileHandler, self).__init__(self._open())

    def _open(self):
        stream = MPIIOStream(self.filename, self.comm, self.mode)
        return stream

    def close(self):
        if self.stream:
            self.stream.close()
            self.stream = None


    def flush(self):
        self.stream.sync()

    def emit(self, record):
        """
        Emit a record.
        We have to override emit, as the logging.StreamHandler has 2 calls
        to 'write'. The first for the message, and the second for the
        terminator. This posed a problem for mpi, where a second process
        could call 'write' in between these two calls and create a
        conjoined log message.
        """
        msg = self.format(record)
        self.stream.write('{}{}'.format(msg, self.terminator))
        self.flush()


def main():
    """
    A sample test run.
    """
    logger = logging.getLogger("node[%i]"%comm.rank)
    logger.setLevel(logging.DEBUG)
    mpi_handler.setFormatter(formatter)
    logger.addHandler(mpi_handler)

    # sample log levels
    logger.debug('debug message')
    logger.info('info message')
    logger.error('error message')
    logger.critical('critical message')


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    mpi_handler = MPIFileHandler("test{}.log".format(comm.rank))
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    main()
