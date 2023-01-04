import mpi4py
mpi4py.rc.recv_mprobe = False
from mpi4py import MPI
import logging
import pickle
logger = logging.getLogger('dna.util.mpi_utils')
logger.addHandler(logging.NullHandler())

"""
Utils to help with sending information across processes in the DNA storage systems.
"""

def communicate_objects(objects,mpi):
    logger.info("Rank {} has {} objects at beginning of communication".format(mpi.rank,len(objects)))
    chunked_objects=[]
    chunk_size = len(objects)//mpi.size
    if chunk_size==0 and mpi.rank==0: chunk_size=len(objects) 
    leftover_objects = len(objects)%mpi.size
    if mpi.rank==0:
        for index, i in enumerate(range(0,len(objects),chunk_size)):
            end_point = i+chunk_size
            if index==mpi.size-1: end_point+=leftover_objects
            chunked_objects.append(objects[i:end_point])
            if index==mpi.size-1: break
        while len(chunked_objects)<mpi.size: chunked_objects+=[[]]
    logger.info("Rank {} at scatter".format(mpi.rank))
    
    return_objects=mpi.scatter(chunked_objects,root=0)
    logger.info("Rank {} has {} objects after communicate_objects".format(mpi.rank,len(return_objects)))
    return return_objects

def object_scatter(objects,comm,objs_per_transaction=10000): #handle the scattering of a large set of objects to avoid overflow
    if comm.rank==0:
        logger.info("{} Total objects need to be communicated".format(len(objects)))
        if len(objects)>0: logger.info("Example pickle size of object is {}".format(len(pickle.dumps(objects[0]))))
    number_of_scatters = len(objects)//objs_per_transaction
    if len(objects)%objs_per_transaction>0: number_of_scatters+=1
    number_of_scatters=comm.bcast(number_of_scatters,root=0)
    rank_objects=[]
    for i in range(number_of_scatters):
        s=[]
        if comm.rank==0:
            s = objects[i*objs_per_transaction:i*objs_per_transaction+objs_per_transaction]
        rank_objects+=communicate_objects(s,comm)
    logger.info("Rank {} has {} objects after complete scatter".format(comm.rank,len(rank_objects)))
    return rank_objects
        

#handle the gathering of large sets of objects to avoid overflow
def object_gather(objects,comm,objs_per_transaction=10000):
    logger.info("Rank {} communicating {} objects back to rank 0".format(comm.Get_rank(),len(objects)))
    if len(objects)>0: logger.info("Example pickle size of object is {}".format(len(pickle.dumps(objects[0]))))
    objects_per_rank = objs_per_transaction//comm.size
    total_transactions = len(objects)//objects_per_rank
    if len(objects)%objects_per_rank>0: total_transactions+=1
    logger.info("Rank {} total_transactions {} ".format(comm.rank,total_transactions))        
    number_of_transactions=comm.allreduce(total_transactions,MPI.MAX)
    return_objects=[]
    if  number_of_transactions==0: return [] #nothing to actually send 
    for i in range(number_of_transactions):
        objects_to_send = objects[i*objects_per_rank:i*objects_per_rank+objects_per_rank]
        logger.info("Rank {} sending {} objects on iteration {}".format(comm.rank,len(objects_to_send),i))
        gathered_objects=comm.gather(objects_to_send,root=0)
        logger.info("Finished object gather iteration {}".format(i))
        if comm.rank==0:
            gathered_objects = [_ for sublist in gathered_objects for _ in sublist] 
            return_objects+=gathered_objects
        logger.info("Rank {} has {} objects after gather communication".format(comm.Get_rank(),len(return_objects)))
    return return_objects

    
