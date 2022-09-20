from mpi4py import MPI
import logging
logger = logging.getLogger('dna.util.mpi_utils')
logger.addHandler(logging.NullHandler())

"""
Utils to help with sending information across processes in the DNA storage systems.
"""

def communicate_strands(strands,mpi):
    logger.info("Rank {} has {} strands at beginning of communication".format(mpi.rank,len(strands)))
    chunked_strands=[]
    chunk_size = len(strands)//mpi.size
    if chunk_size==0 and mpi.rank==0: chunk_size=len(strands) 
    leftover_strands = len(strands)%mpi.size
    if mpi.rank==0:
        for index, i in enumerate(range(0,len(strands),chunk_size)):
            end_point = i+chunk_size
            if index==mpi.size-1: end_point+=leftover_strands
            chunked_strands.append(strands[i:end_point])
            if index==mpi.size-1: break
        while len(chunked_strands)<mpi.size: chunked_strands+=[[]]
    logger.info("Rank {} at scatter".format(mpi.rank))
    return_strands=mpi.scatter(chunked_strands,root=0)
    logger.info("Rank {} has {} strands after communicate_strands".format(mpi.rank,len(return_strands)))
    return return_strands

def strand_scatter(strands,comm): #handle the scattering of a large set of strands to avoid overflow
    if comm.rank==0:logger.info("{} Total strands need to be communicated".format(len(strands)))
    strand_amount_per_transaction=100000 #100k strands has worked in the past, can make smaller if needed
    number_of_scatters = len(strands)//strand_amount_per_transaction
    if len(strands)%strand_amount_per_transaction>0: number_of_scatters+=1
    number_of_scatters=comm.bcast(number_of_scatters,root=0)
    rank_strands=[]
    for i in range(number_of_scatters):
        s=[]
        if comm.rank==0:
            s = strands[i*strand_amount_per_transaction:i*strand_amount_per_transaction+strand_amount_per_transaction]
        rank_strands+=communicate_strands(s,comm)
    logger.info("Rank {} has {} strands after complete scatter".format(comm.rank,len(rank_strands)))
    return rank_strands
        

#handle the gathering of large sets of strands to avoid overflow
def strand_gather(strands,comm):
    logger.info("Rank {} communicating {} strands back to rank 0".format(comm.Get_rank(),len(strands)))
    strand_amount_per_transaction=100000
    strands_per_rank = strand_amount_per_transaction//comm.size
    total_transactions = len(strands)//strands_per_rank
    if len(strands)%strands_per_rank>0: total_transactions+=1
    number_of_transactions=comm.allreduce(total_transactions,MPI.MAX)
    return_strands=[]
    assert number_of_transactions>0
    for i in range(number_of_transactions):
        strands_to_send = strands[i*strands_per_rank:i*strands_per_rank+strands_per_rank]
        logger.info("Rank {} sending {} strands on iteration {}".format(comm.rank,len(strands_to_send),i))
        gathered_strands=comm.gather(strands_to_send,root=0)
        logger.info("Finished strand gather iteration {}".format(i))
        if comm.rank==0:
            gathered_strands = [_ for sublist in gathered_strands for _ in sublist] 
            return_strands+=gathered_strands
    logger.info("Rank {} has {} strands after gather communication".format(comm.Get_rank(),len(return_strands)))
    return return_strands

    
