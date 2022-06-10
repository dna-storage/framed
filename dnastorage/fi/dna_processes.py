'''
Filename: dna_processes.py
Description: functions to carry out some basic functions to transform DNA strands
'''



def T7_Transcription(dna_strand):
    t7 = "TAATACGACTCACTATAG"
    if not t7 in dna_strand:
        return dna_strand
    index = dna_strand.find(t7)
    assert index != -1
    return dna_strand[index+len(t7):]

def dna_process_nop(dna_strand):
    return dna_strand

    
