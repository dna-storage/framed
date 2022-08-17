'''
This file allows for easy customization of different codebooks to be constructed from others for use in decoders, e.g. perhaps taking certain codebook elements 
for testing purposes, using a smaller codebook, etc.
'''


from dnastorage.codec.commafreecodec import cfc_all
import os



def load_text_codebook(path):
    return_list=[]
    with open(path,'r') as codebook:
        for code in codebook:
            x = code.strip()
            if len(x)>0: return_list.append(x)
    return return_list


#-------------Miscellaneous Books------------------------------

def NONE_BOOK(): #used as a null option
    return None

#------------End Miscellaneous Books -------------------------

#----------------Codebook Dictionaries for use with CW hedges-------------------------
def CFC_ALL():#using entire cfc codebook
    codebook={}
    for i in range(0,len(cfc_all)): 
        codebook[i]=cfc_all[i]
    return codebook

def CFC_2(): #using only 2 codewords, should be 8 HD
    codebook = {
        0: cfc_all[5],
        1: cfc_all[54]
    }
    return codebook

def CFC_DUMMY():
    codebook = {
        0:"TAG",
        1:"AGCT",
        2:"AGGCT",
        3:"ACCGT"
    }
    return codebook


def ED_2_L8(): #length 8, edit distance 2
    path = "/tuck_data/kvolkel/dnastorage/dnastorage/codec/codebooks/ed_codebook_cw_size_8_ed_2.txt"
    if not os.path.exists(path):
        assert 0 and "Path to codebook not found"
    
    return {x:y for x,y in enumerate(load_text_codebook(path))}

def ED_4_L8():
    path = "/tuck_data/kvolkel/dnastorage/dnastorage/codec/codebooks/ed_codebook_cw_size_8_ed_4.txt"
    if not os.path.exists(path):
        assert 0 and "Path to codebook not found"
    return {x:y for x,y in enumerate(load_text_codebook(path))}


def ED_5_L8():
    path = "/tuck_data/kvolkel/dnastorage/dnastorage/codec/codebooks/ed_codebook_cw_size_8_ed_5.txt"
    if not os.path.exists(path):
        assert 0 and "Path to codebook not found"
    return {x:y for x,y in enumerate(load_text_codebook(path))}


def ED_5_L11():
    path = "/tuck_data/kvolkel/dnastorage/dnastorage/codec/codebooks/ed_codebook_cw_size_11_ed_5.txt"
    if not os.path.exists(path):
        assert 0 and "Path to codebook not found"
    return {x:y for x,y in enumerate(load_text_codebook(path))}
#-------------------------End Codebook Diction------------------------------------------


#------------------Synchronization Arrays to be used as Sync code books for CW Hedges------------

def TEST_SYNC_BOOK():
    syncbook=["TAATATAAC"]
    return syncbook

#----------------End Synchronization Arrays-----------------------------------------
