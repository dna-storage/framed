'''
This file allows for easy customization of different codebooks to be constructed from others for use in decoders, e.g. perhaps taking certain codebook elements 
for testing purposes, using a smaller codebook, etc.
'''


from dnastorage.codec.commafreecodec import cfc_all


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
        0:"CCCC",
        1:"GGGG",
        2:"TTTT",
        3:"AAAA"
    }
    return codebook


def TEST_SYNC_BOOK():
    syncbook=["AAAAAA",
              "GGGGGG"
    ]
    return syncbook
