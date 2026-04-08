#import exceptions
from dnastorage.exceptions import *

#import all pipelines that have been defined 
from dnastorage.arch.builder import *

# DO NOT ALTER ENTRIES IN THIS TABLE, BUT YOU MAY ADD NEW ONES
# ALL CHANGES NEED TO BE THOROUGHLY TESTED FOR BACKWARDS COMPATIBILITY

FileSystemFormats = {
    #------ Pipelines
    0x0703 : [0x0703,208,15,"BasicHedges","Basic Hedges implementation with just primers",
              Basic_Hedges_Pipeline,Basic_Hedges_Pipeline],
    0x0704 : [0x0704,0,0,"ReedSolomon_Base4_Pipeline","Reed solomon with base encoding and clustering algorithms",
              ReedSolomon_Base4_Pipeline,ReedSolomon_Base4_Pipeline],
    0x0705 : [0x0705,0,0,"Fountain_Base4_Pipeline","LT fountain code with base4 encoding and clustering algorithms",
              Fountain_Base4_Pipeline,Fountain_Base4_Pipeline],
    0x0706 : [0x0706,208,15,"Fountain_Hedges_Pipeline","LT fountain code with Hedges inner code",
              Fountain_Hedges_Pipeline,Fountain_Hedges_Pipeline],
    0x0707 : [0x0707,0,0,"ReedSolomon_Base4_FileLevelFountain_Pipeline",
              "Reed-Solomon + Base4 with file-level LT fountain erasure protection",
              ReedSolomon_Base4_FileLevelFountain_Pipeline,ReedSolomon_Base4_FileLevelFountain_Pipeline],
    0x0708 : [0x0708,0,0,"Fountain_Base4_FileLevelFountain_Pipeline",
              "LT fountain + Base4 with file-level LT fountain erasure protection",
              Fountain_Base4_FileLevelFountain_Pipeline,Fountain_Base4_FileLevelFountain_Pipeline],

}


def file_system_formats():
    return [ v[3] for k,v in FileSystemFormats.items() ]

_abbrevFileSystemDict = { v[3] : v for k,v in FileSystemFormats.items() }

def file_system_format_description(formatid):
    return FileSystemFormats[formatid][4]

def file_system_format_packetsize(formatid):
    return FileSystemFormats[formatid][2]

def file_system_encoder(formatid):
    return FileSystemFormats[formatid][5]

def file_system_decoder(formatid):
    return FileSystemFormats[formatid][6]

def file_system_encoder_by_abbrev(ab):
    return _abbrevFileSystemDict[ab][5]

def file_system_decoder_by_abbrev(ab):
    return _abbrevFileSystemDict[ab][6]

def file_system_formatid_by_abbrev(ab):
    return _abbrevFileSystemDict[ab][0]

def file_system_abbrev(id):
    return FileSystemFormats[id][3]


