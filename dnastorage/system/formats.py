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


