from dnastorage.codec import dense
from dnastorage.codec import norepeatscodec
from dnastorage.codec import commafreecodec
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec.rscodec import *
from dnastorage.arch.strand import *

def ENC_FSMD_200(pf, primer5, primer3, bIndex=0):
    b = commafreecodec.CommaFreeCodec(19,None,1)
    p = StrandPrimers(primer5, primer3, b)
    enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=16,e_inner=2,k_index=1,k_outer=10,e_outer=1,minIndex=bIndex)
    return enc

def DEC_FSMD_200(pf, primer5, primer3, bindex=0):
    b = commafreecodec.CommaFreeCodec(19,None,1)
    p = StrandPrimers(primer5, primer3, b)
    dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=16,e_inner=2,k_index=1,k_outer=10,e_outer=1,minIndex=bindex)
    return dec

def ENC_RS_CFC8_200(pf, primer5, primer3, bIndex=0):
    b = commafreecodec.CommaFreeCodec(20,None,2)
    p = StrandPrimers(primer5, primer3, b)
    enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=16,e_inner=2,k_index=2,minIndex=bIndex)
    return enc

def DEC_RS_CFC8_200(pf, primer5, primer3, bindex=0):
    b = commafreecodec.CommaFreeCodec(20,None,2)
    p = StrandPrimers(primer5, primer3, b)
    dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=16,e_inner=2,k_index=2,minIndex=bindex)
    return dec

def ENC_XOR_ROT_200(pf, primer5, primer3):
    h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
    p = StrandPrimers(primer5, primer3, h)
    pf.packetSize = 20
    enc = EncodeXORStrand(pf,p)
    return enc
    
def DEC_XOR_ROT_200(pf, primer5, primer3):
    h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
    p = StrandPrimers(primer5, primer3, h)
    pf.packetSize = 20
    dec = DecodeXORStrand(pf,p)
    return dec

def ENC_RS_ROT_200(pf, primer5, primer3, bindex=0):
    h = huffman.RotateCodec(huffman.HuffmanCodec(26,None,16,144))
    p = StrandPrimers(primer5, primer3, h)
    enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=24,e_inner=2,k_index=2,minIndex=bindex)
    # h = huffman.RotateCodec(huffman.HuffmanCodec(11,None,11,99))
    # p = StrandPrimers(primer5, primer3, h)
    # enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=9,e_inner=2,k_index=2,minIndex=bindex)
    return enc

def DEC_RS_ROT_200(pf, primer5, primer3, bindex=0):
    h = huffman.RotateCodec(huffman.HuffmanCodec(26,None,16,144))
    p = StrandPrimers(primer5, primer3, h)
    dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=24,e_inner=2,k_index=2,minIndex=bindex)
    # h = huffman.RotateCodec(huffman.HuffmanCodec(11,None,11,99))
    # p = StrandPrimers(primer5, primer3, h)
    # dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=9,e_inner=2,k_index=2,minIndex=bindex)
    return dec

def ENC_Goldman_200(pf, primer5, primer3):
    h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
    p = StrandPrimers(primer5, primer3, h)
    pf.packetSize = 5
    enc = EncodeGoldmanStrand(pf,4,p)
    return enc

def DEC_Goldman_200(pf, primer5, primer3):
    h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
    p = StrandPrimers(primer5, primer3, h)
    pf.packetSize = 5
    dec = DecodeGoldmanStrand(pf,4,p)
    return dec


# DO NOT ALTER ENTRIES IN THIS TABLE, BUT YOU MAY ADD NEW ONES
# ALL CHANGES NEED TO BE THOROUGHLY TESTED FOR BACKWARDS COMPATIBILITY

FileSystemFormats = {
    # KEY      KEY    LEN  PacketSize, Abbrev.   Description                   ENCODER       DECODR
    0x0010 : [0x0010, 200, 16, "FSMD", "File system meta-data format", ENC_FSMD_200, DEC_FSMD_200 ],
    0x0020 : [0x0020, 200, 16, "RS+CFC8", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_200, DEC_RS_CFC8_200 ],
    #0x0100 : [0x0100, 200, "Dense", "Dense encoding", ENC_Dense_200, DEC_Dense_200 ],
    0x0110 : [0x0110, 200, 5, "Goldman", "Goldman with 4-way repetition", ENC_Goldman_200, DEC_Goldman_200],
    0x0120 : [0x0120, 200, 20, "XOR+ROT", "XOR architecture with Rotating Huffman Encoding",
              ENC_XOR_ROT_200, DEC_XOR_ROT_200],
    0x0030 : [0x0030, 200, 11, "RS+ROT", "Inner/Outer RS with Rotating Huffman encoding", ENC_RS_ROT_200,
              DEC_RS_ROT_200],
    0x1000 : [0x1000, 200, 20, "Segmented", "Segmented file format to support Preview", None, None]
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



