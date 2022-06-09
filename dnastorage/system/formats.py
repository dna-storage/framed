from dnastorage.codec import dense
from dnastorage.codec import norepeatscodec
from dnastorage.codec import commafreecodec
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.arch.strand import *
from dnastorage.arch.builder import customize_RS_CFC8, customize_RS_CFC8_pipeline
from dnastorage.arch.builder import build_overhang_bitstring_strand, SDC_pipeline, build_hedges_rs
from dnastorage.exceptions import *


# FIXME: All of these formats need to implement the LayeredCodec.

def ENC_OH_BITSTRING_XXX(*args, **kwargs): #encoding used to experimentally evaluate overhang construction
    assert len(args) >= 3
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    strand_length = kwargs.get('strand_length',512)
    num_overhangs = kwargs.get('num_overhangs',3)
    bits_per_block = kwargs.get('bits_per_block',1)    
    enc = build_overhang_bitstring_strand(True,pf,primer5,primer3,strand_length,num_overhangs,bits_per_block)
    return enc

def DEC_OH_BITSTRING_XXX(pf, primer5,primer3,strand_length,num_overhangs,bits_per_block):
    assert 0 #not implemented
    
def ENC_FSMD_200(*args, **kwargs):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",NoTolerance())
    withCut = kwargs.get("withCut",None)
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,0,2,15,90,policy,withCut=withCut,outerECCStrands=20,minIndex=bIndex)
    return enc    

def DEC_FSMD_200(pf, primer5, primer3, bIndex=0,policy=AllowAll(),withCut=None):
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,0,2,15,90,policy,withCut=withCut,outerECCStrands=20,minIndex=bIndex)
    return dec

def ENC_FSMD_WCUT_160(pf, primer5, primer3, bIndex=0, policy=NoTolerance(),withCut=None):
    withCut="AGGTACCA"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,0,2,10,120,policy,withCut=withCut,outerECCStrands=20,minIndex=bIndex)
    return enc    

def DEC_FSMD_WCUT_160(pf, primer5, primer3, bIndex=0,policy=AllowAll(),withCut=None):
    withCut="AGGTACCA"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,0,2,10,120,policy,withCut=withCut,outerECCStrands=20,minIndex=bIndex)
    return dec

def ENC_RS_CFC8_200(*args, **kwargs):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",NoTolerance())
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,2,2,15,185*15,policy,minIndex=bIndex)
    return enc    

def DEC_RS_CFC8_200(*args, **kwargs):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",AllowAll())
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,2,2,15,185*15,policy,minIndex=bIndex)
    return dec    


def ENC_RS_CFC8_RE1_160(*args, **kwargs):#
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",NoTolerance())   
    withCut="AGGTACCA"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,
                            withCut=withCut)
    return enc    

def DEC_RS_CFC8_RE1_160(*args, **kwargs):#(pf, primer5, primer3, bIndex=0, policy=AllowAll()):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",AllowAll())
    withCut="AGGTACCA"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,\
                            withCut=withCut)
    return dec    

def ENC_RS_CFC8_RE2_160(*args, **kwargs):#(pf, primer5, primer3, bIndex=0, policy=NoTolerance()):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",NoTolerance())   
    withCut="CCTGCAGG"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,
                            withCut=withCut)
    return enc    

def DEC_RS_CFC8_RE2_160(*args, **kwargs):#(pf, primer5, primer3, bIndex=0, policy=AllowAll()):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",AllowAll())
    withCut="CCTGCAGG"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,\
                            withCut=withCut)
    return dec    

def ENC_RS_CFC8_RE3_160(*args, **kwargs):#(pf, primer5, primer3, bIndex=0, policy=NoTolerance()):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",NoTolerance())   
    withCut="GCGGCCGC"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,
                            withCut=withCut)
    return enc    

def DEC_RS_CFC8_RE3_160(*args, **kwargs):#(pf, primer5, primer3, bIndex=0, policy=AllowAll()):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",AllowAll())
    withCut="GCGGCCGC"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,\
                            withCut=withCut)
    return dec    

def ENC_RS_CFC8_RE4_160(*args, **kwargs):#(pf, primer5, primer3, bIndex=0, policy=NoTolerance()):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",NoTolerance())   
    withCut="GTTTAAAC"
    enc = customize_RS_CFC8(True,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,
                            withCut=withCut)
    return enc    

def DEC_RS_CFC8_RE4_160(*args, **kwargs):#(pf, primer5, primer3, bIndex=0, policy=AllowAll()):
    pf = args[0]
    primer5 = args[1]
    primer3 = args[2]
    bIndex = kwargs.get("bIndex",0)
    policy = kwargs.get("policy",AllowAll())
    withCut="GTTTAAAC"
    dec = customize_RS_CFC8(False,pf,primer5,primer3,1,1,3,9,185*9,policy,minIndex=bIndex,\
                            withCut=withCut)
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
    assert False
    
    #h = huffman.RotateCodec(huffman.HuffmanCodec(26,None,16,144))
    #p = StrandPrimers(primer5, primer3, h)
    #enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=24,e_inner=2,k_index=2,minIndex=bindex)
    # h = huffman.RotateCodec(huffman.HuffmanCodec(11,None,11,99))
    # p = StrandPrimers(primer5, primer3, h)
    # enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=9,e_inner=2,k_index=2,minIndex=bindex)
    return enc

def DEC_RS_ROT_200(pf, primer5, primer3, bindex=0):
    assert False
    #h = huffman.RotateCodec(huffman.HuffmanCodec(26,None,16,144))
    #p = StrandPrimers(primer5, primer3, h)
    #dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=24,e_inner=2,k_index=2,minIndex=bindex)
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


#KV: Added pipelines
def PIPE_Hedges_RS(pf,**kwargs):    
    #Pipeline implementation of Hedges and RS outer code, customization values are embedded in kwargs
    primer5 = kwargs.get("primer5",'A'*20)
    primer3 = kwargs.get("primer3",'A'*20)    
    pipe = build_hedges_rs(pf,\
                        blockSizeInBytes=180*8,\
                        strandSizeInBytes=8,\
                        outerECCStrands=255-180,\
                        hedges_rate=1.0/4,\
                        hedges_pad=8,\
                        hedges_prev_bits=8,\
                        dna_length=300,\
                        primer5=primer5,\
                        primer3=primer3,\
                        title='Hedges+RS')
    return pipe

#KV: Added pipelines
def CUSTOM_PIPE_Hedges_RS(pf,**kwargs):
    title = kwargs.get("title","CustomHedges+RS")
    #Pipeline implementation of Hedges and RS outer code, customization values are embedded in kwargs
    pipe = build_hedges_rs(pf,title=title,**kwargs)
    return pipe


def PIPE_250_FSMD(pf,**kwargs):
    blockSizeInBytes=90
    outerECCStrands=40
    strandSizeInBytes=15
    #pipeline implementation fo FSMD
    pipe = customize_RS_CFC8_pipeline(pf,**kwargs,innerECC=2,strandSizeInBytes=strandSizeInBytes,\
                                      blockSizeInBytes=blockSizeInBytes,outerECCStrands=outerECCStrands,\
                                      dna_length=250,title="FSMD-Pipe")
    return pipe

def PIPE_RS_CFC8(pf,**kwargs):
    #pipe = customize_RS_CFC8_pipeline(pf,outerECC=37,innerECC=3,dna_length=208,**kwargs)
    primer5 = kwargs.get("primer5",'A'*20)
    primer3 = kwargs.get("primer3",'A'*20)    
    pipe = customize_RS_CFC8_pipeline(pf,\
                                      innerECC=3,\
                                      blockSizeInBytes=150*15,\
                                      strandSizeInBytes=15,\
                                      outerECCStrands=255-150,\
                                      dna_length=208,\
                                      title='Pipe-RS+CFC8',\
                                      primer5=primer5,\
                                      primer3=primer3)
    return pipe


def CUSTOM_PIPE_RS_CFC8(pf,**kwargs):
    title = kwargs.get("title","CustomPipe-RS+CFC8")
    kwargs2 = kwargs.copy()
    kwargs2.pop("title","")
    #pipe = customize_RS_CFC8_pipeline(pf,outerECC=37,innerECC=3,dna_length=208,**kwargs)        
    pipe = customize_RS_CFC8_pipeline(pf,title=title,**kwargs2)
    return pipe

def SDC_PIPE(pf,**kwargs):
    # only let primers through to influence design
    primer5 = kwargs.get("primer5",'A'*20)
    primer3 = kwargs.get("primer3",'A'*20)
    pipe = SDC_pipeline(pf,\
                        blockSizeInBytes=180*6,\
                        strandSizeInBytes=6,\
                        outerECCStrands=255-180,\
                        #hedges_rate=1.0/4,\
                        hedges_pad=8,\
                        hedges_prev_bits=8,\
                        dna_length=300,\
                        title="SDC",\
                        primer5=primer5,\
                        primer3=primer3\
                        )

    #pipe = customize_RS_CFC8_pipeline(pf,outerECC=37,innerECC=3,dna_length=208,**kwargs)
    return pipe

def CUSTOM_SDC_PIPE(pf,**kwargs):
    title = kwargs.get("title","CustomSDC")
    return SDC_pipeline(pf,title=title,**kwargs)


# DO NOT ALTER ENTRIES IN THIS TABLE, BUT YOU MAY ADD NEW ONES
# ALL CHANGES NEED TO BE THOROUGHLY TESTED FOR BACKWARDS COMPATIBILITY

FileSystemFormats = {
    #------ Meta-data
    # KEY      KEY    LEN  PacketSize, Abbrev.   Description                   ENCODER       DECODR
    0x0010 : [0x0010, 200, 90, "FSMD", "File system meta-data format", ENC_FSMD_200, DEC_FSMD_200 ],
    0x0011 : [0x0011, 160, 120, "FSMD-1", "File system meta-data format with cut",
              ENC_FSMD_WCUT_160, DEC_FSMD_WCUT_160 ],
    0x0012 : [0x0012, 250, 200, "FSMD-Pipe", "File system meta-data format pipelined",
              PIPE_250_FSMD,PIPE_250_FSMD],

    #------ ^^^ Meta-data only above here

    #------ First codecs used in experiments based on CFC and RS, LayeredCodec
    0x0020 : [0x0020, 200, 16, "RS+CFC8", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_200, DEC_RS_CFC8_200 ],
    #0x0100 : [0x0100, 200, "Dense", "Dense encoding", ENC_Dense_200, DEC_Dense_200 ],
    #0x0110 : [0x0110, 200, 5, "Goldman", "Goldman with 4-way repetition", ENC_Goldman_200, DEC_Goldman_200],
    #0x0030 : [0x0030, 200, 11, "RS+ROT", "Inner/Outer RS with Rotating Huffman encoding", ENC_RS_ROT_200,
    #          DEC_RS_ROT_200],
    #0x0040 : [0x0040, 200, 20, "XOR+ROT", "XOR architecture with Rotating Huffman Encoding",
    #          ENC_XOR_ROT_200, DEC_XOR_ROT_200],

    #------ Pipelines
    0x0100 : [0x0100, 208, 15, "Pipe-RS+CFC8","RS+CFC8 implemented with pipelines",PIPE_RS_CFC8,PIPE_RS_CFC8],
    0x0400 : [0x0400, 300, 0, "Hedges+RS",
              "Hedges with ReedSolomon Outer Code, pipeline implementation, fully customizable",
              PIPE_Hedges_RS,PIPE_Hedges_RS],    
    0x0500 : [0x0500, 208, 15, "SDC","Final format for SDC experiments",SDC_PIPE,SDC_PIPE],
    #------ ^^ Above are hardcoded, below are flexible VV    
    0x0700 : [0x0700, 500, 15, "CustomPipe-RS+CFC8","Customizable RS+CFC8 implemented with pipelines",
              CUSTOM_PIPE_RS_CFC8, CUSTOM_PIPE_RS_CFC8],
    0x0701 : [0x0701, 208, 15, "CustomSDC","Support for SDC experiments",\
              CUSTOM_SDC_PIPE,CUSTOM_SDC_PIPE],
    0x0702 : [0x0702, 300, 0, "CustomHedges+RS",
              "Hedges with ReedSolomon Outer Code, pipeline implementation, fully customizable",
              CUSTOM_PIPE_Hedges_RS,CUSTOM_PIPE_Hedges_RS],    
    
    #------ Segmented
    0x2000 : [0x2000, 200, 20, "Segmented", "Segmented file format to support Preview", None, None],
    0x2021 : [0x2021, 160, 9, "RS+CFC8+RE1", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_RE1_160, DEC_RS_CFC8_RE1_160 ],
    0x2022 : [0x2022, 160, 9, "RS+CFC8+RE2", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_RE2_160, DEC_RS_CFC8_RE2_160 ],
    0x2023 : [0x2023, 160, 9, "RS+CFC8+RE3", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_RE3_160, DEC_RS_CFC8_RE3_160 ],
    0x2024 : [0x2024, 160, 9, "RS+CFC8+RE4", "Reed-Solomon coded with Comma-free codewords",
              ENC_RS_CFC8_RE4_160, DEC_RS_CFC8_RE4_160 ],

    #------ Overhang
    0x3000 : [0x3000, "Variable", "Variable", "OH_BITSTRING_XXX",
              "Experimental encoding to evaluate overhang construction",
              ENC_OH_BITSTRING_XXX,DEC_OH_BITSTRING_XXX],
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


