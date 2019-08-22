from dnastorage.codec import dense
from dnastorage.codec import norepeatscodec
from dnastorage.codec.commafreecodec import *
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec.deprecated.rscodec import *
from dnastorage.codec.phys import *
from dnastorage.codec.strand import *
from dnastorage.codec.block import *
from dnastorage.codec.LayeredCodec import *
from dnastorage.arch.strand import *
from dnastorage.exceptions import *

def available_file_architectures():
    return ["UW+MSv1","Illinois",'Goldman','Fountain','Binary','Dense','NRDense','RS+CFC8','RS+ROT','RS+NRD', 'Sys']


def customize_RS_CFC8(is_enc,pf,primer5,primer3,intraBlockIndex=1,\
                      interBlockIndex=2,innerECC=2,strandSizeInBytes=15,\
                      blockSizeInBytes=15*185,Policy=None,\
                      withCut=None,outerECCStrands=None,minIndex=0):
    assert blockSizeInBytes % strandSizeInBytes == 0
    payload=strandSizeInBytes
    blockStrands = blockSizeInBytes / strandSizeInBytes
    if outerECCStrands is None:
        outerECCStrands = 255-blockStrands # strands
    else:
        assert outerECCStrands + blockStrands <= 255
    assert blockStrands + outerECCStrands < 256

    index = intraBlockIndex + interBlockIndex

    blockCodec = ReedSolomonOuterCodec(packetSize=blockSizeInBytes,\
                                       errorSymbols=outerECCStrands,payloadSize=payload,Policy=Policy)

    blockToStrand = BlockToStrand(payload,(blockStrands+outerECCStrands)*payload,Policy=Policy,\
                                  intraIndexSize=intraBlockIndex,\
                                  interIndexSize=interBlockIndex,filterZeroes=True)
    
    # take index into account
    # better to just say number of error symbols
    strandCodec = ReedSolomonInnerCodec(innerECC,Policy=Policy)

    codewords = CommaFreeCodewords(payload+innerECC+index,Policy=Policy)

    if not (withCut is None):
        cut = InsertMidSequence(withCut)
    else:
        cut = None
    
    pre = PrependSequence(primer5,isPrimer=True,Policy=Policy,CodecObj=cut)
    app = AppendSequence(reverse_complement(primer3), CodecObj=pre, isPrimer=True, \
                         Policy=Policy)    
    physCodec = app

    if is_enc==True:
        enc = LayeredEncoder(pf,blockSizeInBytes=blockSizeInBytes,\
                             strandSizeInBytes=strandSizeInBytes,\
                         blockCodec=blockCodec,\
                         blockToStrandCodec=blockToStrand,\
                         strandCodec=strandCodec,\
                         strandToCodewordCodec=codewords,\
                         codewordToPhysCodec=CombineCodewords(),\
                             physCodec=physCodec,Policy=Policy,minIndex=minIndex)
        return enc
    else:
        dec = LayeredDecoder(pf,blockSizeInBytes=blockSizeInBytes,\
                             strandSizeInBytes=strandSizeInBytes,\
                             blockIndexSize=interBlockIndex,\
                             blockCodec=blockCodec,\
                             strandCodec=strandCodec,\
                             physCodec=physCodec,\
                             physToStrandCodec=codewords,\
                             strandToBlockCodec=blockToStrand,Policy=Policy,minIndex=minIndex)
        return dec
        

def build_RS_CFC8(is_encoder,pf,primer5,primer3):
    
    pol = AllowAll()

    intraBlockIndex = 1
    interBlockIndex = 2
    inner_ecc = 2
    payload=15
    baseBlock = 185
    outerECC = 255-baseBlock # strands
    assert baseBlock + outerECC < 256

    index = intraBlockIndex + interBlockIndex

    blockSize = payload*baseBlock # 185 normal strands
    
    blockCodec = ReedSolomonOuterCodec(packetSize=blockSize,\
                                       errorSymbols=70,payloadSize=payload,Policy=pol)

    blockToStrand = BlockToStrand(payload,(baseBlock+outerECC)*payload,Policy=pol,\
                                  intraIndexSize=intraBlockIndex,\
                                  interIndexSize=interBlockIndex)
    
    # take index into account
    # better to just say number of error symbols
    strandCodec = ReedSolomonInnerCodec(inner_ecc,Policy=pol)

    codewords = CommaFreeCodewords(payload+inner_ecc+index,Policy=pol)
    pre = PrependSequence(primer5,isPrimer=True,Policy=pol)
    app = AppendSequence(reverse_complement(primer3), CodecObj=pre, isPrimer=True)    
    physCodec = app

    if is_encoder==True:
        enc = LayeredEncoder(pf,blockSizeInBytes=blockSize,strandSizeInBytes=payload,\
                             blockCodec=blockCodec,\
                             blockToStrandCodec=blockToStrand,\
                             strandCodec=strandCodec,\
                             strandToCodewordCodec=codewords,\
                             codewordToPhysCodec=CombineCodewords(),\
                             physCodec=physCodec,Policy=pol)
        return enc
    else:
        dec = LayeredDecoder(pf,blockSizeInBytes=blockSize,strandSizeInBytes=payload,
                             blockCodec=blockCodec,\
                             strandCodec=strandCodec,\
                             physCodec=physCodec,\
                             physToStrandCodec=codewords,\
                             strandToBlockCodec=blockToStrand,Policy=pol)
        return dec
    


def build_encode_architecture(arch, pf, primer5, primer3):
    if arch == "UW+MSv1":
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        ## build packetized file encode
        pf.packetSize = 20
        enc = EncodeXORStrand(pf,p)
        return enc

    elif arch == "Illinois":
        illini = illinois.IllinoisCodec(primer5,150,Checksum())
        p = StrandPrimers(primer5, primer3, illini)
        pf.packetSize = 149
        enc = EncodeNaiveStrand(pf,p)
        return enc

    elif arch == 'Goldman':
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 5
        enc = EncodeGoldmanStrand(pf,4,p)
        return enc

    elif arch == 'Fountain':
        #assert False and "Not fully implemented"
        h = huffman.RotateCodec(huffman.HuffmanCodec(23,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 22
        enc = EncodeFountainStrand(pf,1.5,p)
        return enc

    elif arch == 'Binary':
        b = binary.BinaryRotateCodec(binary.BinaryCodec(Checksum()))
        #b = binary.BinaryCodec(Checksum())
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 17
        enc = EncodeNaiveStrand(pf,p)
        return enc

    elif arch == 'Dense':
        #b = dense.RandomizeDenseCodec(dense.DenseCodec(Checksum()))
        b = dense.DenseCodec(Checksum())
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 34
        enc = EncodeNaiveStrand(pf,p)
        return enc

    elif arch == 'NRDense':
        b = norepeatscodec.NoRepeatCodec(Checksum(2))
        #b = dense.RandomizeDenseCodec(dense.DenseCodec(Checksum()))
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 28
        enc = EncodeNaiveStrand(pf,p)
        return enc

    elif arch == 'RS+NRD':
        b = norepeatscodec.NoRepeatCodec()
        #b = dense.RandomizeDenseCodec(dense.DenseCodec(Checksum()))
        p = StrandPrimers(primer5, primer3, b)
        #pf.packetSize = 28
        enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=16,e_inner=8,k_index=4)
        return enc

    elif arch == 'RS+CFC8':
        return build_RS_CFC8(True,pf,primer5,primer3)
        #b = commafreecodec.CommaFreeCodec(13,None,2)
        #p = StrandPrimers(primer5, primer3, b)
        #enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        #return enc

    elif arch == 'FSMD':
        b = CommaFreeCodec(19,None,1)
        p = StrandPrimers(primer5, primer3, b)
        enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=16,e_inner=2,k_index=1,e_outer=1)
        return enc

    elif arch == 'RS+ROT':
        h = huffman.RotateCodec(huffman.HuffmanCodec(11,None,11,99))
        p = StrandPrimers(primer5, primer3, h)
        enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        return enc

def build_decode_architecture(arch, pf, primer5, primer3, fountain_table=None):
    if arch == "UW+MSv1":
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 20
        dec = DecodeXORStrand(pf,p)
        return dec

    elif arch == "Illinois":
        illini = illinois.IllinoisCodec(primer5,150,Checksum())
        p = StrandPrimers(primer5, primer3, illini)
        pf.packetSize = 149
        enc = DecodeNaiveStrand(pf,p)
        return enc

    elif arch == 'Goldman':
        h = huffman.RotateCodec(huffman.HuffmanCodec(21,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 5
        dec = DecodeGoldmanStrand(pf,4,p)
        return dec

    elif arch == 'Binary':
        b = binary.BinaryRotateCodec(binary.BinaryCodec(Checksum()))
        #b = binary.BinaryCodec(Checksum())
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 17
        enc = DecodeNaiveStrand(pf,p)
        return enc

    elif arch == 'Dense':
        b = dense.DenseCodec(Checksum())
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 34
        dec = DecodeNaiveStrand(pf,p)
        return dec

    elif arch == 'NRDense':
        b = norepeatscodec.NoRepeatCodec(Checksum(2,CheckLength(28)))
        p = StrandPrimers(primer5, primer3, b)
        pf.packetSize = 28
        dec = DecodeNaiveStrand(pf,p)
        return dec

    elif arch == 'Fountain':
        #assert False and "Not fully implemented"
        h = huffman.RotateCodec(huffman.HuffmanCodec(23,Checksum()))
        p = StrandPrimers(primer5, primer3, h)
        pf.packetSize = 22
        enc = DecodeFountainStrand(pf,fountain_table,p)
        return enc

    elif arch == 'NCState':
        assert False and "Not fully implemented"
        return None

    elif arch == 'RS+NRD':
        b = norepeatscodec.NoRepeatCodec()
        #b = dense.RandomizeDenseCodec(dense.DenseCodec(Checksum()))
        p = StrandPrimers(primer5, primer3, b)
        #pf.packetSize = 28
        dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=16,e_inner=8,k_index=4)
        return dec

    elif arch == 'RS+CFC8':
        return build_RS_CFC8(False,pf,primer5,primer3)
        #print pf.packetSize
        #b = commafreecodec.CommaFreeCodec(13,None,2)
        #p = StrandPrimers(primer5, primer3, b)
        #dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        #return dec
        
    elif arch == 'RS+ROT':
        h = huffman.RotateCodec(huffman.HuffmanCodec(11,None,11,99))
        p = StrandPrimers(primer5, primer3, h)
        enc = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        return enc

    elif arch == 'FSMD':
        b = CommaFreeCodec(19,None,1)
        p = StrandPrimers(primer5, primer3, b)
        enc = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=16,e_inner=2,k_index=1,e_outer=1)
        return enc
