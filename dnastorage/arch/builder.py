from dnastorage.codec import dense
from dnastorage.codec import norepeatscodec
from dnastorage.codec.commafreecodec import *
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec import PipeLine as pipeline
from dnastorage.codec.deprecated.rscodec import *
from dnastorage.codec.phys import *
from dnastorage.codec.strand import *
from dnastorage.codec.block import *
from dnastorage.codec.LayeredCodec import *
from dnastorage.arch.strand import *
from dnastorage.exceptions import *
from dnastorage.codec.binarystringcodec import *
from dnastorage.codec.consolidation import *

from dnastorage.fi.probes import *


def available_file_architectures():
    return ["UW+MSv1","Illinois",'Goldman','Fountain','Binary','Dense','NRDense','RS+CFC8','RS+ROT','RS+NRD', 'Sys', 'OverhangStrand']


def customize_RS_CFC8_pipeline(pf,**kwargs):
    print(kwargs)
    blockSizeInBytes=kwargs.get("blockSizeInBytes",150)
    strandSizeInBytes=kwargs.get("strandSizeInBytes",15)
    primer5 = kwargs.get("primer5","")
    primer3 = kwargs.get("primer3","")
    innerECC = kwargs.get("innerECC",0)
    outerECCStrands = kwargs.get("outerECCStrands",0)
    magic_strand = kwargs.get("magic","")
    cut = kwargs.get("cut","")
    fault_injection= kwargs.get("fi",False)
    upper_strand_length = kwargs.get("dna_length",200)
    index_bytes=kwargs.get("index_bytes",0)
    index_bit_set=kwargs.get("index_bit_set",None)

    
    #create the components we are gonna use
    rsOuter = ReedSolomonOuterPipeline(blockSizeInBytes//strandSizeInBytes,outerECCStrands)
    commafree = CommaFreeCodecPipeline(numberBytes=index_bytes+innerECC+strandSizeInBytes)
    rsInner = ReedSolomonInnerCodecPipeline(innerECC)
    magic = PrependSequencePipeline(magic_strand)
    p5 = PrependSequencePipeline(primer5)
    p3 = AppendSequencePipeline(reverse_complement(primer3))
    consolidator = SimpleMajorityVote()

    if fault_injection is False:
        return pipeline.PipeLine((rsOuter,rsInner,commafree,p3,magic,p5),consolidator,blockSizeInBytes,strandSizeInBytes,upper_strand_length,1,packetizedfile=pf,index_bytes=index_bytes
                                 ,index_bit_set=index_bit_set)
    else:
        innerECCprobe = CodewordErrorRateProbe(probe_name="RSInner")
        commafreeprobe = CodewordErrorRateProbe(probe_name="CommaFree")
        return pipeline.PipeLine((rsOuter,innerECCProbe,rsInner,commafreeprobe,commafree,p3,magic,p5),consolidator,
                                 blockSizeInBytes,strandSizeInBytes,upper_strand_length,1,packetizedFile=pf,index_bytes=index_bytes,index_bit_set=index_bit_set)

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


def build_overhang_bitstring_strand(is_enc,pf,primer5,primer3,strand_length,num_overhangs,bits_per_block): #this encoding will build up a representation of the recurisve tree assuming single bit codewords,num overhangs controls substrand sharing and controls the amount of work that can be completed each reaction
    #strand length is passed in terms of bytes per strand
    assert (strand_length*8)%bits_per_block==0 #make sure that the amount of data in a strand is a multiple of the number of bits in a building block
    pol=AllowAll()
    if is_enc==True:
        intraBlockIndex=4 #we want to get up to large files to study sharability of substrands, therefore 4 bytes for indexing
        interBlockIndex=0 #no block indexes all strand indexes
        index=intraBlockIndex+interBlockIndex
        #This is a simple implementation with no error correction built in, simply reads in a strandlength worth of data and then chooses the sequence of overhangs to use
        #resultant strand from the encoder will be overhangs (of length 4), and a 2 length GG or TT sequence to represent 1's and 0's, this can be replaced later to represent actual code words
        blockCodec=DoNothingOuterCodec(packetSize=strand_length,
                                       payloadSize=strand_length,
                                       Policy=pol)
        blockToStrand = BlockToStrand(strand_length,strand_length,Policy=pol,\
                                  intraIndexSize=intraBlockIndex,\
                                  interIndexSize=interBlockIndex)
        strandCodec=ReedSolomonInnerCodec(0,Policy=pol)#inner ecc set to 0, we need to make a long ECC for inner strand protection of long strands
        codewords=BinaryStringCodec(strand_length+index,bits_per_block,Policy=pol)#converts to a binary string
        overhangs=InsertOverhangs(num_overhangs,Policy=pol,CodecObj=codewords) #overhangs inserts overhang sequences in between codewords        
        pre = PrependSequence(primer5,isPrimer=True,Policy=pol)
        app = AppendSequence(reverse_complement(primer3), CodecObj=pre, isPrimer=True)
        physCodec= app #append and prepend primer strings  
        enc=LayeredEncoder(pf,blockSizeInBytes=strand_length,strandSizeInBytes=strand_length,\
                           blockCodec=blockCodec,\
                           blockToStrandCodec=blockToStrand,\
                           strandCodec=strandCodec,\
                           strandToCodewordCodec=overhangs,\
                           codewordToPhysCodec=CombineCodewords(),\
                           physCodec=physCodec,Policy=pol)
        return enc #result will be a a set of binary strings with overhang tags in between


def build_encode_architecture(arch, pf, primer5, primer3,num_overhangs=None):
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
    
    elif arch == 'OverhangBitStringStrand':
        return build_overhang_bitstring_strand(True,pf,primer5,primer3,strand_length=2000,num_overhangs=num_overhangs,bits_per_block=1)#encoding algorithm creates an encoding file with will be interpreted by a compiler to create instructions that will create the strand 
    
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
