from dnastorage.codec import dense
from dnastorage.codec import norepeatscodec
from dnastorage.codec import commafreecodec
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec.rscodec import *
from dnastorage.arch.strand import *

def available_file_architectures():
    return ["UW+MSv1","Illinois",'Goldman','Fountain','Binary','Dense','NRDense','RS+CFC8','RS+ROT','RS+NRD']

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
        b = commafreecodec.CommaFreeCodec(13,None,2)
        p = StrandPrimers(primer5, primer3, b)
        enc = ReedSolomonInnerOuterEncoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
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
        #print pf.packetSize
        b = commafreecodec.CommaFreeCodec(13,None,2)
        p = StrandPrimers(primer5, primer3, b)
        dec = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        return dec
    elif arch == 'RS+ROT':
        h = huffman.RotateCodec(huffman.HuffmanCodec(11,None,11,99))
        p = StrandPrimers(primer5, primer3, h)
        enc = ReedSolomonInnerOuterDecoder(pf,p,k_datastrand=9,e_inner=2,k_index=2)
        return enc
