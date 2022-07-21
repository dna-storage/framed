from dnastorage.codec import dense
from dnastorage.codec import norepeatscodec
from dnastorage.codec.commafreecodec import *
from dnastorage.codec import illinois
from dnastorage.codec import binary
from dnastorage.codec import huffman
from dnastorage.codec import fountain
from dnastorage.codec import PipeLine as pipeline
import dnastorage.codec.codebooks as codebooks
from dnastorage.codec.deprecated.rscodec import *
from dnastorage.codec.phys import *
from dnastorage.codec.strand import *
from dnastorage.codec.block import *
from dnastorage.codec.LayeredCodec import *
from dnastorage.arch.strand import *
from dnastorage.exceptions import *
from dnastorage.codec.binarystringcodec import *
from dnastorage.codec.consolidation import *
from dnastorage.codec.hedges import *
from dnastorage.codec.cwhedges import *

from dnastorage.fi.probes import *


import logging
logger = logging.getLogger("dnastorage.arch.builder")
logger.addHandler(logging.NullHandler())



def available_file_architectures():
    return ["UW+MSv1","Illinois",'Goldman','Fountain','Binary','Dense','NRDense','RS+CFC8','RS+ROT','RS+NRD', 'Sys', 'OverhangStrand']

def check_required(required, **kwargs):
    missing = []
    for r in required:
        if r not in kwargs:
            missing.append(r)
    title = kwargs.get("title","")
    verb = "is"
    if len(missing) > 1:
        verb = "are"
        missing_s = ", ".join(missing)
    else:
        missing_s = "".join(missing)
    if len(missing) >= 1:
        print ("{} {} missing but required to build a {}.".format(missing_s,verb,title))
        print ("Warning: Hard coded assumptions may not match expectations.")

        
def RS_Codeword_hedges_pipeline(pf,**kwargs):
    blockSizeInBytes=kwargs.get("blockSizeInBytes",150*15)
    strandSizeInBytes=kwargs.get("strandSizeInBytes",15)
    primer5 = kwargs.get("primer5","")
    primer3 = kwargs.get("primer3","")
    outerECCStrands = kwargs.get("outerECCStrands",255-blockSizeInBytes//strandSizeInBytes)
    codebook_func = getattr(codebooks,kwargs.get("codebook","CFC_ALL"))
    syncbook_func = getattr(codebooks,kwargs.get("syncbook","NONE_BOOK"))
    sync_period =  kwargs.get("sync_period",0)
    fault_injection= kwargs.get("fi",False)
    upper_strand_length = kwargs.get("dna_length",208)
    pipeline_title=kwargs.get("title","anonymous_pipeline")
    barcode = kwargs.get("barcode",tuple())
    
    #create the components we are gonna use
    rsOuter = ReedSolomonOuterPipeline(blockSizeInBytes//strandSizeInBytes,outerECCStrands)
    
    #set up the comma free codebook
    codebook = codebook_func()
    syncbook = syncbook_func()
    
    commafree = CodewordHedgesPipeline(codebook,syncbook,sync_period) #use hedges-like decoding method
    
    p5 = PrependSequencePipeline(primer5,handler="align")
    p3 = AppendSequencePipeline(reverse_complement(primer3),handler="align")
    consolidator = SimpleMajorityVote()
    
    if fault_injection is False:
        return pipeline.PipeLine((rsOuter,commafree,p3,p5),consolidator,blockSizeInBytes,strandSizeInBytes,upper_strand_length,1,packetizedfile=pf,
                                 barcode=barcode)
    else:
        innerECCprobe = CodewordErrorRateProbe(probe_name="{}::RSInner".format(pipeline_title))
        commafreeprobe = CodewordErrorRateProbe(probe_name="{}::CommaFree".format(pipeline_title))
        return pipeline.PipeLine((rsOuter,innerECCprobe,commafreeprobe,commafree,p3,p5),consolidator,
                                 blockSizeInBytes,strandSizeInBytes,upper_strand_length,1,packetizedfile=pf,barcode=barcode)

    
def customize_RS_CFC8_pipeline(pf,**kwargs):    
    required = ["blockSizeInBytes","strandSizeInBytes","innerECC",\
                "dna_length"]
    check_required(required,**kwargs)
    
    blockSizeInBytes=kwargs.get("blockSizeInBytes",150*15)
    strandSizeInBytes=kwargs.get("strandSizeInBytes",15)
    primer5 = kwargs.get("primer5","")
    primer3 = kwargs.get("primer3","")
    innerECC = kwargs.get("innerECC",0)
    outerECCStrands = kwargs.get("outerECCStrands",255-blockSizeInBytes//strandSizeInBytes)
    
    cut = kwargs.get("cut","")
    fault_injection= kwargs.get("fi",False)
    upper_strand_length = kwargs.get("dna_length",208)
    pipeline_title=kwargs.get("title","anonymous_pipeline")
    barcode = kwargs.get("barcode",tuple())
    
    #create the components we are gonna use
    rsOuter = ReedSolomonOuterPipeline(blockSizeInBytes//strandSizeInBytes,outerECCStrands)
    commafree = CommaFreeCodecPipeline(numberBytes=innerECC+strandSizeInBytes)
    rsInner = ReedSolomonInnerCodecPipeline(innerECC)
    p5 = PrependSequencePipeline(primer5,handler="align")
    p3 = AppendSequencePipeline(reverse_complement(primer3),handler="align")
    consolidator = SimpleMajorityVote()

    if fault_injection is False:
        return pipeline.PipeLine((rsOuter,rsInner,commafree,p3,p5),consolidator,blockSizeInBytes,strandSizeInBytes,upper_strand_length,1,packetizedfile=pf,
                                 barcode=barcode)
    else:
        innerECCprobe = CodewordErrorRateProbe(probe_name="{}::RSInner".format(pipeline_title))
        commafreeprobe = CodewordErrorRateProbe(probe_name="{}::CommaFree".format(pipeline_title))
        return pipeline.PipeLine((rsOuter,innerECCprobe,rsInner,commafreeprobe,commafree,p3,p5),consolidator,
                                 blockSizeInBytes,strandSizeInBytes,upper_strand_length,1,packetizedfile=pf,barcode=barcode)

    
def SDC_pipeline(pf,**kwargs):
    required = ["blockSizeInBytes","strandSizeInBytes","hedges_rate",\
                "dna_length"]
    check_required(required,**kwargs) 
        
    fault_injection= kwargs.get("fi",False)
    blockSizeInBytes=kwargs.get("blockSizeInBytes",180*15)
    strandSizeInBytes=kwargs.get("strandSizeInBytes",15)
    primer5 = kwargs.get("primer5",'A'*20)
    primer3 =kwargs.get("primer3",'A'*20)
    t7_seq = kwargs.get("T7","CGACTAATACGACTCACTATAGC")
    rt_pcr_seq =kwargs.get("RT-PCR","CGCTAGCTCTAGAGATCTAG")
    check_primers = kwargs.get("check_primers",False)
    other_strands=kwargs.get("other_strands",[])
    
    hedges_rate = kwargs.get("hedges_rate",1/2.)
    # pad_bits and prev_bits should match by default:
    hedges_pad_bits=kwargs.get("hedges_pad",8)
    hedges_previous = kwargs.get("hedge_prev_bits",8)

    if "outerECCStrands" not in kwargs and "blockSizeInBytes" in kwargs:
        outerECCStrands = 255 - blockSizeInBytes//strandSizeInBytes
    else:
        outerECCStrands = kwargs.get("outerECCStrands",255-180)
        
    upper_strand_length = kwargs.get("dna_length",300)
    pipeline_title=kwargs.get("title","")
    barcode = kwargs.get("barcode",tuple())
    
    #Error correction components
    if "outerECCdivisor" not in kwargs:
        rsOuter = ReedSolomonOuterPipeline(blockSizeInBytes//strandSizeInBytes,outerECCStrands)
    elif "outerECCdivisor" in kwargs:
        rsOuter = ReedSolomonOuterPipeline(kwargs["outerECCdivisor"],outerECCStrands)
  
    hedges = FastHedgesPipeline(rate=hedges_rate,pad_bits=hedges_pad_bits,prev_bits=hedges_previous)
    crc = CRC8()
    
    #components related to DNA functionality
    p5 = PrependSequencePipeline(primer5,ignore=True)
    t7 = PrependSequencePipeline(t7_seq,ignore=True)
    rt_pcr = PrependSequencePipeline(rt_pcr_seq,handler="align",search_range=70)
    p3 = AppendSequencePipeline(reverse_complement(primer3),handler="align",search_range=20)

    consolidator = SimpleMajorityVote()

    out_pipeline = (rsOuter,)
    inner_pipeline = (crc,hedges)
    DNA_pipeline = (p3,rt_pcr,t7,p5)

    if fault_injection: #some counters for data collection
        index_probe = IndexDistribution(probe_name=pipeline_title,prefix_to_match=barcode)
        hedges_probe = CodewordErrorRateProbe(probe_name="{}::hedges".format(pipeline_title))
        inner_pipeline = (index_probe,crc,hedges_probe,hedges)
        dna_counter_probe = FilteredDNACounter(probe_name=pipeline_title)
        DNA_pipeline=(dna_counter_probe,)+DNA_pipeline

    if check_primers: #checks data strands for matches in 
        primer_check_probe = StrandCheckProbe(strands=[primer5,t7_seq,rt_pcr_seq,primer3]+other_strands) 
        DNA_pipeline=DNA_pipeline+(primer_check_probe,)

    return pipeline.PipeLine(out_pipeline+inner_pipeline+DNA_pipeline,consolidator,blockSizeInBytes,strandSizeInBytes,upper_strand_length,1,packetizedfile=pf,
                            barcode=barcode)



def Basic_Hedges_Pipeline(pf,**kwargs):
    required = ["blockSizeInBytes","strandSizeInBytes","hedges_rate",\
                "dna_length", "crc_type", "reverse_payload"]
    check_required(required,**kwargs) 
        
    fault_injection= kwargs.get("fi",False)
    blockSizeInBytes=kwargs.get("blockSizeInBytes",180*15)
    strandSizeInBytes=kwargs.get("strandSizeInBytes",15)
    primer5 = kwargs.get("primer5",'A'*20) #basic 5' primer region
    primer3 =kwargs.get("primer3",'A'*20) #basic 3' primer region
    crc_type = kwargs.get("crc_type","strand") #location where to put the CRC check
    reverse_payload = kwargs.get("reverse_payload",False) #should the payload be reversed
    upper_strand_length = kwargs.get("dna_length",300) #if strands longer than this value, throw exceptions
    pipeline_title=kwargs.get("title","") #name for the pipeline
    barcode = kwargs.get("barcode",tuple()) #specific barcode for the pipeline

    hedges_rate = kwargs.get("hedges_rate",1/2.) #rate of hedges encode/decode
    # pad_bits and prev_bits should match by default:
    hedges_pad_bits=kwargs.get("hedges_pad",8)
    hedges_previous = kwargs.get("hedge_prev_bits",8)
    try_reverse = kwargs.get("try_reverse",False) #should the reverse be tried, for hedges decoding, set this if primer3/primer5 are not used
    

    if "outerECCStrands" not in kwargs and "blockSizeInBytes" in kwargs:
        outerECCStrands = 255 - blockSizeInBytes//strandSizeInBytes
    else:
        outerECCStrands = kwargs.get("outerECCStrands",255-180)
  
    
    #Error correction components
    if "outerECCdivisor" not in kwargs:
        rsOuter = ReedSolomonOuterPipeline(blockSizeInBytes//strandSizeInBytes,outerECCStrands)
    elif "outerECCdivisor" in kwargs:
        rsOuter = ReedSolomonOuterPipeline(kwargs["outerECCdivisor"],outerECCStrands)
  
    hedges = FastHedgesPipeline(rate=hedges_rate,pad_bits=hedges_pad_bits,prev_bits=hedges_previous,try_reverse=try_reverse)

    if crc_type=="strand": crc = CRC8()
    elif crc_type=="index":
        logger.info("BasicHedges: Using Index CRC")
        crc = CRC8_Index()
    else: assert 0 and "Invalid CRC selection"
        
    #components related to DNA functionality
    p5 = PrependSequencePipeline(primer5,ignore=False,handler="align",search_range=30)
    p3 = PrependSequencePipeline(primer3,ignore=False,handler="align",search_range=30)

    consolidator = SimpleMajorityVote()
    out_pipeline = (rsOuter,)
    inner_pipeline = (crc,hedges)
    
    if reverse_payload:
        logger.info("BasicHedges: Using reverse after payload DNA")
        r = ReversePipeline()
        DNA_pipeline = (r,p3,p5)
    else:
        DNA_pipeline = (p3,p5)
        
    if fault_injection: #some counters for data collection
        index_probe = IndexDistribution(probe_name=pipeline_title,prefix_to_match=barcode)
        hedges_probe = CodewordErrorRateProbe(probe_name="{}::hedges".format(pipeline_title))
        inner_pipeline = (index_probe,crc,hedges_probe,hedges)
        dna_counter_probe = FilteredDNACounter(probe_name=pipeline_title)
        DNA_pipeline=(dna_counter_probe,)+DNA_pipeline

   
    return pipeline.PipeLine(out_pipeline+inner_pipeline+DNA_pipeline,consolidator,blockSizeInBytes,strandSizeInBytes,upper_strand_length,1,packetizedfile=pf,
                            barcode=barcode)


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
