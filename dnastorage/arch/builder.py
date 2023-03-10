from dnastorage.codec import PipeLine as pipeline #import pipeline backend
from dnastorage.exceptions import * #import special exceptions
from dnastorage.fi.probes import * #probes import


#muscle alignment import
from dnastorage.alignment.muscle import * 

#imports for consolidators
from dnastorage.codec.DNAConsolidatemodels import *
from dnastorage.cluster.lsh import *
from dnastorage.cluster.ideal_cluster import *
from dnastorage.codec.consolidation import *

#import codec components used in the inner code 
from dnastorage.codec.phys import *
from dnastorage.codec.strand import *
from dnastorage.codec.hedges import *

#import outer code
from dnastorage.codec.block import *


import logging
logger = logging.getLogger("dnastorage.arch.builder")
logger.addHandler(logging.NullHandler())


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

def ReedSolomon_Base4_Pipeline(pf,**kwargs):
    required = ["blockSizeInBytes","strandSizeInBytes","hedges_rate",\
                "dna_length", "crc_type", "reverse_payload"]
    check_required(required,**kwargs)
    index_bytes = kwargs.get("index_bytes",None) #optional, sets index bytes to be constant, if the constant number is less than the actual, encoder will fail
    primer5 = kwargs.get("primer5",'A'*20) #basic 5' primer region
    primer3 =kwargs.get("primer3",'A'*20) #basic 3' primer region
    #components related to DNA functionality
    p5 = PrependSequencePipeline(primer5,ignore=False,handler="align",search_range=100)
    p3 = AppendSequencePipeline(primer3,ignore=False,handler="align",search_range=100)
    #set up some fault injection/sequencing options
    fault_injection= kwargs.get("fi",False)
    #set strand size and block size bytes
    blockSizeInBytes=kwargs.get("blockSizeInBytes",180*15)
    strandSizeInBytes=kwargs.get("strandSizeInBytes",15)
    #set some pipeline characteristics
    pipeline_title=kwargs.get("title","") #name for the pipeline
    barcode = kwargs.get("barcode",tuple()) #specific barcode for the pipeline
    #set inner RS parameters
    inner_ECC=kwargs.get("inner_ECC",0)

    if "packeted_inner_strand_size" in kwargs: #allow for packeting to reduce parameter waste 
        inner_ECC,strandSizeInBytes=kwargs["packeted_inner_strand_size"]


    using_DNA_consolidator = kwargs.get("using_DNA_consolidator","lsh")
    lsh_m_sigs = kwargs.get("lsh_m",50) #number of signatures to make per strand
    lsh_kmer = kwargs.get("lsh_k",5) #kmer length to build a signature out of
    lsh_sim = kwargs.get("lsh_sim",0.5) #simularity of strands
    lsh_sig_samples = kwargs.get("lsh_sig_samples",4) #samples to make on signatures
    lsh_sample_length = kwargs.get("lsh_sample_length",100000)  #length of strand to consider when hashing
    align_num_strands = kwargs.get("align_num_strands",15)
    cw_consolidator = SimpleMajorityVote()
    if using_DNA_consolidator:
        if using_DNA_consolidator=="lsh":
            cluster = LocalitySensitiveHashCluster(lsh_m_sigs,lsh_kmer,lsh_sig_samples,int(1/(lsh_sim**lsh_sig_samples)),lsh_sample_length)
        if using_DNA_consolidator=="ideal":
            cluster = IdealCluster()
        align = MuscleAlign(align_num_strands)
        dna_consolidator = BasicDNAClusterModel(cluster,align,name=pipeline_title)
  

    if "outerECCStrands" not in kwargs and "blockSizeInBytes" in kwargs:
        outerECCStrands = 255 - blockSizeInBytes//strandSizeInBytes
    else:
        outerECCStrands = kwargs.get("outerECCStrands",255-180)

    if outerECCStrands>0:
        #Error correction components
        if "outerECCdivisor" not in kwargs:
            rsOuter = ReedSolomonOuterPipeline(blockSizeInBytes//strandSizeInBytes,outerECCStrands)
        elif "outerECCdivisor" in kwargs:
            rsOuter = ReedSolomonOuterPipeline(kwargs["outerECCdivisor"],outerECCStrands)
        out_pipeline = (rsOuter,)
    else:
        out_pipeline=(BaseOuterCodec(int(math.ceil(blockSizeInBytes/strandSizeInBytes))),)
    
    randomize = RandomizePayloadPipeline()
    base4codec = Base4TranscodePipeline()
    RS_inner = ReedSolomonInnerCodecPipeline(inner_ECC)    
    inner_pipeline = (randomize,RS_inner,base4codec)
    
    if fault_injection: #some counters for data collection
        index_probe = IndexDistribution(probe_name=pipeline_title,prefix_to_match=barcode)
        RS_probe = CodewordErrorRateProbe(probe_name="{}::inner_rs".format(pipeline_title))
        Base4Probe = CodewordErrorRateProbe(probe_name="{}::base4".format(pipeline_title))
        RandomizeProbe = CodewordErrorRateProbe(probe_name="{}::randomize".format(pipeline_title))
        dna_hook_probe = HookProbe("dna_strand",RS_probe.name)
        DNA_error_probe = DNAErrorProbe(probe_name=pipeline_title)
        post_cluster_DNA_probe = DNAErrorProbe(probe_name="{}::post_cluster".format(pipeline_title))
        post_cluster_DNA_probe.dna_attr = dna_hook_probe.name 
        inner_pipeline = (index_probe,RandomizeProbe,randomize,RS_probe,RS_inner,Base4Probe,post_cluster_DNA_probe,base4codec)
        DNA_pipeline =(DNA_error_probe,dna_hook_probe,p5,p3)
    upper_strand_length=400
    return pipeline.PipeLine(out_pipeline+inner_pipeline+DNA_pipeline,blockSizeInBytes,strandSizeInBytes,upper_strand_length,1,packetizedfile=pf,
                             barcode=barcode,dna_consolidator=dna_consolidator,cw_consolidator=cw_consolidator,constant_index_bytes=index_bytes)

def Basic_Hedges_Pipeline(pf,**kwargs):
    required = ["blockSizeInBytes","strandSizeInBytes","hedges_rate",\
                "dna_length", "crc_type", "reverse_payload"]
    check_required(required,**kwargs) 

    #set up some fault injection/sequencing options
    fault_injection= kwargs.get("fi",False)
    sequencing_run=kwargs.get("seq",False)

    #set strand size and block size bytes
    blockSizeInBytes=kwargs.get("blockSizeInBytes",180*15)
    strandSizeInBytes=kwargs.get("strandSizeInBytes",15)
    index_bytes = kwargs.get("index_bytes",None) #optional, sets index bytes to be constant, if the constant number is less than the actual, encoder will fail

    #primers to append/prepend to the strand
    primer5 = kwargs.get("primer5",'A'*20) #basic 5' primer region
    primer3 =kwargs.get("primer3",'A'*20) #basic 3' primer region
    crc_type = kwargs.get("crc_type","strand") #location where to put the CRC check
    reverse_payload = kwargs.get("reverse_payload",False) #should the payload be reversed
    
    #set lengths for checking strands that are encoded/decoded
    dna_length = kwargs.get("dna_length",300) #if strands longer than this value, throw exceptions (encode)
    filter_upper_length = kwargs.get("filter_upper_length",float('inf'))
    filter_lower_length = kwargs.get("filter_lower_length",float('-inf'))

    #set some pipeline characteristics
    pipeline_title=kwargs.get("title","") #name for the pipeline
    barcode = kwargs.get("barcode",tuple()) #specific barcode for the pipeline

    #set hedges parameters
    hedges_rate = kwargs.get("hedges_rate",1/2.) #rate of hedges encode/decode
    hedges_pad_bits=kwargs.get("hedges_pad",8)
    hedges_previous = kwargs.get("hedge_prev_bits",8)
    hedges_guesses = kwargs.get("hedges_guesses",100000)
    try_reverse = kwargs.get("try_reverse",False) #should the reverse be tried, for hedges decoding, set this if primer3/primer5 are not used
    

    #consolidator options
    using_DNA_consolidator = kwargs.get("using_DNA_consolidator",False)
    lsh_m_sigs = kwargs.get("lsh_m",50) #number of signatures to make per strand
    lsh_kmer = kwargs.get("lsh_k",5) #kmer length to build a signature out of
    lsh_sim = kwargs.get("lsh_sim",0.5) #simularity of strands
    lsh_sig_samples = kwargs.get("lsh_sig_samples",4) #samples to make on signatures
    lsh_sample_length = kwargs.get("lsh_sample_length",100000)  #length of strand to consider when hashing
    align_num_strands = kwargs.get("align_num_strands",15)
    cw_consolidator = SimpleMajorityVote()


    if "packeted_inner_strand_size" in kwargs: #allow for packeting to reduce parameter waste 
        hedges_rate,strandSizeInBytes=kwargs["packeted_inner_strand_size"]
    
    if using_DNA_consolidator:
        if using_DNA_consolidator=="lsh":
            cluster = LocalitySensitiveHashCluster(lsh_m_sigs,lsh_kmer,lsh_sig_samples,int(1/(lsh_sim**lsh_sig_samples)),lsh_sample_length)
        if using_DNA_consolidator=="ideal":
            cluster = IdealCluster()
        align = MuscleAlign(align_num_strands)
        dna_consolidator = BasicDNAClusterModel(cluster,align,name=pipeline_title)
    else:
        dna_consolidator=None

    if "outerECCStrands" not in kwargs and "blockSizeInBytes" in kwargs:
        outerECCStrands = 255 - blockSizeInBytes//strandSizeInBytes
    else:
        outerECCStrands = kwargs.get("outerECCStrands",255-180)

    if outerECCStrands>0:
        #Error correction components
        if "outerECCdivisor" not in kwargs:
            rsOuter = ReedSolomonOuterPipeline(blockSizeInBytes//strandSizeInBytes,outerECCStrands)
        elif "outerECCdivisor" in kwargs:
            rsOuter = ReedSolomonOuterPipeline(kwargs["outerECCdivisor"],outerECCStrands)
        out_pipeline = (rsOuter,)
    else:
        out_pipeline=(BaseOuterCodec(int(math.ceil(blockSizeInBytes/strandSizeInBytes))),)

  
    hedges = FastHedgesPipeline(rate=hedges_rate,pad_bits=hedges_pad_bits,prev_bits=hedges_previous,try_reverse=try_reverse,guess_limit=hedges_guesses)

    if crc_type=="strand": crc = CRC8()
    elif crc_type=="index":
        logger.info("BasicHedges: Using Index CRC")
        crc = CRC8_Index()
    else: assert 0 and "Invalid CRC selection"
        
    #components related to DNA functionality
    p5 = PrependSequencePipeline(primer5,ignore=False,handler="align",search_range=100)
    p3 = AppendSequencePipeline(primer3,ignore=False,handler="align",search_range=100)
    length_filter = DNALengthFilterPipeline(filter_lower_length,filter_upper_length)
    

    using_randomize = kwargs.get("randomize",False)
    randomize = RandomizePayloadPipeline()
    if using_randomize: inner_pipeline = (crc,randomize,hedges)
    else: inner_pipeline=(crc,hedges)
    
    if reverse_payload:
        logger.info("BasicHedges: Using reverse after payload DNA")
        r = ReversePipeline()
        DNA_pipeline = (r,p3,p5)
    else:
        DNA_pipeline = (p3,p5)
    
    DNA_pipeline = (length_filter,)+DNA_pipeline #add filter to end of DNAtoDNA decoding
    
    if fault_injection: #some counters for data collection
        post_cluster_DNA_probe = DNAErrorProbe(probe_name="{}::post_cluster".format(pipeline_title))            
        index_probe = IndexDistribution(probe_name=pipeline_title,prefix_to_match=barcode)
        hedges_probe = CodewordErrorRateProbe(probe_name="{}::hedges".format(pipeline_title))
        if not using_randomize: inner_pipeline = (index_probe,crc,hedges_probe,post_cluster_DNA_probe,hedges)
        else: inner_pipeline = (index_probe,crc,randomize,hedges_probe,post_cluster_DNA_probe,hedges)
        dna_counter_probe = FilteredDNACounter(probe_name=pipeline_title)
        DNA_pipeline=(dna_counter_probe,)+DNA_pipeline
        if sequencing_run:
            dna_hook_probe = HookProbe("dna_strand",hedges_probe.name)
            post_cluster_DNA_probe.dna_attr = dna_hook_probe.name 
            length_filter.alignment_name=dna_hook_probe.name
            hedges_probe.dna_attr=dna_hook_probe.name
            DNA_error_probe = DNAErrorProbe(probe_name=pipeline_title)
            DNA_pipeline =(DNA_error_probe,dna_hook_probe)+DNA_pipeline
        
    return pipeline.PipeLine(out_pipeline+inner_pipeline+DNA_pipeline,blockSizeInBytes,strandSizeInBytes,dna_length,1,packetizedfile=pf,
                             barcode=barcode,cw_consolidator=cw_consolidator,dna_consolidator=dna_consolidator,constant_index_bytes=index_bytes)

