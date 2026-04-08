from dnastorage.codec.base import *
from math import log, ceil
from dnastorage.codec import base_conversion
from dnastorage.codec.reedsolomon.rs import ReedSolomon,get_reed_solomon,ReedSolomonError
from dnastorage.codec.fountain.lt import LTCode
from collections import Counter
from dnastorage.strand_representation import *
from dnastorage.codec_types import *
import logging
logger = logging.getLogger('dna.storage.codec.block')
logger.addHandler(logging.NullHandler())


class ReedSolomonOuterPipeline(BaseOuterCodec):
    def __init__(self,packet_divisor,parity_packets,c_exp=8,OuterCodecObj=None,Policy=None):
        BaseOuterCodec.__init__(self,packet_divisor,OuterCodecObj=OuterCodecObj,Policy=Policy)
        self._rs = get_reed_solomon(c_exp=c_exp)
        self._parity_packets=parity_packets
        assert(packet_divisor+parity_packets<=self._rs.field_charac)
    def _encode(self,packets):
        parity_packets=[]
        #initialize parity_packets
        for i in range(0,self._parity_packets):
            new= []
            parity_packets.append(new)
        #encode a set of packets
        for i in range(0,len(packets[0])):
            strand_set=[]
            for j in range(0,len(packets)):
                strand_set.append(packets[j][i])
            #use strand set to calculate a message
            assert(len(strand_set)>0)
            strand_length = len(strand_set[0].codewords)
            for i in range(0,self._parity_packets): parity_packets[i].append(BaseDNA(codewords=[])) #make a new strand in each parity packet
            #go column by column, calculating a RS message and splitting the ECC across strands in each packet
            for k in range(0,strand_length):
                message=[]
                for j in range(0,len(strand_set)):
                    message.append(strand_set[j].codewords[k]) #get bytes
                assert(len(message)>0)
                mesecc = self._rs.rs_encode_msg(message, self._parity_packets)
                ecc = mesecc[len(message):]
                for ei, e in enumerate(ecc):
                    parity_packets[ei][-1].codewords.append(e)
            for p in parity_packets:
                assert len(p[-1].codewords)==strand_length
        return packets+parity_packets
    
    def _decode(self,packets):
        #decode a set of packets 
        #Need to get packets into the right order to unroll the encoding, basically will get messages the same way as encoding
        assert len(packets)>0
        for i in range(0,len(packets[0])):
            strands=[]
            for j in range(0,self._total_sub_packets):
                strands.append(packets[j][i])
            #now get messages
            assert len(strands)>0
            for byte_index in range(0,len(strands[0].codewords)):
                message=[]
                for strand_index in range(0,len(strands)):
                    message.append(strands[strand_index].codewords[byte_index])
                assert len(message)==self._total_sub_packets
                try:
                    #perform correction
                    erasures = [i for i in range(0,len(message)) if message[i]==None]
                    corrected_message, corrected_ecc = \
                        self._rs.rs_correct_msg(message, \
                                                self._parity_packets, \
                                                erase_pos=erasures)
                except ReedSolomonError as e:
                    #just gonna assume we can't fix this, hopefully someone else will
                    #stats.inc("RSOuterCodec::ReedSolomonError")
                    wr_e = DNAReedSolomonOuterCodeError(msg=\
                             "RSOuterCodec found error at sub-packet index".format(strands[0].index_ints[:self._level]))
                    corrected_message = message
                except Exception as e:
                    if self._Policy.allow(e):
                        corrected_message = message                        
                        pass
                    else:
                        raise e

                #write the corrected message back into strand classes
                for s,m in zip(strands,corrected_message):
                    s.codewords[byte_index] = m
        return packets


class FountainOuterPipeline(BaseOuterCodec):
    """
    Outer erasure code based on a systematic LT (Luby Transform) fountain code.

    Each sub-packet corresponds to one LT source symbol (a list of strands at a
    given sub-packet index).  Encoding produces `parity_packets` additional parity
    sub-packets whose strand codewords are XOR combinations of source strands
    chosen via the robust soliton distribution.  Decoding runs a belief-propagation
    peeling decoder to recover any erased source sub-packets.

    Unlike Reed-Solomon this code is not bounded by GF field size, so it works for
    any number of source sub-packets.

    Parameters
    ----------
    packet_divisor : int
        Passed directly to BaseOuterCodec; controls how the strand packet is split
        into sub-packets.
    parity_packets : int
        Number of LT parity sub-packets to generate.
    seed : int
        PRNG seed for deterministic Tanner graph construction.  Must match between
        encoder and decoder (stored in the encoded header).
    """

    def __init__(self, packet_divisor, parity_packets, seed=42, OuterCodecObj=None, Policy=None):
        BaseOuterCodec.__init__(self, packet_divisor, OuterCodecObj=OuterCodecObj, Policy=Policy)
        self._parity_packets = parity_packets
        self._seed = seed
        # LTCode is constructed lazily in _encode/_decode once k is known.
        self._lt = None

    def _get_lt(self, k):
        """Return (or create) the LTCode instance for k source symbols."""
        if self._lt is None or self._lt.k != k:
            self._lt = LTCode(k, seed=self._seed)
        return self._lt

    # ------------------------------------------------------------------
    # Header serialisation (extends BaseOuterCodec header)
    # ------------------------------------------------------------------

    def _encode_header(self):
        buf = BaseOuterCodec._encode_header(self)
        buf += convertIntToBytes(self._seed, 4)
        buf += convertIntToBytes(self._parity_packets, 2)
        return buf

    def _decode_header(self, buff):
        buff = BaseOuterCodec._decode_header(self, buff)
        self._seed = convertBytesToInt(buff[0:4])
        self._parity_packets = convertBytesToInt(buff[4:6])
        self._lt = None  # reset so it will be rebuilt with new seed
        return buff[6:]

    # ------------------------------------------------------------------
    # Encoding
    # ------------------------------------------------------------------

    def _encode(self, packets):
        """
        Generate parity sub-packets and append them to `packets`.

        Each sub-packet is a list of BaseDNA strands.  The LT code treats one
        sub-packet as one source symbol.  Codeword bytes are XORed column-by-column
        across strands (matching the RS pattern so sub-packet sizes stay consistent).
        """
        num_data = len(packets)
        lt = self._get_lt(num_data)
        strand_count = len(packets[0])
        strand_length = len(packets[0][0].codewords)

        # Build parity sub-packet skeletons.  Codewords start empty and are filled
        # in below when each strand position is processed.
        parity_packets = [[BaseDNA(codewords=[]) for _ in range(strand_count)]
                          for _ in range(self._parity_packets)]

        # Process one strand position at a time (column across sub-packets)
        for strand_pos in range(strand_count):
            # Collect one bytearray per source sub-packet for this strand_pos
            source_symbols = [
                bytearray(packets[sp][strand_pos].codewords)
                for sp in range(num_data)
            ]

            # Ask LT code to produce parity symbols
            parity_syms = lt.encode(source_symbols, self._parity_packets)

            for pi, encoded_ba, _nbrs in parity_syms:
                parity_packets[pi][strand_pos].codewords = list(encoded_ba)

        # Verify lengths are consistent
        for pp in parity_packets:
            for s in pp:
                assert len(s.codewords) == strand_length

        return packets + parity_packets

    # ------------------------------------------------------------------
    # Decoding
    # ------------------------------------------------------------------

    def _decode(self, packets):
        """
        Recover erased source sub-packets using the LT peeling decoder.

        `packets` has length `_total_sub_packets` (data + parity).  Erased strands
        carry None in their codeword slots.  The decoder writes corrected bytes back
        in-place and returns the full packet list unchanged.
        """
        assert len(packets) > 0
        num_data = self._num_data_sub_packets
        lt = self._get_lt(num_data)
        strand_count = len(packets[0])
        strand_length = len(packets[0][0].codewords)

        # Process one strand position at a time
        for strand_pos in range(strand_count):
            # Build the received-symbol list for the LT decoder.
            # Format: (symbol_index, bytearray_or_None, [source_neighbors])
            received = []

            # Source sub-packets (indices 0..num_data-1)
            for sp in range(num_data):
                strand = packets[sp][strand_pos]
                cws = strand.codewords
                # A sub-packet strand is erased if ANY byte is None
                if any(b is None for b in cws):
                    received.append((sp, None, []))
                else:
                    received.append((sp, bytearray(cws), []))

            # Parity sub-packets (indices num_data..total-1)
            for pi in range(self._parity_packets):
                sp = num_data + pi
                if sp >= len(packets):
                    break
                strand = packets[sp][strand_pos]
                cws = strand.codewords
                nbrs = lt.neighbors(pi)
                if any(b is None for b in cws):
                    received.append((sp, None, nbrs))
                else:
                    received.append((sp, bytearray(cws), nbrs))

            # Run the peeling decoder
            try:
                decoded = lt.decode(received, strand_length)
            except Exception as e:
                if self._Policy.allow(e):
                    decoded = [None] * num_data
                else:
                    raise e

            if decoded is None:
                decoded = [None] * num_data

            # Write recovered bytes back into source sub-packet strands
            for sp in range(num_data):
                strand = packets[sp][strand_pos]
                recovered = decoded[sp]
                if recovered is not None:
                    strand.codewords = list(recovered)
                # If still None leave codewords as-is (already None-filled by base class)

        return packets

