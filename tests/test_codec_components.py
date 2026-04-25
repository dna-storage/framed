# -*- coding: utf-8 -*-
"""
tests/test_codec_components.py

Unit tests for individual dnastorage codec components.

Each class is tested in isolation so failures point directly to the
offending component rather than to an integration issue.
"""

import random
from io import BytesIO

import pytest

from dnastorage.strand_representation import BaseDNA
from dnastorage.codec.base_conversion import (
    convertIntToBytes,
    convertBytesToInt,
    pack_bits_to_bytes,
    unpack_bytes_to_indexes,
)
from dnastorage.codec.strand import (
    CRC8,
    CRC8_Index,
    RandomizePayloadPipeline,
    Base4TranscodePipeline,
    ReedSolomonInnerCodecPipeline,
)
from dnastorage.codec.block import ReedSolomonOuterPipeline
from dnastorage.codec.consolidation import SimpleMajorityVote
from dnastorage.util.packetizedfile import WritePacketizedFilestream, ReadPacketizedFilestream


# ---------------------------------------------------------------------------
# BaseDNA
# ---------------------------------------------------------------------------

class TestBaseDNA:
    def test_default_construction(self):
        s = BaseDNA()
        assert s.dna_strand == ""
        assert s.codewords == []
        assert s.index_ints == ()
        assert s.index_bytes == 1

    def test_construction_with_args(self):
        s = BaseDNA(dna_strand="ACGT", codewords=[1, 2, 3],
                    index_ints=(5,), index_bytes=2)
        assert s.dna_strand == "ACGT"
        assert s.codewords == [1, 2, 3]
        assert s.index_ints == (5,)
        assert s.index_bytes == 2

    def test_setters(self):
        s = BaseDNA()
        s.dna_strand = "TTTT"
        s.codewords = [10, 20]
        s.index_ints = (1, 2)
        s.index_bytes = 3
        assert s.dna_strand == "TTTT"
        assert s.codewords == [10, 20]
        assert s.index_ints == (1, 2)
        assert s.index_bytes == 3

    def test_dynamic_attribute_assignment(self):
        # Pipeline sets arbitrary attributes on BaseDNA; ensure that works.
        s = BaseDNA()
        s.before_decode = "ACGT"
        s.is_reversed = False
        assert s.before_decode == "ACGT"
        assert s.is_reversed is False


# ---------------------------------------------------------------------------
# base_conversion utilities
# ---------------------------------------------------------------------------

class TestBaseConversionUtils:
    def test_int_bytes_roundtrip(self):
        for n in [0, 1, 127, 255, 1000, 65535]:
            for nb in [1, 2, 3, 4]:
                if n < (256 ** nb):
                    encoded = convertIntToBytes(n, nb)
                    assert convertBytesToInt(encoded) == n

    def test_pack_unpack_roundtrip_multifield(self):
        index_set = (3, 7, 1)
        bit_sizes = (4, 4, 4)
        packed = pack_bits_to_bytes(index_set, bit_sizes)
        unpacked = unpack_bytes_to_indexes(packed, bit_sizes)
        assert tuple(unpacked) == index_set

    def test_pack_unpack_single_byte(self):
        for v in range(256):
            packed = pack_bits_to_bytes((v,), (8,))
            unpacked = unpack_bytes_to_indexes(packed, (8,))
            assert tuple(unpacked) == (v,)

    def test_pack_unpack_16bit(self):
        for v in [0, 1, 100, 255, 1000, 65535]:
            packed = pack_bits_to_bytes((v,), (16,))
            unpacked = unpack_bytes_to_indexes(packed, (16,))
            assert tuple(unpacked) == (v,)

    def test_pack_unpack_zero_width_field(self):
        # A 0-bit field always has value 0.
        packed = pack_bits_to_bytes((0, 5), (0, 5))
        unpacked = unpack_bytes_to_indexes(packed, (0, 5))
        assert tuple(unpacked) == (0, 5)

    def test_pipeline_index_bit_set_roundtrip(self):
        # Representative bit sizes produced by PipeLine._encode_pipeline.
        # packet_idx=1 bit, sub_packet_idx=5 bits (for 20 sub-packets)
        for pkt in range(2):
            for sub in range(20):
                packed = pack_bits_to_bytes((pkt, sub), (1, 5))
                unpacked = unpack_bytes_to_indexes(packed, (1, 5))
                assert tuple(unpacked) == (pkt, sub)


# ---------------------------------------------------------------------------
# CRC8
# ---------------------------------------------------------------------------

class TestCRC8:
    def test_roundtrip(self):
        codewords = [5, 6, 7, 8, 9, 214]
        s = BaseDNA(codewords=list(codewords))
        crc = CRC8()
        crc.encode(s)
        assert len(s.codewords) == len(codewords) + 1
        crc.decode(s)
        assert s.codewords == codewords

    def test_detects_corruption(self):
        codewords = [5, 6, 7, 8, 9, 214]
        s = BaseDNA(codewords=list(codewords))
        crc = CRC8()
        crc.encode(s)
        s.codewords[2] ^= 0xFF     # flip a byte
        crc.decode(s)
        assert all(v is None for v in s.codewords)

    def test_none_codewords_detected(self):
        # If codewords contain None the decode must mark the strand as invalid.
        # Encode a valid strand first, then inject None into the received codewords.
        codewords = [5, 6, 7, 8, 9, 214]
        s = BaseDNA(codewords=list(codewords))
        crc = CRC8()
        crc.encode(s)
        s.codewords[2] = None   # simulate a dropout in received data
        crc.decode(s)
        assert all(v is None for v in s.codewords)

    def test_all_zeros(self):
        codewords = [0] * 16
        s = BaseDNA(codewords=list(codewords))
        crc = CRC8()
        crc.encode(s)
        crc.decode(s)
        assert s.codewords == codewords

    def test_all_255(self):
        codewords = [255] * 8
        s = BaseDNA(codewords=list(codewords))
        crc = CRC8()
        crc.encode(s)
        crc.decode(s)
        assert s.codewords == codewords


# ---------------------------------------------------------------------------
# CRC8_Index
# ---------------------------------------------------------------------------

class TestCRC8Index:
    def test_roundtrip(self):
        codewords = [5, 6, 7, 8, 9, 214]
        s = BaseDNA(codewords=list(codewords), index_bytes=2)
        crc = CRC8_Index()
        crc.encode(s)
        # CRC inserted right after index bytes → length increases by 1
        assert len(s.codewords) == len(codewords) + 1
        crc.decode(s)
        assert s.codewords == codewords

    def test_index_byte_corruption_detected(self):
        codewords = [5, 6, 7, 8, 9, 214]
        s = BaseDNA(codewords=list(codewords), index_bytes=2)
        crc = CRC8_Index()
        crc.encode(s)
        s.codewords[0] ^= 0xFF   # corrupt first index byte
        crc.decode(s)
        assert all(v is None for v in s.codewords)

    def test_crc_position(self):
        # CRC must be inserted at index_bytes position, not at the end.
        codewords = [10, 20, 30, 40]
        s = BaseDNA(codewords=list(codewords), index_bytes=2)
        crc = CRC8_Index()
        crc.encode(s)
        # After encode: [10, 20, <crc_byte>, 30, 40]
        assert s.codewords[0] == 10
        assert s.codewords[1] == 20
        assert s.codewords[3] == 30
        assert s.codewords[4] == 40


# ---------------------------------------------------------------------------
# RandomizePayloadPipeline
# ---------------------------------------------------------------------------

class TestRandomizePayloadPipeline:
    def test_roundtrip(self):
        codewords = list(range(1, 17))   # 16 bytes
        s = BaseDNA(codewords=list(codewords), index_bytes=3)
        rng = RandomizePayloadPipeline()
        rng.encode(s)
        # Index bytes must be unchanged
        assert s.codewords[:3] == codewords[:3]
        # Payload should be randomized
        rng.decode(s)
        assert s.codewords == codewords

    def test_deterministic_from_index(self):
        # Same index bytes → identical randomization applied
        codewords = [10, 20, 30, 40, 50]
        s1 = BaseDNA(codewords=list(codewords), index_bytes=2)
        s2 = BaseDNA(codewords=list(codewords), index_bytes=2)
        rng = RandomizePayloadPipeline()
        rng.encode(s1)
        rng.encode(s2)
        assert s1.codewords == s2.codewords

    def test_different_indices_differ(self):
        # Different index bytes → different randomization
        codewords = [10, 20, 30, 40, 50, 60, 70, 80]
        s1 = BaseDNA(codewords=list(codewords), index_bytes=2)
        s2 = BaseDNA(codewords=[0, 1] + codewords[2:], index_bytes=2)
        rng = RandomizePayloadPipeline()
        rng.encode(s1)
        rng.encode(s2)
        assert s1.codewords[2:] != s2.codewords[2:]

    def test_none_passthrough(self):
        # If codewords contain None the pipeline must return without crashing.
        s = BaseDNA(codewords=[None, None, 5, 10], index_bytes=2)
        rng = RandomizePayloadPipeline()
        rng.decode(s)   # must not raise


# ---------------------------------------------------------------------------
# Base4TranscodePipeline
# ---------------------------------------------------------------------------

class TestBase4TranscodePipeline:
    def test_roundtrip(self):
        codewords = list(range(16))
        s = BaseDNA(codewords=list(codewords))
        b4 = Base4TranscodePipeline()
        b4.encode(s)
        assert all(c in "ACGT" for c in s.dna_strand)
        b4.decode(s)
        assert s.codewords == codewords

    def test_dna_length(self):
        # Each byte encodes to exactly 4 DNA bases in base-4.
        codewords = [0, 255, 128, 64]
        s = BaseDNA(codewords=list(codewords))
        b4 = Base4TranscodePipeline()
        b4.encode(s)
        assert len(s.dna_strand) == 4 * len(codewords)

    def test_all_zeros(self):
        codewords = [0] * 8
        s = BaseDNA(codewords=list(codewords))
        b4 = Base4TranscodePipeline()
        b4.encode(s)
        b4.decode(s)
        assert s.codewords == codewords

    def test_all_255(self):
        codewords = [255] * 8
        s = BaseDNA(codewords=list(codewords))
        b4 = Base4TranscodePipeline()
        b4.encode(s)
        b4.decode(s)
        assert s.codewords == codewords

    def test_header_roundtrip(self):
        codewords = list(range(12))
        s = BaseDNA(codewords=list(codewords))
        b4 = Base4TranscodePipeline()
        b4.encode(s)
        hdr = b4.encode_header()
        b4_new = Base4TranscodePipeline()
        b4_new.decode_header(hdr)
        assert b4_new.num_codewords == b4.num_codewords == 12


# ---------------------------------------------------------------------------
# ReedSolomonInnerCodecPipeline
# ---------------------------------------------------------------------------

class TestReedSolomonInnerCodecPipeline:
    def test_roundtrip_no_errors(self):
        codewords = list(range(20))
        s = BaseDNA(codewords=list(codewords))
        rs = ReedSolomonInnerCodecPipeline(numberECCBytes=4)
        rs.encode(s)
        assert len(s.codewords) == len(codewords) + 4
        rs.decode(s)
        assert s.codewords == codewords

    def test_corrects_two_errors(self):
        # 4 ECC bytes → can correct 2 errors (t = n_ecc / 2)
        codewords = list(range(20))
        s = BaseDNA(codewords=list(codewords))
        rs = ReedSolomonInnerCodecPipeline(numberECCBytes=4)
        rs.encode(s)
        s.codewords[3] ^= 0xFF
        s.codewords[10] ^= 0xFF
        rs.decode(s)
        assert s.codewords == codewords

    def test_zero_ecc_passthrough(self):
        codewords = list(range(10))
        s = BaseDNA(codewords=list(codewords))
        rs = ReedSolomonInnerCodecPipeline(numberECCBytes=0)
        rs.encode(s)
        assert s.codewords == codewords
        rs.decode(s)
        assert s.codewords == codewords

    def test_varying_ecc_sizes(self):
        for ecc in [2, 6, 8]:
            codewords = list(range(16))
            s = BaseDNA(codewords=list(codewords))
            rs = ReedSolomonInnerCodecPipeline(numberECCBytes=ecc)
            rs.encode(s)
            assert len(s.codewords) == 16 + ecc
            rs.decode(s)
            assert s.codewords == codewords


# ---------------------------------------------------------------------------
# ReedSolomonOuterPipeline (_encode / _decode layer)
# ---------------------------------------------------------------------------

def _make_sub_packets(num_pkts, strands_per_pkt, strand_len, seed=0):
    """Return list-of-lists of BaseDNA with deterministic codewords."""
    rng = random.Random(seed)
    return [
        [BaseDNA(codewords=[rng.randint(0, 255) for _ in range(strand_len)],
                 index_ints=(0, sp, si))
         for si in range(strands_per_pkt)]
        for sp in range(num_pkts)
    ]


class TestReedSolomonOuterPipeline:
    def test_roundtrip_no_erasures(self):
        k, parity, strand_len = 8, 10, 12
        packets = _make_sub_packets(k, strands_per_pkt=1, strand_len=strand_len)
        original = [[list(s.codewords) for s in pkt] for pkt in packets]

        rs = ReedSolomonOuterPipeline(k, parity)
        rs._num_data_sub_packets = k
        rs._total_sub_packets = k + parity
        rs._level = 1

        all_pkts = rs._encode(packets)
        assert len(all_pkts) == k + parity

        decoded = rs._decode(all_pkts)
        for sp in range(k):
            for si, s in enumerate(decoded[sp]):
                assert s.codewords == original[sp][si]

    def test_single_erasure_recovery(self):
        k, parity = 6, 8
        packets = _make_sub_packets(k, strands_per_pkt=1, strand_len=8, seed=42)
        original = [[list(s.codewords) for s in pkt] for pkt in packets]

        rs = ReedSolomonOuterPipeline(k, parity)
        rs._num_data_sub_packets = k
        rs._total_sub_packets = k + parity
        rs._level = 1

        all_pkts = rs._encode(packets)
        # Erase sub-packet 2
        for s in all_pkts[2]:
            s.codewords = [None] * len(s.codewords)

        decoded = rs._decode(all_pkts)
        for sp in range(k):
            for si, s in enumerate(decoded[sp]):
                assert s.codewords == original[sp][si], \
                    f"Mismatch at sub-packet {sp}, strand {si}"

    def test_header_roundtrip(self):
        k, parity = 5, 6
        rs_enc = ReedSolomonOuterPipeline(k, parity)
        rs_enc._num_data_sub_packets = k
        rs_enc._total_sub_packets = k + parity
        rs_enc._level = 1
        rs_enc._index_bits = 4
        rs_enc._zero_range = tuple()

        hdr = rs_enc._encode_header()
        rs_dec = ReedSolomonOuterPipeline(1, 1)   # placeholder values
        rs_dec._num_data_sub_packets = 1
        rs_dec._total_sub_packets = 1
        rs_dec._level = 1
        rs_dec._index_bits = 0
        rs_dec._zero_range = tuple()
        rs_dec._decode_header(hdr)

        assert rs_dec._parity_packets == parity
        assert rs_dec._num_data_sub_packets == k
        assert rs_dec._total_sub_packets == k + parity


# ---------------------------------------------------------------------------
# SimpleMajorityVote
# ---------------------------------------------------------------------------

class TestSimpleMajorityVote:
    def test_single_strand_passthrough(self):
        s = BaseDNA(codewords=[1, 2, 3], index_ints=(0,))
        result = SimpleMajorityVote().decode([s])
        assert len(result) == 1
        assert result[0].codewords == [1, 2, 3]

    def test_majority_vote_three_strands(self):
        s1 = BaseDNA(codewords=[1, 2, 3], index_ints=(0,))
        s2 = BaseDNA(codewords=[1, 2, 3], index_ints=(0,))
        s3 = BaseDNA(codewords=[1, 9, 3], index_ints=(0,))   # byte 1 wrong
        result = SimpleMajorityVote().decode([s1, s2, s3])
        assert len(result) == 1
        assert result[0].codewords == [1, 2, 3]

    def test_two_groups_returned(self):
        s0 = BaseDNA(codewords=[10, 20], index_ints=(0,))
        s1 = BaseDNA(codewords=[30, 40], index_ints=(1,))
        result = SimpleMajorityVote().decode([s0, s1])
        indices = {tuple(r.index_ints) for r in result}
        assert (0,) in indices
        assert (1,) in indices


# ---------------------------------------------------------------------------
# WritePacketizedFilestream / ReadPacketizedFilestream
# ---------------------------------------------------------------------------

class TestWritePacketizedFilestream:
    def test_basic_fill_and_write(self):
        buf = BytesIO()
        pf = WritePacketizedFilestream(buf, 100, 10)
        for i in range(pf.numberOfPackets):
            pf[i] = bytearray([i % 256] * 10)
        assert pf.complete
        pf.write()
        buf.seek(0)
        assert len(buf.read()) == 100

    def test_incomplete_not_complete(self):
        buf = BytesIO()
        pf = WritePacketizedFilestream(buf, 50, 10)
        pf[0] = bytearray(10)    # only 1 of 5 packets
        assert not pf.complete

    def test_missing_keys_reported(self):
        buf = BytesIO()
        pf = WritePacketizedFilestream(buf, 50, 10)
        pf[0] = bytearray(10)
        missing = pf.getMissingKeys()
        assert 0 not in missing
        for k in range(1, 5):
            assert k in missing

    def test_out_of_range_key_ignored(self):
        buf = BytesIO()
        pf = WritePacketizedFilestream(buf, 50, 10)
        pf[999] = bytearray(10)  # way out of range
        assert not pf.complete

    def test_last_packet_truncated(self):
        """Data of 25 bytes with 10-byte packets: last packet is only 5 bytes."""
        buf = BytesIO()
        pf = WritePacketizedFilestream(buf, 25, 10)
        for i in range(pf.numberOfPackets):
            pf[i] = bytearray([i + 1] * 10)
        pf.write()
        buf.seek(0)
        data = buf.read()
        assert len(data) == 25
        # last 5 bytes should come from packet 2 (value = 3)
        assert data[-5:] == bytes([3] * 5)


class TestReadPacketizedFilestream:
    def test_reads_packets(self):
        data = bytes(range(30))
        buf = BytesIO(data)
        pf = ReadPacketizedFilestream(buf)
        pf.packetSize = 10
        packets = [pf.read() for _ in range(3)]
        assert b"".join(packets) == data

    def test_number_of_packets(self):
        data = bytes(range(25))
        buf = BytesIO(data)
        pf = ReadPacketizedFilestream(buf)
        pf.packetSize = 10
        assert pf.numberOfPackets == 3    # ceil(25 / 10)

    def test_number_of_packets_exact(self):
        data = bytes(range(30))
        buf = BytesIO(data)
        pf = ReadPacketizedFilestream(buf)
        pf.packetSize = 10
        assert pf.numberOfPackets == 3    # 30 / 10 exactly

    def test_iteration(self):
        data = bytes(range(40))
        buf = BytesIO(data)
        pf = ReadPacketizedFilestream(buf)
        pf.packetSize = 10
        packets = list(pf)
        assert len(packets) == 4
        assert b"".join(packets) == data

    def test_pad_last_packet_to_packet_size(self):
        """ReadPacketizedFilestream pads the last read to packetSize."""
        buf = BytesIO(bytes(range(15)))
        pf = ReadPacketizedFilestream(buf)
        pf.packetSize = 10
        _ = pf.read()           # packet 0: bytes 0-9
        pkt1 = pf.read()        # packet 1: bytes 10-14 + 5 pad zeros
        assert len(pkt1) == 10
        assert list(pkt1[:5]) == list(range(10, 15))
        assert list(pkt1[5:]) == [0] * 5


# ---------------------------------------------------------------------------
# formats.py registry
# ---------------------------------------------------------------------------

class TestFormatsRegistry:
    def test_formats_not_empty(self):
        from dnastorage.system.formats import file_system_formats
        assert len(file_system_formats()) > 0

    def test_known_format_ids_present(self):
        from dnastorage.system.formats import FileSystemFormats
        assert 0x0704 in FileSystemFormats   # ReedSolomon_Base4_Pipeline
        assert 0x0705 in FileSystemFormats   # Fountain_Base4_Pipeline
        assert 0x0703 in FileSystemFormats   # Basic_Hedges_Pipeline

    def test_encoder_by_abbrev(self):
        from dnastorage.system.formats import file_system_encoder_by_abbrev
        from dnastorage.arch.builder import ReedSolomon_Base4_Pipeline
        enc = file_system_encoder_by_abbrev("ReedSolomon_Base4_Pipeline")
        assert enc is ReedSolomon_Base4_Pipeline

    def test_formatid_by_abbrev(self):
        from dnastorage.system.formats import file_system_formatid_by_abbrev
        assert file_system_formatid_by_abbrev("ReedSolomon_Base4_Pipeline") == 0x0704
        assert file_system_formatid_by_abbrev("Fountain_Base4_Pipeline") == 0x0705

    def test_abbrev_by_id(self):
        from dnastorage.system.formats import file_system_abbrev
        assert file_system_abbrev(0x0704) == "ReedSolomon_Base4_Pipeline"
        assert file_system_abbrev(0x0705) == "Fountain_Base4_Pipeline"

    def test_encoder_decoder_same_function(self):
        from dnastorage.system.formats import (
            FileSystemFormats, file_system_encoder, file_system_decoder,
        )
        for fid in FileSystemFormats:
            assert file_system_encoder(fid) is file_system_decoder(fid)
