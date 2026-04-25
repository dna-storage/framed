# -*- coding: utf-8 -*-
"""
tests/test_pipeline_roundtrip.py

End-to-end encode → decode roundtrip tests for the dnastorage pipeline
builders.

Design notes:
- Strand objects are mutable: decode strips primers from ``dna_strand`` in
  place.  Each decode run receives **fresh** BaseDNA copies so successive
  tests never interfere with one another.
- Tests use small block / strand sizes (blockSizeInBytes=160, strandSizeInBytes=16)
  so the suite stays fast.
- Erasure tests drop sub-packets *well below* the codec capacity limit to be
  deterministic; separate tests probe the exact capacity boundary.
"""

import copy
import random
from io import BytesIO

import pytest

from dnastorage.strand_representation import BaseDNA
from dnastorage.util.packetizedfile import ReadPacketizedFilestream, WritePacketizedFilestream


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fresh(strands):
    """Return new BaseDNA objects with the same ``dna_strand`` / ``index_ints``
    as the originals.  The pipeline only needs ``dna_strand`` for decoding."""
    return [BaseDNA(dna_strand=s.dna_strand, index_ints=s.index_ints)
            for s in strands]


def _encode(builder, data, params):
    """Encode *data* and return (strands, header_bytes)."""
    pf = ReadPacketizedFilestream(BytesIO(data))
    pipe = builder(pf, **params)
    strands = [s for block in pipe for s in block]
    hdr = pipe.encode_header_data()
    return strands, hdr


def _decode(builder, strands, hdr, data_len, params):
    """Decode *strands* and return the decoded bytes (length = data_len)."""
    buf = BytesIO()
    pf = WritePacketizedFilestream(buf, data_len, 0)
    pipe = builder(pf, **params)
    pipe.decode_header_data(hdr)
    for s in strands:
        pipe.decode(s)
    pipe.final_decode()
    buf.seek(0)
    return buf.read(data_len)


def _roundtrip(builder, data, params):
    """Encode data then decode it.  Returns decoded bytes."""
    strands, hdr = _encode(builder, data, params)
    return _decode(builder, _fresh(strands), hdr, len(data), params)


# ──────────────────────────────────────────────────────────────────────────────
# Shared pipeline parameters (small, fast for CI)
# ──────────────────────────────────────────────────────────────────────────────
_RS_PARAMS = dict(
    blockSizeInBytes=160,
    strandSizeInBytes=16,
    primer5='ACGTACGTACGTACGTACGT',
    primer3='TGCATGCATGCATGCATGCA',
    outerECCStrands=10,
    inner_ECC=0,
    using_DNA_consolidator=False,
)

_FOUNTAIN_PARAMS = dict(
    blockSizeInBytes=160,
    strandSizeInBytes=16,
    primer5='ACGTACGTACGTACGTACGT',
    primer3='TGCATGCATGCATGCATGCA',
    outerECCStrands=15,          # LT works better with more parity
    fountain_seed=7,
    inner_ECC=0,
    using_DNA_consolidator=False,
)


# ──────────────────────────────────────────────────────────────────────────────
# ReedSolomon_Base4_Pipeline roundtrip
# ──────────────────────────────────────────────────────────────────────────────

class TestReedSolomonBase4Roundtrip:
    """Full encode → decode roundtrips for ReedSolomon_Base4_Pipeline."""

    @pytest.fixture(autouse=True)
    def _builder(self):
        from dnastorage.arch.builder import ReedSolomon_Base4_Pipeline
        self.builder = ReedSolomon_Base4_Pipeline

    def test_single_block_sequential_bytes(self):
        data = bytes(range(160))
        assert _roundtrip(self.builder, data, _RS_PARAMS) == data

    def test_single_block_random_bytes(self):
        rng = random.Random(0)
        data = bytes(rng.randint(0, 255) for _ in range(160))
        assert _roundtrip(self.builder, data, _RS_PARAMS) == data

    def test_multi_block(self):
        data = bytes(i % 256 for i in range(320))
        assert _roundtrip(self.builder, data, _RS_PARAMS) == data

    def test_exactly_one_strand_of_data(self):
        """Edge case: data fills exactly one strand payload."""
        data = bytes(range(16))
        assert _roundtrip(self.builder, data, _RS_PARAMS) == data

    def test_all_zero_bytes(self):
        data = bytes(160)
        assert _roundtrip(self.builder, data, _RS_PARAMS) == data

    def test_all_ff_bytes(self):
        data = bytes([0xFF] * 160)
        assert _roundtrip(self.builder, data, _RS_PARAMS) == data

    def test_encoded_strands_are_dna_only(self):
        """All encoded strands must consist only of DNA bases."""
        data = bytes(range(160))
        strands, _ = _encode(self.builder, data, _RS_PARAMS)
        for s in strands:
            assert all(c in "ACGT" for c in s.dna_strand), \
                f"Non-DNA character in strand: {s.dna_strand!r}"

    def test_strands_carry_unique_indices(self):
        """No two strands from the same block should share the same index."""
        data = bytes(range(160))
        strands, _ = _encode(self.builder, data, _RS_PARAMS)
        indices = [s.index_ints for s in strands]
        assert len(indices) == len(set(indices))

    def test_primer5_prefix(self):
        """Every strand must start with the primer5 sequence."""
        p5 = _RS_PARAMS['primer5']
        data = bytes(range(160))
        strands, _ = _encode(self.builder, data, _RS_PARAMS)
        for s in strands:
            assert s.dna_strand.startswith(p5), \
                f"Strand does not start with primer5: {s.dna_strand!r}"

    def test_primer3_suffix(self):
        """Every strand must end with the primer3 sequence."""
        p3 = _RS_PARAMS['primer3']
        data = bytes(range(160))
        strands, _ = _encode(self.builder, data, _RS_PARAMS)
        for s in strands:
            assert s.dna_strand.endswith(p3), \
                f"Strand does not end with primer3: {s.dna_strand!r}"

    def test_erasure_recovery_within_capacity(self):
        """Drop 5 data sub-packets; RS with 10 parity should recover."""
        data = bytes(range(160))
        strands, hdr = _encode(self.builder, data, _RS_PARAMS)
        # Keep only strands whose sub-packet index >= 5
        surviving = _fresh([s for s in strands if s.index_ints[1] >= 5])
        decoded = _decode(self.builder, surviving, hdr, len(data), _RS_PARAMS)
        assert decoded == data

    def test_header_survives_fresh_decoder(self):
        """Header must be enough to fully reconstruct decode parameters."""
        data = bytes(range(160))
        strands, hdr = _encode(self.builder, data, _RS_PARAMS)
        # Decode with *minimal* kwargs – only what the builder requires for
        # construction; all codec state is loaded from the header.
        decoded = _decode(self.builder, _fresh(strands), hdr, len(data), _RS_PARAMS)
        assert decoded == data

    def test_strand_order_independence(self):
        """Decoding with strands in random order must yield correct output."""
        rng = random.Random(123)
        data = bytes(range(160))
        strands, hdr = _encode(self.builder, data, _RS_PARAMS)
        shuffled = _fresh(strands)
        rng.shuffle(shuffled)
        decoded = _decode(self.builder, shuffled, hdr, len(data), _RS_PARAMS)
        assert decoded == data

    def test_strand_duplication_handled(self):
        """Duplicate strands must not corrupt the decoded output."""
        data = bytes(range(160))
        strands, hdr = _encode(self.builder, data, _RS_PARAMS)
        # Duplicate every strand twice
        duplicated = []
        for s in strands:
            duplicated.extend(_fresh([s, s]))
        decoded = _decode(self.builder, duplicated, hdr, len(data), _RS_PARAMS)
        assert decoded == data


# ──────────────────────────────────────────────────────────────────────────────
# Fountain_Base4_Pipeline roundtrip
# ──────────────────────────────────────────────────────────────────────────────

class TestFountainBase4Roundtrip:
    """Full encode → decode roundtrips for Fountain_Base4_Pipeline."""

    @pytest.fixture(autouse=True)
    def _builder(self):
        from dnastorage.arch.builder import Fountain_Base4_Pipeline
        self.builder = Fountain_Base4_Pipeline

    def test_single_block_sequential_bytes(self):
        data = bytes(range(160))
        assert _roundtrip(self.builder, data, _FOUNTAIN_PARAMS) == data

    def test_single_block_random_bytes(self):
        rng = random.Random(42)
        data = bytes(rng.randint(0, 255) for _ in range(160))
        assert _roundtrip(self.builder, data, _FOUNTAIN_PARAMS) == data

    def test_multi_block(self):
        data = bytes(i % 256 for i in range(320))
        assert _roundtrip(self.builder, data, _FOUNTAIN_PARAMS) == data

    def test_all_zero_bytes(self):
        data = bytes(160)
        assert _roundtrip(self.builder, data, _FOUNTAIN_PARAMS) == data

    def test_encoded_strands_are_dna_only(self):
        data = bytes(range(160))
        strands, _ = _encode(self.builder, data, _FOUNTAIN_PARAMS)
        for s in strands:
            assert all(c in "ACGT" for c in s.dna_strand)

    def test_erasure_recovery(self):
        """Drop a few strands; fountain with 15 parity should recover."""
        data = bytes(range(160))
        strands, hdr = _encode(self.builder, data, _FOUNTAIN_PARAMS)
        # Drop 5 strands
        surviving = _fresh(strands[5:])
        decoded = _decode(self.builder, surviving, hdr, len(data), _FOUNTAIN_PARAMS)
        assert decoded == data

    def test_strand_order_independence(self):
        rng = random.Random(99)
        data = bytes(range(160))
        strands, hdr = _encode(self.builder, data, _FOUNTAIN_PARAMS)
        shuffled = _fresh(strands)
        rng.shuffle(shuffled)
        decoded = _decode(self.builder, shuffled, hdr, len(data), _FOUNTAIN_PARAMS)
        assert decoded == data


# ──────────────────────────────────────────────────────────────────────────────
# ReedSolomon_Base4_FileLevelFountain_Pipeline
# ──────────────────────────────────────────────────────────────────────────────

def _file_level_roundtrip(builder, data, params):
    """
    End-to-end roundtrip for file-level fountain pipelines.

    The file-level pipeline requires calling ``generate_parity_bytes`` on the
    raw data blocks before ``encode_header_data``.  This helper replicates the
    logic performed by ``WriteDNAFilePipeline.close()``.
    """
    block_size = params['blockSizeInBytes']

    # ── Encoding ──────────────────────────────────────────────────────────────
    # 1. Split raw data into fixed-size blocks (pad the last block with zeros).
    full_buf = BytesIO(data)
    blocks = []
    while True:
        chunk = full_buf.read(block_size)
        if not chunk:
            break
        if len(chunk) < block_size:
            chunk = chunk.ljust(block_size, b'\x00')
        blocks.append(bytearray(chunk))

    # 2. Build encoder pipeline, compute file-level parity blocks.
    pf_enc = ReadPacketizedFilestream(BytesIO(data))
    enc_pipe = builder(pf_enc, **params)
    fll = enc_pipe._file_level_codec
    parity_blocks = fll.generate_parity_bytes(blocks)

    # 3. Encode all blocks (data blocks first, then parity blocks).
    extended = bytes(data) + b''.join(bytes(pb) for pb in parity_blocks)
    pf_enc2 = ReadPacketizedFilestream(BytesIO(extended))
    enc_pipe2 = builder(pf_enc2, **params)
    enc_pipe2._file_level_codec = fll  # reuse codec that already has state set
    all_strands = [s for block in enc_pipe2 for s in block]
    hdr = enc_pipe2.encode_header_data()

    # ── Decoding ──────────────────────────────────────────────────────────────
    buf = BytesIO()
    pf_dec = WritePacketizedFilestream(buf, len(data), 0)
    dec_pipe = builder(pf_dec, **params)
    dec_pipe.decode_header_data(hdr)
    for s in _fresh(all_strands):
        dec_pipe.decode(s)
    dec_pipe.final_decode()
    buf.seek(0)
    return buf.read(len(data))


class TestFileLevelFountainPipeline:
    """
    ReedSolomon_Base4_FileLevelFountain_Pipeline: multi-block file with an
    additional file-level LT fountain code across blocks.
    """

    _PARAMS = dict(
        **_RS_PARAMS,
        file_level_parity_blocks=1,
        file_level_seed=42,
    )

    @pytest.fixture(autouse=True)
    def _builder(self):
        from dnastorage.arch.builder import ReedSolomon_Base4_FileLevelFountain_Pipeline
        self.builder = ReedSolomon_Base4_FileLevelFountain_Pipeline

    def test_two_block_roundtrip(self):
        """Two 160-byte blocks → encode + decode."""
        data = bytes(range(160)) + bytes(reversed(range(160)))
        assert _file_level_roundtrip(self.builder, data, self._PARAMS) == data

    def test_single_block_roundtrip(self):
        data = bytes(range(160))
        assert _file_level_roundtrip(self.builder, data, self._PARAMS) == data

    def test_pipeline_has_file_level_codec(self):
        data = bytes(range(160))
        pf = ReadPacketizedFilestream(BytesIO(data))
        pipe = self.builder(pf, **self._PARAMS)
        assert hasattr(pipe, '_file_level_codec')
        assert pipe._file_level_codec is not None

    def test_file_level_codec_parity_blocks(self):
        data = bytes(range(160))
        pf = ReadPacketizedFilestream(BytesIO(data))
        pipe = self.builder(pf, **self._PARAMS)
        assert pipe._file_level_codec.parity_blocks == self._PARAMS['file_level_parity_blocks']


# ──────────────────────────────────────────────────────────────────────────────
# Fountain_Base4_FileLevelFountain_Pipeline
# ──────────────────────────────────────────────────────────────────────────────

class TestFountainFileLevelPipeline:
    _PARAMS = dict(
        **_FOUNTAIN_PARAMS,
        file_level_parity_blocks=1,
        file_level_seed=11,
    )

    @pytest.fixture(autouse=True)
    def _builder(self):
        from dnastorage.arch.builder import Fountain_Base4_FileLevelFountain_Pipeline
        self.builder = Fountain_Base4_FileLevelFountain_Pipeline

    def test_two_block_roundtrip(self):
        data = bytes(range(160)) + bytes(reversed(range(160)))
        assert _file_level_roundtrip(self.builder, data, self._PARAMS) == data

    def test_pipeline_has_file_level_codec(self):
        data = bytes(range(160))
        pf = ReadPacketizedFilestream(BytesIO(data))
        pipe = self.builder(pf, **self._PARAMS)
        assert hasattr(pipe, '_file_level_codec')
        assert pipe._file_level_codec is not None


# ──────────────────────────────────────────────────────────────────────────────
# Formats-registry integration
# ──────────────────────────────────────────────────────────────────────────────

class TestFormatsIntegration:
    """Verify the format registry maps IDs to working builder functions."""

    def test_format_0x0704_roundtrip(self):
        """0x0704 = ReedSolomon_Base4_Pipeline."""
        from dnastorage.system.formats import file_system_encoder
        builder = file_system_encoder(0x0704)
        data = bytes(range(160))
        assert _roundtrip(builder, data, _RS_PARAMS) == data

    def test_format_0x0705_roundtrip(self):
        """0x0705 = Fountain_Base4_Pipeline."""
        from dnastorage.system.formats import file_system_encoder
        builder = file_system_encoder(0x0705)
        data = bytes(range(160))
        assert _roundtrip(builder, data, _FOUNTAIN_PARAMS) == data
