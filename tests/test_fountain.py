# -*- coding: utf-8 -*-
"""
tests/test_fountain.py

Unit and integration tests for the LT fountain code outer codec.

Run with:  python -m pytest tests/test_fountain.py -v
"""

import sys
import os
import random
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from dnastorage.codec.fountain.lt import LTCode
from dnastorage.codec.block import FountainOuterPipeline
from dnastorage.strand_representation import BaseDNA
from dnastorage.codec.base_conversion import convertIntToBytes, convertBytesToInt


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_strands(num_packets, strands_per_packet, strand_len=8, seed=0):
    """Return a list-of-lists of BaseDNA with random codewords."""
    rng = random.Random(seed)
    packets = []
    for sp in range(num_packets):
        pkt = []
        for si in range(strands_per_packet):
            cws = [rng.randint(0, 255) for _ in range(strand_len)]
            strand = BaseDNA(codewords=cws, index_ints=(0, sp, si))
            pkt.append(strand)
        packets.append(pkt)
    return packets


def erase_packets(packets, erase_indices):
    """Set all codewords in each listed packet index to None (simulate dropout)."""
    for i in erase_indices:
        for strand in packets[i]:
            strand.codewords = [None] * len(strand.codewords)
    return packets


# ---------------------------------------------------------------------------
# LTCode unit tests
# ---------------------------------------------------------------------------

class TestLTCode(unittest.TestCase):

    def _make_symbols(self, k, sym_len=10, seed=7):
        rng = random.Random(seed)
        return [bytearray(rng.randint(0, 255) for _ in range(sym_len)) for _ in range(k)]

    def test_roundtrip_no_erasures(self):
        """Encode + decode with no erasures should recover all source symbols."""
        k = 10
        num_parity = 15
        sym_len = 12
        source = self._make_symbols(k, sym_len)
        lt = LTCode(k, seed=99)

        parity_syms = lt.encode(source, num_parity)
        self.assertEqual(len(parity_syms), num_parity)

        # Build received list: all source + all parity present
        received = [(i, bytearray(source[i]), []) for i in range(k)]
        for pi, enc_ba, nbrs in parity_syms:
            received.append((k + pi, enc_ba, nbrs))

        decoded = lt.decode(received, sym_len)
        self.assertIsNotNone(decoded)
        for i in range(k):
            self.assertEqual(decoded[i], source[i], f"Mismatch at symbol {i}")

    def test_roundtrip_with_erasures(self):
        """Recover erased source symbols using parity."""
        k = 20
        num_parity = 30
        sym_len = 15
        source = self._make_symbols(k, sym_len, seed=42)
        lt = LTCode(k, seed=7)

        parity_syms = lt.encode(source, num_parity)

        # Erase roughly a third of source symbols
        erase_set = set(range(0, k, 3))  # indices 0,3,6,...
        received = []
        for i in range(k):
            if i in erase_set:
                received.append((i, None, []))
            else:
                received.append((i, bytearray(source[i]), []))
        for pi, enc_ba, nbrs in parity_syms:
            received.append((k + pi, enc_ba, nbrs))

        decoded = lt.decode(received, sym_len)
        self.assertIsNotNone(decoded)
        recovered = 0
        for i in range(k):
            if decoded[i] is not None:
                self.assertEqual(decoded[i], source[i], f"Wrong value at symbol {i}")
                recovered += 1
        # With 30 parity for 20 sources and only ~7 erasures, expect full recovery
        self.assertEqual(recovered, k, "Expected full recovery with sufficient parity")

    def test_degree_distribution_valid(self):
        """Neighbour lists must contain distinct indices in [0, k)."""
        k = 50
        lt = LTCode(k, seed=123)
        for pi in range(100):
            nbrs = lt.neighbors(pi)
            self.assertGreater(len(nbrs), 0)
            self.assertEqual(len(nbrs), len(set(nbrs)), "Duplicate neighbours")
            for idx in nbrs:
                self.assertGreaterEqual(idx, 0)
                self.assertLess(idx, k)

    def test_deterministic_neighbors(self):
        """Same seed must produce identical neighbour lists."""
        k = 15
        lt1 = LTCode(k, seed=55)
        lt2 = LTCode(k, seed=55)
        for pi in range(20):
            self.assertEqual(lt1.neighbors(pi), lt2.neighbors(pi))

    def test_different_seeds_differ(self):
        """Different seeds should (almost always) produce different neighbour lists."""
        k = 15
        lt1 = LTCode(k, seed=1)
        lt2 = LTCode(k, seed=2)
        diffs = sum(1 for pi in range(20) if lt1.neighbors(pi) != lt2.neighbors(pi))
        self.assertGreater(diffs, 0, "Different seeds unexpectedly produced identical graphs")

    def test_single_source_symbol(self):
        """Edge case: k=1 should still work (degree must be 1)."""
        k = 1
        sym_len = 4
        source = [bytearray([10, 20, 30, 40])]
        lt = LTCode(k, seed=0)
        parity = lt.encode(source, 3)
        received = [(0, bytearray(source[0]), [])]
        for pi, enc_ba, nbrs in parity:
            received.append((k + pi, enc_ba, nbrs))
        decoded = lt.decode(received, sym_len)
        self.assertEqual(decoded[0], source[0])


# ---------------------------------------------------------------------------
# FountainOuterPipeline integration tests
# ---------------------------------------------------------------------------

class TestFountainOuterPipeline(unittest.TestCase):

    def _run_pipeline(self, num_data_pkts, strands_per_pkt, strand_len,
                      parity_pkts, erase_indices, seed=42):
        """
        Full encode→(erase)→decode cycle via FountainOuterPipeline._encode/_decode.
        Returns (original_codewords, decoded_codewords) for data sub-packets only.
        """
        # Build FountainOuterPipeline manually (bypass PipeLine infrastructure)
        fp = FountainOuterPipeline(num_data_pkts, parity_pkts, seed=seed)
        # Manually set attrs that BaseOuterCodec.encode would normally set:
        fp._num_data_sub_packets = num_data_pkts
        fp._total_sub_packets = num_data_pkts + parity_pkts
        fp._level = 1

        packets = make_strands(num_data_pkts, strands_per_pkt, strand_len)

        # Save originals before encoding
        original = [
            [list(strand.codewords) for strand in pkt]
            for pkt in packets
        ]

        # Encode
        all_packets = fp._encode(packets)
        self.assertEqual(len(all_packets), num_data_pkts + parity_pkts)

        # Erase some packets
        all_packets = erase_packets(all_packets, erase_indices)

        # Decode
        decoded_packets = fp._decode(all_packets)

        # Extract codewords from data sub-packets
        decoded_cws = [
            [list(strand.codewords) for strand in decoded_packets[sp]]
            for sp in range(num_data_pkts)
        ]
        return original, decoded_cws

    def test_no_erasure_roundtrip(self):
        """Encode + decode with no erasures: data sub-packets must be unchanged."""
        orig, dec = self._run_pipeline(
            num_data_pkts=8, strands_per_pkt=3, strand_len=10,
            parity_pkts=12, erase_indices=[]
        )
        self.assertEqual(orig, dec)

    def test_single_erasure_recovery(self):
        """Erasing one data sub-packet should be recoverable."""
        orig, dec = self._run_pipeline(
            num_data_pkts=10, strands_per_pkt=4, strand_len=8,
            parity_pkts=15, erase_indices=[3]
        )
        self.assertEqual(orig, dec)

    def test_multiple_erasure_recovery(self):
        """Erase several data sub-packets; expect full recovery with ample parity."""
        orig, dec = self._run_pipeline(
            num_data_pkts=12, strands_per_pkt=3, strand_len=6,
            parity_pkts=18, erase_indices=[0, 5, 11]
        )
        self.assertEqual(orig, dec)

    def test_header_serialisation_roundtrip(self):
        """_encode_header / _decode_header must preserve seed and parity count."""
        fp_enc = FountainOuterPipeline(10, 5, seed=12345)
        fp_enc._num_data_sub_packets = 10
        fp_enc._total_sub_packets = 15
        fp_enc._level = 1
        fp_enc._index_bits = 4

        header_bytes = fp_enc._encode_header()

        fp_dec = FountainOuterPipeline(1, 1, seed=0)  # placeholder values
        fp_dec._num_data_sub_packets = 1
        fp_dec._total_sub_packets = 1
        fp_dec._level = 1
        fp_dec._index_bits = 0
        remaining = fp_dec._decode_header(header_bytes)

        self.assertEqual(fp_dec._seed, fp_enc._seed)
        self.assertEqual(fp_dec._parity_packets, fp_enc._parity_packets)
        self.assertEqual(remaining, [])

    def test_parity_only_erasure(self):
        """Erasing only parity sub-packets should still allow full data recovery."""
        num_data = 8
        parity = 10
        orig, dec = self._run_pipeline(
            num_data_pkts=num_data, strands_per_pkt=2, strand_len=5,
            parity_pkts=parity, erase_indices=list(range(num_data, num_data + parity))
        )
        self.assertEqual(orig, dec)

    def test_large_block(self):
        """Validate the key LT advantage: block sizes exceeding the RS GF(2^8) limit of
        255 symbols.  With 300 source sub-packets and 150 parity sub-packets, erasing
        10 source sub-packets should achieve full recovery with high probability."""
        orig, dec = self._run_pipeline(
            num_data_pkts=300, strands_per_pkt=2, strand_len=4,
            parity_pkts=150, erase_indices=list(range(0, 10)), seed=42
        )
        self.assertEqual(orig, dec)


if __name__ == '__main__':
    unittest.main()
