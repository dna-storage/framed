"""
filelevel.py - File-level LT fountain code for cross-block erasure protection.

This module provides :class:`FileLevelFountainCodec`, which applies an LT fountain
code *across* all encoded blocks of a file rather than within a single block.  The
result is a set of "parity blocks" that can recover entire missing data blocks, giving
protection against loss of contiguous regions of the DNA pool.

Relationship to the existing pipeline
--------------------------------------
The existing :class:`~dnastorage.codec.block.FountainOuterPipeline` (and
``ReedSolomonOuterPipeline``) operate *within* a single packet/block: they create
parity *sub-packets* from the strands that make up one block.  This class operates
one level higher: it treats each fully-encoded block as a single opaque byte string
and applies XOR-based LT fountain coding across those byte strings.

How it is integrated
---------------------
* **Encoding** (``WriteDNAFilePipeline.close``): before the pipeline encodes anything,
  the raw block bytes are read from the write buffer, parity block bytes are computed,
  and appended to the buffer.  The pipeline then encodes all blocks (data + parity) in
  the normal streaming fashion.  Parity blocks are assigned packet indices
  ``B, B+1, …, B+r-1`` (where *B* is the number of data blocks).

* **Decoding** (``ReadDNAFilePipeline.__init__``): after ``pipe.final_decode()`` the
  ``WritePacketizedFilestream`` may still have missing data-block keys.  Parity blocks
  (at out-of-range indices *B+* ) were silently ignored by the packetised file but
  captured via a hook set on ``PipeLine``.  The LT peeling decoder reconstructs
  missing data blocks from those parity bytes.

Header serialisation
---------------------
``encode_header`` / ``decode_header`` persist the four parameters needed to reproduce
the Tanner graph and interpret the recovered bytes:

    seed          (4 bytes) – LTCode PRNG seed
    parity_blocks (2 bytes) – number of parity blocks appended during encode
    num_data_blocks (2 bytes) – number of source (data) blocks
    block_size    (2 bytes) – bytes per block (after zero-padding)

These 10 bytes are appended to ``PipeLine.encode_header_data()`` and consumed at the
*end* of ``PipeLine.decode_header_data()``.
"""

from typing import Optional

from dnastorage.codec.fountain.lt import LTCode
from dnastorage.codec.base_conversion import convertIntToBytes, convertBytesToInt

import logging
logger = logging.getLogger('dna.storage.codec.filelevel')
logger.addHandler(logging.NullHandler())


class FileLevelFountainCodec:
    """
    File-level LT fountain code that operates across all blocks in a file.

    Parameters
    ----------
    parity_blocks : int
        Number of parity blocks to generate (appended after the data blocks).
    seed : int
        PRNG seed for the deterministic LT Tanner graph.  Must be identical at
        encoder and decoder; it is stored in the encoded header.
    """

    #: Byte length of the serialised header produced by :meth:`encode_header`.
    HEADER_SIZE = 10  # seed(4) + parity_blocks(2) + num_data_blocks(2) + block_size(2)

    def __init__(self, parity_blocks: int, seed: int = 42) -> None:
        self.parity_blocks = parity_blocks
        self.seed = seed
        # Set during encode (generate_parity_bytes) or restored from header during decode.
        self._num_data_blocks: Optional[int] = None
        self._block_size: Optional[int] = None

    # ------------------------------------------------------------------
    # Header serialisation
    # ------------------------------------------------------------------

    def encode_header(self):
        """Return 10 serialised header bytes.  Must be called *after*
        :meth:`generate_parity_bytes` so that ``_num_data_blocks`` and
        ``_block_size`` are populated."""
        assert self._num_data_blocks is not None, \
            "FileLevelFountainCodec.encode_header called before generate_parity_bytes"
        assert self._block_size is not None, \
            "FileLevelFountainCodec.encode_header called before generate_parity_bytes"
        buf = convertIntToBytes(self.seed, 4)
        buf += convertIntToBytes(self.parity_blocks, 2)
        buf += convertIntToBytes(self._num_data_blocks, 2)
        buf += convertIntToBytes(self._block_size, 2)
        return buf

    def decode_header(self, buff):
        """Parse 10 header bytes from *buff*, update ``self`` in-place, and return
        the remaining bytes."""
        self.seed = convertBytesToInt(buff[0:4])
        self.parity_blocks = convertBytesToInt(buff[4:6])
        self._num_data_blocks = convertBytesToInt(buff[6:8])
        self._block_size = convertBytesToInt(buff[8:10])
        return buff[10:]

    # ------------------------------------------------------------------
    # Encoding
    # ------------------------------------------------------------------

    def generate_parity_bytes(self, blocks):
        """
        Compute ``parity_blocks`` parity byte arrays from *blocks*.

        Parameters
        ----------
        blocks : list of :class:`bytearray`
            One entry per data block; **all must have the same length** (pad the last
            block with zero bytes if necessary before calling).

        Returns
        -------
        list of :class:`bytearray`
            One parity bytearray per parity block, each the same length as *blocks[0]*.
        """
        n = len(blocks)
        if n == 0:
            return []
        block_size = len(blocks[0])
        self._num_data_blocks = n
        self._block_size = block_size

        lt = LTCode(n, seed=self.seed)
        parity_tuples = lt.encode(blocks, self.parity_blocks)
        # parity_tuples: list of (parity_index, bytearray, neighbor_list)
        return [enc_ba for _, enc_ba, _ in parity_tuples]

    # ------------------------------------------------------------------
    # Decoding
    # ------------------------------------------------------------------

    def recover_missing(self, all_block_data, num_data_blocks, block_size):
        """
        Recover missing data blocks using the LT peeling decoder.

        Parameters
        ----------
        all_block_data : dict {int -> bytes-like}
            Mapping from packet index to block byte content.  Data blocks have
            indices ``0 .. num_data_blocks-1``; parity blocks have indices
            ``num_data_blocks .. num_data_blocks + parity_blocks - 1``.  Entries
            may be absent or ``None`` for completely missing blocks.
        num_data_blocks : int
            Number of source (data) blocks.
        block_size : int
            Byte length of each block symbol passed to the LT decoder.

        Returns
        -------
        dict {int -> bytearray}
            Recovered data blocks, keyed by their original packet index.
            Only includes blocks that were *missing* from *all_block_data* and
            could be reconstructed.
        """
        lt = LTCode(num_data_blocks, seed=self.seed)

        received = []

        # Source symbols (data blocks 0 .. num_data_blocks-1)
        for i in range(num_data_blocks):
            data = all_block_data.get(i)
            if data is not None:
                received.append((i, bytearray(data), []))
            else:
                received.append((i, None, []))

        # Parity symbols (indices num_data_blocks .. num_data_blocks+parity_blocks-1)
        for pi in range(self.parity_blocks):
            parity_idx = num_data_blocks + pi
            data = all_block_data.get(parity_idx)
            if data is not None:
                nbrs = lt.neighbors(pi)
                received.append((parity_idx, bytearray(data), nbrs))
            # Missing parity blocks are simply omitted from the received list.

        decoded = lt.decode(received, block_size)

        result = {}
        if decoded is None:
            return result

        for i in range(num_data_blocks):
            if all_block_data.get(i) is None and decoded[i] is not None:
                result[i] = bytearray(decoded[i])

        recovered = len(result)
        if recovered:
            logger.info("FileLevelFountainCodec: recovered %d missing block(s)", recovered)
        return result
