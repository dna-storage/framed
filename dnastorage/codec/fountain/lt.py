"""
lt.py - LT (Luby Transform) fountain code implementation.

This module provides a self-contained systematic LT code suitable for use as an
outer erasure code in the DNA storage pipeline.

Design decisions
----------------
* **Systematic**: output symbols 0..k-1 are the original source symbols unchanged.
  Only symbols k..k+r-1 are XOR-combinations, matching the BaseOuterCodec assumption
  that parity strands come after data strands.
* **Robust Soliton distribution** (Shokrollahi 2006) for good overhead/recovery tradeoffs.
* **Deterministic PRNG**: a seeded `random.Random` instance reproduces the exact same
  Tanner graph at both encoder and decoder — no side-channel needed.
* **Belief-propagation peeling decoder**: O(k log k) average complexity, well-suited for
  typical DNA block sizes (50–2000 strands).
* Symbols are byte-arrays (`bytearray`).  XOR is applied element-wise.
"""

import random
import math
from typing import List, Optional, Tuple


def _robust_soliton_distribution(k: int, c: float = 0.1, delta: float = 0.05) -> List[float]:
    """
    Compute the Robust Soliton degree-probability distribution over [1..k].

    Parameters
    ----------
    k : number of source symbols
    c : tuning constant (0 < c < 1); larger values increase overhead but improve recovery
    delta : desired failure probability (0 < delta < 1)

    Returns a list P where P[d-1] is the probability of degree d (1-indexed list, length k).
    """
    # Ideal Soliton
    rho = [0.0] * k
    rho[0] = 1.0 / k          # d=1
    for d in range(2, k + 1):
        rho[d - 1] = 1.0 / (d * (d - 1))

    # Spike parameter R
    R = c * math.log(k / delta) * math.sqrt(k)
    R = max(1.0, R)  # guard against degenerate small k

    # Robust extra term tau
    tau = [0.0] * k
    thresh = int(math.floor(k / R))
    for d in range(1, thresh):          # d = 1 .. floor(k/R)-1
        tau[d - 1] = (R / k) * (1.0 / d)
    if thresh >= 1 and thresh <= k:
        tau[thresh - 1] = (R / k) * math.log(R / delta)

    # Combine and normalise
    mu = [rho[i] + tau[i] for i in range(k)]
    total = sum(mu)
    mu = [v / total for v in mu]
    return mu


def _cumulative(dist: List[float]) -> List[float]:
    """Return cumulative distribution list for fast sampling."""
    cdf = []
    running = 0.0
    for p in dist:
        running += p
        cdf.append(running)
    return cdf


def _sample_degree(cdf: List[float], rng: random.Random) -> int:
    """Return a degree (1-indexed integer) sampled from the given CDF."""
    u = rng.random()
    # bisect manually to avoid extra imports
    lo, hi = 0, len(cdf)
    while lo < hi:
        mid = (lo + hi) // 2
        if cdf[mid] < u:
            lo = mid + 1
        else:
            hi = mid
    return lo + 1  # 1-indexed degree


def _symbol_neighbors(parity_index: int, k: int, degree: int, rng: random.Random) -> List[int]:
    """Sample `degree` distinct source-symbol indices in [0, k) without replacement."""
    return rng.sample(range(k), degree)


class LTCode:
    """
    Systematic LT fountain code.

    Parameters
    ----------
    k : number of source symbols
    seed : integer seed for the PRNG (must be the same at encoder and decoder)
    c, delta : robust soliton tuning parameters
    """

    def __init__(self, k: int, seed: int = 42, c: float = 0.1, delta: float = 0.05) -> None:
        if k < 1:
            raise ValueError("k must be >= 1")
        self.k = k
        self.seed = seed
        self.c = c
        self.delta = delta
        self._dist = _robust_soliton_distribution(k, c, delta)
        self._cdf = _cumulative(self._dist)

    def _make_rng(self, parity_index: int) -> random.Random:
        """Create a fresh RNG seeded deterministically for parity symbol `parity_index`.

        The multiplicative constant 0x9E3779B9 is the 32-bit truncation of the
        golden-ratio hash multiplier (2^32 / phi).  Multiplying the parity index by
        it before XORing with the base seed ensures that adjacent parity indices
        produce well-dispersed seeds (Fibonacci / Knuth multiplicative hashing).
        """
        return random.Random(self.seed ^ (parity_index * 0x9E3779B9 & 0xFFFFFFFF))

    def neighbors(self, parity_index: int) -> List[int]:
        """Return the list of source indices XORed to produce parity symbol `parity_index`."""
        rng = self._make_rng(parity_index)
        degree = _sample_degree(self._cdf, rng)
        degree = min(degree, self.k)
        return _symbol_neighbors(parity_index, self.k, degree, rng)

    # ------------------------------------------------------------------
    # Encoding
    # ------------------------------------------------------------------

    def encode(
        self,
        source_symbols: List[bytearray],
        num_parity: int,
    ) -> List[Tuple[int, bytearray, List[int]]]:
        """
        Encode *only* the parity symbols (systematic part is implicit — callers
        keep the original source symbols unchanged).

        Parameters
        ----------
        source_symbols : list of k bytearray objects (all the same length)
        num_parity : number of additional parity symbols to produce

        Returns
        -------
        List of (parity_index, encoded_bytearray, [source_indices]) tuples,
        one per parity symbol.  `parity_index` is 0-based and refers to the
        position *after* the k source symbols, i.e. total position = k + parity_index.
        """
        if len(source_symbols) != self.k:
            raise ValueError(
                f"Expected {self.k} source symbols, got {len(source_symbols)}"
            )
        sym_len = len(source_symbols[0])
        result: List[Tuple[int, bytearray, List[int]]] = []

        for pi in range(num_parity):
            nbrs = self.neighbors(pi)
            encoded = bytearray(sym_len)
            for idx in nbrs:
                s = source_symbols[idx]
                for byte_pos in range(sym_len):
                    encoded[byte_pos] ^= s[byte_pos]
            result.append((pi, encoded, nbrs))

        return result

    # ------------------------------------------------------------------
    # Decoding
    # ------------------------------------------------------------------

    def decode(
        self,
        received: List[Tuple[int, Optional[bytearray], List[int]]],
        sym_len: int,
    ) -> Optional[List[Optional[bytearray]]]:
        """
        Recover source symbols using a belief-propagation (peeling) decoder.

        Parameters
        ----------
        received : list of (symbol_index, bytearray_or_None, [source_indices]) tuples.
            * symbol_index in [0, k+r): 0..k-1 are source, k..k+r-1 are parity.
            * bytearray_or_None is the received byte content; None means erased.
            * source_indices is the neighbor list (empty [] for source symbols).
        sym_len : length of each symbol in bytes

        Returns
        -------
        List of k bytearray objects (source symbols), with None entries where
        recovery failed, or None if the decoder could not start.
        """
        k = self.k

        # Initialise source-symbol slots
        decoded: List[Optional[bytearray]] = [None] * k

        # Populate known source symbols from received list
        for sym_idx, data, nbrs in received:
            if sym_idx < k and data is not None:
                decoded[sym_idx] = bytearray(data)

        # Build adjacency: for each parity node, track remaining unknown neighbors
        # parity_nodes: list of [remaining_data, remaining_unknown_nbrs]
        parity_nodes: List[Optional[List]] = []
        for sym_idx, data, nbrs in received:
            if sym_idx >= k and data is not None:
                remaining_data = bytearray(data)
                # Remove already-known source symbols
                unknown_nbrs = []
                for nbr in nbrs:
                    if decoded[nbr] is not None:
                        # XOR out the known symbol
                        for bp in range(sym_len):
                            remaining_data[bp] ^= decoded[nbr][bp]
                    else:
                        unknown_nbrs.append(nbr)
                parity_nodes.append([remaining_data, unknown_nbrs])

        # Peeling loop
        changed = True
        while changed:
            changed = False
            for node in parity_nodes:
                if node is None:
                    continue
                remaining_data, unknown_nbrs = node

                # Remove any neighbors that became known since last iteration
                still_unknown = []
                for nbr in unknown_nbrs:
                    if decoded[nbr] is not None:
                        for bp in range(sym_len):
                            remaining_data[bp] ^= decoded[nbr][bp]
                        changed = True
                    else:
                        still_unknown.append(nbr)
                node[1] = still_unknown

                # Degree-1 check: sole unknown neighbor can be recovered
                if len(still_unknown) == 1:
                    target = still_unknown[0]
                    if decoded[target] is None:
                        decoded[target] = bytearray(remaining_data)
                        changed = True
                        node[1] = []  # this parity node is consumed

        return decoded
