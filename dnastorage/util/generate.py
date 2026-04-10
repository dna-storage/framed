"""Pure-Python fallback for the ``dnastorage.util.generate`` C++ extension.

The compiled extension (built from ``random_int.cpp``) takes priority over this
module when it is present; Python always prefers a compiled ``.so``/``.pyd`` to a
same-named ``.py`` file in the same package directory.  This file is therefore
only used in environments where the extension has not been compiled, providing
the same public interface so that all callers continue to work without changes.
"""

import random as _random

_rng = _random.Random()


def seed():
    """Return a random seed value (mirrors the C++ extension behaviour).

    The C++ implementation uses ``std::random_device`` to produce a 32-bit
    (4-byte) unsigned integer seed, so we do the same here.
    """
    import os
    raw = os.urandom(4)  # 4 bytes → 32-bit integer, matching the C++ extension
    return int.from_bytes(raw, "little")


def set_seed(value):
    """Seed the internal random-number generator."""
    _rng.seed(value)


def rand():
    """Return a uniform random float in [0, 1]."""
    return _rng.random()


def rand_in_range(lower, upper):
    """Return a uniform random integer in [lower, upper] (inclusive)."""
    return _rng.randint(lower, upper)
