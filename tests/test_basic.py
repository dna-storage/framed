# -*- coding: utf-8 -*-

from .context import dnastorage
from .context import nupack

from dnastorage.codec.base import *
from dnastorage.codec.binary import *
from dnastorage.codec.huffman import *
from dnastorage.codec.dense import *
from dnastorage.codec.base_conversion import *

import unittest

class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_convert_with_exclusions(self):
        primer = "GTCTCGTGGGCTCGG"
        for i in range(3**5):
            s = encodeWithExclusion(i,10,primer)
            assert i == decodeWithExclusion(s,primer)

    def test_encoding(self):
        for i in range(2**8):
            x = convertFromBase(2,convertBase(2,i,8))
            assert x == i
            x = convertFromBase(3,convertBase(3,i,8))
            assert x == i

    def test_rotated_encoding(self):
        for i in range(2**8):
            assert binary_unrotate_decode(binary_rotate_encode(convertBase(2,i,8)))==convertBase(2,i,8)

    def test_huffman(self):
        s = "12345678"
        assert s == huffman_decode(huffman_encode(s))
        s = 'AGCAGCAGC'
        assert s == rotate_decode(rotate_encode(s))

    def test_dense(self):
        s = "12345678"
        assert s == dense_decode(dense_encode(s))


    def test_reedsolomon(self):
        pass



if __name__ == '__main__':
    unittest.main()
