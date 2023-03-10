# -*- coding: utf-8 -*-

from .context import dnastorage

import unittest

from dnastorage.primer.design import *

class AdvancedTestSuite(unittest.TestCase):
    """Advanced test cases."""

    def test_design_rules_unmet(self):
        s = "AAAAAAAAAAAAAAAAAAAG"
        dr = build_standard_design_rules([],False)
        assert dr.check(s) == False

    def test_design_rules_met(self):
        s = "TACCAACATTGCCGCAACTG"
        dr = build_standard_design_rules([],True)
        assert dr.check(s) == True



if __name__ == '__main__':
    unittest.main()
