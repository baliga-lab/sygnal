#!/usr/bin/env python3

import sys
import unittest
import xmlrunner

import sygnal

class StatsTest(unittest.TestCase):

    def test_correlation(self):
        a1 = [1.0, 2.0, 1.5]
        a2 = [3.0, 1.0, 1.2]
        coeff, pval = sygnal.correlation(a1, a2)
        self.assertAlmostEquals(coeff, -0.907841299003)
        self.assertAlmostEquals(pval, 0.275458676635)

    def test_ttest(self):
        a1 = [1.0, 2.0, 1.5]
        a2 = [3.0, 1.0, 1.2]
        tstat, pval = sygnal.ttest(a1, a2)
        self.assertAlmostEquals(tstat, -0.33409177193234346)
        self.assertAlmostEquals(pval, 0.75509203857458196)

class CompareMiRNANamesTest(unittest.TestCase):

    def test_compare_identical(self):
        self.assertTrue(sygnal.compareMiRNANames('miRNA', 'miRNA'))

    def test_compare_b_is3p(self):
        self.assertTrue(sygnal.compareMiRNANames('hsa-let-7f-3p', 'hsa-let-7f-3p-1-3p'))
        self.assertTrue(sygnal.compareMiRNANames('hsa-let-7f-3p', 'hsa-let-7f-3pb-1-3p'))

    def test_compare_b_is3p_mismatch(self):
        self.assertFalse(sygnal.compareMiRNANames('hsa-let-7f-3p', 'hsa-let-7f-3pp-1-3p'))

    def test_compare_b_is5p(self):
        self.assertTrue(sygnal.compareMiRNANames('hsa-let-7f', 'hsa-let-7f-2-5p'))
        self.assertTrue(sygnal.compareMiRNANames('hsa-let-7f', 'hsa-let-7fa-2-5p'))

    def test_compare_b_is5p_mismatch(self):
        self.assertFalse(sygnal.compareMiRNANames('hsa-let-7f', 'hsa-let-7fp-2-5p'))

    def test_compare_a_is3p(self):
        self.assertTrue(sygnal.compareMiRNANames('hsa-let-7f-3p-1-3p', 'hsa-let-7f-3p'))
        self.assertTrue(sygnal.compareMiRNANames('hsa-let-7f-3pb-1-3p', 'hsa-let-7f-3p'))

    def test_compare_a_is3p_mismatch(self):
        self.assertFalse(sygnal.compareMiRNANames('hsa-let-7f-3pp-1-3p', 'hsa-let-7f-3p'))

    def test_compare_a_is5p(self):
        self.assertTrue(sygnal.compareMiRNANames('hsa-let-7f-2-5p', 'hsa-let-7f'))
        self.assertTrue(sygnal.compareMiRNANames('hsa-let-7fa-2-5p', 'hsa-let-7f'))

    def test_compare_a_is5p_mismatch(self):
        self.assertFalse(sygnal.compareMiRNANames('hsa-let-7fp-2-5p', 'hsa-let-7f'))


if __name__ == '__main__':
    suite = [unittest.TestLoader().loadTestsFromTestCase(StatsTest),
             unittest.TestLoader().loadTestsFromTestCase(CompareMiRNANamesTest)]
    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(suite))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suite))
