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


if __name__ == '__main__':
    suite = [unittest.TestLoader().loadTestsFromTestCase(StatsTest)]
    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(suite))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suite))
