import sys
from numpy.testing import *
sys.path.insert(0, "../")

from coordinate import *


class TestCoordinate(TestCase):

    def test_lambert(self):
        alat = 36
        alon = -84
        cm = -81
        ps = 30
        pn = 40
        ol = 35
        output = (-269462.2055244348, 114585.91618948989)
        res = gd2lambert(alon, alat, cm, ps, pn, ol)
        assert_almost_equal(res, output, decimal=6)
        assert_almost_equal([alon, alat], lambert2gd(output[0], output[1], cm, ps, pn, ol), decimal=6)
