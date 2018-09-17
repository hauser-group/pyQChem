import unittest
from os.path import dirname, join
import pyqchem as qc
import numpy as _np

class TestOutputClassesForce(unittest.TestCase):

    def setUp(self):
        self.path = join(dirname(__file__), "files")
        self._expected_grad_c2h6 = 0.01302
        self._expected_grad_vec_c2h6 = _np.array(
            [[ 1.30202e-02, -1.36000e-05, -9.20000e-06],
            [-1.30202e-02,  1.36000e-05,  9.20000e-06],
            [-4.92410e-03,  6.41380e-03, -9.36410e-03],
            [-4.92560e-03,  4.91470e-03,  1.02363e-02],
            [-4.91750e-03, -1.13134e-02, -8.64400e-04],
            [ 4.91750e-03,  1.13134e-02,  8.64400e-04],
            [ 4.92560e-03, -4.91470e-03, -1.02363e-02],
            [ 4.92410e-03, -6.41380e-03,  9.36410e-03]]).T

        self._expected_grad_h2 = 2.800E+03
        self._expected_grad_vec_h2 = _np.array(
            [[    0.       ,    -0.       ],
            [    0.       ,    -0.       ],
            [ 2800.1617645, -2800.1617645]])

    def test_ethane(self):
        """ Simple TestCase to ensure correct handling of the
        somewhat strange printing format of QChem gradients"""

        # Test for QChem version 5.0
        job = qc.read(join(self.path, "force_c2h6_5.0.out"))
        self.assertEqual(job.force.gradient,
            self._expected_grad_c2h6)
        _np.testing.assert_array_equal(job.force.gradient_vector,
            self._expected_grad_vec_c2h6)

    def test_h2(self):
        """ Hard TestCase as the large numerical values in
        the gradient are more difficult to parse. """

        # Test for QChem version 5.0
        job = qc.read(join(self.path, "force_h2_5.0.out"))
        self.assertEqual(job.force.gradient,
            self._expected_grad_h2)
        _np.testing.assert_array_equal(job.force.gradient_vector,
            self._expected_grad_vec_h2)

if __name__ == '__main__':
    unittest.main()
