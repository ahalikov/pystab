
import unittest

from pystab.mechanics import solve_slae_by_gauss
from sympy.core.symbol import symbols

class  SolveSlaeByGaussTestCase(unittest.TestCase):
    def test_1(self):
        x,y,z = symbols('x y z')
        eqns = [2*x + y - z - 8, -3*x - y + 2*z + 11, -2*x + y + 2*z + 3]
        vars = [x, y, z]
        self.assertEqual(solve_slae_by_gauss(eqns, vars), {x: 2, y: 3, z: -2});

if __name__ == '__main__':
    unittest.main()

