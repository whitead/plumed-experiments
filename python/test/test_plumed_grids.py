import unittest
import sys
import numpy as np
import textwrap

sys.path.append('..')

EPSILON = 0.000000001

from plumed_grids import *

class TestPlumedGrid(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_init(self):
        g = Grid()
        
    def test_add_cv(self):
        g = Grid()
        g.add_cv("Coordination number", 0, 3, 5)
        g.add_cv("Coordination number", 0, 3, 5)
        self.assertEqual(np.shape(g.pot), (5,5))
        self.assertEqual(g.nbins, (5,5))

    def test_clone(self):
        g = Grid()
        g.add_cv("Coordination number", 0, 3, 5)
        g.add_cv("Coordination number", 0, 3, 5)
        h = g.clone()
        self.assertFalse(h is g)
        self.assertEqual(h.__str__(), g.__str__())
        self.assertTrue(np.all(h.pot[:] == g.pot[:]))

    def test_to_index(self):
        g = Grid()
        g.add_cv("Distance", 0, 100, 100)
        self.assertEqual(g.to_index(0,0) , 0)
        self.assertEqual(g.to_index(49.999,0) , 49)
        self.assertEqual(g.to_index(99,0) , 99)
        self.assertEqual(g.to_index(101,0) , 99)

        g.add_cv("Distance", 0, 10, 10, True)
        self.assertEqual(g.to_index(-0.1,1) , 9)
        self.assertEqual(g.to_index(-1,1) , 9)
        self.assertEqual(g.to_index(-1.01,1) , 8)
        self.assertEqual(g.to_index(10,1) , 0)
        self.assertEqual(g.to_index(11,1) , 1)

    def test_np_to_index(self):

        #this one was failing previously
        g = Grid()
        g.add_cv("Distance", 0, 2, 2)
        g.add_cv("Distance", 0, 2, 2)
        self.assertTrue(np.all(g.np_to_index([0,1]) == (0,1)))


        #test w/o periodicity
        g = Grid()
        g.add_cv("Distance", 0, 100, 100)
        g.add_cv("Distance", 0, 100, 100)
        g.add_cv("Distance", 0, 100, 100)

        self.assertTrue(np.all(g.np_to_index([0,1,2]) == (0,1,2)))

        #compare both methods
        self.assertEqual(g.np_to_index([92,21.3,27.2])[1], g.to_index(21.3, 1))

        #with periodidicity
        g = Grid()
        g.add_cv("Distance", 0, 10, 10)
        g.add_cv("Distance", 0, 10, 10, True)
        g.add_cv("Distance", 0, 10, 10)

        self.assertEqual(g.np_to_index([0,1,2]) , (0,1,2))
        self.assertEqual(g.np_to_index([0,-0.5,2]) , (0,9,2))



                

    def test_set_bin_number(self):
        g = Grid()
        g.add_cv("Distance", 0, 10, 10)        
        g.set_value([0.5], 1)
        g.set_bin_number(100)
        self.assertTrue(g.pot[0]-1 < EPSILON)
        self.assertTrue(np.abs(g.pot[9]) < 0.2)

    def test_load_data(self):
        g = Grid()
        g.add_cv("Distance", 0, 2, 2)
        g.add_cv("Distance", 0, 2, 2)
        with open('test', 'w') as f:
            f.write(textwrap.dedent('''
            0 0 1
            0 1 2
            1 0 3
            1 1 4'''))            
        g.load_data('test')
        self.assertEqual(g.pot[0,0],1)
        self.assertEqual(g.pot[1,0],3)
        self.assertEqual(g.pot[1,1] ,4)
        self.assertEqual(g.get_value( [1.5,1.5] ) ,4)


    def test_load_data_expand(self):
        g = Grid()
        g.add_cv("Distance", 0, 1, 1)
        g.add_cv("Distance", 0, 1, 1)
        with open('test', 'w') as f:
            f.write(textwrap.dedent('''
            0 0 2
            0 1 4
            1 0 4
            1 1 2'''))            
        g.load_data('test')
        self.assertEqual(g.nbins , (1,1))
        print g.pot[0,0]
        self.assertTrue(abs(g.pot[0,0] - 4) < EPSILON)

    def test_load_plumed_1d(self):
        
    

if __name__ == '__main__':
    unittest.main()
