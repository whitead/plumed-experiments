from math import ceil, floor, log
import numpy as np
import copy

class Grid:
    """PLUMED grid class
    
    Member variables: 
      nbins: number of bins for each cv
      min: min for each cv
      max: max for each cv
      periodic: logical, indicating periodicity for each cv
      ncv: number of cv
      dx: grid spacing for each cv
      pot: array for potential      
      types: type for each cv
    """
    cv_type_map = {"Distance":1,
                    "Minimum distance":2,
                    "Coordination number":3,
                    "Angle":4,
                    "Torsion":5,
                    "Alpha-beta similarity":6,
                    "Hydrogen bonds":7,
                    "Dipole":8,
                    "Radius of gyration":11,
                    "Dihedral correlation":16,
                    "Interfacial water":20,
                    "Path collective variable S":30,
                    "Path collective variable Z":31,
                    "Absolute position":32,
                    "Electrostatic potential":33,
                    "Puckering coordinates":34,
                    "Energy":35,
                    "Helix loops":36,
                    "Alpha helix rmsd":37,
                    "Antiparallel beta rmsd":38,
                    "Parallel beta rmsd":39,
                    "PCA projection":42,
                    "Contact Map":45,
                    "SPRINT":55}

    def __init__(self):
        """Create a grid. Add individual CVs using the add_cv method"""
        self.nbins = []
        self.min = []
        self.max = []
        self.min = []
        self.dx = []
        self.types = []
        self.nbins = []     
        self.periodic = []
        self.pot = None
        self.ncv = 0

    def add_cv(self, name, dx, min, max, periodic=False):
        if(self.pot is not None):
            raise ValueError("You must add all CVs before adding values")
        self.nbins.append(int(ceil((max - min) / dx)))
        self.min.append(min)
        self.max.append(max)
        self.periodic.append(periodic)
        self.types.append(Grid.cv_type_map[name])
        self.dx.append(float(dx))
        self.ncv += 1
               
    def add_value(self, x, v):
        if(self.pot is None):
            self.pot = np.zeros(self.nbins)
        if(len(x) != self.ncv):
            raise ValueError("Dimension of given x vector does not match grid dimension!")
        index = [0 for xi in x]
        for i, xi in enumerate(x):
            assert xi >= self.min[i] and xi <= self.max[i],"Mesh point is not within grid dimension {}: {}, [{}, {}]".format(i, xi, self.min[i], self.max[i])
            index[i] = max(0, min(self.nbins[i] - 1, int(floor( (xi - self.min[i]) / self.dx[i]) )))
        self.pot[tuple(index)] += v

    def set_value(self, x, v):
        if(self.pot is None):
            self.pot = np.zeros(self.nbins)
        if(len(x) != self.ncv):
            raise ValueError("Dimension of given x vector does not match grid dimension!")
        index = [0 for xi in x]
        for i, xi in enumerate(x):
            assert xi >= self.min[i] and xi <= self.max[i],"Mesh point is not within grid dimension {}: {}, [{}, {}]".format(i, xi, self.min[i], self.max[i])
            index[i] = max(0, min(self.nbins[i] - 1, int(floor( (xi - self.min[i]) / self.dx[i]) )))
        self.pot[tuple(index)] = v

            
    def get_value(self, x, v):
        if(self.pot is None):
            self.pot = np.zeros(self.nbins)
        index = [0 for xi in x]
        for i, xi in enumerate(x):
            assert xi >= self.min[i] and xi <= self.max[i],"Mesh point is not within grid dimension {}: {}, [{}, {}]".format(i, xi, self.min[i], self.max[i])
            index[i] = max(0, min(self.nbins[i] - 1, int(floor( (xi - self.min[i]) / self.dx[i]) )))
        return self.pot[tuple(index)]


    def _print_header_array(self, name, array, output):
        output.write('#! {} '.format(name))
        for a in array:
            output.write('{} '.format(a))
        output.write('\n')

    @staticmethod
    def _prepend_emit(array, element):
        array_copy = copy.copy(array)
        array_copy.insert(0, element)
        return array_copy
        
    def _enumerate_grid(self, fxn, dim=None, indices=[]):
        if(dim is None):
            dim = self.ncv - 1
        if(dim > 0):
            for i in range(self.nbins[dim]):
                self._enumerate_grid(fxn, 
                               dim - 1, 
                               Grid._prepend_emit(indices, i))
                
        else:
            for i in range(self.nbins[dim]):            
                fxn(Grid._prepend_emit(indices, i))

    def _print_grid(self, indices, output):
        for i,j in enumerate(indices):
            output.write('{0:05} '.format(j * self.dx[i] + self.min[i]))        
        output.write('{0:05}\n'.format(self.pot[tuple(indices)]))
            

    def write(self, output):
        output.write('#! FORCE 0\n')
        output.write('#! NVAR {}\n'.format(self.ncv))
        self._print_header_array('TYPE', self.types, output)

        #Some kind of weird Plumed convention
        mod_bins = copy.copy(self.nbins)        
        mod_max = copy.copy(self.max)
        for i,p in enumerate(self.periodic):
            mod_bins[i] -= 1 if p else 0
            mod_max[i] -= dx[i] if p else 0                

        self._print_header_array('BIN', mod_bins, output)
        self._print_header_array('MIN', self.min, output)
        self._print_header_array('MAX', mod_max, output)
        self._print_header_array('PBC', [0 if x else 1 for x in self.periodic], output)
        self._enumerate_grid(lambda x: self._print_grid(x, output))
        

    def add_png_to_grid(self, filename):
        if(self.ncv != 2):
            raise ValueError("This method only makes sense on 2D grids")
        from pylab import imread, imshow, gray, mean
        a = imread(filename) # read to RGB file
        gray_scale = mean(a,2) # convert to grayscale
        csum = np.sum(gray_scale)        
        if(np.shape(gray_scale)[0] != self.nbins[0] or np.shape(gray_scale)[1] != self.nbins[1]):
            raise ValueError("Your image must be exactly the same number of pixels as your bin number: {}".format(self.nbins))

        if(self.pot is None):
            self.pot = np.zeros(self.nbins)
        for i in range(np.shape(gray_scale)[0]):
            for j in range(np.shape(gray_scale)[1]):
                self.pot[i,j] += log(gray_scale[i,j]) - log(csum)
                                                                

def test():
    import sys
    g = Grid()
    g.add_cv("Distance", 0.5, 0, 16, True)
    g.add_cv("Distance", 0.5, 0, 16, True)
    g.add_png_to_grid("circle.png")
    g.write(sys.stdout)

if __name__ == "__main__":
    test()
