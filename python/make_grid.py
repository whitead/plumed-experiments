from math import ceil, floor, log, exp
import numpy as np
from scipy.integrate import simps
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
        self.min = []
        self.max = []
        self.min = []
        self.types = []
        self.periodic = []
        self.pot = None


    @property
    def dims(self):
        return len(self.min)
    
    @property
    def nbins(self):
        if(self.pot is None):
            return ()
        return np.shape(self.pot)
    
    @property
    def ncv(self):
        if(self.pot is None):
            return 0
        return len(np.shape(self.pot))

    @property
    def dx(self):
        if(self.pot is None):
            return 0
        return [float(max - min) / nb for max,min,nb in zip(self.max, self.min, self.nbins)]

    def add_cv(self, name, min, max, bin_number, periodic=False):
        self.min.append(min)
        self.max.append(max)
        self.periodic.append(periodic)
        self.types.append(Grid.cv_type_map[name])
        if(self.pot is None):
            self.pot = np.zeros(bin_number)
        else:
            self.pot = np.resize(self.pot, self.nbins + (int(bin_number),))

    def to_index(self, x, i):
        return max(0, min(self.nbins[i] - 1, int(floor( (x - self.min[i]) / self.dx[i]) )))
               
    def add_value(self, x, v):
        if(len(x) != self.ncv):
            raise ValueError("Dimension of given x vector does not match grid dimension!")
        index = [0 for xi in x]
        for i, xi in enumerate(x):
            assert xi >= self.min[i] and xi <= self.max[i],"Mesh point is not within grid dimension {}: {}, [{}, {}]".format(i, xi, self.min[i], self.max[i])
            index[i] = self.to_index(xi, i)
        self.pot[tuple(index)] += v

    def set_value(self, x, v):
        if(len(x) != self.ncv):
            raise ValueError("Dimension of given x vector does not match grid dimension!")
        index = [0 for xi in x]
        for i, xi in enumerate(x):
            assert xi >= self.min[i] and xi <= self.max[i],"Mesh point is not within grid dimension {}: {}, [{}, {}]".format(i, xi, self.min[i], self.max[i])
            index[i] = self.to_index(xi, i)
        self.pot[tuple(index)] = v

            
    def get_value(self, x, v):
        index = [0 for xi in x]
        for i, xi in enumerate(x):
            assert xi >= self.min[i] and xi <= self.max[i],"Mesh point is not within grid dimension {}: {}, [{}, {}]".format(i, xi, self.min[i], self.max[i])
            index[i] = self.to_index(xi, i)
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

    def resize(self, new_shape):
        from scipy.ndimage.interpolation import zoom
        zoom_factor = np.array(new_shape, dtype='float') / self.nbins
        self.pot = zoom(self.pot, zoom_factor, prefilter=True, mode='nearest')
                        
    def _enumerate_grid(self, fxn, dim=None, indices=[], end_fxn=None):
        if(dim is None):
            dim = self.ncv - 1
        if(dim > 0):
            for i in range(self.nbins[dim]):
                self._enumerate_grid(fxn, 
                                     dim - 1, 
                                     Grid._prepend_emit(indices, i), end_fxn)
            if(end_fxn is not None):
                self._enumerate_grid(fxn, 
                                     dim - 1, 
                                     Grid._prepend_emit(indices, self.nbins[dim]), end_fxn)

                
        else:
            #check if we are at an end
            if(end_fxn is not None):
                #check if any index is at an end
                if(reduce(lambda x,y: x or y, [x == y for x,y in zip(indices,self.nbins)])):
                    for i in range(self.nbins[dim] + 1):
                        end_fxn(Grid._prepend_emit(indices, i))
                else:
                    for i in range(self.nbins[dim]):
                        fxn(Grid._prepend_emit(indices, i))
                    end_fxn(Grid._prepend_emit(indices, self.nbins[dim]))
            else:
                for i in range(self.nbins[dim]):
                    fxn(Grid._prepend_emit(indices, i))



    def _print_grid(self, indices, output):
        for i,j in enumerate(indices):
            output.write('{:8} '.format(j * self.dx[i] + self.min[i]))        
        output.write('{:08}\n'.format(self.pot[tuple(indices)]))

    def _print_grid_end(self, indices, output):
        for i,j in enumerate(indices):
            output.write('{:8} '.format(j * self.dx[i] + self.min[i]))
        #copy the last bin to the boundary
        indices = [i if i < nb else nb - 1 for i,nb in zip(indices,self.nbins)]
        output.write('{:08}\n'.format(self.pot[tuple(indices)]))


    def normalize(self):
        interval = 250
        grids = [np.arange(min, max, dx) for min,max,dx in zip(self.min, self.max, self.dx)]
        mesh = np.meshgrid(*grids)
        Z = np.exp(self.pot)
        grids.reverse()
        for g in grids:
            Z = simps(Z, g)
        self.pot -= np.log(Z)

    

    def write(self, output):
        output.write('#! FORCE 0\n')
        output.write('#! NVAR {}\n'.format(self.ncv))
        self._print_header_array('TYPE', self.types, output)

        #Some kind of weird Plumed convention
        mod_bins = list(self.nbins)        
        mod_max = copy.copy(self.max)
        for i,p in enumerate(self.periodic):
            mod_bins[i] -= 1 if p else 0
            mod_max[i] -= self.dx[i] if p else 0
        
        self._print_header_array('BIN', np.array(mod_bins), output)
        self._print_header_array('MIN', self.min, output)
        self._print_header_array('MAX', mod_max, output)
        self._print_header_array('PBC', [1 if x else 0 for x in self.periodic], output)
        self._enumerate_grid(lambda x: self._print_grid(x, output), end_fxn=lambda x: self._print_grid_end(x,output))
        

    def add_png_to_grid(self, filename):
        if(self.ncv != 2):
            raise ValueError("This method only makes sense on 2D grids. Grid is currently {} dimension".format(self.ncv))
        from pylab import imread, imshow, gray, mean
        a = imread(filename) # read to RGB file
        if(np.shape(a)[2] == 4):
            gray_scale = mean(a[:,:,0:2],2) * a[:,:,3] # convert to grayscale by multiplication with alpha channel
        else:
            gray_scale = mean(a,2)# convert to gray scale with meana

            
        gray_scale = gray_scale.astype(np.float64)
        too_small = exp(-10)
        gray_scale[np.where(gray_scale < too_small)] = too_small
        self.resize(np.shape(gray_scale))
        self.pot += np.log(gray_scale)
        self.normalize()

def test():
    import sys
    import matplotlib.pyplot as plt
    g = Grid()
    g.add_cv("Distance", 0, 12, 16)
    g.add_cv("Distance", 0, 12, 16)
    g.add_png_to_grid("circle.png")
    g.resize((256,256))
    plt.imshow(np.exp(g.pot), interpolation='none', cmap='gray')
    plt.savefig("circle_out.png")
    g.write(sys.stdout)    

    g = Grid()
    g.add_cv("Distance", 0, 10, 128, True)
    g.add_cv("Distance", 0, 10, 128, True)
    g.add_png_to_grid("uc_shield.png")
    g.resize((128,128))
    plt.imshow(np.exp(g.pot), interpolation='none', cmap='gray')
    plt.savefig("uc_out.png")



if __name__ == "__main__":
    test()
