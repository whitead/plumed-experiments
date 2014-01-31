from math import ceil, floor, log, exp
import numpy as np
from scipy.integrate import simps
import copy

class Grid(object):
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
        self._clear()

    def _clear(self):
        self.min = []
        self.max = []
        self.types = []
        self.periodic = []
        self.pot = None
        self.meshgrid = None
        self._nbins = None


    def clone(self):
        g = Grid()
        g.min = copy.copy(self.min)
        g.max = copy.copy(self.max)
        g.types = copy.copy(self.types)
        g.periodic = copy.copy(self.periodic)
        g.pot = np.copy(self.pot)
        return g


    def load_data(self, filename, reset_bounds=True):
        """Read data line by line from a file and load it into the potential"""
        data = np.genfromtxt(filename)
        if(np.shape(data)[1] != self.dims + 1):
            raise ValueError("Incorrect number of dimensions in file {}".fomrat(filename))

        if(reset_bounds):
            for x in data:
                self.min = [min(m,xi) for xi,m in zip(x, self.min)]
                self.max = [max(m,xi) for xi,m in zip(x, self.max)]
            print 'Reset min to {} and max to {}'.format(self.min, self.max) 
            

        #make sure bins are large enough
        old_bins = self.nbins
        self.set_bin_number([np.shape(data)[0] for x in self.nbin])
        for x in data:
            indexs = self.np_to_index(x[:-1])
            self.pot[tuple(indexs)] += x[-1]

        self.set_bin_number(old_bins)
        

    def read_plumed_grid(self, filename):

        import re
        
        #I'll ignore the force for now
        with open(filename, 'r') as f:
            line = f.readline()
            while(line.startswith('#!')):
                if(line.find('TYPE') != -1):
                    self.types = [int(x) for x in re.findall(r'\d{1,}', line)]
                if(line.find('MIN') != -1):
                    self.min = [float(x) for x in re.findall(r'-*\d+\.*\d*', line)]
                if(line.find('MAX') != -1):
                    self.max = [float(x) for x in re.findall(r'-*\d+\.*\d*', line)]
                if(line.find('BIN') != -1):
                    bins = [int(x) + 1 for x in re.findall(r'\d{1,}', line)]
                if(line.find('NVAR') != -1):
                    ncv = [int(x) for x in re.findall(r'\d{1,}', line)]
                    ncv = ncv[0]
                if(line.find('PBC') != -1):
                    self.periodic = [int(x) == 1 for x in re.findall(r'\d{1,}', line)]
                line = f.readline()
        
        #undo that thing that metadynamics does for periodicity 
        for i,p in enumerate(self.periodic):
            if(p):
                bins[i] -= 1

        #now load data
        data = np.genfromtxt(filename)
        assert np.shape(data)[0] == reduce(lambda x,y: x * y, bins, 1), "Number of lines in grid does not match stated bin size: read {}, bins = {} => {}".format(np.shape(data)[0], bins, reduce(lambda x,y: x * y, bins, 1))
        self.pot = data[:,ncv]
        
        self.pot = np.reshape(self.pot, bins)
      

        if(ncv > 1):
            #reflect
            self.pot = self.pot[::-1]

            #rotate
            self.pot = np.rot90(self.pot, 3)
  

                
    def __str__(self):
        return "{} dimension Grid object from {} to {} with {} bins. Periodic = {}, Types = {}".format(self.dims, self.min, self.max, self.nbins, self.periodic, self.types)

    @property
    def dims(self):
        return len(self.min)
    
    @property
    def nbins(self):
        if(self.pot is None):
            return ()
        if(self._nbins is not None):
            return self._nbins
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
        self.meshgrid = None
        self.periodic.append(periodic)
        if(type(name) == type("")):
            self.types.append(Grid.cv_type_map[name])
        else:
            self.types.append(name)
        if(self.pot is None):
            self.pot = np.zeros(bin_number)
        else:
            self.pot = np.resize(self.pot, self.nbins + (int(bin_number),))

    def add_margin(self, percent,pretty=True):
        """Add a margin around the potential without affecting bin number""" 
        assert len(percent) == self.ncv
        length_diff = [(y - x) * (s) for x,y,s in zip(self.min, self.max, percent)]
        new_min = [x - l / 2 for x,l in zip(self.min, length_diff)]
        new_max = [x + l / 2 for x,l in zip(self.max, length_diff)]
        if(pretty):
            new_min = [round(x * 10) / 10. for x in new_min]
            new_max = [round(x * 10) / 10. for x in new_max]
        self.set_min(new_min)
        self.set_max(new_max)
        print 'Set min to {} and max to {}'.format(self.min, self.max) 

        
            
    def set_min(self,min):
        #Change the mins. Fills with previous boundaries
        g = self.clone()
        self._clear()
        for t,m,x,b,p in zip(g.types, min, g.max, g.nbins, g.periodic):
            self.add_cv(t,m,x,b,p)
        self.add(g)


    def set_max(self,max):
        #Change the maxs. Fills with previous boundaries
        g = self.clone()
        self._clear()
        for t,m,x,b,p in zip(g.types, g.min, max, g.nbins, g.periodic):
            self.add_cv(t,m,x,b,p)
        self.add(g)
        #set values to g minimum
        self.pot[np.where(self.pot == 0)] = np.min(g.pot)


    def rescale(self, scale):
        assert len(scale) == self.ncv
        length_diff = [(y - x) * (s - 1) for x,y,s in zip(self.min, self.max, scale)]
        self.min = [x - l / 2 for x,l in zip(self.min, length_diff)]
        self.max = [x + l / 2 for x,l in zip(self.max, length_diff)]

    def to_index(self, x, i):
        return max(0, min(self.nbins[i] - 1, int(floor( (x - self.min[i]) / self.dx[i]))))

    def np_to_index(self, x):
        return np.fmax(np.zeros(np.shape(x)), np.fmin(np.array(self.nbins) - 1, np.floor( (x - np.array(self.min)) / np.array(self.dx))))


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

    def get_value(self, x):
        if(len(x) != self.ncv):
            raise ValueError("Dimension of given x vector does not match grid dimension!")
        index = [0 for xi in x]
        for i, xi in enumerate(x):
            assert xi >= self.min[i] and xi <= self.max[i],"Mesh point is not within grid dimension {}: {}, [{}, {}]".format(i, xi, self.min[i], self.max[i])
            index[i] = self.to_index(xi, i)
        return self.pot[tuple(index)]

    def add(self, other_grid):
        if(self.meshgrid is None):
            self.meshgrid = np.meshgrid(*[np.arange(min, max, dx) for min,max,dx in zip(self.min, self.max, self.dx)], indexing='ij')
        for x in np.nditer(self.meshgrid):
            indexo = other_grid.np_to_index(x)
            indexs = self.np_to_index(x)
            self.pot[tuple(indexs)] += other_grid.pot[tuple(indexo)]
        

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

    def set_bin_number(self, new_shape):
        from scipy.ndimage.interpolation import zoom
        zoom_factor = np.array(new_shape, dtype='float') / self.nbins
        self.pot = zoom(self.pot, zoom_factor, prefilter=True, mode='nearest')
        self.meshgrid = None
                        
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
                if(reduce(lambda x,y: x or y, [x == y for x,y in zip(indices,self.nbins)], False)):
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
            output.write('{: 10.8f} '.format(j * self.dx[i] + self.min[i]))        
        output.write('{: 10.8f}\n'.format(self.pot[tuple(indices)]))

    def _print_grid_end(self, indices, output):
        for i,j in enumerate(indices):
            output.write('{: 10.8f} '.format(j * self.dx[i] + self.min[i]))
        #copy the last bin to the boundary
        indices = [i if i < nb else nb - 1 for i,nb in zip(indices,self.nbins)]
        output.write('{: 10.8f}\n'.format(self.pot[tuple(indices)]))

    
    def plot_2d(self, filename, cmap='gist_earth', resolution=None, axis=(0,1)):
        assert self.dims >= 2
        import matplotlib.pyplot as plt
        old_bins = self.nbins
        if(self.dims > 2):
            #integreate along non-plotting axis
#            raise NotImplementedError()
            data = self.pot[:,:,self.nbins[2] / 2]
        else:
            data = self.pot
        if(resolution is not None):
            self.set_bin_number([resolution if x in axis else self.nbins[x] for x in range(self.dims)])
        plt.imshow(data, interpolation='none', cmap=cmap, extent=[self.min[axis[0]], self.max[axis[0]],self.max[axis[1]],self.min[axis[1]]])
        if(resolution is not None):
            self.set_bin_number(old_bins)
        plt.savefig(filename)


    def normalize(self):
        #make sure we don't have gigantic numbers to start
        self.pot -= np.max(self.pot)
        grids = [np.arange(min, max, dx) for min,max,dx in zip(self.min, self.max, self.dx)]
        Z = np.exp(-self.pot)
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

        #swap in 
        self._nbins = mod_bins
        self.max, mod_max = mod_max, self.max        
                
        self._print_header_array('BIN', np.shape(self.pot), output)
        self._print_header_array('MIN', self.min, output)
        self._print_header_array('MAX', mod_max, output)
        self._print_header_array('PBC', [1 if x else 0 for x in self.periodic], output)
        self._enumerate_grid(lambda x: self._print_grid(x, output), end_fxn=lambda x: self._print_grid_end(x,output))

        #swap out 
        self._nbins = None
        self.max, mod_max = mod_max, self.max        
        

    def add_png_to_grid(self, filename, invert=False):
        if(self.ncv != 2):
            raise ValueError("This method only makes sense on 2D grids. Grid is currently {} dimension".format(self.ncv))
        from pylab import imread, imshow, gray, mean
        a = imread(filename) # read to RGB file
        if(len(np.shape(a)) == 2):
            gray_scale = a
        elif(np.shape(a)[2] == 4):
            gray_scale = mean(a[:,:,0:2],2) * a[:,:,3] # convert to grayscale by multiplication with alpha channel
        else:
            gray_scale = mean(a,2)# convert to gray scale with meana

        if(invert):
            gray_scale = 1 - gray_scale
            
        gray_scale = gray_scale.astype(np.float64)
        #too_small = exp(-10)
        #gray_scale[np.where(gray_scale < too_small)] = too_small        
        #self.pot += np.log(gray_scale)
        old_bins = self.nbins
        self.set_bin_number(np.shape(gray_scale))
        self.pot += gray_scale
        self.normalize()
        self.set_bin_number(old_bins)

def test():
    import sys
    import matplotlib.pyplot as plt
    g = Grid()
    g.add_cv("Distance", 0, 8, 8)
    g.add_cv("Distance", 0, 8, 8)
    g.add_png_to_grid("circle.png")
    g.set_bin_number((512,512))
    plt.imshow(np.exp(g.pot), interpolation='none', cmap='gray')
    plt.savefig("circle_out.png", res=300)


    g = Grid()
    g.add_cv("Absolute position", 0, 10, 128, False)
    g.add_cv("Absolute position", 0, 10, 128, False)
    g.add_png_to_grid("uc_shield.png", invert=True)
    g.set_bin_number((512,512))
    g.set_min([-1,-1])
    g.set_max([11,11])
    plt.imshow(np.exp(g.pot), interpolation='none', cmap='gray')
    plt.savefig("uc_out.png")
    g.write(sys.stdout)    

def load_plumed_grid(filename):
    import matplotlib.pyplot as plt
    g = Grid()
    g.read_plumed_grid(filename)
    g.plot_2d("plot.png")

def plot_free_energy(bias_grid, target, boltzmann_factor=2, bias_factor=1):
    g = Grid()
    g.read_plumed_grid(bias_grid)
    g.pot *= (bias_factor) / (bias_factor - 1)
    g.normalize()

    t = Grid()
    t.read_plumed_grid(target)
    t.pot *= boltzmann_factor    
    t.normalize()

    g.write(open("bias_before_add.grid", "w"))
    t.write(open("target_before_add.grid", "w"))
    g.add(t)
    g.pot *= -1.
    g.plot_2d("free_energy.png")
    g.write(open("free_energy.grid", "w"))
    
    

if __name__ == "__main__":
    #test()
    import sys
    load_grid(sys.argv[1])
