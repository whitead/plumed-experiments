from plumed_grids import *
import sys


g = Grid()

g.add_cv("Absolute position", -1.5, 1.5, 32, False)
g.add_cv("Absolute position", -1.5, 1.5, 32, False)
g.add_cv("Absolute position", -1.5, 1.5, 32, False)

g.load_data(sys.argv[1])
g.pot *= 10
g.rescale([3 / 1.5, 3 / 1.5, 3 / 1.5])
g.set_min([-4, -4, -4])
g.set_max([4, 4, 4])
g.normalize()
g.write(open("target.grid", 'w'))
