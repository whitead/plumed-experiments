from plumed_grids import *
import sys


g = Grid()

g.add_cv("Absolute position", -1., 1., 64, False)
g.add_cv("Absolute position", -1, 1, 64, False)
g.add_cv("Absolute position", 0, 2, 64, False)
g.load_data(sys.argv[1])
g.pot *= 2
g.normalize()
#g.rescale([2, 2,2])
g.add_margin([0.25, 0.25, 0.25])
g.plot_2d("plot.png")
g.write(open("target.grid", 'w'))
