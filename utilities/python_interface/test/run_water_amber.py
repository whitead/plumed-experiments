#!/usr/bin/env python

import sys
sys.path.append('../pylib')
from plumedlib import *

# This is data for a water molecule
nats = 3
mass = [16.,1.,1.]
charge = [0.,0.,0.]
pbc = 0                         # No periodic boundary conditions
box = [0.,0.,0.]
tstep = 0.5                     # The step size in fs.
filename = "Data/plumed.in"          # PLUMED input file name

init_metadyn (nats, mass, tstep, filename, pbc, charge)

# The positions 
pos = []
pos.append([ -2.190,   4.536, -11.667])
pos.append([ -2.609,   4.484, -12.574])
pos.append([ -2.534,   5.435, -11.356])

# The velocities 
vel = []
vel.append([  0.000,   0.000,   0.000])
vel.append([  0.000,   0.000,   0.000])
vel.append([  0.000,   0.000,   0.000])


for istep in xrange(15) :
        # Keep increasing the x-coordinate of atom 1 (O) by factor, each timestep.
        factor = 0.001
        pos[0][0] += factor
        ene, force = cv_calculation(istep, pos, vel, box)
        print 'ene: ',ene
        print 'force[0]: ',force[0]


