#!/usr/bin/env python

import ase, numpy, time, shutil, sys
from ase.md import npt
sys.path.append('../pylib')
from ase_calculators import PlumedSettings, Plumed, set_maxwell_boltzmann_distribution

# Set up the plumed calculator (containing EMT)
emtjob = ase.calculators.emt.EMT()
plumedsettings = PlumedSettings( 'Data/plumed.in' )
forcejob = Plumed( settings=plumedsettings, job=emtjob )

# This is data for a water molecule
labels = ['O','H','H']
pos = [] 
pos.append([ -2.190,   4.536, -11.667])
pos.append([ -2.609,   4.484, -12.574])
pos.append([ -2.534,   5.435, -11.356])
ASEatoms = ase.Atoms(symbols=labels,positions=pos)
# Assign random velocities to the atoms
set_maxwell_boltzmann_distribution(ASEatoms,300*ase.units.kB,seed=1)
# Give the calculator to the atoms object
ASEatoms.set_calculator(forcejob)

# Choose the integration method
# Here we choose Nose Hoover, at room temperature
ts = 0.5 #timestep
# Room temperature in eV
temp = 1./40.
# Room pressure (1 atm)
pressure = 101325 * ase.units.Pascal
ttime = 25*ase.units.fs
pfact=None
ASEmd = npt.NPT(ASEatoms,ts*ase.units.fs,temp,pressure,ttime,pfact)

# PRINTING
# Print coordinates and energy at t=0
outfile = open('ENERGY.xyz','w')
nats = len(labels)
coords = ASEatoms.get_positions()
ASEatoms.get_forces()
ene = ASEatoms.get_potential_energy()
outfile.write('%i\n'%(nats))
outfile.write('Energy (cycle %4i): %16.5f eV\n'%(0,ene))
for l,c in zip(labels,coords) :
        outfile.write('%6s %16.5f %16.5f %16.5f\n'%(l,c[0],c[1],c[2]))

# Start MD, and print coordinates and energy at each timestep
for i in range(1000) :
        starttime = time.time()
        outfile.write('%i\n'%(nats))
        ASEmd.run(1)
        coords = ASEatoms.get_positions()
        ene = ASEatoms.get_potential_energy()
        outfile.write('Energy (cycle %4i): %16.5f eV\n'%(i+1,ene))
        for l,c in zip(labels,coords) :
                outfile.write('%6s %16.5f %16.5f %16.5f\n'%(l,c[0],c[1],c[2]))
        endtime = time.time()
        print 'cycle %4i'%(i+1)
        print 'tottime = %16.5f'%(endtime-starttime)

outfile.close()

#shutil.move('ENERGY.xyz','Outputs/ASEnvt_plumed_emt_water.xyz')
#shutil.move('PLUMED/COLVAR','Outputs/ASEnvt_plumed_emt_water.col')
#shutil.move('PLUMED/HILLS','Outputs/ASEnvt_plumed_emt_water.hil')
#shutil.rmtree('PLUMED')
