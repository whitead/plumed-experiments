#!/usr/bin/env python
"""
 Defines the Plumed class.

 This module defines the class Plumed, which
 can be used 
 as a calculator object in  ASE molecular dynamics.

 @author:       Rosa Bulo
 @organization: Vrije Universiteit Amsterdam
 @contact:      bulo@few.vu.nl
"""

import os, sys, re, math, commands, time, ase, numpy
eV_to_kcal = ase.units.mol/ase.units.kcal

# The following method is stolen from ASE, and changed to include seeding.
def set_maxwell_boltzmann_distribution(atoms, temp, seed=None):
        """Sets the momenta to a Maxwell-Boltzmann distribution."""
        masses = atoms.get_masses()
        random_number_generator = numpy.random.mtrand.RandomState()
        if seed != None :
                random_number_generator.seed(seed)
        xi = random_number_generator.standard_normal((len(masses),3))
        momenta = xi * numpy.sqrt(masses * temp)[:,numpy.newaxis]
        atoms.set_momenta(momenta)

class PlumedSettings (object) :
        """
        Class for representing the settings of a Plumed calculator

        It contains the following variables:

        inputblock:
                A text block containing the PLUMED input file. If such a
                file is there, no further settings need to be specified.
        """

        def __init__ (self, filename) :
                """
                Constructs a PlumedSettings object

                @param filename: The filename of a Plumed input file
                """
                self.inputblock = self.read_inputblock(filename)
 
        def read_inputblock (self, filename) :
                """
                Returns an input text block, by reading from a file

                @oaram filename: The path to a PLUMED input file
                """
                infile = open(filename)
                input = infile.read()
                if len(input) > 0 :
                        inputblock = input

                return inputblock

class Plumed (object) :
        """
        Class for representing a Plumed calculator

        This class which can be used 
        as a calculator object in  ASE molecular dynamics.

        It contains the following variables:

        settings:
                PlumedSettings object containing the settings (input file
                specifications) for the Plumed calculator
        
        dir:                    
                The directory in which the job is being set up (setup_job()). This is
                the root directory for the calculation, which creates and manages
                several subdirectories.
                                
        dirname:        
                The name of the subdirectory in which the actual calculation will
                take place. It's name is related to the name of the job, and this
                directory is created in setup_job()
                                
        forces:
                A list containing forces for all atoms (in eV/Angstrom)

        energy: 
                The potential energy of the system, from the last calculation
                (in eV)

        virial: 
                The virial of the system from the last calculation
                (in eV/Angstrom)
                
        oldcoords:      
                The coordinates from the atoms object as they are at the
                end of a force computation. 

        npt:    
                Boolean set if the forces will be used in a constant pressure simulation
                (get stress should returns something in this case).
        """

        def __init__ (self, settings, job, counter=0) :
                """
                Creates a Plumed object

                @param settings:
                        PlumedForceSettings object containing the settings (input file
                        specifications) for the PlumedForceJob.

                @param job:
                        Another ASE calculator object.

                @param counter:
                        Counts the stepnumber (the number of times this routine is called).
                """
                self.dirname = 'PLUMED'
                self.dir = None
                self.forces = []
                self.oldcoords = []
                self.energy = None
                self.virial = None
                self.npt = False

                self.settings = settings
                self.job = job
                self.counter = counter

        def setup_job (self,atoms) :
                """
                Gets the job read to run

                @param atoms: An ASE Atoms object
                """
                self.dir = os.getcwd()
                out = commands.getoutput('mkdir %s'%(self.dirname))
                os.chdir(self.dirname)
                self.write_moldata(atoms)
                os.chdir(self.dir)
               
                # Start setting up the subjob 
                # Not sure I need this
                if self.npt :
                        self.job.npt = True

        def write_moldata (self, atoms) :
                """
                Writes input and molecule files
                """
                import plumedlib 

                self.write_input()

                # PLUMED initialization
                masses = atoms.get_masses()
                nats = len(masses)
                # The forcejob does not know the timestep size.
                # I guess if I just specify it as 1, I will get
                # the output in units of dt.
                tstep = 1.
                filename = 'plumed.in'
                # Set the periodic boundary conditions based on the box size in 
                # the ASE atoms object
                pbc = True
                cell = atoms.get_cell()
                default_cell = numpy.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                if (cell == default_cell).all().all() == True :
                        pbc = False
                charge = atoms.get_charges()
                plumedlib.init_metadyn (nats, masses, tstep, filename, pbc, charge)

        def write_input (self) :
                """
                Writes the PLUMED input file to disc
                """
                input = self.get_input_block()
                pyfile = open('plumed.in','w')
                pyfile.write(input)
                pyfile.close()

        def run (self,atoms) :
                """
                Runs the PLUMED job (produces PLUMED forces)

                @param atoms: An ASE Atoms object
                """
                import plumedlib

                if self.dir == None :
                        self.setup_job(atoms)

                os.chdir(self.dirname)

                #Now run the subjob
                # Update the positions (all jobs have a molecule attribute)
                # In case this is a pure ASE calculator, we need an atoms object.
                jobforces = self.job.get_forces(atoms)
                jobenergy = self.job.get_potential_energy(atoms)
                # Convert to kcal/mol
                jobforces *= eV_to_kcal
                jobenergy *= eV_to_kcal

                if self.npt :
                        self.virial = self.job.get_stress()

                istep = self.counter
                pos = atoms.get_positions()
                vel = atoms.get_velocities()
                # The box MUST be orthorombic!
                cell = atoms.get_cell()
                if cell[0][1] != 0 or cell[0][2] != 0. :
                        print 'The cell is not orthrombic!'
                        sys.exit(0)
                elif cell[1][0] != 0 or cell[1][2] != 0. :
                        print 'The cell is not orthrombic!'
                        sys.exit(0)
                elif cell[2][0] != 0 or cell[2][1] != 0. :
                        print 'The cell is not orthrombic!'
                        sys.exit(0)
                box = cell[0][0],cell[1][1],cell[2][2]
                # The energy is actually not changed by PLUMED...
                # However, it is in version 1.3, so I should be careful not to use that changed energy
                # further.
                energy, forces = plumedlib.cv_calculation(istep, pos, vel, box, jobforces, jobenergy)
                endrun = time.time()
                #print 'timings: ',endrun-startrun

                # Convert back to eV
                self.forces = forces / eV_to_kcal
                # The energy of the calculator is the energy of its job, not that spit out by the PLUMED 
                # routine.
                self.energy = jobenergy / eV_to_kcal

                self.oldcoords = pos
                self.counter += 1

                os.chdir(self.dir)

        def get_input_block (self) :
                """
                Creates the text of an input file, and returns it

                Currently I am using the text from an existing input file.
                Later I will create my own default settings.
                """
                block = self.settings.inputblock
                return block

        def get_forces (self, atoms) :
                """
                Returns the forces for the atoms object
                
                It only actually runs the calculation if the coordinates have changed 
                since the last call. Otherwise it just uses the forces that
                are stored in self.forces.
        
                @param atoms: 
                        If called by ASE, this is an ASE atoms object, containing
                        the coordinates.
                """
                startrun = time.time()
                newcoords = atoms.get_positions()
                # Now, check if ASE has changed the coordinates
                if len(self.oldcoords) == 0 :
                        self.run(atoms)
                elif not (newcoords == self.oldcoords).all().all() == True :
                        self.run(atoms)
                return self.forces

        def get_potential_energy (self, atoms=None) :
                """
                Returns the potential energy in eV (ASE units)
                """
                return self.energy

        def get_stress (self, atoms=None) :
                """
                Returns the virial computed in self.job

                @note: Sometimes produces a 0 tensor!
                """
                virial = self.virial

                if virial == None :
                        forces = self.forces
                        virial = 0.
                        for f,r in zip(forces,self.oldcoords) :
                                dot_product = 0.
                                for v1,v2 in zip(f,r) :
                                        dot_product += v1 * v2
                                virial += dot_product
                        virial /= 3.

                return virial
