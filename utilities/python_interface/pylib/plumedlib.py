#!/usr/bin/env python

import ctypes, numpy

def atoi (numstring) :
        """
        Converts a string (if it holds a number) to an integer.

        @param numstring: A string holding a number (i.e. '100')

        @note: 
                I included this jusr for testing. It was the 
                simplest function in the library.
        """
        plumedlib = ctypes.cdll.LoadLibrary('libplumed.so')
        num = plumedlib.plumed_atoi(numstring)
        return num

def init_metadyn (nats, mass, tstep, filename, pbc=False, charge=None) :
        '''
        Calls the PLUMED library to set up all global data for the PLUMED run

        @param nats:            The total number of atoms in the system
        @param mass:            A list of masses for all atoms
        @param tstep:           The timestep in fs
        @param filename:        The PLUMED input filename
        @param pbc:             A logical specifying if periodic boundary conditions apply
        @param charge:          A list of charges for all atoms (not needed to all calculations)
        '''

        # Create the string pointer for the PLUMED input filename
        metainp = ctypes.c_char_p(filename)

        # Now create the C arrays for the atom data
        array_of_doubles = ctypes.c_double * nats
        
        c_mass = array_of_doubles()
        for i,m in enumerate(mass) : c_mass[i] = ctypes.c_double(m)
        address_mass = ctypes.byref(c_mass)
        
        c_charge = array_of_doubles()
        for i,c in enumerate(charge) :c_charge[i] = ctypes.c_double(c)
        address_charge = ctypes.byref(c_charge)
      
        # Create the C array for the periodic box 
        array_of_3_doubles = ctypes.c_double * 3
 
        # Create the double pointers
        # For amber, a timestep in ps is expected...
        c_tstep = ctypes.c_double(tstep/1000)
        address_tstep = ctypes.byref(c_tstep)
       
        # Create the integer pointers
        c_nats = ctypes.c_int(nats)
        address_nats = ctypes.byref(c_nats)

        pbc_int = 0
        if pbc : pbc_int = 1
        c_pbc = ctypes.c_int(pbc_int)
        address_pbc = ctypes.byref(c_pbc)

        # This is the unit of energy with respect to kcal/mol
        # But it is set to 1 for amber inside the code (by the
        # precopiler).
        eunit = 0
        c_eunit = ctypes.c_int(eunit)
        address_eunit = ctypes.byref(c_eunit)
        

        # Call the library method
        plumedlib = ctypes.cdll.LoadLibrary('libplumed.so')

        # Sets up the global data for the PLUMED runs.
        # The last argument is added because the routine can also be called in Fortran.
        plumedlib.init_metadyn_(address_nats, address_tstep, address_mass, address_charge, address_pbc, \
                        address_eunit, metainp, len(filename))

        return


def cv_calculation (istep, pos, vel, box=[0.,0.,0.], force=None, energy=0.) :
        """
        Computes the forces and energy coming from the PLUMED rare events job.

        @param pos: A list of the atom positions
        @type pos: list of lists of 3 floats

        @param vel: A list of the atom velocities 
        @type vel: list of lists of 3 floats

        @param box: The size of the orthorombic box
        @type box: list of 3 floats

        @param istep: The timestep
        @type istep: int

        @param force: 
                The forces on the system (without the PLUMED additions)
                These are needed in case of the energy CV.
        @type force: list of lists, or None

        @param energy: The potential energy of the system (used for the ENERGY cv)
        @type energy: float

        @returns: Tuple containing forst the PLUMED energy, and then a list of forces.
        """

        # Create the box array
        array_of_3_doubles = ctypes.c_double * 3

        c_box = array_of_3_doubles()
        for i,l in enumerate(box) : c_box[i] = ctypes.c_double(l)
        address_box = ctypes.byref(c_box)

        # Create the position array

        # Put the position data in a C acceptable form
        nats = len(pos)
        array_of_doubles = ctypes.c_double * (3 * nats)

        c_pos = array_of_doubles()
        counter = 0
        for i in xrange(nats) :
                for j in range(3) :
                        c_pos[counter] = ctypes.c_double(pos[i][j])
                        counter += 1
        address_pos = ctypes.byref(c_pos)

        v_pos = array_of_doubles()
        counter = 0
        for i in xrange(nats) :
                for j in range(3) :
                        v_pos[counter] = ctypes.c_double(vel[i][j])
                        counter += 1
	address_vel = ctypes.byref(v_pos)
      
        c_force = array_of_doubles()
        if force == None : 
                # Now, declare an empty force array and an energy variable
                address_force = ctypes.byref(c_force)
        else :
                counter = 0
                for i in xrange(nats) :
                        for j in range(3) :
                                c_force[counter] = ctypes.c_double(force[i][j])
                                counter += 1
                address_force = ctypes.byref(c_force)
        
        c_ene = ctypes.c_double(energy)
        address_ene = ctypes.byref(c_ene)

        c_istep = ctypes.c_int(istep)
        address_istep = ctypes.byref(c_istep)

        plumedlib = ctypes.cdll.LoadLibrary('libplumed.so')
        plumedlib.meta_force_calculation_(address_box, address_istep, address_pos, address_vel, address_pos,
                                        address_force, address_force, address_force, address_ene)

        # Convert the force array to a more intuitive form ([[x1,y1,z1],[x2,y2,z2],..])
        # It also gets converted to a numpy array object.
        force = numpy.array(c_force)
        #print force
        #force = force.reshape(nats,3)

        # Return the python type value for the energy (not the ctypes object c_ene)
        ene = c_ene.value

        return ene, force


