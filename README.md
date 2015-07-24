Plumed-Experiment
=============

A fork of PLUMED 1.3.0 that has experiment directed metadynamics and
experiment directed simulation implementations. Also some small
additions such as McGovern-De Pablo boundary corrected hills, global
tempering, and thresholding. Written by Andrew White and James Dama.

New Options
====

The new options are described here instead of the manual. Please see
the manual for how to use the usual options in Plumed 1.3.


Experiment Directed Simulation
----

**Theory**

Experiment directed simulation applies an adaptively fit linear bias
to a CV so that it matches a set value. If a non-NVE ensemble is being
used, just adding time to allow the bias to settle is sufficient to
remove the effect of adaptively finding the bias. If energy
conservation is essential (NVE), first find the bias in NVT and then
fix it to the `average` value from the NVT simulation. This is a brief
overview, please see **Efficient and minimal method to bias molecular
simulations with experimental data **. AD White, GA
Voth. *J. Chem. Theory Comput.* **2014**, *10 (8)*, pp 189-194 for
complete information.

**Input**

To do experiment directed simulation with an adaptive bias, here's
what the input should look like:

```
#Sets the temperature, random number seed, update period/stride, output filename and gives a list of CVs to bias
EDS STRIDE 500 SIMTEMP 300 SEED 4902398 FILENAME EDS_OUT CV LIST 1 2
#This is where the desired average value for the CVs is set
EDS CV CENTERS 0.2 2.5
#Define two CVs that will be biased
TORSION LIST 5 7 9 15 
TORSION LIST 7 9 15 17
```

Choosing the `STRIDE` is discussed in more depth in the White and Voth
paper, but essentially it should be on the order of the correlation
time of you CV. If the `FILENAME` is not given, the EDS bias will be
written out to `EDS_OUT`. The `CV LIST` command specifies which CVs
are being biased. Additional options not shown here are setting the
`EDS CV RANGES`, which allows setting the magnitude of the first bias
change. If the bias is constantly increasing or constantly decreasing,
it means it's having no effect on the system and its magnitude should
be increased. In the example above, you could add a new line such as

```
EDS CV RANGES 10 10
```

By default the CV range is 1 kT. You may also restart a bias by adding
`RESTART EDS_OUT_RES` to the first `EDS` line. The restart is not
exact if the simulation ends in the middle of an EDS period, however
the change is negligible except for the first ~5 periods of EDS.

To create an input file with a fixed bias:

```
EDS SIMTEMP 300 CV LIST 1 2
EDS CV CENTERS 0.2 2.5
EDS CV CONSTANTS 331.0 10.5
TORSION LIST 5 7 9 15 
TORSION LIST 7 9 15 17

```

The extra line is the `EDS CV CONSTANTS`, where the bias constants are
defined. You may also add `RAMP 500` to the `EDS` line, which will
cause the simulation to ramp to the set bias over the first 500 steps
from 0.


**Output File**

The `EDS_OUT` file contains the step in the first column, followed by
the bias coefficients. The bias force applied to the CV is `-a`, where
`a` is the coefficient in the EDS column and `a CV` is the bias
energy, where CV is the value of the CV.



Experiment Directed Metadynamics
----

EDS is for matching average values. EDM is for
matching whole PMFs. It is described in:

**Designing Free Energy Surfaces that Match Experimental Data with Metadynamics**. AD White, JF Dama, GA Voth. *J. Chem. Theory Comput.* **2015**, *11 (6)*, pp 2451-2460


To use EDM, you must have a grid containing the PMF you would like to
match in the plumed grid format. The units of the PMF should be kT, so
that to convert to a probability `e^(-f(x))` would be applied to the
grid where `f(x)` is the value in the grid. See the
[plumed_grids](https://github.com/whitead/plumed_grids) for making
such grids. EDM uses metadynamics, so please read the manual on metadynamics and
see one of the many papers written on the topic.

**Input**

The input below would match the phi/psi angles on alanine dipeptide to the given grid.

```
#These are metadynamics parameters
HILLS HEIGHT 0.2 W_STRIDE 60
WELLTEMPERED SIMTEMP 300 BIASFACTOR 10

#This is the extra part added for EDM
TARGET_DISTRIBUTION FILENAME target.dat

#The usual PLUMED definition of CVs
PRINT W_STRIDE 10
TORSION LIST 5 7 9 15 SIGMA 0.35
TORSION LIST 7 9 15 17 SIGMA 0.35

#CVS must be gridded to use EDM 
GRID CV 1 MIN -pi MAX pi NBIN 300 PBC 
GRID CV 2 MIN -pi MAX pi NBIN 300 PBC
```

**Restrictions**

Your target should always span the CV space. See
[electronic-dance-music](https://github.com/whitead/electronic-dance-music)
for an implementation that relaxes this requirement. If you're doing
targeted on a non-periodic CV, you should have a method that prevents
the CV from leaving the targeted space. This may be done by just
making the target have high free energy outside of a specific
region. You should use McGovern-DePablo hills as well and turn off
splines. Here's an example of a non-periodic CV:

```
GRID CV 1 MIN 0 MAX 10 NBIN 1024
NOSPLINE
INTERVAL CV 1 LOWER_LIMIT 0 UPPER_LIMIT 10
MCGDP_HILLS CV 1 LOWER_BOUND 0 UPPER_BOUND 10
```

CP2K Notes
=========

There is a patch for the 2.5 branch of cp2k included in the source.
