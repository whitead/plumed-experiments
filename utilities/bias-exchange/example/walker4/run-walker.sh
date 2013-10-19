#!/bin/bash

# this script must be executable (e.g.: chmod a+x run-walker.sh)

# note that CPMD+Plumed is launched in foreground (no & at the end of the line):
/home/fabio/CPMD-3.15.1/cpmd.x input . -plumed >> output

# this creates and empty file called READY,
# telling the bias-exchange tool that the MD run has finished:
touch READY

# note that bias-exchange.sh is launched in background:
../bias-exchange.sh &

