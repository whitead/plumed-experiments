#!/bin/bash

#==============================================================================#
#                 SCRIPT FOR BIAS-EXCHANGE SIMULATIONS                         #
#                 compatible with Plumed 1.3                                   #
#                                                                              #
#                 version 1.0 - 7 Nov 2011                                     #
#                 Fabio Pietrucci  (fabio.pietrucci@gmail.com)                 #
#                 Katsumasa Kamiya (kka2masa@gmail.com)                        #
#                                                                              #
#  HOW TO USE: (see also the Plumed 1.3 manual and directory "example")        #
#                                                                              #
#  1) compile exchange-tool.f90 with the command "make"                        #
#     (modify the Makefile if necessary)                                       #
#                                                                              #
#  2) modify the 4 variables in the section here below to suit your setup      #
#                                                                              #
#  3) create NWALKER directories called walker1, walker2, etc.                 #
#                                                                              #
#  4) in each walker-directory prepare all the necessary input for             #
#     a Plumed simulation on one replica, including the command                #
#     HILLS_LABEL [label_for_this_walker]                                      #
#     into the Plumed input file                                               #
#                                                                              #
#  5) in each walker-directory start the simulation,                           #
#     using an executable script called run-walker.sh of this form:            #
#                                                                              #
#       [command to run MD+Plumed, not in background]                          #
#       touch READY                                                            #
#       ../bias-exchange.sh &                                                  #
#                                                                              #
#  6) at this point, all walkers are running: when they will all               #
#     have finished, the present script will take care of                      #
#     performing bias exchanges and re-launching the simulation                #
#     for each walker using the corresponding run-walker.sh script;            #
#     details will be written in bias-exchange.log                             #
#                                                                              #
#==============================================================================#

################################################################################
############ PLEASE MODIFY THE FOLLOWING VARIABLES #############################

# base directory, where this script and exchange-tool.x are located:
BASE_DIR="/home/fabio/Plumed1.3/utilities/bias-exchange/example"

# number of walkers:
NWALKER=4

# k_B*T in units of the MD code:
KBT=0.00095

# maximum number of simulations
NSIMULATIONS=100

############ END ###############################################################
################################################################################

cd $BASE_DIR

# check if all walkers finished their last run
NREADY=0
for ((iw=1;iw<=NWALKER;iw++))
do
  if [ -e walker${iw}/READY ]
  then
    ((NREADY++))
  fi
done
echo "###### NREADY = $NREADY / $NWALKER" >> bias-exchange.log

# if all walkers finished, try to exchange
if [ $NREADY -eq $NWALKER ]
then
  echo "###### START OF BIAS-EXCHANGE ATTEMPT (`date`)" >> bias-exchange.log
  # current configuration
  echo -n "###### CONFIGURATION: " >> bias-exchange.log
  for ((iw=1;iw<=NWALKER;iw++))
  do
    echo -n "${iw}=`grep ACTIVE walker${iw}/COLVAR | awk '{label=$NF}END{print label}'` " >> bias-exchange.log
  done
  echo "" >> bias-exchange.log
  # generate a list of pairs
  ./exchange-tool.x -pairs $NWALKER &> pair-list
  # catch possible problems
  if [ $? -ne 0 ]
  then
    ( echo "*** ERROR calling exchange-tool.x:"
      cat pair-list
      echo "Exiting." ) >> bias-exchange.log
    exit
  fi
  # NPAIRS is equal to NWALKER/2 (truncated integer)
  NPAIRS=$(( NWALKER/2 ))
  # try to perform an exchange for each pair
  for ((ip=1;ip<=NPAIRS;ip++))
  do
    # first walker of this pair
    iw1=`awk -v ip=$ip '{if(NR==ip)print $1}' pair-list`
    # second walker of this pair
    iw2=`awk -v ip=$ip '{if(NR==ip)print $2}' pair-list`
    # evaluate the Metropolis condition for this pair
    ./exchange-tool.x -try $iw1 $iw2 $KBT &> accept
    # catch possible problems
    if [ $? -ne 0 ]
    then
      ( echo "*** ERROR calling exchange-tool.x:"
        cat accept
        echo "Exiting." ) >> bias-exchange.log
      exit
    fi
    # if accepted, do the exchange
    if [ `grep -c ACCEPT accept` == 1 ]
    then
      if [ -e walker${iw1}/HILLS -a -e walker${iw1}/plumed.dat -a -e walker${iw2}/HILLS -a -e walker${iw2}/plumed.dat ]
      then
        # swap the HILLS and plumed input files
        mv walker${iw1}/HILLS         tmpHILLS
        mv walker${iw1}/plumed.dat    tmpplumed.dat
        mv walker${iw2}/HILLS         walker${iw1}
        mv walker${iw2}/plumed.dat    walker${iw1}
        mv tmpHILLS                   walker${iw2}/HILLS
        mv tmpplumed.dat              walker${iw2}/plumed.dat
      else
        ( echo "*** ERROR: not all HILLS and plumed.dat files are in place." 
          echo "Exiting." ) >> bias-exchange.log
        exit
      fi
    fi
    cat accept >> bias-exchange.log
    mv accept accept.old.$ip
  done
  echo "###### END OF BIAS-EXCHANGE ATTEMPT (`date`)" >> bias-exchange.log

  # now we can restart the simulation for each walker, 
  # if the max number of simulations has not been reached
  lastsim=`grep -c STEP bias-exchange.log`
  ((lastsim++))
  if [ $lastsim -le $NSIMULATIONS ]
  then
    echo "###### RESTARTING ALL WALKERS: STEP $lastsim" >> bias-exchange.log
    for ((iw=1;iw<=NWALKER;iw++))
    do 
      cd walker${iw}
      rm READY
      ./run-walker.sh < /dev/null &
      cd ${BASE_DIR}
    done
  else
    echo "###### REACHED THE MAXIMUM NUMBER OF SIMULATIONS = $NSIMULATIONS : THE END !" >> bias-exchange.log
  fi

fi
###### end of bias-exchange script ############################################
