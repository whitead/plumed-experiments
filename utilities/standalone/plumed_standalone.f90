!/*
!*******************************************************************************
!*                                                                             *
!*                                PLUMED                                       *
!*   A Portable Plugin for Free Energy Calculations with Molecular Dynamics    *
!*                              VERSION 1.3                                    *
!*                                                                             *
!*******************************************************************************
!*
!*  
!*  Copyright (c) 2008-2011 The PLUMED team.
!*  See http://www.plumed-code.org for more information. 
!*
!*  This file is part of PLUMED.
!*
!*  PLUMED is free software: you can redistribute it and/or modify
!*  it under the terms of the GNU Lesser General Public License as 
!*  published by the Free Software Foundation, either version 3 of 
!*  the License, or (at your option) any later version.
!*
!*  PLUMED is distributed in the hope that it will be useful,
!*  but WITHOUT ANY WARRANTY; without even the implied warranty of
!*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*  GNU Lesser General Public License for more details.
!*
!*  You should have received a copy of the GNU Lesser General
!*  Public License along with PLUMED.  
!*  If not, see <http://www.gnu.org/licenses/>.
!*
!*  For more info, see:  http://www.plumed-code.org
!*  or subscribe to plumed-users@googlegroups.com
!*
!*/
 PROGRAM standalone 

 IMPLICIT none

 logical      :: go1,go3,go4
 character*300:: wq_char, line
 character*45 :: coordname 
 character*45 :: filename
 character*8  :: str,trash  
 integer      :: i,j,m,numatom,kk
 integer      :: argcount,IARGC,unit,iostat
 integer      :: box_logic
 real*8       :: cell(6),time
 real*8       :: box(3)
 real*8, allocatable  :: pos(:),force(:)
 real*8, allocatable      :: occu(:),charge(:)
 integer  nstep
 real*8   ampli,boltz,tstep,ene
 logical  ampli_on
 

 go1=.false.
 go3=.false.
 coordname=""
 argcount = IARGC()
 box_logic=0
 box=0.d0

 do i=1,argcount

   call getarg(i,wq_char)
   if (INDEX(wq_char,'-coord').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)coordname
     go1=.true.
   endif
   if (INDEX(wq_char,'-plumed').NE.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)filename
     go3=.true.
   endif
 enddo

 if(.not.(go1.and.go3)) then 
  write(*,*)
  write(*,*)"USAGE : "
  write(*,*)"./standalone -coord XYZFILE  -plumed PLUMED_input" 
  write(*,*)"  "
  write(*,*)"  "
  write(*,*)"-plumed     : input file for PLUMED" 
  write(*,*)"-coord      : COORDINATES and other data see below" 
  write(*,*)
  write(*,*)"CELL        : provide box dimension IN COORD" 
  write(*,*)"BOLTZ       : need boltzmann constant in the energy units" 
  write(*,*)"AMPLI      : a factor to convert the native Angstrom CVS" 
  write(*,*)"              as the RMSD (whose template is a PDB in Ang)"
  write(*,*)"              e.g. if my internal CV are in bohr then     "
  write(*,*)"              AMPLI 1.889725989 ;"
  write(*,*)
  write(*,*)"TIME   0.2 0"
  write(*,*)"AMPLI  18.9 "
  write(*,*)"BOLTZ  0.26786879 "
  write(*,*)"CELL  20.3    56.7    78.2   "
  write(*,*)" .      .      .       .           ."
  write(*,*)" .      .      .       .           ."
  write(*,*)" .      .      .       .           ."
  write(*,*)"X[i]   Y[i]   Z[i]   MASS[i]   CHARGE[i]  "
  write(*,*)" .      .      .       .           ."
  write(*,*)" .      .      .       .           ."
  write(*,*)" .      .      .       .           ."
  write(*,*)" "
  write(*,*)" "
  write(*,*)" "
  write(*,*)" and gives back: plumed_forces.dat "
  write(*,*)" ENERGY_FROM_RESTR_AND_HILLS  "
  write(*,*)" .      .      .      "
  write(*,*)" .      .      .      "
  write(*,*)" .      .      .      "
  write(*,*)" .      .      .      "
  write(*,*)"FX[i]  FY[i]  FZ[i]   "
  write(*,*)" .      .      .      "
  write(*,*)" .      .      .      "
  write(*,*)" .      .      .      "
  write(*,*)
  write(*,*)
  stop
 endif
 write(*,*) "-----------------------------------------"
 write(*,*) "-          PLUMED - STANDALONE          -"
 write(*,*) "COORD  input : ",coordname
 write(*,*) "PLUMED input : ",filename
 write(*,*)

! Reading number of atoms in PDB file
 open (unit=10,file=coordname,form="formATTED")
 numatom=0
 do
  read(10,'(a80)',iostat=iostat)line
  if(iostat/=0) exit
  read(line,*)str
  if((index(str,"TIME").eq.0).and.(index(str,"BOX").eq.0).and.&
   &(index(str,"AMPLI").eq.0).and.(index(str,"BOLTZ").eq.0)) numatom=numatom+1
 enddo
 rewind(10)
 write(*,*) "The number of atoms in coord is:",numatom

! allocating stuff
 allocate(occu(numatom))
 allocate(charge(numatom))
 allocate(pos(numatom*3))
 allocate(force(numatom*3))

!! read pdb
 go4=.false.
 kk=0 
 do
  read (10,'(a80)',iostat=iostat) line
!  write (*,*) line
  if(iostat/=0) exit
  if(index(line,"TIME").ne.0)then 
   read (line,*)trash,tstep,nstep
!write (line,*)trash,tstep,nstep
   go4=.true.
  elseif (index(line,"BOX").ne.0)then
   read (line,*)trash,box(1),box(2),box(3)
   box_logic=1; 
  elseif (index(line,"BOLTZ").ne.0)then
   read (line,*)trash,boltz
  elseif (index(line,"AMPLI").ne.0)then
   read (line,*)trash,ampli
  else 
   kk=kk+1
   read (line,*)pos(kk),pos(kk+numatom),pos(kk+2*numatom),occu(kk),charge(kk)
!write(*,*) pos(kk),pos(kk+numatom),pos(kk+2*numatom),occu(kk),charge(kk)
  endif
 enddo
 write(*,*) "-----------------------------------------"
 if(.not.go4) then
 write(*,*) "ABORTING: need time indicators...   TIME timestep nstep"
 stop
 endif 
!! initializing PLUMED 
 call init_metadyn(numatom, occu, charge, box_logic, box, tstep, nstep , boltz, ampli , trim(filename)//char(0))
 call cv_calculation_standalone(box,pos,force,ene)

 OPEN(unit=34,FILE="plumed_forces.dat",FORM="FORMATTED") 
 write(34,*)ene
 do i=1,numatom
    write(34,*)force(i),force(i+numatom),force(i+2*numatom)
 enddo 
 close(34)
 stop
 end
