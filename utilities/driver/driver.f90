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
 PROGRAM driver 

 IMPLICIT none

 logical                  :: go1, go2, go3, go4, do_ext, stop1
 logical                  :: dcd_has_cell, cell_input
 character*1              :: ics
 character*80             :: wq_char, line
 character*2048           :: filename, pdbout, pdbname, dcdname, colvarname
 character*8              :: str  
 integer                  :: i, j, m, kk, pbc_on
 integer                  :: ntap, flag, nredux, nfull, ncv, fncv
 integer                  :: tot_frames, ex_frames
 integer                  :: argcount, IARGC, unit, iostat
 real*8                   :: cell(6), cell_fix(6)
 real*8, allocatable      :: pos_full(:), cv(:), inf(:), sup(:)
 real*8, allocatable      :: occu_full(:), charge_full(:)
 character*5              :: num_ascii
! PDB stuff
 character*6, allocatable :: name(:)
 character*5, allocatable :: atom(:)
 character*3, allocatable :: residue(:)
 character*1, allocatable :: icode(:), altloc(:)
 character*4, allocatable :: numres(:)
 integer, allocatable     :: num(:)
 real*4, allocatable      :: x(:), y(:), z(:)
 real*8, allocatable      :: occu(:), charge(:)
! DCD stuff 
 character*4              :: car4
 character*80             :: car(10)
 integer                  :: nstart, nsanc, nset, ntitle
 integer                  :: charm, namin, i5(5), i9(9)
 real*4                   :: delta
! time 
 real *8                  :: dt 

11   format(a6,a5,1x,a4,a1,a3,1x,a1,a4,a1,3x,f8.3,f8.3,f8.3,f6.2,f6.2)
 
!! initialization
 tot_frames = 0
 ex_frames  = 0
 go1        = .false.
 go2        = .false.
 go3        = .false.
 go4        = .false.
 do_ext     = .false.
 dcdname    = ""
 pdbname    = ""
 pdbout     = ""
 colvarname = "COLVAR"
 argcount   = IARGC()
 ncv        = 1
 pbc_on     = 1 
 cell_input = .false.
 cell_fix   = 0.d0
 cell       = 0.d0
 dt         = 1.d0

 write(*,*)
 write(*,*) "           --> DRIVER <-- "
 write(*,*)


 do i=1,argcount
   call getarg(i,wq_char)
   if (index(wq_char,'-pdb').ne.0)then
     call getarg(i+1,wq_char)
     pdbname=wq_char
     go1=.true.
   endif
   if (index(wq_char,'-dcd').ne.0)then
     call getarg(i+1,wq_char)
     dcdname=wq_char
     go2=.true.
   endif
   if (index(wq_char,'-plumed').ne.0)then
     call getarg(i+1,wq_char)
     filename=wq_char
     go3=.true.
   endif
   if (index(wq_char,'-ncv').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)ncv
     go4=.true.
   endif
   if (index(wq_char,'-out').ne.0)then
     call getarg(i+1,wq_char)
     pdbout=wq_char
   endif
   if (index(wq_char,'-colvar').ne.0)then
     call getarg(i+1,wq_char)
     colvarname=wq_char
   endif
   if (index(wq_char,'-nopbc').ne.0) pbc_on = 0 
   if (index(wq_char,'-cell').ne.0)then
    do j=1,3
     call getarg(i+j,wq_char)
     read(wq_char,*)cell_fix(j)
    enddo
    cell_input =  .true. 
   endif
   if (index(wq_char,'-dt').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*)dt
   endif
 enddo

 allocate(inf(ncv),sup(ncv))
 inf=-1.d40
 sup=1.d40

 i=0
 if (go4) then
  do while(i<argcount)
      i=i+1
      call getarg(i,wq_char) 

      if (index(wq_char,'-interval').ne.0)then
       do_ext = .true.
       do j=1,ncv
        i=i+1
        call getarg(i,wq_char)
        read(wq_char,*) inf(j)
        i=i+1 
        call getarg(i,wq_char)
        read(wq_char,*) sup(j)
       enddo
       if(pdbout=="") pdbout="extract.pdb"
      endif
  enddo
 endif

 if( (do_ext).and.(.not.go4) )then
   write(*,*)"! ERROR - cannot do extract frames unless you also provide -ncv"
   stop1=.true. 
 else
   stop1=.false.
 end if

 if((.not.(go1.and.go2.and.go3)).or.(stop1)) then                      !.and.go4))) then 
  write(*,*)
  write(*,*)"USAGE : "
  write(*,*)"./driver -pdb PDB_FILE -dcd DCD_FILE -plumed PLUMED_input -ncv #" 
  write(*,*)"     (-interval min1 max1 min2 max2 -out OUT_FILE -nopbc -cell a b c)"
  write(*,*) 
  write(*,*)"-pdb      : pdb (one frame) connected to dcd file" 
  write(*,*)"-dcd      : trajectory file" 
  write(*,*)"-plumed   : input file for PLUMED" 
  write(*,*)"-ncv      : number of collective variables                               (optional)"
  write(*,*)"-out      : pdb output filename                                          (optional)" 
!  write(*,*)"-colvar   : colvar output filename [default COLVAR]                      (optional)" 
  write(*,*)"-interval : extract frames with CV in this interval (needs -ncv also)    (optional)" 
  write(*,*)"-nopbc    : don't apply pbc                                              (optional)"
  write(*,*)"-cell     : provide fixed box dimension in Angstrom for orthorhombic PBC (optional)" 
  write(*,*)"-dt       : time interval for dcd frames                                 (optional)"
  write(*,*)
  stop
 endif

 write(*,*)
 write(*,'(a25,a18)') "PDB file name         :: ",trim(pdbname)
 write(*,'(a25,a18)') "DCD file name         :: ",trim(dcdname)
 write(*,'(a25,a18)') "PLUMED input file     :: ",trim(filename)
 write(*,'(a25,a18)') "COLVAR output file     :: ",trim(colvarname)
 write(*,'(a25,13x,i5)') "Number of CVs         :: ",ncv
 write(*,'(a25,13x,f12.6)') "dt for dcd read       :: ",dt
 if(do_ext) then
  write(*,'(a25,a18)') "Output PDB file name  :: ",trim(pdbout)
  write(*,*)
  write(*,'(a8)') "INTERVAL"
  do i=1,ncv
   write(*,*) " CV ",i," :: ",inf(i),sup(i)
  enddo
 endif
 write(*,*)

!! open PDB output
 if(pdbout.ne."") then
  open(17,file=pdbout,form="formatted")
 endif

! cv allocation
 allocate(cv(ncv))

! Reading number of atoms in PDB file
 open (unit=10,file=pdbname,form="formatted")
 nredux = 0 
 do
  read(10,'(a80)',iostat=iostat)line
  if(iostat/=0) exit
  read(line,*)str
  if(index(str,"ATOM").ne.0) nredux = nredux + 1
 enddo
 rewind(10)
 write(*,'(a57,i8)') "The number of atoms in the PDB file is               ::  ", nredux

! allocating stuff
 allocate(name(nredux))
 allocate(atom(nredux))
 allocate(residue(nredux))
 allocate(num(nredux))
 allocate(numres(nredux))
 allocate(x(nredux))
 allocate(y(nredux))
 allocate(z(nredux))
 allocate(occu(nredux))
 allocate(charge(nredux))
 allocate(altloc(nredux))
 allocate(icode(nredux))

!! read PDB 
 kk=0 
 nfull = 0
 do
  read (10,'(a80)',iostat=iostat) line
  if(iostat/=0) exit
  if(index(line,"ATOM").ne.0) then
   kk=kk+1
   read (line,11) name(kk),num_ascii,atom(kk),altloc(kk),residue(kk),ics,numres(kk),icode(kk),x(kk),y(kk),z(kk),occu(kk),charge(kk)
   !! serial contains alphabetic chars for PDBs with > 100k atoms
   read (num_ascii,'(i5)',iostat=iostat) num(kk)
   if(iostat/=0) num(kk)=kk
   if(num(kk).gt.nfull) nfull = num(kk)
  endif
 enddo

 write(*,'(a57,i8)') "The total number of atoms in your system is AT LEAST ::  ", nfull
 if(nredux>nfull) write(*,*) "WARNING: it seems your pdb file has more than one frame: only the first is used"

 allocate(pos_full(3*nfull))
 allocate(occu_full(nfull))
 allocate(charge_full(nfull))
 occu_full   = 0.d0
 charge_full = 0.d0
 do i=1, nredux
  occu_full(num(i))   = occu(i)
  if(occu(i).le.0.d0)then
     write(*,*) "****SEEMS YOUR OCCUPANCY IS SET TO ZERO!!! **** "
     write(*,*) "****THE OCCUPANCY COLUMN IS USED FOR MASSES IN THE DRIVER ****"
     write(*,*) "****PLEASE SET IT TO A VALUE THAT IS LARGER THAN ZERO ******" 
     stop
  endif
  charge_full(num(i)) = charge(i)
 enddo

!!! Reading DCD file    
 open (unit=11,file=dcdname,form="unformatted") 

! Header
 read(11) car4, nset, nstart, nsanc, i5, namin, delta, i9, charm
 read(11) ntitle, (car(i),i=1,ntitle)
 read(11) ntap
 
!! printout info from header
 write(*,'(a57,i8)') "Total number of FRAMES declared in DCD header        ::  ", nset
 write(*,*)

!! check CELL stuff
 if((pbc_on==0).and.cell_input) write(*,'(a65)') "**** Warning !! PBC are OFF. Discarding cell information in INPUT"
 if(pbc_on==1.and.cell_input) write(*,'(a39,3f12.6)') "PBC are ON. FIXED ORTHORHOMBIC CELL :: ",(cell_fix(i),i=1,3)
 if(i9(1)==1) then
   dcd_has_cell = .true.
   if((pbc_on==0).and.(.not.cell_input)) write(*,'(a67)') "**** Warning !! PBC are OFF. Discarding cell information in the DCD"
   if(pbc_on==1.and.(.not.cell_input)) write(*,'(a39)') "PBC are ON. Using cell info from DCD..."
 endif
 if((i9(1)==0)) then 
   dcd_has_cell = .false.
   if(pbc_on==1.and.(.not.cell_input)) then 
    write(*,'(a74)') "**** Error !! PBC are ON but no cell info neither in the DCD nor in INPUT."
    write(*,'(a69)') "Please specify explicitly -cell a b c or turn OFF the PBC with -nopbc" 
    stop
   endif
 endif
 write(*,*)
!! initializing PLUMED 
 call init_metadyn(nfull, occu_full, charge_full, pbc_on, cell_fix, dt,trim(filename)//char(0), fncv)
!! Get ncv from plumed if it is not specified in input
 if (.not.go4) then 
    deallocate(inf,sup,cv)
    ncv=fncv
    allocate(inf(ncv),sup(ncv),cv(ncv))
    inf=-1.d40
    sup=1.d40
    write(*,*) "Found", ncv, "collective coordinates in plumed file"
 end if
 

!! looping on DCD file
 do

  if(dcd_has_cell) then
   read(11,iostat=iostat) cell(1),cell(4),cell(2),cell(5),cell(6),cell(3)
   if(cell_input) cell = cell_fix
   if(iostat/=0) exit
  endif
  read(11,iostat=iostat)(x(m),m=1,nfull)
  if(iostat/=0) exit
  read(11,iostat=iostat)(y(m),m=1,nfull)
  if(iostat/=0) exit
  read(11,iostat=iostat)(z(m),m=1,nfull)
  if(iostat/=0) exit

  pos_full   = 0.d0 
  tot_frames = tot_frames + 1

  if(dcd_has_cell.and.(.not.cell_input).and.pbc_on==1) then
!! check if FAKE cell info 
   if(cell(1)==1..and.cell(2)==1.and.cell(3)==1.) then
    write(*,'(a75)') "**** Error !! Cell info in the DCD are FAKE. Please specify -cell in INPUT"
    write(*,*) 
    stop
   endif
!! check if ORTHORHOMBIC. (TONI) Orthorhombic can be either (0.,0.,0.)
!! [Charmm; NAMD>2.5] or (90.,90.,90.) [older].  See
!! http://mdanalysis.googlecode.com/svn/trunk/src/dcd/dcd.c
   if(maxval(cell(4:6))>0.000000001) then
      if(maxval(cell(4:6)-90.)>0.000000001) then
         write(*,*) "**** Error !! Only orthorhombic cells are supported in this release."  
         write(*,*) "**** Angles in DCD: ", cell(4), cell(5), cell(6)
    write(*,*)
    stop
   endif
   endif
!! printout CELL 
   if(mod(tot_frames,100)==0) write(*,'(a9,i8,2x,a8,6f8.3)') "FRAME :: ",tot_frames,"CELL :: ",cell
  endif

!! calling PLUMED 
  do i=1, nfull
   pos_full(num(i))=dble(x(i))
   pos_full(num(i)+nfull)=dble(y(i))
   pos_full(num(i)+2*nfull)=dble(z(i))
  enddo 
 
  call cv_calculation(cell,pos_full,ncv,cv)
  write(*,*)"CVs are",cv(:)

!! in case of extraction
  if(do_ext) then 

   flag=0
   do i=1,ncv
    if((cv(i).gt.inf(i)).and.(cv(i).lt.sup(i))) flag=flag+1
   enddo 
 
   if(flag.eq.ncv) then
       ex_frames=ex_frames+1
        write(17,*) "REMARK FRAME # ::", tot_frames
        write(17,*) "REMARK CVs ::",(cv(m),m=1,ncv)
        write(17,'(a10,i4)') "MODEL     ",ex_frames-1
        do i=1,nfull
          write (17,11) name(i),num(i),atom(i),altloc(i),residue(i),ics,numres(i),icode(i),x(i),y(i),z(i),occu(i),charge(i)
        enddo
        write(17,'(a3)')"TER"
        write(17,'(a6)')"ENDMDL"
   endif
  
  endif

 enddo !! end loop on DCD

 close(11)

 write(*,*)
 if(do_ext) write(*,'(a57,i8)') "Total number of FRAMES extracted                     ::  ",ex_frames
 write(*,'(a57,i8)')   "Total number of FRAMES in the DCD                    ::  ",tot_frames
 write(*,*)

 end
