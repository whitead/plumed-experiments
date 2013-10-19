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
program reweight
implicit none

integer              :: i, ii, j, jj, is, ip, ix
integer              :: iostat,argcount, IARGC
integer              :: ncv, nvar, nhisto, nneigh, ntimeout
integer              :: ngrid_var, ngrid_ncv, ngrid_fes
integer              :: ngrid_mini, ngrid_rest, ngrid_ext
integer              :: ix_var, ix_cv, nframe, ix_fes
integer              :: ngrid_redux, ix_var_redux
integer              :: flag, gtab,  dp2index
integer              :: ngauss, nfes, hills_stride, nreject
integer              :: ngrid_histo, bias_nvar, line
integer              :: ix_ext, next, next2
integer              :: header, rewtype, nfix

integer, allocatable :: ix_to_var(:,:), ix_to_minigrid(:,:)
integer, allocatable :: grid(:), idfes(:), var_nD(:), cv_nD(:)
integer, allocatable :: grid_fes(:), fes_nD(:), idfix(:)
integer, allocatable :: ix_to_reweight(:,:), mini_grid(:)
integer, allocatable :: ix_to_histo(:,:), histo_grid(:)
integer, allocatable :: grid_bias(:), grid_ext(:), idext(:), ext_nD(:) 
integer, allocatable :: ix_to_reweight_count(:), full_to_redux(:)

real*8               :: time, hills_W, biasfactor
real*8               :: pi, dp2, xx, expo, dp2cutoff
real*8               :: evol, beta, V_punto_ave, sumP_s
real*8               :: memusage, histo_norm, maxV_t
real*8               :: histo_ratio, seed, Vpot 

real*8, allocatable  :: V_t(:), hills_s0(:), hills_sigma(:)
real*8, allocatable  :: hills_sigold(:)
real*8, allocatable  :: V_punto(:), P_s(:), tab_exp(:)
real*8, allocatable  :: colvar_s0(:), colvar_lbox(:), colvar_dx(:)
real*8, allocatable  :: colvar_min(:), colvar_max(:)
real*8, allocatable  :: fix_min(:), fix_max(:)
real*8, allocatable  :: dp(:), histo(:), fes(:), dp_histo(:)
real*8, allocatable  :: zz(:), yy(:), mini_dx(:), mini_lbox(:)
real*8, allocatable  :: histo_sigma(:), histo_dx(:), histo_lbox(:)
real*8, allocatable  :: bias_min(:), bias_max(:), hh(:)
real*8, allocatable  :: ext_min(:), ext_max(:), V_ext(:)

character*80         :: aux, wq_char, colvar, hills
character*80         :: fes_name, complex_name, numb_to_char
character*80         :: bias_name, stringa, str1, str2
character*80         :: ext_name

logical              :: go(7), go_ok, do_period
logical              :: read_bias, do_ext, nofinal 
logical              :: kjoule, welltemp, fix, histout, gauss
logical, allocatable :: period(:), histo_on(:)

write(*,*)
write(*,*) "           --> REWEIGHT <-- "
write(*,*)

! initialization
go        = .false.
go_ok     = .true.
do_period = .false.
fix       = .false.
read_bias = .false.
do_ext    = .false.
histout   = .false.
gauss     = .false.
welltemp  = .false.
kjoule    = .false.
nofinal   = .false.
!pi        = dacos(-1.d0)
pi        = 3.14159265d0
nhisto    = 1
nneigh    = 1
ncv       = 0 
nvar      = 1
next      = 0
nfix      = 0
ntimeout  = 0 
nreject   = 0
dp2cutoff = 6.25d0
fes_name  = "fes.dat"
seed      = 1.d0
rewtype   = 0

argcount = IARGC()
do i=1,argcount
   call getarg(i,wq_char)
   if (index(wq_char,'-colvar ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) colvar 
     go(1)=.true.
   endif
   if (index(wq_char,'-hills ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) hills 
     go(2)=.true.
   endif
   if (index(wq_char,'-ncv ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) ncv 
     go(3)=.true.
   endif
   if (index(wq_char,'-nvar ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) nvar
     go(4)=.true.
   endif
   if (index(wq_char,'-stride ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) hills_stride 
     go(5)=.true.
   endif 
   if (index(wq_char,'-nhisto ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) nhisto 
   endif
   if (index(wq_char,'-nreject ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) nreject
   endif
   if (index(wq_char,'-temp ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) beta 
     go(6)=.true.
   endif
   if (index(wq_char,'-neigh ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) nneigh 
   endif
   if (index(wq_char,'-timeout ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) ntimeout
   endif
   if (index(wq_char,'-welltemp ').ne.0)then
     welltemp=.true. 
   endif
   if (index(wq_char,'-kjoule ').ne.0)then
     kjoule=.true.
   endif
   if (index(wq_char,'-out ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) fes_name 
   endif
   if (index(wq_char,'-histout ').ne.0)then
     histout=.true.
   endif
   if (index(wq_char,'-gauss ').ne.0)then
     gauss=.true.
     call getarg(i+1,wq_char)
     read(wq_char,*) histo_ratio 
   endif 
   if (index(wq_char,'-seed ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) seed 
   endif
   if (index(wq_char,'-read_bias ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) bias_name
     read_bias=.true. 
   endif
   if (index(wq_char,'-ext_file ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) ext_name
     do_ext=.true.
   endif
   if (index(wq_char,'-nofinal ').ne.0)then
     nofinal=.true.
   endif
   if (index(wq_char,'-rewtype ').ne.0)then
     call getarg(i+1,wq_char)
     read(wq_char,*) rewtype 
   endif
enddo

! check ncv and nvar
if(nvar.lt.ncv) then
  write(*,*) "ERROR! NVAR must be greater than or equal to NCV"
  write(*,*)
  stop
endif

!! in case ncv=0 no HILLS file expected
if(ncv==0) then
 go(5)=.true.
 go(2)=.true.
 nofinal=.true.
endif

allocate(period(nvar),idfes(nvar),grid(nvar),idext(nvar))
allocate(fix_min(nvar), fix_max(nvar), idfix(nvar))
if(ncv/=0) allocate(mini_grid(ncv))

period = .false.
idfes  = 0
idext  = 0
idfix  = 0
grid   = 50 

i=0
do while(i<argcount)
   i=i+1
   call getarg(i,wq_char)
   if (index(wq_char,'-pi ').ne.0) then
       do_period=.true.
       do ix=i+1,argcount
         call getarg(ix,wq_char)
         if(index(wq_char,'-') == 1 ) then 
           i=ix-1
           exit
         else
           read(wq_char,*)ip
           period(ip)=.true.
         endif
       enddo
       cycle 
   endif 
   if (index(wq_char,'-fes ').ne.0) then
       nfes=0
       go(7)=.true.
       do ix=i+1,argcount
         call getarg(ix,wq_char)
         if(index(wq_char,'-') == 1 ) then
           i=ix-1
           exit
         else
           nfes=nfes+1
           READ(wq_char,*)idfes(nfes)
         endif
       enddo
       cycle 
   endif 
   if (index(wq_char,'-fix ').ne.0) then
       nfix=0
       do ix=i+1,argcount
         call getarg(ix,wq_char)
         if(index(wq_char,'-') == 1 ) then
           i=ix-1
           exit
         else
           nfix=nfix+1
           READ(wq_char,*)idfix(nfix)
         endif
       enddo
       cycle
   endif
   if (index(wq_char,'-ext ').ne.0) then
       next=0
       do ix=i+1,argcount
         call getarg(ix,wq_char)
         if(index(wq_char,'-') == 1 ) then
           i=ix-1
           exit
         else
           next=next+1
           READ(wq_char,*) idext(next)
         endif
       enddo
       cycle 
   endif
   if (index(wq_char,'-ngrid ').ne.0) then
    do j=1,nvar
     call getarg(i+j,wq_char)
     read(wq_char,*) grid(j)
    enddo
    i=i+nvar
    cycle
   endif
enddo

i=0
do while(i<argcount)
   i=i+1
   call getarg(i,wq_char)
   if (index(wq_char,'-bound ').ne.0) then
    do j=1,nfix
     call getarg(i+2*j-1,wq_char)
     read(wq_char,*) fix_min(idfix(j))
     call getarg(i+2*j,wq_char)
     read(wq_char,*) fix_max(idfix(j))
    enddo
    i=i+2*nfix
    cycle
   endif
enddo 


do i=1,7
 if(.not.go(i)) go_ok = .false.
enddo

if(.not.go_ok)then
    write(*,*)"** USAGE examples: "
    write(*,*)"1) To use the algorithm described in M. Bonomi, A. Barducci and M. Parrinello, J. Comput. Chem. 30 (2009) pp. 1615:"
    write(*,*)"    reweight -colvar COLVAR -hills HILLS -ncv 1 -nvar 2 -stride 2 -fes 1 2 -temp 300 -welltemp"
    write(*,*)  
    write(*,*)"2) For simulations with a fixed external potential, i.e. umbrella-sampling-like reweighting:"
    write(*,*)"    reweight -colvar COLVAR -ncv 0 -nvar 2 -fes 1 2 -temp 300 -ext 1 2 -ext_file bias.dat"
    write(*,*)
    write(*,*)"3) For unbiased simulations:"
    write(*,*)"    reweight -colvar COLVAR -ncv 0 -nvar 2 -fes 1 2 -temp 300"
    write(*,*)
    write(*,*)"** COMPLETE LIST OF ARGUMENTS:"
    write(*,*)"-hills      : HILLS  filename                           (ex. -hills HILLS)" 
    write(*,*)"-colvar     : COLVAR filename                           (ex. -colvar CV.dat)" 
    write(*,*)"-out        : FES    filename                           (default: -out fes.dat)"
    write(*,*)"-ncv        : # of variables in HILLS                   (ex. -ncv 1)" 
    write(*,*)"-nvar       : # of variables in COLVAR                  (ex. -nvar 2)" 
    write(*,*)"-stride     : ratio between COLVAR and HILLS stride     (ex. -stride 2)" 
    write(*,*)"-fes        : ID of the variables for FES in output     (ex. -fes 1 2)"
    write(*,*)"-temp       : temperature in Kelvin                     (ex. -temp 300.0)"
    write(*,*)"-ngrid      : histogram grid dimension                  (ex. -ngrid 50 100)"                          
    write(*,*)"-nreject    : discard initial steps                     (default: -nreject 0)"
    write(*,*)"-timeout    : stride for FES printout                   (ex. -timeout 1000)"
    write(*,*)"-pi         : ID of variables with [-pi:pi] periodicity (ex. -pi 1 2)"
    write(*,*)"-nhisto     : keep data for histogram every # steps     (default: -nhisto 1)"
    write(*,*)"-fix        : ID of variables with fix boundaries       (ex. -fix 3)"
    write(*,*)"-bound      : fixed boundaries                          (ex. -bound 10.0 12.0)"
    write(*,*)"-gauss      : Gaussian sigma for histogram              (ex. -gauss 1.0)"
    write(*,*)"-read_bias  : starting metad bias filename              (ex. -read_bias bias.dat)"
    write(*,*)"-ext        : ID of variables with external potential   (ex. -ext 1 2)"
    write(*,*)"-ext_file   : external potential filename               (ex. -ext_file bias.dat)"
    write(*,*)
    write(*,*)"** AND FLAGS:"
    write(*,*)"-welltemp   : well-tempered metadynamics"
    write(*,*)"-kjoule     : energy in kjoule/mol                      (default in kcal/mol)"
    write(*,*)"-histout    : print histogram rather than fes"
    write(*,*)
    stop
endif 

write(*,'(1x,a31,a16)') "COLVAR filename             :: ", trim(colvar)
if(ncv/=0) write(*,'(1x,a31,a16)') "HILLS  filename             :: ", trim(hills)
write(*,'(1x,a31,a16)') "FES    filename             :: ", trim(fes_name)
if(read_bias) write(*,'(1x,a31,a16)') "Starting metad bias file    :: ", trim(bias_name)
if(do_ext)    write(*,'(1x,a31,a16)') "External potential file     :: ", trim(ext_name)
if(ncv/=0) write(*,'(1x,a31,10x,i6)') "Number of meta CVs          :: ", ncv
write(*,'(1x,a31,10x,i6)') "Number of variables         :: ", nvar
if(ncv/=0) write(*,'(1x,a31,10x,i6)') "Metadynamics stride         :: ", hills_stride
write(*,'(1x,a31,10x,f6.2)') "Temperature                 :: ", beta
write(*,*) "Grid dimension              :: ", (grid(i),i=1,nvar)
write(*,'(1x,a31,10x,i6)') "Histo stride                :: ", nhisto 
if(ncv/=0) write(*,'(1x,a31,10x,i6)') "Neighbours                  :: ", nneigh
write(*,'(1x,a31,10x,i6)') "Reject                      :: ", nreject
!write(*,'(1x,a31,10x,f6.4)') "Seed                        :: ", seed
if(ntimeout.ne.0) write(*,'(1x,a31,10x,i6)') "Timeout                     :: ", ntimeout
if(ncv/=0) then
 if(welltemp) then 
  write(*,'(1x,a31,14x,a2)')"Well-tempered metad         :: ","on"
 else
  write(*,'(1x,a31,13x,a3)')"Well-tempered metad         :: ","off"
 endif
endif

!! check nfes and nvar 
if(nfes.gt.nvar) then
 write(*,*)
 write(*,*) "ERROR! The number of variables for FES calculation must be lower than or equal to NVAR"
 write(*,*)
 stop
endif
!! check next and nvar 
if(do_ext.and.(next.gt.nvar)) then
 write(*,*)
 write(*,*) "ERROR! The number of variables for external potential must be lower than or equal to NVAR"
 write(*,*)
 stop
endif
!! check read bias and meta on
if(read_bias.and.(ncv==0)) then
 write(*,*)
 write(*,*) "ERROR! You must switch on metadynamics (ncv>0) to read bias"
 write(*,*)
 stop
endif
!! check welltemp
if((.not.welltemp).and.(ncv/=0))  then
 write(*,*)
 write(*,*) "ERROR! Reweighting is supported only for well-tempered metadynamics."
 write(*,*)
 stop
endif
!! check period and fix
if(nfix.gt.nvar) then
 write(*,*)
 write(*,*) "ERROR! The number of fixed variables must be lower than or equal to NVAR"
 write(*,*)
 stop
endif
do i=1,nfix
 if(period(idfix(i))) then
  write(*,*) "ERROR! You cannot fix boundaries for periodic variables"
  stop
 endif 
enddo

!! printout
write(*,*) "Output FES CVs              :: ", (idfes(i),i=1,nfes) 
if(nfix.gt.0) write(*,*) "Fix boundaries CVs          :: ", (idfix(i),i=1,nfix)
if(do_ext) then 
 write(*,*) "External potential on       :: ", (idext(i),i=1,next)
endif
!! temperature unit of measure
write(*,*) 
if(kjoule) then
 write(*,*) "** USING kjoule/mol for energy "
 beta = 1.d0/(0.0083145112119d0*beta)
else
 write(*,*) "** USING kcal/mol for energy " 
 beta = 1.d0/(0.001987191d0*beta)
endif
!! histogram?
if(histout) then
 write(*,*) 
 write(*,*) "** Writing histogram in output instead of fes" 
endif

!! open files
open(unit=100,file=colvar,status="unknown")
if(ncv/=0) open(unit=200,file=hills,status="unknown") 

! if the variable is not periodic
! add one point to the grid
do i=1,nvar
 if(.not.period(i)) grid(i) = grid(i) + 1
enddo

! total grid dimension
ngrid_var  = product(grid(1:nvar))

if(ncv/=0) then
 ngrid_ncv  = product(grid(1:ncv))
! meta allocation
 allocate(V_t(ngrid_ncv),dp(ncv))
 allocate(V_punto(ngrid_ncv))
 allocate(P_s(ngrid_ncv))
 allocate(hills_s0(ncv),hills_sigma(ncv))
 allocate(hills_sigold(ncv))
 allocate(cv_nD(ncv))
 allocate(mini_dx(ncv),yy(ncv),mini_lbox(ncv))
! and initialization
 V_t          = 0.d0
endif

!! allocation
allocate(colvar_min(nvar), colvar_max(nvar), colvar_s0(nvar))
allocate(colvar_lbox(nvar), colvar_dx(nvar), var_nD(nvar))
!! initialization
colvar_min   =  1.d30
colvar_max   = -1.d30 

!! final FES arrays
ngrid_fes=product(grid(idfes(1:nfes)))
allocate(fes(ngrid_fes),grid_fes(nfes), fes_nD(nfes))
grid_fes(1:nfes)=grid(idfes(1:nfes))

!! tabulating exponential function
gtab = 1000000
allocate(tab_exp(gtab))
do i=1,gtab
 expo = dp2cutoff/dble(gtab)*dble(i-1)
 tab_exp(i) = dexp(-expo)
enddo

!! calculate memusage 
memusage=0.d0
! ngrid_ncv arrays
if(ncv/=0) memusage=3.d0*dble(ngrid_ncv)*8.d0
! ngrid_fes
memusage=memusage+dble(ngrid_fes)*8.d0
! tab exp
memusage=memusage+dble(gtab)*8.d0

!! scanning COLVAR for boundaries
write(*,*)
write(*,*) "** SCANNING ",trim(colvar)," for BOUNDARIES"
nframe = 0
do
 read(100,*,iostat=iostat) time, (colvar_s0(i),i=1,nvar) 
 if(iostat/=0)exit
 nframe=nframe+1
 do i=1,nvar
  if(colvar_s0(i).le.colvar_min(i)) colvar_min(i)=colvar_s0(i)
  if(colvar_s0(i).ge.colvar_max(i)) colvar_max(i)=colvar_s0(i) 
 enddo
enddo
rewind(100)

if(ncv/=0) then
 ngauss = 0
!! and HILLS 
 do
  read(200,*,iostat=iostat) time, (hills_s0(i),i=1,ncv)
  if(iostat/=0)exit
  ngauss=ngauss+1
 enddo
 rewind(200)
endif

!! in case of period set min and max to -pi pi
if(do_period) write(*,*) "-> PERIODIC COLLECTIVE VARIABLES"
do i=1,nvar
  if(period(i)) then 
     write(*,'(a9,i2,a46)') "    CV # ",i, " is periodic: overwriting boundaries to -pi,pi"
     colvar_min(i)= -pi
     colvar_max(i)= pi
  endif 
enddo 

!! check for boundaries and fix
 do i=1, nfix 
  j=idfix(i)
  if((fix_min(j).gt.colvar_min(j)).or.(fix_max(j).lt.colvar_max(j))) then
   write(*,*) "ERROR! Fixed boundaries must be larger than the natural ones..."
   do jj=1, nvar
    write(*,'(a7,i2,a4,f15.6,2x,f15.6)') " var # ",jj," :: ",colvar_min(jj),colvar_max(jj)
   enddo
   write(*,*)
   stop
  endif
  colvar_min(j) = fix_min(j)
  colvar_max(j) = fix_max(j)
 enddo

!! final printout
 write(*,*)
 write(*,*) "** BOUNDARIES" 
 do i=1, nvar
  write(*,'(a7,i2,a4,f15.6,2x,f15.6)') " var # ",i," :: ",colvar_min(i),colvar_max(i)
 enddo
 write(*,*)

!! printout and check COLVAR/HILLS
write(*,*) "Total number of Colvar    :: ", nframe
if(ncv/=0) write(*,*) "Total number of Gaussians :: ", ngauss
write(*,*)

if((dble(nframe*nneigh)/dble(ngauss).gt.dble(hills_stride)).and.(ncv/=0)) then
 write(*,*) "ERROR! Mismatch between ",trim(colvar), " and ",trim(hills)," file."
 write(*,*)
 stop
endif

!! setting box stuff 
do i=1,nvar
 colvar_lbox(i) = abs(colvar_max(i) - colvar_min(i))
 if(period(i)) colvar_dx(i)   = colvar_lbox(i) / grid(i)
 if(.not.period(i)) then 
  colvar_dx(i)   = colvar_lbox(i) / (grid(i) - 1)
  colvar_max(i)  = colvar_max(i)  + colvar_dx(i)
  colvar_lbox(i) = colvar_lbox(i) + colvar_dx(i)  
 endif
enddo

if(ncv/=0)  mini_dx(1:ncv) = colvar_dx(1:ncv) 

!! if histogram with gaussians
if(gauss) then
 write(*,*) "** DOING HISTOGRAM WITH GAUSSIANS"
 write(*,*)
 allocate(histo_sigma(nvar), histo_grid(nvar), hh(nvar))
 allocate(dp_histo(nvar), histo_dx(nvar), histo_lbox(nvar))
 histo_sigma(:) = histo_ratio*colvar_dx(:)
 histo_dx(:)    = colvar_dx(:)
 histo_lbox(:)  = dsqrt(2.d0*dp2cutoff)*histo_sigma(:)
 histo_grid(:)  = int(2.d0*histo_lbox(:)/histo_dx(:))+1 
 ngrid_histo    = product(histo_grid(:))
 histo_norm     = dexp(-nvar/2.d0*dlog(2.d0*pi))
 allocate(ix_to_histo(ngrid_histo,nvar))
 memusage=memusage+4.d0*dble(ngrid_histo*nvar)
 do i=1,ngrid_histo
  ix_var=i
  call recover_index(nvar,histo_grid,var_nD,ix_var)
  ix_to_histo(i,:) = var_nD(:)
 enddo
endif

write(*,*) "** CREATING INDEX MAP"
write(*,*)

allocate(histo_on(ngrid_var), full_to_redux(ngrid_var))
ngrid_redux   = 0
histo_on      = .false.
full_to_redux = 0
!! add full_to_redux to memory
memusage=memusage+dble(ngrid_var)*4.d0

if(ncv/=0) then
  allocate(ix_to_reweight_count(ngrid_ncv))
  ix_to_reweight_count = 0
  memusage=memusage+dble(ngrid_ncv)*4.d0
endif

do is=1,nframe
 read(100,*) time, (colvar_s0(i),i=1,nvar)
 if(.not.gauss) then 
  var_nD(:)=int((colvar_s0(:)-colvar_min(:))/colvar_dx(:))+1
  call get_index(nvar,grid,var_nD,ix_var)
  if(.not.histo_on(ix_var)) then
   histo_on(ix_var)=.true.
   ngrid_redux=ngrid_redux+1
   full_to_redux(ix_var)=ngrid_redux
   if(ncv/=0) then
    call get_index(ncv,grid,var_nD(1:ncv),ix_cv)
    ix_to_reweight_count(ix_cv)=ix_to_reweight_count(ix_cv)+1
   endif
  endif
 else
  do ii=1,ngrid_histo
   flag = 0
   hh(:) = colvar_s0(:) - histo_lbox(:) + histo_dx(:)*(ix_to_histo(ii,:)-1) 
   do i=1,nvar
    if(period(i)) hh(i)=hh(i)-colvar_lbox(i)*dnint(hh(i)/colvar_lbox(i))
    if((.not.period(i)).and.((hh(i).lt.colvar_min(i)).or.(hh(i).ge.colvar_max(i)))) flag=1
   enddo
   if(flag==1) cycle
   var_nD(:) = int((hh(:)-colvar_min(:))/colvar_dx(:))+1
   call get_index(nvar,grid,var_nD,ix_var)
   if(.not.histo_on(ix_var)) then
    histo_on(ix_var)=.true.
    ngrid_redux=ngrid_redux+1
    full_to_redux(ix_var)=ngrid_redux
    if(ncv/=0) then
     call get_index(ncv,grid,var_nD(1:ncv),ix_cv)
     ix_to_reweight_count(ix_cv)=ix_to_reweight_count(ix_cv)+1
    endif
   endif
  enddo
 endif
enddo
rewind(100)

allocate(histo(ngrid_redux))
allocate(ix_to_var(ngrid_redux,nvar))
histo     = 0.d0 
ix_to_var = 0
! add histo and ix_to_var memory
memusage=memusage+dble(ngrid_redux)*8.d0+dble(ngrid_redux*nvar)*4.d0

if(ncv/=0) then
 allocate(ix_to_reweight(ngrid_ncv,maxval(ix_to_reweight_count(:)))) 
 ix_to_reweight       = 0
 ix_to_reweight_count = 0
 memusage=memusage+dble(ngrid_ncv*maxval(ix_to_reweight_count(:)))*4.d0
endif

histo_on             = .false.

do is=1,nframe
 read(100,*) time, (colvar_s0(i),i=1,nvar)
 if(.not.gauss) then 
  var_nD(:)=int((colvar_s0(:)-colvar_min(:))/colvar_dx(:))+1
  call get_index(nvar,grid,var_nD,ix_var)
  ix_to_var(full_to_redux(ix_var),:)=var_nD(:)
  if(.not.histo_on(ix_var).and.ncv/=0) then
   histo_on(ix_var)=.true.
   call get_index(ncv,grid,var_nD(1:ncv),ix_cv)
   ix_to_reweight_count(ix_cv)=ix_to_reweight_count(ix_cv)+1 
   ix_to_reweight(ix_cv,ix_to_reweight_count(ix_cv))=full_to_redux(ix_var) 
  endif
 else
  do ii=1,ngrid_histo
   flag = 0
   hh(:) = colvar_s0(:) - histo_lbox(:) + histo_dx(:)*(ix_to_histo(ii,:)-1)
   do i=1,nvar
    if(period(i)) hh(i)=hh(i)-colvar_lbox(i)*dnint(hh(i)/colvar_lbox(i))
    if((.not.period(i)).and.((hh(i).lt.colvar_min(i)).or.(hh(i).ge.colvar_max(i)))) flag=1
   enddo
   if(flag==1) cycle
   var_nD(:) = int((hh(:)-colvar_min(:))/colvar_dx(:))+1
   call get_index(nvar,grid,var_nD,ix_var)
   ix_to_var(full_to_redux(ix_var),:)=var_nD(:)
   if(.not.histo_on(ix_var).and.ncv/=0) then
    histo_on(ix_var)=.true.
    call get_index(ncv,grid,var_nD(1:ncv),ix_cv)
    ix_to_reweight_count(ix_cv)=ix_to_reweight_count(ix_cv)+1 
    ix_to_reweight(ix_cv,ix_to_reweight_count(ix_cv))=full_to_redux(ix_var)
   endif
  enddo
 endif
enddo

deallocate(histo_on)
rewind(100)

!! reading initial bias
if(read_bias) then
 write(*,*) "** Reading starting metad bias from ",trim(bias_name)
 open(unit=400,file=bias_name,status="unknown")

 header = 0
!! first look for NVAR
 read(400,'(a80)') stringa
 do
   read(stringa,*,iostat=iostat) str1,str2
   if(iostat/=0) exit
   if (trim(str1)=="#") then 
    write(*,*) "COMMENT: ",trim(stringa)
   else if (trim(str1)=="#!".and.trim(str2)=="NVAR") then
    read(stringa,*,iostat=iostat) str1,str2,bias_nvar
    if(iostat/=0) exit
    allocate(grid_bias(bias_nvar), bias_min(bias_nvar), bias_max(bias_nvar))
    header = 1
   else if (trim(str1)/="#".and.trim(str1)/="#!") then 
    exit
   endif
100   read(400,'(a80)') stringa
   if(len_trim(stringa)==0) goto 100
 enddo

 if(header/=1) then
  write(*,*) "ERROR: something missing or wrong in the header file"
  stop
 endif

!! now read all the other stuff
 rewind(400)
 read(400,'(a80)') stringa
 do
   read(stringa,*,iostat=iostat) str1,str2
   if (iostat/=0) exit
   if (trim(str1)=="#") then 
    write(*,*) "COMMENT: ",trim(stringa)
   else if (trim(str1)=="#!") then 
    if (trim(str2)=="BIN")  then
     read(stringa,*,iostat=iostat) str1,str2, (grid_bias(i),i=1,bias_nvar)
     if (iostat/=0) exit
     header=header+1
    else if (trim(str2)=="MIN")  then
     read(stringa,*,iostat=iostat) str1,str2, (bias_min(i),i=1,bias_nvar)
     if (iostat/=0) exit
     header=header+1
    else if (trim(str2)=="MAX")  then
     read(stringa,*,iostat=iostat) str1,str2, (bias_max(i),i=1,bias_nvar)
     if (iostat/=0) exit
     header=header+1
    endif
   else if(trim(str1)/="#".and.trim(str1)/="#!") then 
    exit
   endif
200   read(400,'(a80)') stringa
   if(len_trim(stringa)==0) goto 200 
 enddo

 if(header/=4) then
  write(*,*) "ERROR: something missing or wrong in the header file"
  stop
 endif

! write stuff
 write(*,*) "   NVAR ::", bias_nvar
 write(*,*) "   BIN  ::", (grid_bias(i),i=1,bias_nvar)
 write(*,*) "   MAX  ::", (bias_min(i), i=1,bias_nvar)
 write(*,*) "   MIN  ::", (bias_max(i), i=1,bias_nvar)
 write(*,*)

!! checking 
 if(bias_nvar/=ncv) then
  write(*,*) "ERROR: NVAR on bias file inconsistent with -ncv!"
  stop
 endif
!! if not periodic adjust bin and boundaries
 do i=1,ncv
  if(.not.period(i)) then
   grid_bias(i) = grid_bias(i) + 1
   bias_max(i)  = bias_max(i) + colvar_dx(i)
  endif
 enddo
!! checking 
 do i=1,bias_nvar
  if(grid_bias(i)/=grid(i)) then
   write(*,*) "ERROR: BIN on bias file inconsistent with -ngrid"
   stop
  endif
  if(dabs(bias_min(i)-colvar_min(i))>0.0001) then
   write(*,*) "ERROR: MIN on bias file inconsistent with boundaries in input"
   stop
  endif
  if(dabs(bias_max(i)-colvar_max(i))>0.0001) then
   write(*,*) "ERROR: MAX on bias file inconsistent with boundaries in input"
   stop
  endif
 enddo

!! reading the bias potential
 line=0
 do
  line=line+1
  read(stringa,*) (yy(i),i=1,ncv), Vpot 
  cv_nD(:) = int((yy(1:ncv)+colvar_dx(1:ncv)/2.d0-colvar_min(1:ncv))/colvar_dx(1:ncv))+1
  call get_index(ncv,grid,cv_nD,ix_cv)
 
  if(ix_cv.ne.line) write(*,*) "WARNING: the grid on file is not in the usual PLUMED format"

  V_t(ix_cv) = Vpot

300  read(400,'(a80)',iostat=iostat) stringa
  if(iostat/=0) exit
  if(len_trim(stringa)==0) goto 300 
 enddo

 if(line/=ngrid_ncv) then
  write(*,*) "ERROR: inconsistency between actual and declared bias dimension"
  stop
 endif

endif !! if read an initial bias

!! reading an external potential
if(do_ext) then
 write(*,*) "** Reading external potential from ",trim(ext_name)
 open(unit=500,file=ext_name,status="unknown")

 header = 0
!! first look for NVAR
 read(500,'(a80)') stringa
 do
   read(stringa,*,iostat=iostat) str1,str2
   if(iostat/=0) exit
   if (trim(str1)=="#") then 
    write(*,*) "COMMENT: ",trim(stringa)
   else if (trim(str1)=="#!".and.trim(str2)=="NVAR") then
     read(stringa,*,iostat=iostat) str1,str2,next2
     if(iostat/=0) exit
     allocate(grid_ext(next2), ext_min(next2), ext_max(next2))
     header = 1
   else if (trim(str1)/="#".and.trim(str1)/="#!") then
    exit 
   endif
400  read(500,'(a80)') stringa
   if(len_trim(stringa)==0) goto 400
 enddo

 if(header/=1) then
  write(*,*) "ERROR: something missing or wrong in the header file"
  stop
 endif

!! now read all the other stuff
 rewind(500)
 read(500,'(a80)') stringa
 do
   read(stringa,*,iostat=iostat) str1,str2
   if (iostat/=0) exit
   if (trim(str1)=="#") then 
    write(*,*) "COMMENT: ",trim(stringa)
   else if (trim(str1)=="#!") then 
    if (trim(str2)=="BIN") then
     read(stringa,*,iostat=iostat) str1,str2, (grid_ext(i),i=1,next2)
     if (iostat/=0) exit
     header=header+1
    else if (trim(str2)=="MIN") then
     read(stringa,*,iostat=iostat) str1,str2, (ext_min(i),i=1,next2)
     if (iostat/=0) exit
     header=header+1
    else if (trim(str2)=="MAX") then 
     read(stringa,*,iostat=iostat) str1,str2, (ext_max(i),i=1,next2)
     if (iostat/=0) exit
     header=header+1
    endif
   else if (trim(str1)/="#".and.trim(str1)/="#!") then 
    exit
   endif
500   read(500,'(a80)') stringa
   if(len_trim(stringa)==0) goto 500
 enddo

 if(header/=4) then
  write(*,*) "ERROR: something missing or wrong in the header file"
  stop
 endif

! write stuff
 write(*,*) "   NVAR ::", next2
 write(*,*) "   BIN  ::", (grid_ext(i),i=1,next2)
 write(*,*) "   MAX  ::", (ext_min(i), i=1,next2)
 write(*,*) "   MIN  ::", (ext_max(i), i=1,next2)
 write(*,*) 

! checking stuff
 if(next2/=next) then
  write(*,*) "ERROR: NVAR on external potential file inconsistent with -ext!"
  stop
 endif
!! if not periodic adjust bin and boundaries
 do i=1,next
  if(.not.period(idext(i))) then
   grid_ext(i)  = grid_ext(i) + 1
   ext_max(i)   = ext_max(i) + colvar_dx(idext(i))
  endif
 enddo
!! checking 
 do i=1,next
  if(grid_ext(i)/=grid(idext(i))) then
   write(*,*) "ERROR: BIN on external potential file inconsistent with -ngrid"
   stop
  endif
  if(dabs(ext_min(i)-colvar_min(idext(i)))>0.0001) then
   write(*,*) "ERROR: MIN on external potential file inconsistent with boundaries in input"
   stop
  endif
  if(dabs(ext_max(i)-colvar_max(idext(i)))>0.0001) then
   write(*,*) "ERROR: MAX on external potential file inconsistent with boundaries in input"
   stop
  endif
 enddo

 ngrid_ext=product(grid(idext(1:next)))
 allocate(V_ext(ngrid_ext), zz(next))
 allocate(ext_nD(next))

 V_ext = 0.d0

!! adding memory info
 memusage=memusage+dble(ngrid_ext)*8.d0

!! reading the external potential
 line=0
 do
  line=line+1
  read(stringa,*) (zz(i),i=1,next), Vpot 
  ext_nD(1:next) = int((zz(1:next)+colvar_dx(idext(1:next))/2.d0-colvar_min(idext(1:next)))/colvar_dx(idext(1:next)))+1
  call get_index(next,grid_ext,ext_nD,ix_ext) 

  if(ix_ext/=line) write(*,*) "WARNING: external potential on file not in the usual PLUMED format "

  V_ext(ix_ext) = Vpot

600  read(500,'(a80)',iostat=iostat) stringa
  if(iostat/=0) exit
  if(len_trim(stringa)==0) goto 600 
 enddo

 if(line/=ngrid_ext) then
  write(*,*) "ERROR: inconsistency between actual and declared external potential dimension", line, ngrid_ext
  stop
 endif

 close(500)

!! renormalizing external potential
V_ext = V_ext - maxval(V_ext(:))

endif !! if read an external potential 

!! adjusting ntimeout
if(ntimeout==0) ntimeout = nframe+1

!! reading COLVAR (and HILLS) file
do is=1,nframe

  if(modulo(is,10000)==0) write(*,'(a21,i10)') " READING COLVAR # :: ",is

!! reading from COLVAR file
  read(100,*) time, (colvar_s0(i),i=1,nvar)

!! time for updating the histogram??
  if(modulo(is,nhisto)==0.and.(is.gt.nreject)) then
!! calculating index in multi dimensional histogram and in MONO
    if(.not.gauss) then
     var_nD(:)=int((colvar_s0(:)-colvar_min(:))/colvar_dx(:))+1
     call get_index(nvar,grid,var_nD,ix_var)
     ix_var_redux  = full_to_redux(ix_var)
     histo(ix_var_redux) = histo(ix_var_redux) + seed 
    else
     do ii=1,ngrid_histo

      flag = 0
      hh(:) = colvar_s0(:) - histo_lbox(:) + histo_dx(:)*(ix_to_histo(ii,:)-1) 

      do i=1,nvar
       if(period(i)) hh(i)=hh(i)-colvar_lbox(i)*dnint(hh(i)/colvar_lbox(i))
       if((.not.period(i)).and.((hh(i).lt.colvar_min(i)).or.(hh(i).ge.colvar_max(i)))) flag=1
      enddo

      if(flag==1) cycle

      var_nD(:) = int((hh(:)-colvar_min(:))/colvar_dx(:))+1
      call get_index(nvar,grid,var_nD,ix_var)
      ix_var_redux  = full_to_redux(ix_var)
      if(ncv/=0) call get_index(ncv,grid,var_nD(1:ncv),ix_cv)

      hh(:)=colvar_min(:)+colvar_dx(:)*(var_nD(:)-1)
      dp_histo(:)=hh(:)-colvar_s0(:)
      do i=1,nvar
       if(period(i)) dp_histo(i)=dp_histo(i)-colvar_lbox(i)*dnint(dp_histo(i)/colvar_lbox(i))
      enddo

      dp_histo(:)=dp_histo(:)/histo_sigma(:)
      dp2  = 0.5d0*sum(dp_histo(:)**2)

      if(dp2.lt.dp2cutoff) then
       dp2index              = int(dp2*gtab/dp2cutoff)+1
       expo                  = histo_norm*tab_exp(dp2index)
       histo(ix_var_redux)   = histo(ix_var_redux) + expo 
      endif
     enddo

    endif
  endif

!! time for reading a new hill?
  if((ncv/=0).and.(modulo(is,hills_stride)==0)) then
  
   do jj=1,nneigh

    read(200,*) time,(hills_s0(i),i=1,ncv),(hills_sigma(i),i=1,ncv),hills_W,biasfactor
!! if welltempered rescale the Gaussian height back to have the bias
    if(welltemp) hills_W = hills_W * (biasfactor-1.d0) / biasfactor

!! check if minigrid should be updated
    if(sum(abs((hills_sigold(:)-hills_sigma(:))/hills_sigold(:)))/ncv.gt.0.05.or.(is==1)) then
     if(is.eq.hills_stride)  write(*,*) "** CREATING MINIGRID"
     if(is.ne.hills_stride)  write(*,*) "** UPDATING MINIGRID at t=",is
     mini_lbox(:) = dsqrt(2.d0*dp2cutoff)*hills_sigma(:)
     mini_grid(:) = int(2.d0*mini_lbox(:)/mini_dx(:))+1 
     ngrid_mini   = product(mini_grid(:))
     hills_sigold = hills_sigma
     write(*,'(a23,f8.2,a4)') "   TOTAL MEMORY USAGE: ",(memusage+dble(ngrid_mini)*8.d0)/1024.d0/1024.d0," MB "     
     write(*,*)
!! create index vector
     if(allocated(ix_to_minigrid)) deallocate(ix_to_minigrid)
     allocate(ix_to_minigrid(ngrid_mini,ncv))
     do i=1,ngrid_mini
      ix_cv=i
      call recover_index(ncv,mini_grid,cv_nD,ix_cv)
      ix_to_minigrid(i,:) = cv_nD(:)
     enddo
    endif
    
!! Estimate of P(s,t) 
    if(rewtype==0) then
     maxV_t      = maxval(V_t(:))
!! This is an approximation. This expression holds only for t->infinity
     P_s(:)      = dexp(beta*(V_t(:)-maxV_t)/(biasfactor-1.d0))
     sumP_s      = sum(P_s(:))
    elseif(rewtype==1) then
!! This is also an approximation.
     do i=1,ngrid_ncv
      P_s(i)     = sum(histo(ix_to_reweight(i,1:ix_to_reweight_count(i))))
     enddo
     sumP_s      = sum(P_s(:)) 
    endif 

!! put some stuff to zero 
    V_punto     = 0.d0
    V_punto_ave = 0.d0
   
    do ii=1,ngrid_mini

     flag=0

     yy(:)=hills_s0(:) - mini_lbox(:) + mini_dx(:)*(ix_to_minigrid(ii,:)-1)

     do i=1,ncv
      if(period(i)) yy(i)=yy(i)-colvar_lbox(i)*dnint(yy(i)/colvar_lbox(i))
      if((.not.period(i)).and.((yy(i).lt.colvar_min(i)).or.(yy(i).ge.colvar_max(i)))) flag=1
     enddo

     if(flag==1) cycle

     cv_nD(:) = int((yy(:)-colvar_min(1:ncv))/colvar_dx(1:ncv))+1
     call get_index(ncv,grid,cv_nD,ix_cv) 

     yy(:)=colvar_min(1:ncv)+colvar_dx(1:ncv)*(cv_nD(:)-1)
     dp(:)=yy(:)-hills_s0(:)
     do i=1,ncv
      if(period(i)) dp(i)=dp(i)-colvar_lbox(i)*dnint(dp(i)/colvar_lbox(i))  
     enddo
     dp(:)=dp(:)/hills_sigma(:)

     dp2  = 0.5d0*sum(dp(:)**2)

     if(dp2.lt.dp2cutoff) then 
       dp2index        = int(dp2*gtab/dp2cutoff)+1
       expo            = hills_W*tab_exp(dp2index)
       V_t(ix_cv)      = V_t(ix_cv)+expo
       if((is.gt.nreject).and.(rewtype/=2)) then
        V_punto(ix_cv) = expo
        V_punto_ave    = V_punto_ave + expo*P_s(ix_cv)/sumP_s
        evol           = dexp(-beta*expo)
        do i=1,ix_to_reweight_count(ix_cv)
         ix=ix_to_reweight(ix_cv,i) 
         histo(ix)=histo(ix) * evol
        enddo
       endif !! reweighting
     endif

    enddo ! end cycle on minigrid

    if((is.gt.nreject).and.(rewtype/=2)) histo(:) = histo(:) * dexp(beta*V_punto_ave) 

   enddo !! end cycle on neighbour
    
  endif !! end time for Gaussian deposition 

!! time for writing intermediate results?
 if(modulo(is,ntimeout)==0.and.(is.gt.nreject)) then

  if(ncv/=0) maxV_t = maxval(V_t(:))

  write(numb_to_char,'(i15)')is/ntimeout
  numb_to_char=(adjustl(numb_to_char))
  complex_name=trim(fes_name)//'.'//trim(numb_to_char)
  open(2000,file=complex_name,status="unknown")

  fes = 0.d0

  do ii=1,ngrid_redux
     fes_nD(1:nfes)=ix_to_var(ii,idfes(1:nfes))
     call get_index(nfes,grid_fes,fes_nD,ix_fes)
     if(ncv/=0) then
      cv_nD(1:ncv)=ix_to_var(ii,1:ncv)
      call get_index(ncv,grid,cv_nD,ix_cv)
     endif
     if(do_ext.and.(.not.nofinal)) then
      ext_nD(1:next)=ix_to_var(ii,idext(1:next))
      call get_index(next,grid_ext,ext_nD,ix_ext)
      fes(ix_fes)=fes(ix_fes)+histo(ii) * dexp(beta*(V_t(ix_cv)-maxV_t+V_ext(ix_ext)))
     else if(do_ext.and.nofinal) then
      ext_nD(1:next)=ix_to_var(ii,idext(1:next))
      call get_index(next,grid_ext,ext_nD,ix_ext)
      fes(ix_fes)=fes(ix_fes)+histo(ii) * dexp(beta*V_ext(ix_ext)) 
     else if((.not.do_ext).and.nofinal) then
      fes(ix_fes)=fes(ix_fes)+histo(ii)
     else
      fes(ix_fes)=fes(ix_fes)+histo(ii) * dexp(beta*(V_t(ix_cv)-maxV_t))
     endif
  enddo

!! calculating fes
  if(histout) then
   fes = fes/sum(fes(:)) 
  else
   fes = -1.d0/beta*dlog(fes) - minval(-1.d0/beta*dlog(fes(:)))
  endif

  do i=1,ngrid_fes
    ix_fes=i
    call recover_index(nfes,grid_fes,fes_nD,ix_fes)
    write(2000,*) (colvar_min(idfes(ix))+colvar_dx(idfes(ix))*(fes_nD(ix)-1), ix=1, nfes), fes(i) 
    if(fes_nD(1)==grid_fes(1)) write(2000,*)
  enddo

  close(2000)

 endif !! end timeout if

enddo !! end of COLVAR reading

write(*,*)
write(*,*) "** WRITING FINAL FES" 
write(*,*)

open(unit=300,file=fes_name,status="unknown")

if(ncv/=0) V_t = V_t - maxval(V_t(:))

!! final printout
fes = 0.d0

do ii=1,ngrid_redux
  fes_nD(1:nfes)=ix_to_var(ii,idfes(1:nfes))
  call get_index(nfes,grid_fes,fes_nD,ix_fes)
  if(ncv/=0) then
   cv_nD(1:ncv)=ix_to_var(ii,1:ncv)
   call get_index(ncv,grid,cv_nD,ix_cv)
  endif
  if(do_ext.and.(.not.nofinal)) then
   ext_nD(1:next)=ix_to_var(ii,idext(1:next))
   call get_index(next,grid_ext,ext_nD,ix_ext)
   fes(ix_fes)=fes(ix_fes)+histo(ii) * dexp(beta*(V_t(ix_cv)+V_ext(ix_ext)))
  else if(do_ext.and.nofinal) then
   ext_nD(1:next)=ix_to_var(ii,idext(1:next))
   call get_index(next,grid_ext,ext_nD,ix_ext)
   fes(ix_fes)=fes(ix_fes)+histo(ii) * dexp(beta*V_ext(ix_ext))                    
  else if((.not.do_ext).and.nofinal) then
   fes(ix_fes)=fes(ix_fes)+histo(ii)
  else
   fes(ix_fes)=fes(ix_fes)+histo(ii) * dexp(beta*V_t(ix_cv))
  endif
enddo

!! calculating fes
if(histout) then
 fes = fes/sum(fes(:))
else
 fes = -1.d0/beta*dlog(fes) - minval(-1.d0/beta*dlog(fes(:)))
endif

do i=1,ngrid_fes
  ix_fes=i
  call recover_index(nfes,grid_fes,fes_nD,ix_fes)
  write(300,*) (colvar_min(idfes(ix))+colvar_dx(idfes(ix))*(fes_nD(ix)-1), ix=1, nfes), fes(i) 
  if(fes_nD(1)==grid_fes(1)) write(300,*)
enddo

close(300)

end 

subroutine get_index(n,grid,indexnD,index1D)
implicit none
integer  :: i, n, index1D
integer  :: indexnD(n), grid(n)

 index1D = indexnD(1) 
 if(n.gt.1) then 
  do i=2,n
   index1D = index1D+(indexnD(i)-1)*product(grid(1:i-1))
  enddo
 endif

 return
end

subroutine recover_index(n,grid,indexnD,index1D)
implicit none
integer  :: i, n
integer  :: indexnD(n), grid(n), index1D

 if(n.eq.1) then 
  indexnD(1) = index1D
 else
  do i=n,1, -1
   indexnD(i) = index1D / product(grid(1:i-1)) + 1
   index1D    = modulo(index1D, product(grid(1:i-1)))
   if(index1d==0) then
    indexnD(i) = indexnD(i) - 1
    index1d=product(grid(1:i-1))
   endif
  enddo
 endif
 
 return
end
