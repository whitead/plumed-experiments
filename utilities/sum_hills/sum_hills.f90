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
PROGRAM fes_compute
  
  IMPLICIT NONE

  INTERFACE
     RECURSIVE SUBROUTINE fes_compute_low(idim,nn,fes,gauss,ind,ind0,nfes,ndim,ngauss,ngrid,wws,period)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: idim, ngauss, ndim, nfes
       INTEGER, DIMENSION(:) :: nn, ind, ind0
       INTEGER, DIMENSION(:), POINTER :: ngrid
       DOUBLE PRECISION, DIMENSION(:), POINTER   :: fes
       DOUBLE PRECISION, DIMENSION(:,:), POINTER :: gauss
       DOUBLE PRECISION :: wws
       LOGICAL, DIMENSION (:), POINTER :: period
     END SUBROUTINE fes_compute_low

     RECURSIVE SUBROUTINE fes_write(idim, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw, kt)
       IMPLICIT NONE
       INTEGER, INTENT(in) :: idim, ndim, ndw
       INTEGER, DIMENSION(:), POINTER :: pos
       INTEGER, DIMENSION(:), POINTER :: ngrid
       DOUBLE PRECISION, DIMENSION(:), POINTER   :: fes, dp_grid, x0
       DOUBLE PRECISION :: kt
     END SUBROUTINE fes_write

     RECURSIVE SUBROUTINE fes_only_write(idim, fes,  pos,  ndim, ngrid, ndw, kt)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ndw
       INTEGER, INTENT(in) :: idim, ndim
       INTEGER, DIMENSION(:), POINTER :: ngrid
       INTEGER, DIMENSION(:), POINTER :: pos
       DOUBLE PRECISION :: kt
       
       DOUBLE PRECISION, DIMENSION(:), POINTER   :: fes
     END SUBROUTINE fes_only_write
  END INTERFACE

  INTEGER :: i, ix ,id, it, ip
  INTEGER :: ndim, ndw, nt_p, nt, nh, nf, nwr, ncount, naver
  INTEGER :: argcount, ngauss, coor, iargc, nfes
  INTEGER, POINTER :: idw(:)
  INTEGER, POINTER :: i_map(:)
  INTEGER, POINTER :: ngrid(:), tmp(:)
  INTEGER, POINTER :: ind(:), inds(:), nn(:,:), nn_max(:)

  DOUBLE PRECISION :: dum, ss
  CHARACTER(LEN=2000) :: hillsline
  DOUBLE PRECISION, POINTER :: x0w(:),xfw(:)
  DOUBLE PRECISION, POINTER :: ds_cut(:)
  DOUBLE PRECISION, POINTER :: dp_cut(:)
  DOUBLE PRECISION, POINTER :: ss0(:,:)
  DOUBLE PRECISION, POINTER :: mindelta(:)
  DOUBLE PRECISION, POINTER :: delta_s(:,:)
  DOUBLE PRECISION, POINTER :: gauss(:,:),fes(:)
  DOUBLE PRECISION, POINTER :: ww(:), wws(:)
  DOUBLE PRECISION, POINTER :: dp_grid(:),x0(:),xf(:), tmpr(:)
  DOUBLE PRECISION :: dp2,  eps_cut
  DOUBLE PRECISION :: pi

  CHARACTER(LEN=80) :: wq_char, file, out

  LOGICAL :: fix, lstride
  LOGICAL :: l_grid, l_dp, l_grad
  LOGICAL, DIMENSION (:), POINTER :: period, l_tmp
  LOGICAL  :: do_pi, do2_pi

  DOUBLE PRECISION :: kt

  LOGICAL          :: lwellt
  DOUBLE PRECISION :: biasfact

  INTEGER :: rank,npe
  DOUBLE PRECISION, POINTER :: fes_reduce(:)

  CALL PARALLEL_INIT(rank,npe)

!!$ Initialization of the variables
  NDIM    = 2
  NDW     = 2
  nt_p    = 9999999
  eps_cut = -1.d0
  nt      = 0
  kt      = 0.5969d0 ! kT @ 300K in kcal/mol
  pi      = DACOS(-1.d0)

  file    = 'HILLS'
  out     = 'fes.dat'

  fix     = .FALSE.
  lstride = .FALSE.
  l_grid  = .FALSE.
  l_dp    = .FALSE.
  do_pi   = .FALSE.
  do2_pi  = .FALSE.
  l_grad  = .FALSE.
  lwellt  = .FALSE.
  biasfact= -1.d0
  naver   = 0

  argcount = IARGC()
  IF(argcount==0)THEN
     WRITE(6,*)' USAGE:'
     WRITE(6,*)' sum_hills.x -file HILLS -out fes.dat -ndim 3 -ndw 1 2 -kt 0.6 -ngrid 100 100 100'
     WRITE(6,*)' [-ndim 3         ] (number of collective variables NCV)'
     WRITE(6,*)' [-ndw 1 ...      ] (CVs for the free energy surface)'
     WRITE(6,*)' [-ngrid 50 ...   ] (mesh dimension. DEFAULT :: 100)'
     WRITE(6,*)' [-dp ...         ] (size of the mesh of the output free energy)'
     WRITE(6,*)' [-fix 1.1 ...    ] (define the region for the FES, if omitted this is automatically calculated)'
     WRITE(6,*)' [-stride 10      ] (how often the FES is written)'
     WRITE(6,*)' [-cutoff_e 1.e-6 ] (the hills are cutoff at 1.e-6)'
     WRITE(6,*)' [-cutoff_s 6.25  ] (the hills are cutoff at 6.25 standard deviations from the center)'
     WRITE(6,*)' [-2pi x          ] ([0;2pi] periodicity on the x CV, if -fix is not used 2pi is assumed)'
     WRITE(6,*)' [-pi x           ] ([-pi;pi] periodicity on the x CV, if -fix is not used 2pi is assumed)'
     WRITE(6,*)' [-kt 0.6         ] (kT in the energy units)'
     WRITE(6,*)' [-grad           ] (apply periodicity using degrees)'
     WRITE(6,*)' [-bias <biasfact>] (writing output the bias for a well tempered metadynamics run)'
     WRITE(6,*)' [-file HILLLS    ] (input file)'
     WRITE(6,*)' [-out  fes.dat   ] (output file)'
     WRITE(6,*)' [-hills nhills   ] (number of gaussians that are read)'
     WRITE(6,*)' [-aver 100       ] (average in time the bias potential over the last 100 hills)'
     STOP
  ENDIF
  
  DO i=1,argcount
     CALL GETARG(i,wq_char)
 
     IF (INDEX(wq_char,'-file').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)file
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-out').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)out
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-hills').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)nt
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-aver').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)naver
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-ndim').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)NDIM
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-bias').NE.0)THEN
        IF(i<argcount)THEN
           CALL GETARG(i+1,wq_char)
           IF(wq_char(1:1)=="-")THEN
              lwellt=.TRUE.
           ELSE
              READ(wq_char,*)biasfact
           ENDIF
        ELSE
           lwellt=.TRUE.
        ENDIF
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-kt').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)KT
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-stride').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)NT_P
        lstride=.TRUE.
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-cutoff_e').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)eps_cut
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-cutoff_s').NE.0)THEN
        CALL GETARG(i+1,wq_char)
        READ(wq_char,*)dum
        ds_cut=dum
        CYCLE
     ENDIF

     IF (INDEX(wq_char,'-grad').NE.0)THEN
        l_grad=.TRUE.
        CYCLE
     ENDIF
  END DO

  ALLOCATE(ngrid(ndim),dp_grid(ndim),idw(ndim),period(ndim))
  ALLOCATE(mindelta(ndim))

  DO i=1,ndim
     idw(i)=i
     period(i)=.FALSE.
  END DO

  i=0
  DO WHILE(i<argcount)
     i=i+1
     CALL GETARG(i,wq_char)

     IF (INDEX(wq_char,'-ndw ').NE.0)THEN
       ndw=0
       do ix=1,argcount-i
         i=i+1
         CALL GETARG(i,wq_char)
         if(INDEX(wq_char,'-') == 1 ) then
           exit
         else
           ndw=ndw+1
           READ(wq_char,*)idw(ndw)
         endif
       enddo
     ENDIF

     IF (INDEX(wq_char,'-ngrid ').NE.0)THEN
       l_grid=.TRUE.
       do ix=1,ndim
         i=i+1
         CALL GETARG(i,wq_char)
         READ(wq_char,*)ngrid(ix)
       enddo
     ENDIF

     IF (INDEX(wq_char,'-dp ').NE.0)THEN
       do ix=1,ndim
         i=i+1
         CALL GETARG(i,wq_char)
         READ(wq_char,*)dp_grid(ix)
         l_dp=.TRUE.
         l_grid=.FALSE.
       enddo
     ENDIF

     IF (INDEX(wq_char,'-fix').NE.0)THEN
        fix=.TRUE.
        ALLOCATE(x0w(ndw))
        ALLOCATE(xfw(ndw))
        DO id=1,NDW
          i=i+1
          CALL GETARG(i,wq_char)
          READ(wq_char,*)x0w(id)
          i=i+1
          CALL GETARG(i,wq_char)
          READ(wq_char,*)xfw(id)
        ENDDO
     ENDIF

     IF (INDEX(wq_char,'-2pi ').NE.0)THEN
       do2_pi=.TRUE.
       do ix=i+1,argcount
         i=i+1
         CALL GETARG(i,wq_char)
         if(INDEX(wq_char,'-') == 1 ) then
           i=i-1
           exit
         else
           read(wq_char,*)ip
           period(ip)=.TRUE.
         endif
       enddo
     ENDIF

     IF (INDEX(wq_char,'-pi ').NE.0)THEN
       do_pi=.TRUE.
       do ix=i+1,argcount
         i=i+1
         CALL GETARG(i,wq_char)
         if(INDEX(wq_char,'-') == 1 ) then
           i=i-1
           exit
         else
           read(wq_char,*)ip
           period(ip)=.TRUE.
         endif
       enddo
     ENDIF
  ENDDO

!!$  defines the order of the collectiv var.: 
!!$  first the "wanted" ones, then the others

  ALLOCATE(i_map(ndim))
  i_map = 0
  DO id=1,NDW
     i_map(idw(id))=id
  ENDDO
  ix=NDW
  DO id=1,NDIM
     IF(i_map(id)==0)THEN
        ix=ix+1
        i_map(id)=ix
     ENDIF
  ENDDO
  i_map=ndim-i_map+1

  IF(l_grid) THEN
     ALLOCATE(tmp(ndim))
     tmp=ngrid
     DO i=1,ndim
        ngrid(i_map(i))=tmp(i)
     END DO
     DEALLOCATE(tmp)
  ELSE
     ngrid=100
  END IF
  
  IF(any(period)) then
     ALLOCATE(l_tmp(ndim))
     l_tmp=period
     DO i=1,ndim
        period(i_map(i))=l_tmp(i)
     END DO
     DEALLOCATE(l_tmp)
  ELSE
     period=.FALSE.
  ENDIF

  WRITE(*,'(a15,5x,a15)'     )"FILENAME     ::",adjustl(trim(file))
  OPEN(10,file=file,status='old')
  IF( nt == 0 ) THEN
     DO WHILE (.TRUE.)
        READ(10,'(a2000)',END=100,ERR=100) hillsline
        IF(hillsline(1:1)=="#") cycle
        READ(hillsline,*,END=100,ERR=100)(dum,i=1,2+2*ndim)
        if( abs(dum) < 1.d-6 ) cycle
        nt=nt+1
     END DO
100 REWIND(10)
  ENDIF

  ALLOCATE( x0(NDIM)         )
  ALLOCATE( xf(NDIM)         )
  ALLOCATE( ind(NDIM)        )
  ALLOCATE( inds(NDIM)       )
  ALLOCATE( nn_max(NDIM)     )
  ALLOCATE( ds_cut(NDIM)     )
  ALLOCATE( dp_cut(NDIM)     )
  ds_cut=-1.d0

  ALLOCATE( ss0(NDIM,NT)     )
  ALLOCATE( delta_s(NDIM,NT) )
  ALLOCATE( ww(NT)           )
  ALLOCATE( wwS(NT)          )
  ALLOCATE( nn(NDIM,nt)      )

  ss0=0.d0
  delta_s=0.d0
  ww=0.d0
  i=1
  mindelta=10000.d0
  DO WHILE ( i<=nt )
    READ(10,'(a2000)',END=200,ERR=200) hillsline
    IF(hillsline(1:1)=="#") cycle
    IF(lwellt)THEN
      READ(hillsline,*,END=200,ERR=200)dum,(ss0(i_map(id),i),id=1,NDIM),&
                  (delta_s(i_map(id),i),id=1,NDIM),ww(i),biasfact
    ELSE
      READ(hillsline,*,END=200,ERR=200)dum,(ss0(i_map(id),i),id=1,NDIM),&
                  (delta_s(i_map(id),i),id=1,NDIM),ww(i)
    ENDIF
    DO id=1,ndim
      IF(mindelta(id)>delta_s(id,i))mindelta(id)=delta_s(id,i)
    ENDDO
    IF( ABS(ww(i)) < 1.d-6 ) CYCLE
    i=i+1
  END DO
  CLOSE(10)

  IF(biasfact>0.d0)THEN
    ww = ww * (biasfact-1.d0) / biasfact
    IF(out=='fes.dat') out='bias.dat'
  ENDIF

! this implements the time average over the last naver hills 
! (assuming they are equally spaced in time)  
  IF(naver>0)THEN
    DO it=NT-naver+1,NT
      ww(it)=ww(it)*dble(NT-it+1)/dble(naver)
    ENDDO
  ENDIF

  DO id=1,NDIM
     x0(id)=1000000
     xf(id)=-1000000
  ENDDO

  IF(fix) THEN
    DO it=1,NT
      DO id=1,NDIM-NDW
        x0(id)=MIN(x0(id),ss0(id,it)-3.*delta_s(id,it))
        xf(id)=MAX(xf(id),ss0(id,it)+3.*delta_s(id,it))
      ENDDO
    ENDDO
    it=0
    DO id=NDIM,NDIM-NDW+1,-1
       it=it+1
       x0(id)=x0w(it)
       xf(id)=xfw(it)
    ENDDO
  ELSE
    DO it=1,NT
      DO id=ndim,1,-1
        if ( period(id) ) then
          if(do_pi)then
            if(l_grad)then
              x0(id)=-180.d0
              xf(id)= 180.d0
            else
              x0(id)=-pi
              xf(id)= pi
            endif
          elseif(do2_pi)then
            if(l_grad)then
              x0(id)= 0.d0
              xf(id)= 360.d0
            else
              x0(id)= 0.d0
              xf(id)= 2.d0*pi
            endif
          endif
        else
          x0(id)=MIN(x0(id),ss0(id,it)-3.*delta_s(id,it))
          xf(id)=MAX(xf(id),ss0(id,it)+3.*delta_s(id,it))
        endif
      ENDDO
    ENDDO
  ENDIF
  
  IF(l_dp)THEN
     ALLOCATE(tmpr(ndim))
     tmpr=dp_grid
     DO i=1,ndim
        dp_grid(i_map(i))=tmpr(i)
     END DO
     DEALLOCATE(tmpr)
     ngrid=INT((xf-x0)/dp_grid)+1
  ELSE
     dp_grid=(xf-x0)/DBLE(NGRID-1)
  END IF
  
  ip=0
  DO id=ndim,1,-1
    ip=ip+1
    IF(dp_grid(id)>1.25d0*mindelta(id))THEN
      WRITE(*,'("*** WARNING GRID SIZE FOR CV",i2," MAY BE TOO SMALL ***")')ip
      WRITE(*,'("     MINDELTA = ",g15.5)')mindelta(id)
      WRITE(*,'("           DP = ",g15.5)')dp_grid(id)
      WRITE(*,'("ACCEPTABLE DP = ",g15.5)')mindelta(id)*1.25d0
    ENDIF
  ENDDO

  DO i=1,ndim
    IF(period(i)) ngrid(i)=ngrid(i)-1
  ENDDO

  WRITE(*,'(a15,5x,i15)'          )"NDIM         ::",NDIM
  WRITE(*,'(a15,5x,i15)'          )"NDW          ::",NDW
  WRITE(*,'(a15,5x,i15)'          )"HILLS        ::",NT
  if(biasfact>0.d0)WRITE(*,'(a15,5x,g15.5)')"BIASFACTOR   ::",biasfact
  it=ndim
  ip=1
  DO i=1,ndim
     WRITE(*,'(a9,i3,a3,5x,2(g15.5,2x))')"RANGE VAR",ip," ::",x0(it),xf(it)
     it=it-1
     ip=ip+1
  END DO
  WRITE(*,'(a15,5x,10(i15,2x))')     "NGRID        ::",(NGRID(id),id=ndim,1,-1)
  WRITE(*,'(a15,5x,10(g15.5,2x))')  "DX           ::",(dp_grid(id),id=ndim,1,-1)
  WRITE(*,'(a15,5x,10(g15.5,2x))')  "MIN DELTA_S  ::",(minval(delta_s(id,:)),id=ndim,1,-1)

  if(eps_cut<=0.d0 .and. ds_cut(1)<=0.d0)then
    ds_cut=6.25
    WRITE(*,'(a15,5x,10g15.5)' )"CUTOFF (STD) ::",ds_cut(1)

    nn_max = 0
    DO i = 1, nt
       nn(:,i) = INT(delta_s(:,i)*ds_cut/dp_grid)
       wws(i)  = ww(i)/abs(ww(i))
       ww(i)   = abs(ww(i))**(1.d0/DBLE(NDIM))
    END DO

  elseif(ds_cut(1)>0.d0)then
    WRITE(*,'(a15,6x,10g10.5)' )"CUTOFF (STD) ::",ds_cut

    nn_max = 0
    DO i = 1, nt
       nn(:,i) = INT(delta_s(:,i)*ds_cut/dp_grid)
       wws(i)  = ww(i)/abs(ww(i))
       ww(i)   = abs(ww(i))**(1.d0/DBLE(NDIM))
    END DO

  elseif(eps_cut>0.d0)then
    WRITE(*,'(a15,5x,10g10.5)' )"CUTOFF (ERG) ::",eps_cut

    nn_max = 0
    DO i = 1, nt
       if(abs(ww(i)) < 1d-6) cycle
       dp_cut  = SQRT(LOG(ABS(ww(i))/eps_cut)*2.d0)*delta_s(:,i) 
       nn(:,i) = INT(dp_cut/dp_grid)
       wws(i)  = ww(i)/abs(ww(i))
       ww(i)   = abs(ww(i))**(1.d0/DBLE(NDIM))
    END DO

  endif

  nn_max = MAXVAL(nn,DIM=2)
  ngauss = MAXVAL(nn_max) * 2 + 1
  nfes   = PRODUCT(ngrid)

  WRITE(*,'(a15,5x,10g15.5)' )"KbT          ::",KT

  ALLOCATE(gauss(-MAXVAL(nn_max):MAXVAL(nn_max),NDIM))
  ALLOCATE(fes(nfes))
  IF(npe>1) THEN
    ALLOCATE(fes_reduce(nfes))
  ENDIF
  fes=0.d0

  nh=1
  nf=MIN(nh+nt_p-1,nt)
  
  IF (lstride)THEN 
     nwr=nt_p
  ELSE
     nwr=INT(nt/10)+1
  END IF

  ncount = 0

  WRITE(*,'(a22)') "==>  Adding Hills  <=="
  stride : DO WHILE (nh <= nt)

     hills : DO it=nh+rank,nf,npe
        
        if(abs(ww(it)) < 1d-6) cycle
        ind=INT((ss0(:,it)-x0)/dp_grid) + 1
        gauss=0.d0
        
        DO i=1,ndim
           coor = ind(i) - nn(i,it) - 1
           ss  = x0(i) + coor * dp_grid(i) - dp_grid(i)
           DO ip=-nn(i,it),nn(i,it)
              coor = coor + 1
              ss = ss + dp_grid(i)
              IF ( coor .GT. ngrid(i) .and. .not. period(i) ) CYCLE
              IF ( coor .LT. 1        .and. .not. period(i) ) CYCLE
              dp2=((ss-ss0(i,it))/delta_s(i,it))**2
              gauss(ip,i)=ww(it)*EXP(-0.5*dp2)
!              WRITE(*,*)i,ip, gauss(ip,i), ind, coor, ss, x0, ss0(i,it)
           END DO
        END DO
        inds = ind
        CALL fes_compute_low(NDIM,nn(:,it),fes,gauss,ind,inds,nfes,ndim,ngauss,ngrid,wws(it),period)

        IF(.NOT. lstride .AND. MOD(it,nwr)==0)THEN
           WRITE(6,'(a13,i4,a2)') "Hills done ::",INT(10*ANINT(10.*it/nt))," %"
        ELSEIF(.NOT. lstride .AND. it==nt)THEN
           WRITE(6,'(a13,i4,a2)') "Hills done ::",INT(10*ANINT(10.*it/nt))," %"
        END IF

     END DO hills

     IF (lstride) THEN
        ncount = ncount+1
        WRITE(*,'(a13,i5," |-| Hills up to ",i6)') "Done frame ::",ncount,nf
        IF (out=='bias.dat') THEN
          IF(ncount<10) THEN
             WRITE(file,'("bias.dat.",i1)')ncount
          ELSEIF(ncount<100) THEN
             WRITE(file,'("bias.dat.",i2)')ncount
          ELSEIF(ncount<1000) THEN
             WRITE(file,'("bias.dat.",i3)')ncount
          ELSEIF(ncount<10000) THEN
             WRITE(file,'("bias.dat.",i4)')ncount
          ELSEIF(ncount<100000) THEN
             WRITE(file,'("bias.dat.",i5)')ncount
          ELSE
             STOP "Too many..."
          END IF
        ELSE 
          IF(ncount<10) THEN
             WRITE(file,'("fes.dat.",i1)')ncount
          ELSEIF(ncount<100) THEN
             WRITE(file,'("fes.dat.",i2)')ncount
          ELSEIF(ncount<1000) THEN
             WRITE(file,'("fes.dat.",i3)')ncount
          ELSEIF(ncount<10000) THEN
             WRITE(file,'("fes.dat.",i4)')ncount
          ELSEIF(ncount<100000) THEN
             WRITE(file,'("fes.dat.",i5)')ncount
          ELSE
             STOP "Too many..."
          END IF
        END IF
        OPEN(123,file=file)
        ind   = 1
        IF(npe>1) THEN
          CALL PARALLEL_SUM(fes,fes_reduce,size(fes))
          IF(rank==0) CALL fes_write(ndim, fes_reduce,  ind, ndim, ngrid, dp_grid, x0, ndw, kt)
        ELSE
          CALL fes_write(ndim, fes,  ind, ndim, ngrid, dp_grid, x0, ndw, kt)
        END IF
        !CALL fes_only_write(ndim, fes,  ind,  ndim, ngrid, ndw, kt)
        CLOSE(123)

     END IF

     nh=nh+nt_p
     nf=MIN(nh+nt_p-1,nt)
     
  END DO stride
  DEALLOCATE(gauss)

  WRITE(*,'(a22)') "==> Writing output <=="
  OPEN(123,file=out)
  ix=0
!  WRITE(123,'(10g12.5)')(NGRID(id),id=ndim,ndim-ndw+1,-1),ix
  ind   = 1   
  IF(npe>1) THEN
    CALL PARALLEL_SUM(fes,fes_reduce,size(fes))
    IF(rank==0) CALL fes_write(ndim, fes_reduce,  ind, ndim, ngrid, dp_grid, x0, ndw, kt)
  ELSE
    CALL fes_write(ndim, fes,  ind, ndim, ngrid, dp_grid, x0, ndw, kt)
  END IF
  CLOSE(123)
!  WRITE(77,'(g15.9)')fes

  CALL PARALLEL_FINALIZE()

  STOP "NORMAL END"

200 STOP "ERROR IN READING THE HILLS FILE, NOT ENOUGH GAUSSIANS"
    
END PROGRAM fes_compute

RECURSIVE SUBROUTINE fes_compute_low(idim, nn, fes, gauss, ind, ind0, nfes, ndim, ngauss, ngrid,wws,period)
  IMPLICIT NONE
  INTEGER :: i,pnt,j,k
  INTEGER, INTENT(in) :: idim, nfes, ngauss, ndim
  INTEGER, DIMENSION(:) :: nn, ind, ind0
  INTEGER, DIMENSION(:), POINTER :: pos, pos0, ll, ngrid

  DOUBLE PRECISION, DIMENSION(:), POINTER   :: fes
  DOUBLE PRECISION, DIMENSION(:,:), POINTER :: gauss
  DOUBLE PRECISION :: prod, wws
  LOGICAL, DIMENSION (:), POINTER :: period

  ALLOCATE(pos0(ndim))
  ALLOCATE(pos(ndim))
  ALLOCATE(ll(ndim))
  k=nn(idim)
  pos=ind

  DO i = -k, k
    pos(idim)=ind(idim)+i
    IF ( pos(idim) .GT. ngrid(idim) .and. .not. period(idim) ) CYCLE
    IF ( pos(idim) .LT. 1           .and. .not. period(idim) ) CYCLE
    IF(idim/=1) THEN
      CALL fes_compute_low(idim-1, nn, fes, gauss, pos, ind0, nfes, ndim, ngauss, ngrid,wws,period)
    ELSE
      pos0=pos
      do j=1,ndim
        if( period(j) )then
          IF (pos(j) .GT. ngrid(j)) then
            pos0(j)=pos(j)-ngrid(j)
          endif
          IF (pos(j) .LT. 1       ) then
            pos0(j)=pos(j)+ngrid(j)
          endif
        endif
      enddo
      pnt=point(pos0)
      ll=pos-ind0
      prod=1.
      DO j=1,ndim
        prod=prod*gauss(ll(j),j)
      END DO
!      print*,"GAUSS",prod,pos,pos0,ll,(gauss(ll(j),j),j=1,ndim)
      fes(pnt)=fes(pnt)+prod*wws
    END IF
  END DO
  DEALLOCATE(pos0)
  DEALLOCATE(pos)
  DEALLOCATE(ll)

  RETURN
  
CONTAINS
  
  INTEGER FUNCTION point (pos) RESULT(pnt)
    
    INTEGER, DIMENSION (:) :: pos
    INTEGER :: i
    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO
    
  END FUNCTION point
  
END SUBROUTINE fes_compute_low

RECURSIVE SUBROUTINE fes_write(idim, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw, kt)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ndw
  INTEGER :: i, pnt, id, dimval, ip
  INTEGER, INTENT(in) :: idim, ndim
  INTEGER, DIMENSION(:), POINTER :: pos, ngrid
  
  DOUBLE PRECISION, DIMENSION(:), POINTER   :: xx
  DOUBLE PRECISION, DIMENSION(:), POINTER   :: fes, dp_grid, x0
  DOUBLE PRECISION :: erg, kt

  ALLOCATE(xx(ndim))
  xx = x0
  DO i = 1,ngrid(idim)
     pos(idim)=i
     IF(idim/=ndim-ndw+1) THEN
        CALL fes_write(idim-1, fes,  pos,  ndim, ngrid, dp_grid, x0, ndw, kt)
     ELSE        
        pnt=point(pos)
        xx = x0 + dp_grid * ( pos - 1 )
        dimval = PRODUCT(NGRID(1:NDIM-NDW)) !**(NDIM-NDW)
        erg=kt * log( sum( exp( -fes(pnt:pnt+dimval-1)/kt ) ) )
        if (kt>1.e-5) then
          erg=0.d0
          do ip=pnt,pnt+dimval-1
            if(fes(ip)>1.e-6) then
              erg=erg+exp( fes(ip)/kt )
            endif
          enddo
          if(erg>1e-6) then
            erg = -kt * log(erg)
          else
            erg=0.d0
          endif
        else
          erg=MINVAL(-fes(pnt:pnt+dimval-1))
        endif
        WRITE(123,'(10f20.10)')(xx(id),id=ndim,ndim-ndw+1,-1),erg

     END IF
     
  END DO
  WRITE(123,*)
  DEALLOCATE(xx)
  RETURN
  
CONTAINS
  
  INTEGER FUNCTION point (pos) RESULT(pnt)
    
    INTEGER, DIMENSION (:) :: pos
    INTEGER :: i
    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO
    
  END FUNCTION point
  
END SUBROUTINE FES_WRITE

RECURSIVE SUBROUTINE fes_only_write(idim, fes,  pos,  ndim, ngrid, ndw, kt)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ndw
  INTEGER :: i, pnt, dimval, ip
  INTEGER, INTENT(in) :: idim, ndim
  INTEGER, DIMENSION(:), POINTER :: pos, ngrid
  
  DOUBLE PRECISION, DIMENSION(:), POINTER   :: fes
  DOUBLE PRECISION :: erg, kt

  DO i = 1,ngrid(idim)
     pos(idim)=i
     IF(idim/=ndim-ndw+1) THEN
        CALL fes_only_write(idim-1, fes,  pos,  ndim, ngrid, ndw, kt)
     ELSE        
        pnt=point(pos)
        dimval = PRODUCT(NGRID(1:NDIM-NDW)) !NGRID**(NDIM-NDW)
        if (kt>1.e-5) then
          erg=0.d0
          do ip=pnt,pnt+dimval-1
            if(fes(ip)>1.e-6) then
              erg=erg+exp( fes(ip)/kt )
            endif
          enddo
          if(erg>1e-6) then
            erg = -kt * log(erg)
          else
            erg=0.d0
          endif
        else
          erg=MINVAL(-fes(pnt:pnt+dimval-1))
        endif
        WRITE(123,'(1f12.5)')erg
     END IF
     
  END DO
  WRITE(123,*)
  RETURN
  
CONTAINS
  
  INTEGER FUNCTION point (pos) RESULT(pnt)
    
    INTEGER, DIMENSION (:) :: pos
    INTEGER :: i
    pnt=pos(1)
    DO i=2,ndim
       pnt=pnt+(pos(i)-1) * PRODUCT(ngrid(1:i-1))
    END DO
    
  END FUNCTION point
  
END SUBROUTINE FES_ONLY_WRITE
