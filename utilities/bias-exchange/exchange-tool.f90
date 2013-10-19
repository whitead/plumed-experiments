!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This tool for the script bias-exchange.sh generates a random list
!  of pairs to exchange, and evaluates the exchange probability.
!  Compatible with Plumed 1.3
!
!  version 1.0 - 7 Nov 2011
!  Fabio Pietrucci  (fabio.pietrucci@gmail.com)                 
!  Katsumasa Kamiya (kka2masa@gmail.com)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program exchange_tool
implicit none
!
integer :: argcount,i,nwalker,ind1,ind2
character :: word*80
real*8 :: kBT
! check if there are command-line arguments
argcount=IARGC()
if(argcount==0)then
  write(*,*) 'ERROR: exchange-tool.x invoked without arguments !!!'
  write(*,*) 'USAGE:'
  write(*,*) '  exchange-tool.x -pairs [NWALKER]'
  write(*,*) '    generates a list of random pairs'
  write(*,*) '  exchange-tool.x -try [ind1] [ind2] [kBT]'
  write(*,*) '    tries an exchange between walkers ind1 and ind2'
  write(*,*) 'Exiting.'
endif
! initialize rundom numbers based on the clock
call init_random
! parse command line
do i=1,argcount
  call GETARG(i,word)
  ! generate list of pairs
  if (index(word,'-pairs').ne.0)then
    call GETARG(i+1,word)
    read(word,*) nwalker
    call generate_pairs(nwalker)
    stop
  endif
  if (index(word,'-try').ne.0)THEN
    call GETARG(i+1,word)
    read(word,*) ind1 
    call GETARG(i+2,word)
    read(word,*) ind2 
    call GETARG(i+3,word)
    read(word,*) kBT
    call try_exchange(ind1,ind2,kBT)
    stop
  endif
  if ((index(word,'-pairs').eq.0).and.(index(word,'-try').eq.0))then
    write(*,*) 'ERROR: exchange-tool.x invoked with unknown arguments !!!'
    write(*,*) 'Exiting.'
    stop
  endif
enddo
end program exchange_tool
!
!
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine try_exchange(ind1,ind2,kBT)
    integer :: ind1,ind2,ipair,ind(2),nl,iact(2),nhills(2),nh,nword,i,iacttmp
    character :: filename*50,line*2000,linetmp*2000,answer*10
    real*8 :: t,cv(2,1000),hillcv,sigma,height,ds,kBT,r,VG(2,2),delta,prob
    ind(1)=ind1
    ind(2)=ind2
    ! here the program is verbose
    write(*,'(a,2i4)') "exchange-tool: pair = ",ind(1),ind(2)
    ! read COLVAR
    do ipair=1,2
      if(ind(ipair)<10) write(filename,'(a,i1,a)') "walker",ind(ipair),"/COLVAR"
      if(ind(ipair)>=10.and.ind(ipair)<100) &
      & write(filename,'(a,i2,a)') "walker",ind(ipair),"/COLVAR"
      if(ind(ipair)>=100.and.ind(ipair)<1000) &
      & write(filename,'(a,i3,a)') "walker",ind(ipair),"/COLVAR"
      if(ind(ipair)>=1000)then
        write(*,'(a)') "exchange-tool: ERROR!!! TOO MANY WALKERS !!!"
        stop
      endif 
      open(20,file=filename,status="old")  
      call count_lines(20,nl,iact(ipair)) ! this also finds which CV is active
      do i=1,nl-1
        read(20,*)
      enddo
      read(20,'(A2000)') line
      linetmp=line
      call line2words(linetmp,nword)
      read(line,*) t,cv(ipair,1:nword-1)
      close(20)
      write(*,'(3a,i4,a,f14.6)') "exchange-tool: colvar file = ",trim(filename),&
      & "  active CV = ",iact(ipair),"  value = ",cv(ipair,iact(ipair))
    enddo
    !
    ! read HILLS and sum them
    VG(1:2,1:2)=0.d0
    do ipair=1,2
      if(ind(ipair)<10) write(filename,'(a,i1,a)') "walker",ind(ipair),"/HILLS"
      if(ind(ipair)>=10.and.ind(ipair)<100) &
      & write(filename,'(a,i2,a)') "walker",ind(ipair),"/HILLS"
      if(ind(ipair)>=100.and.ind(ipair)<1000) &
      & write(filename,'(a,i3,a)') "walker",ind(ipair),"/HILLS"
      open(20,file=filename,status="old")  
      call count_lines(20,nl,iacttmp)  
      nh=0
      do i=1,nl
        read(20,'(A2000)') line 
        if (line(1:1)=="#") cycle
        nh=nh+1
        read(line,*) t,hillcv,sigma,height
    ! debug ! write(*,*) "hillcv,sigma,height",hillcv,sigma,height
        ds=0.5*((hillcv-cv(ipair,iact(ipair)))/sigma)**2
    ! debug ! write(*,*) "ds for VG(1)",ds
        if(ds<6.25) then
          ! this is the active walker
          VG(ipair,1)=VG(ipair,1)+height*dexp(-ds)
        endif
        ds=0.5*((hillcv-cv(3-ipair,iact(ipair)))/sigma)**2
    ! debug ! write(*,*) "ds for VG(2)",ds
        if(ds<6.25) then
          ! this is the other walker
          VG(ipair,2)=VG(ipair,2)+height*dexp(-ds)
        endif
      enddo
      nhills(ipair)=nh
      close(20)
    enddo
    ! compute Metropolis acceptance
    ! ( VG(i,j): i = file HILLS, j = 1 (self) or 2 (other) )
    write(*,'(a,f14.6,a,f14.6,a,f14.6)') &
    & "exchange-tool: V_G before exchange = ",&
    & VG(1,1)," + ",VG(2,1)," = ",VG(1,1)+VG(2,1)
    write(*,'(a,f14.6,a,f14.6,a,f14.6)') &
    & "exchange-tool: V_G  after exchange = ",&
    & VG(1,2)," + ",VG(2,2)," = ",VG(1,2)+VG(2,2)
    !
    delta=((VG(1,1)+VG(2,1))-(VG(1,2)+VG(2,2)))/kBT
    if(delta<0.d0)then
      prob=dexp(delta)
    else
      prob=1.d0
    endif
    call random_number(r)
    if(r<=prob)then
      answer="ACCEPT"
    else
      answer="REJECT"
    endif
    write(*,'(a,2f10.7,2x,a)') "exchange-tool: prob, random_num = ",&
    & prob,r,answer
  end subroutine try_exchange
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine generate_pairs(nwalker)
    integer :: i,j
    integer :: nwalker,nw
    integer, allocatable :: uw(:), mw(:,:)
    real*8 :: r
    logical :: found
    ! convert to the closest even integer (greater or equal)
    nw=((nwalker+1)/2)*2
    ! generate a random list of pairs
    allocate(uw(nw),mw(nw,nw))
    uw(:)=0
    mw(:,:)=0
    do i=1,nw
      if(uw(i)==1) cycle
      uw(i)=1
      found=.false.
      do while (.not.found)
        call random_number(r)
        j=int(r*dble(nw))+1
        if(j>0.and.j<nw+1.and.uw(j)==0)then
          mw(i,j)=1
          mw(j,i)=1
          uw(j)=1
          found=.true.
        endif
      enddo
    enddo
    do i=1,nw-1
      do j=i+1,nw
        if(mw(i,j)==1 .and. i<=nwalker .and. j<=nwalker) write(*,'(2i3)') i,j
      enddo
    enddo
  end subroutine generate_pairs
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine count_lines(un,nlines,iact)
    integer :: un,nlines,ios,nact,iact,i
    character :: line*2000
    nact=0
    iact=-1
    nlines=0
    rewind(un)
    do
       read(un,'(A2000)',iostat=ios) line
       if (ios/=0) exit
       nlines=nlines+1
       i=index(line,"ACTIVE")
       if(i.gt.0)then
         read(line(i+7:),*) nact,iact
         if(nact.ne.1)then
           write(*,*) "ERROR!! ACTIVE CVS =",nact
           write(*,*) "you need one active CV"
           write(*,*) "you must employ keyword HILLS_LABEL in each plumed.dat"
           stop
         endif
       endif
    enddo
    rewind(un)
    return
  end subroutine count_lines
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine line2words(line,nword)
    implicit none
    character :: line*2000,word(200)*2000,tmp*2000
    integer :: nword
    integer :: i,p,l,ios,iw
    word(:)=' '
    nword=0
    do iw=1,200
       read(line,*,iostat=ios) tmp
       if (ios/=0) exit
       nword=nword+1
       word(nword)=trim(tmp)
       p=index(line,trim(tmp))
       l=len(trim(tmp))
       do i=p,p+l-1
          line(i:i)=' '
       enddo
    enddo
    return
  end subroutine line2words
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  subroutine init_random
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed
    real*8 :: r
    !
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    call SYSTEM_CLOCK(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
    call random_number(r)
    !!!!!!!
    do i=1,100
      call random_number(r)
    enddo
    !!!!!!!
  end subroutine init_random
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


