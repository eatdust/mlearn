program shepmain
  use shepprec
  use shepnd
  use ioml
  use ioshep
  implicit none
 
  integer :: ndim
  integer(ip) :: ndata
  
  integer(ip), dimension(:), allocatable, save :: ictrs
  integer(ip) :: nctrs
  integeR(ip) :: nfits

  integer, parameter :: ndump = 100

  real(fp), dimension(:,:), pointer :: xdata => null()
  real(fp), dimension(:), pointer :: fdata => null()
  real(fp), dimension(:,:), pointer:: xcubes

  real(fp), dimension(:), allocatable :: x
  real(fp), dimension(:), allocatable :: xnmin,xnmax  

  real(fp) :: f, fmin, fmax

  integer :: i,j
  real(fp) :: lnA, sr1, sr2, sr3
  real(fp) :: lnAmin, lnAmax, sr1min, sr1max, sr2min, sr2max
  real(fp) :: sr3min,sr3max

  real(fp) :: rmax


  logical, parameter :: training = .true.



  call read_binned_posterior('mp_sr2ndlog_litecore120mhi_D_posterior_4D_10.dat',fdata,xdata)

  ndata = size(fdata)
  ndim = size(xdata,1)
  if (size(xdata,2).ne.ndata) stop 'internal error'

!3D 
!  nctrs = 90
!  nfits = 300

!!4D
  nctrs = 90 
  nfits = 300

  print *,'ndata= ',ndata
  print *,'ndim= ',ndim
  print *,'nctrs= nfits= ',nctrs,nfits
  fmax = maxval(fdata)
  fmin = minval(fdata)
  print *,'fdata max',fmax
  print *,'fdata min',fmin

  allocate(xnmin(ndim), xnmax(ndim))
  allocate(xcubes(ndim,ndata))

  call posterior_boundaries(xdata,xnmin,xnmax)

  print *,'xnmin',xnmin
  print *,'xnmax',xnmax
  print *

  call cubize_paramspace(xdata,xcubes)


!  call regularize_zerosphere(fdata,xcubes)

  deallocate(xdata)

  if (training) then

     rmax = shepard_maxradius(ndim,ndata,xcubes,fdata,nctrs,nfits)
     call save_shepdata('shepdata.dat',rmax)
     call save_posterior('postcubed.dat',fdata,xcubes)
     call save_boundaries('bounds.dat',xnmin,xnmax,fmin,fmax)
  else

     deallocate(fdata,xcubes)
     call load_posterior('postcubed.dat',fdata,xcubes)
     call load_shepdata('shepdata.dat',rmax)

  endif

  allocate(x(ndim))
  
  lnAmin = xnmin(1)
  lnAmax = xnmax(1)
  sr1min = xnmin(2)
  sr1max = xnmax(2)
  sr2min = xnmin(3)
  sr2max = xnmax(3)
  sr3min = xnmin(4)
  sr3max = xnmax(4)

  call delete_file('output.dat')
  call delete_file('output2.dat')
  call delete_file('output3.dat')
  lnA = 3.09
  sr3 = 0.0


  do i=1,ndump
     sr1 = sr1min + (sr1max-sr1min)*real(i-1,fp)/real(ndump-1,fp)

     do j=1,ndump
        sr2 = sr2min + (sr2max-sr2min)*real(j-1,fp)/real(ndump-1,fp)
        x(1) = (lnA-xnmin(1))/(xnmax(1)-xnmin(1))
        x(2) = (sr1-xnmin(2))/(xnmax(2)-xnmin(2))
        x(3) = (sr2-xnmin(3))/(xnmax(3)-xnmin(3))
        if (ndim.eq.4) x(4) = (sr3-xnmin(4))/(xnmax(4)-xnmin(4))

        f = shepard_eval(ndim,ndata,xcubes,fdata,rmax,x)

        call livewrite('output.dat',sr1,sr2,f)
!        call livewrite('test.dat',sr1,sr2,rbflike_eval(x))
     enddo
  enddo

  sr3=0
  sr2 = 0.04
  do i=1,ndump
     sr1 = sr1min + (sr1max-sr1min)*real(i-1,fp)/real(ndump-1,fp)
     
     do j=1,ndump
        lnA = lnAmin + (lnAmax-lnAmin)*real(j-1,fp)/real(ndump-1,fp)
        x(1) = (lnA-xnmin(1))/(xnmax(1)-xnmin(1))
        x(2) = (sr1-xnmin(2))/(xnmax(2)-xnmin(2))
        x(3) = (sr2-xnmin(3))/(xnmax(3)-xnmin(3))
        if (ndim.eq.4) x(4) = (sr3-xnmin(4))/(xnmax(4)-xnmin(4))

        f = shepard_eval(ndim,ndata,xcubes,fdata,rmax,x)
        call livewrite('output2.dat',sr1,lnA,f)

     enddo

  enddo

  if (ndim.lt.4) stop

  lnA = 3.09
  sr2 = 0.04
  do i=1,ndump
     sr1 = sr1min + (sr1max-sr1min)*real(i-1,fp)/real(ndump-1,fp)
     
     do j=1,ndump
        sr3 = sr3min + (sr3max-sr3min)*real(j-1,fp)/real(ndump-1,fp)
        x(1) = (lnA-xnmin(1))/(xnmax(1)-xnmin(1))
        x(2) = (sr1-xnmin(2))/(xnmax(2)-xnmin(2))
        x(3) = (sr2-xnmin(3))/(xnmax(3)-xnmin(3))
        if (ndim.eq.4) x(4) = (sr3-xnmin(4))/(xnmax(4)-xnmin(4))

        f = shepard_eval(ndim,ndata,xcubes,fdata,rmax,x)

        call livewrite('output3.dat',sr1,sr3,f)
!        call livewrite('test.dat',sr1,sr2,rbflike_eval(x))
     enddo

  enddo

  
  deallocate(xcubes,fdata)
  
  deallocate(xnmin,xnmax)
  deallocate(x)

contains


  subroutine vecdump(name,vec)
    implicit none
    character(*) :: name
    integer(ip) :: i,j,npts,ncol   
    real(fp), dimension(:,:) :: vec
    
    npts=size(vec,2)
    ncol = size(vec,1)
      
    if (size(vec,2).ne.npts) then
       write(*,*)'inoutfile: WARNING: vectors length differ'
       read(*,*)
    endif
    
    write(*,*)'__write: save in ',name
    open(10,file=name,status='unknown')
    
    
    do j=1,npts      
       write(10,100) (vec(i,j),i=1,ncol)
    enddo
    
    close(10)
    
100 format(50(ES25.16E3))      

  end subroutine vecdump
  

  subroutine regularize_zerosphere(f,xcubes)
    implicit none
    real(fp), dimension(:), intent(inout) :: f
    real(fp), dimension(:,:), intent(in) :: xcubes

    integer(ip) :: ndim
    real(fp) :: distance
    real(fp) :: ZeroData,ZeroCenter,ZeroEdge

    
    integer(ip) :: i,j

    ndim = size(xcubes,1)
    ndata = size(xcubes,2)

    if (size(f,1).ne.ndata) stop 'regularize_zerosphere: wrong size!'

    ZeroData = minval(f)
    ZeroCenter = ZeroData
    ZeroEdge = ZeroCenter - 1._fp

    write(*,*)'regularize_zerosphere:     '
    write(*,*)'renormalizing Zeros=       ',ZeroData
    write(*,*)'to: ZeroCenter= ZeroEdge=  ',ZeroCenter,ZeroEdge
    write(*,*)

    do i=1,ndata
       distance = sum((xcubes(:,i)-0.5_fp)**2)
       if (f(i).eq.ZeroData) then
          f(i) = ZeroCenter + (ZeroEdge-ZeroCenter) &
               *4._fp*distance/real(ndim,fp)
       end if
    enddo

  end subroutine regularize_zerosphere


end program shepmain
