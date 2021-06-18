program rbfmain
  use rbfprec
  use rbfnd
  use rbflearn, only : learn_with_fitscale_rbfngp, learn_with_rbfngp,learn_with_rbf
  use rbflearn, only : ptr_rbf_func
  use ioml
  use iocmb
  use iorbf
  implicit none
 
  integer :: ndim
  integer(ip) :: ndata
  
  integer(ip), dimension(:), allocatable, save :: ictrs
  integer(ip), save :: nctrs = 0

  integer, parameter :: ndump = 100

  real(fp), dimension(:,:), pointer :: xdata => null()
  real(fp), dimension(:), pointer :: fdata => null()

  real(fp), save, dimension(:,:), pointer :: xctrs => null()
  real(fp), save, dimension(:), pointer :: weights => null()
  real(fp) :: scale, sigma

  real(fp), dimension(:), allocatable :: x
  real(fp), dimension(:), allocatable :: xnmin,xnmax  
  real(fp), dimension(:,:), allocatable :: xcubes
  real(fp) :: f, fmin, fmax

  integer :: i,j
  real(fp) :: lnA, sr1, sr2, sr3
  real(fp) :: lnAmin, lnAmax, sr1min, sr1max, sr2min, sr2max
  real(fp) :: sr3min,sr3max


  logical, parameter :: training = .true.


  call read_binned_posterior('test_posterior.dat',fdata,xdata)

  ndata = size(fdata)
  ndim = size(xdata,1)
  if (size(xdata,2).ne.ndata) stop 'internal error'

  allocate(ictrs(ndim))
  
!  ictrs = (/4,4,4,11/)
!  ictrs = (/12,18,12,12/)
!  ictrs = (/12,20,16/)
  ictrs = (/6,6,6,6/)
  
  nctrs = product(ictrs)
!  nctrs = 1400

  print *,'ndata= ',ndata
  print *,'ndim= ',ndim
  print *,'nctrs= ',nctrs
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

     call learn_with_fitscale_rbfngp(ndim,ndata,ictrs,xcubes,fdata)

     call save_boundaries('bounds.dat',xnmin,xnmax,fmin,fmax)

     print *,'done!'

  endif

  call load_weights('weights.dat',scale,weights)
  call load_centres('centres.dat',xctrs)
  if (size(weights,1).ne.size(xctrs,2)) stop 'weights/centres mismatch!'
  nctrs = size(weights,1)
  
  allocate(x(ndim))
  
  lnAmin = xnmin(1)
  lnAmax = xnmax(1)
  sr1min = xnmin(2)
  sr1max = xnmax(2)
  sr2min = xnmin(3)
  sr2max = xnmax(3)
  if (ndim.gt.3) then
     sr3min = xnmin(4)
     sr3max = xnmax(4)
  endif
     
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
        x(4) = (sr3-xnmin(4))/(xnmax(4)-xnmin(4))

        f = rbf_svd_eval(ndim,nctrs,scale,ptr_rbf_func,xctrs,weights,x)
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
        x(4) = (sr3-xnmin(4))/(xnmax(4)-xnmin(4))

        f = rbf_svd_eval(ndim,nctrs,scale,ptr_rbf_func,xctrs,weights,x)
        call livewrite('output2.dat',sr1,lnA,f)
!        call livewrite('test.dat',sr1,sr2,rbflike_eval(x))
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
        x(4) = (sr3-xnmin(4))/(xnmax(4)-xnmin(4))

        f = rbf_svd_eval(ndim,nctrs,scale,ptr_rbf_func,xctrs,weights,x)
!        call livewrite('output3.dat',sr2,sr3,f)
        call livewrite('output3.dat',sr1,sr3,f)
!        call livewrite('test.dat',sr1,sr2,rbflike_eval(x))
     enddo

  enddo

  
  deallocate(xcubes,fdata)
  deallocate(weights)
  
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


end program rbfmain
