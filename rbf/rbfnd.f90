module rbfnd  
  use rbfprec
  use rbfuncs
  implicit none

  logical, parameter :: display = .true.

  private

  public rbf_random_centers, rbf_grid_centers, rbf_subgrid_centers
  public rbf_datagrid_centers

  public rbf_svd_weights, rbf_svd_eval
  public rbf_multiquad, rbf_inverse_multiquad, rbf_thin_plate, rbf_gaussian
  public rbf_polyharmonic_one, rbf_polyharmonic_three, rbf_polyharmonic_five
  public rbf_polyharmonic_two, rbf_polyharmonic_four, rbf_polyharmonic_six
  public rbf_polyharmonic_log_two, rbf_polyharmonic_log_four
  
contains


!returns a set of random position in the ndim-dimensional unit cube
  subroutine rbf_random_centers(xctrs)
    implicit none
    real(fp), dimension(:,:), intent(out) :: xctrs

    integer(ip), parameter :: iseed = 14
    logical, save :: seeded = .false.

    if (.not.seeded) then
       call cte_random_seed(iseed)
       seeded = .true.
    endif
    
    call random_number(xctrs)

  end subroutine rbf_random_centers

  
  
!return of set of grid points centered and equally spaced in all dimensions
  subroutine rbf_grid_centers(xctrs,ictrs)
    implicit none
    real(fp), dimension(:,:), intent(out) :: xctrs
    integer(ip), dimension(:), intent(in), optional :: ictrs

    integer(ip) :: ndim, nctrs
    integer(ip) :: npts, q
    integer(ip), dimension(:), allocatable :: isize, ivec

    ndim = size(xctrs,1)
    nctrs = size(xctrs,2)
    

    allocate(isize(ndim), ivec(ndim))

    if (present(ictrs)) then
       if (size(ictrs,1).ne.ndim) stop 'rbf_grid_centers: ictrs mismatch!'
       isize = ictrs
       write(*,*)'rbf_grid_centers:'
       write(*,*)'ictrs= ',ictrs
    else
       npts = nint(real(nctrs,fp)**(1._fp/real(ndim,fp)))

       if (npts**ndim.ne.nctrs) then
          write(*,*)'nctrs= ndim= npts= ',nctrs,ndim,npts
          stop 'rbf_grid_centers: cannot get equally spaced grid!'
       endif
       isize = npts
    endif

    do q=1,nctrs
       ivec = unflatten_indices(ndim,isize,q)
       xctrs(:,q) = real(ivec-1,fp)/real(isize-1,fp)
    enddo
    
    deallocate(isize, ivec)

  end subroutine rbf_grid_centers



!compress the grid to a subgrid given by xmin xmax [0,1]  
  subroutine rbf_subgrid_centers(xctrs,ictrs,xmin,xmax)
    use rbfprec, only : unflatten_indices
    implicit none
    real(fp), dimension(:,:), intent(out) :: xctrs
    integer(ip), dimension(:), intent(in), optional :: ictrs
    real(fp), dimension(:), intent(in) :: xmin, xmax
    
    integer(ip) :: ndim,i

    ndim = size(xctrs,1)

    if ((size(xmin,1).ne.ndim).or.(size(xmax,1).ne.ndim)) then
       stop 'rbf_subgrid_centers: inconsistent input dimensions'
    endif
    
    if ((any(xmin.gt.1._fp)) .or. (any(xmin.lt.0._fp))) then
       write(*,*)'rbf_subgrid_centers:'
       stop 'xmin not in [0,1]'
    endif

    if ((any(xmax.gt.1._fp)) .or. (any(xmax.lt.0._fp))) then
       write(*,*)'rbf_subgrid_centers:'
       stop 'xmmax not in [0,1]'
    endif
    
    call rbf_grid_centers(xctrs,ictrs)

    do i=1,ndim
       xctrs(i,:) = xmin(i) + (xmax(i)-xmin(i))*xctrs(i,:)
    enddo
    

  end subroutine rbf_subgrid_centers

  

   
!Cut of set of grid points equally spaced in all dimensions by
!discarding region set in the function rbf_datagrid_discard()
  subroutine rbf_datagrid_centers(xctrs,ictrs,xcubes,fdata,fcut)
    use ioml, only : allwrite, delete_file
    implicit none
    real(fp), dimension(:,:), allocatable, intent(out) :: xctrs
    integer(ip), dimension(:), intent(in) :: ictrs
    real(fp), dimension(:,:), intent(in) :: xcubes
    real(fp), dimension(:), intent(in) :: fdata

    real(fp), intent(in), optional :: fcut
    
    integer(ip) :: ndim, ndata, nctrs
    integer(ip) :: nskip, i, q, count

    real(fp) :: fskip, fval
    real(fp), dimension(:,:), allocatable :: xgrid
    real(fp), dimension(:), allocatable :: fgrid
    integer(ip), dimension(:), allocatable :: fcount

    logical, parameter :: dumpngp = .true.
    
        
    ndim = size(xcubes,1)
    ndata = size(xcubes,2)

    if ((size(fdata,1).ne.ndata).or.(size(ictrs,1).ne.ndim)) then
       stop 'rbf_data_centers: array size pb!'
    endif

    if (present(fcut)) then
       fskip = fcut
    else
       fskip = 0._fp
    endif
   
    

    allocate(fgrid(product(ictrs)))
    allocate(fcount(product(ictrs)))
    
    fcount = 0
    fgrid = 0._fp

    
    do i=1,ndata
       call nearest_grid_point(ndim,ictrs,xcubes(:,i),fdata(i),fgrid,fcount)
    enddo


    allocate(xgrid(ndim,product(ictrs)))
    call rbf_grid_centers(xgrid,ictrs)

    
    count = 0
    do q=1,product(ictrs)
       fval = 1._fp/sqrt(real(fcount(q),fp))
       if (rbf_datagrid_discard(fval,fskip)) cycle
       count = count + 1
    enddo

    
    nctrs = count

    allocate(xctrs(ndim,nctrs))
    

    count = 0
    do q=1,product(ictrs)
       fval = 1._fp/sqrt(real(fcount(q),fp))
       if (rbf_datagrid_discard(fval,fskip)) cycle
       count = count + 1
       if (count.gt.nctrs) stop 'rbf_datagrid_centers: count screwed!'
       xctrs(1:ndim,count) = xgrid(1:ndim,q)
    enddo

    if (dumpngp) then
       call delete_file('datagrid.dat')
       call allwrite('datagrid.dat',xgrid(1,:),xgrid(2,:),fgrid(:))
    endif

    deallocate(fcount)
    deallocate(fgrid)
    deallocate(xgrid)

  end subroutine rbf_datagrid_centers


  

  function rbf_datagrid_discard(fval,fskip)
    implicit none
    logical :: rbf_datagrid_discard
    real(fp), intent(in) :: fval, fskip
    
    rbf_datagrid_discard = .false.

    rbf_datagrid_discard = (fval.gt.fskip) 

  end function rbf_datagrid_discard
  
  

  

  subroutine clock_random_seed()
    integer(ip) :: i, n, clock
    integer(ip), dimension(:), allocatable :: seed
    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(PUT = seed)
    deallocate(seed)
  end subroutine clock_random_seed


  subroutine cte_random_seed(seedcte)
    integer(ip) :: i, n, seedcte
    integer(ip), dimension(:), allocatable :: seed
    call random_seed(size = n)
    allocate(seed(n))
    
    seed = seedcte  + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
  end subroutine cte_random_seed



  subroutine rbf_svd_weights(ndim,ndata,nctrs,scale,phi,xdata,fdata &
       ,xctrs,weights)
    implicit none

    integer(ip), intent(in) :: ndim, ndata, nctrs
    real(fp), intent(in) :: scale
    real(fp), dimension(ndim,ndata), intent(in) :: xdata
    real(fp), dimension(ndata), intent(in) :: fdata
    real(fp), dimension(ndim,nctrs), intent(in) :: xctrs
    real(fp), dimension(nctrs), intent(out) :: weights
    
    integer(ip) :: i,j

    real(fp), dimension(nctrs) :: r,v
    real(fp), dimension(ndata,nctrs) :: a

    include 'rbfnd.h'

    do i = 1, ndata
       do j = 1,nctrs
          r(j) = sqrt ( sum ( ( xdata(1:ndim,i) - xctrs(1:ndim,j) )**2 ) )
       end do
       call phi(nctrs,r,scale,v)
       a(i,1:nctrs) = v(1:nctrs)
    end do

    call solve_svd (ndata,nctrs,a,fdata,weights)

  end subroutine rbf_svd_weights



  function rbf_svd_eval(ndim,nctrs,scale,phi,xctrs,weights,x)
    implicit none
    real(fp) :: rbf_svd_eval
    integer(ip), intent(in) :: ndim, nctrs
    real(fp), dimension(ndim,nctrs), intent(in) :: xctrs
    real(fp), intent(in) :: scale
    real(fp), dimension(nctrs), intent(in) :: weights
    real(fp), dimension(ndim), intent(in) :: x
    
    real(fp), dimension(nctrs) :: r,v
    integer(ip) :: j

    include 'rbfnd.h'

    do j = 1, nctrs
       r(j) = sqrt ( sum ( ( x(1:ndim) - xctrs(1:ndim,j) )**2 ) )
    enddo

    call phi(nctrs,r,scale,v)

    rbf_svd_eval  = dot_product(v,weights)

  end function rbf_svd_eval
  


  subroutine solve_svd ( m, n, a, b, x )

    !*****************************************************************************80
    !
    !! solves a linear system A*x=b using the SVD.
    !
    !  Discussion:
    !
    !    When the system is determined, the solution is the solution in the
    !    ordinary sense, and A*x = b.
    !
    !    When the system is overdetermined, the solution minimizes the
    !    L2 norm of the residual ||A*x-b||.
    !
    !    When the system is underdetermined, ||A*x-b|| should be zero, and
    !    the solution is the solution of minimum L2 norm, ||x||.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 June 2012
    !
    !  From:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( ip ) M, N, the number of rows and columns
    !    in the matrix A.
    !
    !    Input, real ( fp ) A(M,N), the matrix.
    !
    !    Input, real ( fp ) B(M), the right hand side.
    !
    !    Output, real ( fp ) X(N), the solution.
    !
    implicit none

    integer :: m
    integer :: n

    real ( fp ) a(m,n)
!    real ( fp ) a_pseudo(n,m)
    real ( fp ) b(m)
    integer ( ip ) i,j
    integer(ip) :: info
    integer(ip) :: lda
    integer(ip):: ldu
    integer(ip) :: ldv
    integer(ip) :: lwork, liwork
    character jobu, jobv, jobz
    real ( fp ) s(m,n)
    real ( fp ) sp(n,m)
    real ( fp ) sdiag(min(m,n))
    real ( fp ) u(m,m)
    real ( fp ) v(n,n)
    real ( fp), dimension(:), allocatable :: work
    integer, dimension(:), allocatable :: iwork
    real ( fp ) x(n)
  
    jobu = 'A'
    jobv = 'A'
    jobz = 'A'

    lda = m
    ldu = m
    ldv = n


!for A
!the man is screwed    lwork = 3*min(m,n)*min(m,n) + max(max(m,n),4*min(m,n)*min(m,n)+4*min(m,n))
    lwork=  4*MIN(m,n)*MIN(m,n) + MAX(m,n) + 9*MIN(m,n)
    liwork = 8*min(m,n)


    allocate (work(1:lwork))
    allocate(iwork(1:liwork))

!    call dgesvd ( jobu, jobv, m, n, a, lda, sdiag, u, ldu, v, ldv, work, &
!    lwork, info )

    if (display) then
       write(*,*)'solve_svd: calling lapack SVD...'
    endif

    call dgesdd(jobz, m, n, a, lda, sdiag, u, ldu, v, ldv, work, lwork, &
         iwork, info )

    if (display) then
       write(*,*)'solve_svd: svd done!'
       write(*,*)
    endif
    
    v = transpose(v)

    deallocate(work)
    deallocate(iwork)

    if ( info /= 0 ) then
       write(*,*)'info= ',info
       stop 'solve_svd: lapack error '
    end if

    forall(i=1:n,j=1:m)
       s(j,i) = 0._fp
       sp(i,j) = 0._fp
    end forall

    do i = 1, min ( m, n )
       s(i,i) = sdiag(i)
    end do

!    sp(1:n,1:m) = 0.0D+00
    do i = 1, min ( m, n )
       if ( s(i,i) /= 0.0D+00 ) then
          sp(i,i) = 1.0D+00 / s(i,i)
       end if
    end do

!    a_pseudo(1:n,1:m) = matmul ( v(1:n,1:n), &
!         matmul ( sp(1:n,1:m), transpose ( u(1:m,1:m) ) ) )
  
!    x(1:n) = matmul ( a_pseudo(1:n,1:m), b(1:m) )

    x = matmul( matmul ( v, &
         matmul ( sp, transpose ( u ) ) ), b )

  end subroutine solve_svd



end module rbfnd
