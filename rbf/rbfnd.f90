module rbfnd  
  use rbfprec
  implicit none

  logical, parameter :: display = .true.

  private

  public rbf_random_centers, rbf_grid_centers, rbf_svd_weights, rbf_svd_eval
  public rbf_multiquad, rbf_inverse_multiquad, rbf_thin_plate, rbf_gaussian
  public rbf_polyharmonic_one, rbf_polyharmonic_three, rbf_polyharmonic_five
  public rbf_polyharmonic_two, rbf_polyharmonic_four, rbf_polyharmonic_six

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


  subroutine rbf_grid_centers(xctrs,ictrs)
    use rbfprec, only : unflatten_indices
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
          stop 'rbf_grid_centers: cannot equally grid!'
       endif
       isize = npts
    endif

    do q=1,nctrs
       ivec = unflatten_indices(ndim,isize,q)
       xctrs(:,q) = real(ivec-1,fp)/real(isize-1,fp)
    enddo
    
    deallocate(isize, ivec)

  end subroutine rbf_grid_centers


  


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
       write(*,*)
       write(*,*)'solve_svd: calling lapack SVD...'
    endif

    call dgesdd(jobz, m, n, a, lda, sdiag, u, ldu, v, ldv, work, lwork, &
         iwork, info )

    if (display) write(*,*)'solve_svd: svd done!'

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


  
  subroutine rbf_multiquad ( n, r, r0, v )
    implicit none

    integer ( ip ) n

    real ( fp ) r(n)
    real ( fp ) r0
    real ( fp ) v(n)

    v(1:n) = sqrt ( r(1:n)**2 + r0**2 )
  end subroutine rbf_multiquad


  subroutine rbf_inverse_multiquad ( n, r, r0, v )
       implicit none

    integer ( ip ) n

    real ( fp ) r(n)
    real ( fp ) r0
    real ( fp ) v(n)

    v = 1.0D+00 / sqrt ( r**2 + r0**2 )

  end subroutine rbf_inverse_multiquad


  subroutine rbf_thin_plate ( n, r, r0, v )
    implicit none

    integer ( ip ) n  
    integer ( ip ) i
    real ( fp ) r(n)
    real ( fp ) r0
    real ( fp ) v(n)


    do i = 1, n
       if ( r(i) .le. 0.0D+00 ) then
          v(i) = 0.0D+00       
       else
          v(i) = r(i)**2 * log ( r(i) / r0 )
       end if
    end do
     
  end subroutine rbf_thin_plate


  subroutine rbf_gaussian ( n, r, r0, v )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)

    v(1:n) = exp ( - 0.5D+00 * r(1:n)**2 / r0**2 )

  end subroutine rbf_gaussian




  subroutine rbf_polyharmonic_log ( n, r, r0, v, k )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0, ru
    real(fp) :: v(n)
    integer(ip) :: k,i

    do i = 1, n
       ru = r(i)/r0
       if ( ru .le. 1._fp ) then
          v(i) = r(i)**(k-1) * log(ru**r(i))
       else
          v(i) = r(i)**k * log ( ru )
       end if
    end do

  end subroutine rbf_polyharmonic_log


 

  subroutine rbf_polyharmonic ( n, r, r0, v, k )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)
    integer(ip) :: k

    v(1:n) = (r(1:n)/r0)**k

  end subroutine rbf_polyharmonic



  subroutine rbf_polyharmonic_one ( n, r, r0, v )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)

    call rbf_polyharmonic(n, r, r0, v,1_ip)

  end subroutine rbf_polyharmonic_one


  subroutine rbf_polyharmonic_two(n, r, r0, v )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)

    call rbf_polyharmonic_log(n, r, r0, v,2_ip)

  end subroutine rbf_polyharmonic_two



  subroutine rbf_polyharmonic_three ( n, r, r0, v )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)

    call rbf_polyharmonic(n, r, r0, v, 3_ip)

  end subroutine rbf_polyharmonic_three


  subroutine rbf_polyharmonic_four(n, r, r0, v )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)

    call rbf_polyharmonic_log(n, r, r0, v,4_ip)

  end subroutine rbf_polyharmonic_four



  subroutine rbf_polyharmonic_five ( n, r, r0, v )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)

    call rbf_polyharmonic(n, r, r0, v, 5_ip)
    
  end subroutine rbf_polyharmonic_five


  subroutine rbf_polyharmonic_six(n, r, r0, v )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)

    call rbf_polyharmonic_log(n, r, r0, v,6_ip)

  end subroutine rbf_polyharmonic_six

  
 

end module rbfnd
