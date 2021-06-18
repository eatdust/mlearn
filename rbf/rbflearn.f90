module rbflearn
  use rbfprec
  use ioml
  use rbfuncs
  implicit none

  private

  type fitshare
     real(fp), dimension(:), allocatable :: paramin, paramax
     real(fp), dimension(:,:), allocatable :: xgrid, xcubes
     real(fp), dimension(:,:), allocatable :: xctrs
     real(fp), dimension(:), allocatable :: weights
     real(fp), dimension(:), allocatable :: fgrid, fdata
  end type fitshare
  
  
!to share data between routines here only
  type(fitshare), save :: fit
!$omp threadprivate(fit)


!  procedure(rbf_func), pointer :: ptr_rbf_func => rbf_polyharmonic_two
!  procedure(rbf_func), pointer :: ptr_rbf_func => rbf_gaussian
   procedure(rbf_func), pointer :: ptr_rbf_func => rbf_inverse_monomial_two

  public :: ptr_rbf_func
  public :: learn_with_fitscale_rbfngp, learn_with_rbfngp, learn_with_rbf


contains


    subroutine learn_with_fitscale_rbfngp(ndim,ndata,igrid,xcubes,fdata)
    use rbfprec, only : nearest_grid_point
    use rbfnd, only : rbf_grid_centers
    use iorbf, only : save_weights, save_centres
    use ioml, only : delete_file, allwrite
    implicit none
    integer, intent(in) :: ndim, ndata
    integer, dimension(ndim), intent(in) :: igrid
    real(fp), dimension(ndim,ndata), intent(in) :: xcubes
    real(fp), dimension(ndata), intent(in) :: fdata

    integer :: ngrid
    real(fp), dimension(:,:), allocatable :: xgrid
    real(fp), dimension(:), allocatable :: fgrid
    integer, dimension(:), allocatable :: fcount

!    real(fp), dimension(ndata) :: compressfdata
    
    integer :: nctrs
    integer, dimension(ndim) :: ictrs

    real(fp) :: scale
    real(fp), dimension(:,:), allocatable :: xctrs
    real(fp), dimension(:), allocatable :: weights

    integer :: i,j
    

    ngrid = product(igrid)

    allocate(xgrid(ndim,ngrid))
    call rbf_grid_centers(xgrid,igrid)
    
    allocate(fgrid(ngrid))
    allocate(fcount(ngrid))

    fcount = 0
    fgrid = 0._fp

    
    do i=1,ndata
       call nearest_grid_point(ndim,igrid,xcubes(:,i),fdata(i),fgrid,fcount)
    enddo

!    call delete_file('datagrid_rbfngp.dat')
!    call allwrite('datagrid_rbfngp.dat',xgrid(1,:),xgrid(2,:),fgrid(:))


!by nature, the fit finds the scale that renders the network closest
!to the data, here the NGP data. Let's pick the same number of control
!points than grid points, for balance between smoothing and
!overfitting
    ictrs = igrid
    nctrs = product(ictrs)

    allocate(xctrs(ndim,nctrs))
    allocate(weights(nctrs))


    call rbf_grid_centers(xctrs,ictrs)

    write(*,*)'learn_with_fit_scale_rbfngp:'
    write(*,*)'ngrid= nctrs= ',ngrid, nctrs
    write(*,*)'searching for best scale...'
    
    call fit_scale_rbfngp(ndim,ngrid,nctrs,ndata,xgrid,fgrid,xctrs,xcubes,fdata,weights,scale)

    write(*,*)'best scale= ',scale
        
    call save_weights('weights.dat',scale,weights)
    call save_centres('centres.dat',xctrs)
    call save_centres('cubes.dat',xcubes)

    deallocate(xgrid,fgrid,fcount)
    deallocate(xctrs,weights)
    
  end subroutine learn_with_fitscale_rbfngp

  
   

  subroutine fit_scale_rbfngp(ndim,ngrid,nctrs,ndata,xgrid,fgrid,xctrs,xcubes,fdata,weights,scale)
    use lm, only : lmdif1, compactify, decompactify
    implicit none
    integer, intent(in) :: ndim, ngrid, nctrs, ndata
    real(fp), dimension(ndim,ngrid), intent(in) :: xgrid
    real(fp), dimension(ngrid), intent(in) :: fgrid
    real(fp), dimension(ndim,nctrs), intent(in) :: xctrs
    real(fp), dimension(ndim,ndata), intent(in) :: xcubes
    real(fp), dimension(ndata), intent(in) :: fdata
    
    real(fp), dimension(nctrs), intent(out) :: weights
    real(fp), intent(out) :: scale

   
    integer :: info

    integer, parameter :: ns = 1
    integer, dimension(ns) :: iwb
    real(fp), dimension(ns) :: s,sinf,smin,smax

    integer :: npoints
    real(fp), dimension(:), allocatable :: rest

    real(fp), parameter :: tol=1d-4

        
    s(1) = 1.0_fp
       
    smin(1) = 0.1_fp
    smax(1) = 2._fp

    sinf = decompactify(s,smin,smax,ns)

    allocate(fit%paramin(ns),fit%paramax(ns))
    allocate(fit%xgrid(ndim,ngrid))
    allocate(fit%xcubes(ndim,ndata))
    allocate(fit%xctrs(ndim,nctrs))
    allocate(fit%fgrid(ngrid))
    allocate(fit%fdata(ndata))
    allocate(fit%weights(nctrs))

    fit%paramin = smin
    fit%paramax = smax
    fit%xgrid = xgrid
    fit%fgrid = fgrid
    fit%xctrs = xctrs
    fit%xcubes = xcubes
    fit%fdata = fdata

    npoints = ndata

    allocate(rest(npoints))
    
    call lmdif1(fcn_scale_rbfngp,ndata,ns,sinf,rest,tol,info,iwb)

    deallocate(rest)
    
    s = compactify(sinf,smin,smax,ns)

    if (info.ge.4) then
       write(*,*)'fit_scale_rbfngp:'
       write(*,*)'info= ',info
       stop
    endif

    weights = fit%weights
    scale = s(1) / real(nctrs,fp)**(1._fp/real(ndim,fp))

    deallocate(fit%paramin,fit%paramax)
    deallocate(fit%xgrid)
    deallocate(fit%xcubes)
    deallocate(fit%xctrs)
    deallocate(fit%fgrid)
    deallocate(fit%fdata)
    deallocate(fit%weights)
    
    
  end subroutine fit_scale_rbfngp



!we minimize the difference between the machine learned ngp and the
!unbinned data (the ngp are just an intermediate smoothing method)
  subroutine fcn_scale_rbfngp(npoints, npar, sinf, fvec, iflag)
    use lm, only : compactify
    use rbfnd, only : rbf_svd_eval, rbf_svd_weights
    implicit none    
    integer, intent(in) :: npoints,npar
    real(fp), dimension(:), intent(in) :: sinf
    real(fp), dimension(:), intent(inout) :: fvec
    integer, intent(inout) :: iflag
   
    integer :: i, ndim, ngrid, nctrs, ndata
    real(fp), dimension(npar) :: s
    real(fp) :: scale
    
    if (npar.ne.1) stop 'fcn_scale_rbfngp: npar is not one!'
    
    s = compactify(sinf,fit%paramin,fit%paramax,npar)
    
    ndim = size(fit%xgrid,1)
    ngrid = size(fit%xgrid,2)
    nctrs = size(fit%weights,1)
    ndata = size(fit%xcubes,2)
               
    scale = s(1)/real(nctrs,fp)**(1._fp/real(ndim,fp))

    call rbf_svd_weights(ndim,ngrid,nctrs,scale,ptr_rbf_func &
         ,fit%xgrid,fit%fgrid,fit%xctrs,fit%weights)

    write(*,*)'test scale= ',scale
    
    do i=1,ndata
       fvec(i) = rbf_svd_eval(ndim,nctrs,scale, ptr_rbf_func &
            ,fit%xctrs,fit%weights,fit%xcubes(:,i)) - fit%fdata(i)
    enddo

    
    
  end subroutine fcn_scale_rbfngp
  



   subroutine learn_with_rbfngp(ndim,ndata,igrid,xcubes,fdata)
    use rbfprec, only : nearest_grid_point
    use rbfnd, only : rbf_grid_centers, rbf_svd_weights
    use iorbf, only : save_weights, save_centres
    use ioml, only : delete_file, allwrite
    implicit none
    integer, intent(in) :: ndim, ndata
    integer, dimension(ndim),intent(in) :: igrid
    real(fp), dimension(ndim,ndata), intent(in) :: xcubes
    real(fp), dimension(ndata), intent(in) :: fdata

    integer :: ngrid
    real(fp), dimension(:,:), allocatable :: xgrid
    real(fp), dimension(:), allocatable :: fgrid
    integer, dimension(:), allocatable :: fcount

    
    integer :: nctrs
    integer, dimension(ndim) :: ictrs

    real(fp) :: scale
    real(fp), dimension(:,:), allocatable :: xctrs
    real(fp), dimension(:), allocatable :: weights

    integer :: i,j
    

    ngrid = product(igrid)

    allocate(xgrid(ndim,ngrid))
    call rbf_grid_centers(xgrid,igrid)
    
    allocate(fgrid(ngrid))
    allocate(fcount(ngrid))

    fcount = 0
    fgrid = 0._fp

    do i=1,ndata
       call nearest_grid_point(ndim,igrid,xcubes(:,i),fdata(i),fgrid,fcount)
    enddo

    call delete_file('datagrid_rbfngp.dat')
    call allwrite('datagrid_rbfngp.dat',xgrid(1,:),xgrid(2,:),fgrid(:))

    ictrs = igrid
    nctrs = product(ictrs)
    allocate(xctrs(ndim,nctrs))
    allocate(weights(nctrs))


    call rbf_grid_centers(xctrs,ictrs)


    scale = 0.9_fp
    scale = scale/real(nctrs,fp)**(1._fp/real(ndim,fp))
    
    write(*,*)'learn_with_rbfngp:'
    write(*,*)'ngrid= nctrs= ',ngrid, nctrs
    write(*,*)'scale= ',scale    
    
    call rbf_svd_weights(ndim,ngrid,nctrs,scale,ptr_rbf_func &
         ,xgrid,fgrid,xctrs,weights)
    
    call save_weights('weights.dat',scale,weights)
    call save_centres('centres.dat',xctrs)
    call save_centres('cubes.dat',xcubes)

    deallocate(xgrid,fgrid,fcount)
    deallocate(xctrs,weights)
    
  end subroutine learn_with_rbfngp


  



  
  subroutine learn_with_rbf(ndim,ndata,ictrs,xcubes,fdata)
    use rbfnd, only : rbf_grid_centers, rbf_svd_weights
    use iorbf, only : save_weights, save_centres
    implicit none
    integer, intent(in) :: ndim, ndata
    integer, dimension(ndim), intent(in) :: ictrs
    real(fp), dimension(ndim,ndata), intent(in) :: xcubes
    real(fp), dimension(ndata), intent(in) :: fdata

    integer :: nctrs


    real(fp) :: scale
    real(fp), dimension(:,:), allocatable :: xctrs
    real(fp), dimension(:), allocatable :: weights

    integer :: i,j

    nctrs = product(ictrs)
    allocate(xctrs(ndim,nctrs))
    
    call rbf_grid_centers(xctrs,ictrs)


    
    write(*,*)'learn_with_rbf:'
    write(*,*)'nctrs= ',nctrs
    
    allocate(weights(nctrs))


    scale = 1.10_fp/real(nctrs,fp)**(1._fp/real(ndim,fp))


    call rbf_svd_weights(ndim,ndata,nctrs,scale,ptr_rbf_func &
         ,xcubes,fdata,xctrs,weights)

    call save_weights('weights.dat',scale,weights)
    call save_centres('centres.dat',xctrs)
    call save_centres('cubes.dat',xcubes)
    

    deallocate(xctrs,weights)
    
  end subroutine learn_with_rbf
  




end module rbflearn
