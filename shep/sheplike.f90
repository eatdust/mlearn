module sheplike
  use shepprec
  use qshepmdata
  implicit none


  private
   
 
  integer(ip), save :: ndim, ndata
  real(fp), save :: rmax
  real(fp), save, dimension(:,:), pointer :: xcubes => null()
  real(fp), save, dimension(:), pointer :: fdata => null()
  real(fp), save, dimension(:), pointer :: xpmin => null(), xpmax=>null()
  real(fp), save :: fmin = huge(1._fp)
  real(fp), save :: fmax = -huge(1._fp)

  logical, parameter :: display = .false.

  public initialize_shep_like, free_shep
  public sheplike_eval, uncubize_shepparams, cubize_shepparams
  public check_shep, get_shep_ndim, get_shep_ndata
  public cutmin_shepparams, cutmax_shepparams
  public get_shep_xpmin, get_shep_xpmax
  public get_shep_fmin, get_shep_fmax

contains


  function check_shep()
    implicit none
    logical :: check_shep
    
    check_shep = associated(xcubes) .and. associated(fdata) &
         .and. associated(xpmin) .and. associated(xpmax)

  end function check_shep


  subroutine free_shep()
    implicit none
    

    if (.not.check_shep()) stop 'free_shep: data not allocated!'

    if (associated(xcubes)) deallocate(xcubes)
    xcubes => null()
    if (associated(fdata)) deallocate(fdata)
    fdata => null()
    if (associated(xpmin)) deallocate(xpmin)
    xpmin => null()
    if (associated(xpmin)) deallocate(xpmax)
    xpmax => null()

  end subroutine free_shep


  function get_shep_ndim()
    implicit none
    integer :: get_shep_ndim

    if (.not.check_shep()) then
       get_shep_ndim = 0
       return
    endif

    get_shep_ndim = ndim

  end function get_shep_ndim


  function get_shep_ndata()
    implicit none
    integer :: get_shep_ndata

    if (.not.check_shep()) then
       get_shep_ndata = 0
       return
    endif

    get_shep_ndata = ndata

  end function get_shep_ndata


  function get_shep_xpmin(i)
    implicit none
    integer, intent(in) :: i
    real(fp) :: get_shep_xpmin

    if (.not.associated(xpmin)) then
       stop 'get_shep_xpmin: not allocated!'
    elseif ((i.gt.ndim).or.(i.le.0)) then
       stop 'get_shep_xpmin: i out of range!'
    endif

    get_shep_xpmin = xpmin(i)

  end function get_shep_xpmin


  function get_shep_xpmax(i)
    implicit none
    integer, intent(in) :: i
    real(fp) :: get_shep_xpmax

    if (.not.associated(xpmax)) then
       stop 'get_shep_xpmax: not allocated!'
    elseif ((i.gt.ndim).or.(i.le.0)) then
       stop 'get_shep_xpmax: i out of range!'
    endif

    get_shep_xpmax = xpmax(i)

  end function get_shep_xpmax


  function get_shep_fmin()
    implicit none
    real(fp) :: get_shep_fmin

    if (fmin.eq.huge(1._fp)) then
       stop 'get_shep_fmin: not initialized!'
    endif
    
    get_shep_fmin = fmin

  end function get_shep_fmin


  function get_shep_fmax()
    implicit none
    real(fp) :: get_shep_fmax

    if (fmax.eq.-huge(1._fp)) then
       stop 'get_shep_fmax: not initialized!'
    endif
    
    get_shep_fmax = fmax

  end function get_shep_fmax



  subroutine initialize_shep_like(fileshep, filepost, filebounds)
    use ioshep, only : load_shepdata, load_posterior
    use iond, only :  read_boundaries
    implicit none   
    character(len=*), intent(in) :: fileshep, filepost, filebounds
    
    call load_shepdata(fileshep,rmax)
    call load_posterior(filepost,fdata,xcubes)
    call read_boundaries(filebounds,xpmin,xpmax,fmin,fmax)

    if (size(fdata,1).ne.size(xcubes,2)) then
       stop 'initialize_shep_like: size mismatch!'
    endif

    ndim = size(xcubes,1)
    ndata = size(fdata,1)

  end subroutine initialize_shep_like



  function sheplike_eval(x)
    use shepnd, only : shepard_eval
    implicit none
    real(fp) :: sheplike_eval
    real(fp), dimension(:), intent(in) :: x
    if (size(x,1).ne.ndim) stop 'sheplike_eval: x dim screwed!'

    if (any(x.gt.1._fp).or.any(x.lt.0._fp)) stop 'sheplike_eval: uncubed input!'

    sheplike_eval = shepard_eval(ndim,ndata,xcubes,fdata,rmax,x)
    
  end function sheplike_eval




  function cutmin_shepparams(ndim,icut,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim,icut
    real(fp), dimension(ndim) :: cutmin_shepparams
    real(fp), dimension(ndim), intent(in) :: uncubed

    integer(ip) :: i

    if (icut.gt.ndim) stop 'cutmin_shepparams: icut > ndim!'

    cutmin_shepparams = uncubed

    if (uncubed(icut).lt.xpmin(icut)) then
       cutmin_shepparams(icut) = xpmin(icut)
       if (display) then
          write(*,*)'cutmin_shepparams: ',icut,uncubed(icut),xpmin(icut)
       endif
    endif    

  end function cutmin_shepparams



  function cutmax_shepparams(ndim,icut,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim,icut
    real(fp), dimension(ndim) :: cutmax_shepparams
    real(fp), dimension(ndim), intent(in) :: uncubed

    integer(ip) :: i

    if (icut.gt.ndim) stop 'cutmax_shepparams: icut > ndim!'

    cutmax_shepparams = uncubed

    if (uncubed(icut).gt.xpmax(icut)) then
       cutmax_shepparams(icut) = xpmax(icut)
       if (display) then
          write(*,*)'cutmax_shepparams: ',icut,uncubed(icut),xpmax(icut)
       endif
    endif    

  end function cutmax_shepparams




  function cubize_shepparams(ndim,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim
    real(fp), dimension(ndim) :: cubize_shepparams
    real(fp), dimension(ndim), intent(in) :: uncubed

    if (ndim.ne.size(uncubed,1)) then
       stop 'uncubize_shepparams: sizes do not match!'
    endif

    cubize_shepparams = (uncubed - xpmin)/(xpmax-xpmin)

  end function cubize_shepparams


  function uncubize_shepparams(ndim,cubed)
    implicit none
    integer(ip), intent(in) :: ndim
    real(fp), dimension(ndim) :: uncubize_shepparams
    real(fp), dimension(ndim), intent(in) :: cubed

    if (ndim.ne.size(cubed,1)) then
       stop 'uncubize_shepparams: sizes do not match!'
    endif

    uncubize_shepparams = xpmin + (xpmax-xpmin)*cubed

  end function uncubize_shepparams
 


end module sheplike
