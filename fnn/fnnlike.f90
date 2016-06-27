module fnnlike
  use fnnprec
  use fann
  implicit none


  private
   
 
  integer(ip), save :: ndim, nout

  real(fp), save, dimension(:), pointer :: xpmin => null(), xpmax=>null()
  real(fp), save :: fmin = huge(1._fp), fmax = -huge(1._fp)

  type(C_PTR) :: ann


  logical, parameter :: display = .true.

!if the neural has been trained on the ln(like), or directly on the
!like
  logical, parameter :: loglearned = .false.

  public initialize_fnn_like, free_fnn
  public fnnlike_eval, uncubize_fnnparams, cubize_fnnparams
  public check_fnn, get_fnn_ndim, get_fnn_nout
  public cutmin_fnnparams, cutmax_fnnparams
  public get_fnn_xpmin, get_fnn_xpmax
  public get_fnn_fmin, get_fnn_fmax

contains


  function check_fnn()
    implicit none
    logical :: check_fnn
    
    check_fnn = C_ASSOCIATED(ann) .and. associated(xpmin) .and. associated(xpmax)

  end function check_fnn


  subroutine free_fnn()
    implicit none    
    if (.not.check_fnn()) stop 'free_fnn: data not allocated!'

    if (associated(xpmin)) deallocate(xpmin)
    xpmin => null()
    if (associated(xpmin)) deallocate(xpmax)
    xpmax => null()
    if (C_ASSOCIATED(ann)) call fann_destroy(ann)

  end subroutine free_fnn



  function get_fnn_nout()
    implicit none
    integer :: get_fnn_nout

    if (.not.check_fnn()) then
       get_fnn_nout = 0
       return
    endif

    get_fnn_nout = nout

  end function get_fnn_nout



  function get_fnn_ndim()
    implicit none
    integer :: get_fnn_ndim

    if (.not.check_fnn()) then
       get_fnn_ndim = 0
       return
    endif

    get_fnn_ndim = ndim

  end function get_fnn_ndim




  function get_fnn_xpmin(i)
    implicit none
    integer, intent(in) :: i
    real(fp) :: get_fnn_xpmin

    if (.not.associated(xpmin)) then
       stop 'get_fnn_xpmin: not allocated!'
    elseif ((i.gt.ndim).or.(i.le.0)) then
       stop 'get_fnn_xpmin: i out of range!'
    endif

    get_fnn_xpmin = xpmin(i)

  end function get_fnn_xpmin


  function get_fnn_xpmax(i)
    implicit none
    integer, intent(in) :: i
    real(fp) :: get_fnn_xpmax

    if (.not.associated(xpmax)) then
       stop 'get_fnn_xpmax: not allocated!'
    elseif ((i.gt.ndim).or.(i.le.0)) then
       stop 'get_fnn_xpmax: i out of range!'
    endif

    get_fnn_xpmax = xpmax(i)

  end function get_fnn_xpmax


  function get_fnn_fmin()
    implicit none
    real(fp) :: get_fnn_fmin

    if (fmin.eq.huge(1._fp)) then
       stop 'get_fnn_fmin: not initialized!'
    endif
    
    if (loglearned) then
       get_fnn_fmin = fmin
    else
       get_fnn_fmin = log(fmin)
    endif

  end function get_fnn_fmin


  function get_fnn_fmax()
    implicit none
    real(fp) :: get_fnn_fmax

    if (fmax.eq.-huge(1._fp)) then
       stop 'get_fnn_fmax: not initialized!'
    endif
    
    if (loglearned) then
       get_fnn_fmax = fmax
    else
       get_fnn_fmax = log(fmax)
    endif
    
  end function get_fnn_fmax



  subroutine initialize_fnn_like(filefnn, filebounds)   
    use ioml, only :  read_boundaries
    implicit none   
    character(len=*), intent(in) :: filefnn, filebounds
    
    ann = fann_create_from_file(f_c_string(filefnn))

    call read_boundaries(filebounds,xpmin,xpmax,fmin,fmax)
    
    ndim = int(fann_get_num_input(ann),ip)
    nout = int(fann_get_num_output(ann),ip)

    write(*,*)
    call fann_print_connections(ann)
    write(*,*)


  end subroutine initialize_fnn_like



  function fnnlike_eval(x)
    implicit none
    real(fp) :: fnnlike_eval
    real(fp), dimension(nout) :: feval
    real(fp), dimension(:), intent(in) :: x

    if (size(x,1).ne.ndim) stop 'fnnlike_eval: x dim screwed!'

    if (any(x.gt.1._fp).or.any(x.lt.0._fp)) stop 'fnnlike_eval: uncubed input!'

    feval = real(f_fann_run(ann,real(x,FANN_TYPE)),fp)

    if (loglearned) then
       fnnlike_eval = feval(1) * (fmax-fmin) + fmin
    else
       fnnlike_eval = log( feval(1) * (fmax-fmin) + fmin )
    endif

  end function fnnlike_eval




  function cutmin_fnnparams(ndim,icut,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim,icut
    real(fp), dimension(ndim) :: cutmin_fnnparams
    real(fp), dimension(ndim), intent(in) :: uncubed

    integer(ip) :: i

    if (icut.gt.ndim) stop 'cutmin_fnnparams: icut > ndim!'

    cutmin_fnnparams = uncubed

    if (uncubed(icut).lt.xpmin(icut)) then
       cutmin_fnnparams(icut) = xpmin(icut)
       if (display) then
          write(*,*)'cutmin_fnnparams: ',icut,uncubed(icut),xpmin(icut)
       endif
    endif    

  end function cutmin_fnnparams



  function cutmax_fnnparams(ndim,icut,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim,icut
    real(fp), dimension(ndim) :: cutmax_fnnparams
    real(fp), dimension(ndim), intent(in) :: uncubed

    integer(ip) :: i

    if (icut.gt.ndim) stop 'cutmax_fnnparams: icut > ndim!'

    cutmax_fnnparams = uncubed

    if (uncubed(icut).gt.xpmax(icut)) then
       cutmax_fnnparams(icut) = xpmax(icut)
       if (display) then
          write(*,*)'cutmax_fnnparams: ',icut,uncubed(icut),xpmax(icut)
       endif
    endif    

  end function cutmax_fnnparams




  function cubize_fnnparams(ndim,uncubed)
    implicit none
    integer(ip), intent(in) :: ndim
    real(fp), dimension(ndim) :: cubize_fnnparams
    real(fp), dimension(ndim), intent(in) :: uncubed

    if (ndim.ne.size(uncubed,1)) then
       stop 'uncubize_fnnparams: sizes do not match!'
    endif

    cubize_fnnparams = (uncubed - xpmin)/(xpmax-xpmin)

  end function cubize_fnnparams


  function uncubize_fnnparams(ndim,cubed)
    implicit none
    integer(ip), intent(in) :: ndim
    real(fp), dimension(ndim) :: uncubize_fnnparams
    real(fp), dimension(ndim), intent(in) :: cubed

    if (ndim.ne.size(cubed,1)) then
       stop 'uncubize_fnnparams: sizes do not match!'
    endif

    uncubize_fnnparams = xpmin + (xpmax-xpmin)*cubed

  end function uncubize_fnnparams
 


end module fnnlike
