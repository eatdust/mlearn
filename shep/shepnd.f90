module shepnd
  use shepprec
  implicit none

  logical, parameter :: display = .true.

  private

  public shepard_maxradius, shepard_eval

contains

  function shepard_num_cells(ndim,ndata)
    implicit none
    integer(ip) :: shepard_num_cells
    integer(ip), intent(in) :: ndim, ndata
    
    shepard_num_cells = nint((real(ndata,fp)/3._fp)**(1._fp/real(ndim,fp)))

!    shepard_num_cells = shepard_num_cells - 5

  end function shepard_num_cells



  function shepard_maxradius(ndim,ndata,xdata,fdata,nctrs,nq)
    use qshepmd_mod, only : qshepm
    implicit none
    real(fp) :: shepard_maxradius
    integer(ip), intent(in) :: ndim, ndata, nctrs,nq
   
    real(fp), dimension(ndim,ndata), intent(in) :: xdata
    real(fp), dimension(ndata), intent(in) :: fdata
    
    real(fp) :: rmax

    integer :: ncells,info


    ncells = shepard_num_cells(ndim,ndata)
    
    if (display) then
       write(*,*)'shepard_maxradius: ncells= ',ncells
    end if

    call qshepm(ndim,ndata,xdata,fdata,nq,nctrs,ncells,rmax,info)

    select case (info)

    case (0)
       if (display) then
          write(*,*)'shepmd_maxradius:'
          write(*,*)'rmax= ',rmax
       endif
    case (1)
       write(*,*)'shepmd_maxradius: input out of range'
       stop
    case (2)
       write(*,*)'shepmd_maxradius: duplicate nodes found'
       stop
    case (3)
       write(*,*)'shepmd_maxradius: nodes lies in an affine subspace!'
       stop
    case (4)
       write(*,*)'shepmd_maxradius: memory allocation error!'
       stop
    case default
       stop 'internal error, RTFM!'
    end select

    shepard_maxradius = rmax

  end function shepard_maxradius



  function shepard_eval(ndim,ndata,xdata,fdata,rmax,x)
    use qshepmd_mod, only : qsmval
    implicit none
    real(fp) :: shepard_eval
    integer(ip), intent(in) :: ndim, ndata
    real(fp), dimension(ndim,ndata), intent(in) :: xdata
    real(fp), dimension(ndata), intent(in) :: fdata
    real(fp), intent(in) :: rmax
    real(fp), dimension(ndim), intent(in) :: x

    integer(ip) :: ncells

    ncells = shepard_num_cells(ndim,ndata)

    shepard_eval = qsmval(ndim,ndata,x,xdata,fdata,ncells,rmax)

  end function shepard_eval
  
 

end module shepnd
