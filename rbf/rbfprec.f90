!   This file is part of mlearn.
!
!   Copyright (C) 2021 C. Ringeval
!   
!   mlearn is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   mlearn is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Foobar.  If not, see <https://www.gnu.org/licenses/>.

module rbfprec
  use mlprec

  implicit none

  public

  
contains


  subroutine nearest_grid_point(ndim,isize,x,fvalue,fngp,fcount)
    implicit none    
    integer(ip), intent(in) :: ndim
    integer(ip), dimension(ndim), intent(in) :: isize
    real(fp), dimension(ndim), intent(in) :: x
    real(fp), intent(in) :: fvalue
    real(fp), dimension(product(isize)), intent(inout) :: fngp
    integer(ip), dimension(product(isize)), intent(inout), optional :: fcount

    integer(ip), dimension(ndim) :: ivec
    integer(ip) :: n,q
    

    ivec(:) = 1 + int(x(:)*real(isize(:)-1,fp) + 0.5_fp)

    
    q = flatten_indices(ndim,isize,ivec)

!do a recursive mean    
    if (present(fcount)) then

       fcount(q) = fcount(q) + 1

       fngp(q) = (fvalue + fngp(q)*real(fcount(q)-1,fp))/real(fcount(q),fp)

    else
!if not, we integrate       

       fngp(q) = fvalue + fngp(q)

    end if
    
  end subroutine nearest_grid_point



  
  function flatten_indices(ndim,isize,ivec)
    implicit none
    integer(ip), intent(in) :: ndim
    integer(ip), intent(in) :: isize(ndim)
    integer(ip), dimension(ndim), intent(in) :: ivec

    integer(ip) :: flatten_indices

    integer(ip) :: q,j

    q = ivec(1)
    do j=2,ndim
       q = q + product(isize(1:j-1))*(ivec(j)-1)
    enddo

    flatten_indices = q

  end function flatten_indices


  function unflatten_indices(ndim,isize,q)
    implicit none
    integer(ip), intent(in) :: ndim
    integer(ip), dimension(ndim) :: unflatten_indices
    integer(ip), intent(in) :: isize(ndim)
    integer(ip), intent(in) :: q

    integer(ip), dimension(ndim) :: ivec
    integer(ip) :: sum
    integer(ip) :: k

    if (ndim.eq.1) then
       ivec(1) = q
       return
    endif

    ivec(ndim) = (q-1)/product(isize(1:ndim-1)) + 1
    sum = 0

    do k=ndim-1,1,-1
       sum = sum + (ivec(k+1)-1) * product(isize(1:k))
       if (k.gt.1) then
          ivec(k) = (q-1 - sum)/product(isize(1:k-1)) + 1
       else
          ivec(k) = (q-1 - sum) + 1
       endif
    enddo

    unflatten_indices = ivec

  end function unflatten_indices

end module rbfprec
