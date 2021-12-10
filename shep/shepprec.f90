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

module shepprec
  use mlprec
end module shepprec


module real_precision
  use shepprec, only : fp
  integer, parameter :: r8 = fp
end module real_precision
