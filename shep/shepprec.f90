module shepprec
  use prec
end module shepprec


module real_precision
  use shepprec, only : fp
  integer, parameter :: r8 = fp
end module real_precision
