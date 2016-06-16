  interface
     subroutine phi( N, R, scale, V )
       use rbfprec
       implicit none
       integer (ip) :: n
       real (fp) :: r(n)
       real (fp) :: scale
       real (fp) :: v(n)
     end subroutine phi
  end interface
