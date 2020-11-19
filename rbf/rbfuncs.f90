module rbfuncs
  use rbfprec
  implicit none

  public

contains


  subroutine rbf_inverse_monomial(n, r, r0, v, k)
    implicit none
    integer(ip) :: n
    real(fp), dimension(n) :: r
    real(fp) :: r0
    real(fp), dimension(n) :: v
    real(fp) :: k
    
    v = 1._fp/(1._fp + (r/r0)**k)
    
  end subroutine rbf_inverse_monomial


  subroutine rbf_inverse_monomial_two(n, r, r0, v)
    implicit none
    integer(ip) :: n
    real(fp), dimension(n) :: r
    real(fp) :: r0
    real(fp), dimension(n) :: v
    real(fp) :: k

    call rbf_inverse_monomial(n,r,r0,v,2._fdp)

  end subroutine rbf_inverse_monomial_two
  
  
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


  
  subroutine rbf_polyharmonic_log_two(n, r, r0, v )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)

    call rbf_polyharmonic_log(n, r, r0, v,2_ip)

  end subroutine rbf_polyharmonic_log_two


  
  subroutine rbf_polyharmonic_log_four(n, r, r0, v )
    implicit none
    integer(ip) :: n
    real(fp) :: r(n)
    real(fp) :: r0
    real(fp) :: v(n)

    call rbf_polyharmonic_log(n, r, r0, v,4_ip)

  end subroutine rbf_polyharmonic_log_four



  
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

    call rbf_polyharmonic(n, r, r0, v,2_ip)

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

    call rbf_polyharmonic(n, r, r0, v,4_ip)

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

    call rbf_polyharmonic(n, r, r0, v,6_ip)

  end subroutine rbf_polyharmonic_six

  
 
end module rbfuncs