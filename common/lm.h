INTERFACE
  SUBROUTINE fcn1(m, n, x, fvec, iflag)
    use mlprec, only : fp
    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: m, n
    REAL (fp), INTENT(IN)      :: x(:)
    REAL (fp), INTENT(IN OUT)  :: fvec(:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn1



  SUBROUTINE fcn2(m, n, x, fvec, fjac, iflag)
    use mlprec, only : fp
    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: m, n
    REAL (fp), INTENT(IN)      :: x(:)
    REAL (fp), INTENT(IN OUT)  :: fvec(:)
    REAL (fp), INTENT(OUT)     :: fjac(:,:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn2

END INTERFACE
