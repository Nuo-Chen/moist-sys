!--- Subroutine Thomas for Pressure Poisson Equation ----
subroutine Thomas(x,a,b,c,r,nxph,nyp,m);
    use, intrinsic :: iso_c_binding
    implicit none
        !Allocatable
        complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: fac
        complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: cwrk
        !Parameters:
        !Integers
        integer(C_INT) :: i
        !In:
        !Integers
        integer(C_INT), intent(in) :: nxph,nyp,m
        !Real
        real(C_DOUBLE), dimension(nxph,nyp,m), intent(in) :: a,b,c
        !Out:
        !Complex
        complex(C_DOUBLE_COMPLEX), dimension(nxph,nyp,m), intent(inout) :: r
        complex(C_DOUBLE_COMPLEX), dimension(nxph,nyp,m), intent(out) :: x
        !---------------------------------------------------
    
        allocate(fac(nxph,nyp)); allocate(cwrk(nxph,nyp,m));
    
        fac = 1d0/b(:,:,1);
        cwrk(:,:,1) = c(:,:,1) * fac;
        r(:,:,1)    = r(:,:,1) * fac;
        do i = 2,m
            fac = 1d0/( b(:,:,i) - a(:,:,i)*cwrk(:,:,i-1) );
            cwrk(:,:,i) = c(:,:,i)*fac;
            r(:,:,i) = ( r(:,:,i) - r(:,:,i-1)*a(:,:,i) )*fac;
        end do
        x(:,:,m) = r(:,:,m);
        do i = m-1,1,-1
            x(:,:,i) = r(:,:,i) - cwrk(:,:,i)*x(:,:,i+1);
        end do
        deallocate(fac,cwrk);
        return;
    end subroutine Thomas

!--- Subroutine Thomas for Pressure Poisson Equation ----
subroutine Thomas(x,a,b,c,r,nxph,nyp,m)
    use, intrinsic :: iso_c_binding
    implicit none

    !Parameters:
    !Integers
    integer :: i
            !Allocatable
    complex :: fac(nxph,nyp)
    complex:: cwrk(nxph,nyp,m)
    !In:
    !Integers
    integer, intent(in) :: nxph,nyp,m
    !Real
    real, dimension(nxph,nyp,m), intent(in) :: a,b,c
    !Out:
    !Complex
    complex, dimension(nxph,nyp,m), intent(inout) :: r
    complex, dimension(nxph,nyp,m), intent(out) :: x
    
    !---------------------------------------------------

    !allocate(fac(nxph,nyp)); allocate(cwrk(nxph,nyp,m));

    fac = 1d0/b(:,:,1);
    cwrk(:,:,1) = c(:,:,1) * fac;
    r(:,:,1)    = r(:,:,1) * fac;
    do i = 2,m
        fac = 1d0/( b(:,:,i) - a(:,:,i)*cwrk(:,:,i-1) );
        cwrk(:,:,i) = c(:,:,i)*fac;
        r(:,:,i) = ( r(:,:,i) - r(:,:,i-1)*a(:,:,i) )*fac;
    end do
    x(:,:,m) = r(:,:,m);
    do i = m-1,1,-1
        x(:,:,i) = r(:,:,i) - cwrk(:,:,i)*x(:,:,i+1);
    end do
    !deallocate(fac,cwrk);
    return
end subroutine Thomas

SUBROUTINE THOMAS_D(x, xd, a, b, c, r, rd, nxph, nyp, m)
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
  !Parameters:
  !Integers
    INTEGER :: i
    INTEGER, INTENT(IN) :: nxph, nyp, m
  !In:
  !Integers
  !Allocatable
    COMPLEX :: fac(nxph, nyp)
    COMPLEX :: cwrk(nxph, nyp, m)
  !Real
    REAL, DIMENSION(nxph, nyp, m), INTENT(IN) :: a, b, c
  !Out:
  !Complex
    COMPLEX, DIMENSION(nxph, nyp, m), INTENT(INOUT) :: r
    COMPLEX, DIMENSION(nxph, nyp, m), INTENT(INOUT) :: rd
    COMPLEX, DIMENSION(nxph, nyp, m), INTENT(OUT) :: x
    COMPLEX, DIMENSION(nxph, nyp, m), INTENT(OUT) :: xd
    fac = 1d0/b(:, :, 1)
  !---------------------------------------------------
  !allocate(fac(nxph,nyp)); allocate(cwrk(nxph,nyp,m));
    cwrk(:, :, 1) = c(:, :, 1)*fac
    rd(:, :, 1) = fac*rd(:, :, 1)
    r(:, :, 1) = r(:, :, 1)*fac
    DO i=2,m
      fac = 1d0/(b(:, :, i)-a(:, :, i)*cwrk(:, :, i-1))
      cwrk(:, :, i) = c(:, :, i)*fac
      rd(:, :, i) = fac*(rd(:, :, i)-a(:, :, i)*rd(:, :, i-1))
      r(:, :, i) = (r(:, :, i)-r(:, :, i-1)*a(:, :, i))*fac
    END DO
    xd = (0.0,0.0)
    xd(:, :, m) = rd(:, :, m)
    x(:, :, m) = r(:, :, m)
    DO i=m-1,1,-1
      xd(:, :, i) = rd(:, :, i) - cwrk(:, :, i)*xd(:, :, i+1)
      x(:, :, i) = r(:, :, i) - cwrk(:, :, i)*x(:, :, i+1)
    END DO
  !deallocate(fac,cwrk);
    RETURN
  END SUBROUTINE THOMAS_D

  SUBROUTINE THOMAS_B(x, xb, a, b, c, r, rb, nxph, nyp, m)
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
  !Parameters:
  !Integers
    INTEGER :: i
    INTEGER, INTENT(IN) :: nxph, nyp, m
  !In:
  !Integers
  !Allocatable
    COMPLEX :: fac(nxph, nyp)
    COMPLEX :: cwrk(nxph, nyp, m)
  !Real
    REAL, DIMENSION(nxph, nyp, m), INTENT(IN) :: a, b, c
  !Out:
  !Complex
    COMPLEX, DIMENSION(nxph, nyp, m), INTENT(INOUT) :: r
    COMPLEX, DIMENSION(nxph, nyp, m), INTENT(INOUT) :: rb
    COMPLEX, DIMENSION(nxph, nyp, m) :: x
    COMPLEX, DIMENSION(nxph, nyp, m) :: xb
  !---------------------------------------------------
  !allocate(fac(nxph,nyp)); allocate(cwrk(nxph,nyp,m));
    fac = 1d0/b(:, :, 1)
    cwrk(:, :, 1) = c(:, :, 1)*fac
    DO i=2,m
      CALL PUSHCOMPLEX8ARRAY(fac, nxph*nyp)
      fac = 1d0/(b(:, :, i)-a(:, :, i)*cwrk(:, :, i-1))
      cwrk(:, :, i) = c(:, :, i)*fac
    END DO
    DO i=1,m-1,1
      rb(:, :, i) = rb(:, :, i) + xb(:, :, i)
      xb(:, :, i+1) = xb(:, :, i+1) + CONJG(-cwrk(:, :, i))*xb(:, :, i)
      xb(:, :, i) = (0.0,0.0)
    END DO
    rb(:, :, m) = rb(:, :, m) + xb(:, :, m)
    xb(:, :, m) = (0.0,0.0)
    DO i=m,2,-1
      rb(:, :, i-1) = rb(:, :, i-1) + CONJG(-(a(:, :, i)*fac))*rb(:, :, i)
      rb(:, :, i) = CONJG(fac)*rb(:, :, i)
      CALL POPCOMPLEX8ARRAY(fac, nxph*nyp)
    END DO
    fac = 1d0/b(:, :, 1)
    rb(:, :, 1) = CONJG(fac)*rb(:, :, 1)
  END SUBROUTINE THOMAS_B