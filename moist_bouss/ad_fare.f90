SUBROUTINE FARE_B()
    IMPLICIT NONE
  !*** end of while **************  
  ! interface
  !     !This is qv(z) or qvs (z)
  !     pure subroutine FQv(qv,m,dz,qv0)
  !     use, intrinsic :: iso_c_binding
  !     implicit none
  !         integer(C_INT) :: ix,iy,iz
  !         real(C_DOUBLE) :: zk,k0,k1,k2,k3,k4, qvs, pz
  !         integer(C_INT), intent(in) :: m
  !         real(C_DOUBLE), intent(in) :: dz, qv0
  !         real(C_DOUBLE), dimension(m),intent(out) :: qv
  !     end subroutine FQv
  !     !This is dqv/dz for the initial data
  !     pure subroutine Fdqvdz(dqvdz,m,dz,qv0)
  !     use, intrinsic :: iso_c_binding
  !     implicit none
  !         integer(C_INT) :: iz
  !         real(C_DOUBLE) :: k0,z_k,k1,k2,k3,k4, qvs, pz
  !         integer(C_INT), intent(in) :: m
  !         real(C_DOUBLE), intent(in) :: dz, qv0
  !         real(C_DOUBLE), dimension(m),intent(out) :: dqvdz
  !     end subroutine Fdqvdz
  ! end interface
  !Integer
    INTEGER :: it, ix, iy, iz
    INTEGER, SAVE :: nx=128
    INTEGER, SAVE :: ny=128
    INTEGER, SAVE :: m=101
    INTEGER, SAVE :: nxh=65
    INTEGER, SAVE :: nxp=86
    INTEGER, SAVE :: nxph=43
    INTEGER, SAVE :: nyp=86
    INTEGER, SAVE :: nyph=43
    REAL, SAVE :: pi=3.14
    REAL, SAVE :: pib=0.0
    REAL :: us, ts, ths
    REAL :: lx, ly, lz, dx, dy, dz
    REAL :: l, vt, tau, qv0
    REAL :: a_squall, f
    INTEGER :: ls, qs, nn
    REAL :: tfinal, im, eps, epsbar, nz, f_star, g_star, b_star
    REAL :: drx, dry, drz, zk
  !Real
    REAL :: max_theta, max_pert
    REAL :: xi, yj, da, x_c, y_c, z_c, r_c, ampl_bubble
    REAL, DIMENSION(128, 128, 101) :: u, v, w, theta, thetar, qt, qr, qv
    REAL, DIMENSION(101) :: qvini, rc, mr, ubg, vbg, dqvdz, qvs, theta_bar&
  & , u_bar, v_bar, wrk1d, qvs0
    INTEGER :: ti, mindt
    REAL :: dt0, dt, cfl, n_theta
    REAL, DIMENSION(2) :: nu
    COMPLEX, DIMENSION(128, 128) :: in, wrk
    COMPLEX, DIMENSION(65, 128) :: out
  ! complex, dimension(43,86) :: wrk;
    COMPLEX, DIMENSION(101) :: tau_z
  ! complex, dimension(43,86,101) :: uhat, vhat, what, ThetaRhat, qthat
  ! complex, dimension(43,86,101) :: u1hat, v1hat, w1hat, ThetaR1hat, qt1hat
  ! complex, dimension(43,86,101) :: fu1hat, fv1hat, fw1hat, fThetaR1hat, fqt1hat
    COMPLEX, DIMENSION(128, 128, 101) :: uhat, vhat, what, thetarhat, &
  & qthat
    COMPLEX, DIMENSION(128, 128, 101) :: u1hat, v1hat, w1hat, thetar1hat, &
  & qt1hat
    COMPLEX, DIMENSION(128, 128, 101) :: fuhat, fvhat, fwhat, fthetarhat, &
  & fqthat
    COMPLEX, DIMENSION(128, 128, 101) :: fu1hat, fv1hat, fw1hat, &
  & fthetar1hat, fqt1hat
    COMPLEX, DIMENSION(43, 86, 101) :: kx, ky, kk
    COMPLEX, DIMENSION(43, 86, 101) :: e_nu_1, e_nu_2, e_nu_3
    COMPLEX, DIMENSION(43, 86, 101) :: e_nu_1uv, e_nu_2uv, e_nu_3uv
    COMPLEX, DIMENSION(43, 86, 101) :: e_nu_1w, e_nu_2w, e_nu_3w
    COMPLEX, DIMENSION(43, 86, 100) :: a, b, c, rhat, phat
    INTRINSIC SIN
    INTRINSIC EXP
    INTRINSIC SQRT
    INTRINSIC COS
    INTRINSIC SUM
    INTRINSIC ABS
    INTRINSIC MAXVAL
    INTRINSIC MAX
    INTRINSIC MIN
    REAL :: y1
    REAL :: x1
    REAL :: y2
    REAL, DIMENSION(4:m-1) :: abs0
    REAL, DIMENSION(128, 128, 101) :: abs1
    REAL :: max1
    REAL, DIMENSION(128, 128, 101) :: abs2
    REAL, DIMENSION(128, 128, 101) :: abs3
    REAL, DIMENSION(128, 128, 101) :: abs4
    LOGICAL, DIMENSION(128, 128) :: mask
    LOGICAL, DIMENSION(128, 128) :: mask0
    LOGICAL, DIMENSION(m-4) :: mask1
    REAL :: result1
    REAL, DIMENSION(128, 128, 101) :: arg1
    LOGICAL, DIMENSION(128, 128, 101) :: mask2
    LOGICAL, DIMENSION(128, 128, 101) :: mask3
    LOGICAL, DIMENSION(128, 128, 101) :: mask4
    LOGICAL, DIMENSION(128, 128, 101) :: mask5
    LOGICAL, DIMENSION(128, 128) :: mask6
    LOGICAL, DIMENSION(128, 128) :: mask7
  !----------------PARAMETERS----------------------------------------------
  !Grid size:
  ! nx  = 128; ny  = 128; m  = 100+1;
  ! nxh = nx/2+1; nxp = nx/3; nxp = 2*nxp+1; nxph = (nxp-1)/2+1;
  ! nyp = ny/3; nyp = 2*nyp+1; nyph = (nyp-1)/2+1;
  !Quantity scales:
  !In m/s
  !In hours
  !In Kelvin
  !In kms
  !In g/kg
  !*******DOMAIN SIZE*******
  !Non-dimentionalization
  !Rain fall velocity
  !Non-dimentionalization
  !Latent heat
  !In 10^{6}*J*kg^{-1}
  !Non-dimentionalization prefector epsilon^{-1} L = epsilon^2 * L^d * \theta_0 ([\Theta]*c_p*T_0)^{-1}
  !Pendiente: Cambiar a LcpTheta0
  !Relaxation time
  ! In hours
  !Non-dimensionalization
  !Turned off
  !Qvs at surface
  !In g/kg
  !Non-dimensionalization 
  !In g/kg
  !Max perturbation in theta
  !In Kelvin
  !Non-dimensionalization
  !Parameter for the zonal vel background
  !In m/s
  !Turned off
  !Coriolis parameter
  !Final time
  !In hours
  !Non-dimensionalization
  !***************
  !***************
  !***************
  !7.9461;
  !Wave numbers are multiple of those
  !!!!Data for Bubble:
  !rx_c = 10; ry_c = 10; rz_c = 1;
  !rx_c = rx_c/Ls; ry_c = ry_c/Ls; rz_c = rz_c/Ls; 
   100 CONTINUE
   110 CONTINUE
   120 CONTINUE
  END SUBROUTINE FARE_B