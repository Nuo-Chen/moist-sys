

pure subroutine FQv(qv,m,dz,qv0)
use, intrinsic :: iso_c_binding
implicit none
	!Parameters
	!Integers
	integer(C_INT) :: ix,iy,iz
	!Real
	real(C_DOUBLE) :: zk,k0,k1,k2,k3,k4, pz
	!In:
	!Integers
	integer(C_INT), intent(in) :: m
	!Real
	real(C_DOUBLE), intent(in) :: dz, qv0
	!Out
	!Real
	real(C_DOUBLE), dimension(m),intent(out) :: qv
!-------------------------------------------------------------------------------------	
	k0 = 18.04d0; k1 = 3.27d0; k2 = 0.1d0; k3 = 0.1804d0; k4 = 3.48d0; !Pendiente: definir cada k
	do iz=1,m
		zk = (iz-0.5d0)*dz;
		pz = (1d0-k1*log(1d0+k2*zk))**k4;

		qv(iz) = (qv0/pz)*exp( -k0*( 1d0/((1d0-k1*log(1d0+k2*zk))*(1d0+k2*zk))-1d0 ) );
	end do
end subroutine FQv

pure subroutine Fdqvdz(dqvdz,m,dz,qv0)
use, intrinsic :: iso_c_binding
implicit none
	!Parameters
	!Integers
	integer(C_INT) :: iz
	!Real
	real(C_DOUBLE) :: zk, k0, k1,k2,k3,k4, qv, pz
	!In:
	!Integers
	integer(C_INT), intent(in) :: m
	!Real
	real(C_DOUBLE), intent(in) :: dz, qv0
	!Out:
	!Real
	real(C_DOUBLE), dimension(m),intent(out) :: dqvdz
!-------------------------------------------------------------------------------------	
	k0 = 18.04d0; k1 = 3.27d0; k2 = 0.1d0; k3 = 0.1804d0; k4 = 3.48d0;
	do iz=1,m
		zk = (iz-0.5d0)*dz;
		pz = (1d0-k1*log(1d0+k2*zk))**k4;

		qv = (qv0/pz)*exp( -k0*( 1d0/((1d0-k1*log(1d0+k2*zk))*(1d0+k2*zk))-1d0 ) );
		dqvdz(iz) = qv*( k1*k2*k4/((1d0-k1*log(1d0+k2*zk))*(1d0+k2*zk)) &
		+ k0*k2*(1d0-k1-k1*log(1d0+k2*zk))/( ((1d0-k1*log(1d0+k2*zk))**2d0)*((1d0+k2*zk)**2d0) ) );
	end do
end subroutine Fdqvdz
