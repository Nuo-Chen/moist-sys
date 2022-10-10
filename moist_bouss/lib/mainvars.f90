!-----------------------------MODULE--------------------------

module mainvars
	use, intrinsic :: iso_c_binding
	!Parameters:
	!Character
	character(len=30) :: filename;
	!Integer
    integer, parameter :: long=selected_real_kind(16,99)
	integer(C_INT) :: ix,iy,iz,nx,nxh,nxp,nxph,ny,nyp,nyph,m,nn
	integer(C_INT) :: nthreads
	type(C_PTR) :: planf,planr
	integer(C_INT), dimension(:), allocatable :: seed
	integer(C_INT) :: fileindex, iret, n
	!Parameters:
    real(long), parameter ::pi=3.14159265358979323846_long
	!Real:
	real(C_DOUBLE) :: dt,dx,dy,dz,zk,Ti,tau,a_squall,dt0,N_theta,CFL,mindt
	real(C_DOUBLE) :: eps,f,epsbar,nu(2),Nz,L,vt
	real(C_DOUBLE) :: savedata,savedata3D, movies2D,movies3D
	real(C_DOUBLE) :: Lx,Ly,Lz,drx,dry,drz,qv0,qvs0, Tfinal, pertrand(100000000)
	real(C_DOUBLE) :: Us, Ts, Ths, Ls, qs, f_star, g_star, B_star
	!Complex
	complex(C_DOUBLE_COMPLEX) :: IM
	!Allocatable
	!Real:
	real(C_DOUBLE) ,dimension(:,:,:), allocatable ::  u,v,w,Theta,ThetaR,qv,qt,qr,kx,ky,kk,a,b,c
	real(C_DOUBLE) ,dimension(:,:,:), allocatable :: e_nu_1,e_nu_2,e_nu_3, e_nu_1uv,e_nu_2uv,e_nu_3uv,  e_nu_1w,e_nu_2w,e_nu_3w
	real(C_DOUBLE) ,dimension(:,:), allocatable :: in
	real(C_DOUBLE) ,dimension(:), allocatable :: qvini
	real(C_DOUBLE) ,dimension(:), allocatable :: dqvdz, RC, Mr, ubg, vbg, qvs, tau_z
	!Complex:
	complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: phat,rhat,uhat,vhat,what,ThetaRhat,qthat
	complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: fuhat,fvhat,fwhat,fThetaRhat,fqthat
	complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: u1hat,v1hat,w1hat,ThetaR1hat,qt1hat
	complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: fu1hat,fv1hat,fw1hat,fThetaR1hat,fqt1hat
	complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: out, wrk
end module mainvars

!COMENTED CODE:

!----------------------------

	!!Cold pool:
	!br = 1d0;
	!bcx = Lx/2d0;
	
	!do ix=1,nx
	!	xi = (ix-1)*dx;
	!	da = sqrt((xi-bcx)**2d0);
	!	do iz=1,m-1
	!		zk = (iz-0.5)*dz;
	!		if (da .le. br .and. zk .ge. 0d0 .and. zk .le. 0.2d0) then
	!			do iy=1,ny
	!				!Cold pool, uniform in y:	
	!				Theta(ix,iy,iz) = -(1d0/3d0)*cos(pi*da/br*0.5d0)*(10d0*(zk-0d0))**2d0*(10d0*(zk-0.2d0))**2d0;
	!				!Now add the perturbation (up to +- 10% of its value from above)
	!				Theta(ix,iy,iz) = Theta(ix,iy,iz)*(1d0+2d0*(pertrand( (ix-1)*ny*(m-1)+(iy-1)*(m-1)+iz )-0.5d0)*0.1d0) !*max_pert;
	!			end do
	!		end if
	!	end do
	!end do

!----------------------------

	!open(unit=1,file='U3DT12.dat',form="formatted",status="old",action="read")
	!do iz=1,m
	!	do iy=1,ny
	!		do ix=1,nx
	!			read(1,*) u(ix,iy,iz);
	!		end do
	!	end do
	!end do
	!close(1)
	!u = u/Us;

	!open(unit=1,file='V3DT12.dat',form="formatted",status="old",action="read")
	!do iz=1,m
	!	do iy=1,ny
	!		do ix=1,nx
	!			read(1,*) v(ix,iy,iz);
	!		end do
	!	end do
	!end do
	!close(1)
	!v = v/Us;

	!open(unit=1,file='W3DT12.dat',form="formatted",status="old",action="read")
	!do iz=1,m
	!	do iy=1,ny
	!		do ix=1,nx
	!			read(1,*) w(ix,iy,iz);
	!		end do
	!	end do
	!end do
	!close(1)
	!w = w/Us;

	!open(unit=1,file='Theta3DT12.dat',form="formatted",status="old",action="read")
	!do iz=1,m
	!	do iy=1,ny
	!		do ix=1,nx
	!			read(1,*) Theta(ix,iy,iz);
	!		end do
	!	end do
	!end do
	!close(1)
	!Theta = Theta/Ths;

	!open(unit=1,file='Qv3DT12.dat',form="formatted",status="old",action="read")
	!do iz=1,m
	!	do iy=1,ny
	!		do ix=1,nx
	!			read(1,*) qv(ix,iy,iz);
	!		end do
	!	end do
	!end do
	!close(1)
	!qv = qv/qs;

	!open(unit=1,file='Qr3DT12.dat',form="formatted",status="old",action="read")
	!do iz=1,m
	!	do iy=1,ny
	!		do ix=1,nx
	!			read(1,*) qr(ix,iy,iz);
	!		end do
	!	end do
	!end do
	!close(1)
	!qr = qr/Us;

