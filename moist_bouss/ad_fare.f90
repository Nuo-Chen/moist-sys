! SUBROUTINE FARE_B()
! 	implicit none
! 	integer :: it, ix, iy, iz
! 	integer :: nx = 128
! 	integer :: ny = 128
! 	integer :: m = 101;
! 	integer :: nxh = 65
! 	integer :: nxp = 86
! 	integer ::  nxph =43
! 	integer :: nyp = 86
! 	integer :: nyph = 43
! 	real :: pi = 3.14

! 	real :: Us, Ts, Ths
! 	real :: Lx, Ly, Lz, dx, dy, dz
! 	real :: L, vt, tau, qv0
! 	real :: a_squall, f
! 	integer :: Ls, qs, nn
! 	real :: Tfinal, IM, eps, epsbar, Nz, f_star, g_star, B_star
! 	real :: drx, dry, drz, zk

! 	!Real
! 	real :: max_theta, max_pert
! 	real :: xi,yj,da,x_c,y_c,z_c,r_c, ampl_bubble

! 	real, dimension(128,128,101) :: u,v,w, Theta, ThetaR, qt, qr, qv;
! 	real, dimension(101) :: qvini, RC, Mr, ubg, vbg, dqvdz, qvs, theta_bar, u_bar, v_bar, wrk1D, qvs0

! 	integer :: Ti, mindt
! 	real :: dt0, dt, CFL, N_theta
! 	real, dimension(2) :: nu

! 	complex, dimension(128,128) :: in, wrk
! 	complex, dimension(65,128) :: out;
! 	! complex, dimension(43,86) :: wrk;
! 	complex, dimension(101) :: tau_z;
! 	! complex, dimension(43,86,101) :: uhat, vhat, what, ThetaRhat, qthat
! 	! complex, dimension(43,86,101) :: u1hat, v1hat, w1hat, ThetaR1hat, qt1hat
! 	! complex, dimension(43,86,101) :: fu1hat, fv1hat, fw1hat, fThetaR1hat, fqt1hat
! 	complex, dimension(128,128,101) :: uhat, vhat, what, ThetaRhat, qthat
! 	complex, dimension(128,128,101) :: u1hat, v1hat, w1hat, ThetaR1hat, qt1hat
! 	complex, dimension(128,128,101) :: fuhat, fvhat, fwhat, fThetaRhat, fqthat
! 	complex, dimension(128,128,101) :: fu1hat, fv1hat, fw1hat, fThetaR1hat, fqt1hat

! 	complex, dimension(43,86,101) :: kx, ky, kk
! 	complex, dimension(43,86,101) :: e_nu_1, e_nu_2, e_nu_3
! 	complex, dimension(43,86,101) :: e_nu_1uv, e_nu_2uv, e_nu_3uv
! 	complex, dimension(43,86,101) :: e_nu_1w, e_nu_2w, e_nu_3w
! 	complex, dimension(43,86,100) :: a, b, c, rhat, phat


! 	Tfinal = 5; !In hours
! 	Tfinal = Tfinal/Ts; !Non-dimensionalization

! 	Ti = Tfinal + 5*dt 
! 	do while (Ti>0)
! 		adj_Runge_Kutta_loop:  DO rk_step = rk_order, 1, -1
! 		IF ( rk_order == 3 ) THEN ! third order Runge-Kutta

! 			IF ( rk_step == 1) THEN
! 			dt_rk = grid%dt/3.
! 			dts_rk = dt_rk
! 			number_of_small_timesteps = 1
! 			ELSE IF (rk_step == 2) THEN
! 			dt_rk  = 0.5*grid%dt
! 			dts_rk = dts
! 			number_of_small_timesteps = num_sound_steps/2
! 			ELSE
! 			dt_rk = grid%dt
! 			dts_rk = dts
! 			number_of_small_timesteps = num_sound_steps
! 			ENDIF
! 		END IF adj_Runge_Kutta_loop
		
! 		Ti = Ti - dt

! 	end do adj_Runge_Kutta_loop


! END SUBROUTINE FARE_B


subroutine simplefare(u, uhat, u1hat, fuhat, fu1hat, rhat, phat, wrk, in, out)
	implicit none 

	integer :: it, ix, iy, iz

	integer :: nx = 128
	integer :: ny = 128
	integer :: m = 101;
	! integer :: nxh = 65
	! integer :: nxp = 86
	integer ::  nxph = 128!43
	integer :: nyp = 128 !86
	integer :: nyph = 64 !43
	! real :: pi = 3.14

	real :: Us !, Ts, Ths
	real :: Lx, Ly, Lz, dx, dy, dz
	real :: f, Ls, Ts
	! real :: L, vt, tau, qv0
	! real :: a_squall, f
	! integer :: Ls, qs, 
	integer :: nn = 128*128
	real :: Tfinal, IM!, eps, epsbar, Nz, f_star, g_star, B_star
	real :: drx, dry, drz, zk

	!Real
	! real :: max_theta, max_pert
	! real :: xi,yj,da,x_c,y_c,z_c,r_c, ampl_bubble

	real, dimension(128,128,101) :: u ! ,v,w, Theta, ThetaR, qt, qr, qv;
	! real, dimension(101) :: qvini, RC, Mr, ubg, vbg, dqvdz, qvs, theta_bar, u_bar, v_bar, wrk1D, qvs0
	real, dimension(101) :: u_bar, ubg

	real :: Ti, mindt
	real :: dt0, dt, CFL!, N_theta
	real, dimension(2) :: nu

	real,  dimension(128,128) :: in, wrk, out 
	real,  dimension(101) :: tau_z;
	real,  dimension(128,128,101) :: uhat!, vhat, what, ThetaRhat, qthat
	real,  dimension(128,128,101) :: u1hat!, v1hat, w1hat, ThetaR1hat, qt1hat
	real,  dimension(128,128,101) :: fuhat!, fvhat, fwhat, fThetaRhat, fqthat
	real,  dimension(128,128,101) :: fu1hat!, fv1hat, fw1hat, fThetaR1hat, fqt1hat

	real,  dimension(128,128,101) :: kx, ky, kk
	!real,  dimension(43,86,101) :: e_nu_1, e_nu_2, e_nu_3
	real,  dimension(128,128,101) :: e_nu_1uv, e_nu_2uv, e_nu_3uv
	!real,  dimension(43,86,101) :: e_nu_1w, e_nu_2w, e_nu_3w
	real,  dimension(128,128,100) :: a, b, c, rhat, phat

	!Quantity scales:
	Us = 100/9; !In m/s
	Ls = 10; !In kms
	Ts = 0.25
	Lx = 128; Ly = 128; Lz = 15;
	Lx = Lx/Ls; Ly = Ly/Ls; Lz = Lz/Ls; !Non-dimentionalization
	dx = Lx/nx; dy = Ly/ny; dz = Lz/(m-1);
	f = sin(25* 3.14 /180);
	!Final time
	Tfinal = 5; !In hours
	Tfinal = Tfinal/Ts; !Non-dimensionalization
	!***************
	IM = 1!(0,1)
	tau_z = 0

	!Wave numbers are multiple of those
	drx = 2*3.14/Lx; dry = 2*3.14/Ly; drz = 2*3.14/Lz;

	do ix = 1, nxph
		kx(ix,:,:) = ix-1;
		do iy = 1, nyph			
			ky(ix,iy,:) = iy-1;
		end do
		do iy = nyph+1, nyp
			ky(ix,iy,:) = iy-nyp-1;
		end do
	end do
	kx = drx*kx; ky = dry*ky; kk = kx*kx+ky*ky;

	!!!!!!!!!!!!!!!Pressure solver !!!!
    !Three diagonals for Thomas' method
    !Main diagonal (entries from 1 to m-1), to compute variables in the staggered grid, level at p
    b = - ( kk(:,:,:m-1) + 1/dz**2);
    b(:,:,2:m-2) = b(:,:,2:m-2) - 1/dz**2;
    b(1,1,1) = 1/dz**2; ! 1; !Pendiente
    b(1,1,m-1) = -1/dz**2;
    !Lower diagonal, a_2, a_3, ... a_{m-1}
    a = 1/dz**2;
    a(:,:,1) = 0;
    !Upper diagonal: c_1, c_2 ... c_{m-2}
    c = 1/dz**2;
    c(:,:,m-1) = 0;
    c(1,1,1) = 0;

	u = 0
	ubg = 0

	do iz=1,m 
		in = u(:,:,iz)
		wrk = in 
		uhat(:,:,iz) = wrk/nn 
	end do 

	uhat(1,1,:) = ubg

	do iz=1,m 
		wrk = uhat(:,:,iz)
		out = wrk
		in = out 
		u(:,:,iz) = in
	end do 

!----- Time iteration start----------------------

	Ti  = 0; dt0 = 2.5*10**(-3);

	nu = 0;
	nu(1) = 4*10**(-5)/1; !0.1/(((pi/dx)**2)**2 * dt0); !
	nu(2) = 4*10**(-4)/1; !0.005*(dz**2)/dt0; !

	CFL = 0.9; mindt = 1; dt=dt0;
	do while ( Ti <= Tfinal+5*dt )

		dt = dt0; 
		do iz=1,m
			e_nu_1uv(:,:,iz) = exp(dt/3*(-nu(1)*kk(:,:,iz)**2 - nu(2)*2/dz**2 - tau_z(iz) )); 
			e_nu_1uv(1,1,iz) = exp(dt/3*(- tau_z(iz) )); 
			e_nu_2uv(:,:,iz) = e_nu_1uv(:,:,iz)**2;
			e_nu_3uv(:,:,iz) = e_nu_1uv(:,:,iz)**3;
		end do
 
            !---------------------------------------------------
            !---- RK3 for momentum equations start ---------
            !****- RK1  start ********************************
		! call RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vt,
		! uhat,vhat,what,ThetaRhat,qthat,u,v,w,ThetaR,qt,&
		! fuhat,fvhat,fwhat,fThetaRhat,fqthat,RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)
		fuhat = uhat + u  ! RKs1
			
		u1hat = (uhat + dt/3*fuhat)*e_nu_1uv;

		!******- Possion equation for pressure p -*******
		rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*u1hat(:,:,1:m-1)
		phat = 0;
		phat = a+b+c+rhat ! Thomas

		!***- Update veloctiy field ----------------------
		u1hat(:,:,1:m-1) = u1hat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat;

		!***Neumann boundary condition, w1hat(1)=w1hat(1),w1hat(m)=w1hat(m)***
		!--- Get u1,v1,w1 & ThetaR1 ---
		do iz=1,m
			wrk = u1hat(:,:,iz); 
			out = wrk
			in = out
			u(:,:,iz) = in;
		end do			

		!************  RK1 end ***********************************
		!****- RK2  start ********************************
		! call RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vt,
		! u1hat,v1hat,w1hat,ThetaR1hat,qt1hat,u,v,w,ThetaR,qt,&
		! fu1hat,fv1hat,fw1hat,fThetaR1hat,fqt1hat,RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)

		fu1hat = u1hat + u ! RKs2

		u1hat = uhat*e_nu_2uv + 2/3*dt*fu1hat*e_nu_1uv;

		!******- Possion equation for pressure p -*******
		rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*u1hat(:,:,1:m-1)
		phat = 0;
		phat = a+b+c+rhat ! Thomas
		!***- Update veloctiy field ----------------------
		u1hat(:,:,1:m-1) = u1hat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat;

		!***Neumann boundary condition, w1hat(1)=w1hat(1),w1hat(m)=w1hat(m)***
		!--- Get u1,v1,w1 & ThetaR1 ---
		do iz=1,m
			wrk = u1hat(:,:,iz); 
			out = wrk
			in = out
			u(:,:,iz) = in;
		end do
		!------------------------
		!************  RK2 end ***************************
		!****- RK3  start ********************************
		! call RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vt,u1hat,v1hat,w1hat,ThetaR1hat,qt1hat,u,v,w,ThetaR,qt,&
		! fu1hat,fv1hat,fw1hat,fThetaR1hat,fqt1hat,RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)
		fu1hat = u1hat + u ! RKs3

		uhat = uhat*e_nu_3uv + dt/4*fuhat*e_nu_3uv + 3/4*dt*fu1hat*e_nu_1uv;
	
		!******- Possion equation for pressure p -*******
		rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*uhat(:,:,1:m-1)
		phat = 0;
		phat = a+b+c+rhat ! Thomas

		!***- Update veloctiy field ----------------------
		uhat(:,:,1:m-1) = uhat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat;

		!***Neumann boundary condition, what(1)=what(1),what(m)=what(m)***
		!--- Get u1,v1,w1 & ThetaR1 ---
		do iz=1,m
			wrk = u1hat(:,:,iz); 
			out = wrk
			in = out
			u(:,:,iz) = in;
		end do
		!************  RK3 end ***************************
!---------------------------------------------------
		Ti = Ti+dt;
	end do

endsubroutine simplefare

subroutine simple_fare_ad  ! with only u
	implicit none 

	integer :: it, ix, iy, iz

	integer :: nx = 128
	integer :: ny = 128
	integer :: m = 101;
	! integer :: nxh = 65
	! integer :: nxp = 86
	integer ::  nxph = 128!43
	integer :: nyp = 128 !86
	integer :: nyph = 64 !43
	! real :: pi = 3.14

	real :: Us !, Ts, Ths
	real :: Lx, Ly, Lz, dx, dy, dz
	real :: f, Ls, Ts
	! real :: L, vt, tau, qv0
	! real :: a_squall, f
	! integer :: Ls, qs, 
	integer :: nn = 128*128
	real :: Tfinal, IM!, eps, epsbar, Nz, f_star, g_star, B_star
	real :: drx, dry, drz, zk

	!Real
	! real :: max_theta, max_pert
	! real :: xi,yj,da,x_c,y_c,z_c,r_c, ampl_bubble

	real, dimension(128,128,101) :: u ! ,v,w, Theta, ThetaR, qt, qr, qv;
	real, dimension(128,128,101) :: a_u ! ,v,w, Theta, ThetaR, qt, qr, qv;
	! real, dimension(101) :: qvini, RC, Mr, ubg, vbg, dqvdz, qvs, theta_bar, u_bar, v_bar, wrk1D, qvs0
	real, dimension(101) :: u_bar, ubg

	real :: Ti, mindt
	real :: dt0, dt, CFL!, N_theta
	real, dimension(2) :: nu

	real,  dimension(128,128) :: in, wrk, out 
	real,  dimension(101) :: tau_z;
	real,  dimension(128,128,101) :: uhat!, vhat, what, ThetaRhat, qthat
	real,  dimension(128,128,101) :: u1hat!, v1hat, w1hat, ThetaR1hat, qt1hat
	real,  dimension(128,128,101) :: fuhat!, fvhat, fwhat, fThetaRhat, fqthat
	real,  dimension(128,128,101) :: fu1hat!, fv1hat, fw1hat, fThetaR1hat, fqt1hat
	real,  dimension(128,128,101) :: a_uhat!, vhat, what, ThetaRhat, qthat
	real,  dimension(128,128,101) :: a_u1hat!, v1hat, w1hat, ThetaR1hat, qt1hat
	real,  dimension(128,128,101) :: a_fuhat!, fvhat, fwhat, fThetaRhat, fqthat
	real,  dimension(128,128,101) :: a_fu1hat!, fv1hat, fw1hat, fThetaR1hat, fqt1hat

	real,  dimension(128,128,101) :: kx, ky, kk
	!real,  dimension(43,86,101) :: e_nu_1, e_nu_2, e_nu_3
	real,  dimension(128,128,101) :: e_nu_1uv, e_nu_2uv, e_nu_3uv
	!real,  dimension(43,86,101) :: e_nu_1w, e_nu_2w, e_nu_3w
	real,  dimension(128,128,100) :: a, b, c, rhat, phat
	real,  dimension(128,128,100) :: a_rhat, a_phat
	do while ( Ti > 0 )

		!************  RK3 end ***************************
		
		!***Neumann boundary condition, what(1)=what(1),what(m)=what(m)***
        !--- Get u1,v1,w1 & ThetaR1 ---
		do iz = m,1, -1
			a_uhat(:,:,iz) = a_u(:,:,iz)
			! a_in = a_u(:,:,iz)
			! call dfftw_execute_dft_c2r(planr,a_out,a_in); 
			! call padm(nxh,nxph,ny,nyp,nyph,a_wrk,a_out)
			! a_uhat(:,:,iz) = a_wrk
		end do 

		!***- Update veloctiy field ----------------------
		a_phat = a_phat - IM*kx(:,:,1:m-1)*a_uhat

		!******- Possion equation for pressure p -*******
		! a_Thomas: a_phat -> a_rhat
		call a_Thomas(a_phat,a,b,c,a_rhat,nxph,nyp,m-1)
		a_uhat(:,:,1:m-1) = IM*kx(:,:,1:m-1) * a_rhat(:,:,1:m-1)

		a_fu1hat = 3/4 * dt * e_nu_1uv * a_uhat  ! k3: k3_bsb = dt*3*bsb(t+1)/4
		a_u1hat = 0.0   						 ! tmp_wb = 0.0
		! a_RK_flux: a_fu1hat -> a_u1hat
		! call a_RK_flux((nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,
		! 				f_star,g_star,epsbar,L,B_star,nu,vt,
		! 				a_u1hat,a_v1hat, a_w1hat, a_ThetaR1hat, a_qt1hat,
		! 				u, v, w, ThetaR, qt,
		! 				a_fu1hat, a_fv1hat, a_fw1hat, a_fThetaR1hat, a_fqt1hat,
		! 				RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)
		call a_RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,
						f_star,g_star, a_u1hat, u, a_fu1hat, in,out,planf)

		!****- RK3  start ********************************
		!************  RK2 end ***************************

		do iz = m,1, -1
			a_u(:,:,iz) = a_u1hat(:,:,iz)
			! a_in = a_u(:,:,iz)
			! call dfftw_execute_dft_c2r(planr,a_out,a_in); 
			! call padm(nxh,nxph,ny,nyp,nyph,a_wrk,a_out)
			! a_u1hat(:,:,iz) = a_wrk
		end do 

		!***- Update veloctiy field ----------------------
		a_phat = a_phat - IM*kx(:,:,1:m-1)*a_u1hat

		!******- Possion equation for pressure p -*******
		call a_Thomas(a_phat,a,b,c,a_rhat,nxph,nyp,m-1)
		a_u1hat(:,:,1:m-1) = IM*kx(:,:,1:m-1) * a_rhat(:,:,1:m-1)

		a_uhat(:,:,1:m-1) = a_uhat(:,:,1:m-1) + a_u1hat(:,:,1:m-1) !bsb(t) = bsb(t) + bsb(t+1) + tmp_bsb
		a_fu1hat = 2/3*dt *e_nu_1uv * a_u1hat   ! k2:  k2_bsb = dt*2*tmp_bsb/3
		a_u1hat = 0								! tmp_wb = 0.0

		! a_RK_flux: a_fu1hat -> a_u1hat
		call a_RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star, a_u1hat, u, a_fu1hat, in,out,planf)
		!****- RK2  start ********************************
		!************  RK1 end ***********************************

		!***Neumann boundary condition, w1hat(1)=w1hat(1),w1hat(m)=w1hat(m)***
        !--- Get u1,v1,w1 & ThetaR1 ---
		do iz = m,1, -1
			a_u(:,:,iz) = a_u1hat(:,:,iz)
			! a_in = a_u(:,:,iz)
			! call dfftw_execute_dft_c2r(planr,a_out,a_in); 
			! call padm(nxh,nxph,ny,nyp,nyph,a_wrk,a_out)
			! a_u1hat(:,:,iz) = a_wrk
		end do 

		!***- Update veloctiy field ----------------------
		a_phat = a_phat - IM*kx(:,:,1:m-1)*a_u1hat

		!******- Possion equation for pressure p -*******
		! a_Thomas: a_phat -> a_rhat
		call a_Thomas(a_phat,a,b,c,a_rhat,nxph,nyp,m-1)
		a_u1hat(:,:,1:m-1) = IM*kx(:,:,1:m-1) * a_rhat(:,:,1:m-1)

		a_fuhat = dt/3*e_nu_1uv * a_u1hat + dt/4*e_nu_3uv * a_uhat ! k1: k1_bsb = dt*bsb(t+1)/4 + dt*tmp_bsb/3

		a_uhat = a_uhat + a_u1hat  									! bsb(t) = bsb(t) + tmp_bsb
		! a_RK_flux: a_fuhat -> a_uhat
		call a_RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star, a_uhat, u, a_fuhat, in,out,planf)
		!****- RK1  start ********************************

		Ti = Ti - dt
	end do
endsubroutine






