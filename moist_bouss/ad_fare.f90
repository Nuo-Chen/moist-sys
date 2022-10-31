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

		do iz = m,1, -1
			a_wrk = a_ThetaRhat(:,:,iz)/nn
			call unpadm(nxh,nxph,ny,nyp,nyph,a_wrk,a_out)
			call dfftw_execute_dft_r2c(planf,a_out,a_in);
			a_ThetaR(:,:,iz) = a_in 

			a_wrk = a_qthat(:,:,iz)/nn
			call unpadm(nxh,nxph,ny,nyp,nyph,a_wrk,a_out)
			call dfftw_execute_dft_r2c(planf,a_out,a_in);
			a_qt(:,:,iz) = a_in 
		end do 

		a_thetaR = a_thetaR + a_theta
		a_qr = a_qr + L * a_theta + a_qt 
		a_qv = a_qv + a_qt 
		!a_qt = 0
		do iz = m, 1, -1
			where (mask0(:,:)) a_qv(:,:,iz) = 0
			call POPBOOLEANARRAY(mask0, 128**2)
			a_qt(:,:,iz) = a_qt(:,:,iz) + a_qv(:,:,iz)
			a_qr(:,:,iz) = a_qr(:,:,iz) - a_qv(:,:,iz)
			!a_qv(:,:,iz) = 0

			mask(:,:) = qvini(iz) + qt(:,:,iz) .gt. qvs(iz)
			where (mask(:,:))
				a_qt(:,:,iz) = a_qt(:,:,iz) + a_qr(:,:,iz)
				!a_qr(:,:,iz) = 0
			end where 
		end do
		!a_qr = 0
		!a_theta = 0

		!************  RK3 end ***************************
		
		!***Neumann boundary condition, what(1)=what(1),what(m)=what(m)***
        !--- Get u1,v1,w1 & ThetaR1 ---
		do iz = m,1, -1
			! a_u(:,:,iz) = a_uhat(:,:,iz)
			a_in = a_u(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_uhat(:,:,iz) = a_wrk

			a_in = a_v(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_vhat(:,:,iz) = a_wrk

			a_in = a_w(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_what(:,:,iz) = a_wrk

			a_in = a_qt(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_qthat(:,:,iz) = a_wrk

			a_in = a_ThetaR(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_ThetaRhat(:,:,iz) = a_wrk
		end do 


		!***- Update veloctiy field ----------------------
		a_what(:, :, m) = 0.0
		a_what(:, :, 1) = 0.0
		a_phat(:, :, 2:m-1) = a_phat(:, :, 2:m-1) - a_what(:, :, 2:m-1)/dz
		a_phat(:, :, 1:m-2) = a_phat(:, :, 1:m-2) + a_what(:, :, 2:m-1)/dz
		a_phat = a_phat - IM*ky(:, :, 1:m-1) * a_vhat(:, :, 1:m-1) - IM*kx(:, :, 1: m-1) * a_uhat(:, :, 1:m-1)

		!******- Possion equation for pressure p -*******
		! a_Thomas: a_phat -> a_rhat
		call a_Thomas(a_phat,a,b,c,a_rhat,nxph,nyp,m-1)

		a_uhat(:, :, 1:m-1) = a_uhat(:, :, 1:m-1) + IM*kx(:, :, 1:m-1) * a_rhat(:, :, 1:m-1)
		a_vhat(:, :, 1:m-1) = a_vhat(:, :, 1:m-1) + IM*ky(:, :, 1:m-1) * a_rhat(:, :, 1:m-1)
		a_what(:, :, 2:m) =   a_what(:, :, 2:m) + a_rhat(:, :, 1:m-1)/dz
		a_what(:, :, 1:m-1) = a_what(:, :, 1:m-1) - a_rhat(:, :, 1:m-1)/dz
		!a_rhat(:, :, 1:m-1) = 0.0

		a_fqthat = a_fqthat + dt*e_nu_3/4 * a_qthat
		a_fqt1hat = a_fqt1hat + dt*3*e_nu_1/4 * a_qthat
		a_qthat = e_nu_3 * a_qthat
		a_fThetaRhat = a_fThetaRhat + dt*e_nu_3/4 * a_ThetaRhat
		a_fThetaR1hat = a_fThetaR1hat + dt*3*e_nu_1/4 * a_ThetaRhat
		a_ThetaRhat = e_nu_3 * a_ThetaRhat
		a_fwhat = a_fwhat + dt*e_nu_3w/4 * a_what
		a_fw1hat = a_fw1hat + dt*3*e_nu_1w/ 4 * a_what
		a_what = e_nu_3w * a_what
		a_fvhat = a_fvhat + dt*e_nu_3uv/4 * a_vhat
		a_fv1hat = a_fv1hat + dt*3*e_nu_1uv/4 * a_vhat
		a_vhat = e_nu_3uv * a_vhat
		a_fu1hat = a_fu1hat + dt*3*e_nu_1uv/4 * a_uhat   ! k3: k3_bsb = dt*3*bsb(t+1)/4
		a_uhat = e_nu_3uv * a_uhat

		a_u1hat = 0.0   						 ! tmp_wb = 0.0
		a_v1hat = 0.0   
		a_w1hat = 0.0   
		a_qt1hat = 0.0   
		a_ThetaR1hat = 0.0   
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
			a_in = a_u(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out)
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_u1hat(:,:,iz) = a_wrk

			a_in = a_v(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out)
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_v1hat(:,:,iz) = a_wrk

			a_in = a_w(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out)
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_w1hat(:,:,iz) = a_wrk

			a_in = a_qt(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out)
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_qt1hat(:,:,iz) = a_wrk

			a_in = a_ThetaR(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out)
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_ThetaR1hat(:,:,iz) = a_wrk
		end do 

		!***- Update veloctiy field ---------------------

		a_w1hat(:, :, m) = 0.0
		a_w1hat(:, :, 1) = 0.0
		a_phat(:, :, 2:m-1) = a_phat(:, :, 2:m-1) - a_w1hat(:, :, 2:m-1)/dz
		a_phat(:, :, 1:m-2) = a_phat(:, :, 1:m-2) + a_w1hat(:, :, 2:m-1)/dz
		a_phat = a_phat - IM*ky(:, :, 1:m-1) * a_v1hat(:, :, 1:m-1) - IM*kx(:, :, 1: m-1) * a_u1hat(:, :, 1:m-1)


		!******- Possion equation for pressure p -*******
		call a_Thomas(a_phat,a,b,c,a_rhat,nxph,nyp,m-1)
		a_u1hat(:, :, 1:m-1) = a_u1hat(:, :, 1:m-1) + IM*kx(:, :, 1:m-1) * a_rhat(:, :, 1:m-1)
		a_v1hat(:, :, 1:m-1) = a_v1hat(:, :, 1:m-1) + IM*ky(:, :, 1:m-1) * a_rhat(:, :, 1:m-1)
		a_w1hat(:, :, 2:m) =   a_w1hat(:, :, 2:m) + a_rhat(:, :, 1:m-1)/dz
		a_w1hat(:, :, 1:m-1) = a_w1hat(:, :, 1:m-1) - a_rhat(:, :, 1:m-1)/dz
		!a_rhat(:, :, 1:m-1) = 0.0

		a_uhat(:,:,1:m-1) = a_uhat(:,:,1:m-1) + a_u1hat(:,:,1:m-1) !bsb(t) = bsb(t) + bsb(t+1) + tmp_bsb
		a_vhat(:,:,1:m-1) = a_vhat(:,:,1:m-1) + a_v1hat(:,:,1:m-1) 
		a_what(:,:,1:m-1) = a_what(:,:,1:m-1) + a_w1hat(:,:,1:m-1) 
		! a_qthat(:,:,1:m-1) = a_qthat(:,:,1:m-1) + a_qt1hat(:,:,1:m-1) 
		! a_ThetaRhat(:,:,1:m-1) = a_ThetaRhat(:,:,1:m-1) + a_ThetaR1hat(:,:,1:m-1) 
		
		!******- Possion equation for pressure p -*******
		! a_Thomas: a_phat -> a_rhat
		call a_Thomas(a_phat,a,b,c,a_rhat,nxph,nyp,m-1)

		a_qthat = a_qthat + e_nu_2* a_qt1hat   ! k2:  k2_bsb = dt*2*tmp_bsb/3
		a_fqt1hat = a_fqt1hat + dt*2*e_nu_1/3 * a_qt1hat
		a_ThetaRhat = a_ThetaRhat + e_nu_2 * a_ThetaR1hat
		a_fThetaR1hat = a_fThetaR1hat + dt*2*e_nu_1/3 * a_ThetaR1hat
		a_what = a_what + e_nu_2w * a_w1hat
		a_fw1hat = a_fw1hat + dt*2*e_nu_1w/3 * a_w1hat
		a_vhat = a_vhat + e_nu_2uv * a_v1hat
		a_fv1hat = fa_v1hat + dt*2*e_nu_1uv/3 * a_v1hat
		a_uhat = a_uhat + e_nu_2uv * a_u1hat
		a_fu1hat = a_fu1hat + dt*2*e_nu_1uv/3 * a_u1hat
		
		a_v1hat = 0.0  ! tmp_wb = 0.0
		a_w1hat = 0.0
		a_qt1hat = 0.0
		a_ThetaR1hat = 0.0
		a_u1hat = 0.0							

		! a_RK_flux: a_fu1hat -> a_u1hat
		call a_RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star, a_u1hat, u, a_fu1hat, in,out,planf)
		!****- RK2  start ********************************
		!************  RK1 end ***********************************

		!***Neumann boundary condition, w1hat(1)=w1hat(1),w1hat(m)=w1hat(m)***
        !--- Get u1,v1,w1 & ThetaR1 ---
		do iz = m,1, -1
			! a_u(:,:,iz) = a_u1hat(:,:,iz)
			a_in = a_u(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_u1hat(:,:,iz) = a_wrk

			a_in = a_v(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_v1hat(:,:,iz) = a_wrk

			a_in = a_w(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_w1hat(:,:,iz) = a_wrk

			a_in = a_qt(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_qt1hat(:,:,iz) = a_wrk

			a_in = a_ThetaR(:,:,iz)
			call dfftw_execute_dft_c2r(planr,a_in,a_out); 
			call padm(nxh,nxph,ny,nyp,nyph,a_out,a_wrk)
			a_ThetaR1hat(:,:,iz) = a_wrk
		end do 

		!***- Update veloctiy field ----------------------

		a_w1hat(:, :, m) = 0.0
		a_w1hat(:, :, 1) = 0.0
		a_phat(:, :, 2:m-1) = a_phat(:, :, 2:m-1) - a_w1hat(:, :, 2:m-1)/dz
		a_phat(:, :, 1:m-2) = a_phat(:, :, 1:m-2) + a_w1hat(:, :, 2:m-1)/dz
		a_phat = a_phat - IM*ky(:, :, 1:m-1) * a_v1hat(:, :, 1:m-1) - IM*kx(:, :, 1: m-1) * a_u1hat(:, :, 1:m-1)
		
		!******- Possion equation for pressure p -*******
		! a_Thomas: a_phat -> a_rhat
		call a_Thomas(a_phat,a,b,c,a_rhat,nxph,nyp,m-1)

		a_u1hat(:, :, 1:m-1) = a_u1hat(:, :, 1:m-1) + IM*kx(:, :, 1:m-1) * a_rhat(:, :, 1:m-1)
		a_v1hat(:, :, 1:m-1) = a_v1hat(:, :, 1:m-1) + IM*ky(:, :, 1:m-1) * a_rhat(:, :, 1:m-1)
		a_w1hat(:, :, 2:m) =   a_w1hat(:, :, 2:m) + a_rhat(:, :, 1:m-1)/dz
		a_w1hat(:, :, 1:m-1) = a_w1hat(:, :, 1:m-1) - a_rhat(:, :, 1:m-1)/dz
		!a_rhat(:, :, 1:m-1) = 0.0

		a_qthat = a_qthat + e_nu_1 * a_qt1hat  ! k1: k1_bsb = dt*bsb(t+1)/4 + dt*tmp_bsb/3
		a_fqthat = a_fqthat + dt*e_nu_1/3 * a_qt1hat
		a_ThetaRhat = a_ThetaRhat + e_nu_1 * a_ThetaR1hat
		a_fThetaRhat = a_fThetaRhat + dt*e_nu_1/3 * a_ThetaR1hat
		a_what = a_what + e_nu_1w * a_w1hat
		a_fwhat = a_fwhat + dt*e_nu_1w/3 * a_w1hat
		a_vhat = a_vhat + e_nu_1uv * a_v1hat
		a_fvhat = a_fvhat + dt*e_nu_1uv/3 * a_v1hat
		a_uhat = a_uhat + e_nu_1uv * a_u1hat
		a_fuhat = a_fuhat + dt*e_nu_1uv3 * a_u1hat

		a_uhat = a_uhat + a_u1hat  									! bsb(t) = bsb(t) + tmp_bsb
		a_vhat = a_vhat + a_v1hat
		a_what = a_what + a_w1hat
		a_qthat = a_qthat + a_qt1hat
		a_ThetaRhat = a_ThetaRhat + a_ThetaR1hat
		
		! a_RK_flux: a_fuhat -> a_uhat
		call a_RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star, a_uhat, u, a_fuhat, in,out,planf)
		!****- RK1  start ********************************

		Ti = Ti - dt
	end do
endsubroutine





subroutine update_q(Theta, ThetaR, qt, qr, qv)
	implicit none
	real, dimension(128,128,101) :: Theta, ThetaR, qt, qr, qv
	real, dimension(101) :: qvini, qvs
	real :: L
	integer :: iz, m
	logical, dimension(128, 128) :: mask, mask0
	qvini = 0
	qvs = 0
	qr = 0
	m = 101
	L = 1005
	do iz=1,m
		! where ( qvini(iz)+qt(:,:,iz) > qvs(iz) ) qr(:,:,iz) = qvini(iz)+qt(:,:,iz)-qvs(iz);
		! qv(:,:,iz) = qt(:,:,iz)-qr(:,:,iz);
		! where ( qv(:,:,iz)+qvini(iz) < 0 ) qv(:,:,iz) = -qvini(iz);
		mask(:, :) = qvini(iz) + qt(:,:,iz) .gt. qvs(iz)
		where (mask(:,:))  qr(:, :, iz) = qvini(iz) + qt(:, :, iz) - qvs(iz)
		qv(:,:,iz) = qt(:,:,iz)-qr(:,:,iz)
		call PUSHBOOLEANARRAY(mask0, 128**2)
		mask0(:, :) = qv(:, :, iz) + qvini(iz) .LT. 0
	end do
	qt = qv+qr
	Theta = ThetaR+L*qr

	a_thetaR = a_thetaR + a_theta
	a_qr = a_qr + L * a_theta + a_qt 
	a_qv = a_qv + a_qt 
	a_qt = 0
	do iz = m, 1, -1
		where (mask0(:,:)) a_qv(:,:,iz) = 0
		call POPBOOLEANARRAY(mask0, 128**2)
		a_qt(:,:,iz) = a_qt(:,:,iz) + a_qv(:,:,iz)
		a_qr(:,:,iz) = a_qr(:,:,iz) - a_qv(:,:,iz)
		a_qv(:,:,iz) = 0

		mask(:,:) = qvini(iz) + qt(:,:,iz) .gt. qvs(iz)
		where (mask(:,:))
			a_qt(:,:,iz) = a_qt(:,:,iz) + a_qr(:,:,iz)
			a_qr(:,:,iz) = 0
		end where 
	end do
	a_qr = 0
	a_theta = 0
		
endsubroutine update_q

subroutine update_w(uhat, vhat, what, phat, ThetaRhat, qthat, uhat, fvhat, fwhat, fThetaRhat, fqthat, fu1hat, fv1hat, fw1hat, fThetaR1hat, fqt1hat)
	implicit none
	real, dimension(128,128,101) :: uhat, vhat, what, phat, ThetaRhat, qthat, rhat
	real, dimension(128,128,101) :: fuhat, fvhat, fwhat, fThetaRhat, fqthat
	real, dimension(128,128,101) :: fu1hat, fv1hat, fw1hat, fThetaR1hat, fqt1hat
	real, dimension(128,128,101) :: u1hat, v1hat, w1hat, ThetaR1hat, qt1hat

	integer, dimension(128, 128, 101) :: kx, ky
	integer :: IM, dz, m, dt
	integer :: e_nu_3uv, e_nu_3w, e_nu_3
	integer :: e_nu_1uv, e_nu_1w, e_nu_1

	kx = 8
	ky = 8
	IM = 1
	dz = 10
	m = 101
	dt = 1
	e_nu_3uv = e_nu_3w = e_nu_3 = 1 
	e_nu_1uv = e_nu_1w = e_nu_1 = 1
	e_nu_2uv = e_nu_2w = e_nu_2 = 1
	! rk 1
	u1hat = uhat*e_nu_1uv + dt/3*fuhat*e_nu_1uv
	v1hat = vhat*e_nu_1uv + dt/3*fvhat*e_nu_1uv
	w1hat = what*e_nu_1w + dt/3*fwhat*e_nu_1w		
	ThetaR1hat = ThetaRhat*e_nu_1 + dt/3*fThetaRhat*e_nu_1
	qt1hat = qthat*e_nu_1 + dt/3*fqthat*e_nu_1
	! rk2
	u1hat = uhat*e_nu_2uv + 2/3*dt*fu1hat*e_nu_1uv;
	v1hat = vhat*e_nu_2uv + 2/3*dt*fv1hat*e_nu_1uv;
	w1hat = what*e_nu_2w + 2/3*dt*fw1hat*e_nu_1w;
	ThetaR1hat = ThetaRhat*e_nu_2 + 2/3*dt*fThetaR1hat*e_nu_1;
	qt1hat = qthat*e_nu_2 + 2/3*dt*fqt1hat*e_nu_1;
	! rk3
	uhat = uhat*e_nu_3uv + dt/4*fuhat*e_nu_3uv + 3/4*dt*fu1hat*e_nu_1uv;
	vhat = vhat*e_nu_3uv + dt/4*fvhat*e_nu_3uv + 3/4*dt*fv1hat*e_nu_1uv;
	what = what*e_nu_3w + dt/4*fwhat*e_nu_3w + 3/4*dt*fw1hat*e_nu_1w;
	ThetaRhat = ThetaRhat*e_nu_3 + dt/4*fThetaRhat*e_nu_3 + 3/4*dt*fThetaR1hat*e_nu_1;
	qthat = qthat*e_nu_3 + dt/4*fqthat*e_nu_3 + 3/4*dt*fqt1hat*e_nu_1;

	!******- Possion equation for pressure p -*******
	rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*uhat(:,:,1:m-1)+IM*ky(:,:,1:m-1)*vhat(:,:,1:m-1)+(what(:,:,2:m)-what(:,:,1:m-1))/dz;
	phat = 0;

	phat = rhat ! Thomas

	uhat(:,:,1:m-1) = uhat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat
	vhat(:,:,1:m-1) = vhat(:,:,1:m-1) - IM*ky(:,:,1:m-1)*phat
	what(:,:,2:m-1) = what(:,:,2:m-1) - (phat(:,:,2:m-1)-phat(:,:,1:m-2))/dz
	what(:,:,1) = 0
	what(:,:,m) = 0

	qthat = qthat + e_nu_1*qt1hat
	fqthat = fqthat + dt*e_nu_1*qt1hat/3
	ThetaRhat = ThetaRhat + e_nu_1*ThetaR1hat
	fThetaRhat = fThetaRhat + dt*e_nu_1*ThetaR1hat/3
	what = what + e_nu_1w*w1hat
	fwhat = fwhat + dt*e_nu_1w*w1hat/3
	vhat = vhat + e_nu_1uv*v1hat
	fvhat = fvhat + dt*e_nu_1uv*v1hat/3
	uhat = uhat + e_nu_1uv*u1hat
	fuhat = fuhat + dt*e_nu_1uv*u1hat/3
	v1hat = 0.0
	w1hat = 0.0
	qt1hat = 0.0
	ThetaR1hat = 0.0
	u1hat = 0.0

	qthat = qthat + e_nu_2*qt1hat
	fqt1hat = fqt1hat + dt*2*e_nu_1*qt1hat/3
	ThetaRhat = ThetaRhat + e_nu_2*ThetaR1hat
	fThetaR1hat = fThetaR1hat + dt*2*e_nu_1*ThetaR1hat/3
	what = what + e_nu_2w*w1hat
	fw1hat = fw1hat + dt*2*e_nu_1w*w1hat/3
	vhat = vhat + e_nu_2uv*v1hat
	fv1hat = fv1hat + dt*2*e_nu_1uv*v1hat/3
	uhat = uhat + e_nu_2uv*u1hat
	fu1hat = fu1hat + dt*2*e_nu_1uv*u1hat/3
	v1hat = 0.0
	w1hat = 0.0
	qt1hat = 0.0
	ThetaR1hat = 0.0
	u1hat = 0.0

	fqthat = fqthat + dt*e_nu_3*qthat/4
	fqt1hat = fqt1hat + dt*3*e_nu_1*qthat/4
	qthat = e_nu_3*qthat
	fThetaRhat = fThetaRhat + dt*e_nu_3*ThetaRhat/4
	fThetaR1hat = fThetaR1hat + dt*3*e_nu_1*ThetaRhat/4
	ThetaRhat = e_nu_3*ThetaRhat
	fwhat = fwhat + dt*e_nu_3w*what/4
	fw1hat = fw1hat + dt*3*e_nu_1w*what/4
	what = e_nu_3w*what
	fvhat = fvhat + dt*e_nu_3uv*vhat/4
	fv1hat = fv1hat + dt*3*e_nu_1uv*vhat/4
	vhat = e_nu_3uv*vhat
	fu1hat = fu1hat + dt*3*e_nu_1uv*uhat/4
	uhat = e_nu_3uv*uhat

endsubroutine