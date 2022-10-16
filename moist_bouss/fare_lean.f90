
!ifort lib/*.f90 FARE.f90 -o FARE -lpthread -lfftw3_threads -lfftw3 -lm
!ifort lib/*.f90 FARE.f90 -o FARE -L/usr/local/lib -lfftw3 -lfftw3_threads -lm -C 

program	FARE
    use, intrinsic :: iso_c_binding
    use FFTW3
    use mainsubroutines
    use mainvars
    use mpads
    implicit none
    interface
        !This is qv(z) or qvs (z)
        pure subroutine FQv(qv,m,dz,qv0)
        use, intrinsic :: iso_c_binding
        implicit none
            integer(C_INT) :: ix,iy,iz
            real(C_DOUBLE) :: zk,k0,k1,k2,k3,k4, qvs, pz
            integer(C_INT), intent(in) :: m
            real(C_DOUBLE), intent(in) :: dz, qv0
            real(C_DOUBLE), dimension(m),intent(out) :: qv
        end subroutine FQv
        !This is dqv/dz for the initial data
        pure subroutine Fdqvdz(dqvdz,m,dz,qv0)
        use, intrinsic :: iso_c_binding
        implicit none
            integer(C_INT) :: iz
            real(C_DOUBLE) :: k0,z_k,k1,k2,k3,k4, qvs, pz
            integer(C_INT), intent(in) :: m
            real(C_DOUBLE), intent(in) :: dz, qv0
            real(C_DOUBLE), dimension(m),intent(out) :: dqvdz
        end subroutine Fdqvdz
    end interface
    
        !Integer
        integer(C_INT) :: it
        !Real
        real(C_DOUBLE) :: max_theta, max_pert
        real(C_DOUBLE) :: xi,yj,da,x_c,y_c,z_c,r_c, ampl_bubble
        !Allocatable
        !Real
        real(C_DOUBLE) ,dimension(:), allocatable :: theta_bar, u_bar, v_bar, wrk1D
    !----------------PARAMETERS----------------------------------------------
    
        !Grid size:
        nx  = 128; ny  = 128; m  = 100+1;  ! m: vertial levels
        nxh = nx/2+1; nxp = nx/3; nxp = 2*nxp+1; nxph = (nxp-1)/2+1;
        nyp = ny/3; nyp = 2*nyp+1; nyph = (nyp-1)/2+1;
    
        !Quantity scales:
        Us = 100d0/9d0; !In m/s
        Ts = 15d0/60d0; !In hours
        Ths = 3.0d0; !In Kelvin
        Ls = 10d0; !In kms
        qs = 10d0; !In g/kg
    
        !*******DOMAIN SIZE*******
        Lx = 128d0; Ly = 128d0; Lz = 15d0;
        Lx = Lx/Ls; Ly = Ly/Ls; Lz = Lz/Ls; !Non-dimentionalization
        dx = Lx/nx; dy = Ly/ny; dz = Lz/(m-1);
    
        !Rain fall velocity
        vt = 0d0;
        vt = vt/Us; !Non-dimentionalization
    
        !Latent heat
        L = 2.5d0; !In 10^{6}*J*kg^{-1}
        L = L/3d0*10d0; !Non-dimentionalization prefector epsilon^{-1} L = epsilon^2 * L^d * \theta_0 ([\Theta]*c_p*T_0)^{-1}
        !Pendiente: Cambiar a LcpTheta0
    
        !Relaxation time   for forcing in squall line, term in momentum back to velocity from jet, or remoist the model lower boundary
        ! tau = 4d0; ! In hours
        ! tau = tau/Ts; !Non-dimensionalization
        tau = 0d0; !Turned off
    
        !Qvs at surface
        qvs0 = 20d0; !In g/kg
        qvs0 = qvs0/qs; !Non-dimensionalization 
        
        qv0 = 18d0; !In g/kg
        qv0 = qv0/qs;
    
        !Max perturbation in theta
        max_pert = 0.1d0*0d0; !In Kelvin
        max_pert = max_pert/Ths; !Non-dimensionalization
    
        !Parameter for the zonal vel background
        ! a_squall = 10d0; !In m/s
        ! a_squall = a_squall/Us;
        a_squall = a_squall*0d0; !Turned off
    
        !Coriolis parameter
        f = sin(25d0*pi/180d0);
    
        !Final time
        Tfinal = 5d0; !In hours
        Tfinal = Tfinal/Ts; !Non-dimensionalization
    
        !***************
        IM = (0d0,1d0)
    
        !***************
        eps = 0.1d0;
        epsbar = 0.622d0;
        Nz = 1d0;
    
        !***************
        
        f_star = eps*f;
        g_star = 10d0; !7.9461d0;
        B_star = 10d0;
        
        !Wave numbers are multiple of those
        drx = 2d0*pi/Lx; dry = 2d0*pi/Ly; drz = 2d0*pi/Lz;

        !!!!Data for Bubble:
        x_c = Lx/2d0; y_c = Ly/2d0; z_c = 2d0/Ls; 
        !rx_c = 10d0; ry_c = 10d0; rz_c = 1;
        !rx_c = rx_c/Ls; ry_c = ry_c/Ls; rz_c = rz_c/Ls; 
        ampl_bubble = 8d0/4d0/qs;   ! gaussian shaped bubble, how strong the heating/moistening inside
        r_c = 5d0;
        r_c = r_c/Ls;
    !-------Allocation----------------------
    
        ! !Real
        ! allocate(u(nx,ny,m)); allocate(v(nx,ny,m)); allocate(w(nx,ny,m)); allocate(Theta(nx,ny,m));	
        ! allocate(ThetaR(nx,ny,m)); allocate(qv(nx,ny,m)); allocate(qr(nx,ny,m)); allocate(qt(nx,ny,m));
        ! allocate(qvini(m)); allocate(RC(m)); allocate(Mr(m)); allocate(ubg(m)); allocate(vbg(m)); 
        ! allocate(dqvdz(m)); allocate(qvs(m));
        ! allocate(theta_bar(m)); allocate(u_bar(m)); allocate(v_bar(m));
        ! allocate(wrk1D(m));
        
        ! !Complex
        ! allocate(in(nx,ny)); allocate(out(nxh,ny));
        ! allocate(wrk(nxph,nyp)); allocate(tau_z(m));
        ! allocate(uhat(nxph,nyp,m)); allocate(vhat(nxph,nyp,m)); allocate(what(nxph,nyp,m)); 
        ! allocate(ThetaRhat(nxph,nyp,m)); allocate(qthat(nxph,nyp,m));
        ! allocate(u1hat(nxph,nyp,m)); allocate(v1hat(nxph,nyp,m)); allocate(w1hat(nxph,nyp,m)); 
        ! allocate(ThetaR1hat(nxph,nyp,m)); allocate(qt1hat(nxph,nyp,m));
        ! allocate(fuhat(nxph,nyp,m)); allocate(fvhat(nxph,nyp,m)); allocate(fwhat(nxph,nyp,m)); 
        ! allocate(fThetaRhat(nxph,nyp,m)); allocate(fqthat(nxph,nyp,m));
        ! allocate(fu1hat(nxph,nyp,m)); allocate(fv1hat(nxph,nyp,m)); allocate(fw1hat(nxph,nyp,m)); 
        ! allocate(fThetaR1hat(nxph,nyp,m)); allocate(fqt1hat(nxph,nyp,m));
        ! allocate(kx(nxph,nyp,m)); allocate(ky(nxph,nyp,m)); allocate(kk(nxph,nyp,m));
        ! allocate(e_nu_1(nxph,nyp,m)); allocate(e_nu_2(nxph,nyp,m)); allocate(e_nu_3(nxph,nyp,m));
        ! allocate(e_nu_1uv(nxph,nyp,m)); allocate(e_nu_2uv(nxph,nyp,m)); allocate(e_nu_3uv(nxph,nyp,m));
        ! allocate(e_nu_1w(nxph,nyp,m)); allocate(e_nu_2w(nxph,nyp,m)); allocate(e_nu_3w(nxph,nyp,m));
        ! allocate(a(nxph,nyp,m-1));  allocate(b(nxph,nyp,m-1));  allocate(c(nxph,nyp,m-1)); 
        ! allocate(rhat(nxph,nyp,m-1)); allocate(phat(nxph,nyp,m-1));
    
    !-------------------FFTW PACKAGE-----------------------------------------
    
        call dfftw_init_threads(iret)
        !print*, 'iret=',iret
        nthreads = 4;
          call dfftw_plan_with_nthreads(nthreads)
    
        planf = fftw_plan_dft_r2c_2d(nx,ny,in,out,FFTW_ESTIMATE); !FFTW_PATIENT); !FFTW_MEASURE); !
        planr = fftw_plan_dft_c2r_2d(nx,ny,out,in,FFTW_ESTIMATE); !FFTW_PATIENT); ! !FFTW_MEASURE); !
    
    !---------------GENERATING RANDOM NUMBERS FOR THE PERTURBATION-----------------------------------
    
    
    !---------------FUNCTIONS-----------------------------------
    
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
        
        call FQv(qvs,m,dz,qvs0);
        call FQv(qvini,m,dz,qv0);
        call Fdqvdz(dqvdz,m,dz,qv0);
        ! qvs, qvini, dqvdz are constants?
        
        Theta = 0d0;
    
        !Define the radiactive cooling rate
        RC = 0d0; !Turned off
        !do iz=1,m
        !	zk = (iz-0.5d0)*dz;
        !	RC(iz) = -exp(-0.2d0*zk)*zk*(1.5d0-zk)**2d0;
        !	RC(iz) = RC(iz)*(1d0/0.45d0)*(1d0/288d0)*25d0; !Profile similar to Majda-Xing, max about 25 K/day
        !end do 
        
        !Define the moistening proportional to dqvdz
        Mr = 0d0; !Turned off
        !do iz=1,m
        !	Mr(iz) = 0.1d0*RC(iz)*dqvdz(iz)*0d0; !Turned off
        !end do
    
        !!!!!!!!!!!!!!!!!!----Damping in w: ----
        tau_z = 0d0;
        do iz=1,m
            zk=(iz-1d0)*dz;
            tau_z(iz) = (1d0/(1d0/Ts))*exp( 50d0*( (1d0/Lz)**2d0-(1d0/(zk-Lz))**2d0 ) );
        end do
        tau_z = 0d0; !Turned off
        
        ubg = 0d0; vbg = 0d0; !Turned off
        !do iz=1,m
        !	zk = (iz-0.5d0)*dz;
        !	if ( zk < 1.2d0 ) then
        !		ubg(iz) = a_squall*(cos(pi*zk/1.2d0)-cos(2d0*pi*zk/1.2d0));
        !	else
        !		ubg(iz) = -2d0*a_squall;
        !	end if
        !
        !	!if ( zk .ge. 1.1d0 ) then
        !	!	vbg(iz) = 4d0/Us;
        !	!else if ( zk .ge. 0.5d0 ) then
        !	!	vbg(iz)=(4d0/Us)*(zk-0.5d0)/(1.1d0-0.5d0);
        !	!end if
        !end do
    
        !!Smoothing out
        !do it = 1,10
        !	v_bar(1) = vbg(2);
        !	do iz=1,m
        !		v_bar(iz) = (vbg(iz-1)+vbg(iz)+vbg(iz+1))/3d0;
        !	end do
        !	v_bar(m) = vbg(m-1);
        !	vbg = v_bar;
        !end do
        
        !!!!!!!!!!!!!!!Pressure solver !!!!
        !Three diagonals for Thomas' method
        !Main diagonal (entries from 1 to m-1), to compute variables in the staggered grid, level at p
        b = - ( kk + 1d0/dz**2d0);
        b(:,:,2:m-2) = b(:,:,2:m-2) - 1d0/dz**2d0;
        b(1,1,1) = 1d0/dz**2d0; ! 1d0; !Pendiente
        b(1,1,m-1) = -1d0/dz**2d0;
        !Lower diagonal, a_2, a_3, ... a_{m-1}
        a = 1d0/dz**2d0;
        a(:,:,1) = 0d0;
        !Upper diagonal: c_1, c_2 ... c_{m-2}
        c = 1d0/dz**2d0;
        c(:,:,m-1) = 0d0;
        c(1,1,1) = 0d0;
    
    !--------------INITIAL DATA-----------------------------------------------
    
        u = 0d0; v = 0d0; w = 0d0; Theta = 0d0; qv = 0d0; qr = 0d0;
        !do iz=1,m
        !	u(:,:,iz) = ubg(iz);
        !	v(:,:,iz) = vbg(iz);
        !end do
        
        !---------Perturbation in temperature-------------------
        Theta = 0d0; !Turned off
        !do iz=1,m-1
        !	zk = (iz-0.5d0)*dz;
        !	do iy=1,ny
        !		do ix=1,nx
        !			!Low level perturbation of theta +- 0.1K
        !			if ( zk >= 0d0 .and. zk <= 0.2d0 ) then
        !				Theta(ix,iy,iz) = Theta(ix,iy,iz)+2d0*(pertrand( (ix-1)*ny*(m-1)+(iy-1)*(m-1)+iz )-0.5d0)*max_pert;
        !			end if
        !		end do	
        !	end do
        !end do
    
        
        !-------Moist bubble
        qt = 0d0;
        do iy = 1,ny
            yj = (iy-1)*dy;
            do ix = 1,nx
                xi = (ix-1)*dx;
                da = sqrt((xi-x_c)**2d0+(yj-y_c)**2d0);
                do iz=1,m-1			
                    zk = (iz-0.5)*dz; !For level at u
                    if ( da <= r_c .and. zk >= z_c-0.1d0 .and. zk <= z_c+0.1d0 ) then		
                        qt(ix,iy,iz) = ampl_bubble*cos(pi*da/(2d0*r_c))*(10d0*(zk-z_c-0.1d0))**2d0*(10d0*(zk-z_c+0.1d0))**2d0; !*(10d0*(zk-z_c-0.3d0))**2d0;
                        !if ( qvini(iz)+qv(ix,iy,iz) > qvs(iz) ) then
                        !	print*,'Reduce amplitude'
                        !	stop
                        !end if
                    end if			
                end do	
            end do
        end do
        !print*,maxval(qv(:,ny/2,:))*qs, ampl_bubble*qs, r_c*Ls, pi
        !stop
        
        qr = 0d0;
        do iz=1,m
            where ( qvini(iz)+qt(:,:,iz) > qvs(iz) ) qr(:,:,iz) = qvini(iz)+qt(:,:,iz)-qvs(iz);
            qv(:,:,iz) = qt(:,:,iz)-qr(:,:,iz);
            where ( qv(:,:,iz)+qvini(iz) < 0d0 ) qv(:,:,iz) = -qvini(iz); !Positivity of qvini+qv
        end do
        qt = qv+qr;
        ThetaR = Theta-L*qr;
    
    !----Pass to Fourier space---------------------------------------------------------
    
        nn = nx*ny
        do iz=1,m
            in = u(:,:,iz); 
            call dfftw_execute_dft_r2c(planf,in,out); 
            call unpadm(nxh,nxph,ny,nyp,nyph,out,wrk)
            uhat(:,:,iz) = wrk/nn;
            !--------------------------------
            in = v(:,:,iz); 
            call dfftw_execute_dft_r2c(planf,in,out);
            call unpadm(nxh,nxph,ny,nyp,nyph,out,wrk)
            vhat(:,:,iz) = wrk/nn;
            !--------------------------------
            in = w(:,:,iz); 
            call dfftw_execute_dft_r2c(planf,in,out);
            call unpadm(nxh,nxph,ny,nyp,nyph,out,wrk)
            what(:,:,iz) = wrk/nn;
            !--------------------------------
            in = ThetaR(:,:,iz); 
            call dfftw_execute_dft_r2c(planf,in,out);
            call unpadm(nxh,nxph,ny,nyp,nyph,out,wrk)
            ThetaRhat(:,:,iz) = wrk/nn;
            !--------------------------------
            in = qt(:,:,iz); 
            call dfftw_execute_dft_r2c(planf,in,out);
            call unpadm(nxh,nxph,ny,nyp,nyph,out,wrk)
            qthat(:,:,iz) = wrk/nn;
            !--------------------------------	
        end do
    
        uhat(1,1,:) = ubg;
        vhat(1,1,:) = vbg;
    
        do iz=1,m
            wrk = uhat(:,:,iz); 
            call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
            call dfftw_execute_dft_c2r(planr,out,in);
            u(:,:,iz) = in;
            !-------------------------
            wrk = vhat(:,:,iz); 
            call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
            call dfftw_execute_dft_c2r(planr,out,in);
            v(:,:,iz) = in;
            !-------------------------
        end do
    
    
    !----- Time iteration start----------------------
    
        Ti  = 0d0; dt0 = 2.5d0*10d0**(-3d0);
    
        nu = 0d0;
        nu(1) = 4d0*10d0**(-5d0)/1d0; !0.1d0/(((pi/dx)**2d0)**2d0 * dt0); !
        nu(2) = 4d0*10d0**(-4d0)/1d0; !0.005d0*(dz**2d0)/dt0; !
    
        CFL = 0.9d0; mindt = 1d0; max_theta=0d0; dt=dt0;
    
        !print*,'started'
        do while ( Ti <= Tfinal+5d0*dt )
    
            !---------------------------------------------------
    
            do iz=1,m
                theta_bar(iz) = sum(Theta(:,:,iz))/(nx*ny);
            end do
    
            N_theta = ( 81d0+8.1d0*maxval( abs( ( theta_bar(4:m-1)-theta_bar(2:m-3) )/(2d0*dz) ) ) )**0.5d0;
    
            dt = dt0; 
            dt = min(dt, CFL/max( maxval( sqrt(u*u/(dx**2d0)+v*v/(dy**2d0)+w*w/(dz**2d0)) ) , maxval(abs(u)/dx+abs(v)/dy+abs(w)/dz) , 1d0) ); 
            dt = min( dt, pi/(10d0*N_theta) );
            e_nu_1 = exp(dt/3d0*(-nu(1)*kk**2d0 - nu(2)*2d0/dz**2d0 ));
            e_nu_2 = e_nu_1**2d0;
            e_nu_3 = e_nu_1**3d0;
            do iz=1,m
                e_nu_1uv(:,:,iz) = exp(dt/3d0*(-nu(1)*kk(:,:,iz)**2d0 - nu(2)*2d0/dz**2d0 - tau_z(iz) )); !Damping in all modes
                e_nu_1uv(1,1,iz) = exp(dt/3d0*(- tau_z(iz) )); !1d0; !Removing the viscosity in the bar quantities and adding damping
                e_nu_2uv(:,:,iz) = e_nu_1uv(:,:,iz)**2d0;
                e_nu_3uv(:,:,iz) = e_nu_1uv(:,:,iz)**3d0;
            end do
            do iz=1,m
                e_nu_1w(:,:,iz) = exp(dt/3d0*(-nu(1)*kk(:,:,iz)**2d0 - nu(2)*2d0/dz**2d0-tau_z(iz) )); !Damping in all modes
                e_nu_1w(1,1,iz) = exp(dt/3d0*(- nu(2)*2d0/dz**2d0-tau_z(iz) ));
                e_nu_2w(:,:,iz) = e_nu_1w(:,:,iz)**2d0;
                e_nu_3w(:,:,iz) = e_nu_1w(:,:,iz)**3d0;
            end do
    
            max_theta = maxval(abs(Theta)); !max(max_theta,maxval(abs(Theta)));
            print*,"T=",Ti*Ts,' hours'
            print*, 'max_theta=', max_theta*Ths,' kelvin'
    
            !---------------------------------------------------
            !---- RK3 for momentum equations start ---------
            !  0  |  0   0   0
            ! 1/3 | 1/3  0   0
            ! 2/3 | 0   2/3  0 
            ! -----------------
            !     | 1/4  0  3/4
            ! Huen's third order method
            !****- RK1  start ********************************
            call RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vt,
            uhat,vhat,what,ThetaRhat,qthat,u,v,w,ThetaR,qt,&
            fuhat,fvhat,fwhat,fThetaRhat,fqthat,RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)
                u1hat = (uhat + dt/3d0*fuhat)*e_nu_1uv;
                v1hat = (vhat + dt/3d0*fvhat)*e_nu_1uv;		
                w1hat = (what + dt/3d0*fwhat)*e_nu_1w;		
                ThetaR1hat = (ThetaRhat + dt/3d0*fThetaRhat)*e_nu_1;
                qt1hat = (qthat + dt/3d0*fqthat)*e_nu_1;
    
                !******- Possion equation for pressure p -*******
                rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*u1hat(:,:,1:m-1)+IM*ky(:,:,1:m-1)*v1hat(:,:,1:m-1)+(w1hat(:,:,2:m)-w1hat(:,:,1:m-1))/dz;
                phat = 0d0;
                call Thomas(phat,a,b,c,rhat,nxph,nyp,m-1);
    
                !***- Update veloctiy field ----------------------
                u1hat(:,:,1:m-1) = u1hat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat;
                v1hat(:,:,1:m-1) = v1hat(:,:,1:m-1) - IM*ky(:,:,1:m-1)*phat;
                w1hat(:,:,2:m-1) = w1hat(:,:,2:m-1) - (phat(:,:,2:m-1)-phat(:,:,1:m-2))/dz;
                w1hat(:,:,1) = 0d0;
                w1hat(:,:,m) = 0d0;
    
                !***Neumann boundary condition, w1hat(1)=w1hat(1),w1hat(m)=w1hat(m)***
                !--- Get u1,v1,w1 & ThetaR1 ---
                do iz=1,m
                    wrk = u1hat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute(planr); 
                    u(:,:,iz) = in;
                    !-------------------------
                    wrk = v1hat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute(planr); 
                    v(:,:,iz) = in;
                    !-------------------------
                    wrk = w1hat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute(planr); 
                    w(:,:,iz) = in;
                    !-------------------------
                    wrk = ThetaR1hat(:,:,iz);
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out) 
                    call dfftw_execute(planr); 
                    ThetaR(:,:,iz) = in;
                    !-------------------------
                    wrk = qt1hat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute(planr); 
                    qt(:,:,iz) = in;
                end do			
    
            !************  RK1 end ***********************************
            !****- RK2  start ********************************
                  call RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vt,
                  u1hat,v1hat,w1hat,ThetaR1hat,qt1hat,u,v,w,ThetaR,qt,&
            fu1hat,fv1hat,fw1hat,fThetaR1hat,fqt1hat,RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)
                u1hat = uhat*e_nu_2uv + 2d0/3d0*dt*fu1hat*e_nu_1uv;
                v1hat = vhat*e_nu_2uv + 2d0/3d0*dt*fv1hat*e_nu_1uv;
                w1hat = what*e_nu_2w + 2d0/3d0*dt*fw1hat*e_nu_1w;
                ThetaR1hat = ThetaRhat*e_nu_2 + 2d0/3d0*dt*fThetaR1hat*e_nu_1;
                qt1hat = qthat*e_nu_2 + 2d0/3d0*dt*fqt1hat*e_nu_1;
    
                !******- Possion equation for pressure p -*******
                rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*u1hat(:,:,1:m-1)+IM*ky(:,:,1:m-1)*v1hat(:,:,1:m-1)+(w1hat(:,:,2:m)-w1hat(:,:,1:m-1))/dz;
                phat = 0d0;
                call Thomas(phat,a,b,c,rhat,nxph,nyp,m-1);			
                !*************************************************
                !***- Update veloctiy field ----------------------
                u1hat(:,:,1:m-1) = u1hat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat;
                v1hat(:,:,1:m-1) = v1hat(:,:,1:m-1) - IM*ky(:,:,1:m-1)*phat;
                w1hat(:,:,2:m-1) = w1hat(:,:,2:m-1) - (phat(:,:,2:m-1)-phat(:,:,1:m-2))/dz;
                w1hat(:,:,1) = 0d0;
                w1hat(:,:,m) = 0d0;
    
                !***Neumann boundary condition, w1hat(1)=w1hat(1),w1hat(m)=w1hat(m)***
                !--- Get u1,v1,w1 & ThetaR1 ---
                do iz=1,m
                    wrk = u1hat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in);; 
                    u(:,:,iz) = in;
                    !-------------------------
                    wrk = v1hat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in); 
                    v(:,:,iz) = in;
                    !-------------------------
                    wrk = w1hat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in);
                    w(:,:,iz) = in;
                    !-------------------------
                    wrk = ThetaR1hat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in);
                    ThetaR(:,:,iz) = in;
                    !-------------------------
                    wrk = qt1hat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in);
                    qt(:,:,iz) = in;
                end do	
                !------------------------
            !************  RK2 end ***************************
            !****- RK3  start ********************************
                 call RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vt,u1hat,v1hat,w1hat,ThetaR1hat,qt1hat,u,v,w,ThetaR,qt,&
            fu1hat,fv1hat,fw1hat,fThetaR1hat,fqt1hat,RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)
                uhat = uhat*e_nu_3uv + dt/4d0*fuhat*e_nu_3uv + 3d0/4d0*dt*fu1hat*e_nu_1uv;
                vhat = vhat*e_nu_3uv + dt/4d0*fvhat*e_nu_3uv + 3d0/4d0*dt*fv1hat*e_nu_1uv;
                what = what*e_nu_3w + dt/4d0*fwhat*e_nu_3w + 3d0/4d0*dt*fw1hat*e_nu_1w;
                ThetaRhat = ThetaRhat*e_nu_3 + dt/4d0*fThetaRhat*e_nu_3 + 3d0/4d0*dt*fThetaR1hat*e_nu_1;
                qthat = qthat*e_nu_3 + dt/4d0*fqthat*e_nu_3 + 3d0/4d0*dt*fqt1hat*e_nu_1;
            
                !******- Possion equation for pressure p -*******
                rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*uhat(:,:,1:m-1)+IM*ky(:,:,1:m-1)*vhat(:,:,1:m-1)+(what(:,:,2:m)-what(:,:,1:m-1))/dz;
                phat = 0d0;
                call Thomas(phat,a,b,c,rhat,nxph,nyp,m-1);	
        
                !***- Update veloctiy field ----------------------
                uhat(:,:,1:m-1) = uhat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat;
                vhat(:,:,1:m-1) = vhat(:,:,1:m-1) - IM*ky(:,:,1:m-1)*phat;
                what(:,:,2:m-1) = what(:,:,2:m-1) - (phat(:,:,2:m-1)-phat(:,:,1:m-2))/dz;
                what(:,:,1) = 0d0;
                what(:,:,m) = 0d0;
    
                !***Neumann boundary condition, what(1)=what(1),what(m)=what(m)***
                !--- Get u1,v1,w1 & ThetaR1 ---
                do iz=1,m
                    wrk = uhat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in); 
                    u(:,:,iz) = in;
                    !-------------------------
                    wrk = vhat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in);
                    v(:,:,iz) = in;
                    !-------------------------
                    wrk = what(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in);
                    w(:,:,iz) = in;
                    !-------------------------
                    wrk = ThetaRhat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in);
                    ThetaR(:,:,iz) = in;
                    !-------------------------
                    wrk = qthat(:,:,iz); 
                    call padm(nxh,nxph,ny,nyp,nyph,wrk,out)
                    call dfftw_execute_dft_c2r(planr,out,in);
                    qt(:,:,iz) = in;
                end do
            !************  RK3 end ***************************
            !Now we need to impose positivity in qv+qvini
            qr = 0d0;
            do iz=1,m
                where ( qvini(iz)+qt(:,:,iz) > qvs(iz) ) qr(:,:,iz) = qvini(iz)+qt(:,:,iz)-qvs(iz);
                qv(:,:,iz) = qt(:,:,iz)-qr(:,:,iz);
                where ( qv(:,:,iz)+qvini(iz) < 0d0 ) qv(:,:,iz) = -qvini(iz);
            end do
            qt = qv+qr;
            Theta = ThetaR+L*qr;
    
            !Now we need to recompute theta,qr,qv to be in agreement with the modification
            !print*,'min qv tot=',minval(qvini+qv)
            !We still need to pass ThetaR, qt to Fourier space
            do iz=1,m
                !--------------------------------
                in = ThetaR(:,:,iz); 
                call dfftw_execute_dft_r2c(planf,in,out);
                call unpadm(nxh,nxph,ny,nyp,nyph,out,wrk)
                ThetaRhat(:,:,iz) = wrk/nn;
                !--------------------------------
                in = qt(:,:,iz); 
                call dfftw_execute_dft_r2c(planf,in,out);
                call unpadm(nxh,nxph,ny,nyp,nyph,out,wrk)
                qthat(:,:,iz) = wrk/nn;
            end do
    !---------------------------------------------------
            Ti = Ti+dt;
        end do !*** end of while **************  
    
    !------------------Destroy FFTW plans------------------------------------------------
        
        call fftw_destroy_plan(planf); !fftw_estimate); !fftw_patient); !fftw_measure);
        call fftw_destroy_plan(planr); !fftw_estimate); !fftw_patient); !fftw_measure);
        call fftw_cleanup_threads();
    
    !------------------------------------------------------------------
        
        deallocate(out,phat,rhat,uhat,vhat,what,ThetaRhat,qthat,fu1hat,fv1hat,fw1hat,fThetaR1hat,fqt1hat)
        deallocate(fuhat,fvhat,fwhat,fThetaRhat,fqthat,u1hat,v1hat,w1hat,ThetaR1hat,qt1hat)
        deallocate(in,u,v,w,Theta,ThetaR,qt,qv,qr,qvini,kx,ky,kk,a,b,c)
        deallocate(e_nu_1,e_nu_2,e_nu_3,e_nu_1uv,e_nu_2uv,e_nu_3uv,e_nu_1w,e_nu_2w,e_nu_3w,dqvdz,RC,Mr,ubg,vbg,qvs,theta_bar,wrk)
        deallocate(u_bar,v_bar,wrk1D)
    
        close(10); close(11); close(12); close(13); close(14); close(15); close(16); close(17);
        close(24); close(25);
    end program FARE
    
    
    
    
program rk3 
    implicit none 
    real, dimension(100) :: w, bu, bs 
    real :: k1_w, k1_bu, k1_bs
    real :: k2_w, k2_bu, k2_bs
    real :: k3_w, k3_bu, k3_bs
    real :: tmp_w, tmp_bu, tmp_bs
    real :: dt, Ns, Nu 
    integer :: t
    Nu = 0.01
    Ns = 0.005
    dt = 0.1
    w(1) = 0.1
    bu(1) = 5
    bs(1) = 1
    do t = 1, 100-1
        call rhs(w(t), bu(t), bs(t), k1_w, k1_bu, k1_bs, dt, Ns, Nu)

        tmp_w = w(t) + dt * k1_w / 3
        tmp_bu = bu(t) + dt * k1_bu / 3
        tmp_bs = bs(t) + dt * k1_bs / 3
        call rhs(tmp_w, tmp_bu, tmp_bs, k2_w, k2_bu, k2_bs, dt, Ns, Nu)

        tmp_w = w(t) + dt * k2_w * 2/3
        tmp_bu = bu(t) + dt * k2_bu * 2/3
        tmp_bs = bs(t) + dt * k2_bs * 2/3
        call rhs(tmp_w, tmp_bu, tmp_bs, k3_w, k3_bu, k3_bs, dt, Ns, Nu)

        w(t+1) = w(t) + dt* 1/4 * k1_w + dt* 3/4 * k3_w
        bu(t+1) = bu(t) + dt* 1/4 * k1_bu + dt* 3/4 * k3_bu
        bs(t+1) = bs(t) + dt* 1/4 * k1_bs + dt* 3/4 * k3_bs
    end do
end program

subroutine rk3_ad
    DO t=99,1,-1
        CALL POPREAL4(bs(t+1))
        k3_bsb = dt*3*bsb(t+1)/4
        CALL POPREAL4(bu(t+1))
        k3_bub = dt*3*bub(t+1)/4
        k3_wb = dt*3*wb(t+1)/4
        tmp_wb = 0.0
        tmp_bub = 0.0
        tmp_bsb = 0.0
        CALL RHS_B(tmp_w, tmp_wb, tmp_bu, tmp_bub, tmp_bs, tmp_bsb, k3_w, &
        &        k3_wb, k3_bu, k3_bub, k3_bs, k3_bsb, dt, ns, nu)
        bsb(t) = bsb(t) + bsb(t+1) + tmp_bsb
        bub(t) = bub(t) + bub(t+1) + tmp_bub
        wb(t) = wb(t) + wb(t+1) + tmp_wb
        CALL POPREAL4(tmp_bs)
        k2_bsb = dt*2*tmp_bsb/3
        CALL POPREAL4(tmp_bu)
        k2_bub = dt*2*tmp_bub/3
        k2_wb = dt*2*tmp_wb/3
        tmp_wb = 0.0
        tmp_bub = 0.0
        tmp_bsb = 0.0
        CALL RHS_B(tmp_w, tmp_wb, tmp_bu, tmp_bub, tmp_bs, tmp_bsb, k2_w, &
        &        k2_wb, k2_bu, k2_bub, k2_bs, k2_bsb, dt, ns, nu)
        k1_bsb = dt*bsb(t+1)/4 + dt*tmp_bsb/3
        bsb(t+1) = 0.0
        k1_bub = dt*bub(t+1)/4 + dt*tmp_bub/3
        bub(t+1) = 0.0
        k1_wb = dt*wb(t+1)/4 + dt*tmp_wb/3
        wb(t+1) = 0.0
        CALL POPREAL4(tmp_bs)
        bsb(t) = bsb(t) + tmp_bsb
        CALL POPREAL4(tmp_bu)
        bub(t) = bub(t) + tmp_bub
        wb(t) = wb(t) + tmp_wb
        CALL RHS_B(w(t), wb(t), bu(t), bub(t), bs(t), bsb(t), k1_w, k1_wb, &
        &        k1_bu, k1_bub, k1_bs, k1_bsb, dt, ns, nu)
    END DO
        bsb(1) = 0.0
        bub(1) = 0.0
        wb(1) = 0.0
endsubroutine rk3_ad 
    
SUBROUTINE RHS_B(wr, wrb, bur, burb, bsr, bsrb, dwr, dwrb, dbur, dburb, &
    & dbsr, dbsrb, dt, ns, nu)
      IMPLICIT NONE
      REAL, INTENT(IN) :: wr, bur, bsr
      REAL :: wrb, burb, bsrb
      REAL :: dwr, dbur, dbsr
      REAL :: dwrb, dburb, dbsrb
      REAL :: dt
      INTEGER :: branch
      REAL :: ns
      REAL :: nu
      IF (bsr .LT. bur) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      wrb = wrb - ns**2*dbsrb - nu**2*dburb
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        burb = burb + dwrb
      ELSE
        bsrb = bsrb + dwrb
      END IF
END SUBROUTINE RHS_B



subroutine rhs(wr, bur, bsr, dwr, dbur, dbsr, dt, Ns, Nu)
	implicit none
	real, intent(in) :: wr, bur, bsr
	real, intent(out) :: dwr, dbur, dbsr
	real :: dt 

	IF (bsr .LT. bur) THEN
		dwr = bur
	ELSE
		dwr = bsr
	END IF
	dbur = -(Nu*Nu*wr)
	dbsr = -(Ns*Ns*wr)
endsubroutine rhs