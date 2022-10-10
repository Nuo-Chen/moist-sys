subroutine FARE
    ! use, intrinsic :: iso_c_binding
    ! use FFTW3
    ! use mainsubroutines
    ! use mainvars
    ! use mpads
    implicit none
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
    integer :: it
    integer :: nx, ny, m, nxh, nxp, nxph, nyp, nyph

    real :: Us, Ts, Ths
    real :: Lx, Ly, Lz, dx, dy, dz
    real :: L, vt, tau, qvs0, qv0
    real :: a_squall, f
    integer :: Ls, qs
    real :: Tfinal, IM, eps, epsbar, Nz, f_star, g_star, B_star
    real :: drx, dry, drz

    !Real
    real :: max_theta, max_pert
    real :: xi,yj,da,x_c,y_c,z_c,r_c, ampl_bubble
    real ,dimension(:), allocatable :: theta_bar, u_bar, v_bar, wrk1D

!----------------PARAMETERS----------------------------------------------
    !Grid size:
    nx  = 128; ny  = 128; m  = 100+1;
    nxh = nx/2+1; nxp = nx/3; nxp = 2*nxp+1; nxph = (nxp-1)/2+1;
    nyp = ny/3; nyp = 2*nyp+1; nyph = (nyp-1)/2+1;

    !Quantity scales:
    Us = 100/9; !In m/s
    Ts = 15/60; !In hours
    Ths = 3.0; !In Kelvin
    Ls = 10; !In kms
    qs = 10; !In g/kg

    !*******DOMAIN SIZE*******
    Lx = 128; Ly = 128; Lz = 15;
    Lx = Lx/Ls; Ly = Ly/Ls; Lz = Lz/Ls; !Non-dimentionalization
    dx = Lx/nx; dy = Ly/ny; dz = Lz/(m-1);

    !Rain fall velocity
    vt = 5;
    vt = vt/Us; !Non-dimentionalization

    !Latent heat
    L = 2.5; !In 10^{6}*J*kg^{-1}
    L = L/3*10; !Non-dimentionalization prefector epsilon^{-1} L = epsilon^2 * L^d * \theta_0 ([\Theta]*c_p*T_0)^{-1}
    !Pendiente: Cambiar a LcpTheta0

    !Relaxation time
    tau = 4; ! In hours
    tau = tau/Ts; !Non-dimensionalization
    tau = 0; !Turned off

    !Qvs at surface
    qvs0 = 20; !In g/kg
    qvs0 = qvs0/qs; !Non-dimensionalization 
    
    qv0 = 18; !In g/kg
    qv0 = qv0/qs;

    !Max perturbation in theta
    max_pert = 0.1*0; !In Kelvin
    max_pert = max_pert/Ths; !Non-dimensionalization

    !Parameter for the zonal vel background
    a_squall = 10; !In m/s
    a_squall = a_squall/Us;
    a_squall = a_squall*0; !Turned off

    !Coriolis parameter
    f = sin(25* pi /180);

    !Final time
    Tfinal = 5; !In hours
    Tfinal = Tfinal/Ts; !Non-dimensionalization

    !***************
    IM = (0,1)

    !***************
    eps = 0.1;
    epsbar = 0.622;
    Nz = 1;

    !***************
    
    f_star = eps*f;
    g_star = 10; !7.9461;
    B_star = 10;
    
    !Wave numbers are multiple of those
    drx = 2*pi/Lx; dry = 2*pi/Ly; drz = 2*pi/Lz;

    !!!!Data for Bubble:
    x_c = Lx/2; y_c = Ly/2; z_c = 2/Ls; 
    !rx_c = 10; ry_c = 10; rz_c = 1;
    !rx_c = rx_c/Ls; ry_c = ry_c/Ls; rz_c = rz_c/Ls; 
    ampl_bubble = 8/4/qs;
    r_c = 5;
    r_c = r_c/Ls;

!-------Allocation----------------------

	!Real
	allocate(u(nx,ny,m)); allocate(v(nx,ny,m)); allocate(w(nx,ny,m)); allocate(Theta(nx,ny,m));	
	allocate(ThetaR(nx,ny,m)); allocate(qv(nx,ny,m)); allocate(qr(nx,ny,m)); allocate(qt(nx,ny,m));
	allocate(qvini(m)); allocate(RC(m)); allocate(Mr(m)); allocate(ubg(m)); allocate(vbg(m)); 
	allocate(dqvdz(m)); allocate(qvs(m));
	allocate(theta_bar(m)); allocate(u_bar(m)); allocate(v_bar(m));
	allocate(wrk1D(m));
	
	!Complex
	allocate(in(nx,ny)); allocate(out(nxh,ny));
	allocate(wrk(nxph,nyp)); allocate(tau_z(m));
	allocate(uhat(nxph,nyp,m)); allocate(vhat(nxph,nyp,m)); allocate(what(nxph,nyp,m)); 
	allocate(ThetaRhat(nxph,nyp,m)); allocate(qthat(nxph,nyp,m));
	allocate(u1hat(nxph,nyp,m)); allocate(v1hat(nxph,nyp,m)); allocate(w1hat(nxph,nyp,m)); 
	allocate(ThetaR1hat(nxph,nyp,m)); allocate(qt1hat(nxph,nyp,m));
	allocate(fuhat(nxph,nyp,m)); allocate(fvhat(nxph,nyp,m)); allocate(fwhat(nxph,nyp,m)); 
	allocate(fThetaRhat(nxph,nyp,m)); allocate(fqthat(nxph,nyp,m));
	allocate(fu1hat(nxph,nyp,m)); allocate(fv1hat(nxph,nyp,m)); allocate(fw1hat(nxph,nyp,m)); 
	allocate(fThetaR1hat(nxph,nyp,m)); allocate(fqt1hat(nxph,nyp,m));
	allocate(kx(nxph,nyp,m)); allocate(ky(nxph,nyp,m)); allocate(kk(nxph,nyp,m));
	allocate(e_nu_1(nxph,nyp,m)); allocate(e_nu_2(nxph,nyp,m)); allocate(e_nu_3(nxph,nyp,m));
	allocate(e_nu_1uv(nxph,nyp,m)); allocate(e_nu_2uv(nxph,nyp,m)); allocate(e_nu_3uv(nxph,nyp,m));
	allocate(e_nu_1w(nxph,nyp,m)); allocate(e_nu_2w(nxph,nyp,m)); allocate(e_nu_3w(nxph,nyp,m));
	allocate(a(nxph,nyp,m-1));  allocate(b(nxph,nyp,m-1));  allocate(c(nxph,nyp,m-1)); 
	allocate(rhat(nxph,nyp,m-1)); allocate(phat(nxph,nyp,m-1));

    real, dimension()
!-------------------FFTW PACKAGE-----------------------------------------

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

    do iz=1,m
        qvs0(m) = 18.04 + 3.27 + 0.1  + 0.1804+3.48
        qvini(m) = 18.04 + 3.27 + 0.1  + 0.1804+3.48
        dqvdz = 18.04 + 3.27 + 0.1 
    end do
    ! call FQv(qvs,m,dz,qvs0);
    ! call FQv(qvini,m,dz,qv0);
    ! call Fdqvdz(dqvdz,m,dz,qv0);
    
    Theta = 0;

    !Define the radiactive cooling rate
    RC = 0; !Turned off
    !do iz=1,m
    !	zk = (iz-0.5)*dz;
    !	RC(iz) = -exp(-0.2*zk)*zk*(1.5-zk)**2;
    !	RC(iz) = RC(iz)*(1/0.45)*(1/288)*25; !Profile similar to Majda-Xing, max about 25 K/day
    !end do 
    
    !Define the moistening proportional to dqvdz
    Mr = 0; !Turned off
    !do iz=1,m
    !	Mr(iz) = 0.1*RC(iz)*dqvdz(iz)*0; !Turned off
    !end do

    !!!!!!!!!!!!!!!!!!----Damping in w: ----
    tau_z = 0;
    do iz=1,m
        zk=(iz-1)*dz;
        tau_z(iz) = (1/(1/Ts))*exp( 50*( (1/Lz)**2-(1/(zk-Lz))**2 ) );
    end do

    tau_z = 0; !Turned off
    
    ubg = 0; vbg = 0; !Turned off
     
    !!!!!!!!!!!!!!!Pressure solver !!!!
    !Three diagonals for Thomas' method
    !Main diagonal (entries from 1 to m-1), to compute variables in the staggered grid, level at p
    b = - ( kk + 1/dz**2);
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

!--------------INITIAL DATA-----------------------------------------------

    u = 0; v = 0; w = 0; Theta = 0; qv = 0; qr = 0;
    
    !---------Perturbation in temperature-------------------
    Theta = 0; !Turned off

    
    !-------Moist bubble
    qt = 0;
    do iy = 1,ny
        yj = (iy-1)*dy;
        do ix = 1,nx
            xi = (ix-1)*dx;
            da = sqrt((xi-x_c)**2+(yj-y_c)**2);
            do iz=1,m-1			
                zk = (iz-0.5)*dz; !For level at u
                if ( da <= r_c .and. zk >= z_c-0.1 .and. zk <= z_c+0.1 ) then		
                    qt(ix,iy,iz) = ampl_bubble*cos(pi*da/(2*r_c))*(10*(zk-z_c-0.1))**2*(10*(zk-z_c+0.1))**2; !*(10*(zk-z_c-0.3))**2;
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
    
    qr = 0;
    do iz=1,m
        where ( qvini(iz)+qt(:,:,iz) > qvs(iz) ) qr(:,:,iz) = qvini(iz)+qt(:,:,iz)-qvs(iz);
        qv(:,:,iz) = qt(:,:,iz)-qr(:,:,iz);
        where ( qv(:,:,iz)+qvini(iz) < 0 ) qv(:,:,iz) = -qvini(iz); !Positivity of qvini+qv
    end do
    qt = qv+qr;
    ThetaR = Theta-L*qr;

!----Pass to Fourier space---------------------------------------------------------

    nn = nx*ny
    do iz=1,m
        in = u(:,:,iz); 
         
        wrk = in
        uhat(:,:,iz) = wrk/nn;
        !--------------------------------
        in = v(:,:,iz); 
        
        wrk = in
        vhat(:,:,iz) = wrk/nn;
        !--------------------------------
        in = w(:,:,iz); 
        
        wrk = in
        what(:,:,iz) = wrk/nn;
        !--------------------------------
        in = ThetaR(:,:,iz); 
        
        wrk = in
        ThetaRhat(:,:,iz) = wrk/nn;
        !--------------------------------
        in = qt(:,:,iz); 
        
        wrk = in
        qthat(:,:,iz) = wrk/nn;
        !--------------------------------	
    end do

    uhat(1,1,:) = ubg;
    vhat(1,1,:) = vbg;

    do iz=1,m
        wrk = uhat(:,:,iz); 
        wrk = in
        
        u(:,:,iz) = in;
        !-------------------------
        wrk = vhat(:,:,iz); 
        wrk = in
        
        v(:,:,iz) = in;
        !-------------------------
    end do


!----- Time iteration start----------------------

    Ti  = 0; dt0 = 2.5*10**(-3);

    nu = 0;
    nu(1) = 4*10**(-5)/1; !0.1/(((pi/dx)**2)**2 * dt0); !
    nu(2) = 4*10**(-4)/1; !0.005*(dz**2)/dt0; !

    CFL = 0.9; mindt = 1; max_theta=0; dt=dt0;

    !print*,'started'
    do while ( Ti <= Tfinal+5*dt )
        !---------------------------------------------------

        do iz=1,m
            theta_bar(iz) = sum(Theta(:,:,iz))/(nx*ny);
        end do

        N_theta = ( 81+8.1*maxval( abs( ( theta_bar(4:m-1)-theta_bar(2:m-3) )/(2*dz) ) ) )**0.5;

        dt = dt0; 
        dt = min(dt, CFL/max( maxval( sqrt(u*u/(dx**2)+v*v/(dy**2)+w*w/(dz**2)) ) , maxval(abs(u)/dx+abs(v)/dy+abs(w)/dz) , 1) ); 
        dt = min( dt, pi/(10*N_theta) );
        e_nu_1 = exp(dt/3*(-nu(1)*kk**2 - nu(2)*2/dz**2 ));
        e_nu_2 = e_nu_1**2;
        e_nu_3 = e_nu_1**3;
        do iz=1,m
            e_nu_1uv(:,:,iz) = exp(dt/3*(-nu(1)*kk(:,:,iz)**2 - nu(2)*2/dz**2 - tau_z(iz) )); !Damping in all modes
            e_nu_1uv(1,1,iz) = exp(dt/3*(- tau_z(iz) )); !1; !Removing the viscosity in the bar quantities and adding damping
            e_nu_2uv(:,:,iz) = e_nu_1uv(:,:,iz)**2;
            e_nu_3uv(:,:,iz) = e_nu_1uv(:,:,iz)**3;
        end do
        do iz=1,m
            e_nu_1w(:,:,iz) = exp(dt/3*(-nu(1)*kk(:,:,iz)**2 - nu(2)*2/dz**2-tau_z(iz) )); !Damping in all modes
            e_nu_1w(1,1,iz) = exp(dt/3*(- nu(2)*2/dz**2-tau_z(iz) ));
            e_nu_2w(:,:,iz) = e_nu_1w(:,:,iz)**2;
            e_nu_3w(:,:,iz) = e_nu_1w(:,:,iz)**3;
        end do

        max_theta = maxval(abs(Theta)); !max(max_theta,maxval(abs(Theta)));
        print*,"T=",Ti*Ts,' hours'
        print*, 'max_theta=', max_theta*Ths,' kelvin'

        !---------------------------------------------------
        !---- RK3 for momentum equations start ---------
        !****- RK1  start ********************************
        call RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vt,uhat,vhat,what,ThetaRhat,qthat,u,v,w,ThetaR,qt,&
        fuhat,fvhat,fwhat,fThetaRhat,fqthat,RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)
        
        u1hat = (uhat + dt/3*fuhat)*e_nu_1uv;
        v1hat = (vhat + dt/3*fvhat)*e_nu_1uv;		
        w1hat = (what + dt/3*fwhat)*e_nu_1w;		
        ThetaR1hat = (ThetaRhat + dt/3*fThetaRhat)*e_nu_1;
        qt1hat = (qthat + dt/3*fqthat)*e_nu_1;

        !******- Possion equation for pressure p -*******
        rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*u1hat(:,:,1:m-1)+IM*ky(:,:,1:m-1)*v1hat(:,:,1:m-1)+(w1hat(:,:,2:m)-w1hat(:,:,1:m-1))/dz;
        phat = 0;

        call Thomas(phat,a,b,c,rhat,nxph,nyp,m-1);

        !***- Update veloctiy field ----------------------
        u1hat(:,:,1:m-1) = u1hat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat;
        v1hat(:,:,1:m-1) = v1hat(:,:,1:m-1) - IM*ky(:,:,1:m-1)*phat;
        w1hat(:,:,2:m-1) = w1hat(:,:,2:m-1) - (phat(:,:,2:m-1)-phat(:,:,1:m-2))/dz;
        w1hat(:,:,1) = 0;
        w1hat(:,:,m) = 0;

        !***Neumann boundary condition, w1hat(1)=w1hat(1),w1hat(m)=w1hat(m)***
        !--- Get u1,v1,w1 & ThetaR1 ---
        do iz=1,m
            wrk = u1hat(:,:,iz); 
            wrk = in
             
            u(:,:,iz) = in;
            !-------------------------
            wrk = v1hat(:,:,iz); 
            wrk = in
             
            v(:,:,iz) = in;
            !-------------------------
            wrk = w1hat(:,:,iz); 
            wrk = in
             
            w(:,:,iz) = in;
            !-------------------------
            wrk = ThetaR1hat(:,:,iz);
            wrk = in 
             
            ThetaR(:,:,iz) = in;
            !-------------------------
            wrk = qt1hat(:,:,iz); 
            wrk = in
             
            qt(:,:,iz) = in;
        end do			

        !************  RK1 end ***********************************
        !****- RK2  start ********************************
        call RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vt,u1hat,v1hat,w1hat,ThetaR1hat,qt1hat,u,v,w,ThetaR,qt,&
        fu1hat,fv1hat,fw1hat,fThetaR1hat,fqt1hat,RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)
        
        u1hat = uhat*e_nu_2uv + 2/3*dt*fu1hat*e_nu_1uv;
        v1hat = vhat*e_nu_2uv + 2/3*dt*fv1hat*e_nu_1uv;
        w1hat = what*e_nu_2w + 2/3*dt*fw1hat*e_nu_1w;
        ThetaR1hat = ThetaRhat*e_nu_2 + 2/3*dt*fThetaR1hat*e_nu_1;
        qt1hat = qthat*e_nu_2 + 2/3*dt*fqt1hat*e_nu_1;

        !******- Possion equation for pressure p -*******
        rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*u1hat(:,:,1:m-1)+IM*ky(:,:,1:m-1)*v1hat(:,:,1:m-1)+(w1hat(:,:,2:m)-w1hat(:,:,1:m-1))/dz;
        phat = 0;

        call Thomas(phat,a,b,c,rhat,nxph,nyp,m-1);

        !*************************************************
        !***- Update veloctiy field ----------------------
        u1hat(:,:,1:m-1) = u1hat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat;
        v1hat(:,:,1:m-1) = v1hat(:,:,1:m-1) - IM*ky(:,:,1:m-1)*phat;
        w1hat(:,:,2:m-1) = w1hat(:,:,2:m-1) - (phat(:,:,2:m-1)-phat(:,:,1:m-2))/dz;
        w1hat(:,:,1) = 0;
        w1hat(:,:,m) = 0;

        !***Neumann boundary condition, w1hat(1)=w1hat(1),w1hat(m)=w1hat(m)***
        !--- Get u1,v1,w1 & ThetaR1 ---
        do iz=1,m
            wrk = u1hat(:,:,iz); 
            wrk = in
            
            u(:,:,iz) = in;
            !-------------------------
            wrk = v1hat(:,:,iz); 
            wrk = in
             
            v(:,:,iz) = in;
            !-------------------------
            wrk = w1hat(:,:,iz); 
            wrk = in
            
            w(:,:,iz) = in;
            !-------------------------
            wrk = ThetaR1hat(:,:,iz); 
            wrk = in
            
            ThetaR(:,:,iz) = in;
            !-------------------------
            wrk = qt1hat(:,:,iz); 
            wrk = in
            
            qt(:,:,iz) = in;
        end do	
            !------------------------
        !************  RK2 end ***************************
        !****- RK3  start ********************************
        call RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vt,u1hat,v1hat,w1hat,ThetaR1hat,qt1hat,u,v,w,ThetaR,qt,&
        fu1hat,fv1hat,fw1hat,fThetaR1hat,fqt1hat,RC,Mr,ubg,vbg,tau,qvs,qvini,dqvdz,in,out,planf)
        
        uhat = uhat*e_nu_3uv + dt/4*fuhat*e_nu_3uv + 3/4*dt*fu1hat*e_nu_1uv;
        vhat = vhat*e_nu_3uv + dt/4*fvhat*e_nu_3uv + 3/4*dt*fv1hat*e_nu_1uv;
        what = what*e_nu_3w + dt/4*fwhat*e_nu_3w + 3/4*dt*fw1hat*e_nu_1w;
        ThetaRhat = ThetaRhat*e_nu_3 + dt/4*fThetaRhat*e_nu_3 + 3/4*dt*fThetaR1hat*e_nu_1;
        qthat = qthat*e_nu_3 + dt/4*fqthat*e_nu_3 + 3/4*dt*fqt1hat*e_nu_1;
    
        !******- Possion equation for pressure p -*******
        rhat(:,:,1:m-1) = IM*kx(:,:,1:m-1)*uhat(:,:,1:m-1)+IM*ky(:,:,1:m-1)*vhat(:,:,1:m-1)+(what(:,:,2:m)-what(:,:,1:m-1))/dz;
        phat = 0;

        call Thomas(phat,a,b,c,rhat,nxph,nyp,m-1);	

        !***- Update veloctiy field ----------------------
        uhat(:,:,1:m-1) = uhat(:,:,1:m-1) - IM*kx(:,:,1:m-1)*phat;
        vhat(:,:,1:m-1) = vhat(:,:,1:m-1) - IM*ky(:,:,1:m-1)*phat;
        what(:,:,2:m-1) = what(:,:,2:m-1) - (phat(:,:,2:m-1)-phat(:,:,1:m-2))/dz;
        what(:,:,1) = 0;
        what(:,:,m) = 0;

        !***Neumann boundary condition, what(1)=what(1),what(m)=what(m)***
        !--- Get u1,v1,w1 & ThetaR1 ---
        do iz=1,m
            wrk = uhat(:,:,iz); 
            wrk = in
             
            u(:,:,iz) = in;
            !-------------------------
            wrk = vhat(:,:,iz); 
            wrk = in
            
            v(:,:,iz) = in;
            !-------------------------
            wrk = what(:,:,iz); 
            wrk = in
            
            w(:,:,iz) = in;
            !-------------------------
            wrk = ThetaRhat(:,:,iz); 
            wrk = in
            
            ThetaR(:,:,iz) = in;
            !-------------------------
            wrk = qthat(:,:,iz); 
            wrk = in
            
            qt(:,:,iz) = in;
        end do
        !************  RK3 end ***************************
        !Now we need to impose positivity in qv+qvini
        qr = 0;
        do iz=1,m
            where ( qvini(iz)+qt(:,:,iz) > qvs(iz) ) qr(:,:,iz) = qvini(iz)+qt(:,:,iz)-qvs(iz);
            qv(:,:,iz) = qt(:,:,iz)-qr(:,:,iz);
            where ( qv(:,:,iz)+qvini(iz) < 0 ) qv(:,:,iz) = -qvini(iz);
        end do
        qt = qv+qr;
        Theta = ThetaR+L*qr;

        !Now we need to recompute theta,qr,qv to be in agreement with the modification
        !print*,'min qv tot=',minval(qvini+qv)
        !We still need to pass ThetaR, qt to Fourier space
        do iz=1,m
            !--------------------------------
            in = ThetaR(:,:,iz); 
            
            wrk = in
            ThetaRhat(:,:,iz) = wrk/nn;
            !--------------------------------
            in = qt(:,:,iz); 
            
            wrk = in
            qthat(:,:,iz) = wrk/nn;
        end do
!---------------------------------------------------
        Ti = Ti+dt;
    end do !*** end of while **************  

end subroutine FARE





