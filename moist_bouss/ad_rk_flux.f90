subroutine RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vrain, &
	uhat,vhat,what,ThetaRhat,qthat,u,v,w,ThetaR,qt, &
	fu,fv,fw,fThetaR,fqt,RC,Mr,ubg,vbg,tau,qvs,qvini,fqvdz,in1,out1,planf1)
	! use, intrinsic :: iso_c_binding
	! use FFTW3
	! use mpads
	! use varsRK
	implicit none
    !Parameters:
    !Integers
    ! type(C_PTR),intent(in) :: planf1
    integer :: ix,iy,iz,ns
    !Complex
    complex :: IM
    !Real:
    real :: qvbar
    !In:
    !Integers
    integer, intent(in) :: nx,nxh,nxph,ny,nyp,nyph,m
    !Real:
    real, intent(in) :: dx,dz,vrain, tau, f_star, g_star, epsbar, L, B_star
    real, dimension(2),intent(in) :: nu   !*** diffusion parameters 
    real, dimension(nxph,nyp,m), intent(in) :: kx,ky
    real, dimension(m), intent(in) :: qvs, fqvdz
    real, dimension(m), intent(in) :: qvini,RC,Mr,ubg,vbg
    real,dimension(nx,ny,m),intent(in) :: u,v,w,ThetaR, qt,    qr
    real ,dimension(nx,ny),intent(inout) ::  in1
    !Complex
    complex, dimension(nxph,nyp,m), intent(in) :: uhat,vhat,what, qthat, ThetaRhat,    wrkNLhat
    complex, dimension(nxh,ny), intent(inout) :: out1
    complex, dimension(nx,ny, m) :: wrkNL
    complex, dimension(nxph,nyp) :: wrk1

    !Out:
    !Complex:
    complex, dimension(nxph,nyp,m), intent(out) :: fu,fv,fw,fThetaR,fqt

    !---------------------------------------------
    IM = (0,1);
    ns = nx*ny;
    
    !---------------------------------------------
    !Allocation:
    !Real
    ! allocate(qr(nx,ny,m));
    !Complex
    ! allocate(wrkNL(nx,ny,m)); allocate(wrkNLhat(nxph,nyp,m));
    ! allocate(wrk1(nxph,nyp));
    
    !------------------------Equations to solve in non-dimensional units:------------------------
    !d u /dt = fv - (uu)_x- (vu)_y -(wu)_z + nu \nabla^2 u
    !d v /dt = -fu - (uv)_x -(vv)_y -(wv)_z + nu \nabla^2 v
    !d w /dt = -(uw)_x - (vw)_y - (ww)_z +nu \nabla^2 w + g_star theta' + g_star eps_o q_v' - g_star*q_r
    !d theta_r' /dt = -(u theta_r)_x - (v \theta_r)_y -B_star w - V_T L d q_r/dz
    !d qt' /dt = -(uq_t')_x-(vq_t')_y - (wq_t')_z -w d_qvditle(z)/dz -V_T dq_r/dz 
    
    !-----------Boundary conditions-------------------------------------
    !u=0 at the bottom, du/dz=0 on the top
    !**** v=0 at the bottom, dv/dz=0 on the top ****
    !Zero Dirichlet for w
    !Neumann bc for theta_r'
    !Neumann boundary condition for q_r		
    
    !----------- Notes for us-------------------------------------
    !q_t here is q_t', the fluctuation
    !Theta_r is theta_r', the fluctuation
    !Remove viscosity to VSHF (Vertically Shear Horizontal Flows), except for theta_r', q_t'

    !---------Coriolis---------------------------------------
    
    fu = f_star*vhat;
    fv = -f_star*uhat;
    
    !---------Non-linear advection term: - (uu)_x---------------------------------------
    !wrkNL = u*u;
    do iz = 1,m
        in1 = u(:,:,iz)*u(:,:,iz); !wrkNL(:,:,iz); !This is why we call it pseudo-spectral
        
        wrk1 = in1; !de-aliasing
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do
    fu = fu - IM*kx*wrkNLhat;
    
    !---------Non-linear advection terms: - (uv)_x and - (uv)_y---------------------------------------	
    !wrkNL = u*v;
    do iz = 1,m
        in1 = u(:,:,iz)*v(:,:,iz); !wrkNL(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;		
    end do
    fu = fu - IM*ky*wrkNLhat;
    fv = fv - IM*kx*wrkNLhat;
    
    !---------Non-linear advection term: - (vv)_x----------------------------------------	
    !wrkNL = v*v;
    do iz = 1,m
        in1 = v(:,:,iz)*v(:,:,iz); !wrkNL(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;		
    end do	
    fv = fv - IM*ky*wrkNLhat;
    
    !----------Non-linear advection term: - (wu)_z affecting the u equation --------------------------------------	

    !Need to compute u at the w level, due to staggered grid
    !Interpolate u:
    wrkNL(:,:,2:m-1) = 0.5*(u(:,:,1:m-2)+u(:,:,2:m-1)); !*** u on position as w
    !**** No boundary condition, but w=0 ****
    wrkNL(:,:,1) = 0;
    wrkNL(:,:,m) = 0;
    wrkNL = wrkNL*w;
    do iz = 1,m
        in1 = wrkNL(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do	
    
    fu(:,:,2:m-2) = fu(:,:,2:m-2) - (wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,2:m-2))/dz &
                    + nu(2)*(uhat(:,:,3:m-1)+uhat(:,:,1:m-3))/dz**2;
                    
    !**** u=0 at the bottom, du/dz=0 on the top ****
    fu(:,:,1)     = fu(:,:,1) - (wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/dz &
                + nu(2)*(uhat(:,:,2)-uhat(:,:,1))/dz**2; !u(0) = -u(1)
    fu(:,:,m-1)   = fu(:,:,m-1) - (wrkNLhat(:,:,m)-wrkNLhat(:,:,m-1))/dz &
                + nu(2)*(uhat(:,:,m-1)+uhat(:,:,m-2))/dz**2;	!u(m)  = u(m-1)

    fu(1,1,2:m-2) = fu(1,1,2:m-2)-nu(2)*(uhat(1,1,3:m-1)+uhat(1,1,1:m-3))/dz**2; !Remove viscosity to VSHF
    fu(1,1,1) = fu(1,1,1)-nu(2)*(uhat(1,1,2)-uhat(1,1,1))/dz**2; !Remove viscosity to VSHF
    fu(1,1,m-1) = fu(1,1,m-1)-nu(2)*(uhat(1,1,m-1)+uhat(1,1,m-2))/dz**2; !Remove viscosity to VSHF
    
    !**************************************	
    !!Ralaxation terms
    !!(u,v) horizontal velocity component
    !do iz=1,m
    !	fu(1,1,iz) = fu(1,1,iz)-(1/tau)*(uhat(1,1,iz)-ubg(iz));
    !	fv(1,1,iz) = fv(1,1,iz)-(1/tau)*(vhat(1,1,iz)-vbg(iz));
    !end do
    
    !----------Non-linear advection term: - (wu)_x affecting the w equation --------------------------------------	
    
    fw = -IM*kx*wrkNLhat;
    
    !----------Non-linear advection term: - (wv)_z affecting the v equation --------------------------------------	
    !Need to compute v at the w level, due to staggered grid
    !********* Interpolation for v ********
    wrkNL(:,:,2:m-1) = 0.5*(v(:,:,1:m-2)+v(:,:,2:m-1)); !*** v on position as w
    !**** No boundary condition, but w=0 ****
    wrkNL(:,:,1) = 0;
    wrkNL(:,:,m) = 0;
    wrkNL = wrkNL*w;
    !----------------------------------
    do iz = 1,m
        in1 = wrkNL(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do	
    !-----------------------------------
    fv(:,:,2:m-2) = fv(:,:,2:m-2) - (wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,2:m-2))/dz &
                    + nu(2)*(vhat(:,:,3:m-1)+vhat(:,:,1:m-3))/dz**2;
                    
    !**** v=0 at the bottom, dv/dz=0 on the top ****
    fv(:,:,1)     = fv(:,:,1) - (wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/dz &
                    + nu(2)*(vhat(:,:,2)-vhat(:,:,1))/dz**2;
    fv(:,:,m-1)   = fv(:,:,m-1) - (wrkNLhat(:,:,m)-wrkNLhat(:,:,m-1))/dz &
                    + nu(2)*(vhat(:,:,m-1)+vhat(:,:,m-2))/dz**2;	

    fv(1,1,2:m-2) = fv(1,1,2:m-2)-nu(2)*(vhat(1,1,3:m-1)+vhat(1,1,1:m-3))/dz**2; !Remove viscosity to VSHF
    fv(1,1,1) = fv(1,1,1)-nu(2)*(vhat(1,1,2)-vhat(1,1,1))/dz**2; !Remove viscosity to VSHF
    fv(1,1,m-1) = fv(1,1,m-1)-nu(2)*(vhat(1,1,m-1)+vhat(1,1,m-2))/dz**2; !Remove viscosity to VSHF
    
    !----------Non-linear advection term: - (wv)_y affecting the w equation --------------------------------------	
    fw = fw - IM*ky*wrkNLhat;
    
    !---------------------------------------------------------------------------------
    !Here we need to get qr
    qr = 0;
    do iz=1,m
        where ( qvini(iz)+qt(:,:,iz) > qvs(iz) ) qr(:,:,iz) = qvini(iz)+qt(:,:,iz)-qvs(iz);
    end do
    
    !------------q_v '  term affecting the w equation---------------------------------------------------------------------

    do iz = 1,m
        qvbar = sum(qt(:,:,iz)-qr(:,:,iz))/ns; !Horizontal average of qv'
        in1 = qt(:,:,iz)-qr(:,:,iz)-qvbar;
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do
    fw(:,:,2:m-1) = fw(:,:,2:m-1) + g_star*epsbar*0.5*(wrkNLhat(:,:,1:m-2)+wrkNLhat(:,:,2:m-1));
    
    !------------q_r '  term affecting the w equation---------------------------------------------------------------------
    do iz = 1,m
        in1 = qr(:,:,iz);
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do
    fw(:,:,2:m-1) = fw(:,:,2:m-1) - g_star*0.5*(wrkNLhat(:,:,1:m-2)+wrkNLhat(:,:,2:m-1));
    
    
    !------------q_r  term affecting the theta_r equation---------------------------------------------------------------------

    fThetaR = 0;
    fThetaR(:,:,2:m-2) = fThetaR(:,:,2:m-2)-vrain*L*(wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,1:m-3))/(2*dz);
    !**** Neumann boundary condition ****
    fThetaR(:,:,1) = fThetaR(:,:,1)-vrain*L*(wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/(2*dz);
    fThetaR(:,:,m-1) = fThetaR(:,:,m-1)-vrain*L*(wrkNLhat(:,:,m-1)-wrkNLhat(:,:,m-2))/(2*dz);

    !------------q_r  term affecting the q_t equation---------------------------------------------------------------------

    fqt = 0;
    fqt(:,:,2:m-2) = fqt(:,:,2:m-2) + vrain*(wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,1:m-3))/(2*dz);
    !**** Neumann boundary condition for q_r		
    fqt(:,:,1) = fqt(:,:,1) + vrain*(wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/(2*dz); !q_r(0) = q_r(1)
    fqt(:,:,m-1) = fqt(:,:,m-1) + vrain*(wrkNLhat(:,:,m-1)-wrkNLhat(:,:,m-2))/(2*dz); !q_r(m) = q_r(m-1)
    
    !-------------Theta term affecting the w equation--------------------------------------------------------------------
    do iz = 1,m
        in1 = ThetaR(:,:,iz)+L*qr(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do
    fw(:,:,2:m-1) = fw(:,:,2:m-1) + g_star*( 0.5*(wrkNLhat(:,:,1:m-2)+wrkNLhat(:,:,2:m-1)) );
    
    !----------Non-linear advection term: - (ww)_z affecting the w equation --------------------------------------	
    !********** Interpolation for w ***************
    wrkNL(:,:,1:m-1) = 0.5*(w(:,:,1:m-1)+w(:,:,2:m));  !*** w on position as u
    !**** Zero boundary condition, 0.5*(w(m+1)+w(m-1))=0, w(m)=0.0
    wrkNL(:,:,m) = -wrkNL(:,:,m-1);
    wrkNL = wrkNL*wrkNL;

    do iz = 1,m
        in1 = wrkNL(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do	
    
    fw(:,:,2:m-1) = fw(:,:,2:m-1) - (wrkNLhat(:,:,2:m-1)-wrkNLhat(:,:,1:m-2))/dz &
    + nu(2)*(what(:,:,3:m)+what(:,:,1:m-2))/dz**2 ;

    !fw(1,1,2:m-1) = fw(1,1,2:m-1)-nu(2)*(what(1,1,3:m)+what(1,1,1:m-2))/dz**2; !Remove viscosity to VSHF

    !**** Zero boundary condition ******	
    fw(:,:,1) = 0;
    fw(:,:,m) = 0;
    
    !----------Non-linear advection term: - (u theta_r')_x affecting the theta_r equation --------------------------------------	
    !wrkNL = u*ThetaR;
    do iz = 1,m
        in1 = u(:,:,iz)*ThetaR(:,:,iz); !wrkNL(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do
    fThetaR = fThetaR-IM*kx*wrkNLhat;
    
    !----------Non-linear advection term: - (v theta_r')_y affecting the theta_r equation --------------------------------------	
    !wrkNL = v*ThetaR;
    do iz = 1,m
        in1 = v(:,:,iz)*ThetaR(:,:,iz); !wrkNL(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do
    fThetaR = fThetaR - IM*ky*wrkNLhat;
    
    !---------------------------------------------------------------------------------
    !********************* For Fixed Background
    do iz = 1,m-2
        fThetaR(:,:,iz) = fThetaR(:,:,iz)-0.5*(what(:,:,iz)+what(:,:,iz+1))*B_star; 
    end do
    
    !----------Non-linear advection term: - (w theta_r')_z affecting the theta_r equation --------------------------------------	
    !********** Interpolation for ThetaR **************
    wrkNL(:,:,2:m-1) = 0.5*(ThetaR(:,:,1:m-2)+ThetaR(:,:,2:m-1)); !*** ThetaR on position as w
    !****   Neumann boundary condition ****
    wrkNL(:,:,1) = ThetaR(:,:,1);
    wrkNL(:,:,m) = ThetaR(:,:,m-1);
    wrkNL = w*wrkNL;
    !----------------------------------
    do iz = 1,m
        in1 = wrkNL(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do	
    !**************************************	
    
    fThetaR(:,:,2:m-2) = fThetaR(:,:,2:m-2) - (wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,2:m-2))/dz &
    + nu(2)*(ThetaRhat(:,:,3:m-1)+ThetaRhat(:,:,1:m-3))/dz**2;
    !**** Neumann boundary condition ****
    fThetaR(:,:,1) = fThetaR(:,:,1) - (wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/dz &
    + nu(2)*(ThetaRhat(:,:,2)+ThetaRhat(:,:,1))/dz**2;
    fThetaR(:,:,m-1) = fThetaR(:,:,m-1) - (wrkNLhat(:,:,m)-wrkNLhat(:,:,m-1))/dz &
    + nu(2)*(ThetaRhat(:,:,m-1)+ThetaRhat(:,:,m-2))/dz**2;

    !fThetaR(1,1,2:m-2) = fThetaR(1,1,2:m-2)-nu(2)*(ThetaRhat(1,1,3:m-1)+ThetaRhat(1,1,1:m-3))/dz**2; !Remove viscosity to VSHF
    !fThetaR(1,1,1) = fThetaR(1,1,1)-nu(2)*(ThetaRhat(1,1,2)+ThetaRhat(1,1,1))/dz**2; !Remove viscosity to VSHF
    !fThetaR(1,1,m-1) = fThetaR(1,1,m-1)-nu(2)*(ThetaRhat(1,1,m-1)+ThetaRhat(1,1,m-2))/dz**2; !Remove viscosity to VSHF
    
    !----------------------------------
    !Radiative cooling
    do iz=1,m
        fThetaR(1,1,iz) = fThetaR(1,1,iz)+RC(iz); !Only one eigenmode
    end do
    
    !----------Non-linear advection term: - (u q_t')_x affecting the q_t equation --------------------------------------	
    do iz = 1,m
        in1 = u(:,:,iz)*qt(:,:,iz); !wrkNL(:,:,iz); 
        
        wrk1 = in1 
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do		
    fqt = fqt-IM*kx*wrkNLhat;
    
    !----------Non-linear advection term: - (v q_t')_y affecting the q_t equation --------------------------------------	
    do iz = 1,m
        in1 = v(:,:,iz)*qt(:,:,iz); !wrkNL(:,:,iz); 
        
        wrk1 = in1
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do	
    fqt = fqt - IM*ky*wrkNLhat;
    
    !----------Non-linear advection term: - (w q_t')_z affecting the q_t equation --------------------------------------	
    !********** Interpolation for qv **************
    wrkNL(:,:,2:m-1) = 0.5*(qt(:,:,1:m-2)+qt(:,:,2:m-1));!*** qv on position as w
    !****"zero on bottom" Neumann boundary condition **** 
    wrkNL(:,:,1) = qt(:,:,1) ;
    wrkNL(:,:,m) = qt(:,:,m-1);
    wrkNL = w*wrkNL;
    do iz = 1,m
        in1 = wrkNL(:,:,iz); 
        
        wrk1 = in1 
        wrkNLhat(:,:,iz) = wrk1/ns;
    end do	
        
    fqt(:,:,2:m-2) = fqt(:,:,2:m-2) - (wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,2:m-2))/dz &
    + nu(2)*(qthat(:,:,3:m-1)+qthat(:,:,1:m-3))/dz**2;
    !****"Non-zero on bottom" Neumann boundary condition **** 		
    fqt(:,:,1) = fqt(:,:,1) - (wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/dz &
    + nu(2)*(qthat(:,:,2)+qthat(:,:,1))/dz**2;
    fqt(:,:,m-1) = fqt(:,:,m-1) - (wrkNLhat(:,:,m)-wrkNLhat(:,:,m-1))/dz &
    + nu(2)*(qthat(:,:,m-1)+qthat(:,:,m-2))/dz**2;

    !fqt(1,1,2:m-2) = fqt(1,1,2:m-2)-nu(2)*(qthat(1,1,3:m-1)+qthat(1,1,1:m-3))/dz**2; !Remove viscosity to VSHF
    !fqt(1,1,1) = fqt(1,1,1)-nu(2)*(qthat(1,1,2)+qthat(1,1,1))/dz**2; !Remove viscosity to VSHF
    !fqt(1,1,m-1) = fqt(1,1,m-1)-nu(2)*(qthat(1,1,m-1)+qthat(1,1,m-2))/dz**2; !Remove viscosity to VSHF
    
    !----------Background term -w d_qvditle(z)/dz--------------------------------------	
    !********************* For Fixed Background
    do iz = 1,m-2
        fqt(:,:,iz) = fqt(:,:,iz)-0.5*(what(:,:,iz)+what(:,:,iz+1))*fqvdz(iz); 
    end do
    
    !-------------------------------------------------------------------------------------		
    !Moistening
    do iz=1,m
        fqt(1,1,iz) = fqt(1,1,iz)+Mr(iz); !Only the first eigenmode
    end do
    !-------------------------------------------------------------------------------------	
    ! deallocate(qr,wrk1,wrkNL,wrkNLhat);
    return
end subroutine RK_flux	


SUBROUTINE RK_FLUX_D(nx, nxh, nxph, ny, nyp, nyph, m, kx, ky, dx, dz, &
    & f_star, g_star, epsbar, l, b_star, nu, vrain, uhat, vhat, what, &
    & thetarhat, qthat, u, ud, v, vd, w, wd, thetar, thetard, qt, qtd, fu, &
    & fud, fv, fvd, fw, fwd, fthetar, fthetard, fqt, fqtd, rc, mr, ubg, vbg&
    & , tau, qvs, qvini, fqvdz, in1, in1d, out1, planf1)
      IMPLICIT NONE
    !Parameters:
    !Integers
    ! type(C_PTR),intent(in) :: planf1
      INTEGER :: ix, iy, iz, ns
    !Complex
      COMPLEX :: im
    !Real:
      REAL :: qvbar
      REAL :: qvbard
    !In:
    !Integers
      INTEGER, INTENT(IN) :: nx, nxh, nxph, ny, nyp, nyph, m
    !Real:
      REAL, INTENT(IN) :: dx, dz, vrain, tau, f_star, g_star, epsbar, l, &
    & b_star
    !*** diffusion parameters 
      REAL, DIMENSION(2), INTENT(IN) :: nu
      REAL, DIMENSION(nxph, nyp, m), INTENT(IN) :: kx, ky
      REAL, DIMENSION(m), INTENT(IN) :: qvs, fqvdz
      REAL, DIMENSION(m), INTENT(IN) :: qvini, rc, mr, ubg, vbg
      REAL, DIMENSION(nx, ny, m), INTENT(IN) :: u, v, w, thetar, qt, qr
      REAL, DIMENSION(nx, ny, m), INTENT(IN) :: ud, vd, wd, thetard, qtd, &
    & qrd
      REAL, DIMENSION(nx, ny), INTENT(INOUT) :: in1
      REAL, DIMENSION(nx, ny), INTENT(INOUT) :: in1d
    !Complex
      COMPLEX, DIMENSION(nxph, nyp, m), INTENT(IN) :: uhat, vhat, what, &
    & qthat, thetarhat, wrknlhat
      COMPLEX, DIMENSION(nxph, nyp, m), INTENT(IN) :: wrknlhatd
      COMPLEX, DIMENSION(nxh, ny), INTENT(INOUT) :: out1
      COMPLEX, DIMENSION(nx, ny, m) :: wrknl
      COMPLEX, DIMENSION(nx, ny, m) :: wrknld
      COMPLEX, DIMENSION(nxph, nyp) :: wrk1
      COMPLEX, DIMENSION(nxph, nyp) :: wrk1d
    !Out:
    !Complex:
      COMPLEX, DIMENSION(nxph, nyp, m), INTENT(OUT) :: fu, fv, fw, fthetar, &
    & fqt
      COMPLEX, DIMENSION(nxph, nyp, m), INTENT(OUT) :: fud, fvd, fwd, &
    & fthetard, fqtd
      INTRINSIC SUM
      REAL :: temp
      TYPE(UNKNOWNTYPE) :: planf1
    !---------------------------------------------
      im = (0,1)
      ns = nx*ny
    !---------------------------------------------
    !Allocation:
    !Real
    ! allocate(qr(nx,ny,m));
    !Complex
    ! allocate(wrkNL(nx,ny,m)); allocate(wrkNLhat(nxph,nyp,m));
    ! allocate(wrk1(nxph,nyp));
    !------------------------Equations to solve in non-dimensional units:------------------------
    !d u /dt = fv - (uu)_x- (vu)_y -(wu)_z + nu \nabla^2 u
    !d v /dt = -fu - (uv)_x -(vv)_y -(wv)_z + nu \nabla^2 v
    !d w /dt = -(uw)_x - (vw)_y - (ww)_z +nu \nabla^2 w + g_star theta' + g_star eps_o q_v' - g_star*q_r
    !d theta_r' /dt = -(u theta_r)_x - (v \theta_r)_y -B_star w - V_T L d q_r/dz
    !d qt' /dt = -(uq_t')_x-(vq_t')_y - (wq_t')_z -w d_qvditle(z)/dz -V_T dq_r/dz 
    !-----------Boundary conditions-------------------------------------
    !u=0 at the bottom, du/dz=0 on the top
    !**** v=0 at the bottom, dv/dz=0 on the top ****
    !Zero Dirichlet for w
    !Neumann bc for theta_r'
    !Neumann boundary condition for q_r		
    !----------- Notes for us-------------------------------------
    !q_t here is q_t', the fluctuation
    !Theta_r is theta_r', the fluctuation
    !Remove viscosity to VSHF (Vertically Shear Horizontal Flows), except for theta_r', q_t'
    !---------Coriolis---------------------------------------
      fu = f_star*vhat
      fv = -(f_star*uhat)
      wrknlhatd = (0.0,0.0)
    !---------Non-linear advection term: - (uu)_x---------------------------------------
    !wrkNL = u*u;
      DO iz=1,m
    !wrkNL(:,:,iz); !This is why we call it pseudo-spectral
        in1d = 2*u(:, :, iz)*ud(:, :, iz)
        in1 = u(:, :, iz)*u(:, :, iz)
    !de-aliasing
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fud = -(im*kx*wrknlhatd)
      fu = fu - im*kx*wrknlhat
    !---------Non-linear advection terms: - (uv)_x and - (uv)_y---------------------------------------	
    !wrkNL = u*v;
      DO iz=1,m
    !wrkNL(:,:,iz); 
        in1d = v(:, :, iz)*ud(:, :, iz) + u(:, :, iz)*vd(:, :, iz)
        in1 = u(:, :, iz)*v(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fud = fud - im*ky*wrknlhatd
      fu = fu - im*ky*wrknlhat
      fvd = -(im*kx*wrknlhatd)
      fv = fv - im*kx*wrknlhat
    !---------Non-linear advection term: - (vv)_x----------------------------------------	
    !wrkNL = v*v;
      DO iz=1,m
    !wrkNL(:,:,iz); 
        in1d = 2*v(:, :, iz)*vd(:, :, iz)
        in1 = v(:, :, iz)*v(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fvd = fvd - im*ky*wrknlhatd
      fv = fv - im*ky*wrknlhat
    !----------Non-linear advection term: - (wu)_z affecting the u equation --------------------------------------	
    !Need to compute u at the w level, due to staggered grid
    !Interpolate u:
    !*** u on position as w
      wrknld = (0.0,0.0)
      wrknld(:, :, 2:m-1) = 0.5*(ud(:, :, 1:m-2)+ud(:, :, 2:m-1))
      wrknl(:, :, 2:m-1) = 0.5*(u(:, :, 1:m-2)+u(:, :, 2:m-1))
    !**** No boundary condition, but w=0 ****
      wrknld(:, :, 1) = (0.0,0.0)
      wrknl(:, :, 1) = 0
      wrknld(:, :, m) = (0.0,0.0)
      wrknl(:, :, m) = 0
      wrknld = w*wrknld + wrknl*wd
      wrknl = wrknl*w
      DO iz=1,m
        in1d = wrknld(:, :, iz)
        in1 = wrknl(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fud(:, :, 2:m-2) = fud(:, :, 2:m-2) - (wrknlhatd(:, :, 3:m-1)-&
    &   wrknlhatd(:, :, 2:m-2))/dz
      fu(:, :, 2:m-2) = fu(:, :, 2:m-2) - (wrknlhat(:, :, 3:m-1)-wrknlhat(:&
    &   , :, 2:m-2))/dz + nu(2)*(uhat(:, :, 3:m-1)+uhat(:, :, 1:m-3))/dz**2
    !**** u=0 at the bottom, du/dz=0 on the top ****
    !u(0) = -u(1)
      fud(:, :, 1) = fud(:, :, 1) - (wrknlhatd(:, :, 2)-wrknlhatd(:, :, 1))/&
    &   dz
      fu(:, :, 1) = fu(:, :, 1) - (wrknlhat(:, :, 2)-wrknlhat(:, :, 1))/dz +&
    &   nu(2)*(uhat(:, :, 2)-uhat(:, :, 1))/dz**2
    !u(m)  = u(m-1)
      fud(:, :, m-1) = fud(:, :, m-1) - (wrknlhatd(:, :, m)-wrknlhatd(:, :, &
    &   m-1))/dz
      fu(:, :, m-1) = fu(:, :, m-1) - (wrknlhat(:, :, m)-wrknlhat(:, :, m-1)&
    &   )/dz + nu(2)*(uhat(:, :, m-1)+uhat(:, :, m-2))/dz**2
    !Remove viscosity to VSHF
      fu(1, 1, 2:m-2) = fu(1, 1, 2:m-2) - nu(2)*(uhat(1, 1, 3:m-1)+uhat(1, 1&
    &   , 1:m-3))/dz**2
    !Remove viscosity to VSHF
      fu(1, 1, 1) = fu(1, 1, 1) - nu(2)*(uhat(1, 1, 2)-uhat(1, 1, 1))/dz**2
    !Remove viscosity to VSHF
      fu(1, 1, m-1) = fu(1, 1, m-1) - nu(2)*(uhat(1, 1, m-1)+uhat(1, 1, m-2)&
    &   )/dz**2
    !**************************************	
    !!Ralaxation terms
    !!(u,v) horizontal velocity component
    !do iz=1,m
    !	fu(1,1,iz) = fu(1,1,iz)-(1/tau)*(uhat(1,1,iz)-ubg(iz));
    !	fv(1,1,iz) = fv(1,1,iz)-(1/tau)*(vhat(1,1,iz)-vbg(iz));
    !end do
    !----------Non-linear advection term: - (wu)_x affecting the w equation --------------------------------------	
      fwd = -(im*kx*wrknlhatd)
      fw = -(im*kx*wrknlhat)
    !----------Non-linear advection term: - (wv)_z affecting the v equation --------------------------------------	
    !Need to compute v at the w level, due to staggered grid
    !********* Interpolation for v ********
    !*** v on position as w
      wrknld(:, :, 2:m-1) = 0.5*(vd(:, :, 1:m-2)+vd(:, :, 2:m-1))
      wrknl(:, :, 2:m-1) = 0.5*(v(:, :, 1:m-2)+v(:, :, 2:m-1))
    !**** No boundary condition, but w=0 ****
      wrknld(:, :, 1) = (0.0,0.0)
      wrknl(:, :, 1) = 0
      wrknld(:, :, m) = (0.0,0.0)
      wrknl(:, :, m) = 0
      wrknld = w*wrknld + wrknl*wd
      wrknl = wrknl*w
    !----------------------------------
      DO iz=1,m
        in1d = wrknld(:, :, iz)
        in1 = wrknl(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
    !-----------------------------------
      fvd(:, :, 2:m-2) = fvd(:, :, 2:m-2) - (wrknlhatd(:, :, 3:m-1)-&
    &   wrknlhatd(:, :, 2:m-2))/dz
      fv(:, :, 2:m-2) = fv(:, :, 2:m-2) - (wrknlhat(:, :, 3:m-1)-wrknlhat(:&
    &   , :, 2:m-2))/dz + nu(2)*(vhat(:, :, 3:m-1)+vhat(:, :, 1:m-3))/dz**2
    !**** v=0 at the bottom, dv/dz=0 on the top ****
      fvd(:, :, 1) = fvd(:, :, 1) - (wrknlhatd(:, :, 2)-wrknlhatd(:, :, 1))/&
    &   dz
      fv(:, :, 1) = fv(:, :, 1) - (wrknlhat(:, :, 2)-wrknlhat(:, :, 1))/dz +&
    &   nu(2)*(vhat(:, :, 2)-vhat(:, :, 1))/dz**2
      fvd(:, :, m-1) = fvd(:, :, m-1) - (wrknlhatd(:, :, m)-wrknlhatd(:, :, &
    &   m-1))/dz
      fv(:, :, m-1) = fv(:, :, m-1) - (wrknlhat(:, :, m)-wrknlhat(:, :, m-1)&
    &   )/dz + nu(2)*(vhat(:, :, m-1)+vhat(:, :, m-2))/dz**2
    !Remove viscosity to VSHF
      fv(1, 1, 2:m-2) = fv(1, 1, 2:m-2) - nu(2)*(vhat(1, 1, 3:m-1)+vhat(1, 1&
    &   , 1:m-3))/dz**2
    !Remove viscosity to VSHF
      fv(1, 1, 1) = fv(1, 1, 1) - nu(2)*(vhat(1, 1, 2)-vhat(1, 1, 1))/dz**2
    !Remove viscosity to VSHF
      fv(1, 1, m-1) = fv(1, 1, m-1) - nu(2)*(vhat(1, 1, m-1)+vhat(1, 1, m-2)&
    &   )/dz**2
    !----------Non-linear advection term: - (wv)_y affecting the w equation --------------------------------------	
      fwd = fwd - im*ky*wrknlhatd
      fw = fw - im*ky*wrknlhat
    !---------------------------------------------------------------------------------
    !Here we need to get qr
      qr = 0
      qrd = 0.0
      DO iz=1,m
        WHERE (qvini(iz) + qt(:, :, iz) .GT. qvs(iz)) 
          qrd(:, :, iz) = qtd(:, :, iz)
          qr(:, :, iz) = qvini(iz) + qt(:, :, iz) - qvs(iz)
        END WHERE
      END DO
    !------------q_v '  term affecting the w equation---------------------------------------------------------------------
      DO iz=1,m
    !Horizontal average of qv'
        qvbard = SUM(qtd(:, :, iz)-qrd(:, :, iz))/ns
        qvbar = SUM(qt(:, :, iz)-qr(:, :, iz))/ns
        in1d = qtd(:, :, iz) - qrd(:, :, iz) - qvbard
        in1 = qt(:, :, iz) - qr(:, :, iz) - qvbar
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      temp = 0.5*g_star*epsbar
      fwd(:, :, 2:m-1) = fwd(:, :, 2:m-1) + temp*(wrknlhatd(:, :, 1:m-2)+&
    &   wrknlhatd(:, :, 2:m-1))
      fw(:, :, 2:m-1) = fw(:, :, 2:m-1) + temp*(wrknlhat(:, :, 1:m-2)+&
    &   wrknlhat(:, :, 2:m-1))
    !------------q_r '  term affecting the w equation---------------------------------------------------------------------
      DO iz=1,m
        in1d = qrd(:, :, iz)
        in1 = qr(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fwd(:, :, 2:m-1) = fwd(:, :, 2:m-1) - g_star*0.5*(wrknlhatd(:, :, 1:m-&
    &   2)+wrknlhatd(:, :, 2:m-1))
      fw(:, :, 2:m-1) = fw(:, :, 2:m-1) - g_star*0.5*(wrknlhat(:, :, 1:m-2)+&
    &   wrknlhat(:, :, 2:m-1))
    !------------q_r  term affecting the theta_r equation---------------------------------------------------------------------
      fthetar = 0
      fthetard = (0.0,0.0)
      fthetard(:, :, 2:m-2) = -(vrain*l*(wrknlhatd(:, :, 3:m-1)-wrknlhatd(:&
    &   , :, 1:m-3))/(2*dz))
      fthetar(:, :, 2:m-2) = fthetar(:, :, 2:m-2) - vrain*l*(wrknlhat(:, :, &
    &   3:m-1)-wrknlhat(:, :, 1:m-3))/(2*dz)
    !**** Neumann boundary condition ****
      fthetard(:, :, 1) = fthetard(:, :, 1) - vrain*l*(wrknlhatd(:, :, 2)-&
    &   wrknlhatd(:, :, 1))/(2*dz)
      fthetar(:, :, 1) = fthetar(:, :, 1) - vrain*l*(wrknlhat(:, :, 2)-&
    &   wrknlhat(:, :, 1))/(2*dz)
      fthetard(:, :, m-1) = fthetard(:, :, m-1) - vrain*l*(wrknlhatd(:, :, m&
    &   -1)-wrknlhatd(:, :, m-2))/(2*dz)
      fthetar(:, :, m-1) = fthetar(:, :, m-1) - vrain*l*(wrknlhat(:, :, m-1)&
    &   -wrknlhat(:, :, m-2))/(2*dz)
    !------------q_r  term affecting the q_t equation---------------------------------------------------------------------
      fqt = 0
      fqtd = (0.0,0.0)
      fqtd(:, :, 2:m-2) = vrain*(wrknlhatd(:, :, 3:m-1)-wrknlhatd(:, :, 1:m-&
    &   3))/(2*dz)
      fqt(:, :, 2:m-2) = fqt(:, :, 2:m-2) + vrain*(wrknlhat(:, :, 3:m-1)-&
    &   wrknlhat(:, :, 1:m-3))/(2*dz)
    !**** Neumann boundary condition for q_r		
    !q_r(0) = q_r(1)
      fqtd(:, :, 1) = fqtd(:, :, 1) + vrain*(wrknlhatd(:, :, 2)-wrknlhatd(:&
    &   , :, 1))/(2*dz)
      fqt(:, :, 1) = fqt(:, :, 1) + vrain*(wrknlhat(:, :, 2)-wrknlhat(:, :, &
    &   1))/(2*dz)
    !q_r(m) = q_r(m-1)
      fqtd(:, :, m-1) = fqtd(:, :, m-1) + vrain*(wrknlhatd(:, :, m-1)-&
    &   wrknlhatd(:, :, m-2))/(2*dz)
      fqt(:, :, m-1) = fqt(:, :, m-1) + vrain*(wrknlhat(:, :, m-1)-wrknlhat(&
    &   :, :, m-2))/(2*dz)
    !-------------Theta term affecting the w equation--------------------------------------------------------------------
      DO iz=1,m
        in1d = thetard(:, :, iz) + l*qrd(:, :, iz)
        in1 = thetar(:, :, iz) + l*qr(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fwd(:, :, 2:m-1) = fwd(:, :, 2:m-1) + g_star*0.5*(wrknlhatd(:, :, 1:m-&
    &   2)+wrknlhatd(:, :, 2:m-1))
      fw(:, :, 2:m-1) = fw(:, :, 2:m-1) + g_star*(0.5*(wrknlhat(:, :, 1:m-2)&
    &   +wrknlhat(:, :, 2:m-1)))
    !----------Non-linear advection term: - (ww)_z affecting the w equation --------------------------------------	
    !********** Interpolation for w ***************
    !*** w on position as u
      wrknld(:, :, 1:m-1) = 0.5*(wd(:, :, 1:m-1)+wd(:, :, 2:m))
      wrknl(:, :, 1:m-1) = 0.5*(w(:, :, 1:m-1)+w(:, :, 2:m))
    !**** Zero boundary condition, 0.5*(w(m+1)+w(m-1))=0, w(m)=0.0
      wrknld(:, :, m) = -wrknld(:, :, m-1)
      wrknl(:, :, m) = -wrknl(:, :, m-1)
      wrknld = 2*wrknl*wrknld
      wrknl = wrknl*wrknl
      DO iz=1,m
        in1d = wrknld(:, :, iz)
        in1 = wrknl(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fwd(:, :, 2:m-1) = fwd(:, :, 2:m-1) - (wrknlhatd(:, :, 2:m-1)-&
    &   wrknlhatd(:, :, 1:m-2))/dz
      fw(:, :, 2:m-1) = fw(:, :, 2:m-1) - (wrknlhat(:, :, 2:m-1)-wrknlhat(:&
    &   , :, 1:m-2))/dz + nu(2)*(what(:, :, 3:m)+what(:, :, 1:m-2))/dz**2
    !fw(1,1,2:m-1) = fw(1,1,2:m-1)-nu(2)*(what(1,1,3:m)+what(1,1,1:m-2))/dz**2; !Remove viscosity to VSHF
    !**** Zero boundary condition ******	
      fwd(:, :, 1) = (0.0,0.0)
      fw(:, :, 1) = 0
      fwd(:, :, m) = (0.0,0.0)
      fw(:, :, m) = 0
    !----------Non-linear advection term: - (u theta_r')_x affecting the theta_r equation --------------------------------------	
    !wrkNL = u*ThetaR;
      DO iz=1,m
    !wrkNL(:,:,iz); 
        in1d = thetar(:, :, iz)*ud(:, :, iz) + u(:, :, iz)*thetard(:, :, iz)
        in1 = u(:, :, iz)*thetar(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fthetard = fthetard - im*kx*wrknlhatd
      fthetar = fthetar - im*kx*wrknlhat
    !----------Non-linear advection term: - (v theta_r')_y affecting the theta_r equation --------------------------------------	
    !wrkNL = v*ThetaR;
      DO iz=1,m
    !wrkNL(:,:,iz); 
        in1d = thetar(:, :, iz)*vd(:, :, iz) + v(:, :, iz)*thetard(:, :, iz)
        in1 = v(:, :, iz)*thetar(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fthetard = fthetard - im*ky*wrknlhatd
      fthetar = fthetar - im*ky*wrknlhat
    !---------------------------------------------------------------------------------
    !********************* For Fixed Background
      DO iz=1,m-2
        fthetar(:, :, iz) = fthetar(:, :, iz) - 0.5*(what(:, :, iz)+what(:, &
    &     :, iz+1))*b_star
      END DO
    !----------Non-linear advection term: - (w theta_r')_z affecting the theta_r equation --------------------------------------	
    !********** Interpolation for ThetaR **************
    !*** ThetaR on position as w
      wrknld(:, :, 2:m-1) = 0.5*(thetard(:, :, 1:m-2)+thetard(:, :, 2:m-1))
      wrknl(:, :, 2:m-1) = 0.5*(thetar(:, :, 1:m-2)+thetar(:, :, 2:m-1))
    !****   Neumann boundary condition ****
      wrknld(:, :, 1) = thetard(:, :, 1)
      wrknl(:, :, 1) = thetar(:, :, 1)
      wrknld(:, :, m) = thetard(:, :, m-1)
      wrknl(:, :, m) = thetar(:, :, m-1)
      wrknld = wrknl*wd + w*wrknld
      wrknl = w*wrknl
    !----------------------------------
      DO iz=1,m
        in1d = wrknld(:, :, iz)
        in1 = wrknl(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
    !**************************************	
      fthetard(:, :, 2:m-2) = fthetard(:, :, 2:m-2) - (wrknlhatd(:, :, 3:m-1&
    &   )-wrknlhatd(:, :, 2:m-2))/dz
      fthetar(:, :, 2:m-2) = fthetar(:, :, 2:m-2) - (wrknlhat(:, :, 3:m-1)-&
    &   wrknlhat(:, :, 2:m-2))/dz + nu(2)*(thetarhat(:, :, 3:m-1)+thetarhat(&
    &   :, :, 1:m-3))/dz**2
    !**** Neumann boundary condition ****
      fthetard(:, :, 1) = fthetard(:, :, 1) - (wrknlhatd(:, :, 2)-wrknlhatd(&
    &   :, :, 1))/dz
      fthetar(:, :, 1) = fthetar(:, :, 1) - (wrknlhat(:, :, 2)-wrknlhat(:, :&
    &   , 1))/dz + nu(2)*(thetarhat(:, :, 2)+thetarhat(:, :, 1))/dz**2
      fthetard(:, :, m-1) = fthetard(:, :, m-1) - (wrknlhatd(:, :, m)-&
    &   wrknlhatd(:, :, m-1))/dz
      fthetar(:, :, m-1) = fthetar(:, :, m-1) - (wrknlhat(:, :, m)-wrknlhat(&
    &   :, :, m-1))/dz + nu(2)*(thetarhat(:, :, m-1)+thetarhat(:, :, m-2))/&
    &   dz**2
    !fThetaR(1,1,2:m-2) = fThetaR(1,1,2:m-2)-nu(2)*(ThetaRhat(1,1,3:m-1)+ThetaRhat(1,1,1:m-3))/dz**2; !Remove viscosity to VSHF
    !fThetaR(1,1,1) = fThetaR(1,1,1)-nu(2)*(ThetaRhat(1,1,2)+ThetaRhat(1,1,1))/dz**2; !Remove viscosity to VSHF
    !fThetaR(1,1,m-1) = fThetaR(1,1,m-1)-nu(2)*(ThetaRhat(1,1,m-1)+ThetaRhat(1,1,m-2))/dz**2; !Remove viscosity to VSHF
    !----------------------------------
    !Radiative cooling
      DO iz=1,m
    !Only one eigenmode
        fthetar(1, 1, iz) = fthetar(1, 1, iz) + rc(iz)
      END DO
    !----------Non-linear advection term: - (u q_t')_x affecting the q_t equation --------------------------------------	
      DO iz=1,m
    !wrkNL(:,:,iz); 
        in1d = qt(:, :, iz)*ud(:, :, iz) + u(:, :, iz)*qtd(:, :, iz)
        in1 = u(:, :, iz)*qt(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fqtd = fqtd - im*kx*wrknlhatd
      fqt = fqt - im*kx*wrknlhat
    !----------Non-linear advection term: - (v q_t')_y affecting the q_t equation --------------------------------------	
      DO iz=1,m
    !wrkNL(:,:,iz); 
        in1d = qt(:, :, iz)*vd(:, :, iz) + v(:, :, iz)*qtd(:, :, iz)
        in1 = v(:, :, iz)*qt(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fqtd = fqtd - im*ky*wrknlhatd
      fqt = fqt - im*ky*wrknlhat
    !----------Non-linear advection term: - (w q_t')_z affecting the q_t equation --------------------------------------	
    !********** Interpolation for qv **************
    !*** qv on position as w
      wrknld(:, :, 2:m-1) = 0.5*(qtd(:, :, 1:m-2)+qtd(:, :, 2:m-1))
      wrknl(:, :, 2:m-1) = 0.5*(qt(:, :, 1:m-2)+qt(:, :, 2:m-1))
    !****"zero on bottom" Neumann boundary condition **** 
      wrknld(:, :, 1) = qtd(:, :, 1)
      wrknl(:, :, 1) = qt(:, :, 1)
      wrknld(:, :, m) = qtd(:, :, m-1)
      wrknl(:, :, m) = qt(:, :, m-1)
      wrknld = wrknl*wd + w*wrknld
      wrknl = w*wrknl
      DO iz=1,m
        in1d = wrknld(:, :, iz)
        in1 = wrknl(:, :, iz)
        wrk1d = in1d
        wrk1 = in1
        wrknlhatd(:, :, iz) = wrk1d/ns
        wrknlhat(:, :, iz) = wrk1/ns
      END DO
      fqtd(:, :, 2:m-2) = fqtd(:, :, 2:m-2) - (wrknlhatd(:, :, 3:m-1)-&
    &   wrknlhatd(:, :, 2:m-2))/dz
      fqt(:, :, 2:m-2) = fqt(:, :, 2:m-2) - (wrknlhat(:, :, 3:m-1)-wrknlhat(&
    &   :, :, 2:m-2))/dz + nu(2)*(qthat(:, :, 3:m-1)+qthat(:, :, 1:m-3))/dz&
    &   **2
    !****"Non-zero on bottom" Neumann boundary condition **** 		
      fqtd(:, :, 1) = fqtd(:, :, 1) - (wrknlhatd(:, :, 2)-wrknlhatd(:, :, 1)&
    &   )/dz
      fqt(:, :, 1) = fqt(:, :, 1) - (wrknlhat(:, :, 2)-wrknlhat(:, :, 1))/dz&
    &   + nu(2)*(qthat(:, :, 2)+qthat(:, :, 1))/dz**2
      fqtd(:, :, m-1) = fqtd(:, :, m-1) - (wrknlhatd(:, :, m)-wrknlhatd(:, :&
    &   , m-1))/dz
      fqt(:, :, m-1) = fqt(:, :, m-1) - (wrknlhat(:, :, m)-wrknlhat(:, :, m-&
    &   1))/dz + nu(2)*(qthat(:, :, m-1)+qthat(:, :, m-2))/dz**2
    !fqt(1,1,2:m-2) = fqt(1,1,2:m-2)-nu(2)*(qthat(1,1,3:m-1)+qthat(1,1,1:m-3))/dz**2; !Remove viscosity to VSHF
    !fqt(1,1,1) = fqt(1,1,1)-nu(2)*(qthat(1,1,2)+qthat(1,1,1))/dz**2; !Remove viscosity to VSHF
    !fqt(1,1,m-1) = fqt(1,1,m-1)-nu(2)*(qthat(1,1,m-1)+qthat(1,1,m-2))/dz**2; !Remove viscosity to VSHF
    !----------Background term -w d_qvditle(z)/dz--------------------------------------	
    !********************* For Fixed Background
      DO iz=1,m-2
        fqt(:, :, iz) = fqt(:, :, iz) - 0.5*(what(:, :, iz)+what(:, :, iz+1)&
    &     )*fqvdz(iz)
      END DO
    !-------------------------------------------------------------------------------------		
    !Moistening
      DO iz=1,m
    !Only the first eigenmode
        fqt(1, 1, iz) = fqt(1, 1, iz) + mr(iz)
      END DO
    !-------------------------------------------------------------------------------------	
    ! deallocate(qr,wrk1,wrkNL,wrkNLhat);
      RETURN
END SUBROUTINE RK_FLUX_D


SUBROUTINE RK_FLUX_B(nx, nxh, nxph, ny, nyp, nyph, m, kx, ky, dx, dz, &
    & f_star, g_star, epsbar, l, b_star, nu, vrain, uhat, vhat, what, &
    & thetarhat, qthat, u, ub, v, vb, w, wb, thetar, thetarb, qt, qtb, fu, &
    & fub, fv, fvb, fw, fwb, fthetar, fthetarb, fqt, fqtb, rc, mr, ubg, vbg&
    & , tau, qvs, qvini, fqvdz, in1, in1b, out1, planf1)
      IMPLICIT NONE
    !Parameters:
    !Integers
    ! type(C_PTR),intent(in) :: planf1
      INTEGER :: ix, iy, iz, ns
    !Complex
      COMPLEX :: im
    !Real:
      REAL :: qvbar
      REAL :: qvbarb
    !In:
    !Integers
      INTEGER, INTENT(IN) :: nx, nxh, nxph, ny, nyp, nyph, m
    !Real:
      REAL, INTENT(IN) :: dx, dz, vrain, tau, f_star, g_star, epsbar, l, &
    & b_star
    !*** diffusion parameters 
      REAL, DIMENSION(2), INTENT(IN) :: nu
      REAL, DIMENSION(nxph, nyp, m), INTENT(IN) :: kx, ky
      REAL, DIMENSION(m), INTENT(IN) :: qvs, fqvdz
      REAL, DIMENSION(m), INTENT(IN) :: qvini, rc, mr, ubg, vbg
      REAL, DIMENSION(nx, ny, m), INTENT(IN) :: u, v, w, thetar, qt, qr
      REAL, DIMENSION(nx, ny, m) :: ub, vb, wb, thetarb, qtb, qrb
      REAL, DIMENSION(nx, ny), INTENT(INOUT) :: in1
      REAL, DIMENSION(nx, ny), INTENT(INOUT) :: in1b
    !Complex
      COMPLEX, DIMENSION(nxph, nyp, m), INTENT(IN) :: uhat, vhat, what, &
    & qthat, thetarhat, wrknlhat
      COMPLEX, DIMENSION(nxph, nyp, m) :: wrknlhatb
      COMPLEX, DIMENSION(nxh, ny), INTENT(INOUT) :: out1
      COMPLEX, DIMENSION(nx, ny, m) :: wrknl
      COMPLEX, DIMENSION(nx, ny, m) :: wrknlb
      COMPLEX, DIMENSION(nxph, nyp) :: wrk1
      COMPLEX, DIMENSION(nxph, nyp) :: wrk1b
    !Out:
    !Complex:
      COMPLEX, DIMENSION(nxph, nyp, m) :: fu, fv, fw, fthetar, fqt
      COMPLEX, DIMENSION(nxph, nyp, m) :: fub, fvb, fwb, fthetarb, fqtb
      INTRINSIC SUM
      LOGICAL, DIMENSION(nx, ny) :: mask
      REAL :: temp
      COMPLEX, DIMENSION(nxph, nyp, m-3) :: tempb
      COMPLEX, DIMENSION(nxph, nyp) :: tempb0
      TYPE(UNKNOWNTYPE) :: planf1
    !---------------------------------------------
      im = (0,1)
      ns = nx*ny
    !---------------------------------------------
    !Allocation:
    !Real
    ! allocate(qr(nx,ny,m));
    !Complex
    ! allocate(wrkNL(nx,ny,m)); allocate(wrkNLhat(nxph,nyp,m));
    ! allocate(wrk1(nxph,nyp));
    !------------------------Equations to solve in non-dimensional units:------------------------
    !d u /dt = fv - (uu)_x- (vu)_y -(wu)_z + nu \nabla^2 u
    !d v /dt = -fu - (uv)_x -(vv)_y -(wv)_z + nu \nabla^2 v
    !d w /dt = -(uw)_x - (vw)_y - (ww)_z +nu \nabla^2 w + g_star theta' + g_star eps_o q_v' - g_star*q_r
    !d theta_r' /dt = -(u theta_r)_x - (v \theta_r)_y -B_star w - V_T L d q_r/dz
    !d qt' /dt = -(uq_t')_x-(vq_t')_y - (wq_t')_z -w d_qvditle(z)/dz -V_T dq_r/dz 
    !-----------Boundary conditions-------------------------------------
    !u=0 at the bottom, du/dz=0 on the top
    !**** v=0 at the bottom, dv/dz=0 on the top ****
    !Zero Dirichlet for w
    !Neumann bc for theta_r'
    !Neumann boundary condition for q_r		
    !----------- Notes for us-------------------------------------
    !q_t here is q_t', the fluctuation
    !Theta_r is theta_r', the fluctuation
    !Remove viscosity to VSHF (Vertically Shear Horizontal Flows), except for theta_r', q_t'
    !---------Coriolis---------------------------------------
    !----------Non-linear advection term: - (wu)_z affecting the u equation --------------------------------------	
    !Need to compute u at the w level, due to staggered grid
    !Interpolate u:
    !*** u on position as w
      wrknl(:, :, 2:m-1) = 0.5*(u(:, :, 1:m-2)+u(:, :, 2:m-1))
    !**** No boundary condition, but w=0 ****
      wrknl(:, :, 1) = 0
      wrknl(:, :, m) = 0
      CALL PUSHCOMPLEX8ARRAY(wrknl, nx*ny*m)
      wrknl = wrknl*w
    !**** u=0 at the bottom, du/dz=0 on the top ****
    !u(0) = -u(1)
    !u(m)  = u(m-1)
    !Remove viscosity to VSHF
    !Remove viscosity to VSHF
    !Remove viscosity to VSHF
    !**************************************	
    !!Ralaxation terms
    !!(u,v) horizontal velocity component
    !do iz=1,m
    !	fu(1,1,iz) = fu(1,1,iz)-(1/tau)*(uhat(1,1,iz)-ubg(iz));
    !	fv(1,1,iz) = fv(1,1,iz)-(1/tau)*(vhat(1,1,iz)-vbg(iz));
    !end do
    !----------Non-linear advection term: - (wu)_x affecting the w equation --------------------------------------	
    !----------Non-linear advection term: - (wv)_z affecting the v equation --------------------------------------	
    !Need to compute v at the w level, due to staggered grid
    !********* Interpolation for v ********
    !*** v on position as w
      wrknl(:, :, 2:m-1) = 0.5*(v(:, :, 1:m-2)+v(:, :, 2:m-1))
    !**** No boundary condition, but w=0 ****
      wrknl(:, :, 1) = 0
      wrknl(:, :, m) = 0
      CALL PUSHCOMPLEX8ARRAY(wrknl, nx*ny*m)
      wrknl = wrknl*w
    !-----------------------------------
    !**** v=0 at the bottom, dv/dz=0 on the top ****
    !Remove viscosity to VSHF
    !Remove viscosity to VSHF
    !Remove viscosity to VSHF
    !----------Non-linear advection term: - (wv)_y affecting the w equation --------------------------------------	
    !---------------------------------------------------------------------------------
    !Here we need to get qr
      qr = 0
      DO iz=1,m
        mask(:, :) = qvini(iz) + qt(:, :, iz) .GT. qvs(iz)
        WHERE (mask(:, :)) qr(:, :, iz) = qvini(iz) + qt(:, :, iz) - qvs(iz)
      END DO
    !----------Non-linear advection term: - (ww)_z affecting the w equation --------------------------------------	
    !********** Interpolation for w ***************
    !*** w on position as u
      wrknl(:, :, 1:m-1) = 0.5*(w(:, :, 1:m-1)+w(:, :, 2:m))
    !**** Zero boundary condition, 0.5*(w(m+1)+w(m-1))=0, w(m)=0.0
      wrknl(:, :, m) = -wrknl(:, :, m-1)
      CALL PUSHCOMPLEX8ARRAY(wrknl, nx*ny*m)
      wrknl = wrknl*wrknl
    !----------Non-linear advection term: - (w theta_r')_z affecting the theta_r equation --------------------------------------	
    !********** Interpolation for ThetaR **************
    !*** ThetaR on position as w
      wrknl(:, :, 2:m-1) = 0.5*(thetar(:, :, 1:m-2)+thetar(:, :, 2:m-1))
    !****   Neumann boundary condition ****
      wrknl(:, :, 1) = thetar(:, :, 1)
      wrknl(:, :, m) = thetar(:, :, m-1)
      CALL PUSHCOMPLEX8ARRAY(wrknl, nx*ny*m)
      wrknl = w*wrknl
    !----------Non-linear advection term: - (w q_t')_z affecting the q_t equation --------------------------------------	
    !********** Interpolation for qv **************
    !*** qv on position as w
      wrknl(:, :, 2:m-1) = 0.5*(qt(:, :, 1:m-2)+qt(:, :, 2:m-1))
    !****"zero on bottom" Neumann boundary condition **** 
      wrknl(:, :, 1) = qt(:, :, 1)
      wrknl(:, :, m) = qt(:, :, m-1)
      wrknlhatb = (0.0,0.0)
      wrknlhatb(:, :, m) = wrknlhatb(:, :, m) - fqtb(:, :, m-1)/dz
      wrknlhatb(:, :, m-1) = wrknlhatb(:, :, m-1) + fqtb(:, :, m-1)/dz
      wrknlhatb(:, :, 2) = wrknlhatb(:, :, 2) - fqtb(:, :, 1)/dz
      wrknlhatb(:, :, 1) = wrknlhatb(:, :, 1) + fqtb(:, :, 1)/dz
      wrknlhatb(:, :, 3:m-1) = wrknlhatb(:, :, 3:m-1) - fqtb(:, :, 2:m-2)/dz
      wrknlhatb(:, :, 2:m-2) = wrknlhatb(:, :, 2:m-2) + fqtb(:, :, 2:m-2)/dz
      wrknlb = (0.0,0.0)
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        wrknlb(:, :, iz) = wrknlb(:, :, iz) + in1b
      END DO
      wb = 0.0
      wb = REAL(CONJG(wrknl)*wrknlb)
      wrknlb = w*wrknlb
      qtb = 0.0
      qtb(:, :, m-1) = qtb(:, :, m-1) + REAL(wrknlb(:, :, m))
      wrknlb(:, :, m) = (0.0,0.0)
      qtb(:, :, 1) = qtb(:, :, 1) + REAL(wrknlb(:, :, 1))
      wrknlb(:, :, 1) = (0.0,0.0)
      qtb(:, :, 1:m-2) = qtb(:, :, 1:m-2) + REAL(0.5*wrknlb(:, :, 2:m-1))
      qtb(:, :, 2:m-1) = qtb(:, :, 2:m-1) + REAL(0.5*wrknlb(:, :, 2:m-1))
      wrknlb(:, :, 2:m-1) = (0.0,0.0)
      wrknlhatb = wrknlhatb + CONJG(-(im*ky))*fqtb
      vb = 0.0
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        vb(:, :, iz) = vb(:, :, iz) + qt(:, :, iz)*in1b
        qtb(:, :, iz) = qtb(:, :, iz) + v(:, :, iz)*in1b
      END DO
      wrknlhatb = wrknlhatb + CONJG(-(im*kx))*fqtb
      ub = 0.0
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        ub(:, :, iz) = ub(:, :, iz) + qt(:, :, iz)*in1b
        qtb(:, :, iz) = qtb(:, :, iz) + u(:, :, iz)*in1b
      END DO
      wrknlhatb(:, :, m) = wrknlhatb(:, :, m) - fthetarb(:, :, m-1)/dz
      wrknlhatb(:, :, m-1) = wrknlhatb(:, :, m-1) + fthetarb(:, :, m-1)/dz
      wrknlhatb(:, :, 2) = wrknlhatb(:, :, 2) - fthetarb(:, :, 1)/dz
      wrknlhatb(:, :, 1) = wrknlhatb(:, :, 1) + fthetarb(:, :, 1)/dz
      wrknlhatb(:, :, 3:m-1) = wrknlhatb(:, :, 3:m-1) - fthetarb(:, :, 2:m-2&
    &   )/dz
      wrknlhatb(:, :, 2:m-2) = wrknlhatb(:, :, 2:m-2) + fthetarb(:, :, 2:m-2&
    &   )/dz
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        wrknlb(:, :, iz) = wrknlb(:, :, iz) + in1b
      END DO
      CALL POPCOMPLEX8ARRAY(wrknl, nx*ny*m)
      wb = wb + REAL(CONJG(wrknl)*wrknlb)
      wrknlb = w*wrknlb
      thetarb = 0.0
      thetarb(:, :, m-1) = thetarb(:, :, m-1) + REAL(wrknlb(:, :, m))
      wrknlb(:, :, m) = (0.0,0.0)
      thetarb(:, :, 1) = thetarb(:, :, 1) + REAL(wrknlb(:, :, 1))
      wrknlb(:, :, 1) = (0.0,0.0)
      thetarb(:, :, 1:m-2) = thetarb(:, :, 1:m-2) + REAL(0.5*wrknlb(:, :, 2:&
    &   m-1))
      thetarb(:, :, 2:m-1) = thetarb(:, :, 2:m-1) + REAL(0.5*wrknlb(:, :, 2:&
    &   m-1))
      wrknlb(:, :, 2:m-1) = (0.0,0.0)
      wrknlhatb = wrknlhatb + CONJG(-(im*ky))*fthetarb
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        vb(:, :, iz) = vb(:, :, iz) + thetar(:, :, iz)*in1b
        thetarb(:, :, iz) = thetarb(:, :, iz) + v(:, :, iz)*in1b
      END DO
      wrknlhatb = wrknlhatb + CONJG(-(im*kx))*fthetarb
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        ub(:, :, iz) = ub(:, :, iz) + thetar(:, :, iz)*in1b
        thetarb(:, :, iz) = thetarb(:, :, iz) + u(:, :, iz)*in1b
      END DO
      fwb(:, :, m) = (0.0,0.0)
      fwb(:, :, 1) = (0.0,0.0)
      wrknlhatb(:, :, 2:m-1) = wrknlhatb(:, :, 2:m-1) - fwb(:, :, 2:m-1)/dz
      wrknlhatb(:, :, 1:m-2) = wrknlhatb(:, :, 1:m-2) + fwb(:, :, 2:m-1)/dz
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        wrknlb(:, :, iz) = wrknlb(:, :, iz) + in1b
      END DO
      CALL POPCOMPLEX8ARRAY(wrknl, nx*ny*m)
      wrknlb = CONJG(2*wrknl)*wrknlb
      wrknlb(:, :, m-1) = wrknlb(:, :, m-1) - wrknlb(:, :, m)
      wrknlb(:, :, m) = (0.0,0.0)
      wb(:, :, 1:m-1) = wb(:, :, 1:m-1) + REAL(0.5*wrknlb(:, :, 1:m-1))
      wb(:, :, 2:m) = wb(:, :, 2:m) + REAL(0.5*wrknlb(:, :, 1:m-1))
      wrknlb(:, :, 1:m-1) = (0.0,0.0)
      wrknlhatb(:, :, 1:m-2) = wrknlhatb(:, :, 1:m-2) + g_star*0.5*fwb(:, :&
    &   , 2:m-1)
      wrknlhatb(:, :, 2:m-1) = wrknlhatb(:, :, 2:m-1) + g_star*0.5*fwb(:, :&
    &   , 2:m-1)
      qrb = 0.0
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        thetarb(:, :, iz) = thetarb(:, :, iz) + in1b
        qrb(:, :, iz) = qrb(:, :, iz) + l*in1b
      END DO
      tempb0 = vrain*fqtb(:, :, m-1)/(2*dz)
      wrknlhatb(:, :, m-1) = wrknlhatb(:, :, m-1) + tempb0
      wrknlhatb(:, :, m-2) = wrknlhatb(:, :, m-2) - tempb0
      tempb0 = vrain*fqtb(:, :, 1)/(2*dz)
      wrknlhatb(:, :, 2) = wrknlhatb(:, :, 2) + tempb0
      wrknlhatb(:, :, 1) = wrknlhatb(:, :, 1) - tempb0
      wrknlhatb(:, :, 3:m-1) = wrknlhatb(:, :, 3:m-1) + vrain*fqtb(:, :, 2:m&
    &   -2)/(2*dz)
      wrknlhatb(:, :, 1:m-3) = wrknlhatb(:, :, 1:m-3) - vrain*fqtb(:, :, 2:m&
    &   -2)/(2*dz)
      tempb0 = -(vrain*l*fthetarb(:, :, m-1)/(2*dz))
      wrknlhatb(:, :, m-1) = wrknlhatb(:, :, m-1) + tempb0
      wrknlhatb(:, :, m-2) = wrknlhatb(:, :, m-2) - tempb0
      tempb0 = -(vrain*l*fthetarb(:, :, 1)/(2*dz))
      wrknlhatb(:, :, 2) = wrknlhatb(:, :, 2) + tempb0
      wrknlhatb(:, :, 1) = wrknlhatb(:, :, 1) - tempb0
      tempb = -(vrain*l*fthetarb(:, :, 2:m-2)/(2*dz))
      wrknlhatb(:, :, 3:m-1) = wrknlhatb(:, :, 3:m-1) + tempb
      wrknlhatb(:, :, 1:m-3) = wrknlhatb(:, :, 1:m-3) - tempb
      wrknlhatb(:, :, 1:m-2) = wrknlhatb(:, :, 1:m-2) - g_star*0.5*fwb(:, :&
    &   , 2:m-1)
      wrknlhatb(:, :, 2:m-1) = wrknlhatb(:, :, 2:m-1) - g_star*0.5*fwb(:, :&
    &   , 2:m-1)
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        qrb(:, :, iz) = qrb(:, :, iz) + in1b
      END DO
      temp = 0.5*g_star*epsbar
      wrknlhatb(:, :, 1:m-2) = wrknlhatb(:, :, 1:m-2) + temp*fwb(:, :, 2:m-1&
    &   )
      wrknlhatb(:, :, 2:m-1) = wrknlhatb(:, :, 2:m-1) + temp*fwb(:, :, 2:m-1&
    &   )
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        qvbarb = -SUM(in1b)
        qtb(:, :, iz) = qtb(:, :, iz) + in1b + qvbarb/ns
        qrb(:, :, iz) = qrb(:, :, iz) - in1b - qvbarb/ns
      END DO
      DO iz=m,1,-1
        mask(:, :) = qvini(iz) + qt(:, :, iz) .GT. qvs(iz)
        WHERE (mask(:, :)) 
          qtb(:, :, iz) = qtb(:, :, iz) + qrb(:, :, iz)
          qrb(:, :, iz) = 0.0
        END WHERE
      END DO
      wrknlhatb = wrknlhatb + CONJG(-(im*ky))*fwb
      wrknlhatb(:, :, m) = wrknlhatb(:, :, m) - fvb(:, :, m-1)/dz
      wrknlhatb(:, :, m-1) = wrknlhatb(:, :, m-1) + fvb(:, :, m-1)/dz
      wrknlhatb(:, :, 2) = wrknlhatb(:, :, 2) - fvb(:, :, 1)/dz
      wrknlhatb(:, :, 1) = wrknlhatb(:, :, 1) + fvb(:, :, 1)/dz
      wrknlhatb(:, :, 3:m-1) = wrknlhatb(:, :, 3:m-1) - fvb(:, :, 2:m-2)/dz
      wrknlhatb(:, :, 2:m-2) = wrknlhatb(:, :, 2:m-2) + fvb(:, :, 2:m-2)/dz
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        wrknlb(:, :, iz) = wrknlb(:, :, iz) + in1b
      END DO
      CALL POPCOMPLEX8ARRAY(wrknl, nx*ny*m)
      wb = wb + REAL(CONJG(wrknl)*wrknlb)
      wrknlb = w*wrknlb
      wrknlb(:, :, m) = (0.0,0.0)
      wrknlb(:, :, 1) = (0.0,0.0)
      vb(:, :, 1:m-2) = vb(:, :, 1:m-2) + REAL(0.5*wrknlb(:, :, 2:m-1))
      vb(:, :, 2:m-1) = vb(:, :, 2:m-1) + REAL(0.5*wrknlb(:, :, 2:m-1))
      wrknlb(:, :, 2:m-1) = (0.0,0.0)
      wrknlhatb = wrknlhatb + CONJG(-(im*kx))*fwb
      wrknlhatb(:, :, m) = wrknlhatb(:, :, m) - fub(:, :, m-1)/dz
      wrknlhatb(:, :, m-1) = wrknlhatb(:, :, m-1) + fub(:, :, m-1)/dz
      wrknlhatb(:, :, 2) = wrknlhatb(:, :, 2) - fub(:, :, 1)/dz
      wrknlhatb(:, :, 1) = wrknlhatb(:, :, 1) + fub(:, :, 1)/dz
      wrknlhatb(:, :, 3:m-1) = wrknlhatb(:, :, 3:m-1) - fub(:, :, 2:m-2)/dz
      wrknlhatb(:, :, 2:m-2) = wrknlhatb(:, :, 2:m-2) + fub(:, :, 2:m-2)/dz
      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        wrknlb(:, :, iz) = wrknlb(:, :, iz) + in1b
      END DO
      CALL POPCOMPLEX8ARRAY(wrknl, nx*ny*m)
      wb = wb + REAL(CONJG(wrknl)*wrknlb)
      wrknlb = w*wrknlb
      wrknlb(:, :, m) = (0.0,0.0)
      wrknlb(:, :, 1) = (0.0,0.0)


      ub(:, :, 1:m-2) = ub(:, :, 1:m-2) + REAL(0.5*wrknlb(:, :, 2:m-1))
      ub(:, :, 2:m-1) = ub(:, :, 2:m-1) + REAL(0.5*wrknlb(:, :, 2:m-1))
      wrknlhatb = wrknlhatb + CONJG(-(im*ky))*fvb

      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        vb(:, :, iz) = vb(:, :, iz) + 2*v(:, :, iz)*in1b
      END DO
      wrknlhatb = wrknlhatb + CONJG(-(im*kx))*fvb + CONJG(-(im*ky))*fub

      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        ub(:, :, iz) = ub(:, :, iz) + v(:, :, iz)*in1b
        vb(:, :, iz) = vb(:, :, iz) + u(:, :, iz)*in1b
      END DO
      wrknlhatb = wrknlhatb + CONJG(-(im*kx))*fub

      DO iz=m,1,-1
        wrk1b = (0.0,0.0)
        wrk1b = wrknlhatb(:, :, iz)/ns
        wrknlhatb(:, :, iz) = (0.0,0.0)
        in1b = 0.0
        in1b = REAL(wrk1b)
        ub(:, :, iz) = ub(:, :, iz) + 2*u(:, :, iz)*in1b
      END DO
    !   fthetarb = (0.0,0.0)
    !   fqtb = (0.0,0.0)
    !   fub = (0.0,0.0)
    !   fvb = (0.0,0.0)
    !   fwb = (0.0,0.0)
END SUBROUTINE RK_FLUX_B