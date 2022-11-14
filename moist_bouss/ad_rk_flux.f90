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
  type(C_PTR), intent(in) :: planf1
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
  real,dimension(nx,ny,m),intent(in) :: u,v,w,ThetaR, qt
  real, dimension(nx,ny,m) :: qr
  real ,dimension(nx,ny),intent(inout) ::  in1
  !Complex
  complex, dimension(nxph,nyp,m), intent(in) :: uhat,vhat,what, qthat, ThetaRhat
  complex, dimension(nxh,ny), intent(inout) :: out1
  complex, dimension(nx,ny, m) :: wrkNL
  complex, dimension(nxph,nyp) :: wrk1
  complex, dimension(nxph,nyp,m) :: wrkNLhat

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

!-------------------------------------------------------------------------------------		
!-------------------------------------------------------------------------------------		
!-------------------------------------------------------------------------------------		
!-------------------------------------------------------------------------------------		

SUBROUTINE RK_FLUX_D(nx, nxh, nxph, ny, nyp, nyph, m, kx, ky, dx, dz, &
& f_star, g_star, epsbar, l, b_star, nu, vrain, uhat, vhat, what, &
& ThetaRhat, qthat, u, ud, v, vd, w, wd, ThetaR, ThetaRd, qt, qtd, fu, &
& fud, fv, fvd, fw, fwd, fThetaR, fThetaRd, fqt, fqtd, rc, mr, ubg, vbg&
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
  REAL, DIMENSION(nx, ny, m), INTENT(IN) :: u, v, w, ThetaR, qt
  REAL, DIMENSION(nx, ny, m), INTENT(IN) :: ud, vd, wd, ThetaRd, qtd
  REAL, DIMENSION(nx, ny), INTENT(INOUT) :: in1
  REAL, DIMENSION(nx, ny), INTENT(INOUT) :: in1d
  REAL, DIMENSION(nx, ny, m) :: qr, qrd
!Complex
  COMPLEX, DIMENSION(nxph, nyp, m), INTENT(IN) :: uhat, vhat, what, &
& qthat, ThetaRhat
COMPLEX, DIMENSION(nxph, nyp, m) :: uhatd, vhatd, whatd, &
& qthatd, ThetaRhatd
  COMPLEX, DIMENSION(nxh, ny), INTENT(INOUT) :: out1
  COMPLEX, DIMENSION(nx, ny, m) :: wrkNL
  COMPLEX, DIMENSION(nx, ny, m) :: wrkNLd
  COMPLEX, DIMENSION(nxph, nyp) :: wrk1
  COMPLEX, DIMENSION(nxph, nyp) :: wrk1d
!Out:
!Complex:
  COMPLEX, DIMENSION(nxph, nyp, m), INTENT(OUT) :: fu, fv, fw, fThetaR, &
& fqt
  COMPLEX, DIMENSION(nxph, nyp, m), INTENT(OUT) :: fud, fvd, fwd, &
& fThetaRd, fqtd
  complex, dimension(nxph,nyp,m) :: wrkNLhat, wrkNLhatd, out1d
  INTRINSIC SUM
  REAL :: temp
  TYPE(C_PTR) :: planf1
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

!---------Coriolis---------------------------------------
  fud = f_star*vhatd
  fu = f_star*vhat
  fvd = -(f_star*uhatd)
  fv = -(f_star*uhat)

  !wrkNLhatd = (0.0,0.0)
  wrkNLhatd = 0
  wrkNLhat = 0
!---------Non-linear advection term: - (uu)_x---------------------------------------
!wrkNL = u*u;
  DO iz=1,m
!wrkNL(:,:,iz); !This is why we call it pseudo-spectral
    in1d = 2*u(:, :, iz)*ud(:, :, iz)    ! du'/dt <- - d/dx(2*u*u')
    in1 = u(:, :, iz)*u(:, :, iz)
!de-aliasing
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fud = fud -(im*kx*wrkNLhatd)               ! du'/dt <- - d/dx(2*u*u')
  fu = fu - im*kx*wrkNLhat
!---------Non-linear advection terms: - (uv)_x and - (uv)_y---------------------------------------	
!wrkNL = u*v;
  DO iz=1,m
!wrkNL(:,:,iz); 
    in1d = v(:, :, iz)*ud(:, :, iz) + u(:, :, iz)*vd(:, :, iz)  
    in1 = u(:, :, iz)*v(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fud = fud - im*ky*wrkNLhatd  ! dv'/dt <- - d/dy(u*v' + u'*v)
  fu = fu - im*ky*wrkNLhat
  fvd = fvd -(im*kx*wrkNLhatd)     ! du'/dt <- - d/dx(u*v' + u'*v)
  fv = fv - im*kx*wrkNLhat
!---------Non-linear advection term: - (vv)_x----------------------------------------	
!wrkNL = v*v;
  DO iz=1,m
!wrkNL(:,:,iz); 
    in1d = 2*v(:, :, iz)*vd(:, :, iz)    ! dv'/dt <- - d/dy(2*v*v')
    in1 = v(:, :, iz)*v(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fvd = fvd - im*ky*wrkNLhatd            ! dv'/dt <- - d/dy(2*v*v')
  fv = fv - im*ky*wrkNLhat
!----------Non-linear advection term: - (wu)_z affecting the u equation --------------------------------------	
!Need to compute u at the w level, due to staggered grid
!Interpolate u:
!*** u on position as w
  ! wrkNLd = (0.0,0.0)
  wrkNLd = 0
  wrkNLd(:, :, 2:m-1) = 0.5*(ud(:, :, 1:m-2)+ud(:, :, 2:m-1))
  wrkNL(:, :, 2:m-1) = 0.5*(u(:, :, 1:m-2)+u(:, :, 2:m-1))
!**** No boundary condition, but w=0 ****
  ! wrkNLd(:, :, 1) = (0.0,0.0)
  ! wrkNL(:, :, 1) = 0
  ! wrkNLd(:, :, m) = (0.0,0.0)
  ! wrkNL(:, :, m) = 0
  wrkNLd(:, :, 1) = 0
  wrkNL(:, :, 1) = 0
  wrkNLd(:, :, m) = 0
  wrkNL(:, :, m) = 0
  wrkNLd = w*wrkNLd + wrkNL*wd      ! du'/dt <- - d/dz(w*u' + w'*u)
  wrkNL = wrkNL*w
  DO iz=1,m
    in1d = wrkNLd(:, :, iz)
    in1 = wrkNL(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fud(:, :, 2:m-2) = fud(:, :, 2:m-2) - (wrkNLhatd(:, :, 3:m-1)-&
&   wrkNLhatd(:, :, 2:m-2))/dz                ! du'/dt <- - d/dz(w*u' + w'*u)
  fu(:, :, 2:m-2) = fu(:, :, 2:m-2) - (wrkNLhat(:, :, 3:m-1)-wrkNLhat(:&
&   , :, 2:m-2))/dz + nu(2)*(uhat(:, :, 3:m-1)+uhat(:, :, 1:m-3))/dz**2
!**** u=0 at the bottom, du/dz=0 on the top ****
!u(0) = -u(1)
  fud(:, :, 1) = fud(:, :, 1) - (wrkNLhatd(:, :, 2)-wrkNLhatd(:, :, 1))/&
&   dz
  fu(:, :, 1) = fu(:, :, 1) - (wrkNLhat(:, :, 2)-wrkNLhat(:, :, 1))/dz +&
&   nu(2)*(uhat(:, :, 2)-uhat(:, :, 1))/dz**2
!u(m)  = u(m-1)
  fud(:, :, m-1) = fud(:, :, m-1) - (wrkNLhatd(:, :, m)-wrkNLhatd(:, :, &
&   m-1))/dz
  fu(:, :, m-1) = fu(:, :, m-1) - (wrkNLhat(:, :, m)-wrkNLhat(:, :, m-1)&
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
  fwd = -(im*kx*wrkNLhatd)       ! dw'/dt <- - d/dx(w*u' + w'*u)
  fw = -(im*kx*wrkNLhat)
!----------Non-linear advection term: - (wv)_z affecting the v equation --------------------------------------	
!Need to compute v at the w level, due to staggered grid
!********* Interpolation for v ********
!*** v on position as w
  wrkNLd(:, :, 2:m-1) = 0.5*(vd(:, :, 1:m-2)+vd(:, :, 2:m-1))
  wrkNL(:, :, 2:m-1) = 0.5*(v(:, :, 1:m-2)+v(:, :, 2:m-1))
!**** No boundary condition, but w=0 ****
  wrkNLd(:, :, 1) = 0
  wrkNL(:, :, 1) = 0
  wrkNLd(:, :, m) = 0
  wrkNL(:, :, m) = 0
  wrkNLd = w*wrkNLd + wrkNL*wd      ! dv'/dt <- - d/dz(w*v' + w'*v)
  wrkNL = wrkNL*w
!----------------------------------
  DO iz=1,m
    in1d = wrkNLd(:, :, iz)
    in1 = wrkNL(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
!-----------------------------------
  fvd(:, :, 2:m-2) = fvd(:, :, 2:m-2) - (wrkNLhatd(:, :, 3:m-1)-&
&   wrkNLhatd(:, :, 2:m-2))/dz      ! dv'/dt <- - d/dz(w*v' + w'*v)
  fv(:, :, 2:m-2) = fv(:, :, 2:m-2) - (wrkNLhat(:, :, 3:m-1)-wrkNLhat(:&
&   , :, 2:m-2))/dz + nu(2)*(vhat(:, :, 3:m-1)+vhat(:, :, 1:m-3))/dz**2
!**** v=0 at the bottom, dv/dz=0 on the top ****
  fvd(:, :, 1) = fvd(:, :, 1) - (wrkNLhatd(:, :, 2)-wrkNLhatd(:, :, 1))/&
&   dz
  fv(:, :, 1) = fv(:, :, 1) - (wrkNLhat(:, :, 2)-wrkNLhat(:, :, 1))/dz +&
&   nu(2)*(vhat(:, :, 2)-vhat(:, :, 1))/dz**2
  fvd(:, :, m-1) = fvd(:, :, m-1) - (wrkNLhatd(:, :, m)-wrkNLhatd(:, :, &
&   m-1))/dz
  fv(:, :, m-1) = fv(:, :, m-1) - (wrkNLhat(:, :, m)-wrkNLhat(:, :, m-1)&
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
  fwd = fwd - im*ky*wrkNLhatd      ! dw'/dt <- - d/dy(w*v' + w'*v)
  fw = fw - im*ky*wrkNLhat
!---------------------------------------------------------------------------------
!Here we need to get qr
! ?? if qt (and theta_r) are perturbation, does the equation need to be linearized?
! ?? if qt = qtd and theta_r = theta_rd
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
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  temp = 0.5*g_star*epsbar
  fwd(:, :, 2:m-1) = fwd(:, :, 2:m-1) + temp*(wrkNLhatd(:, :, 1:m-2)+&
&   wrkNLhatd(:, :, 2:m-1))    !d w /dt <- g_star eps_o q_v'
  fw(:, :, 2:m-1) = fw(:, :, 2:m-1) + temp*(wrkNLhat(:, :, 1:m-2)+&
&   wrkNLhat(:, :, 2:m-1))

!------------q_r '  term affecting the w equation---------------------------------------------------------------------
  DO iz=1,m
    in1d = qrd(:, :, iz)
    in1 = qr(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fwd(:, :, 2:m-1) = fwd(:, :, 2:m-1) - g_star*0.5*(wrkNLhatd(:, :, 1:m-&
&   2)+wrkNLhatd(:, :, 2:m-1))   !d w /dt <- - g_star*q_r
  fw(:, :, 2:m-1) = fw(:, :, 2:m-1) - g_star*0.5*(wrkNLhat(:, :, 1:m-2)+&
&   wrkNLhat(:, :, 2:m-1))

!------------q_r  term affecting the theta_r equation---------------------------------------------------------------------
  fThetaR = 0
  fThetaRd = (0.0,0.0)
  fThetaRd(:, :, 2:m-2) = -(vrain*l*(wrkNLhatd(:, :, 3:m-1)-wrkNLhatd(:&
&   , :, 1:m-3))/(2*dz))  ! d theta_r' /dt <- - V_T L dqr/dz (vrain * l * dqr/dz)
  fThetaR(:, :, 2:m-2) = fThetaR(:, :, 2:m-2) - vrain*l*(wrkNLhat(:, :, &
&   3:m-1)-wrkNLhat(:, :, 1:m-3))/(2*dz)

!**** Neumann boundary condition ****
  fThetaRd(:, :, 1) = fThetaRd(:, :, 1) - vrain*l*(wrkNLhatd(:, :, 2)-&
&   wrkNLhatd(:, :, 1))/(2*dz)
  fThetaR(:, :, 1) = fThetaR(:, :, 1) - vrain*l*(wrkNLhat(:, :, 2)-&
&   wrkNLhat(:, :, 1))/(2*dz)
  fThetaRd(:, :, m-1) = fThetaRd(:, :, m-1) - vrain*l*(wrkNLhatd(:, :, m&
&   -1)-wrkNLhatd(:, :, m-2))/(2*dz)
  fThetaR(:, :, m-1) = fThetaR(:, :, m-1) - vrain*l*(wrkNLhat(:, :, m-1)&
&   -wrkNLhat(:, :, m-2))/(2*dz)

!------------q_r  term affecting the q_t equation---------------------------------------------------------------------
  fqt = 0
  fqtd = (0.0,0.0)
  fqtd(:, :, 2:m-2) = vrain*(wrkNLhatd(:, :, 3:m-1)-wrkNLhatd(:, :, 1:m-&
&   3))/(2*dz)             ! d theta_r' /dt <- - V_T dqr/dz (vrain * dqr/dz)
  fqt(:, :, 2:m-2) = fqt(:, :, 2:m-2) + vrain*(wrkNLhat(:, :, 3:m-1)-&
&   wrkNLhat(:, :, 1:m-3))/(2*dz)

!**** Neumann boundary condition for q_r		
!q_r(0) = q_r(1)
  fqtd(:, :, 1) = fqtd(:, :, 1) + vrain*(wrkNLhatd(:, :, 2)-wrkNLhatd(:&
&   , :, 1))/(2*dz)
  fqt(:, :, 1) = fqt(:, :, 1) + vrain*(wrkNLhat(:, :, 2)-wrkNLhat(:, :, &
&   1))/(2*dz)
!q_r(m) = q_r(m-1)
  fqtd(:, :, m-1) = fqtd(:, :, m-1) + vrain*(wrkNLhatd(:, :, m-1)-&
&   wrkNLhatd(:, :, m-2))/(2*dz)
  fqt(:, :, m-1) = fqt(:, :, m-1) + vrain*(wrkNLhat(:, :, m-1)-wrkNLhat(&
&   :, :, m-2))/(2*dz)

!-------------Theta term affecting the w equation--------------------------------------------------------------------
  DO iz=1,m
    in1d = ThetaRd(:, :, iz) + l*qrd(:, :, iz)
    in1 = ThetaR(:, :, iz) + l*qr(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fwd(:, :, 2:m-1) = fwd(:, :, 2:m-1) + g_star*0.5*(wrkNLhatd(:, :, 1:m-&
&   2)+wrkNLhatd(:, :, 2:m-1))      ! dw'/dt <- g_star (theta' + L*qr')
  fw(:, :, 2:m-1) = fw(:, :, 2:m-1) + g_star*(0.5*(wrkNLhat(:, :, 1:m-2)&
&   +wrkNLhat(:, :, 2:m-1)))

!----------Non-linear advection term: - (ww)_z affecting the w equation --------------------------------------	
!********** Interpolation for w ***************
!*** w on position as u
  wrkNLd(:, :, 1:m-1) = 0.5*(wd(:, :, 1:m-1)+wd(:, :, 2:m))
  wrkNL(:, :, 1:m-1) = 0.5*(w(:, :, 1:m-1)+w(:, :, 2:m))

!**** Zero boundary condition, 0.5*(w(m+1)+w(m-1))=0, w(m)=0.0
  wrkNLd(:, :, m) = -wrkNLd(:, :, m-1)
  wrkNL(:, :, m) = -wrkNL(:, :, m-1)
  wrkNLd = 2*wrkNL*wrkNLd         ! dw'/dt <- d/dz(2w*w')
  wrkNL = wrkNL*wrkNL
  DO iz=1,m
    in1d = wrkNLd(:, :, iz)
    in1 = wrkNL(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fwd(:, :, 2:m-1) = fwd(:, :, 2:m-1) - (wrkNLhatd(:, :, 2:m-1)-&
&   wrkNLhatd(:, :, 1:m-2))/dz         ! dw'/dt <- d/dz(2w*w')
  fw(:, :, 2:m-1) = fw(:, :, 2:m-1) - (wrkNLhat(:, :, 2:m-1)-wrkNLhat(:&
&   , :, 1:m-2))/dz + nu(2)*(what(:, :, 3:m)+what(:, :, 1:m-2))/dz**2
!fw(1,1,2:m-1) = fw(1,1,2:m-1)-nu(2)*(what(1,1,3:m)+what(1,1,1:m-2))/dz**2; !Remove viscosity to VSHF

!**** Zero boundary condition ******	
  !fwd(:, :, 1) = (0.0,0.0)
  fw(:, :, 1) = 0
  fwd(:, :, 1) = 0
  !fwd(:, :, m) = (0.0,0.0)
  fw(:, :, m) = 0
  fwd(:, :, m) = 0

!----------Non-linear advection term: - (u theta_r')_x affecting the theta_r equation --------------------------------------	
!wrkNL = u*ThetaR;
  DO iz=1,m
!wrkNL(:,:,iz); 
    in1d = ThetaR(:, :, iz)*ud(:, :, iz) + u(:, :, iz)*ThetaRd(:, :, iz)
    in1 = u(:, :, iz)*ThetaR(:, :, iz)  !dtheta_r'/dt <- -d/dx(u*theta_r' + u'*theta_r)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fThetaRd = fThetaRd - im*kx*wrkNLhatd  !dtheta_r'/dt <- -d/dx(u*theta_r' + u'*theta_r)
  fThetaR = fThetaR - im*kx*wrkNLhat

!----------Non-linear advection term: - (v theta_r')_y affecting the theta_r equation --------------------------------------	
!wrkNL = v*ThetaR;
  DO iz=1,m
!wrkNL(:,:,iz); 
    in1d = ThetaR(:, :, iz)*vd(:, :, iz) + v(:, :, iz)*ThetaRd(:, :, iz)
    in1 = v(:, :, iz)*ThetaR(:, :, iz)   !dtheta_r'/dt <- -d/dy(v*theta_r' + v'*theta_r)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fThetaRd = fThetaRd - im*ky*wrkNLhatd     !dtheta_r'/dt <- -d/dy(v*theta_r' + v'*theta_r)
  fThetaR = fThetaR - im*ky*wrkNLhat
!---------------------------------------------------------------------------------
!********************* For Fixed Background
  DO iz=1,m-2
    fThetaRd(:, :, iz) = fThetaRd(:, :, iz) - b_star*0.5*(whatd(:, :, iz&
&     )+whatd(:, :, iz+1))
    fThetaR(:, :, iz) = fThetaR(:, :, iz) - 0.5*(what(:, :, iz)+what(:, &
&     :, iz+1))*b_star
  END DO
!----------Non-linear advection term: - (w theta_r')_z affecting the theta_r equation --------------------------------------	
!********** Interpolation for ThetaR **************
!*** ThetaR on position as w
  wrkNLd(:, :, 2:m-1) = 0.5*(ThetaRd(:, :, 1:m-2)+ThetaRd(:, :, 2:m-1))
  wrkNL(:, :, 2:m-1) = 0.5*(ThetaR(:, :, 1:m-2)+ThetaR(:, :, 2:m-1))
!****   Neumann boundary condition ****
  wrkNLd(:, :, 1) = ThetaRd(:, :, 1)
  wrkNL(:, :, 1) = ThetaR(:, :, 1)
  wrkNLd(:, :, m) = ThetaRd(:, :, m-1)
  wrkNL(:, :, m) = ThetaR(:, :, m-1)
  wrkNLd = wrkNL*wd + w*wrkNLd              !dThetaR'/dt <- -d/dz(w*theta_r' + w'*theta_r)
  wrkNL = w*wrkNL
!----------------------------------
  DO iz=1,m
    in1d = wrkNLd(:, :, iz)
    in1 = wrkNL(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
!**************************************	
  fThetaRd(:, :, 2:m-2) = fThetaRd(:, :, 2:m-2) - (wrkNLhatd(:, :, 3:m-1&
&   )-wrkNLhatd(:, :, 2:m-2))/dz               !dThetaR'/dt <- -d/dz(w*theta_r' + w'*theta_r)
  fThetaR(:, :, 2:m-2) = fThetaR(:, :, 2:m-2) - (wrkNLhat(:, :, 3:m-1)-&
&   wrkNLhat(:, :, 2:m-2))/dz + nu(2)*(ThetaRhat(:, :, 3:m-1)+ThetaRhat(&
&   :, :, 1:m-3))/dz**2

!**** Neumann boundary condition ****
  fThetaRd(:, :, 1) = fThetaRd(:, :, 1) - (wrkNLhatd(:, :, 2)-wrkNLhatd(&
&   :, :, 1))/dz
  fThetaR(:, :, 1) = fThetaR(:, :, 1) - (wrkNLhat(:, :, 2)-wrkNLhat(:, :&
&   , 1))/dz + nu(2)*(ThetaRhat(:, :, 2)+ThetaRhat(:, :, 1))/dz**2
  fThetaRd(:, :, m-1) = fThetaRd(:, :, m-1) - (wrkNLhatd(:, :, m)-&
&   wrkNLhatd(:, :, m-1))/dz
  fThetaR(:, :, m-1) = fThetaR(:, :, m-1) - (wrkNLhat(:, :, m)-wrkNLhat(&
&   :, :, m-1))/dz + nu(2)*(ThetaRhat(:, :, m-1)+ThetaRhat(:, :, m-2))/&
&   dz**2
!fThetaR(1,1,2:m-2) = fThetaR(1,1,2:m-2)-nu(2)*(ThetaRhat(1,1,3:m-1)+ThetaRhat(1,1,1:m-3))/dz**2; !Remove viscosity to VSHF
!fThetaR(1,1,1) = fThetaR(1,1,1)-nu(2)*(ThetaRhat(1,1,2)+ThetaRhat(1,1,1))/dz**2; !Remove viscosity to VSHF
!fThetaR(1,1,m-1) = fThetaR(1,1,m-1)-nu(2)*(ThetaRhat(1,1,m-1)+ThetaRhat(1,1,m-2))/dz**2; !Remove viscosity to VSHF
!----------------------------------
!Radiative cooling
  DO iz=1,m
!Only one eigenmode
    fThetaR(1, 1, iz) = fThetaR(1, 1, iz) + rc(iz)
  END DO
!----------Non-linear advection term: - (u q_t')_x affecting the q_t equation --------------------------------------	
  DO iz=1,m
!wrkNL(:,:,iz); 
    in1d = qt(:, :, iz)*ud(:, :, iz) + u(:, :, iz)*qtd(:, :, iz)
    in1 = u(:, :, iz)*qt(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fqtd = fqtd - im*kx*wrkNLhatd     !dqt'/dt <- -d/dx(u*qt' + u'*qt)
  fqt = fqt - im*kx*wrkNLhat
!----------Non-linear advection term: - (v q_t')_y affecting the q_t equation --------------------------------------	
  DO iz=1,m
!wrkNL(:,:,iz); 
    in1d = qt(:, :, iz)*vd(:, :, iz) + v(:, :, iz)*qtd(:, :, iz)
    in1 = v(:, :, iz)*qt(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fqtd = fqtd - im*ky*wrkNLhatd     !dqt'/dt <- -d/dy(v*qt' + v'*qt)
  fqt = fqt - im*ky*wrkNLhat

!----------Non-linear advection term: - (w q_t')_z affecting the q_t equation --------------------------------------	
!********** Interpolation for qv **************
!*** qv on position as w
  wrkNLd(:, :, 2:m-1) = 0.5*(qtd(:, :, 1:m-2)+qtd(:, :, 2:m-1))
  wrkNL(:, :, 2:m-1) = 0.5*(qt(:, :, 1:m-2)+qt(:, :, 2:m-1))
!****"zero on bottom" Neumann boundary condition **** 
  wrkNLd(:, :, 1) = qtd(:, :, 1)
  wrkNL(:, :, 1) = qt(:, :, 1)
  wrkNLd(:, :, m) = qtd(:, :, m-1)
  wrkNL(:, :, m) = qt(:, :, m-1)
  wrkNLd = wrkNL*wd + w*wrkNLd     !dqt'/dt <- -d/dz(w*qt' + w'*qt)
  wrkNL = w*wrkNL
  DO iz=1,m
    in1d = wrkNLd(:, :, iz)
    in1 = wrkNL(:, :, iz)
    ! wrk1d = in1d
    ! wrk1 = in1
    call dfftw_execute_dft_r2c(planf1,in1,out1)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
    call dfftw_execute_dft_r2c(planf1,in1d,out1d)
    call unpadm(nxh,nxph,ny,nyp,nyph,out1d,wrk1d) 
    wrkNLhatd(:, :, iz) = wrk1d/ns
    wrkNLhat(:, :, iz) = wrk1/ns
  END DO
  fqtd(:, :, 2:m-2) = fqtd(:, :, 2:m-2) - (wrkNLhatd(:, :, 3:m-1)-&
&   wrkNLhatd(:, :, 2:m-2))/dz     !dqt'/dt <- -d/dz(w*qt' + w'*qt)
  fqt(:, :, 2:m-2) = fqt(:, :, 2:m-2) - (wrkNLhat(:, :, 3:m-1)-wrkNLhat(&
&   :, :, 2:m-2))/dz + nu(2)*(qthat(:, :, 3:m-1)+qthat(:, :, 1:m-3))/dz&
&   **2

!****"Non-zero on bottom" Neumann boundary condition **** 		
  fqtd(:, :, 1) = fqtd(:, :, 1) - (wrkNLhatd(:, :, 2)-wrkNLhatd(:, :, 1)&
&   )/dz
  fqt(:, :, 1) = fqt(:, :, 1) - (wrkNLhat(:, :, 2)-wrkNLhat(:, :, 1))/dz&
&   + nu(2)*(qthat(:, :, 2)+qthat(:, :, 1))/dz**2
  fqtd(:, :, m-1) = fqtd(:, :, m-1) - (wrkNLhatd(:, :, m)-wrkNLhatd(:, :&
&   , m-1))/dz
  fqt(:, :, m-1) = fqt(:, :, m-1) - (wrkNLhat(:, :, m)-wrkNLhat(:, :, m-&
&   1))/dz + nu(2)*(qthat(:, :, m-1)+qthat(:, :, m-2))/dz**2
!fqt(1,1,2:m-2) = fqt(1,1,2:m-2)-nu(2)*(qthat(1,1,3:m-1)+qthat(1,1,1:m-3))/dz**2; !Remove viscosity to VSHF
!fqt(1,1,1) = fqt(1,1,1)-nu(2)*(qthat(1,1,2)+qthat(1,1,1))/dz**2; !Remove viscosity to VSHF
!fqt(1,1,m-1) = fqt(1,1,m-1)-nu(2)*(qthat(1,1,m-1)+qthat(1,1,m-2))/dz**2; !Remove viscosity to VSHF

!----------Background term -w d_qvditle(z)/dz--------------------------------------	
!********************* For Fixed Background
  DO iz=1,m-2
    fqtd(:, :, iz) = fqtd(:, :, iz) - 0.5*(fqvdz(iz)*(whatd(:, :, iz)+&
&     whatd(:, :, iz+1)))
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

!-------------------------------------------------------------------------------------		
!-------------------------------------------------------------------------------------		
!-------------------------------------------------------------------------------------		
!-------------------------------------------------------------------------------------	

SUBROUTINE RK_FLUX_B(nx, nxh, nxph, ny, nyp, nyph, m, kx, ky, dx, dz, &
& f_star, g_star, epsbar, l, b_star, nu, vrain, uhat, vhat, what, &
& ThetaRhat, qthat, u, ub, v, vb, w, wb, ThetaR, ThetaRb, qt, qtb, fu, &
& fub, fv, fvb, fw, fwb, fThetaR, fThetaRb, fqt, fqtb, rc, mr, ubg, vbg&
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
  REAL, DIMENSION(nx, ny, m), INTENT(IN) :: u, v, w, ThetaR, qt
  REAL, DIMENSION(nx, ny, m) :: ub, vb, wb, ThetaRb, qtb, qrb
  REAL, DIMENSION(nx, ny, m) :: qr
  REAL, DIMENSION(nx, ny), INTENT(INOUT) :: in1
  REAL, DIMENSION(nx, ny), INTENT(INOUT) :: in1b
!Complex
  COMPLEX, DIMENSION(nxph, nyp, m), INTENT(IN) :: uhat, vhat, what, &
& qthat, ThetaRhat
  COMPLEX, DIMENSION(nxph, nyp, m) :: wrkNLhat, wrkNLhatb
  COMPLEX, DIMENSION(nxh, ny), INTENT(INOUT) :: out1
  COMPLEX, DIMENSION(nx, ny, m) :: wrkNL
  COMPLEX, DIMENSION(nx, ny, m) :: wrkNLb
  COMPLEX, DIMENSION(nxph, nyp) :: wrk1
  COMPLEX, DIMENSION(nxph, nyp) :: wrk1b, out1b
!Out:
!Complex:
  COMPLEX, DIMENSION(nxph, nyp, m) :: fu, fv, fw, fThetaR, fqt
  COMPLEX, DIMENSION(nxph, nyp, m) :: fub, fvb, fwb, fThetaRb, fqtb
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
  wrkNL(:, :, 2:m-1) = 0.5*(u(:, :, 1:m-2)+u(:, :, 2:m-1))
!**** No boundary condition, but w=0 ****
  wrkNL(:, :, 1) = 0
  wrkNL(:, :, m) = 0
  CALL PUSHCOMPLEX8ARRAY(wrkNL, nx*ny*m)
  wrkNL = wrkNL*w
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
  wrkNL(:, :, 2:m-1) = 0.5*(v(:, :, 1:m-2)+v(:, :, 2:m-1))
!**** No boundary condition, but w=0 ****
  wrkNL(:, :, 1) = 0
  wrkNL(:, :, m) = 0
  CALL PUSHCOMPLEX8ARRAY(wrkNL, nx*ny*m)
  wrkNL = wrkNL*w
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
  wrkNL(:, :, 1:m-1) = 0.5*(w(:, :, 1:m-1)+w(:, :, 2:m))
!**** Zero boundary condition, 0.5*(w(m+1)+w(m-1))=0, w(m)=0.0
  wrkNL(:, :, m) = -wrkNL(:, :, m-1)
  CALL PUSHCOMPLEX8ARRAY(wrkNL, nx*ny*m)
  wrkNL = wrkNL*wrkNL
!----------Non-linear advection term: - (w theta_r')_z affecting the theta_r equation --------------------------------------	
!********** Interpolation for ThetaR **************
!*** ThetaR on position as w
  wrkNL(:, :, 2:m-1) = 0.5*(ThetaR(:, :, 1:m-2)+ThetaR(:, :, 2:m-1))
!****   Neumann boundary condition ****
  wrkNL(:, :, 1) = ThetaR(:, :, 1)
  wrkNL(:, :, m) = ThetaR(:, :, m-1)
  CALL PUSHCOMPLEX8ARRAY(wrkNL, nx*ny*m)
  wrkNL = w*wrkNL
!----------Non-linear advection term: - (w q_t')_z affecting the q_t equation --------------------------------------	
!********** Interpolation for qv **************
!*** qv on position as w
  wrkNL(:, :, 2:m-1) = 0.5*(qt(:, :, 1:m-2)+qt(:, :, 2:m-1))
!****"zero on bottom" Neumann boundary condition **** 
  wrkNL(:, :, 1) = qt(:, :, 1)
  wrkNL(:, :, m) = qt(:, :, m-1)
  DO iz=m-2,1,-1
    whatb(:, :, iz) = whatb(:, :, iz) - fqvdz(iz)*0.5*fqtb(:, :, iz)
    whatb(:, :, iz+1) = whatb(:, :, iz+1) - fqvdz(iz)*0.5*fqtb(:, :, iz)
!     fqvdzb(iz) = fqvdzb(iz) - SUM((what(:, :, iz)+what(:, :, iz+1))*fqtb&
! &     (:, :, iz))*0.5
  END DO
  ! wrkNLhatb = (0.0,0.0)
  wrkNLhatb(:, :, m) = wrkNLhatb(:, :, m) - fqtb(:, :, m-1)/dz
  wrkNLhatb(:, :, m-1) = wrkNLhatb(:, :, m-1) + fqtb(:, :, m-1)/dz
  wrkNLhatb(:, :, 2) = wrkNLhatb(:, :, 2) - fqtb(:, :, 1)/dz
  wrkNLhatb(:, :, 1) = wrkNLhatb(:, :, 1) + fqtb(:, :, 1)/dz
  wrkNLhatb(:, :, 3:m-1) = wrkNLhatb(:, :, 3:m-1) - fqtb(:, :, 2:m-2)/dz
  wrkNLhatb(:, :, 2:m-2) = wrkNLhatb(:, :, 2:m-2) + fqtb(:, :, 2:m-2)/dz
  ! wrkNLb = (0.0,0.0)
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    wrkNLb(:,:,iz) = in1b
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = 0
    ! in1b = in1b + wrk1b
    ! wrkNLb(:, :, iz) = wrkNLb(:, :, iz) + in1b
  END DO
  wb = wb + wrkNL*wrkNLb
      
  !----------Non-linear advection term: - (w q_t')_z affecting the q_t equation --------------------------------------      
  !dqt'/dt <- -d/dz(w*qt' + w'*qt)
  ! a_w = d/dz(qt*a_fqt)
  ! a_qt = d/dz(w*a_fqt)
  wrkNLb = w*wrkNLb
  qtb = 0.0
  qtb(:, :, m-1) = qtb(:, :, m-1) + wrkNLb(:, :, m)
  wrkNLb(:, :, m) = (0.0,0.0)
  qtb(:, :, 1) = qtb(:, :, 1) + REAL(wrkNLb(:, :, 1))
  wrkNLb(:, :, 1) = (0.0,0.0)
  qtb(:, :, 1:m-2) = qtb(:, :, 1:m-2) + 0.5*wrkNLb(:, :, 2:m-1)
  qtb(:, :, 2:m-1) = qtb(:, :, 2:m-1) + 0.5*wrkNLb(:, :, 2:m-1)
  wrkNLb(:, :, 2:m-1) = (0.0,0.0)
  wrkNLhatb = wrkNLhatb - im*ky*fqtb

  !----------Non-linear advection term: - (v q_t')_y affecting the q_t equation --------------------------------------      
  !dqt'/dt <- -d/dy(v*qt' + v'*qt)
  ! a_v = d/dy(qt*a_v)
  ! a_qt = -d/dy(v*a_qt)
  vb = 0.0
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = 0
    ! in1b = in1b + wrk1b
    vb(:, :, iz) = vb(:, :, iz) + qt(:, :, iz)*in1b
    qtb(:, :, iz) = qtb(:, :, iz) + v(:, :, iz)*in1b
    in1b = 0.0
  END DO
  wrkNLhatb = wrkNLhatb - im*kx*fqtb
  
  !----------Non-linear advection term: - (u q_t')_x affecting the q_t equation --------------------------------------      
  !dqt'/dt <- -d/dx(u*qt' + u'*qt)
  ! a_u = d/dx(qt*a_u)
  ! a_qt = d/dx(u*a_qt)
  ub = 0.0
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = in1b + wrk1b
    ub(:, :, iz) = ub(:, :, iz) + qt(:, :, iz)*in1b
    qtb(:, :, iz) = qtb(:, :, iz) + u(:, :, iz)*in1b
    in1b = 0.0
  END DO
  
  !----------Non-linear advection term: - (w theta_r')_z affecting the theta_r equation --------------------------------------
  !dThetaR'/dt <- -d/dz(w*theta_r' + w'*theta_r)
  ! a_w = d/dz(ThetaR*a_fThetaR) ? 
  ! a_ThetaR = d/dz(w*a_fThetaR) ?
  wrkNLhatb(:, :, m) = wrkNLhatb(:, :, m) - fThetaRb(:, :, m-1)/dz
  wrkNLhatb(:, :, m-1) = wrkNLhatb(:, :, m-1) + fThetaRb(:, :, m-1)/dz
  wrkNLhatb(:, :, 2) = wrkNLhatb(:, :, 2) - fThetaRb(:, :, 1)/dz
  wrkNLhatb(:, :, 1) = wrkNLhatb(:, :, 1) + fThetaRb(:, :, 1)/dz
  wrkNLhatb(:, :, 3:m-1) = wrkNLhatb(:, :, 3:m-1) - fThetaRb(:, :, 2:m-2&
&   )/dz
  wrkNLhatb(:, :, 2:m-2) = wrkNLhatb(:, :, 2:m-2) + fThetaRb(:, :, 2:m-2&
&   )/dz
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    wrkNLb(:,:,iz) = in1b
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = 0
    ! in1b = in1b + wrk1b
    ! wrkNLb(:, :, iz) = wrkNLb(:, :, iz) + in1b
  END DO
  CALL POPCOMPLEX8ARRAY(wrkNL, nx*ny*m)
  wb = wb + wrkNL*wrkNLb
  wrkNLb = w*wrkNLb
  ThetaRb = 0.0
  ThetaRb(:, :, m-1) = ThetaRb(:, :, m-1) + wrkNLb(:, :, m)
  wrkNLb(:, :, m) = (0.0,0.0)
  ThetaRb(:, :, 1) = ThetaRb(:, :, 1) + wrkNLb(:, :, 1)
  wrkNLb(:, :, 1) = (0.0,0.0)
  ThetaRb(:, :, 1:m-2) = ThetaRb(:, :, 1:m-2) + 0.5*wrkNLb(:, :, 2:m-1)
  ThetaRb(:, :, 2:m-1) = ThetaRb(:, :, 2:m-1) + 0.5*wrkNLb(:, :, 2:m-1)
  wrkNLb(:, :, 2:m-1) = (0.0,0.0)
  DO iz=m-2,1,-1
    whatb(:, :, iz) = whatb(:, :, iz) - b_star*0.5*fThetaRb(:, :, iz)
    whatb(:, :, iz+1) = whatb(:, :, iz+1)- b_star*0.5*fThetaRb(:, :, iz)
  END DO
  
  !----------Non-linear advection term: - (v theta_r')_y affecting the theta_r equation --------------------------------------
  !dtheta_r'/dt <- -d/dy(v*theta_r' + v'*theta_r)
  ! a_v = d/dy(ThetaR * a_fThetaR)
  ! a_ThetaR = d/dy(v * a_fThetaR)
  wrkNLhatb = wrkNLhatb - im*ky*fThetaRb
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    wrkNLb(:,:,iz) = in1b
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = 0
    ! in1b = in1b + wrk1b
    ! wrkNLb(:, :, iz) = wrkNLb(:, :, iz) + in1b
    vb(:, :, iz) = vb(:, :, iz) + ThetaR(:, :, iz)*in1b
    ThetaRb(:, :, iz) = ThetaRb(:, :, iz) + v(:, :, iz)*in1b
  END DO

  !----------Non-linear advection term: - (u theta_r')_x affecting the theta_r equation --------------------------------------
  !dtheta_r'/dt <- -d/dx(u*theta_r' + u'*theta_r)
  ! a_u = d/dx(ThetaR * a_fThetaR)
  ! a_ThetaR = d/dx(u * a_fThetaR)
  wrkNLhatb = wrkNLhatb - im*kx*fThetaRb
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = in1b + wrk1b
    ub(:, :, iz) = ub(:, :, iz) + ThetaR(:, :, iz)*in1b
    ThetaRb(:, :, iz) = ThetaRb(:, :, iz) + u(:, :, iz)*in1b
  END DO

  !----------Non-linear advection term: - (ww)_z affecting the w equation --------------------------------------
  ! dw'/dt <- -d/dz(2w*w')
  ! a_w = d/dz (2*w * a_fw)
  fwb(:, :, m) = 0
  fwb(:, :, 1) = 0
  wrkNLhatb(:, :, 2:m-1) = wrkNLhatb(:, :, 2:m-1) - fwb(:, :, 2:m-1)/dz
  wrkNLhatb(:, :, 1:m-2) = wrkNLhatb(:, :, 1:m-2) + fwb(:, :, 2:m-1)/dz
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    wrkNLb(:,:,iz) = in1b
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = 0
    ! in1b = in1b + wrk1b
    ! wrkNLb(:, :, iz) = wrkNLb(:, :, iz) + in1b
  END DO

  !-------------Theta term affecting the w equation-------------------------------------------------
  ! dw'/dt <- g_star (theta' + L*qr')  ?
  ! a_w = 
  CALL POPCOMPLEX8ARRAY(wrkNL, nx*ny*m)
  wrkNLb = 2*wrkNL*wrkNLb
  wrkNLb(:, :, m-1) = wrkNLb(:, :, m-1) - wrkNLb(:, :, m)
  wrkNLb(:, :, m) = (0.0,0.0)
  wb(:, :, 1:m-1) = wb(:, :, 1:m-1) + REAL(0.5*wrkNLb(:, :, 1:m-1))
  wb(:, :, 2:m) = wb(:, :, 2:m) + REAL(0.5*wrkNLb(:, :, 1:m-1))
  wrkNLb(:, :, 1:m-1) = (0.0,0.0)
  ! a_ThetaR = g_star * 0.5? * a_fw
  wrkNLhatb(:, :, 1:m-2) = wrkNLhatb(:, :, 1:m-2) + g_star*0.5*fwb(:, :&
&   , 2:m-1)
  wrkNLhatb(:, :, 2:m-1) = wrkNLhatb(:, :, 2:m-1) + g_star*0.5*fwb(:, :&
&   , 2:m-1)

!------------q_r  terms -------------------------------------------------
! a_qr from dw/dt
! a_qr = g_star * a_fw  ?
! a_theta = g_star * a_fw
  qrb = 0.0
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = in1b + wrk1b
    ThetaRb(:, :, iz) = ThetaRb(:, :, iz) + in1b
    qrb(:, :, iz) = qrb(:, :, iz) + l*in1b
  END DO    
  !------------q_r  term affecting the q_t equation-------------------------------------------------
  ! dqt' /dt <- - V_T*dqr/dz 
  ! a_qr <- V_T* d/dz (a_fqt) [a_qr = wrkNLhatb]
  tempb0 = vrain*fqtb(:, :, m-1)/(2*dz)
  wrkNLhatb(:, :, m-1) = wrkNLhatb(:, :, m-1) + tempb0
  wrkNLhatb(:, :, m-2) = wrkNLhatb(:, :, m-2) - tempb0
  tempb0 = vrain*fqtb(:, :, 1)/(2*dz)
  wrkNLhatb(:, :, 2) = wrkNLhatb(:, :, 2) + tempb0
  wrkNLhatb(:, :, 1) = wrkNLhatb(:, :, 1) - tempb0
  wrkNLhatb(:, :, 3:m-1) = wrkNLhatb(:, :, 3:m-1) + vrain*fqtb(:, :, 2:m&
&   -2)/(2*dz)
  wrkNLhatb(:, :, 1:m-3) = wrkNLhatb(:, :, 1:m-3) - vrain*fqtb(:, :, 2:m&
&   -2)/(2*dz)
  !------------q_r  term affecting the theta_r equation---------------------------------------------------------------------
  ! d theta_r' /dt <- - V_T L dqr/dz 
  ! a_qr <- V_T * L * d/dz(a_fThetaR)  [a_qr = wrkNLhatb]
  tempb0 = -(vrain*l*fThetaRb(:, :, m-1)/(2*dz))
  wrkNLhatb(:, :, m-1) = wrkNLhatb(:, :, m-1) + tempb0
  wrkNLhatb(:, :, m-2) = wrkNLhatb(:, :, m-2) - tempb0
  tempb0 = -(vrain*l*fThetaRb(:, :, 1)/(2*dz))
  wrkNLhatb(:, :, 2) = wrkNLhatb(:, :, 2) + tempb0
  wrkNLhatb(:, :, 1) = wrkNLhatb(:, :, 1) - tempb0
  tempb = -(vrain*l*fThetaRb(:, :, 2:m-2)/(2*dz))
  wrkNLhatb(:, :, 3:m-1) = wrkNLhatb(:, :, 3:m-1) + tempb
  wrkNLhatb(:, :, 1:m-3) = wrkNLhatb(:, :, 1:m-3) - tempb
  !------------q_r '  term affecting the w equation---------------------------------------------------------------------
  !d w /dt <- - g_star*q_r
  ! a_qr <- - g_star* a_fw  [a_qr = wrkNLhatb]
  wrkNLhatb(:, :, 1:m-2) = wrkNLhatb(:, :, 1:m-2) - g_star*0.5*fwb(:, :&
&   , 2:m-1)
  wrkNLhatb(:, :, 2:m-1) = wrkNLhatb(:, :, 2:m-1) - g_star*0.5*fwb(:, :&
&   , 2:m-1)
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = in1b + wrk1b
    qrb(:, :, iz) = qrb(:, :, iz) + in1b
    ! update a_qr = a_qr + wrkNLhatb
  END DO

  !------------q_v '  term affecting the w equation--------------------------------------------
  !d w /dt <- g_star eps_o q_v'
  temp = 0.5*g_star*epsbar
  wrkNLhatb(:, :, 1:m-2) = wrkNLhatb(:, :, 1:m-2) + temp*fwb(:, :, 2:m-1&
&   )
  wrkNLhatb(:, :, 2:m-1) = wrkNLhatb(:, :, 2:m-1) + temp*fwb(:, :, 2:m-1&
&   )
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = in1b + wrk1b
    qvbarb = -SUM(in1b)
    qtb(:, :, iz) = qtb(:, :, iz) + in1b + qvbarb/ns
    qrb(:, :, iz) = qrb(:, :, iz) - in1b - qvbarb/ns
  END DO

  !Here we need to get qr
  DO iz=m,1,-1
    mask(:, :) = qvini(iz) + qt(:, :, iz) .GT. qvs(iz)
    WHERE (mask(:, :)) 
      qtb(:, :, iz) = qtb(:, :, iz) + qrb(:, :, iz)
      qrb(:, :, iz) = 0.0
    END WHERE
  END DO
  !----------Non-linear advection term: - (wv)_y affecting the w equation ---------------------
  ! dw'/dt <- - d/dy(w*v' + w'*v)
  ! a_v = d/dy(w*a_fw)
  ! a_w = d/dy(v*a_fw)
  !----------Non-linear advection term: - (wv)_z affecting the v equation ---------------------
  ! dv'/dt <- - d/dz(w*v' + w'*v)
  ! a_v = d/dz(w*a_fv)
  ! a_w = d/dz(u*a_fv)
  ! a_w = d/dy(v*a_fw) + d/dz(v*a_fv)
  ! a_v = d/dy(w*a_fw) + d/dz(w*a_fv)
  wrkNLhatb = wrkNLhatb - im*ky*fwb
  wrkNLhatb(:, :, m) = wrkNLhatb(:, :, m) - fvb(:, :, m-1)/dz
  wrkNLhatb(:, :, m-1) = wrkNLhatb(:, :, m-1) + fvb(:, :, m-1)/dz
  wrkNLhatb(:, :, 2) = wrkNLhatb(:, :, 2) - fvb(:, :, 1)/dz
  wrkNLhatb(:, :, 1) = wrkNLhatb(:, :, 1) + fvb(:, :, 1)/dz
  wrkNLhatb(:, :, 3:m-1) = wrkNLhatb(:, :, 3:m-1) - fvb(:, :, 2:m-2)/dz
  wrkNLhatb(:, :, 2:m-2) = wrkNLhatb(:, :, 2:m-2) + fvb(:, :, 2:m-2)/dz
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    wrkNLb(:,:,iz) = in1b
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = 0
    ! in1b = in1b + wrk1b
    ! wrkNLb(:, :, iz) = wrkNLb(:, :, iz) + in1b
  END DO
  CALL POPCOMPLEX8ARRAY(wrkNL, nx*ny*m)
  wb = wb + REAL(CONJG(wrkNL)*wrkNLb)
  wrkNLb = w*wrkNLb
  wrkNLb(:, :, m) = (0.0,0.0)
  wrkNLb(:, :, 1) = (0.0,0.0)
  vb(:, :, 1:m-2) = vb(:, :, 1:m-2) + 0.5*wrkNLb(:, :, 2:m-1)
  vb(:, :, 2:m-1) = vb(:, :, 2:m-1) + 0.5*wrkNLb(:, :, 2:m-1)
  wrkNLb(:, :, 2:m-1) = (0.0,0.0)
  !----------Non-linear advection term: - (wu)_x affecting the w equation ---------------------
  ! dw'/dt <- - d/dx(w*u' + w'*u)
  ! a_u = d/dx(w*a_fw)
  ! a_w = d/dx(u*a_fw)
  !----------Non-linear advection term: - (wu)_z affecting the u equation -------------------
  ! du'/dt <- - d/dz(w*u' + w'*u)
  ! a_u = d/dz(w*a_fu)
  ! a_w = d/dz(u*a_fu)
  ! a_w = d/dx(u*a_fw) + d/dz(u*a_fu)
  ! a_u = d/dx(w*a_fw) + d/dz(w*a_fu)
  wrkNLhatb = wrkNLhatb - im*kx*fwb
  wrkNLhatb(:, :, m) = wrkNLhatb(:, :, m) - fub(:, :, m-1)/dz
  wrkNLhatb(:, :, m-1) = wrkNLhatb(:, :, m-1) + fub(:, :, m-1)/dz
  wrkNLhatb(:, :, 2) = wrkNLhatb(:, :, 2) - fub(:, :, 1)/dz
  wrkNLhatb(:, :, 1) = wrkNLhatb(:, :, 1) + fub(:, :, 1)/dz
  wrkNLhatb(:, :, 3:m-1) = wrkNLhatb(:, :, 3:m-1) - fub(:, :, 2:m-2)/dz
  wrkNLhatb(:, :, 2:m-2) = wrkNLhatb(:, :, 2:m-2) + fub(:, :, 2:m-2)/dz
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    wrkNLb = in1b(:,:,iz)
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = 0
    ! in1b = in1b + wrk1b
    ! wrkNLb(:, :, iz) = wrkNLb(:, :, iz) + in1b
  END DO
  CALL POPCOMPLEX8ARRAY(wrkNL, nx*ny*m)
  wb = wb + wrkNL*wrkNLb
  wrkNLb = w*wrkNLb
  wrkNLb(:, :, m) = (0.0,0.0)
  wrkNLb(:, :, 1) = (0.0,0.0)

  ! a_u = w* (d/dx(*a_fw) + d/dz(*a_fu))
  ub(:, :, 1:m-2) = ub(:, :, 1:m-2) + 0.5*wrkNLb(:, :, 2:m-1)
  ub(:, :, 2:m-1) = ub(:, :, 2:m-1) + 0.5*wrkNLb(:, :, 2:m-1)
  
  ! dv'/dt <- - d/dy(2*v*v')
  ! a_v = d/dy (2*v*a_fv)
  wrkNLhatb = wrkNLhatb - im*ky*fvb
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = in1b + wrk1b
    vb(:, :, iz) = vb(:, :, iz) + 2*v(:, :, iz)*in1b
  END DO

  ! dv'/dt <- - d/dy(u*v' + u'*v)
  ! du'/dt <- - d/dx(u*v' + u'*v)
  ! a_u = d/dy(v* a_fv) + d/dx(v*a_fu)
  ! a_v = d/dy(u* a_fv) + d/dx(u*a_fu)
  wrkNLhatb = wrkNLhatb - im*kx*fvb - im*ky*fub
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = in1b + wrk1b
    ub(:, :, iz) = ub(:, :, iz) + v(:, :, iz)*in1b
    vb(:, :, iz) = vb(:, :, iz) + u(:, :, iz)*in1b
  END DO

  ! du'/dt <- - d/dx(2*u*u')
  ! a_u = d/dx (2*u*a_fu)
  wrkNLhatb = wrkNLhatb - im*kx*fub
  DO iz=m,1,-1
    wrk1b = wrkNLhatb(:,:,iz)/ns;
    call dfftw_execute_dft_r2c(planf1,wrk1b,out1b);
    call unpadm(nxh,nxph,ny,nyp,nyph,out1b,in1b) 
    ! wrk1b = (0.0,0.0)
    ! wrk1b = wrkNLhatb(:, :, iz)/ns
    ! wrkNLhatb(:, :, iz) = (0.0,0.0)
    ! in1b = in1b + wrk1b
    ub(:, :, iz) = ub(:, :, iz) + 2*u(:, :, iz)*in1b
  END DO
  uhatb = uhatb - f_star*fvb
  vhatb = vhatb + f_star*fub
!   fThetaRb = (0.0,0.0)
!   fqtb = (0.0,0.0)
!   fub = (0.0,0.0)
!   fvb = (0.0,0.0)
!   fwb = (0.0,0.0)
END SUBROUTINE RK_FLUX_B