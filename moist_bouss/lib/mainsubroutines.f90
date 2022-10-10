
!-----------------------------MODULE--------------------------

module varsRK
	use, intrinsic :: iso_c_binding
	!Allocatable:
	!Real:
	real(C_DOUBLE), dimension(:,:,:),allocatable :: qr, wrkNL
	!Complex:
	complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: wrk1
	complex(C_DOUBLE_COMPLEX), dimension(:,:,:),allocatable :: wrkNLhat
end module varsRK

!-----------------------------MODULE--------------------------

module varsThomas
	use, intrinsic :: iso_c_binding
	!Allocatable
	complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: fac
	complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable :: cwrk
end module varsThomas

!-----------------------------MODULE--------------------------

!mpads: De-aliasing
module mpads
contains
	!-----subroutine-------
	pure subroutine unpadm(nxh,nxph,ny,nyp,nyph,var,varp)
	use, intrinsic :: iso_c_binding
	implicit none
		!In:
		!Integer
		integer(C_INT), intent(in) :: nxh,nxph,ny,nyp,nyph
		!Complex
		complex(C_DOUBLE_COMPLEX), dimension(nxh,ny), intent(in) :: var
		!Out:
		!Complex
		complex(C_DOUBLE_COMPLEX), dimension(nxph,nyp), intent(out) :: varp
		!---------
		varp(:,1:nyph) = var(1:nxph,1:nyph);
		varp(:,nyph+1:nyp) = var(1:nxph,ny-nyph+2:ny);
	end subroutine unpadm
	!-----subroutine-------
	pure subroutine padm(nxh,nxph,ny,nyp,nyph,varp,var)
	use, intrinsic :: iso_c_binding
	implicit none
		!In:
		!Integer
		integer(C_INT), intent(in) :: nxh,nxph,ny,nyp,nyph
		!Complex
		complex(C_DOUBLE_COMPLEX), dimension(nxph,nyp), intent(in) :: varp
		!Out:
		!Complex
		complex(C_DOUBLE_COMPLEX), dimension(nxh,ny), intent(out) :: var
		!---------
		var=0d0;
		var(1:nxph,1:nyph) = varp(:,1:nyph);
		var(1:nxph,ny-nyph+2:ny) = varp(:,nyph+1:nyp);
	end subroutine padm
end module mpads

!-----------------------------MODULE--------------------------

module mainsubroutines
contains

!--- Subroutine Runge-Kuta

	subroutine RK_flux(nx,nxh,nxph,ny,nyp,nyph,m,kx,ky,dx,dz,f_star,g_star,epsbar,L,B_star,nu,vrain, &
	uhat,vhat,what,ThetaRhat,qthat,u,v,w,ThetaR,qt, &
	fu,fv,fw,fThetaR,fqt,RC,Mr,ubg,vbg,tau,qvs,qvini,fqvdz,in1,out1,planf1)
	use, intrinsic :: iso_c_binding
	use FFTW3
	use mpads
	use varsRK
	implicit none
		!Parameters:
		!Integers
		type(C_PTR),intent(in) :: planf1
		integer(C_INT) :: ix,iy,iz,ns
		!Complex
		complex(C_DOUBLE_COMPLEX) :: IM
		!Real:
		real(C_DOUBLE) :: qvbar
		!In:
		!Integers
		integer, intent(in) :: nx,nxh,nxph,ny,nyp,nyph,m
		!Real:
		real(C_DOUBLE), intent(in) :: dx,dz,vrain, tau, f_star, g_star, epsbar, L, B_star
		real(C_DOUBLE), dimension(2),intent(in) :: nu   !*** diffusion parameters 
		real(C_DOUBLE), dimension(nxph,nyp,m), intent(in) :: kx,ky
		real(C_DOUBLE), dimension(m), intent(in) :: qvs, fqvdz
		real(C_DOUBLE), dimension(m), intent(in) :: qvini,RC,Mr,ubg,vbg
		real(C_DOUBLE),dimension(nx,ny,m),intent(in) :: u,v,w,ThetaR, qt
		real(C_DOUBLE) ,dimension(nx,ny),intent(inout) ::  in1
		!Complex
		complex(C_DOUBLE_COMPLEX), dimension(nxph,nyp,m), intent(in) :: uhat,vhat,what, qthat, ThetaRhat
		complex(C_DOUBLE_COMPLEX), dimension(nxh,ny), intent(inout) :: out1
		!Out:
		!Complex:
		complex(C_DOUBLE_COMPLEX), dimension(nxph,nyp,m), intent(out) :: fu,fv,fw,fThetaR,fqt

		!---------------------------------------------
		IM = (0d0,1d0);
		ns = nx*ny;
		
		!---------------------------------------------
		!Allocation:
		!Real
		allocate(qr(nx,ny,m));
		!Complex
		allocate(wrkNL(nx,ny,m)); allocate(wrkNLhat(nxph,nyp,m));
		allocate(wrk1(nxph,nyp));
		
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
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1); !de-aliasing
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do
		fu = fu - IM*kx*wrkNLhat;
		
		!---------Non-linear advection terms: - (uv)_x and - (uv)_y---------------------------------------	
		!wrkNL = u*v;
		do iz = 1,m
			in1 = u(:,:,iz)*v(:,:,iz); !wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;		
		end do
		fu = fu - IM*ky*wrkNLhat;
		fv = fv - IM*kx*wrkNLhat;
		
		!---------Non-linear advection term: - (vv)_x----------------------------------------	
		!wrkNL = v*v;
		do iz = 1,m
			in1 = v(:,:,iz)*v(:,:,iz); !wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;		
		end do	
		fv = fv - IM*ky*wrkNLhat;
		
		!----------Non-linear advection term: - (wu)_z affecting the u equation --------------------------------------	

		!Need to compute u at the w level, due to staggered grid
		!Interpolate u:
		wrkNL(:,:,2:m-1) = 0.5d0*(u(:,:,1:m-2)+u(:,:,2:m-1)); !*** u on position as w
		!**** No boundary condition, but w=0 ****
		wrkNL(:,:,1) = 0d0;
		wrkNL(:,:,m) = 0d0;
		wrkNL = wrkNL*w;
		do iz = 1,m
			in1 = wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do	
		
		fu(:,:,2:m-2) = fu(:,:,2:m-2) - (wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,2:m-2))/dz &
					 + nu(2)*(uhat(:,:,3:m-1)+uhat(:,:,1:m-3))/dz**2d0;
					 
		!**** u=0 at the bottom, du/dz=0 on the top ****
		fu(:,:,1)     = fu(:,:,1) - (wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/dz &
				 + nu(2)*(uhat(:,:,2)-uhat(:,:,1))/dz**2d0; !u(0) = -u(1)
		fu(:,:,m-1)   = fu(:,:,m-1) - (wrkNLhat(:,:,m)-wrkNLhat(:,:,m-1))/dz &
				 + nu(2)*(uhat(:,:,m-1)+uhat(:,:,m-2))/dz**2d0;	!u(m)  = u(m-1)

		fu(1,1,2:m-2) = fu(1,1,2:m-2)-nu(2)*(uhat(1,1,3:m-1)+uhat(1,1,1:m-3))/dz**2d0; !Remove viscosity to VSHF
		fu(1,1,1) = fu(1,1,1)-nu(2)*(uhat(1,1,2)-uhat(1,1,1))/dz**2d0; !Remove viscosity to VSHF
		fu(1,1,m-1) = fu(1,1,m-1)-nu(2)*(uhat(1,1,m-1)+uhat(1,1,m-2))/dz**2d0; !Remove viscosity to VSHF
		
		!**************************************	
		!!Ralaxation terms
		!!(u,v) horizontal velocity component
		!do iz=1,m
		!	fu(1,1,iz) = fu(1,1,iz)-(1d0/tau)*(uhat(1,1,iz)-ubg(iz));
		!	fv(1,1,iz) = fv(1,1,iz)-(1d0/tau)*(vhat(1,1,iz)-vbg(iz));
		!end do
		
		!----------Non-linear advection term: - (wu)_x affecting the w equation --------------------------------------	
		
		fw = -IM*kx*wrkNLhat;
		
		!----------Non-linear advection term: - (wv)_z affecting the v equation --------------------------------------	
		!Need to compute v at the w level, due to staggered grid
		!********* Interpolation for v ********
		wrkNL(:,:,2:m-1) = 0.5d0*(v(:,:,1:m-2)+v(:,:,2:m-1)); !*** v on position as w
		!**** No boundary condition, but w=0 ****
		wrkNL(:,:,1) = 0d0;
		wrkNL(:,:,m) = 0d0;
		wrkNL = wrkNL*w;
		!----------------------------------
		do iz = 1,m
			in1 = wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do	
		!-----------------------------------
		fv(:,:,2:m-2) = fv(:,:,2:m-2) - (wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,2:m-2))/dz &
					 + nu(2)*(vhat(:,:,3:m-1)+vhat(:,:,1:m-3))/dz**2d0;
					 
		!**** v=0 at the bottom, dv/dz=0 on the top ****
		fv(:,:,1)     = fv(:,:,1) - (wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/dz &
					 + nu(2)*(vhat(:,:,2)-vhat(:,:,1))/dz**2d0;
		fv(:,:,m-1)   = fv(:,:,m-1) - (wrkNLhat(:,:,m)-wrkNLhat(:,:,m-1))/dz &
					 + nu(2)*(vhat(:,:,m-1)+vhat(:,:,m-2))/dz**2d0;	
	
		fv(1,1,2:m-2) = fv(1,1,2:m-2)-nu(2)*(vhat(1,1,3:m-1)+vhat(1,1,1:m-3))/dz**2d0; !Remove viscosity to VSHF
		fv(1,1,1) = fv(1,1,1)-nu(2)*(vhat(1,1,2)-vhat(1,1,1))/dz**2d0; !Remove viscosity to VSHF
		fv(1,1,m-1) = fv(1,1,m-1)-nu(2)*(vhat(1,1,m-1)+vhat(1,1,m-2))/dz**2d0; !Remove viscosity to VSHF
		
		!----------Non-linear advection term: - (wv)_y affecting the w equation --------------------------------------	
		fw = fw - IM*ky*wrkNLhat;
		
		!---------------------------------------------------------------------------------
		!Here we need to get qr
		qr = 0d0;
		do iz=1,m
			where ( qvini(iz)+qt(:,:,iz) > qvs(iz) ) qr(:,:,iz) = qvini(iz)+qt(:,:,iz)-qvs(iz);
		end do
		
		!------------q_v '  term affecting the w equation---------------------------------------------------------------------

		do iz = 1,m
			qvbar = sum(qt(:,:,iz)-qr(:,:,iz))/ns; !Horizontal average of qv'
			in1 = qt(:,:,iz)-qr(:,:,iz)-qvbar;
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do
		fw(:,:,2:m-1) = fw(:,:,2:m-1) + g_star*epsbar*0.5d0*(wrkNLhat(:,:,1:m-2)+wrkNLhat(:,:,2:m-1));
		
		!------------q_r '  term affecting the w equation---------------------------------------------------------------------
		do iz = 1,m
			in1 = qr(:,:,iz);
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do
		fw(:,:,2:m-1) = fw(:,:,2:m-1) - g_star*0.5d0*(wrkNLhat(:,:,1:m-2)+wrkNLhat(:,:,2:m-1));
		
		
		!------------q_r  term affecting the theta_r equation---------------------------------------------------------------------

		fThetaR = 0d0;
		fThetaR(:,:,2:m-2) = fThetaR(:,:,2:m-2)-vrain*L*(wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,1:m-3))/(2d0*dz);
		!**** Neumann boundary condition ****
		fThetaR(:,:,1) = fThetaR(:,:,1)-vrain*L*(wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/(2d0*dz);
		fThetaR(:,:,m-1) = fThetaR(:,:,m-1)-vrain*L*(wrkNLhat(:,:,m-1)-wrkNLhat(:,:,m-2))/(2d0*dz);

		!------------q_r  term affecting the q_t equation---------------------------------------------------------------------

		fqt = 0d0;
		fqt(:,:,2:m-2) = fqt(:,:,2:m-2) + vrain*(wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,1:m-3))/(2d0*dz);
		!**** Neumann boundary condition for q_r		
		fqt(:,:,1) = fqt(:,:,1) + vrain*(wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/(2d0*dz); !q_r(0) = q_r(1)
		fqt(:,:,m-1) = fqt(:,:,m-1) + vrain*(wrkNLhat(:,:,m-1)-wrkNLhat(:,:,m-2))/(2d0*dz); !q_r(m) = q_r(m-1)
		
		!-------------Theta term affecting the w equation--------------------------------------------------------------------
		do iz = 1,m
			in1 = ThetaR(:,:,iz)+L*qr(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do
		fw(:,:,2:m-1) = fw(:,:,2:m-1) + g_star*( 0.5d0*(wrkNLhat(:,:,1:m-2)+wrkNLhat(:,:,2:m-1)) );
		
		!----------Non-linear advection term: - (ww)_z affecting the w equation --------------------------------------	
		!********** Interpolation for w ***************
		wrkNL(:,:,1:m-1) = 0.5d0*(w(:,:,1:m-1)+w(:,:,2:m));  !*** w on position as u
		!**** Zero boundary condition, 0.5*(w(m+1)+w(m-1))=0, w(m)=0.0
		wrkNL(:,:,m) = -wrkNL(:,:,m-1);
		wrkNL = wrkNL*wrkNL;
	
		do iz = 1,m
			in1 = wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do	
		
		fw(:,:,2:m-1) = fw(:,:,2:m-1) - (wrkNLhat(:,:,2:m-1)-wrkNLhat(:,:,1:m-2))/dz &
		+ nu(2)*(what(:,:,3:m)+what(:,:,1:m-2))/dz**2d0 ;
	
		!fw(1,1,2:m-1) = fw(1,1,2:m-1)-nu(2)*(what(1,1,3:m)+what(1,1,1:m-2))/dz**2d0; !Remove viscosity to VSHF
	
		!**** Zero boundary condition ******	
		fw(:,:,1) = 0d0;
		fw(:,:,m) = 0d0;
		
		!----------Non-linear advection term: - (u theta_r')_x affecting the theta_r equation --------------------------------------	
		!wrkNL = u*ThetaR;
		do iz = 1,m
			in1 = u(:,:,iz)*ThetaR(:,:,iz); !wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do
		fThetaR = fThetaR-IM*kx*wrkNLhat;
		
		!----------Non-linear advection term: - (v theta_r')_y affecting the theta_r equation --------------------------------------	
		!wrkNL = v*ThetaR;
		do iz = 1,m
			in1 = v(:,:,iz)*ThetaR(:,:,iz); !wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do
		fThetaR = fThetaR - IM*ky*wrkNLhat;
		
		!---------------------------------------------------------------------------------
		!********************* For Fixed Background
		do iz = 1,m-2
			fThetaR(:,:,iz) = fThetaR(:,:,iz)-0.5d0*(what(:,:,iz)+what(:,:,iz+1))*B_star; 
		end do
		
		!----------Non-linear advection term: - (w theta_r')_z affecting the theta_r equation --------------------------------------	
		!********** Interpolation for ThetaR **************
		wrkNL(:,:,2:m-1) = 0.5d0*(ThetaR(:,:,1:m-2)+ThetaR(:,:,2:m-1)); !*** ThetaR on position as w
		!****   Neumann boundary condition ****
		wrkNL(:,:,1) = ThetaR(:,:,1);
		wrkNL(:,:,m) = ThetaR(:,:,m-1);
		wrkNL = w*wrkNL;
		!----------------------------------
		do iz = 1,m
			in1 = wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do	
		!**************************************	
		
		fThetaR(:,:,2:m-2) = fThetaR(:,:,2:m-2) - (wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,2:m-2))/dz &
		+ nu(2)*(ThetaRhat(:,:,3:m-1)+ThetaRhat(:,:,1:m-3))/dz**2d0;
		!**** Neumann boundary condition ****
		fThetaR(:,:,1) = fThetaR(:,:,1) - (wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/dz &
		+ nu(2)*(ThetaRhat(:,:,2)+ThetaRhat(:,:,1))/dz**2d0;
		fThetaR(:,:,m-1) = fThetaR(:,:,m-1) - (wrkNLhat(:,:,m)-wrkNLhat(:,:,m-1))/dz &
		+ nu(2)*(ThetaRhat(:,:,m-1)+ThetaRhat(:,:,m-2))/dz**2d0;
	
		!fThetaR(1,1,2:m-2) = fThetaR(1,1,2:m-2)-nu(2)*(ThetaRhat(1,1,3:m-1)+ThetaRhat(1,1,1:m-3))/dz**2d0; !Remove viscosity to VSHF
		!fThetaR(1,1,1) = fThetaR(1,1,1)-nu(2)*(ThetaRhat(1,1,2)+ThetaRhat(1,1,1))/dz**2d0; !Remove viscosity to VSHF
		!fThetaR(1,1,m-1) = fThetaR(1,1,m-1)-nu(2)*(ThetaRhat(1,1,m-1)+ThetaRhat(1,1,m-2))/dz**2d0; !Remove viscosity to VSHF
		
		!----------------------------------
		!Radiative cooling
		do iz=1,m
			fThetaR(1,1,iz) = fThetaR(1,1,iz)+RC(iz); !Only one eigenmode
		end do
		
		!----------Non-linear advection term: - (u q_t')_x affecting the q_t equation --------------------------------------	
		do iz = 1,m
			in1 = u(:,:,iz)*qt(:,:,iz); !wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1) 
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do		
		fqt = fqt-IM*kx*wrkNLhat;
		
		!----------Non-linear advection term: - (v q_t')_y affecting the q_t equation --------------------------------------	
		do iz = 1,m
			in1 = v(:,:,iz)*qt(:,:,iz); !wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1)
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do	
		fqt = fqt - IM*ky*wrkNLhat;
		
		!----------Non-linear advection term: - (w q_t')_z affecting the q_t equation --------------------------------------	
		!********** Interpolation for qv **************
		wrkNL(:,:,2:m-1) = 0.5d0*(qt(:,:,1:m-2)+qt(:,:,2:m-1));!*** qv on position as w
		!****"zero on bottom" Neumann boundary condition **** 
		wrkNL(:,:,1) = qt(:,:,1) ;
		wrkNL(:,:,m) = qt(:,:,m-1);
		wrkNL = w*wrkNL;
		do iz = 1,m
			in1 = wrkNL(:,:,iz); 
			call dfftw_execute_dft_r2c(planf1,in1,out1);
			call unpadm(nxh,nxph,ny,nyp,nyph,out1,wrk1) 
			wrkNLhat(:,:,iz) = wrk1/ns;
		end do	
			
		fqt(:,:,2:m-2) = fqt(:,:,2:m-2) - (wrkNLhat(:,:,3:m-1)-wrkNLhat(:,:,2:m-2))/dz &
		+ nu(2)*(qthat(:,:,3:m-1)+qthat(:,:,1:m-3))/dz**2d0;
		!****"Non-zero on bottom" Neumann boundary condition **** 		
		fqt(:,:,1) = fqt(:,:,1) - (wrkNLhat(:,:,2)-wrkNLhat(:,:,1))/dz &
		+ nu(2)*(qthat(:,:,2)+qthat(:,:,1))/dz**2d0;
		fqt(:,:,m-1) = fqt(:,:,m-1) - (wrkNLhat(:,:,m)-wrkNLhat(:,:,m-1))/dz &
		+ nu(2)*(qthat(:,:,m-1)+qthat(:,:,m-2))/dz**2d0;
	
		!fqt(1,1,2:m-2) = fqt(1,1,2:m-2)-nu(2)*(qthat(1,1,3:m-1)+qthat(1,1,1:m-3))/dz**2d0; !Remove viscosity to VSHF
		!fqt(1,1,1) = fqt(1,1,1)-nu(2)*(qthat(1,1,2)+qthat(1,1,1))/dz**2d0; !Remove viscosity to VSHF
		!fqt(1,1,m-1) = fqt(1,1,m-1)-nu(2)*(qthat(1,1,m-1)+qthat(1,1,m-2))/dz**2d0; !Remove viscosity to VSHF
		
		!----------Background term -w d_qvditle(z)/dz--------------------------------------	
		!********************* For Fixed Background
		do iz = 1,m-2
			fqt(:,:,iz) = fqt(:,:,iz)-0.5d0*(what(:,:,iz)+what(:,:,iz+1))*fqvdz(iz); 
		end do
		
		!-------------------------------------------------------------------------------------		
		!Moistening
		do iz=1,m
			fqt(1,1,iz) = fqt(1,1,iz)+Mr(iz); !Only the first eigenmode
		end do
		!-------------------------------------------------------------------------------------	
		deallocate(qr,wrk1,wrkNL,wrkNLhat);
		return;
	end subroutine RK_flux	

	!--- Subroutine Thomas for Pressure Poisson Equation ----
	subroutine Thomas(x,a,b,c,r,nxph,nyp,m);
	use, intrinsic :: iso_c_binding
	use varsThomas
	implicit none
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

end module mainsubroutines
