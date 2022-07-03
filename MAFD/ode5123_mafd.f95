! test code for modified leap frog
!      SUBROUTINE FWD_B(nt, w0, w0b, bu0, bu0b, bs0, bs0b, w, wb, bu, bub &
!                      , bs, bsb, nu, ns, dt)
      SUBROUTINE FWD_B(nt, w, wb, bu, bub, bs, bsb, nu, ns, dt)
      IMPLICIT NONE
      INTEGER nt, t
      REAL dt, nu, ns
      REAL w(nt+1), bu(nt+1), bs(nt+1)
      REAL wb(nt+1), bub(nt+1), bsb(nt+1)
      INTEGER branch
      
!      nu = 0.12
!      ns = 0.012
!      dt = 0.1
!

      
! leap frog first time step      
      t = 1
      write(14,*) w(t),bu(t), bs(t)
      IF (bs(t) .LT. bu(t)) THEN
        w(t+1) = w(t) + bu(t)*dt
        CALL PUSHCONTROL1B(0)
      ELSE
        w(t+1) = w(t) + bs(t)*dt
        CALL PUSHCONTROL1B(1)
      END IF
!
      bu(t+1) = bu(t) - nu*nu*w(t)*dt
      bs(t+1) = bs(t) - ns*ns*w(t)*dt
      CALL PUSHINTEGER4(t)
! leap frog          
      DO t=2,nt
	write(14,*) w(t),bu(t), bs(t)
        IF (bs(t) .LT. bu(t)) THEN
          w(t+1) = w(t-1) + bu(t)*2*dt
          CALL PUSHCONTROL1B(0)
        ELSE
          w(t+1) = w(t-1) + bs(t)*2*dt
          CALL PUSHCONTROL1B(1)
        END IF
!
        bu(t+1) = bu(t-1) - nu*nu*w(t)*2*dt
        bs(t+1) = bs(t-1) - ns*ns*w(t)*2*dt
      ENDDO
!      write(14,*) w(nt),bu(nt), bs(nt)
!
! ------------- ADJOINT PART--------------  
! -----------------------------------------------------------------

! y(n+1) =  \zeta(n)                                -> y(n+1)
!      wb(nt+1) = w0b
!      bub(nt+1) = bu0b 
!      bsb(nt+1) = bs0b

! yn = 2*dt*P^T(n)*y(n+1)
      t = nt
      write(16,*) wb(t+1),bub(t+1), bsb(t+1)
      !bsb(t-1) = bsb(t-1) + bsb(t+1)
      wb(t) = - ns**2*dt*2*bsb(t+1) - nu**2*dt*2*bub(t+1)
      !bub(t-1) = bub(t-1) + bub(t+1)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        !wb(t-1) = wb(t-1) + wb(t+1)
        bub(t) = dt*2*wb(t+1)
      ELSE
        !wb(t-1) = wb(t-1) + wb(t+1)
        bsb(t) = dt*2*wb(t+1)
      END IF

! y(n) = 0.5* y(n) + y(n+1)                          -> y(n)
      wb(nt) = 0.5*wb(nt) + wb(nt+1)
      bub(nt) = 0.5*bub(nt) + bub(nt+1)
      bsb(nt) = 0.5*bsb(nt) + bsb(nt+1)

      write(16,*) wb(nt),bub(nt), bsb(nt)

! y(j) = y(j+2) + 2*dt*P^T(j)*y(j+1)                  -> y(n-1), y(n-2), ..., y(2), y(1?)
! including y(2) = y(4) + 2*dt*P^T(2)*y(3) ?          -> t=n-1, n-2, ..., 2, 1(except for t=1 no need for t-1=t+1)
!                                                   ->  t+1 = nt, nt-1
!                                                   -> when t = nt-1,  wb(t) = wb(t) -2dt*ns^2*bs(t+1) - 2dt*nu^2*bu(t+1)
!                                                                      wb(t) = wb(t+2) -2dt*ns^2*bs(t+1) - 2dt*nu^2*bu(t+1)
!                                                                      wb(t) = wb(t+2) is needed before the loop starts, 
!                                                                      wb(nt-1) = wb(nt+1)  same for bub(nt-1)=bub(nt+1) and bsb(nt-1)=bsb(nt+1)
! zb(t) = zb(t+2) + (zb(t)(=0))
      wb(nt-1) = wb(nt+1)
      bub(nt-1) = bub(nt+1)
      bsb(nt-1) = bsb(nt+1)

! not sure if this is needed?
      wb(t+1) = 0.0
      bub(t+1) = 0.0
      bsb(t+1) = 0.0

! start looping                                       
!      -> solve for y(t) at each loop, not y(t-1), therefore print y(t) at the end of the loop
      DO t=nt-1,2,-1
        
        bsb(t-1) = bsb(t-1) + bsb(t+1)       
        wb(t) = wb(t) - ns**2*dt*2*bsb(t+1) - nu**2*dt*2*bub(t+1)
        bsb(t+1) = 0.0
        bub(t-1) = bub(t-1) + bub(t+1)
        bub(t+1) = 0.0
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          wb(t-1) = wb(t-1) + wb(t+1)
          bub(t) = bub(t) + dt*2*wb(t+1)
          wb(t+1) = 0.0
        ELSE
          wb(t-1) = wb(t-1) + wb(t+1)
          bsb(t) = bsb(t) + dt*2*wb(t+1)
          wb(t+1) = 0.0
        END IF

        write(16,*) wb(t),bub(t), bsb(t)
      ENDDO
      write(6,*) 't = ', t
      
!      CALL POPINTEGER4(t)
      end