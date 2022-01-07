! test code for modified leap frog
      SUBROUTINE FWD_B(nt, w0, w0b, bu0, bu0b, bs0, bs0b, w, wb, bu, bub
     +                 , bs, bsb, nu, ns, dt)
!      SUBROUTINE FWD_B(nt, w, wb, bu, bub, bs, bsb, nu, ns, dt)
      IMPLICIT NONE
      INTEGER nt, t
      REAL dt, nu, ns
      REAL w(nt), bu(nt), bs(nt)
      REAL wb(nt+1), bub(nt+1), bsb(nt+1)
      INTEGER branch
      
      nu = 0.12
      ns = 0.012
      dt = 0.1
C

      
C leap frog first time step      
      t = 1
      write(14,*) w(t),bu(t), bs(t)
      IF (bs(t) .LT. bu(t)) THEN
        w(t+1) = w(t) + bu(t)*dt
        CALL PUSHCONTROL1B(0)
      ELSE
        w(t+1) = w(t) + bs(t)*dt
        CALL PUSHCONTROL1B(1)
      END IF
C
      bu(t+1) = bu(t) - nu*nu*w(t)*dt
      bs(t+1) = bs(t) - ns*ns*w(t)*dt
      CALL PUSHINTEGER4(t)
C leap frog          
      DO t=2,nt-1
	write(14,*) w(t),bu(t), bs(t)
        IF (bs(t) .LT. bu(t)) THEN
          w(t+1) = w(t-1) + bu(t)*2*dt
          CALL PUSHCONTROL1B(0)
        ELSE
          w(t+1) = w(t-1) + bs(t)*2*dt
          CALL PUSHCONTROL1B(1)
        END IF
C
        bu(t+1) = bu(t-1) - nu*nu*w(t)*2*dt
        bs(t+1) = bs(t-1) - ns*ns*w(t)*2*dt
      ENDDO
      write(14,*) w(nt),bu(nt), bs(nt)
C
C adjoint

! y(n+1) =  \zeta(n)
      wb(nt+1) = w0b
      bub(nt+1) = bu0b 
      bsb(nt+1) = bs0b
! yn = 2*dt*P^T(n)*y(n+1)
      t=nt
!        write(16,*) wb(t+1),bub(t+1), bsb(t+1)
!------------- this should be P^T(n) ?--------------     
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
!------------- this should be P^T(n) -------------- 
! y(n) = 0.5* y(n) + y(n+1)
      wb(nt) = 0.5*wb(nt) + wb(nt+1)
      bub(nt) = 0.5*bub(nt) + bub(nt+1)
      bsb(nt) = 0.5*bsb(nt) + bsb(nt+1)

! y(j) = y(j+2) + 2*dt*P^T(j)*y(j+1)
! including y(2) = y(4) + 2*dt*P^T(2)*y(3) ?
! y(2) should be the sensitivity at the first time step? y(n+1) be the last?
      DO t=nt,2,-1
!        write(16,*) wb(t+1),bub(t+1), bsb(t+1)
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
      ENDDO
!      write(6,*) 't+1 = ', t+1
!      write(16,*) wb(t+1),bub(t+1), bsb(t+1)
      
      CALL POPINTEGER4(t)