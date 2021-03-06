C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.16 (develop) - 31 May 2021 11:17
C
C  Differentiation of fwd in reverse (adjoint) mode:
C   gradient     of useful results: w bs bu
C   with respect to varying inputs: w w0 bs bu0 bu bs0
C   RW status of diff variables: w:in-out w0:out bs:in-out bu0:out
C                bu:in-out bs0:out
C234567
C      SUBROUTINE FWD_B(nt, w0, w0b, bu0, bu0b, bs0, bs0b, w, wb, bu, bub
C     +                 , bs, bsb, nu, ns, dt)
      SUBROUTINE FWD_B(nt, w, wb, bu, bub, bs, bsb, nu, ns, dt)
      IMPLICIT NONE
      INTEGER nt, t
      REAL dt, nu, ns
      REAL w(nt), bu(nt), bs(nt)
      REAL wb(nt), bub(nt), bsb(nt)
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
      DO t=nt-1,2,-1
        write(16,*) wb(t+1),bub(t+1), bsb(t+1)
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
      write(6,*) 't+1 = ', t+1
      write(16,*) wb(t+1),bub(t+1), bsb(t+1)
      
      CALL POPINTEGER4(t)
      
      t = 1
      bsb(t) = bsb(t) + bsb(t+1)
      wb(t) = wb(t) - dt*ns**2*bsb(t+1) - dt*nu**2*bub(t+1)
      bsb(t+1) = 0.0
      bub(t) = bub(t) + bub(t+1)
      bub(t+1) = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        wb(t) = wb(t) + wb(t+1)
        bub(t) = bub(t) + dt*wb(t+1)
        wb(t+1) = 0.0
      ELSE
        wb(t) = wb(t) + wb(t+1)
        bsb(t) = bsb(t) + dt*wb(t+1)
        wb(t+1) = 0.0
      END IF
      
      write(16,*) wb(1),bub(1), bsb(1)

      END
