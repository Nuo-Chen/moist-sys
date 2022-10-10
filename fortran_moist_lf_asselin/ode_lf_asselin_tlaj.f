C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.16 (develop) - 31 May 2021 11:17
C
C  Differentiation of fwd in reverse (adjoint) mode:
C   gradient     of useful results: w bs bu
C   with respect to varying inputs: w bs bu
C   RW status of diff variables: w:in-out bs:in-out bu:in-out
C234567
      SUBROUTINE FWD_B(nt, w, wb, bu, bub, bs, bsb, nu, ns, dt)
      IMPLICIT NONE
      INTEGER nt, t
      REAL dt, nu, ns, alpha
      REAL w(nt), bu(nt), bs(nt)
      REAL wb(nt), bub(nt), bsb(nt)
      REAL tempb
      INTEGER branch
      alpha = 0.01
C
      nu = 0.012
      ns = 0.0012
      dt = 0.1
C
      w(1) = 0.1
      bs(1) = 0
      bu(1) = 0
C
C leap frog first time step      
      t = 1
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
C
C leap frog          
      DO t=2,nt-1
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
C Apply Asselin filter 
C
        bu(t) = bu(t) + alpha*(bu(t+1)-2*bu(t)+bu(t-1))
        bs(t) = bs(t) + alpha*(bs(t+1)-2*bs(t)+bs(t-1))
        w(t) = w(t) + alpha*(w(t+1)-2*w(t)+w(t-1))
      ENDDO
      DO t=nt-1,2,-1
        tempb = alpha*wb(t)
        wb(t+1) = wb(t+1) + tempb
        wb(t) = wb(t) - 2*tempb
        wb(t-1) = wb(t-1) + tempb
        tempb = alpha*bsb(t)
        bsb(t+1) = bsb(t+1) + tempb
        bsb(t) = bsb(t) - 2*tempb
        bsb(t-1) = bsb(t-1) + tempb + bsb(t+1)
        tempb = alpha*bub(t)
        bub(t+1) = bub(t+1) + tempb
        bub(t) = bub(t) - 2*tempb
        bub(t-1) = bub(t-1) + tempb + bub(t+1)
        wb(t) = wb(t) - ns**2*dt*2*bsb(t+1) - nu**2*dt*2*bub(t+1)
        bsb(t+1) = 0.0
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
      bub(1) = 0.0
      bsb(1) = 0.0
      wb(1) = 0.0
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE FWD_D(nt, w, wd, bu, bud, bs, bsd, nu, ns, dt)
      IMPLICIT NONE
      INTEGER nt, t
      REAL dt, nu, ns, alpha
      REAL w(nt), bu(nt), bs(nt)
      REAL wd(nt), bud(nt), bsd(nt)
      alpha = 0.01
C
      nu = 0.012
      ns = 0.0012
      dt = 0.1
C
      wd(1) = 0.0
      w(1) = 0.1
      bsd(1) = 0.0
      bs(1) = 0
      bud(1) = 0.0
      bu(1) = 0
C
C leap frog first time step      
      t = 1
      IF (bs(t) .LT. bu(t)) THEN
        wd(t+1) = wd(t) + dt*bud(t)
        w(t+1) = w(t) + bu(t)*dt
      ELSE
        wd(t+1) = wd(t) + dt*bsd(t)
        w(t+1) = w(t) + bs(t)*dt
      END IF
C
      bud(t+1) = bud(t) - nu**2*dt*wd(t)
      bu(t+1) = bu(t) - nu*nu*w(t)*dt
      bsd(t+1) = bsd(t) - ns**2*dt*wd(t)
      bs(t+1) = bs(t) - ns*ns*w(t)*dt
C
C leap frog          
      DO t=2,nt-1
        IF (bs(t) .LT. bu(t)) THEN
          wd(t+1) = wd(t-1) + dt*2*bud(t)
          w(t+1) = w(t-1) + bu(t)*2*dt
        ELSE
          wd(t+1) = wd(t-1) + dt*2*bsd(t)
          w(t+1) = w(t-1) + bs(t)*2*dt
        END IF
C
        bud(t+1) = bud(t-1) - dt*2*nu**2*wd(t)
        bu(t+1) = bu(t-1) - nu*nu*w(t)*2*dt
        bsd(t+1) = bsd(t-1) - dt*2*ns**2*wd(t)
        bs(t+1) = bs(t-1) - ns*ns*w(t)*2*dt
C Apply Asselin filter 
C
        bud(t) = bud(t) + alpha*(bud(t+1)-2*bud(t)+bud(t-1))
        bu(t) = bu(t) + alpha*(bu(t+1)-2*bu(t)+bu(t-1))
        bsd(t) = bsd(t) + alpha*(bsd(t+1)-2*bsd(t)+bsd(t-1))
        bs(t) = bs(t) + alpha*(bs(t+1)-2*bs(t)+bs(t-1))
        wd(t) = wd(t) + alpha*(wd(t+1)-2*wd(t)+wd(t-1))
        w(t) = w(t) + alpha*(w(t+1)-2*w(t)+w(t-1))
      ENDDO
C
      RETURN
      END
