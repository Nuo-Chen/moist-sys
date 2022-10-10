      SUBROUTINE FWD_D(nt, w, wd, bu, bud, bs, bsd, nu, ns, dt, alpha)
      IMPLICIT NONE
      INTEGER nt, t
      REAL dt, nu, ns, alpha
      REAL w(nt), bu(nt), bs(nt)
      REAL wd(nt), bud(nt), bsd(nt)
C      alpha = 0.01
C
C      nu = 0.012
C      ns = 0.0012
C      dt = 0.1
C
C      wd(1) = 0.0
C      w(1) = 0.1
C      bsd(1) = 0.0
C      bs(1) = 0
C      bud(1) = 0.0
C      bu(1) = 0
C
C leap frog first time step      
      t = 1
      write(15,*) wd(t),bud(t), bsd(t)
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
        write(15,*) wd(t),bud(t), bsd(t)
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
      write(15,*) wd(nt),bud(nt), bsd(nt)
C
      RETURN
      END
