c234567
      subroutine fwd(nt, w0, bu0, bs0, w, bu, bs, nu, ns, dt)
      integer nt, t
      real dt, nu, ns
      real w(nt), bu(nt), bs(nt)
      
      
      nu = 0.012
      ns = 0.0012
      dt = 0.1

      w(1) = w0
      bs(1) = bs0
      bu(1) = bu0
! leap frog first time step      
      t = 1
      if (bs(t) .lt. bu(t)) then
              w(t+1) = w(t) + bu(t)*dt
          else
              w(t+1) = w(t) + bs(t)*dt
          end if

      bu(t+1) = bu(t) - nu*nu*w(t)*dt
      bs(t+1) = bs(t) - ns*ns*w(t)*dt
! leap frog          
      do t = 2, nt-1
          if (bs(t) .lt. bu(t)) then
              w(t+1) = w(t-1) + bu(t)*2*dt
          else
              w(t+1) = w(t-1) + bs(t)*2*dt
          end if

          bu(t+1) = bu(t-1) - nu*nu*w(t)*2*dt
          bs(t+1) = bs(t-1) - ns*ns*w(t)*2*dt

      end do
      
      return 
      end
