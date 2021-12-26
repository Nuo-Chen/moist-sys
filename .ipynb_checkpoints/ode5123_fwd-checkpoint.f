c234567
      subroutine fwd(nt, w, bu, bs, nu, ns, dt)
      integer nt, t
      real dt, nu, ns
      real w(nt), bu(nt), bs(nt)
      
      
      nu = 0.012
      ns = 0.0012
      dt = 0.1

      w(1) = 0.1
      bs(1) = 0
      bu(1) = 0
      
      do t = 1, nt-1
          if (bs(t) .lt. bu(t)) then
              w(t+1) = w(t) + bu(t)*dt
          else
              w(t+1) = w(t) + bs(t)*dt
          end if

          bu(t+1) = bu(t) - nu*nu*w(t)*dt
          bs(t+1) = bs(t) - ns*ns*w(t)*dt

      end do
      
      return 
      end