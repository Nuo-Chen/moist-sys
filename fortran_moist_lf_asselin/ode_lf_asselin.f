c234567
      subroutine fwd(nt, w, bu, bs, nu, ns, dt)
      integer nt, t
      real dt, nu, ns, alpha
      real w(nt), bu(nt), bs(nt)
      
      alpha = 0.01 

      nu = 0.012
      ns = 0.0012
      dt = 0.1

      w(1) = 0.1
      bs(1) = 0
      bu(1) = 0

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
	
! Apply Asselin filter 

          bu(t) = bu(t) + alpha*(bu(t+1) - 2*bu(t) + bu(t-1))
          bs(t) = bs(t) + alpha*(bs(t+1) - 2*bs(t) + bs(t-1))
          w(t) = w(t) + alpha*(w(t+1) - 2*w(t) + w(t-1))

      end do
      
      return 
      end