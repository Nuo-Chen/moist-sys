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
! leap frog first time step      
      t = 1
      
      w(t+1) = w(t) + nu*nu*bu(t)*dt + ns*ns*bs(t)*dt
      if (bs(t) .lt. bu(t)) then
              bu(t+1) = bu(t) - w(t)*dt
          else
              bs(t+1) = bs(t) - w(t)*dt
          end if

! leap frog          
      do t = 2, nt-1
          w(t+1) = w(t-1) + nu*nu*bu(t)*2*dt + ns*ns*bs(t)*2*dt
          
          if (bs(t) .lt. bu(t)) then
              bu(t+1) = bu(t) - w(t)*dt*2
          else
              bu(t+1) = bu(t) - w(t)*dt*2
          end if
          
      end do
      
      return 
      end