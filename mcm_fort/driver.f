c234567
      program driver

      IMPLICIT NONE
      INTEGER nt, t, j, it
      parameter(nt = 50000)
      REAL dt, nu, ns, right, left
      REAL w(nt), bu(nt), bs(nt)
      REAL wd(nt), bud(nt), bsd(nt)
      REAL wb(nt), bub(nt), bsb(nt)
      REAL w0, bu0, bs0
      REAL w0d, bu0d, bs0d
      REAL w0b, bu0b, bs0b
      
      nu = 0.012
      ns = 0.0012
      dt = 0.1
C                                                                                                                                 

      w0d = 0.1
      w0 = 0.1
      bs0d = 0.1
      bs0 = 0
      bu0d = 0.1
      bu0 = 0

      call FWD_D(nt, w0, w0d, bu0, bu0d, bs0, bs0d, w, wd, bu, bud
     +                 , bs, bsd, nu, ns, dt)


      wb(nt)= wd(nt)
      bub(nt) = bud(nt)
      bsb(nt) = bsd(nt)

      write(6,*) bud(nt), bsd(nt)

      call FWD_B(nt, w0, w0b, bu0, bu0b, bs0, bs0b, w, wb, bu, bub
     +                 , bs, bsb, nu, ns, dt)

                
      left = wd(nt)*wd(nt)+bud(nt)*bud(nt) + bsd(nt)*bsd(nt)
      right = w0b*w0d + bu0b*bu0d + bs0b*bs0d
      write(6,*) 'left  = ', left
      write(6,*) 'right = ', right
      write(6,*) wb(1), bub(1), bsb(1)
      
      stop
      end
      

CCCC
