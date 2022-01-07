c234567
      program driver

      IMPLICIT NONE
      INTEGER nt, t, j, it
      parameter(nt = 5000)
      REAL dt, nu, ns, right, left
      REAL w(nt+1), bu(nt+1), bs(nt+1)
      REAL wd(nt), bud(nt), bsd(nt)
      REAL wb(nt+1), bub(nt+1), bsb(nt+1)
      REAL w0, bu0, bs0
      REAL w0d, bu0d, bs0d
      REAL w0b, bu0b, bs0b
      
      nu = 0.012
      ns = 0.0012
      dt = 0.1
C                                                                                                                                 

      w0d = 0.01
      w0 = 0.1
      bs0d = 0.02
      bs0 = 0.15
      bu0d = 0.01
      bu0 = 0.1


      w(1) = w0
      bu(1) = bu0
      bs(1) = bs0
      wd(1) = w0d
      bud(1) = bu0d
      bsd(1) = bs0d
      call FWD_D(nt, w, wd, bu, bud, bs, bsd, nu, ns, dt)

      wb(nt+1)= wd(nt)
      bub(nt+1) = bud(nt)
      bsb(nt+1) = bsd(nt)

      write(6,*) bud(nt), bsd(nt)

      call FWD_B(nt, w, wb, bu, bub, bs, bsb, nu, ns, dt)

                
      left = wd(nt)*wd(nt)+bud(nt)*bud(nt) + bsd(nt)*bsd(nt)
      right = wb(1)*w0d + bub(1)*bu0d + bsb(1)*bs0d
      write(6,*) 'left  = ', left
      write(6,*) 'right = ', right
      write(6,*) wb(1), bub(1), bsb(1)
      
      stop
      end
      

CCCC
