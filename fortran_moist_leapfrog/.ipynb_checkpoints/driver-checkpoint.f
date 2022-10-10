C 14 for NLM
C 15 for TLM
C 16 for ADJ
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
      
      nu = 0.12
      ns = 0.012
      dt = 0.1
C                                                                                                                                 

      w0d = 0.0001
      w0 = 0.1
      bs0d = 0.0001
      bs0 = 0.1
      bu0d = 0.0001
      bu0 = 0.5


      w(1) = w0
      bu(1) = bu0
      bs(1) = bs0
      wd(1) = w0d
      bud(1) = bu0d
      bsd(1) = bs0d
      call FWD_D(nt, w, wd, bu, bud, bs, bsd, nu, ns, dt)

      wb(nt)= wd(nt)
      bub(nt) = bud(nt)
      bsb(nt) = bsd(nt)

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
