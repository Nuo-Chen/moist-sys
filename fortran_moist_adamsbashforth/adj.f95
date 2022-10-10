SUBROUTINE FWD_B(nt, w, wb, bu, bub, bs, bsb, nu, ns, dt)
    IMPLICIT NONE
    INTEGER :: nt, t
    REAL :: dt, nu, ns
    REAL :: dwf, dbuf, dbsf
    REAL :: dwfb, dbufb, dbsfb
    REAL :: dwpf, dbupf, dbspf
    REAL :: dwpfb, dbupfb, dbspfb
    REAL :: dwppf, dbuppf, dbsppf
    REAL :: dwppfb, dbuppfb, dbsppfb
    REAL, INTENT(INOUT) :: w(nt), bu(nt), bs(nt)
    REAL, INTENT(INOUT) :: wb(nt), bub(nt), bsb(nt)
    INTEGER :: branch
!    nu = 0.012
!    ns = 0.0012
!    dt = 0.1
!    w(1) = 0.1
!    bs(1) = 0
!    bu(1) = 0
  ! first step
    t = 1
    write(14,*) w(t),bu(t), bs(t)
    IF (bs(t) .LT. bu(t)) THEN
      dwf = bu(t)*dt
      CALL PUSHCONTROL1B(0)
    ELSE
      dwf = bs(t)*dt
      CALL PUSHCONTROL1B(1)
    END IF
    w(t+1) = w(t) + dwf
    dbuf = -(nu*nu*w(t)*dt)
    dbsf = -(ns*ns*w(t)*dt)
    bu(t+1) = bu(t) + dbuf
    bs(t+1) = bs(t) + dbsf
    dwpf = dwf
    dbspf = dbsf
    dbupf = dbuf
    
  ! second step
    t = 2
    write(14,*) w(t),bu(t), bs(t)
    IF (bs(t) .LT. bu(t)) THEN
      dwf = bu(t)*dt
      CALL PUSHCONTROL1B(0)
    ELSE
      dwf = bs(t)*dt
      CALL PUSHCONTROL1B(1)
    END IF
    w(t+1) = w(t) + dwf*1.5 - dwpf*0.5
    dbuf = -(nu*nu*w(t)*dt)
    dbsf = -(ns*ns*w(t)*dt)
    bu(t+1) = bu(t) + dbuf*1.5 - dbupf*0.5
    bs(t+1) = bs(t) + dbsf*1.5 - dbspf*0.5
    dwpf = dwf
    dbspf = dbsf
    dbupf = dbuf
    dwppf = dwpf
    dbsppf = dbspf
    dbuppf = dbupf
    CALL PUSHINTEGER4(t)
  ! leap frog          
    DO t=3,nt-1
      write(14,*) w(t),bu(t), bs(t)
      IF (bs(t) .LT. bu(t)) THEN
        dwf = bu(t)*dt
        CALL PUSHCONTROL1B(0)
      ELSE
        dwf = bs(t)*dt
        CALL PUSHCONTROL1B(1)
      END IF
      w(t+1) = w(t) + dwf*23./12. - dwpf*16./12. + dwppf*5./12.
      dbuf = -(nu*nu*w(t)*dt)
      dbsf = -(ns*ns*w(t)*dt)
      bu(t+1) = bu(t) + dbuf*23./12. - dbupf*16./12. + dbuppf*5./12.
      bs(t+1) = bs(t) + dbsf*23./12. - dbspf*16./12. + dbsppf*5./12.
      dwpf = dwf
      dbspf = dbsf
      dbupf = dbuf
      dwppf = dwpf
      dbsppf = dbspf
      dbuppf = dbupf
    END DO
    write(6, *) 'time nlm: ', t
    write(14,*) w(t),bu(t), bs(t)
  
  
    dbspfb = 0.0
    dbuppfb = 0.0
    dbupfb = 0.0
    dbsppfb = 0.0
    dwpfb = 0.0
    dwppfb = 0.0
  
  
  
    DO t=nt-1,3,-1
  
      write(16,*) wb(t+1),bub(t+1), bsb(t+1)
  
      dbupfb = dbupfb + dbuppfb
      dbspfb = dbspfb + dbsppfb
      dwpfb = dwpfb + dwppfb
      dbufb = dbupfb + 23.*bub(t+1)/12.
      dbsfb = dbspfb + 23.*bsb(t+1)/12.
      dwfb = dwpfb + 23.*wb(t+1)/12.
      bsb(t) = bsb(t) + bsb(t+1)
      dbsppfb = 5.*bsb(t+1)/12.
      dbspfb = -(16.*bsb(t+1)/12.)
      bsb(t+1) = 0.0
      bub(t) = bub(t) + bub(t+1)
      dbuppfb = 5.*bub(t+1)/12.
      dbupfb = -(16.*bub(t+1)/12.)
      bub(t+1) = 0.0
      wb(t) = wb(t) + wb(t+1) - dt*ns**2*dbsfb - dt*nu**2*dbufb
      dwppfb = 5.*wb(t+1)/12.
      dwpfb = -(16.*wb(t+1)/12.)
      wb(t+1) = 0.0
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        bub(t) = bub(t) + dt*dwfb
      ELSE
        bsb(t) = bsb(t) + dt*dwfb
      END IF
    END DO
    write(16,*) wb(t+1),bub(t+1), bsb(t+1)
    CALL POPINTEGER4(t)
  
    t = 2
    dbupfb = dbupfb + dbuppfb
    dbspfb = dbspfb + dbsppfb
    dwpfb = dwpfb + dwppfb
    dbufb = dbupfb + 1.5*bub(t+1)
    dbsfb = dbspfb + 1.5*bsb(t+1)
    dwfb = dwpfb + 1.5*wb(t+1)
    bsb(t) = bsb(t) + bsb(t+1)
    dbspfb = -(0.5*bsb(t+1))
    bsb(t+1) = 0.0
    bub(t) = bub(t) + bub(t+1)
    dbupfb = -(0.5*bub(t+1))
    bub(t+1) = 0.0
    wb(t) = wb(t) + wb(t+1) - dt*ns**2*dbsfb - dt*nu**2*dbufb
    dwpfb = -(0.5*wb(t+1))
    wb(t+1) = 0.0
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      bub(t) = bub(t) + dt*dwfb
    ELSE
      bsb(t) = bsb(t) + dt*dwfb
    END IF
    write(16,*) wb(t),bub(t), bsb(t)
  
    t = 1
    dbufb = dbupfb + bub(t+1)
    dbsfb = dbspfb + bsb(t+1)
    dwfb = dwpfb + wb(t+1)
    bsb(t) = bsb(t) + bsb(t+1)
    bsb(t+1) = 0.0
    bub(t) = bub(t) + bub(t+1)
    bub(t+1) = 0.0
    wb(t) = wb(t) + wb(t+1) - dt*ns**2*dbsfb - dt*nu**2*dbufb
    wb(t+1) = 0.0
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      bub(t) = bub(t) + dt*dwfb
    ELSE
      bsb(t) = bsb(t) + dt*dwfb
    END IF
    
    write(6, *) 'time adj: ', t
    write(16,*) wb(t),bub(t), bsb(t)
  
!    bub(1) = 0.0
!    bsb(1) = 0.0
!    wb(1) = 0.0
    return 
END SUBROUTINE FWD_B
