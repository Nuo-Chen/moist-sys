SUBROUTINE FWD_D(nt, w, wd, bu, bud, bs, bsd, nu, ns, dt)
    IMPLICIT NONE
    INTEGER :: nt, t
    REAL :: dt, nu, ns
    REAL :: dwf, dbuf, dbsf
    REAL :: dwfd, dbufd, dbsfd
    REAL :: dwpf, dbupf, dbspf
    REAL :: dwpfd, dbupfd, dbspfd
    REAL :: dwppf, dbuppf, dbsppf
    REAL :: dwppfd, dbuppfd, dbsppfd
    REAL, INTENT(INOUT) :: w(nt), bu(nt), bs(nt)
    REAL, INTENT(INOUT) :: wd(nt), bud(nt), bsd(nt)
  
  ! first step
    t = 1
    write(15,*) wd(t),bud(t), bsd(t)
    IF (bs(t) .LT. bu(t)) THEN
      dwfd = dt*bud(t)
      dwf = bu(t)*dt
    ELSE
      dwfd = dt*bsd(t)
      dwf = bs(t)*dt
    END IF
    wd(t+1) = wd(t) + dwfd
    w(t+1) = w(t) + dwf
  
    dbufd = -(nu**2*dt*wd(t))
    dbuf = -(nu*nu*w(t)*dt)
    dbsfd = -(ns**2*dt*wd(t))
    dbsf = -(ns*ns*w(t)*dt)
    bud(t+1) = bud(t) + dbufd
    bu(t+1) = bu(t) + dbuf
    bsd(t+1) = bsd(t) + dbsfd
    bs(t+1) = bs(t) + dbsf
  
    dwpfd = dwfd
    dwpf = dwf
    dbspfd = dbsfd
    dbspf = dbsf
    dbupfd = dbufd
    dbupf = dbuf
  
  ! second step
    t = 2
    write(15,*) wd(t),bud(t), bsd(t)
    IF (bs(t) .LT. bu(t)) THEN
      dwfd = dt*bud(t)
      dwf = bu(t)*dt
    ELSE
      dwfd = dt*bsd(t)
      dwf = bs(t)*dt
    END IF
    wd(t+1) = wd(t) + 1.5*dwfd - 0.5*dwpfd
    w(t+1) = w(t) + dwf*1.5 - dwpf*0.5
  
    dbufd = -(nu**2*dt*wd(t))
    dbuf = -(nu*nu*w(t)*dt)
    dbsfd = -(ns**2*dt*wd(t))
    dbsf = -(ns*ns*w(t)*dt)
    bud(t+1) = bud(t) + 1.5*dbufd - 0.5*dbupfd
    bu(t+1) = bu(t) + dbuf*1.5 - dbupf*0.5
    bsd(t+1) = bsd(t) + 1.5*dbsfd - 0.5*dbspfd
    bs(t+1) = bs(t) + dbsf*1.5 - dbspf*0.5
  
    dwpfd = dwfd
    dwpf = dwf
    dbspfd = dbsfd
    dbspf = dbsf
    dbupfd = dbufd
    dbupf = dbuf
    dwppfd = dwpfd
    dwppf = dwpf
    dbsppfd = dbspfd
    dbsppf = dbspf
    dbuppfd = dbupfd
    dbuppf = dbupf
  
  ! leap frog          
    DO t=3,nt-1
      write(15,*) wd(t),bud(t), bsd(t)
      IF (bs(t) .LT. bu(t)) THEN
        dwfd = dt*bud(t)
        dwf = bu(t)*dt
      ELSE
        dwfd = dt*bsd(t)
        dwf = bs(t)*dt
      END IF
      wd(t+1) = wd(t) + 23.*dwfd/12. + 5.*dwppfd/12. - 16.*dwpfd/12.
      w(t+1) = w(t) + dwf*23./12. - dwpf*16./12. + dwppf*5./12.
  
      dbufd = -(nu**2*dt*wd(t))
      dbuf = -(nu*nu*w(t)*dt)
      dbsfd = -(ns**2*dt*wd(t))
      dbsf = -(ns*ns*w(t)*dt)
      bud(t+1) = bud(t) + 23.*dbufd/12. + 5.*dbuppfd/12. - 16.*dbupfd/12.
      bu(t+1) = bu(t) + dbuf*23./12. - dbupf*16./12. + dbuppf*5./12.
      bsd(t+1) = bsd(t) + 23.*dbsfd/12. + 5.*dbsppfd/12. - 16.*dbspfd/12.
      bs(t+1) = bs(t) + dbsf*23./12. - dbspf*16./12. + dbsppf*5./12.
  
      dwpfd = dwfd
      dwpf = dwf
      dbspfd = dbsfd
      dbspf = dbsf
      dbupfd = dbufd
      dbupf = dbuf
      dwppfd = dwpfd
      dwppf = dwpf
      dbsppfd = dbspfd
      dbsppf = dbspf
      dbuppfd = dbupfd
      dbuppf = dbupf
    END DO
    write(15,*) wd(t),bud(t), bsd(t)
    return
  END SUBROUTINE FWD_D