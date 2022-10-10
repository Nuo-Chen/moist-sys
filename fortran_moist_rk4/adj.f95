!  Differentiation of fwd in reverse (adjoint) mode:
!   gradient     of useful results: w bs bu
!   with respect to varying inputs: w bs bu
!   RW status of diff variables: w:in-out bs:in-out bu:in-out
SUBROUTINE FWD_B(nt, w, wb, bu, bub, bs, bsb, nu, ns, dt)
  IMPLICIT NONE
  INTEGER :: nt, t
  REAL :: dt, nu, ns
  REAL :: dbu, dbs, dw
  REAL :: dbub, dbsb, dwb
  REAL :: k1_bu, k1_bs, k1_w
  REAL :: k1_bub, k1_bsb, k1_wb
  REAL :: k2_bu, k2_bs, k2_w
  REAL :: k2_bub, k2_bsb, k2_wb
  REAL :: k3_bu, k3_bs, k3_w
  REAL :: k3_bub, k3_bsb, k3_wb
  REAL :: k4_bu, k4_bs, k4_w
  REAL :: k4_bub, k4_bsb, k4_wb
  REAL :: xtmp_bu, xtmp_bs, xtmp_w
  REAL :: xtmp_bub, xtmp_bsb, xtmp_wb
  REAL, INTENT(INOUT) :: w(nt), bu(nt), bs(nt)
  REAL, INTENT(INOUT) :: wb(nt), bub(nt), bsb(nt)
  REAL :: tempb
  INTEGER :: branch
! rhs 
  ! nu = 0.01
  ! ns = 0.005
  ! dt = 0.1
  ! w(1) = 0.1
  ! bu(1) = 5
  ! bs(1) = 1

  DO t=1,nt-1
    write(14,*) w(t),bu(t), bs(t)
! rhs 
    IF (bs(t) .LT. bu(t)) THEN
      dw = bu(t)
      CALL PUSHCONTROL1B(0)
    ELSE
      dw = bs(t)
      CALL PUSHCONTROL1B(1)
    END IF
    dbu = -(nu*nu*w(t))
    dbs = -(ns*ns*w(t))
! k1=self.rhs(self.xvar)*self.dt
    k1_bu = dbu*dt
    k1_bs = dbs*dt
    k1_w = dw*dt
! xtmp=self.xvar+0.5*k1
    xtmp_w = w(t) + 0.5*k1_w
    xtmp_bu = bu(t) + 0.5*k1_bu
    xtmp_bs = bs(t) + 0.5*k1_bs
! rhs 
    IF (xtmp_bs .LT. xtmp_bu) THEN
      dw = xtmp_bu
      CALL PUSHCONTROL1B(0)
    ELSE
      dw = xtmp_bs
      CALL PUSHCONTROL1B(1)
    END IF
    dbu = -(nu*nu*xtmp_w)
    dbs = -(ns*ns*xtmp_w)
! k2=self.rhs(xtmp)*self.dt
    k2_bu = dbu*dt
    k2_bs = dbs*dt
    k2_w = dw*dt
! xtmp=self.xvar+0.5*k2
    xtmp_w = w(t) + 0.5*k2_w
    xtmp_bu = bu(t) + 0.5*k2_bu
    xtmp_bs = bs(t) + 0.5*k2_bs
! rhs 
    IF (xtmp_bs .LT. xtmp_bu) THEN
      dw = xtmp_bu
      CALL PUSHCONTROL1B(0)
    ELSE
      dw = xtmp_bs
      CALL PUSHCONTROL1B(1)
    END IF
    dbu = -(nu*nu*xtmp_w)
    dbs = -(ns*ns*xtmp_w)
! k3=self.rhs(xtmp)*self.dt
    k3_bu = dbu*dt
    k3_bs = dbs*dt
    k3_w = dw*dt
! xtmp=self.xvar+k2
    xtmp_w = w(t) + 0.5*k2_w
    xtmp_bu = bu(t) + 0.5*k2_bu
    xtmp_bs = bs(t) + 0.5*k2_bs
! rhs 
    IF (xtmp_bs .LT. xtmp_bu) THEN
      dw = xtmp_bu
      CALL PUSHCONTROL1B(0)
    ELSE
      dw = xtmp_bs
      CALL PUSHCONTROL1B(1)
    END IF
    dbu = -(nu*nu*xtmp_w)
    dbs = -(ns*ns*xtmp_w)
! k4=self.rhs(xtmp)*self.dt
    k4_bu = dbu*dt
    k4_bs = dbs*dt
    k4_w = dw*dt
! self.xvar+=(k1+2*k2+2*k3+k4)/6.
    bu(t+1) = bu(t) + (k1_bu+2*k2_bu+2*k3_bu+k4_bu)/6.
    bs(t+1) = bs(t) + (k1_bs+2*k2_bs+2*k3_bs+k4_bs)/6.
    w(t+1) = w(t) + (k1_w+2*k2_w+2*k3_w+k4_w)/6.
  END DO
  write(6,*) 'time nlm = ', t
  write(14,*) w(t),bu(t), bs(t)


!! ######################### begin adjoint ##################
  !! ###########################################
  !! ###########################################
  DO t=nt-1,1,-1

    write(16,*) wb(t+1),bub(t+1), bsb(t+1)

    wb(t) = wb(t) + wb(t+1)
    tempb = wb(t+1)/6.
    wb(t+1) = 0.0
    k1_wb = tempb
    k2_wb = 2*tempb
    k3_wb = 2*tempb
    k4_wb = tempb
    bsb(t) = bsb(t) + bsb(t+1)
    tempb = bsb(t+1)/6.
    bsb(t+1) = 0.0
    k1_bsb = tempb
    k2_bsb = 2*tempb
    k3_bsb = 2*tempb
    k4_bsb = tempb
    bub(t) = bub(t) + bub(t+1)
    tempb = bub(t+1)/6.
    bub(t+1) = 0.0
    k1_bub = tempb
    k2_bub = 2*tempb
    k3_bub = 2*tempb
    k4_bub = tempb
    dwb = dt*k4_wb
    dbsb = dt*k4_bsb
    dbub = dt*k4_bub
    xtmp_wb = -(ns**2*dbsb) - nu**2*dbub
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xtmp_bub = dwb
      xtmp_bsb = 0.0
    ELSE
      xtmp_bsb = dwb
      xtmp_bub = 0.0
    END IF
    bsb(t) = bsb(t) + xtmp_bsb
    k2_bsb = k2_bsb + 0.5*xtmp_bsb
    bub(t) = bub(t) + xtmp_bub
    k2_bub = k2_bub + 0.5*xtmp_bub
    wb(t) = wb(t) + xtmp_wb
    k2_wb = k2_wb + 0.5*xtmp_wb
    dwb = dt*k3_wb
    dbsb = dt*k3_bsb
    dbub = dt*k3_bub
    xtmp_wb = -(ns**2*dbsb) - nu**2*dbub
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xtmp_bub = dwb
      xtmp_bsb = 0.0
    ELSE
      xtmp_bsb = dwb
      xtmp_bub = 0.0
    END IF
    bsb(t) = bsb(t) + xtmp_bsb
    k2_bsb = k2_bsb + 0.5*xtmp_bsb
    bub(t) = bub(t) + xtmp_bub
    k2_bub = k2_bub + 0.5*xtmp_bub
    wb(t) = wb(t) + xtmp_wb
    k2_wb = k2_wb + 0.5*xtmp_wb
    dwb = dt*k2_wb
    dbsb = dt*k2_bsb
    dbub = dt*k2_bub
    xtmp_wb = -(ns**2*dbsb) - nu**2*dbub
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      xtmp_bub = dwb
      xtmp_bsb = 0.0
    ELSE
      xtmp_bsb = dwb
      xtmp_bub = 0.0
    END IF
    bsb(t) = bsb(t) + xtmp_bsb
    k1_bsb = k1_bsb + 0.5*xtmp_bsb
    bub(t) = bub(t) + xtmp_bub
    k1_bub = k1_bub + 0.5*xtmp_bub
    k1_wb = k1_wb + 0.5*xtmp_wb
    dwb = dt*k1_wb
    dbsb = dt*k1_bsb
    dbub = dt*k1_bub
    wb(t) = wb(t) + xtmp_wb - ns**2*dbsb - nu**2*dbub
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      bub(t) = bub(t) + dwb
    ELSE
      bsb(t) = bsb(t) + dwb
    END IF
  END DO

  write(6, *) 'time adj: ', t+1
  ! write(6,*) 'adj t=1: ', wb(1),bub(1), bsb(1)
  write(16,*) wb(1),bub(1), bsb(1)

  ! bsb(1) = 0.0
  ! bub(1) = 0.0
  ! wb(1) = 0.0
END SUBROUTINE FWD_B