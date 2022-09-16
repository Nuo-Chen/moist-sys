subroutine fwd(nt, w, bu, bs, nu, ns, dt)
    implicit none
    integer :: nt, t
    real :: dt, nu, ns
    real :: dbu, dbs, dw 
    real :: k1_bu, k1_bs, k1_w 
    real :: k2_bu, k2_bs, k2_w 
    real :: k3_bu, k3_bs, k3_w 
    real :: k4_bu, k4_bs, k4_w 
    real :: xtmp_bu, xtmp_bs, xtmp_w

    real, intent(inout) :: w(nt), bu(nt), bs(nt)

    ! rhs 
    nu = 0.01
    ns = 0.005
    dt = 0.1

    w(1) = 0.1
    bu(1) = 5
    bs(1) = 1

    do t = 1, nt-1
        ! rhs 
        if (bs(t) .lt. bu(t)) then
            dw = bu(t)
        else
            dw = bs(t)
        end if

        dbu = - nu*nu * w(t)
        dbs = - ns*ns * w(t)
        
        ! k1=self.rhs(self.xvar)*self.dt
        k1_bu = dbu * dt
        k1_bs = dbs * dt
        k1_w = dw * dt

        ! xtmp=self.xvar+0.5*k1
        xtmp_w = w(t) + 0.5 * k1_w     
        xtmp_bu = bu(t) + 0.5 * k1_bu
        xtmp_bs = bs(t) + 0.5 * k1_bs
        
        ! rhs 
        if (xtmp_bs .lt. xtmp_bu) then
            dw = xtmp_bu
        else
            dw = xtmp_bs
        end if

        dbu = - nu*nu * xtmp_w
        dbs = - ns*ns * xtmp_w
        
        ! k2=self.rhs(xtmp)*self.dt
        k2_bu = dbu * dt
        k2_bs = dbs * dt
        k2_w = dw * dt

        ! xtmp=self.xvar+0.5*k2
        xtmp_w = w(t) + 0.5 * k2_w     
        xtmp_bu = bu(t) + 0.5 * k2_bu
        xtmp_bs = bs(t) + 0.5 * k2_bs

        ! rhs 
        if (xtmp_bs .lt. xtmp_bu) then
            dw = xtmp_bu
        else
            dw = xtmp_bs
        end if

        dbu = - nu*nu * xtmp_w
        dbs = - ns*ns * xtmp_w
        
        ! k3=self.rhs(xtmp)*self.dt
        k3_bu = dbu * dt
        k3_bs = dbs * dt
        k3_w = dw * dt

        ! xtmp=self.xvar+k2
        xtmp_w = w(t) + 0.5 * k2_w     
        xtmp_bu = bu(t) + 0.5 * k2_bu
        xtmp_bs = bs(t) + 0.5 * k2_bs
        
        ! rhs 
        if (xtmp_bs .lt. xtmp_bu) then
            dw = xtmp_bu
        else
            dw = xtmp_bs
        end if

        dbu = - nu*nu * xtmp_w
        dbs = - ns*ns * xtmp_w
        
        ! k4=self.rhs(xtmp)*self.dt
        k4_bu = dbu * dt
        k4_bs = dbs * dt
        k4_w = dw * dt
        
        ! self.xvar+=(k1+2*k2+2*k3+k4)/6.
        bu(t+1) = bu(t) + (k1_bu + 2*k2_bu + 2*k3_bu + k4_bu)/6.
        bs(t+1) = bs(t) + (k1_bs + 2*k2_bs + 2*k3_bs + k4_bs)/6.
        w(t+1) = w(t) + (k1_w + 2*k2_w + 2*k3_w + k4_w)/6.
    end do

end subroutine fwd