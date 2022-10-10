subroutine fwd(nt, w, bu, bs, nu, ns, dt)
    implicit none
    integer :: nt, t
    real :: dt, nu, ns
    real :: dwf, dbuf, dbsf
    real :: dwpf, dbupf, dbspf
    real :: dwppf, dbuppf, dbsppf
    real, intent(inout) :: w(nt), bu(nt), bs(nt)


    nu = 0.012
    ns = 0.0012
    dt = 0.1

    w(1) = 0.1
    bs(1) = 0
    bu(1) = 0

    ! first step
    t = 1
    if (bs(t) .lt. bu(t)) then
        dwf = bu(t)*dt
    else
        dwf = bs(t)*dt
    end if
    w(t+1) = w(t) + dwf

    dbuf = - nu*nu*w(t)*dt
    dbsf = - ns*ns*w(t)*dt
    bu(t+1) = bu(t) + dbuf 
    bs(t+1) = bs(t) + dbsf 

    dwpf = dwf
    dbspf = dbsf
    dbupf = dbuf

    ! second step
    t = 2
    if (bs(t) .lt. bu(t)) then
        dwf = bu(t)*dt
    else
        dwf = bs(t)*dt
    end if
    w(t+1) = w(t) + dwf*1.5 - dwpf*0.5

    dbuf = - nu*nu*w(t)*dt
    dbsf = - ns*ns*w(t)*dt
    bu(t+1) = bu(t) + dbuf*1.5 - dbupf*0.5
    bs(t+1) = bs(t) + dbsf*1.5 - dbspf*0.5

    dwppf = dwpf
    dbsppf = dbspf
    dbuppf = dbupf
    dwpf = dwf
    dbspf = dbsf
    dbupf = dbuf
   

    ! leap frog          
    do t = 3, nt-1
        if (bs(t) .lt. bu(t)) then
            dwf = bu(t)*dt
        else
            dwf = bs(t)*dt
        end if
        w(t+1) = w(t) + dwf*23./12. - dwpf*16./12. + dwppf*5./12.

        dbuf = - nu*nu*w(t)*dt
        dbsf = - ns*ns*w(t)*dt
        bu(t+1) = bu(t) + dbuf*23./12. - dbupf*16./12. + dbuppf*5./12.
        bs(t+1) = bs(t) + dbsf*23./12. - dbspf*16./12. + dbsppf*5./12.

        dwppf = dwpf
        dbsppf = dbspf
        dbuppf = dbupf
        dwpf = dwf
        dbspf = dbsf
        dbupf = dbuf
    end do

end subroutine fwd