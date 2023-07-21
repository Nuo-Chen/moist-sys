### RK_flux
in:
    - u,v,w,ThetaR,qt, 
    - uhat,vhat,what,ThetaRhat,qthat,  ! u,v,w,ThetaR,qt in fourier space
out:
    - fu,fv,fw,fThetaR,fqt

intermediate:
    - qr

### ad_rk_flux.f90 adjoint
input : 
    basic trajectory:
        - u     (nx,ny,m)
        - v     (nx,ny,m)
        - w     (nx,ny,m)
        - qvini (m)
        - qvs   (m)
        - qt    (nx,ny,m)
        - fqvdz (m)
    
    adjoint forcing:
        - fqtb      (nxph,nyp,m)
        - fThetaRb  (nxph,nyp,m)
        - fwb       (nxph,nyp,m)
        - fvb       (nxph,nyp,m)
        - fub       (nxph,nyp,m)

output:
        - ub    (nx,ny,m) 
        - vb    (nx,ny,m) 
        - wb    (nx,ny,m) 
        - qtb   (nx,ny,m) 
        - ThetaRb   (nx,ny,m) 

        <!-- - uhatb     (nxph,nyp,m) 
        - vhatb     (nxph,nyp,m) 
        - whatb     (nxph,nyp,m) 
        - ThetaRhat (nxph,nyp,m) = 0 ! only in viscosity term
        - qthat     (nxph,nyp,m) = 0 ! only in viscosity term -->



intermediate:
    - wrkNL     (nx,ny,m)
    - qr        (nx,ny,m)

    - wrkNLb    (nx,ny,m)
    - wrkNLhatb (nxph,nyp,m)
    - qrb       (nx,ny,m)    
    
    - wrk1b     (nxph,nyp)
    - out1b     (nxh,ny)
    - in1b      (nx,ny)

### fare.f90
initial data:
    - u     (nx,ny,m) ;0
    - v     (nx,ny,m) ;0
    - w     (nx,ny,m) ;0
    - Theta (nx,ny,m) ;0
    - qr    (nx,ny,m) ;0

out:
    - u, v, w, Theta, qr, qt

intermediate:
    - qvini (m)       ;0
    - qvs   (m)       ;0
    - dqvdz (m)       ;0
    - tau_z (m)
    - qv    (nx,ny,m) ; <- qr
    - qt    (nx,ny,m) ; <- qv+qr
    - ThetaR (nx,ny,m) ; Theta-L*qr

    - uhatb     (nxph,nyp,m) ; <- u 
    - vhatb     (nxph,nyp,m) ; <- v 
    - whatb     (nxph,nyp,m) ; <- w  
    - ThetaRhat (nxph,nyp,m) ; <- ThetaR 
    - qthat     (nxph,nyp,m) ; <- qt

    padm + c2r
    - uhat -> wrk   -> out  ->   in -> u 
    - (nxph,nyp,m) -> (nxh,ny) -> (nx,ny)

    r2c + unpadm
    u -> in   ->    out   -> wrk -> uhat
    - (nx,ny) -> (nxh,ny) -> (nxph,nyp,m)

coefficients:
    - e_nu_1/2/3    (nxph,nyp,m)
    - e_nu_1/2/3uv  (nxph,nyp,m)
    - e_nu_1/2/3w   (nxph,nyp,m)

    - a,b,c     (nxph,nyp,m-1)

RK flux 1:
    - u,v,w,ThetaR,qt, uhat,vhat,what,ThetaRhat,qthat, 
    ->  - fu,fv,fw,fThetaR,fqt

rk loops: 
    rk1: u, uhat -> fuhat
    u1hat = (uhat + dt/3* fuhat) * e_nu_1uv
    rhat <- uhat
    phat <- rhat
    u1hat = u1hat - IM * kx * phat
                    u1hat -> u
    rk2: u, u1hat -> fu1hat
    u1hat = uhat * e_nu_2uv + 2/3* dt * fu1hat * e_nu_1uv;
    rhat <- uhat1
    phat <- rhat
    u1hat = u1hat - IM * kx * phat
                    u1hat -> u
    rk3: u, u1hat -> fu1hat
    uhat = uhat * e_nu_3uv + dt/4 * fuhat * e_nu_3uv + 3/4 * dt * fu1hat * e_nu_1uv;
    rhat <- uhat
    phat <- rhat
    uhat = uhat - IM * kx * phat
                    uhat -> u
    qv = f(qt,qr)
    qt = qv + qr
    Theta = ThetaR + L * qr
    ThetaRhat <- ThetaR
    qthat <- qt
end loop

### ad_fare.f90
input : 
    basic trajectory:
        - u     (nx,ny,m)
        - v     (nx,ny,m)
        - w     (nx,ny,m)
        - qvini (m)
        - qvs   (m)
        - qt    (nx,ny,m)
        - fqvdz (m)
    
    adjoint forcing:
        - a_qt      (nx,ny,m)
        - a_ThetaR  (nx,ny,m)
        - a_u       (nx,ny,m)
        - a_v       (nx,ny,m)
        - a_w       (nx,ny,m)

a_ThetaR -> a_ThetaRhat
a_qt -> a_qthat
a_u -> a_uhat
a_v -> a_vhat
a_w -> a_what
rk loops: 
    a_ThetaRhat -> a_ThetaR
    a_qthat -> a_qt

    a_ThetaR <- f(a_ThetaR, a_Theta)
    a_qr <- f(a_qr, a_Theta, a_qt) 
    a_qv <- f(a_qv, a_qt)
    a_qr <- f(a_qt, a_qr) 

    a_phat <- f(a_what, a_uhat, a_vhat)
    a_rhat <- a_phat
    a_what, a_uhat, a_vhat <- a_rhat
    a_fuhat = a_fuhat + f(a_uhat * 1/4)
    a_fu1hat = a_fuhat + f(a_uhat * 3/4)
    a_uhat = f(a_uhat)
    a_u1hat = 0

    rk1: u, uhat -> fuhat
    u1hat = (uhat + dt/3* fuhat) * e_nu_1uv
    rhat <- uhat
    phat <- rhat
    u1hat = u1hat - IM * kx * phat
                    u1hat -> u
    rk2: u, u1hat -> fu1hat
    u1hat = uhat * e_nu_2uv + 2/3* dt * fu1hat * e_nu_1uv;
    rhat <- uhat1
    phat <- rhat
    u1hat = u1hat - IM * kx * phat
                    u1hat -> u
    rk3: u, u1hat -> fu1hat
    uhat = uhat * e_nu_3uv + dt/4 * fuhat * e_nu_3uv + 3/4 * dt * fu1hat * e_nu_1uv;
    rhat <- uhat
    phat <- rhat
    uhat = uhat - IM * kx * phat
                    uhat -> u
    qv = f(qt,qr)
    qt = qv + qr
    Theta = ThetaR + L * qr
    
end loop

output: