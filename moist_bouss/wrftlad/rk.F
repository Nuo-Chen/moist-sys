SUBROUTINE rk_update_scalar( scs, sce,                      &
    scalar_1, scalar_2, sc_tend,   &
    advh_t, advz_t,                &
    advect_tend,                   &
    h_tendency, z_tendency,        &
    msftx, msfty, c1, c2,          &
    mu_old, mu_new, mu_base,       &
    rk_step, dt, spec_zone,        &
    config_flags,                  &
    tenddec,                       &
    ids, ide, jds, jde, kds, kde,  &
    ims, ime, jms, jme, kms, kme,  &
    its, ite, jts, jte, kts, kte  )

      IMPLICIT NONE

      !  Input data.

      TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags
      LOGICAL ,                INTENT(IN   ) :: tenddec

      INTEGER ,                INTENT(IN   ) :: scs, sce, rk_step, spec_zone
      INTEGER ,                INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                          ims, ime, jms, jme, kms, kme, &
                          its, ite, jts, jte, kts, kte

      REAL,                    INTENT(IN   ) :: dt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce),                &
      INTENT(INOUT)                                  :: scalar_1,   scalar_2

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce),                &
      INTENT(IN)                                     :: sc_tend

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
      INTENT(IN)                                  :: advect_tend

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme ), INTENT(INOUT) , OPTIONAL :: advh_t,  advz_t ! accumulating for output
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: h_tendency, z_tendency ! from rk_scalar_tend

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN   ) ::  mu_old,  &
                                      mu_new,  &
                                      mu_base, &
                                      msftx,   &
                                      msfty

      REAL, DIMENSION(kms:kme ), INTENT(IN   ) ::  c1, c2

      INTEGER :: i,j,k,im
      REAL    :: sc_middle, msfsq
      REAL, DIMENSION(its:ite, kts:kte, jts:jte  ) :: tendency

      INTEGER :: i_start,i_end,j_start,j_end,k_start,k_end
      INTEGER :: i_start_spc,i_end_spc,j_start_spc,j_end_spc,k_start_spc,k_end_spc

      REAL, DIMENSION(its:ite) :: muold, munew












      !<DESCRIPTION>
      !
      !  rk_scalar_update advances the scalar equation given the time t value
      !  of the scalar and the scalar tendency.  
      !
      !</DESCRIPTION>


      !
      !  set loop limits.

      i_start = its
      i_end   = min(ite,ide-1)
      j_start = jts
      j_end   = min(jte,jde-1)
      k_start = kts
      k_end   = kte-1

      i_start_spc = i_start
      i_end_spc   = i_end
      j_start_spc = j_start
      j_end_spc   = j_end
      k_start_spc = k_start
      k_end_spc   = k_end

      IF( config_flags%nested .or. config_flags%specified ) THEN
      IF( .NOT. config_flags%periodic_x)i_start = max( its,ids+spec_zone )
      IF( .NOT. config_flags%periodic_x)i_end   = min( ite,ide-spec_zone-1 )
      j_start = max( jts,jds+spec_zone )
      j_end   = min( jte,jde-spec_zone-1 )
      k_start = kts
      k_end   = min( kte, kde-1 )
      ENDIF

      IF ( rk_step == 1 ) THEN

      !  replace t-dt values (in scalar_1) with t values scalar_2,
      !  then compute new values by adding tendency to values at t

      DO  im = scs,sce

        !
        DO  j = jts, min(jte,jde-1)
        DO  k = kts, min(kte,kde-1)
        DO  i = its, min(ite,ide-1)
        tendency(i,k,j) = 0.
        ENDDO
        ENDDO
        ENDDO

        DO  j = j_start,j_end
        DO  k = k_start,k_end
        DO  i = i_start,i_end
        ! scalar was coupled with my
        tendency(i,k,j) = advect_tend(i,k,j) * msfty(i,j)
        ENDDO
        ENDDO
        ENDDO

        DO  j = j_start_spc,j_end_spc
        DO  k = k_start_spc,k_end_spc
        DO  i = i_start_spc,i_end_spc
        tendency(i,k,j) = tendency(i,k,j) + sc_tend(i,k,j,im)
        ENDDO
        ENDDO
        ENDDO

        DO  j = jts, min(jte,jde-1)

        !
          DO  i = its, min(ite,ide-1)
          MUOLD(i) = MU_OLD(i,j) + MU_BASE(i,j)
          MUNEW(i) = MU_NEW(i,j) + MU_BASE(i,j)
          ENDDO

          DO  k = kts, min(kte,kde-1)
            DO  i = its, min(ite,ide-1)

            scalar_1(i,k,j,im) = scalar_2(i,k,j,im)
            scalar_2(i,k,j,im) = ((c1(k)*muold(i)+c2(k))*scalar_1(i,k,j,im)   &
                + dt*tendency(i,k,j))/(c1(k)*munew(i)+c2(k))

            ENDDO !i
          ENDDO !k
        ENDDO !j

      ENDDO !im

      ELSE ! rk != 1
      !  just compute new values, scalar_1 already at time t.

      DO  im = scs, sce

      !
        DO  j = jts, min(jte,jde-1)
        DO  k = kts, min(kte,kde-1)
        DO  i = its, min(ite,ide-1)
        tendency(i,k,j) = 0.
        ENDDO
        ENDDO
        ENDDO

        DO  j = j_start,j_end
        DO  k = k_start,k_end
        DO  i = i_start,i_end
        ! scalar was coupled with my
        tendency(i,k,j) = advect_tend(i,k,j) * msfty(i,j)
        ENDDO
        ENDDO
        ENDDO

        DO  j = j_start_spc,j_end_spc
        DO  k = k_start_spc,k_end_spc
        DO  i = i_start_spc,i_end_spc
        tendency(i,k,j) = tendency(i,k,j) + sc_tend(i,k,j,im)
        ENDDO
        ENDDO
        ENDDO

        DO  j = jts, min(jte,jde-1)

        !
          DO  i = its, min(ite,ide-1)
            MUOLD(i) = MU_OLD(i,j) + MU_BASE(i,j)
            MUNEW(i) = MU_NEW(i,j) + MU_BASE(i,j)
          ENDDO

          DO  k = kts, min(kte,kde-1)
          DO  i = its, min(ite,ide-1)

            scalar_2(i,k,j,im) = ((c1(k)*muold(i)+c2(k))*scalar_1(i,k,j,im)   &
                + dt*tendency(i,k,j))/(c1(k)*munew(i)+c2(k))

          ENDDO !i
          ENDDO !k

          ! This is separated from the k/i-loop above for better performance
          IF ( PRESENT(advh_t) .AND. PRESENT(advz_t) ) THEN
            IF(tenddec.and.rk_step.eq.config_flags%rk_ord) THEN
              DO  k = kts, min(kte,kde-1)
              DO  i = its, min(ite,ide-1)

              advh_t(i,k,j) = advh_t(i,k,j) + (dt*h_tendency(i,k,j)* msfty(i,j))/(c1(k)*munew(i)+c2(k))
              advz_t(i,k,j) = advz_t(i,k,j) + (dt*z_tendency(i,k,j)* msfty(i,j))/(c1(k)*munew(i)+c2(k))

              ENDDO !i
              ENDDO !k
            END IF
          END IF
      ENDDO !j

      ENDDO !im

      END IF

END SUBROUTINE rk_update_scalar 
















SUBROUTINE rk_tendency ( config_flags, rk_step,        &
      ru_tend, rv_tend, rw_tend, ph_tend, t_tend,      &
      ru_tendf, rv_tendf, rw_tendf, ph_tendf, t_tendf, &
      mu_tend, u_save, v_save, w_save, ph_save,        &
      t_save, mu_save, RTHFTEN,                        &
      ru, rv, rw, ww, wwE, wwI,                        &
      u, v, w, t, ph,                                  &
      u_old, v_old, w_old, t_old, ph_old,              &
      h_diabatic, phb,t_init,                          &
      mu_old, mu, mut, muu, muv, mub,                  &
      c1h, c2h, c1f, c2f,                              &
      al, ht, alt, p, pb, php, cqu, cqv, cqw,          &
      u_base, v_base, t_base, qv_base, z_base,         &
      msfux, msfuy, msfvx, msfvx_inv,                  &
      msfvy, msftx, msfty,                             &
      clat, f, e, sina, cosa,                          &
      fnm, fnp, rdn, rdnw,                             &
      dt, rdx, rdy, khdif, kvdif, xkmhd, xkhh,         &
      diff_6th_opt, diff_6th_factor,                   &
      adv_opt,                                         &
      dampcoef,zdamp,damp_opt,rad_nudge,               &
      cf1, cf2, cf3, cfn, cfn1, n_moist,               &
      non_hydrostatic, top_lid,                        &
      u_frame, v_frame,                                &
      ids, ide, jds, jde, kds, kde,                    &
      ims, ime, jms, jme, kms, kme,                    &
      its, ite, jts, jte, kts, kte,                    &
      max_vert_cfl, max_horiz_cfl)

IMPLICIT NONE

!  Input data.

TYPE(grid_config_rec_type)    ,           INTENT(IN   ) :: config_flags

INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte

LOGICAL ,               INTENT(IN   ) :: non_hydrostatic, top_lid

INTEGER ,               INTENT(IN   ) :: n_moist, rk_step

REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,              &
                   INTENT(IN   ) :: ru,      &
                                    rv,      &
                                    rw,      &
                                    ww,      &
                                    u,       &
                                    v,       &
                                    w,       &
                                    t,       &
                                    ph,      &
                                    u_old,   &
                                    v_old,   &
                                    w_old,   &
                                    t_old,   &
                                    ph_old,  &
                                    phb,     &
                                    al,      &
                                    alt,     &
                                    p,       &
                                    pb,      &
                                    php,     &
                                    cqu,     &
                                    cqv,     &
                                    t_init,  &
                                    xkmhd,   &
                                    xkhh,    &
                                    h_diabatic

REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,              &
                   INTENT(OUT  ) :: ru_tend, &
                                    rv_tend, &
                                    rw_tend, &
                                    t_tend,  &
                                    ph_tend, &
                                    RTHFTEN, &
                                     u_save, &
                                     v_save, &
                                     w_save, &
                                    ph_save, &
                                     t_save

REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,               &
                   INTENT(INOUT) :: ru_tendf, &
                                    rv_tendf, &
                                    rw_tendf, &
                                    t_tendf,  &
                                    ph_tendf, &
                                    cqw

REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(  OUT) :: mu_tend, &
                                               mu_save

REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfux,   &
                                               msfuy,   &
                                               msfvx,   &
                                               msfvx_inv,   &
                                               msfvy,   &
                                               msftx,   &
                                               msfty,   &
                                               clat,    &
                                               f,       &
                                               e,       &
                                               sina,    &
                                               cosa,    &
                                               mu_old,  &
                                               mu,      &
                                               mut,     &
                                               mub,     &
                                               muu,     &
                                               muv,     &
                                               ht             

REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fnm,     &
                                             fnp,     &
                                             rdn,     &
                                             rdnw,    &
                                             u_base,  &
                                             v_base,  &
                                             t_base,  &
                                             qv_base, &
                                             z_base,  &
                                             c1h,     &
                                             c2h,     &
                                             c1f,     &
                                             c2f

REAL ,                                      INTENT(IN   ) :: rdx,     &
                                           rdy,     &
                                           dt,      &
                                           u_frame, &
                                           v_frame, &
                                           khdif,   &
                                           kvdif
INTEGER, INTENT( IN ) :: diff_6th_opt
REAL,    INTENT( IN ) :: diff_6th_factor
INTEGER, INTENT( IN ) :: adv_opt

INTEGER, INTENT( IN ) :: damp_opt, rad_nudge

REAL, INTENT( IN ) :: zdamp, dampcoef

REAL, INTENT( OUT ) :: max_horiz_cfl
REAL, INTENT( OUT ) :: max_vert_cfl

REAL    :: kdift, khdq, kvdq, cfn, cfn1, cf1, cf2, cf3
INTEGER :: i,j,k
INTEGER :: time_step
INTEGER :: rk_order
REAL    :: dt_step

! IEVA declarations
REAL , DIMENSION( ims:ime, kms:kme, jms:jme ) :: wwE, wwI  
REAL , DIMENSION( ims:ime , jms:jme )         :: mut_old, mut_new                                                              
REAL , DIMENSION( ims:ime , jms:jme )         :: muu_old, muu_new                                                              
REAL , DIMENSION( ims:ime , jms:jme )         :: muv_old, muv_new                                                              
LOGICAL                                       :: ieva

!<DESCRIPTION>
!
!  rk_tendency computes the large-timestep tendency terms in the
!  momentum, thermodynamic (theta), and geopotential equations.  
!  These terms include:
!
!  (1) advection (for u, v, w, theta - calls to advect_u, advect_v,
!                 advect_w, and advact_scalar).
!
!  (2) geopotential equation terms (advection and "gw" - call to rhs_ph).
!
!  (3) buoyancy term in vertical momentum equation (call to pg_buoy_w).
!
!  (4) Coriolis and curvature terms in u,v,w momentum equations
!      (calls to subroutines coriolis, curvature)
!
!  (5) 3D diffusion on coordinate surfaces.
!
!</DESCRIPTION>

CALL zero_tend ( ru_tend,                      &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend ( rv_tend,                      &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend ( rw_tend,                      &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend ( t_tend,                       &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend ( ph_tend,                      &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend ( u_save,                       &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend ( v_save,                       &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend ( w_save,                       &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend ( ph_save,                       &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend ( t_save,                       &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )

CALL zero_tend2d( mu_tend,                  &
ids, ide, jds, jde, 1, 1, &
ims, ime, jms, jme, 1, 1, &
its, ite, jts, jte, 1, 1 )

CALL zero_tend2d( mu_save,                  &
ids, ide, jds, jde, 1, 1, &
ims, ime, jms, jme, 1, 1, &
its, ite, jts, jte, 1, 1 )

CALL nl_get_time_step ( 1, time_step )

rk_order = config_flags%rk_ord
dt_step  = dt / (rk_order - rk_step + 1)   ! needed for calculations using sub-rk step

   
! Code for splitting vertical velocity into ex/im parts

CALL WW_SPLIT(wwE, wwI,                         &
u,  v,  ww,                       &
mut, rdnw, msfty,                 &
c1f, c2f,                         &
rdx, rdy, msfux, msfuy,           &
msfvx, msfvy, dt,                 &
config_flags, rk_step,            &
ids, ide, jds, jde, kds, kde,     &
ims, ime, jms, jme, kms, kme,     &
its, ite, jts, jte, kts, kte )              

! IEVA ON/OFF

ieva = CHK_IEVA( config_flags, rk_step ) 

IF( ieva ) THEN

! Need an estimate of the mut at the original ('n') time level...

CALL calculate_full( mut_old,  mub, mu_old,     &
   ids, ide, jds, jde, 1, 2,  &
   ims, ime, jms, jme, 1, 1,  &
   its, ite, jts, jte, 1, 1 )

! Estimate new temporary column mass
   
CALL calc_mut_new ( u, v, c1h, c2h,                      &
  mut_old, muu, muv, mut_new,          &
  dt, rdx, rdy, msftx, msfty,          &
  msfux, msfuy, msfvx, msfvx_inv,      &
  msfvy, rdnw,                         &
  ids, ide, jds, jde, kds, kde,        &
  ims, ime, jms, jme, kms, kme,        &
  its, ite, jts, jte, kts, kte    )
    
! Create staggered mass values for u and v points

CALL calc_mu_uv_1 ( config_flags,                 &
  mut_old, muu_old, muv_old,    &
  ids, ide, jds, jde, kds, kde, &
  ims, ime, jms, jme, kms, kme, &
  its, ite, jts, jte, kts, kte )
 
CALL calc_mu_uv_1 ( config_flags,                 &
  mut_new, muu_new, muv_new,    &
  ids, ide, jds, jde, kds, kde, &
  ims, ime, jms, jme, kms, kme, &
  its, ite, jts, jte, kts, kte )

ENDIF

! Okay do normal advection now....using the wwE array

IF( (rk_step == rk_order) .and. ( adv_opt == WENO_MOM ) ) THEN

CALL advect_weno_u ( u, u , ru_tend, ru, rv, wwE,  &
     c1h, c2h,                     &
     mut, time_step, config_flags, &
     msfux, msfuy, msfvx, msfvy,   &
     msftx, msfty,                 &
     fnm, fnp, rdx, rdy, rdnw,     &
     ids, ide, jds, jde, kds, kde, &
     ims, ime, jms, jme, kms, kme, &
     its, ite, jts, jte, kts, kte )

ELSE

CALL advect_u ( u, u , ru_tend, ru, rv, wwE,  &
c1h, c2h,                     &
mut, time_step, config_flags, &
msfux, msfuy, msfvx, msfvy,   &
msftx, msfty,                 &
fnm, fnp, rdx, rdy, rdnw,     &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )
ENDIF

IF( ieva ) THEN

CALL advect_u_implicit ( u, u_old, ru_tend, ru, rv, wwI,  &
         c1h, c2h,                        &
         muu_old, muu, muu_new,           &
         config_flags,                    &
         msfux, msfuy, msfvx, msfvy,      &
         msftx, msfty,                    &
         fnm, fnp,                        &
         dt_step,                         &
         rdx, rdy, rdnw,                  &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         its, ite, jts, jte, kts, kte )
ENDIF

IF( (rk_step == rk_order) .and. ( adv_opt == WENO_MOM ) ) THEN

CALL advect_weno_v ( v, v , rv_tend, ru, rv, wwE,  &
      c1h, c2h,                     &
      mut, time_step, config_flags, &
      msfux, msfuy, msfvx, msfvy,   &
      msftx, msfty,                 &
      fnm, fnp, rdx, rdy, rdnw,     &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte )

ELSE

CALL advect_v ( v, v , rv_tend, ru, rv, wwE,  &
c1h, c2h,                     &
mut, time_step, config_flags, &
msfux, msfuy, msfvx, msfvy,   &
msftx, msfty,                 &
fnm, fnp, rdx, rdy, rdnw,     &
ids, ide, jds, jde, kds, kde, &
ims, ime, jms, jme, kms, kme, &
its, ite, jts, jte, kts, kte )
ENDIF

IF( ieva ) THEN

CALL advect_v_implicit ( v, v_old, rv_tend, ru, rv, wwI,  &
         c1h, c2h,                        &
         muv_old, muv, muv_new,           &
         config_flags,                    &
         msfux, msfuy, msfvx, msfvy,      &
         msftx, msfty,                    &
         fnm, fnp,                        &
         dt_step,                         &
         rdx, rdy, rdnw,                  &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         its, ite, jts, jte, kts, kte )
ENDIF

IF (non_hydrostatic) THEN

IF( (rk_step == rk_order) .and. ( adv_opt == WENO_MOM ) ) THEN

CALL advect_weno_w ( w, w, rw_tend, ru, rv, wwE,   &
       c1h, c2h,                     &
       mut, time_step, config_flags, &
       msfux, msfuy, msfvx, msfvy,   &
       msftx, msfty,                 &
       fnm, fnp, rdx, rdy, rdn,      &
       ids, ide, jds, jde, kds, kde, &
       ims, ime, jms, jme, kms, kme, &
       its, ite, jts, jte, kts, kte )

ELSE

CALL advect_w ( w, w, rw_tend, ru, rv, wwE,   &
  c1h, c2h,                     &
  mut, time_step, config_flags, &
  msfux, msfuy, msfvx, msfvy,   &
  msftx, msfty,                 &
  fnm, fnp, rdx, rdy, rdn,      &
  ids, ide, jds, jde, kds, kde, &
  ims, ime, jms, jme, kms, kme, &
  its, ite, jts, jte, kts, kte )
ENDIF

ENDIF  ! non-hydrostatic

!  theta flux divergence

! 11/2016 ERM: Use WENO for theta flux on 3rd RK step if using WENO_SCALAR or WENOPD_SCALAR
! to be consistent with other scalar fluxes
IF(  ( config_flags%scalar_adv_opt == WENO_SCALAR      &
.or. config_flags%scalar_adv_opt == WENOPD_SCALAR    &
.or. config_flags%moist_adv_opt == WENO_SCALAR       &
.or. config_flags%moist_adv_opt == WENOPD_SCALAR )   &
.and. (rk_step == rk_order) ) THEN

! also use weno for monotonic scalar option so that the h_ and z_tendency arrays are not needed

CALL advect_scalar_weno ( t, t, t_tend, ru, rv, wwE,     &
          c1h, c2h, mut, time_step,      &
          config_flags,                  &
          msfux, msfuy, msfvx, msfvy,    &
          msftx, msfty, fnm, fnp,        &
          rdx, rdy, rdnw,                &
          ids, ide, jds, jde, kds, kde,  &
          ims, ime, jms, jme, kms, kme,  &
          its, ite, jts, jte, kts, kte  )
ELSE

CALL advect_scalar ( t, t, t_tend, ru, rv, wwE,    &
     c1h, c2h,                     &
     mut, time_step, config_flags, &
     msfux, msfuy, msfvx, msfvy,   &
     msftx, msfty, fnm, fnp,       &
     rdx, rdy, rdnw,               &
     ids, ide, jds, jde, kds, kde, &
     ims, ime, jms, jme, kms, kme, &
     its, ite, jts, jte, kts, kte )

ENDIF

IF( ieva ) THEN

CALL advect_s_implicit ( t, t_old, t_tend, ru, rv, wwI,  &
         c1h, c2h,                       &
         mut_old, mut, mut_new,          &
         config_flags,                   &
         msfux, msfuy, msfvx, msfvy,     &
         msftx, msfty,                   &
         fnm, fnp,                       &
         dt_step,                        &
         rdx, rdy, rdnw,                 &
         ids, ide, jds, jde, kds, kde,   &
         ims, ime, jms, jme, kms, kme,   &
         its, ite, jts, jte, kts, kte )

ENDIF

IF ( config_flags%cu_physics == GDSCHEME  .OR.     &
config_flags%cu_physics == GFSCHEME  .OR.     &
config_flags%cu_physics == G3SCHEME  .OR.     &
config_flags%cu_physics == NTIEDTKESCHEME )  THEN     ! NTiedtke

! theta advection only:

CALL set_tend( RTHFTEN, t_tend, msfty,          &
ids, ide, jds, jde, kds, kde,    &
ims, ime, jms, jme, kms, kme,    &
its, ite, jts, jte, kts, kte     )

END IF

CALL rhs_ph( ph_tend, u, v, wwE, ph, ph, phb, w, &
mut, muu, muv,                     &
c1f, c2f,                          &
fnm, fnp,                          &
rdnw, cfn, cfn1, rdx, rdy,         &
msfux, msfuy, msfvx,               &
msfvx_inv, msfvy,                  &
msftx, msfty,                      &
non_hydrostatic,                   &
config_flags,                      &
ids, ide, jds, jde, kds, kde,      &
ims, ime, jms, jme, kms, kme,      &
its, ite, jts, jte, kts, kte      )

IF( ieva ) THEN

CALL advect_ph_implicit ( ph, ph_old, ph_tend, phb,       &
             ru, rv, wwE, wwI, w,            & 
             c1f, c2f,                       &
             mut, config_flags,              &
             msfux, msfuy, msfvx, msfvy,     &
             msftx, msfty,                   &
             fnm, fnp,                       &
             dt_step,                        &
             rdx, rdy, rdnw,                 &
             ids, ide, jds, jde, kds, kde,   &
             ims, ime, jms, jme, kms, kme,   &
             its, ite, jts, jte, kts, kte )

IF (non_hydrostatic) THEN

CALL advect_w_implicit ( w, w_old, rw_tend,              &
           ru_tend, rv_tend, ht, wwI,      &
           ph, ph_old, ph_tend,            &
           c1f, c2f, cf1, cf2, cf3,        &
           mut_old, mut, mut_new,          &
           config_flags,                   &
           msfux, msfuy, msfvx, msfvy,     &
           msftx, msfty,                   &
           fnm, fnp,                       &
           dt_step,                        &
           rdx, rdy, rdn,                  &
           ids, ide, jds, jde, kds, kde,   &
           ims, ime, jms, jme, kms, kme,   &
           its, ite, jts, jte, kts, kte )
ENDIF

ENDIF  ! IEVA

CALL horizontal_pressure_gradient( ru_tend,rv_tend,                 &
                  ph,alt,p,pb,al,php,cqu,cqv,     &
                  muu,muv,mu,c1h,c2h,fnm,fnp,rdnw,&
                  cf1,cf2,cf3,cfn,cfn1,           &
                  rdx,rdy,msfux,msfuy,            &
                  msfvx,msfvy,msftx,msfty,        &
                  config_flags, non_hydrostatic,  &
                  top_lid,                        &
                  ids, ide, jds, jde, kds, kde,   &
                  ims, ime, jms, jme, kms, kme,   &
                  its, ite, jts, jte, kts, kte   )

IF (non_hydrostatic) THEN
CALL pg_buoy_w( rw_tend, p, cqw, mu, mub,       &
   c1f,c2f,                        &
   rdnw, rdn, g, msftx, msfty,     &
   ids, ide, jds, jde, kds, kde,   &
   ims, ime, jms, jme, kms, kme,   &
   its, ite, jts, jte, kts, kte   )
ENDIF

CALL w_damp   ( rw_tend, max_vert_cfl,             &
max_horiz_cfl,                    &
u, v, ww, w, mut, c1f, c2f, rdnw, &
rdx, rdy, msfux, msfuy, msfvx,    &
msfvy, dt, config_flags,          &
ids, ide, jds, jde, kds, kde,     &
ims, ime, jms, jme, kms, kme,     &
its, ite, jts, jte, kts, kte     )

IF(config_flags%pert_coriolis) THEN

CALL perturbation_coriolis ( ru, rv, rw,                   &
                ru_tend,  rv_tend,  rw_tend,  &
                config_flags,                 &
                u_base, v_base, z_base,       &
                muu, muv, c1h, c2h, phb, ph,  &
                msftx, msfty, msfux, msfuy,   &
                msfvx, msfvy,                 &
                f, e, sina, cosa, fnm, fnp,   &
                ids, ide, jds, jde, kds, kde, &
                ims, ime, jms, jme, kms, kme, &
                its, ite, jts, jte, kts, kte )
ELSE
CALL coriolis ( ru, rv, rw,                   &
   ru_tend,  rv_tend,  rw_tend,  &
   config_flags,                 &
   msftx, msfty, msfux, msfuy,   &
   msfvx, msfvy,                 &
   f, e, sina, cosa, fnm, fnp,   &
   ids, ide, jds, jde, kds, kde, &
   ims, ime, jms, jme, kms, kme, &
   its, ite, jts, jte, kts, kte )

END IF

CALL curvature ( ru, rv, rw, u, v, w,            &
ru_tend,  rv_tend,  rw_tend,    &
config_flags,                   &
msfux, msfuy, msfvx, msfvy,     &
msftx, msfty,                   &
clat, fnm, fnp, rdx, rdy,       &
ids, ide, jds, jde, kds, kde,   &
ims, ime, jms, jme, kms, kme,   &
its, ite, jts, jte, kts, kte   )

! Damping option added for Held-Suarez test (also uses lw option HELDSUAREZ)

IF (config_flags%ra_lw_physics == HELDSUAREZ) THEN
CALL held_suarez_damp ( ru_tend, rv_tend,               &   
          ru,rv,p,pb,                     &
          ids, ide, jds, jde, kds, kde,   &
          ims, ime, jms, jme, kms, kme,   &
          its, ite, jts, jte, kts, kte   )
END IF

!**************************************************************
!
!  Next, the terms that we integrate only with forward-in-time
!  (evaluate with time t variables).
!
!**************************************************************

forward_step: IF( rk_step == 1 ) THEN

diff_opt1 : IF (config_flags%diff_opt .eq. 1) THEN

CALL horizontal_diffusion ('u', u, ru_tendf, mut,               &
                   c1h, c2h, config_flags,         &
                   msfux, msfuy, msfvx, msfvx_inv, &
                   msfvy,msftx, msfty,             &
                   khdif, xkmhd, rdx, rdy,         &
                   ids, ide, jds, jde, kds, kde,   &
                   ims, ime, jms, jme, kms, kme,   &
                   its, ite, jts, jte, kts, kte   )

CALL horizontal_diffusion ('v', v, rv_tendf, mut,               &
                   c1h, c2h, config_flags,         &
                   msfux, msfuy, msfvx, msfvx_inv, &
                   msfvy,msftx, msfty,             &
                   khdif, xkmhd, rdx, rdy,         &
                   ids, ide, jds, jde, kds, kde,   &
                   ims, ime, jms, jme, kms, kme,   &
                   its, ite, jts, jte, kts, kte   )

CALL horizontal_diffusion ('w', w, rw_tendf, mut,               &
                   c1f, c2f, config_flags,         &
                   msfux, msfuy, msfvx, msfvx_inv, &
                   msfvy,msftx, msfty,             &
                   khdif, xkmhd, rdx, rdy,         &
                   ids, ide, jds, jde, kds, kde,   &
                   ims, ime, jms, jme, kms, kme,   &
                   its, ite, jts, jte, kts, kte   )

khdq = 3.*khdif
CALL horizontal_diffusion_3dmp ( 'm', t, t_tendf, mut,            &
                    c1h, c2h,                        &
                    config_flags, t_init,            &
                    msfux, msfuy, msfvx, msfvx_inv,  &
                    msfvy, msftx, msfty,             &
                    khdq , xkhh, rdx, rdy,           &
                    ids, ide, jds, jde, kds, kde,    &
                    ims, ime, jms, jme, kms, kme,    &
                    its, ite, jts, jte, kts, kte    )

pbl_test : IF (config_flags%bl_pbl_physics .eq. 0) THEN

CALL vertical_diffusion_u ( u, ru_tendf, config_flags,      &
                 u_base, c1h, c2h,               &
                 alt, muu, rdn, rdnw, kvdif,     &
                 ids, ide, jds, jde, kds, kde,   &
                 ims, ime, jms, jme, kms, kme,   &
                 its, ite, jts, jte, kts, kte   )

CALL vertical_diffusion_v ( v, rv_tendf, config_flags,      &
                 v_base, c1h, c2h,               &
                 alt, muv, rdn, rdnw, kvdif,     &
                 ids, ide, jds, jde, kds, kde,   &
                 ims, ime, jms, jme, kms, kme,   &
                 its, ite, jts, jte, kts, kte   )

IF (non_hydrostatic)                                           &
CALL vertical_diffusion ( 'w', w, rw_tendf, config_flags,      &
               c1f, c2f,                            &
               alt, mut, rdn, rdnw, kvdif,          &
               ids, ide, jds, jde, kds, kde,        &
               ims, ime, jms, jme, kms, kme,        &
               its, ite, jts, jte, kts, kte        )

kvdq = 3.*kvdif
CALL vertical_diffusion_3dmp ( t, t_tendf, config_flags, t_init,     &
                    c1h, c2h,                             &
                    alt, mut, rdn, rdnw, kvdq ,           &
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte         )

ENDIF pbl_test

!  Theta tendency computations.

END IF diff_opt1

IF ( diff_6th_opt .NE. 0 ) THEN

      CALL sixth_order_diffusion( 'u', u, ru_tendf, mut, dt,          &
                        config_flags, c1h, c2h,        &
                        diff_6th_opt, diff_6th_factor, &
                        phb, ph,                       &
                        rdx, rdy,                      &
                        msftx, msfty,                  &
                        msfux, msfuy,                  &
                        msfvx, msfvy,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte )

      CALL sixth_order_diffusion( 'v', v, rv_tendf, mut, dt,          &
                        config_flags, c1h, c2h,        &
                        diff_6th_opt, diff_6th_factor, &
                        phb, ph,                       &
                        rdx, rdy,                      &
                        msftx, msfty,                  &
                        msfux, msfuy,                  &
                        msfvx, msfvy,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte )

      IF (non_hydrostatic)                                            &
        CALL sixth_order_diffusion( 'w', w, rw_tendf, mut, dt,          &
                          config_flags, c1f, c2f,        &
                          diff_6th_opt, diff_6th_factor, &
                          phb, ph,                       &
                          rdx, rdy,                      &
                          msftx, msfty,                  &
                          msfux, msfuy,                  &
                          msfvx, msfvy,                  &
                          ids, ide, jds, jde, kds, kde,  &
                          ims, ime, jms, jme, kms, kme,  &
                          its, ite, jts, jte, kts, kte )

      CALL sixth_order_diffusion( 'm', t,  t_tendf, mut, dt,          &
                        config_flags, c1h, c2h,        &
                        diff_6th_opt, diff_6th_factor, &
                        phb, ph,                       &
                        rdx, rdy,                      &
                        msftx, msfty,                  &
                        msfux, msfuy,                  &
                        msfvx, msfvy,                  &
                        ids, ide, jds, jde, kds, kde,  &
                        ims, ime, jms, jme, kms, kme,  &
                        its, ite, jts, jte, kts, kte )

ENDIF

IF( damp_opt .eq. 2 )                                      &
CALL rk_rayleigh_damp( ru_tendf, rv_tendf,              &
         rw_tendf, t_tendf,               &
         u, v, w, t, t_init,              &
         c1h, c2h, c1f, c2f,              &
         mut, muu, muv, ph, phb,          &
         u_base, v_base, t_base, z_base,  &
         dampcoef, zdamp,                 &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         its, ite, jts, jte, kts, kte   )

IF( rad_nudge .eq. 1 )                                     &
CALL theta_relaxation( t_tendf, t, t_init,              &
         mut, c1h, c2h, ph, phb,          &
         t_base, z_base,                  &
         ids, ide, jds, jde, kds, kde,    &
         ims, ime, jms, jme, kms, kme,    &
         its, ite, jts, jte, kts, kte   )

END IF forward_step

END SUBROUTINE rk_tendency