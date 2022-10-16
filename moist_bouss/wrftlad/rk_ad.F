SUBROUTINE a_rk_update_scalar ( scs, sce,                      &
      scalar_1, a_scalar_1, scalar_2, a_scalar_2, sc_tend, a_sc_tend,  &
      advh_t, a_advh_t, advz_t,  a_advz_t,             & 
      advect_tend, a_advect_tend,    &
      h_tendency, a_h_tendency, z_tendency, a_z_tendency,  & 
      msftx, msfty,                  &
      mu_old, a_mu_old, mu_new, a_mu_new, mu_base,  &
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
      INTEGER :: i,j,k,im
      REAL    :: sc_middle, msfsq
      REAL, DIMENSION(its:ite, kts:kte, jts:jte  ) :: tendency

      INTEGER :: i_start,i_end,j_start,j_end,k_start,k_end
      INTEGER :: i_start_spc,i_end_spc,j_start_spc,j_end_spc,k_start_spc,k_end_spc

      REAL, DIMENSION(its:ite) :: muold, r_munew
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce),                &
      INTENT(INOUT)                                  :: a_scalar_1, a_scalar_2

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme , scs:sce),                &
      INTENT(INOUT)                                  :: a_sc_tend

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                &
      INTENT(INOUT)                               :: a_advect_tend

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme ), OPTIONAL :: a_advh_t,  a_advz_t ! accumulating for output
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme ) :: a_h_tendency, a_z_tendency ! from rk_scalar_tend


      REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) ::  a_mu_old,  &
                                 a_mu_new

      REAL, DIMENSION(its:ite) :: a_muold, a_r_munew

      REAL, DIMENSION(its:ite, kts:kte, jts:jte  ) :: a_tendency

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, scs:sce) :: scalar_old



      !<DESCRIPTION> 
      !
      !  Basic states: mu_old, mu_new, advect_tend, sc_tend, scalar_2(rk_step=1), scalar_1(rk_step/=1)
      !
      !</DESCRIPTION>

      !  Initilize local adjoint variables
      a_muold = 0.0
      a_r_munew = 0.0
      a_tendency = 0.0
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

      !
      !
         
      DO  im = sce,scs,-1

         !     Recalculate tendency, for trajectory?
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

         !     Recalculate muold and r_munew, for traj
            DO  i = its, min(ite,ide-1)
            muold(i) = mu_old(i,j) + mu_base(i,j)
            r_munew(i) = 1./(mu_new(i,j) + mu_base(i,j))
            ENDDO

            DO  k = kts, min(kte,kde-1)
               DO  i = its, min(ite,ide-1)
               !       Recalculate scalar_1 (i.e. scalar_old)
               scalar_old(i,k,j,im) = scalar_2(i,k,j,im)
               ! adjoint of 
               ! scalar_1(i,k,j,im) = scalar_2(i,k,j,im)
               ! scalar_2(i,k,j,im) = ((c1(k)*muold(i)+c2(k))*scalar_1(i,k,j,im)   &
               ! + dt*tendency(i,k,j))/(c1(k)*munew(i)+c2(k))

               a_scalar_1(i,k,j,im) = a_scalar_1(i,k,j,im) + muold(i)*r_munew(i) * a_scalar_2(i,k,j,im)
               a_muold(i) = a_muold(i) + scalar_old(i,k,j,im)*r_munew(i) * a_scalar_2(i,k,j,im)
               a_tendency(i,k,j) = a_tendency(i,k,j) + dt*r_munew(i) * a_scalar_2(i,k,j,im)
               a_r_munew(i) = a_r_munew(i) + (muold(i)*scalar_old(i,k,j,im)+dt*tendency(i,k,j)) * a_scalar_2(i,k,j,im)
               a_scalar_2(i,k,j,im) = 0.0

               a_scalar_2(i,k,j,im) = a_scalar_2(i,k,j,im) + a_scalar_1(i,k,j,im)
               a_scalar_1(i,k,j,im) = 0.0
               ENDDO !i
            ENDDO !k

            DO  i = its, min(ite,ide-1)
               a_mu_new(i,j) = a_mu_new(i,j) - a_r_munew(i) / ((mu_new(i,j)+mu_base(i,j))*(mu_new(i,j)+mu_base(i,j)))
               a_r_munew(i) = 0.0

               a_mu_old(i,j) = a_mu_old(i,j) + a_muold(i)
               a_muold(i) = 0.0
            ENDDO

         ENDDO !j

         DO  j = j_start_spc,j_end_spc
         DO  k = k_start_spc,k_end_spc
         DO  i = i_start_spc,i_end_spc
         a_sc_tend(i,k,j,im) = a_sc_tend(i,k,j,im) + a_tendency(i,k,j)
         ENDDO
         ENDDO
         ENDDO

         DO  j = j_start,j_end
         DO  k = k_start,k_end
         DO  i = i_start,i_end
         ! scalar was coupled with my
         a_advect_tend(i,k,j) = a_advect_tend(i,k,j) + msfty(i,j) * a_tendency(i,k,j)
         a_tendency(i,k,j) = 0.0
         ENDDO
         ENDDO
         ENDDO

         DO  j = jts, min(jte,jde-1)
         DO  k = kts, min(kte,kde-1)
         DO  i = its, min(ite,ide-1)
         a_tendency(i,k,j) = 0.
         ENDDO
         ENDDO
         ENDDO

      ENDDO !im

      ELSE ! rk != 1
      !  \

      DO  im = sce, scs, -1

         !     Recalculate tendency
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

         !     Recalculate muold and r_munew
            DO  i = its, min(ite,ide-1)
            muold(i) = mu_old(i,j) + mu_base(i,j)
            r_munew(i) = 1./(mu_new(i,j) + mu_base(i,j))
            ENDDO

            ! This is separated from the k/i-loop above for better performance
            IF ( PRESENT(advh_t) .AND. PRESENT(advz_t) .AND. PRESENT(a_advh_t) .AND. PRESENT(a_advz_t) ) THEN
               IF (tenddec.and.rk_step.eq.config_flags%rk_ord) THEN
                  DO k = kts, min(kte,kde-1)
                     DO i = its, min(ite,ide-1)
                     ! adjoint of 
                     ! advh_t(i,k,j) = advh_t(i,k,j) + (dt*h_tendency(i,k,j)* msfty(i,j))/(c1(k)*munew(i)+c2(k))
                     ! advz_t(i,k,j) = advz_t(i,k,j) + (dt*z_tendency(i,k,j)* msfty(i,j))/(c1(k)*munew(i)+c2(k))
                     a_h_tendency(i,k,j) = a_h_tendency(i,k,j) + dt*msfty(i,j)*r_munew(i)*a_advh_t(i,k,j)
                     a_r_munew(i) = a_r_munew(i) + (dt*h_tendency(i,k,j)* msfty(i,j))*a_advh_t(i,k,j)
                     a_z_tendency(i,k,j) = a_z_tendency(i,k,j) + dt*msfty(i,j)*r_munew(i)*a_advz_t(i,k,j)
                     a_r_munew(i) = a_r_munew(i) + (dt*z_tendency(i,k,j)* msfty(i,j))*a_advz_t(i,k,j)

                     ENDDO
                  ENDDO
               END IF
            END IF

            DO  k = kts, min(kte,kde-1)
            DO  i = its, min(ite,ide-1)
               ! adjoint of 
               ! scalar_2(i,k,j,im) = ((c1(k)*muold(i)+c2(k))*scalar_1(i,k,j,im)   &
               !  + dt*tendency(i,k,j))/(c1(k)*munew(i)+c2(k))

               a_scalar_1(i,k,j,im) = a_scalar_1(i,k,j,im) + muold(i)*r_munew(i) * a_scalar_2(i,k,j,im)
               a_muold(i) = a_muold(i) + scalar_1(i,k,j,im)*r_munew(i) * a_scalar_2(i,k,j,im)
               a_tendency(i,k,j) = a_tendency(i,k,j) + dt*r_munew(i) * a_scalar_2(i,k,j,im)
               a_r_munew(i) = a_r_munew(i) + (muold(i)*scalar_1(i,k,j,im)+dt*tendency(i,k,j)) * a_scalar_2(i,k,j,im)

               a_scalar_2(i,k,j,im) = 0.0

            ENDDO !i
            ENDDO !k

            DO  i = its, min(ite,ide-1)
               ! adjoint of 
               ! MUOLD(i) = MU_OLD(i,j) + MU_BASE(i,j)
               ! MUNEW(i) = MU_NEW(i,j) + MU_BASE(i,j)
               a_mu_new(i,j) = a_mu_new(i,j) - a_r_munew(i) / ((mu_new(i,j)+mu_base(i,j))*(mu_new(i,j)+mu_base(i,j)))
               a_r_munew(i) = 0.0

               a_mu_old(i,j) = a_mu_old(i,j) + a_muold(i)
               a_muold(i) = 0.0
            ENDDO
         ENDDO !j

         DO  j = j_start_spc,j_end_spc
         DO  k = k_start_spc,k_end_spc
         DO  i = i_start_spc,i_end_spc
         a_sc_tend(i,k,j,im) = a_sc_tend(i,k,j,im) + a_tendency(i,k,j)
         ENDDO
         ENDDO
         ENDDO

         DO  j = j_start,j_end
         DO  k = k_start,k_end
         DO  i = i_start,i_end
         ! scalar was coupled with my
         a_advect_tend(i,k,j) = a_advect_tend(i,k,j) + msfty(i,j) * a_tendency(i,k,j)
         a_tendency(i,k,j) = 0.0
         ENDDO
         ENDDO
         ENDDO

         DO  j = jts, min(jte,jde-1)
         DO  k = kts, min(kte,kde-1)
         DO  i = its, min(ite,ide-1)
         a_tendency(i,k,j) = 0.
         ENDDO
         ENDDO
         ENDDO

      ENDDO !im

      END IF

END SUBROUTINE a_rk_update_scalar
















SUBROUTINE a_rk_tendency( config_flags, rk_step,        &
      ru_tend,a_ru_tend,rv_tend, a_rv_tend,rw_tend,a_rw_tend,ph_tend,a_ph_tend,t_tend,a_t_tend,       &
      ru_tendf, a_ru_tendf,rv_tendf,a_rv_tendf,rw_tendf,a_rw_tendf,ph_tendf,a_ph_tendf,t_tendf, a_t_tendf,              &
      mu_tend,a_mu_tend,u_save,a_u_save,v_save,a_v_save,w_save, a_w_save,ph_save,a_ph_save,        &
      t_save,a_t_save,mu_save,a_mu_save,RTHFTEN, a_RTHFTEN,                        &
      ru,a_ru,rv,a_rv,rw,a_rw,ww,a_ww,                        &
      u,a_u,v,a_v,w,a_w,t,a_t,ph,a_ph,                                  &
      u_old,a_u_old,v_old,a_v_old,w_old,a_w_old,t_old,a_t_old,ph_old, a_ph_old,              &
      h_diabatic,a_h_diabatic,phb,t_init,                          &
      mu,a_mu,mut,a_mut,muu,a_muu,muv,a_muv,mub,                  &
      al,a_al,alt,a_alt,p,a_p,pb,php,a_php,cqu,a_cqu,cqv a_cqv,cqw,a_cqw,          &
      u_base, v_base, t_base, qv_base, z_base,         &
      msfux, msfuy, msfvx, msfvx_inv,                  &
      msfvy, msftx, msfty,                             &
      clat, f, e, sina, cosa,                          &
      fnm, fnp, rdn, rdnw,                             &
      dt,rdx,rdy,khdif,kvdif,xkmhd,a_xkmhd,xkhh,a_xkhh, &
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
 
 ! PART I: DECLARATION OF VARIABLES
      IMPLICIT NONE

      ! 
      INTEGER :: K0_ADJ,K1_ADJ,K2_ADJ,K3_ADJ
      TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags

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
      REAL,                    INTENT(IN   ) :: dt
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme) :: a_ru,a_rv,a_rw,a_ww, &
                                                a_u,a_v,a_w,a_t,a_ph,a_u_old,a_v_old, &
                                                a_w_old,a_t_old,a_ph_old,a_al,a_alt,a_p, &
                                                a_php,a_cqu,a_cqv,a_xkmhd,a_xkhh, a_h_diabatic
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
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme) :: a_ru_tend,a_rv_tend, &
                                                a_rw_tend,a_t_tend,a_ph_tend,a_RTHFTEN, &
                                                a_u_save,a_v_save,a_w_save,a_ph_save,a_t_save

      REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) ,               &
                     INTENT(INOUT) :: ru_tendf, &
                                    rv_tendf, &
                                    rw_tendf, &
                                    t_tendf,  &
                                    ph_tendf, &
                                    cqw
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme) :: a_ru_tendf,a_rv_tendf,a_rw_tendf,a_t_tendf,a_ph_tendf,a_cqw
      REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(  OUT) :: mu_tend, &
                                                mu_save
      REAL,DIMENSION(ims:ime,jms:jme) :: a_mu_tend,a_mu_save
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
      REAL,DIMENSION(ims:ime,jms:jme) :: a_mu,a_mut,a_muu,a_muv
      REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fnm,     &
                                             fnp,     &
                                             rdn,     &
                                             rdnw,    &
                                             u_base,  &
                                             v_base,  &
                                             t_base,  &
                                             qv_base, &
                                             z_base,  
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



      CALL nl_get_time_step ( 1, time_step )

      ! Remarked by Ning Pan, 2010-07-30 : Part II is not needed
      !! PART! II: CALCULATIONS OF B. S. TRAJECTORY


      !! temp NING
      forward_step: IF( rk_step == 1 ) THEN

         IF( rad_nudge .eq. 1 )                                     &
            CALL a_theta_relaxation( t_tendf, a_t_tendf, t, a_t, t_init,  &
                                    mut, a_mut, ph, a_ph, phb,       &
                                    t_base, z_base,                  &
                                    ids, ide, jds, jde, kds, kde,    &
                                    ims, ime, jms, jme, kms, kme,    &
                                    its, ite, jts, jte, kts, kte   )

         IF( damp_opt .eq. 2 ) THEN
            CALL a_rk_rayleigh_damp(ru_tendf,a_ru_tendf,rv_tendf,a_rv_tendf,rw_tendf,  &
                                 a_rw_tendf,t_tendf,a_t_tendf,u,a_u,v,a_v,w,a_w,t,a_t,t_init,  &
                                 mut,a_mut,muu,a_muu,muv,a_muv,ph,a_ph,phb,u_base,v_base,t_base,z_base,  &
                                 dampcoef,zdamp,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,  &
                                 kme,its,ite,jts,jte,kts,kte)

         END IF

         IF( diff_6th_opt .NE. 0 ) THEN

            CALL a_sixth_order_diffusion('m',t,a_t,t_tendf,a_t_tendf,mut,a_mut,dt,  &
                                       config_flags,diff_6th_opt,diff_6th_factor,ids,ide,jds,jde,kds,  &
                                       kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

            IF(non_hydrostatic) THEN

               CALL a_sixth_order_diffusion('w',w,a_w,rw_tendf,a_rw_tendf,mut,a_mut,dt,  &
                                          config_flags,diff_6th_opt,diff_6th_factor,ids,ide,jds,jde,kds,  &
                                          kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

            END IF


            CALL a_sixth_order_diffusion('v',v,a_v,rv_tendf,a_rv_tendf,mut,a_mut,dt,  &
                                       config_flags,diff_6th_opt,diff_6th_factor,ids,ide,jds,jde,kds,  &
                                       kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)


            CALL a_sixth_order_diffusion('u',u,a_u,ru_tendf,a_ru_tendf,mut,a_mut,dt,  &
                                       config_flags,diff_6th_opt,diff_6th_factor,ids,ide,jds,jde,kds,  &
                                       kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

         ENDIF ! diff_6th_opt


         diff_opt1: IF(config_flags%diff_opt .eq. 1) THEN
         ! Revised by Ning Pan, 2010-07-30 : reverse the adjoint computation order
         !                                   revise actual arguments
         !                                   remark useless recalculation
            pbl_test: IF(config_flags%bl_pbl_physics .eq. 0) THEN

               kvdq = 3.*kvdif

               CALL a_vertical_diffusion_3dmp(t,a_t,t_tendf,a_t_tendf,config_flags,t_init,  &
                                             alt,a_alt,mut,a_mut,rdn,rdnw,kvdq,ids,ide,jds,jde,kds,kde,  &
                                             ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

               IF(non_hydrostatic) THEN
                  CALL a_vertical_diffusion('w',w,a_w,rw_tendf,a_rw_tendf,config_flags,alt,  &
                                             a_alt,mut,a_mut,rdn,rdnw,kvdif,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,  &
                                             its,ite,jts,jte,kts,kte)
               END IF

               CALL a_vertical_diffusion_v(v,a_v,rv_tendf,a_rv_tendf,config_flags,v_base,  &
                                             alt,a_alt,muv,a_muv,rdn,rdnw,kvdif,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,  &
                                             kme,its,ite,jts,jte,kts,kte)

               CALL a_vertical_diffusion_u(u,a_u,ru_tendf,a_ru_tendf,config_flags,u_base,  &
                                             alt,a_alt,muu,a_muu,rdn,rdnw,kvdif,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,  &
                                             kme,its,ite,jts,jte,kts,kte)

            ENDIF pbl_test

            khdq = 3.*khdif
            CALL a_horizontal_diffusion_3dmp('m',t,a_t,t_tendf,a_t_tendf,mut,a_mut,  &
                                             config_flags,t_init,msfux,msfuy,msfvx,msfvx_inv,msfvy,msftx,msfty,khdq,  &
                                             xkhh,a_xkhh,rdx,rdy,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,  &
                                             ite,jts,jte,kts,kte)

            CALL a_horizontal_diffusion('w',w,a_w,rw_tendf,a_rw_tendf,mut,a_mut,  &
                                          config_flags,msfux,msfuy,msfvx,msfvx_inv,msfvy,msftx,msfty,khdif,xkmhd,a_xkmhd,rdx,  &
                                          rdy,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

            CALL a_horizontal_diffusion('v',v,a_v,rv_tendf,a_rv_tendf,mut,a_mut,  &
                                          config_flags,msfux,msfuy,msfvx,msfvx_inv,msfvy,msftx,msfty,khdif,xkmhd,a_xkmhd,rdx,  &
                                          rdy,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

            CALL a_horizontal_diffusion('u',u,a_u,ru_tendf,a_ru_tendf,mut,a_mut,  &
                                          config_flags,msfux,msfuy,msfvx,msfvx_inv,msfvy,msftx,msfty,khdif,xkmhd,a_xkmhd,rdx,  &
                                          rdy,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

         END IF diff_opt1
      END IF forward_step

      ! IF(config_flags%ra_lw_physics == HELDSUAREZ) THEN
      
      ! END IF



      CALL a_curvature(ru,a_ru,rv,a_rv,rw,a_rw,u,a_u,v,a_v,w,ru_tend,  &
                        a_ru_tend,rv_tend,a_rv_tend,rw_tend,a_rw_tend,config_flags,msfux,msfuy,msfvx,  &
                        msfvy,msftx,msfty,xlat,fnm,fnp,rdx,rdy,ids,ide,jds,jde,kds,kde,ims,ime,jms,  &
                        jme,kms,kme,its,ite,jts,jte,kts,kte)


      IF(config_flags%pert_coriolis) THEN
         CALL a_perturbation_coriolis(ru,a_ru,rv,a_rv,rw,a_rw,ru_tend,a_ru_tend,  &
                                    rv_tend,a_rv_tend,rw_tend,a_rw_tend,config_flags,u_base,v_base,z_base,muu,  &
                                    a_muu,muv,a_muv,phb,ph,a_ph,msftx,msfty,msfux,msfuy,msfvx,msfvy,f,e,sina,cosa,  &
                                    fnm,fnp,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      ELSE
         CALL a_coriolis(ru,a_ru,rv,a_rv,rw,a_rw,ru_tend,a_ru_tend,rv_tend,  &
                  a_rv_tend,rw_tend,a_rw_tend,config_flags,msftx,msfty,msfux,msfuy,msfvx,msfvy,f,e,  &
                  sina,cosa,fnm,fnp,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      END IF ! config_flags%pert_coriolis


      CALL a_w_damp(rw_tend,a_rw_tend,max_vert_cfl,max_horiz_cfl,  &
                     u,a_u,v,a_v,ww,a_ww,w,a_w,mut,a_mut,rdnw,rdx,rdy,msfux,  &
                     msfuy,msfvx,msfvy,dt,config_flags,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,  &
                     its,ite,jts,jte,kts,kte)


      IF(non_hydrostatic) THEN
         CALL a_pg_buoy_w(rw_tend,a_rw_tend,p,a_p,cqw,a_cqw,mu,a_mu,mub,rdnw,rdn,  &
                           g,msftx,msfty,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      ENDIF

      CALL a_horizontal_pressure_gradient(ru_tend,a_ru_tend,rv_tend,a_rv_tend,ph,  &
                                          a_ph,alt,a_alt,p,a_p,pb,al,a_al,php,a_php,cqu,a_cqu,cqv,a_cqv,muu,  &
                                          a_muu,muv,a_muv,mu,a_mu,fnm,fnp,rdnw,cf1,cf2,cf3,cfn,cfn1,rdx,rdy,msfux,msfuy,msfvx,  &
                                          msfvy,msftx,msfty,config_flags,non_hydrostatic,top_lid,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,  &
                                          ite,jts,jte,kts,kte)
      ! IF (non_hydrostatic) THEN  CALL advect_w_implicit ( w, w_old, rw_tend,    
      ! IF( ieva ) THEN CALL advect_ph_implicit ( ph, ph_old, ph_tend, phb, 

      CALL a_rhs_ph(ph_tend,a_ph_tend,u,a_u,v,a_v,ww,a_ww,ph,a_ph,ph,a_ph,  &
                     phb,w,a_w,mut,a_mut,muu,a_muu,muv,a_muv,fnm,fnp,rdnw,cfn,cfn1,rdx,rdy,msfux,  &
                     msfuy,msfvx,msfvx_inv,msfvy,msftx,msfty,non_hydrostatic,config_flags,ids,ide,jds,jde,kds,kde,ims,  &
                     ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)


      IF( config_flags%cu_physics == GDSCHEME  .OR.       &
         config_flags%cu_physics == G3SCHEME  .OR.       &
         config_flags%cu_physics == NTIEDTKESCHEME ) THEN
         
            CALL a_set_tend(RTHFTEN,a_RTHFTEN,t_tend,a_t_tend,msfty,ids,ide,jds,jde,kds,  &
                           kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
      END IF

      ! IF( ieva ) THEN CALL advect_s_implicit ( t, t_old, t_tend, ru, rv, wwI,  &


      !hcl 11/2016 ERM: Use WENO for theta flux on 3rd RK step if using WENO_SCALAR or WENOPD_SCALAR
      ! to be consistent with other scalar fluxes
      IF(  ( config_flags%scalar_adv_opt == WENO_SCALAR           &
            .or. config_flags%scalar_adv_opt == WENOPD_SCALAR   &
            .or. config_flags%moist_adv_opt == WENO_SCALAR      &
            .or. config_flags%moist_adv_opt == WENOPD_SCALAR    &
                        )  .and. (rk_step == 3) ) THEN

         CALL a_advect_scalar_weno(t,a_t,t,a_t,t_tend,a_t_tend,ru,a_ru,rv,a_rv,ww,  &
                              a_ww,mut,time_step,config_flags,msfux,msfuy,msfvx,msfvy,msftx,msfty,fnm,  &
                              fnp,rdx,rdy,rdnw,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

      ELSE

         CALL a_advect_scalar(t,a_t,t,a_t,t_tend,a_t_tend,ru,a_ru,rv,a_rv,ww,  &
                              a_ww,mut,time_step,config_flags,msfux,msfuy,msfvx,msfvy,msftx,msfty,fnm,  &
                              fnp,rdx,rdy,rdnw,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

      ENDIF



      IF(non_hydrostatic) THEN
         IF( (rk_step == 3) .and. ( adv_opt == WENO_MOM ) ) THEN
            CALL a_advect_weno_w ( w, a_w, w, a_w, rw_tend, a_rw_tend, &
                              ru, a_ru, rv, a_rv, ww, a_ww, &
                              mut, time_step, config_flags, &
                              msfux, msfuy, msfvx, msfvy,   &
                              msftx, msfty,                 &
                              fnm, fnp, rdx, rdy, rdn,      &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte )
         ELSE
            CALL a_advect_w(w,a_w,w,a_w,rw_tend,a_rw_tend,ru,a_ru,rv,a_rv,ww,  &
                           a_ww,mut,time_step,config_flags,msfux,msfuy,msfvx,msfvy,msftx,msfty,fnm,  &
                           fnp,rdx,rdy,rdn,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

         END IF
      END IF
      ! IF( ieva ) THEN CALL advect_v_implicit ( v, v_old, rv_tend, ru, rv, wwI,  &


      IF( (rk_step == 3) .and. ( adv_opt == WENO_MOM ) ) THEN

         CALL a_advect_weno_v ( v, a_v, v, a_v, rv_tend, a_rv_tend, &
                                 ru, a_ru, rv, a_rv, ww, a_ww, &
                                 mut, a_mut, time_step, config_flags, &
                                 msfux, msfuy, msfvx, msfvy,   &
                                 msftx, msfty,                 &
                                 fnm, fnp, rdx, rdy, rdnw,     &
                                 ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte )
      ELSE
         CALL a_advect_v(v,a_v,v,a_v,rv_tend,a_rv_tend,ru,a_ru,rv,a_rv,ww,  &
                           a_ww,mut,a_mut,time_step,config_flags,msfux,msfuy,msfvx,msfvy,msftx,msfty,fnm,  &
                           fnp,rdx,rdy,rdnw,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

      ENDIF


      IF( (rk_step == 3) .and. ( adv_opt == WENO_MOM ) ) THEN

         CALL a_advect_weno_u ( u, a_u, u, a_u, ru_tend, a_ru_tend, &
                                 ru, a_ru, rv, a_rv, ww, a_ww, &
                                 mut, a_mut, time_step, config_flags, &
                                 msfux, msfuy, msfvx, msfvy,   &
                                 msftx, msfty,                 &
                                 fnm, fnp, rdx, rdy, rdnw,     &
                                 ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte )
      ELSE
         CALL a_advect_u(u,a_u,u,a_u,ru_tend,a_ru_tend,ru,a_ru,rv,a_rv,ww,  &
                        a_ww,mut,a_mut,time_step,config_flags,msfux,msfuy,msfvx,msfvy,msftx,msfty,fnm,  &
                        fnp,rdx,rdy,rdnw,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

      ENDIF

      ! IF( ieva ) THEN CALL calc_mu_uv_1 ( config_flags,
      ! calc_mu_uv_1 ( config_flags,                
      ! calc_mut_new ( u, v, c1h, c2h,          
      ! calculate_full( mut_old,  mub, mu_old,     &   

      ! WW_SPLIT(wwE, wwI,          

      ! Added by Ning Pan, 2010-07-30
      CALL a_zero_tend2d(a_mu_save,ids,ide,jds,jde,1,1,ims,ime,jms,jme,1,1,  &
      its,ite,jts,jte,1,1)
      CALL a_zero_tend2d(a_mu_tend,ids,ide,jds,jde,1,1,ims,ime,jms,jme,1,1,  &
      its,ite,jts,jte,1,1)
      CALL a_zero_tend(a_t_save,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,  &
      kme,its,ite,jts,jte,kts,kte)
      CALL a_zero_tend(a_ph_save,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,  &
      kms,kme,its,ite,jts,jte,kts,kte)
      CALL a_zero_tend(a_w_save,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,  &
      kme,its,ite,jts,jte,kts,kte)
      CALL a_zero_tend(a_v_save,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,  &
      kme,its,ite,jts,jte,kts,kte)
      CALL a_zero_tend(a_u_save,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,  &
      kme,its,ite,jts,jte,kts,kte)
      CALL a_zero_tend(a_ph_tend,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,  &
      kms,kme,its,ite,jts,jte,kts,kte)
      CALL a_zero_tend(a_t_tend,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,  &
      kme,its,ite,jts,jte,kts,kte)
      CALL a_zero_tend(a_rw_tend,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,  &
      kms,kme,its,ite,jts,jte,kts,kte)
      CALL a_zero_tend(a_rv_tend,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,  &
      kms,kme,its,ite,jts,jte,kts,kte)
      CALL a_zero_tend(a_ru_tend,ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,  &
      kms,kme,its,ite,jts,jte,kts,kte)

END SUBROUTINE a_rk_tendency