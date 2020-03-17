!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_sf_sfclay

 !USE ccpp_kinds, ONLY: kind_phys
USE sf_sfclay

 PUBLIC :: sfclay

 REAL    , PARAMETER ::  VCONVC=1.
 REAL    , PARAMETER ::  CZO=0.0185
 REAL    , PARAMETER ::  OZO=1.59E-5

 REAL,   DIMENSION(0:1000 ),SAVE :: psim_stab,psim_unstab,psih_stab,psih_unstab

CONTAINS

!-------------------------------------------------------------------
   SUBROUTINE SFCLAY(U3D,V3D,T3D,QV3D,P3D,dz8w,                    &
                     CP,G,ROVCP,R,XLV,PSFC,CHS,CHS2,CQS2,CPM,      &
                     ZNT,UST,PBLH,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
                     FM,FH,                                        &
                     XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
                     U10,V10,TH2,T2,Q2,                            &
                     GZ1OZ0,WSPD,BR,ISFFLX,DX,                     &
                     SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                     KARMAN,                                       &
                     P1000mb,                                      &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte,                             &
                     ustm,ck,cka,cd,cda,isftcflx,iz0tlnd,scm_force_flux,   &
                     errmsg,errflg           )
!-------------------------------------------------------------------

      IMPLICIT NONE
!-------------------------------------------------------------------
!   Changes in V3.7 over water surfaces:
!          1. for ZNT/Cd, replacing constant OZO with 0.11*1.5E-5/UST(I)
!             the COARE 3.5 (Edson et al. 2013) formulation is also available
!          2. for VCONV, reducing magnitude by half
!          3. for Ck, replacing Carlson-Boland with COARE 3
!-------------------------------------------------------------------
!-- U3D         3D u-velocity interpolated to theta points (m/s)
!-- V3D         3D v-velocity interpolated to theta points (m/s)
!-- T3D         temperature (K)
!-- QV3D        3D water vapor mixing ratio (Kg/Kg)
!-- P3D         3D pressure (Pa)
!-- dz8w        dz between full levels (m)
!-- CP          heat capacity at constant pressure for dry air (J/kg/K)
!-- G           acceleration due to gravity (m/s^2)
!-- ROVCP       R/CP
!-- R           gas constant for dry air (J/kg/K)
!-- XLV         latent heat of vaporization for water (J/kg)
!-- PSFC        surface pressure (Pa)
!-- ZNT         roughness length (m)
!-- UST         u* in similarity theory (m/s)
!-- USTM        u* in similarity theory (m/s) without vconv correction
!               used to couple with TKE scheme
!-- PBLH        PBL height from previous time (m)
!-- MAVAIL      surface moisture availability (between 0 and 1)
!-- ZOL         z/L height over Monin-Obukhov length
!-- MOL         T* (similarity theory) (K)
!-- REGIME      flag indicating PBL regime (stable, unstable, etc.)
!-- PSIM        similarity stability function for momentum
!-- PSIH        similarity stability function for heat
!-- FM          integrated stability function for momentum
!-- FH          integrated stability function for heat
!-- XLAND       land mask (1 for land, 2 for water)
!-- HFX         upward heat flux at the surface (W/m^2)
!-- QFX         upward moisture flux at the surface (kg/m^2/s)
!-- LH          net upward latent heat flux at surface (W/m^2)
!-- TSK         surface temperature (K)
!-- FLHC        exchange coefficient for heat (W/m^2/K)
!-- FLQC        exchange coefficient for moisture (kg/m^2/s)
!-- CHS         heat/moisture exchange coefficient for LSM (m/s)
!-- QGH         lowest-level saturated mixing ratio
!-- QSFC        ground saturated mixing ratio
!-- U10         diagnostic 10m u wind
!-- V10         diagnostic 10m v wind
!-- TH2         diagnostic 2m theta (K)
!-- T2          diagnostic 2m temperature (K)
!-- Q2          diagnostic 2m mixing ratio (kg/kg)
!-- GZ1OZ0      log(z/z0) where z0 is roughness length
!-- WSPD        wind speed at lowest model level (m/s)
!-- BR          bulk Richardson number in surface layer
!-- ISFFLX      isfflx=1 for surface heat and moisture fluxes
!-- DX          horizontal grid size (m)
!-- SVP1        constant for saturation vapor pressure (kPa)
!-- SVP2        constant for saturation vapor pressure (dimensionless)
!-- SVP3        constant for saturation vapor pressure (K)
!-- SVPT0       constant for saturation vapor pressure (K)
!-- EP1         constant for virtual temperature (R_v/R_d - 1) (dimensionless)
!-- EP2         constant for specific humidity calculation 
!               (R_d/R_v) (dimensionless)
!-- KARMAN      Von Karman constant
!-- ck          enthalpy exchange coeff at 10 meters
!-- cd          momentum exchange coeff at 10 meters
!-- cka         enthalpy exchange coeff at the lowest model level
!-- cda         momentum exchange coeff at the lowest model level
!-- isftcflx    =0, (Charnock and Carlson-Boland); =1, AHW Ck, Cd, =2 Garratt
!-- iz0tlnd     =0 Carlson-Boland, =1 Czil_new
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-------------------------------------------------------------------
      INTEGER,  INTENT(IN )   ::        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte
!                                                               
      INTEGER,  INTENT(IN )   ::        ISFFLX
      REAL,     INTENT(IN )   ::        SVP1,SVP2,SVP3,SVPT0
      REAL,     INTENT(IN )   ::        EP1,EP2,KARMAN
      REAL,     INTENT(IN )   ::        P1000mb
!
      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                           dz8w
                                        
      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                           QV3D, &
                                                              P3D, &
                                                              T3D

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN   )               ::             MAVAIL, &
                                                             PBLH, &
                                                            XLAND, &
                                                              TSK
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(OUT  )               ::                U10, &
                                                              V10, &
                                                              TH2, &
                                                               T2, &
                                                               Q2, &
                                                             QSFC

!
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(OUT  )               ::             REGIME, &
                                                              HFX, &
                                                              QFX, &
                                                               LH, &
                                                          MOL,RMOL
! the following 5 are change to memory size
!
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(OUT  )   ::                 GZ1OZ0,WSPD,BR, &
                                                  PSIM,PSIH,FM,FH

      REAL,     DIMENSION( ims:ime, kms:kme, jms:jme )           , &
                INTENT(IN   )   ::                            U3D, &
                                                              V3D
                                        
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN   )               ::               PSFC

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(OUT  )   ::                            ZNT, &
                                                              ZOL, &
                                                              UST, &
                                                              CPM, &
                                                             CHS2, &
                                                             CQS2, &
                                                              CHS

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(OUT  )   ::                      FLHC,FLQC

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(OUT  )   ::                            QGH
                                    
      REAL,     INTENT(IN   )               ::   CP,G,ROVCP,R,XLV

      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN   )               ::                 DX
 
      REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme )              , &
                INTENT(OUT)     ::                  ck,cka,cd,cda

      REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme )              , &
                INTENT(OUT  )   ::                           USTM

      INTEGER,  OPTIONAL,  INTENT(IN )   ::     ISFTCFLX, IZ0TLND
      INTEGER,  OPTIONAL,  INTENT(IN )   ::     SCM_FORCE_FLUX

      CHARACTER(LEN=*), INTENT(OUT  ) ::                   errmsg
      INTEGER,          INTENT(OUT  ) ::                   errflg

! LOCAL VARS

      REAL,     DIMENSION( its:ite ) ::                       U1D, &
                                                              V1D, &
                                                             QV1D, &
                                                              P1D, &
                                                              T1D

      REAL,     DIMENSION( its:ite ) ::                    dz8w1d

      REAL,     DIMENSION( its:ite ) ::                      DX2D

      INTEGER ::  I,J

!  Local tile-sized arrays for contiguous data for sf_sfclay_run call.

   real,                   dimension(its:ite)      ::     psfc_hv, &
                                                           chs_hv, &
                                                          chs2_hv, &
                                                          cqs2_hv, &
                                                           cpm_hv, &
                                                          pblh_hv, &
                                                          rmol_hv, &
                                                           znt_hv, &
                                                           ust_hv, &
                                                        mavail_hv, &
                                                           zol_hv, &
                                                           mol_hv, &
                                                        regime_hv, &
                                                          psim_hv, &
                                                          psih_hv, &
                                                            fm_hv, &
                                                            fh_hv, &
                                                         xland_hv, &
                                                           hfx_hv, &
                                                           qfx_hv, &
                                                           tsk_hv, &
                                                           u10_hv, &
                                                           v10_hv, &
                                                           th2_hv, &
                                                            t2_hv, &
                                                            q2_hv, &
                                                          flhc_hv, &
                                                          flqc_hv, &
                                                           qgh_hv, &
                                                          qsfc_hv, &
                                                            lh_hv, &
                                                        gz1oz0_hv, &
                                                          wspd_hv, &
                                                            br_hv, &
                                                          ustm_hv, &
                                                            ck_hv, &
                                                           cka_hv, &
                                                            cd_hv, &
                                                           cda_hv
 
 !  Assign input data to local tile-sized arrays.

      DO J=jts,jte

        DO i=its,ite
           DX2D(i)=DX(i,j)
        ENDDO

        DO i=its,ite
          dz8w1d(I) = dz8w(i,1,j)
        ENDDO
   
        DO i=its,ite
           U1D(i) =U3D(i,1,j)
           V1D(i) =V3D(i,1,j)
           QV1D(i)=QV3D(i,1,j)
           P1D(i) =P3D(i,1,j)
           T1D(i) =T3D(i,1,j)
        ENDDO
 
        do i = its, ite
           psfc_hv(i) = psfc(i,j)
           chs_hv(i) = chs(i,j)
           chs2_hv(i) = chs2(i,j)
           cqs2_hv(i) = cqs2(i,j)
           cpm_hv(i) = cpm(i,j)
           pblh_hv(i) = pblh(i,j)
           rmol_hv(i) = rmol(i,j)
           znt_hv(i) = znt(i,j)
           ust_hv(i) = ust(i,j)
           mavail_hv(i) = mavail(i,j)
           zol_hv(i) = zol(i,j)
           mol_hv(i) = mol(i,j)
           regime_hv(i) = regime(i,j)
           psim_hv(i) = psim(i,j)
           psih_hv(i) = psih(i,j)
           fm_hv(i) = fm(i,j)
           fh_hv(i) = fh(i,j)
           xland_hv(i) = xland(i,j)
           hfx_hv(i) = hfx(i,j)
           qfx_hv(i) = qfx(i,j)
           tsk_hv(i) = tsk(i,j)
           u10_hv(i) = u10(i,j)
           v10_hv(i) = v10(i,j)
           th2_hv(i) = th2(i,j)
           t2_hv(i) = t2(i,j)
           q2_hv(i) = q2(i,j)
           flhc_hv(i) = flhc(i,j)
           flqc_hv(i) = flqc(i,j)
           qgh_hv(i) = qgh(i,j)
           qsfc_hv(i) = qsfc(i,j)
           lh_hv(i) = lh(i,j)
           gz1oz0_hv(i) = gz1oz0(i,j)
           wspd_hv(i) = wspd(i,j)
           br_hv(i) = br(i,j)
           ustm_hv(i) = ustm(i,j)
           ck_hv(i) = ck(i,j)
           cka_hv(i) = cka(i,j)
           cd_hv(i) = cd(i,j)
           cda_hv(i) = cda(i,j)
        end do

        !  Sending array starting locations of optional variables may cause
        !  troubles, so we explicitly change the call.

        CALL sf_sfclay_run(U1D,V1D,T1D=t1d,QV1D=qv1d,P1D=p1d,dz8w1d=dz8w1d,                  &
                CP=cp,G=g,ROVCP=rovcp,R=r,XLV=xlv,psfcpa=psfc_hv,CHS=chs_hv,CHS2=chs2_hv, &
                CQS2=cqs2_hv,CPM=cpm_hv,PBLH=pblh_hv,RMOL=rmol_hv,                      &
                ZNT=znt_hv,UST=ust_hv,MAVAIL=mavail_hv,ZOL=zol_hv,                      &
                MOL=mol_hv,REGIME=regime_hv,PSIM=psim_hv,PSIH=psih_hv,                  &
                FM=fm_hv,FH=fh_hv,                                                      &
                XLAND=xland_hv,HFX=hfx_hv,QFX=qfx_hv,TSK=tsk_hv,                        &
                U10=u10_hv,V10=v10_hv,TH2=th2_hv,T2=t2_hv,                              &
                Q2=q2_hv,FLHC=flhc_hv,FLQC=flqc_hv,QGH=qgh_hv,                          &
                QSFC=qsfc_hv,LH=lh_hv,                                                  &
                GZ1OZ0=gz1oz0_hv,WSPD=wspd_hv,BR=br_hv,ISFFLX=isfflx,DX=dx2d,           &
                SVP1=svp1,SVP2=svp2,SVP3=svp3,SVPT0=svpt0,EP1=ep1,EP2=ep2,KARMAN=karman,&
                P1000mb=p1000mb,                                                        &
                its=its,ite=ite,                                                        &
                errmsg=errmsg,errflg=errflg                                             &
#if ( EM_CORE == 1 )
                ,isftcflx=isftcflx,iz0tlnd=iz0tlnd,scm_force_flux=scm_force_flux,       &
                USTM=ustm_hv,CK=ck_hv,CKA=cka_hv,                                       &
                CD=cd_hv,CDA=cda_hv                                                     &
#endif
                     ) 

 !  Assign local data back to full-sized arrays.
      !  Only required for the INTENT(OUT) or INTENT(INOUT) arrays.

        do i = its, ite
           chs(i,j) = chs_hv(i)
           chs2(i,j) = chs2_hv(i)
           cqs2(i,j) = cqs2_hv(i)
           cpm(i,j) = cpm_hv(i)
           rmol(i,j) = rmol_hv(i)
           znt(i,j) = znt_hv(i)
           ust(i,j) = ust_hv(i)
           zol(i,j) = zol_hv(i)
           mol(i,j) = mol_hv(i)
           regime(i,j) = regime_hv(i)
           psim(i,j) = psim_hv(i)
           psih(i,j) = psih_hv(i)
           fm(i,j) = fm_hv(i)
           fh(i,j) = fh_hv(i)
           hfx(i,j) = hfx_hv(i)
           qfx(i,j) = qfx_hv(i)
           u10(i,j) = u10_hv(i)
           v10(i,j) = v10_hv(i)
           th2(i,j) = th2_hv(i)
           t2(i,j) = t2_hv(i)
           q2(i,j) = q2_hv(i)
           flhc(i,j) = flhc_hv(i)
           flqc(i,j) = flqc_hv(i)
           qgh(i,j) = qgh_hv(i)
           qsfc(i,j) = qsfc_hv(i)
           lh(i,j) = lh_hv(i)
           gz1oz0(i,j) = gz1oz0_hv(i)
           wspd(i,j) = wspd_hv(i)
           br(i,j) = br_hv(i)
           ustm(i,j) = ustm_hv(i)
           ck(i,j) = ck_hv(i)
           cka(i,j) = cka_hv(i)
           cd(i,j) = cd_hv(i)
           cda(i,j) = cda_hv(i)
        end do

      ENDDO


   END SUBROUTINE SFCLAY

END MODULE module_sf_sfclay
