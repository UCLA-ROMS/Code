      subroutine lmd_kpp_tile (istr,iend,jstr,jend,       Kv,Kt,Ks,
     &                         ustar, Bo,Bosol, hbl,hbbl, FX,FE,FE1,
     &                         Cr,FC, wrk1,wrk2,
     &                         Gm1,dGm1dS,  Gt1,dGt1dS,  Gs1,dGs1dS,
     &                                                     kbl,kmo)
      implicit none
      integer(kind=4), parameter ::
     &               LLm=435, MMm=660, N=60
      integer(kind=4), parameter ::
     &      NP_XI=8, NP_ETA=32, NSUB_X=1, NSUB_E=1
      integer(kind=4), parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA
      integer(kind=4) ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
      logical west_exchng,  east_exchng
      logical south_exchng, north_exchng
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
     &                , west_exchng,  east_exchng
     &                , south_exchng, north_exchng
      integer(kind=4), parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
     &       , itemp=1
     &       , isalt=2
     &       , ntrc_salt=1
      integer(kind=4), parameter :: ntrc_pas=0
      integer(kind=4), parameter :: itrc_bio=itemp+ntrc_salt+ntrc_pas+1
      integer(kind=4), parameter :: iPO4=itrc_bio,iNO3=iPO4+1, 
     &                           iSIO3=iPO4+2,
     &     iNH4=iPO4+3,
     &     iFE=iPO4+4, iO2=iPO4+5, iDIC=iPO4+6,
     &     iALK=iPO4+7, iDOC=iPO4+8, iSPC=iPO4+9,
     &     iSPCHL=iPO4+10, iSPCACO3=iPO4+11, iDIATC=iPO4+12,
     &     iDIATCHL=iPO4+13, iZOOC=iPO4+14, iSPFE=iPO4+15,
     &     iDIATSI=iPO4+16, iDIATFE=iPO4+17, iDIAZC=iPO4+18,
     &     iDIAZCHL=iPO4+19, iDIAZFE=iPO4+20, iDON=iPO4+21,
     &     iDOFE=iPO4+22, iDOP=iPO4+23
      integer(kind=4), parameter :: iNO2 = iDOP + 1
      integer(kind=4), parameter :: iN2O = iNO2 + 1, iN2 = iN2O + 1
      integer(kind=4), parameter :: iSedOrgC = 1
      integer(kind=4), parameter :: iSedCaCO3 = iSedOrgC + 1
      integer(kind=4), parameter :: iSedSi = iSedCaCO3 + 1
      integer(kind=4), parameter :: NT_sed = 3
      integer(kind=4), parameter :: ntrc_bio = 24
     &     + 3
      integer(kind=4), parameter :: NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio
      integer(kind=4) istr,iend,jstr,jend, i,j,k
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2,0:N) :: Kv, 
     &                               Kt, Ks
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: ustar, 
     &                             Bo, Bosol
     &                                                  , hbl, hbbl
     &                                              , FX, FE, FE1
      real(kind=8), dimension(istr-2:iend+2,0:N) :: Cr,FC, wrk1,wrk2
      real(kind=8), dimension(istr-2:iend+2) ::  Bfsfc_bl,
     &                               Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS
      integer(kind=4), dimension(istr-2:iend+2) :: kbl, kmo
      real(kind=8), parameter ::
     &   Ricr=0.15_8,
     &   Ri_inv=1._8/Ricr,
     &   epssfc=0.1_8,
     &   betaT=-0.2_8,
     &   nubl=0.01_8,
     &   nu0c=0.1_8,
     &   Cv=1.8_8,
     &   C_MO=1._8,
     &   C_Ek=258._8,
     &   Cstar=10._8,
     &   zeta_m=-0.2_8,
     &   a_m=1.257_8,
     &   c_m=8.360_8,
     &   zeta_s=-1.0_8,
     &   a_s=-28.86_8,
     &   c_s=98.96_8
      real(kind=8), parameter :: eps=1.D-20,  r2=0.5_8, r3=1._8/3._8, 
     &                             r4=0.25_8
      real(kind=8) Cg, ustar3, Bfsfc, zscale, zetahat, ws,wm, Kern, 
     &                             Vtc,Vtsq,
     &  sigma, z_bl, Av_bl,dAv_bl, At_bl,dAt_bl,  As_bl,dAs_bl,
     &                       f1,a1,a2,a3,  cff,cff1, cff_up,cff_dn
      real(kind=8)   Kv0, Kt0, Ks0
      real(kind=8) h(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) hinv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) f(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) fomn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grd_h/h /grd_hinv/hinv /grd_f/f /grd_fomn/fomn
      real(kind=8) angler(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grd_angler/angler
      real(kind=8) latr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) lonr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) mycoeff(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grd_latr/latr /grd_lonr/lonr /grd_mycoeff/mycoeff
      real(kind=8) pm(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) pn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dm_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dn_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dm_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dn_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dm_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dn_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dm_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dn_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mtrix_pm/pm     /mtrix_pn/pn
     &       /mtrix_dm_r/dm_r /mtrix_dn_r/dn_r
     &       /mtrix_dm_u/dm_u /mtrix_dn_u/dn_u
     &       /mtrix_dm_v/dm_v /mtrix_dn_v/dn_v
     &       /mtrix_dm_p/dm_p /mtrix_dn_p/dn_p
      real(kind=8) iA_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) iA_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mtrix_iAu/iA_u  /mtrix_iAv/iA_v
      real(kind=8) dmde(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) dndx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mtrix_dmde/dmde   /mtrix_dndx/dndx
      real(kind=8) pmon_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) pnom_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) grdscl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mtrix_pmon_u/pmon_u /mtrix_pnom_v/pnom_v
     &                            /mtrix_grdscl/grdscl
      real(kind=8) rmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) pmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) umask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) vmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mask_r/rmask /mask_p/pmask
     &       /mask_u/umask /mask_v/vmask
      real(kind=8) u(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real(kind=8) v(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real(kind=8) t(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real(kind=8) FlxU(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) FlxV(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) We(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real(kind=8) Wi(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /flx_FU/FlxU /flx_FV/FlxV /flx_We/We /flx_Wi/Wi
      real(kind=8) Hz(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) z_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) z_w(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /grid_zw/z_w /grid_zr/z_r /grid_Hz/Hz
      real(kind=8) sustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) svstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      integer(kind=4) ntrad,nttra,ntprec,ntwnd
      common /blk_nt/ntrad,nttra,ntprec,ntwnd
      integer(kind=4) tra_file_id, rad_file_id, prec_file_id,
     &        wnd_file_id
      common /blk_id/ tra_file_id,rad_file_id,
     &                prec_file_id,wnd_file_id
      real(kind=8) stflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_stflx/stflx
      real(kind=8) ust(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulk_ust/ust
      real(kind=8) tst(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulk_tst/tst
      real(kind=8) qst(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulk_qst/qst
      real(kind=8) prec_scale,tra_scale,wnd_scale,srf_scale
      common /blk_scale/prec_scale,tra_scale,wnd_scale,srf_scale
      real(kind=8) tair(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) qair(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) rain(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) radlw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) radsw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) uwnd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) vwnd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /bulk_tair/tair
     &       /bulk_qair/qair
     &       /bulk_rain/rain
     &       /bulk_radlw/radlw
     &       /bulk_radsw/radsw
     &       /bulk_uwnd/uwnd
     &       /bulk_vwnd/vwnd
      real(kind=8) tairg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) qairg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) raing(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) radlwg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) radswg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) uwndg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) vwndg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulkdat_tairg/tairg
     &       /bulkdat_qairg/qairg
     &       /bulkdat_raing/raing
     &       /bulkdat_radlwg/radlwg
     &       /bulkdat_radswg/radswg
     &       /bulk_uwndg/uwndg
     &       /bulk_vwndg/vwndg
      real(kind=8)    rad_time(2), rad_cycle
      real(kind=8)    tra_time(2), tra_cycle
      real(kind=8)    wnd_time(2), wnd_cycle
      real(kind=8)    prec_time(2),prec_cycle
      real(kind=8)    ztref, zref
      common /bulkdat2/
     &        rad_time, rad_cycle,
     &        tra_time, tra_cycle,
     &        wnd_time, wnd_cycle,
     &        prec_time, prec_cycle
     &        ,ztref, zref
      integer(kind=4) tair_id,uwnd_id,radsw_id,prec_id
      integer(kind=4) qair_id,vwnd_id,radlw_id
      integer(kind=4) itrad,rad_ncycle,rad_rec,rad_tid
      integer(kind=4) ittra,tra_ncycle,tra_rec,tra_tid
      integer(kind=4) itwnd,wnd_ncycle,wnd_rec,wnd_tid
      integer(kind=4) itprec,prec_ncycle,prec_rec,prec_tid
      common /bulkdat1/
     &        tair_id,uwnd_id,radsw_id,prec_id,
     &        qair_id,vwnd_id,radlw_id,
     &        itrad,rad_ncycle,rad_rec,rad_tid,
     &        ittra,tra_ncycle,tra_rec,tra_tid,
     &        itwnd,wnd_ncycle,wnd_rec,wnd_tid,
     &        itprec,prec_ncycle,prec_rec,prec_tid
      real(kind=8) sssg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /sss_dat/sssg
      real(kind=8) sss_cycle, sss_time(2)
      integer(kind=4) sss_ncycle,   sss_rec,  itsss,    ntsss,
     &        sss_file_id,  sss_id,   sss_tid
      common /sssrest_data/ sss_cycle,sss_time,sss_ncycle,
     &  sss_rec,itsss,ntsss,sss_file_id,sss_id, sss_tid
      real(kind=8), parameter :: coeff_ssuv =
     &     0.77_8
      real(kind=8) sustr_blk(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_sustr_blk/sustr_blk
      real(kind=8) svstr_blk(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_svstr_blk/svstr_blk
      real(kind=8) shflx_net(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_net/shflx_net
      real(kind=8) shflx_lat(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_lat/shflx_lat
      real(kind=8) shflx_sen(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_sen/shflx_sen
      real(kind=8) shflx_rad(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_rad/shflx_rad
      real(kind=8) swflx_emp(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_swflx_emp/swflx_emp
      real(kind=8) shflx_wwk(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_wwk/shflx_wwk
      real(kind=8) surf_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_surf_u/surf_u
      real(kind=8) surf_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_surf_v/surf_v
      real(kind=8) time_bulk_diags_his
      common /t_bulk_diags_his/time_bulk_diags_his
      real(kind=8) zeta_bulk_diags_his(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /zeta_bulk_diags_his/zeta_bulk_diags_his
      logical new_bulk_diags_his
      integer(kind=4) n_bulk_diags_his
      common /nc_bulk_diags_his/n_bulk_diags_his,new_bulk_diags_his
      integer(kind=4) nrpf_bulk_diags_his,ncid_bulk_diags_his
     &       , nrec_bulk_diags_his,bulk_diags_hisTstep, bulk_diags_hisZ
     &       , bulk_diags_hisTime
      common /nc_bulk_diags_his/ nrpf_bulk_diags_his
     &     , ncid_bulk_diags_his, nrec_bulk_diags_his
     &     , bulk_diags_hisTstep, bulk_diags_hisZ
     &     , bulk_diags_hisTime
      character(len=80) bulk_diags_his_name
      common /nc_bulk_diags_his/bulk_diags_his_name
      real(kind=8) time_bulk_diags_avg
      common /t_bulk_diags_avg/time_bulk_diags_avg
      real(kind=8) zeta_bulk_diags_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /zeta_bulk_diags_avg/zeta_bulk_diags_avg
      logical new_bulk_diags_avg
      integer(kind=4) nts_bulk_diags_avg, n_bulk_diags_avg
      common /nc_bulk_diags_avg/nts_bulk_diags_avg
     $        ,n_bulk_diags_avg,new_bulk_diags_avg
      integer(kind=4) nrpf_bulk_diags_avg,ncid_bulk_diags_avg
     &       , nrec_bulk_diags_avg,bulk_diags_avgTstep, bulk_diags_avgZ
     &       , bulk_diags_avgTime
      common /nc_bulk_diags_avg/ nrpf_bulk_diags_avg
     &     , ncid_bulk_diags_avg, nrec_bulk_diags_avg
     &     , bulk_diags_avgTstep, bulk_diags_avgZ
     &     , bulk_diags_avgTime
      character(len=80) bulk_diags_avg_name
      common /nc_bulk_diags_avg/bulk_diags_avg_name
         real(kind=8) dust(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
           common /forces_dust/dust
       real(kind=8) dustg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
       common /dustdat_dustg/dustg
        real(kind=8) dustp(2), dust_time(2),dust_cycle, scldqdt
        integer(kind=4) itdust,dust_id,ldustgrd ,dust_ncycle,
     &  dust_rec,dust_tid,dust_file_id,iron_file_id
       common/dustdat/itdust,dust_id,ldustgrd,
     &  dust_ncycle,dust_rec,dust_tid,dust_file_id,iron_file_id
       common/dustdat1/dustp,dust_time,dust_cycle,scldqdt
       real(kind=8) iron(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       common /forces_iron/iron
       real(kind=8) irong(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
       common /irondat_irong/irong
       real(kind=8) ironp(2),iron_time(2),iron_cycle
       integer(kind=4) itiron,iron_id,lirongrd,iron_ncycle,
     &  iron_rec,iron_tid
       common/irondat/ironp,iron_time,iron_cycle
       common/irondat1/itiron,iron_id,lirongrd,
     &  iron_ncycle,iron_rec,iron_tid
      real(kind=8) srflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) srflx_dailyavg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer(kind=4), parameter :: n_srflx_day = 24
      integer(kind=4) :: num_srflx_day
      integer(kind=4), parameter :: MAX_NUM_SRFLX_DAY = 144
      integer(kind=4) iptr_srflx_day(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer(kind=4) iptr_srflx_day_set(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer(kind=4) PARinc_rst_read
      real(kind=8) srflx_day(-1:Lm+2+padd_X,-1:Mm+2+padd_E,
     &                         MAX_NUM_SRFLX_DAY)
      real(kind=8) srflx_frac
      common /forces_srflx/srflx
     &       , srflx_dailyavg
     &     , srflx_day, srflx_frac
      common /i_forces_srflx/ iptr_srflx_day, PARinc_rst_read,
     &     iptr_srflx_day_set,  num_srflx_day
      character(len=100) pco2_atm_file
      common /pco2_atm_file/ pco2_atm_file
      real(kind=8) visc2_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) visc2_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      real(kind=8) diff2(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff2/diff2
      real(kind=8) Akv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real(kind=8) Akt(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N,NT)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real(kind=8) bvf(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_bvf/ bvf
      real(kind=8) hbls(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /kpp_hbl/hbls
      real(kind=8) swr_frac(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /kpp_swr_frac/swr_frac
      real(kind=8) ghat(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /kpp_ghat/ghat
      real(kind=8) hbbls(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /kpp_hbbl/hbbls
      real(kind=8) weight(2,288)
      common /coup_weight/ weight
      real(kind=8) rufrc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) rvfrc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /coup_rufrc/rufrc /coup_rvfrc/rvfrc
      real(kind=8) rhoA(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) rhoS(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /coup_rhoA/rhoA /coup_rhoS/rhoS
      real(kind=8) r_D(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /coup_r_D/r_D
      real(kind=8) Zt_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) DU_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) DV_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) DU_avg2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) DV_avg2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /coup_Zt_avg1/Zt_avg1
     &       /coup_DU_avg1/DU_avg1 /coup_DV_avg1/DV_avg1
     &       /coup_DU_avg2/DU_avg2 /coup_DV_avg2/DV_avg2
      real(kind=4) cpu_time(4)
      real(kind=8) WallClock, time, tdays
      integer(kind=4) proc(2), numthreads, iic, kstp, knew
     &                           , iif, nstp, nnew, nrhs
     &                           , priv_count(16)
      logical synchro_flag, diag_sync
      common /priv_scalars/  WallClock, cpu_time,   proc,
     &         time, tdays, numthreads, iic,  kstp, knew
     &                           , iif, nstp, nnew, nrhs
     &       , priv_count, synchro_flag, diag_sync
C$OMP THREADPRIVATE(/priv_scalars/)
      real(kind=8) start_time, dt, dtfast, time_avg, xl,el, 
     &                          rdrg,rdrg2,Zob,
     &                                                 visc2,gamma2
      common /scalars_main/ start_time, dt, dtfast, time_avg, xl,el,
     &                                 rdrg,rdrg2,Zob, visc2,gamma2
      real(kind=8) rho0, tnu2(NT)
      common /scalars_main/ rho0, tnu2
      real(kind=8) v_sponge
      common /scalars_main/ v_sponge
      real(kind=8) tauM2_in, tauM2_out, attnM2
      common /scalars_main/ tauM2_in, tauM2_out, attnM2
      real(kind=8) tauM3_in, tauM3_out,  tauT_in, tauT_out
      common /scalars_main/ tauM3_in,tauM3_out, tauT_in,tauT_out
      real(kind=8) dSdt,dSdh
      common /scalars_sss/ dSdt,dSdh
      integer(kind=4) ntstart, ntimes, ndtfast, nfast, ninfo, 
     &                           may_day_flag,
     &                                                barr_count(16)
      common /scalars_main/ ntstart, ntimes, ndtfast, nfast, ninfo,
     &                               may_day_flag,    barr_count
      integer(kind=4) forw_start
      common /scalars_main/ forw_start
      real(kind=8), parameter :: pi=3.14159265358979323_8, 
     &                        Eradius=6371315._8,
     &              deg2rad=pi/180._8, rad2deg=180._8/pi, 
     &                         day2sec=86400._8,
     &                   sec2day=1._8/86400._8, Cp=3985._8, 
     &                           vonKar=0.41_8
     &                 , g=9.81_8
      real(kind=8) nmol_cm2_to_mmol_m2
      parameter (nmol_cm2_to_mmol_m2 = 0.01_8)
      integer(kind=4) imin,imax,jmin,jmax
      if (istr.eq.iwest .and. .not.west_exchng) then
        imin=istr
      else
        imin=istr-1
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        imax=iend
      else
        imax=iend+1
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        jmax=jend
      else
        jmax=jend+1
      endif
      Cg=Cstar * vonKar * (c_s*vonKar*epssfc)**(1._8/3._8)
      Vtc=Cv * sqrt(-betaT/(c_s*epssfc)) / (Ricr*vonKar**2)
      call alfabeta_tile (istr,iend,jstr,jend, imin,imax,
     &                             jmin,jmax, Bosol,Bo)
      do j=jmin,jmax
        do i=imin,imax
          Bo(i,j)=g*( Bosol(i,j)*(stflx(i,j,itemp)-srflx(i,j))
     &                              -Bo(i,j)*stflx(i,j,isalt)
     &                                                        )
          Bosol(i,j)=g*Bosol(i,j)*srflx(i,j)
          ustar(i,j)=sqrt(sqrt( sustr(i,j)**2+svstr(i,j)**2) )
          hbl(i,j)=hbls(i,j,nstp)
          kbl(i)=0
          Cr(i,N)=0._8
          Cr(i,0)=0._8
          FC(i,N)=0._8
        enddo
        do k=N-1,1,-1
          do i=imin,imax
            zscale=z_w(i,j,N)-z_r(i,j,k)
            Kern=zscale/(zscale+epssfc*hbl(i,j))
            FC(i,k)=FC(i,k+1) + Kern*(
     &                 0.5_8*( ( u(i  ,j,k+1,nstp)-u(i  ,j,k,nstp)
     &                        +u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp) )**2
     &                      +( v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)
     &                        +v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp) )**2
     &                      )/(Hz(i,j,k)+Hz(i,j,k+1))
     &               -0.5_8*(Hz(i,j,k)+Hz(i,j,k+1))*( Ri_inv*bvf(i,j,k)
     &                                            +C_Ek*f(i,j)*f(i,j)
     &                                                             ))
          enddo
        enddo
        do k=N,1,-1
          do i=imin,imax
            wrk1(i,k)=sqrt(swr_frac(i,j,k)*swr_frac(i,j,k-1))
            zscale=z_w(i,j,N)-z_r(i,j,k)
            Bfsfc=Bo(i,j)+Bosol(i,j)*(1._8-wrk1(i,k))
          if (Bfsfc .lt. 0._8) zscale=min(zscale, hbl(i,j)*epssfc)
          zscale=zscale*rmask(i,j)
          zetahat=vonKar*zscale*Bfsfc
          ustar3=ustar(i,j)**3
          if (zetahat .ge. 0._8) then
            ws=vonKar*ustar(i,j)*ustar3/max(ustar3+5._8*zetahat, 1.D-20)
          elseif (zetahat .gt. zeta_s*ustar3) then
            ws=vonKar*( (ustar3-16._8*zetahat)/ustar(i,j) )**r2
          else
            ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
          endif
            Vtsq=Vtc*ws*sqrt(max(0._8, bvf(i,j,k-1) ))
            Cr(i,k)=FC(i,k)+Vtsq
            if (kbl(i).eq.0 .and. Cr(i,k).lt.0._8) kbl(i)=k
          enddo
        enddo
        do i=imin,imax
          hbl(i,j)=z_w(i,j,N)-z_w(i,j,0) +eps
          if (kbl(i).gt.0) then
            k=kbl(i)
            if (k.eq.N) then
              hbl(i,j)=z_w(i,j,N)-z_r(i,j,N)
            else
              hbl(i,j)=z_w(i,j,N)-( z_r(i,j,k)*Cr(i,k+1)
     &                              -z_r(i,j,k+1)*Cr(i,k)
     &                              )/(Cr(i,k+1)-Cr(i,k))
            endif
          endif
          hbl(i,j)=hbl(i,j)*rmask(i,j)
        enddo
        do i=imin,imax
          kbl(i)=0
          Cr(i,0)=0.D0
          FC(i,0)=1.5D0*FC(i,1)-0.5D0*FC(i,2)
        enddo
        do k=1,N,+1
          do i=imin,imax
            Cr(i,k)=FC(i,k)-FC(i,0)
            if (kbl(i).eq.0 .and. Cr(i,k).gt.0._8) kbl(i)=k
          enddo
        enddo
        do i=imin,imax
          hbbl(i,j)=z_w(i,j,N)-z_w(i,j,0)
          if (kbl(i).gt.0) then
            k=kbl(i)
            if (k.eq.1) then
              hbbl(i,j)=z_r(i,j,1)-z_w(i,j,0)
            else
              hbbl(i,j)=( z_r(i,j,k-1)*Cr(i,k)-z_r(i,j,k)*Cr(i,k-1)
     &                            )/(Cr(i,k)-Cr(i,k-1)) -z_w(i,j,0)
            endif
          endif
          hbbl(i,j)=hbbl(i,j)*rmask(i,j)
        enddo
      enddo
      cff=1.D0/12.D0
      cff1=3.D0/16.D0
      if (istr.eq.iwest .and. .not.west_exchng) then
        do j=jmin,jmax
          hbl(istr-1,j)=hbl(istr,j)
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        do j=jmin,jmax
          hbl(iend+1,j)=hbl(iend,j)
        enddo
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        do i=imin,imax
          hbl(i,jstr-1)=hbl(i,jstr)
        enddo
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        do i=imin,imax
          hbl(i,jend+1)=hbl(i,jend)
        enddo
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        hbl(istr-1,jstr-1)=hbl(istr,jstr)
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        hbl(istr-1,jend+1)=hbl(istr,jend)
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        hbl(iend+1,jstr-1)=hbl(iend,jstr)
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        hbl(iend+1,jend+1)=hbl(iend,jend)
      endif
      do j=jstr-1,jend+1
        do i=istr,iend+1
          FX(i,j)=(hbl(i,j)-hbl(i-1,j))
     &                      *umask(i,j)
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr-1,iend+1
          FE1(i,j)=(hbl(i,j)-hbl(i,j-1))
     &                      *vmask(i,j)
        enddo
        do i=istr,iend
          FE(i,j)=FE1(i,j) + cff*( FX(i+1,j)+FX(i  ,j-1)
     &                            -FX(i  ,j)-FX(i+1,j-1))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend+1
          FX(i,j)=FX(i,j) + cff*( FE1(i,j+1)+FE1(i-1,j  )
     &                           -FE1(i,j  )-FE1(i-1,j+1))
        enddo
        do i=istr,iend
          hbl(i,j)=hbl(i,j) + cff1*( FX(i+1,j)-FX(i,j)
     &                              +FE(i,j+1)-FE(i,j))
          hbl(i,j)=hbl(i,j)*rmask(i,j)
        enddo
      enddo
      cff=1.D0/12.D0
      cff1=3.D0/16.D0
      if (istr.eq.iwest .and. .not.west_exchng) then
        do j=jmin,jmax
          hbbl(istr-1,j)=hbbl(istr,j)
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        do j=jmin,jmax
          hbbl(iend+1,j)=hbbl(iend,j)
        enddo
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        do i=imin,imax
          hbbl(i,jstr-1)=hbbl(i,jstr)
        enddo
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        do i=imin,imax
          hbbl(i,jend+1)=hbbl(i,jend)
        enddo
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        hbbl(istr-1,jstr-1)=hbbl(istr,jstr)
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        hbbl(istr-1,jend+1)=hbbl(istr,jend)
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        hbbl(iend+1,jstr-1)=hbbl(iend,jstr)
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        hbbl(iend+1,jend+1)=hbbl(iend,jend)
      endif
      do j=jstr-1,jend+1
        do i=istr,iend+1
          FX(i,j)=(hbbl(i,j)-hbbl(i-1,j))
     &                      *umask(i,j)
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr-1,iend+1
          FE1(i,j)=(hbbl(i,j)-hbbl(i,j-1))
     &                      *vmask(i,j)
        enddo
        do i=istr,iend
          FE(i,j)=FE1(i,j) + cff*( FX(i+1,j)+FX(i  ,j-1)
     &                            -FX(i  ,j)-FX(i+1,j-1))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend+1
          FX(i,j)=FX(i,j) + cff*( FE1(i,j+1)+FE1(i-1,j  )
     &                           -FE1(i,j  )-FE1(i-1,j+1))
        enddo
        do i=istr,iend
          hbbl(i,j)=hbbl(i,j) + cff1*( FX(i+1,j)-FX(i,j)
     &                              +FE(i,j+1)-FE(i,j))
          hbbl(i,j)=hbbl(i,j)*rmask(i,j)
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          cff=z_w(i,j,N)-z_w(i,j,0)
          if (hbl(i,j)+hbbl(i,j).gt.cff) then
            hbl(i,j)  = cff
            hbbl(i,j) = cff
          endif
        enddo
        do i=istr,iend
          kbl(i)=N
        enddo
        do k=N-1,1,-1
          do i=istr,iend
            if (z_w(i,j,k) .gt. z_w(i,j,N)-hbl(i,j)) kbl(i)=k
          enddo
        enddo
        do i=istr,iend
          k=kbl(i)
          z_bl=z_w(i,j,N)-hbl(i,j)
          zscale=hbl(i,j)
          if (swr_frac(i,j,k-1).gt. 0._8) then
            Bfsfc=Bo(i,j) +Bosol(i,j)*( 1._8 -swr_frac(i,j,k-1)
     &              *swr_frac(i,j,k)*(z_w(i,j,k)-z_w(i,j,k-1))
     &               /( swr_frac(i,j,k  )*(z_w(i,j,k)   -z_bl)
     &                 +swr_frac(i,j,k-1)*(z_bl -z_w(i,j,k-1))
     &                                                      ))
          else
            Bfsfc=Bo(i,j)+Bosol(i,j)
          endif
            if (Bfsfc.lt.0._8) zscale=min(zscale, hbl(i,j)*epssfc)
            zscale=zscale*rmask(i,j)
            zetahat=vonKar*zscale*Bfsfc
            ustar3=ustar(i,j)**3
            if (zetahat .ge. 0._8) then
              wm=vonKar*ustar(i,j)*ustar3/max( ustar3+5._8*zetahat,
     &                                                   1.D-20 )
              ws=wm
            else
              if (zetahat .gt. zeta_m*ustar3) then
                wm=vonKar*( ustar(i,j)*(ustar3-16._8*zetahat) )**r4
              else
                wm=vonKar*(a_m*ustar3-c_m*zetahat)**r3
              endif
              if (zetahat .gt. zeta_s*ustar3) then
                ws=vonKar*( (ustar3-16._8*zetahat)/ustar(i,j) )**r2
              else
                ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
              endif
            endif
          f1=5.0_8 * max(0._8, Bfsfc) * vonKar/(ustar(i,j)**4+eps)
          cff=1._8/(z_w(i,j,k)-z_w(i,j,k-1))
          cff_up=cff*(z_bl -z_w(i,j,k-1))
          cff_dn=cff*(z_w(i,j,k)   -z_bl)
          Av_bl=cff_up*Kv(i,j,k)+cff_dn*Kv(i,j,k-1)
          dAv_bl=cff * (Kv(i,j,k)  -   Kv(i,j,k-1))
          Gm1(i)=Av_bl/(hbl(i,j)*wm+eps)
          dGm1dS(i)=min(0._8, Av_bl*f1-dAv_bl/(wm+eps))
          At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
          dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
          Gt1(i)=At_bl/(hbl(i,j)*ws+eps)
          dGt1dS(i)=min(0._8, At_bl*f1-dAt_bl/(ws+eps))
          As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
          dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
          Gs1(i)=As_bl/(hbl(i,j)*ws+eps)
          dGs1dS(i)=min(0._8, As_bl*f1-dAs_bl/(ws+eps))
          Bfsfc_bl(i)=Bfsfc
        enddo
        do i=istr,iend
          do k=N-1,kbl(i),-1
            Bfsfc=Bfsfc_bl(i)
            zscale=z_w(i,j,N)-z_w(i,j,k)
            if (Bfsfc.lt.0._8) zscale=min(zscale, hbl(i,j)*epssfc)
            zscale=zscale*rmask(i,j)
            zetahat=vonKar*zscale*Bfsfc
            ustar3=ustar(i,j)**3
            if (zetahat .ge. 0._8) then
              wm=vonKar*ustar(i,j)*ustar3/max( ustar3+5._8*zetahat,
     &                                                   1.D-20 )
              ws=wm
            else
              if (zetahat .gt. zeta_m*ustar3) then
                wm=vonKar*( ustar(i,j)*(ustar3-16._8*zetahat) )**r4
              else
                wm=vonKar*(a_m*ustar3-c_m*zetahat)**r3
              endif
              if (zetahat .gt. zeta_s*ustar3) then
                ws=vonKar*( (ustar3-16._8*zetahat)/ustar(i,j) )**r2
              else
                ws=vonKar*(a_s*ustar3-c_s*zetahat)**r3
              endif
            endif
            sigma=(z_w(i,j,N)-z_w(i,j,k))/max(hbl(i,j),eps)
            a1=sigma-2._8
            a2=3._8-2._8*sigma
            a3=sigma-1._8
            if (sigma.lt.0.07D0) then
              cff=0.5_8*(sigma-0.07D0)**2/0.07D0
            else
              cff=0.D0
            endif
            Kv(i,j,k)=wm*hbl(i,j)*( cff + sigma*( 1._8+sigma*(
     &                           a1+a2*Gm1(i)+a3*dGm1dS(i) )))
            Kt(i,j,k)=ws*hbl(i,j)*( cff + sigma*( 1._8+sigma*(
     &                           a1+a2*Gt1(i)+a3*dGt1dS(i) )))
            Ks(i,j,k)=ws*hbl(i,j)*( cff + sigma*( 1._8+sigma*(
     &                           a1+a2*Gs1(i)+a3*dGs1dS(i) )))
            if (Bfsfc .lt. 0._8) then
              ghat(i,j,k)=Cg * sigma*(1._8-sigma)**2
            else
              ghat(i,j,k)=0._8
            endif
          enddo
          do k=kbl(i)-1,1,-1
            ghat(i,j,k)=0._8
            if (bvf(i,j,k).lt.0._8) then
              Kv(i,j,k)=Kv(i,j,k) + nu0c
              Kt(i,j,k)=Kt(i,j,k) + nu0c
              Ks(i,j,k)=Ks(i,j,k) + nu0c
            endif
          enddo
        enddo
        do i=istr,iend
          kbl(i)=N
        enddo
        do k=N-1,1,-1
          do i=istr,iend
            if (z_r(i,j,k)-z_w(i,j,0).gt.hbbl(i,j)) kbl(i)=k
          enddo
        enddo
        do i=istr,iend
          wm=vonKar*sqrt( r_D(i,j)*sqrt( 0.333333333333_8*(
     &                u(i,j,1,nstp)**2 +u(i+1,j,1,nstp)**2
     &                       +u(i,j,1,nstp)*u(i+1,j,1,nstp)
     &               +v(i,j,1,nstp)**2 +v(i,j+1,1,nstp)**2
     &                       +v(i,j,1,nstp)*v(i,j+1,1,nstp)
     &                                                ) ) )
          ws=wm
          k=kbl(i)
          z_bl=z_w(i,j,0)+hbbl(i,j)
          if (z_bl.lt.z_w(i,j,k-1)) k=k-1
          cff=1._8/(z_w(i,j,k)-z_w(i,j,k-1))
          cff_up=cff*(z_bl -z_w(i,j,k-1))
          cff_dn=cff*(z_w(i,j,k)   -z_bl)
          Av_bl=cff_up*Kv(i,j,k)+cff_dn*Kv(i,j,k-1)
          dAv_bl=cff * (Kv(i,j,k)  -   Kv(i,j,k-1))
          Gm1(i)=Av_bl/(hbbl(i,j)*wm+eps)
          dGm1dS(i)=min(0._8, -dAv_bl/(wm+eps))
          At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
          dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
          Gt1(i)=At_bl/(hbbl(i,j)*ws+eps)
          dGt1dS(i)=min(0._8, -dAt_bl/(ws+eps))
          As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
          dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
          Gs1(i)=As_bl/(hbbl(i,j)*ws+eps)
          dGs1dS(i)=min(0._8, -dAs_bl/(ws+eps))
          do k=1,N
            if (k.lt.kbl(i)) then
              sigma=min((z_w(i,j,k)-z_w(i,j,0)+Zob)
     &                       /(hbbl(i,j)+Zob),1._8)
              a1=sigma-2._8
              a2=3._8-2._8*sigma
              a3=sigma-1._8
              Kv0=wm*hbbl(i,j)*( sigma*( 1._8+sigma*(
     &                            a1+a2*Gm1(i)+a3*dGm1dS(i) )))
              Kt0=ws*hbbl(i,j)*( sigma*( 1._8+sigma*(
     &                            a1+a2*Gt1(i)+a3*dGt1dS(i) )))
              Ks0=ws*hbbl(i,j)*( sigma*( 1._8+sigma*(
     &                            a1+a2*Gs1(i)+a3*dGs1dS(i) )))
              z_bl=z_w(i,j,N)-hbl(i,j)
              if (z_w(i,j,k).gt.z_bl) then
                Kv0=max(Kv(i,j,k),Kv0)
                Kt0=max(Kt(i,j,k),Kt0)
                Ks0=max(Ks(i,j,k),Ks0)
              endif
              Kv(i,j,k)=Kv0
              Kt(i,j,k)=Kt0
              Ks(i,j,k)=Ks0
            else
              if (bvf(i,j,k).lt.0._8) then
                z_bl=z_w(i,j,N)-hbl(i,j)
                if (z_w(i,j,k).lt.z_bl) then
                  Kv(i,j,k)=Kv(i,j,k) + nu0c
                  Kt(i,j,k)=Kt(i,j,k) + nu0c
                  Ks(i,j,k)=Ks(i,j,k) + nu0c
                endif
              endif
            endif
          enddo
        enddo
        do i=istr,iend
          if (rmask(i,j).gt.0.5_8) then
            do k=1,N
              Akv(i,j,k)=Kv(i,j,k)
              Akt(i,j,k,itemp)=Kt(i,j,k)
              Akt(i,j,k,isalt)=Ks(i,j,k)
            enddo
          else
            do k=1,N
              Akv(i,j,k)=0._8
              Akt(i,j,k,itemp)=0._8
              Akt(i,j,k,isalt)=0._8
            enddo
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          hbls(i,j,3-nstp)=hbl(i,j)
        enddo
      enddo
      if (istr.eq.iwest .and. .not.west_exchng) then
        do j=jstr,jend
          hbls(istr-1,j,3-nstp)=hbls(istr,j,3-nstp)
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        do j=jstr,jend
          hbls(iend+1,j,3-nstp)=hbls(iend,j,3-nstp)
        enddo
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        do i=istr,iend
          hbls(i,jstr-1,3-nstp)=hbls(i,jstr,3-nstp)
        enddo
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        do i=istr,iend
          hbls(i,jend+1,3-nstp)=hbls(i,jend,3-nstp)
        enddo
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        hbls(istr-1,jstr-1,3-nstp)=hbls(istr,jstr,3-nstp)
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        hbls(istr-1,jend+1,3-nstp)=hbls(istr,jend,3-nstp)
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        hbls(iend+1,jstr-1,3-nstp)=hbls(iend,jstr,3-nstp)
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        hbls(iend+1,jend+1,3-nstp)=hbls(iend,jend,3-nstp)
      endif
      do j=jstr,jend
        do i=istr,iend
          hbbls(i,j,3-nstp)=hbbl(i,j)
        enddo
      enddo
      if (istr.eq.iwest .and. .not.west_exchng) then
        do j=jstr,jend
          hbbls(istr-1,j,3-nstp)=hbbls(istr,j,3-nstp)
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        do j=jstr,jend
          hbbls(iend+1,j,3-nstp)=hbbls(iend,j,3-nstp)
        enddo
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        do i=istr,iend
          hbbls(i,jstr-1,3-nstp)=hbbls(i,jstr,3-nstp)
        enddo
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        do i=istr,iend
          hbbls(i,jend+1,3-nstp)=hbbls(i,jend,3-nstp)
        enddo
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        hbbls(istr-1,jstr-1,3-nstp)=hbbls(istr,jstr,3-nstp)
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        hbbls(istr-1,jend+1,3-nstp)=hbbls(istr,jend,3-nstp)
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        hbbls(iend+1,jstr-1,3-nstp)=hbbls(iend,jstr,3-nstp)
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        hbbls(iend+1,jend+1,3-nstp)=hbbls(iend,jend,3-nstp)
      endif
      call exchange_4_tile ( istr,iend,jstr,jend,
     &                 Akt(-1,-1,0,isalt), N+1,
     &       Akv, N+1,  Akt(-1,-1,0,itemp), N+1,
     &                    hbls(-1,-1,3-nstp), 1
     &                                                 )
      return
      end
      subroutine check_kpp_switches (ierr)
      implicit none
      integer(kind=4) ierr, is,ie, lenstr
      integer(kind=4), parameter ::
     &               LLm=435, MMm=660, N=60
      integer(kind=4), parameter ::
     &      NP_XI=8, NP_ETA=32, NSUB_X=1, NSUB_E=1
      integer(kind=4), parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA
      integer(kind=4) ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
      logical west_exchng,  east_exchng
      logical south_exchng, north_exchng
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
     &                , west_exchng,  east_exchng
     &                , south_exchng, north_exchng
      integer(kind=4), parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
     &       , itemp=1
     &       , isalt=2
     &       , ntrc_salt=1
      integer(kind=4), parameter :: ntrc_pas=0
      integer(kind=4), parameter :: itrc_bio=itemp+ntrc_salt+ntrc_pas+1
      integer(kind=4), parameter :: iPO4=itrc_bio,iNO3=iPO4+1, 
     &                           iSIO3=iPO4+2,
     &     iNH4=iPO4+3,
     &     iFE=iPO4+4, iO2=iPO4+5, iDIC=iPO4+6,
     &     iALK=iPO4+7, iDOC=iPO4+8, iSPC=iPO4+9,
     &     iSPCHL=iPO4+10, iSPCACO3=iPO4+11, iDIATC=iPO4+12,
     &     iDIATCHL=iPO4+13, iZOOC=iPO4+14, iSPFE=iPO4+15,
     &     iDIATSI=iPO4+16, iDIATFE=iPO4+17, iDIAZC=iPO4+18,
     &     iDIAZCHL=iPO4+19, iDIAZFE=iPO4+20, iDON=iPO4+21,
     &     iDOFE=iPO4+22, iDOP=iPO4+23
      integer(kind=4), parameter :: iNO2 = iDOP + 1
      integer(kind=4), parameter :: iN2O = iNO2 + 1, iN2 = iN2O + 1
      integer(kind=4), parameter :: iSedOrgC = 1
      integer(kind=4), parameter :: iSedCaCO3 = iSedOrgC + 1
      integer(kind=4), parameter :: iSedSi = iSedCaCO3 + 1
      integer(kind=4), parameter :: NT_sed = 3
      integer(kind=4), parameter :: ntrc_bio = 24
     &     + 3
      integer(kind=4), parameter :: NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio
      integer(kind=4), parameter :: max_opt_size=2048
      character(len=max_opt_size) cpps, srcs, kwds
      common /strings/ cpps, srcs, kwds
      ie=lenstr(cpps)
      is=ie+2
      ie=is+10
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='<lmd_kpp.F>'
      is=ie+2
      ie=is+16
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='INT_AT_RHO_POINTS'
      is=ie+2
      ie=is+9
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='SMOOTH_HBL'
      is=ie+2
      ie=is+18
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='LIMIT_UNSTABLE_ONLY'
      is=ie+2
      ie=is+13
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='MERGE_OVERLAP'
      return
  99  if (mynode.eq.0) write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      return
      end
