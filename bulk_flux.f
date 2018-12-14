      subroutine bulk_flux_tile (istr,iend,jstr,jend, ierr)
      implicit none
      integer(kind=4)  istr,iend,jstr,jend, ierr,   imin,imax,jmin,jmax
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
       real(kind=8) tracer(-1:Lm+2+padd_X,N,ntrc_bio,2)
        common /tracers/ tracer
        real(kind=8) tracer_sed(-1:Lm+2+padd_X,NT_sed,2)
        common /tracer_sed/ tracer_sed
        real(kind=8) ifrac(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &    press(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
        common /fic_ap/ifrac,press
        real(kind=8) PH_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   pCO2sw(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   pCO2air(-1:Lm+2+padd_X,-1:Mm+2+padd_E) ,
     &   PAR(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   PARinc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
     &   ,PARinc_rst(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       common /time_averaging/
     &       ph_hist,pCO2sw,pCO2air,
     &       PAR,PARinc
     &  , PARinc_rst
        real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E) :: 
     &                     dom_sp_sfc, dom_diat_sfc,
     &       dom_diaz_sfc, dom_sp_int, dom_diat_int, dom_diaz_int,
     &       spchl_int, diatchl_int, diazchl_int
        common /specdom/ dom_sp_sfc, dom_diat_sfc,
     &       dom_diaz_sfc, dom_sp_int, dom_diat_int, dom_diaz_int,
     &       spchl_int, diatchl_int, diazchl_int
       real(kind=8) WS_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   XKW_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   AP_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   SCHMIDT_O2_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   O2SAT_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   FG_O2_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &    SCHMIDT_CO2_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   CO2STAR_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   DCO2STAR_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &    FG_CO2_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   IRON_FLUX_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       real(kind=8) fg_n2o_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       real(kind=8) fg_n2_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
        real(kind=8)
     &    PO4_RESTORE_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    NO3_RESTORE_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    SiO3_RESTORE_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   PO4STAR_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    POC_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    POC_PROD_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   POC_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    CaCO3_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    CaCO3_PROD_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    CaCO3_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    SiO2_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    SiO2_PROD_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    SiO2_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    dust_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
          real(kind=8)  dust_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,
     &                                N),
     &    P_iron_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    P_iron_PROD_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    P_iron_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    graze_sp_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    graze_diat_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    graze_tot_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_loss_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_loss_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    zoo_loss_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_agg_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_agg_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
          real(kind=8)  photoC_sp_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    f_ratio_sp_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    photoC_diat_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    f_ratio_diat_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    tot_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    no3_v_sp_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    nh4_v_sp_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    no3_v_diat_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    nh4_v_diat_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DOC_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DOC_remin_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    Fe_scavenge_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_N_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_Fe_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &     sp_PO4_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_light_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_N_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &     diat_Fe_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_PO4_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_SiO3_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
         real(kind=8) diat_light_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,
     &                                N),
     &   CaCO3_form_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diaz_Nfix_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    graze_diaz_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diaz_loss_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
        real(kind=8) photoC_diaz_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diaz_P_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diaz_Fe_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &     diaz_light_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &     Fe_scavenge_rate_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DON_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DON_remin_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DOFe_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
         real(kind=8)  DOFe_remin_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   DOP_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DOP_remin_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    bSI_form_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    photoFe_diaz_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    photoFe_diat_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    photoFe_sp_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    nitrif_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
         real(kind=8) ammox_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &        denitr_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &        n2o_prod_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &        n2_prod_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &        denitr_sed_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       common /time_averaging1/WS_HIST, XKW_HIST,
     &   AP_HIST, SCHMIDT_O2_HIST, O2SAT_HIST, FG_O2_HIST,
     &    SCHMIDT_CO2_HIST, CO2STAR_HIST,
     &    DCO2STAR_HIST,
     &    FG_CO2_HIST, IRON_FLUX_HIST,
     &    PO4_RESTORE_HIST, NO3_RESTORE_HIST,
     &    SiO3_RESTORE_HIST, PO4STAR_HIST,
     &    POC_FLUX_IN_HIST, POC_PROD_HIST, POC_REMIN_HIST,
     &    CaCO3_FLUX_IN_HIST, CaCO3_PROD_HIST,
     &    CaCO3_REMIN_HIST,  SiO2_FLUX_IN_HIST,
     &    SiO2_PROD_HIST, SiO2_REMIN_HIST, dust_FLUX_IN_HIST,
     &    dust_REMIN_HIST, P_iron_FLUX_IN_HIST,
     &    P_iron_PROD_HIST, P_iron_REMIN_HIST,
     &    graze_sp_HIST, graze_diat_HIST, graze_tot_HIST,
     &    sp_loss_HIST, diat_loss_HIST, zoo_loss_HIST,
     &    sp_agg_HIST, diat_agg_HIST,
     &    photoC_sp_HIST, f_ratio_sp_hist,
     &    photoC_diat_HIST, f_ratio_diat_hist, tot_prod_HIST,
     &    no3_v_sp_hist, nh4_v_sp_hist,
     &    no3_v_diat_hist, nh4_v_diat_hist,
     &    DOC_prod_HIST, DOC_remin_HIST, Fe_scavenge_HIST
     &        , fg_n2o_hist
     &        , fg_n2_hist
     &        , denitr_sed_hist
       common /time_averaging2/
     &    sp_N_lim_HIST, sp_Fe_lim_HIST, sp_PO4_lim_HIST,
     &    sp_light_lim_HIST, diat_N_lim_HIST, diat_Fe_lim_HIST,
     &    diat_PO4_lim_HIST, diat_SiO3_lim_HIST,
     &    diat_light_lim_HIST, CaCO3_form_HIST,
     &    diaz_Nfix_HIST, graze_diaz_HIST, diaz_loss_HIST,
     &     photoC_diaz_HIST, diaz_P_lim_HIST,
     &    diaz_Fe_lim_HIST, diaz_light_lim_HIST,
     &     Fe_scavenge_rate_HIST, DON_prod_HIST,
     &    DON_remin_HIST, DOFe_prod_HIST,
     &    DOFe_remin_HIST, DOP_prod_HIST,
     &    DOP_remin_HIST, bSI_form_HIST,
     &    photoFe_diaz_HIST, photoFe_diat_HIST, photoFe_sp_HIST,
     &    nitrif_HIST
     &    , ammox_hist, denitr_hist, n2o_prod_hist, n2_prod_hist
       real(kind=8) bot_flux_poc_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &      bot_flux_caco3_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &      bot_flux_si_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &      bot_flux_fe_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       common /sed_flux/ bot_flux_poc_hist, bot_flux_caco3_hist,
     &      bot_flux_si_hist,bot_flux_fe_hist
       integer(kind=4)  po4_ind , no3_ind,sio3_ind, nh4_ind,fe_ind,
     & o2_ind, dic_ind,alk_ind,doc_ind,spC_ind,spChl_ind,
     & spCaCO3_ind,diatC_ind,diatChl_ind,zooC_ind,spFe_ind,
     &  diatSi_ind,diatFe_ind,diazC_ind,diazChl_ind, diazFe_ind,
     &  don_ind,dofe_ind,dop_ind
       parameter (po4_ind=1 , no3_ind=2,sio3_ind=3, nh4_ind=4,
     &  fe_ind=5,o2_ind=6, dic_ind=7,alk_ind=8,doc_ind=9,
     &  spC_ind=10,spChl_ind=11, spCaCO3_ind=12,diatC_ind=13,
     &  diatChl_ind=14,zooC_ind=15,spFe_ind=16,
     &  diatSi_ind=17,diatFe_ind=18,diazC_ind=19,
     &  diazChl_ind=20, diazFe_ind=21,
     &  don_ind=22,dofe_ind=23,dop_ind=24)
       integer(kind=4), parameter :: no2_ind = dop_ind + 1
       integer(kind=4), parameter :: n2o_ind = no2_ind + 1
       integer(kind=4), parameter :: n2_ind = n2o_ind + 1
       integer(kind=4), parameter :: sed_poc_ind = 1
       integer(kind=4), parameter :: sed_caco3_ind = 2
       integer(kind=4), parameter :: sed_si_ind = 3
       logical lsource_sink,lflux_gas_o2, lflux_gas_co2,
     &  liron_flux,ldust_flux
        common /ecoflag/lsource_sink,lflux_gas_o2,lflux_gas_co2,
     &   liron_flux,ldust_flux
       logical lrest_po4,lrest_no3,lrest_sio3
       real(kind=8) po4_clim(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   no3_clim(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   sio3_clim(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
           real(kind=8) nutr_rest_time_inv(N)
        common /restore_flag/lrest_po4,lrest_no3,lrest_sio3
        common /restore_clim/po4_clim,
     &      no3_clim,sio3_clim,nutr_rest_time_inv
      real(kind=8) t_sed(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT_sed)
      common /ocean_t_sed/t_sed
        real(kind=8) sinking_particle_POC(6,-1:Lm+2+padd_X,
     &                          -1:Mm+2+padd_E),
     & sinking_particle_P_CaCO3(6,-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     & sinking_particle_P_sio2(6,-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     & sinking_particle_dust(6,-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     & sinking_particle_P_iron(6,-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       real(kind=8) diss(5),gamma(5),mass(5),rhoo(5)
      common /sinking_part/sinking_particle_POC,
     &  sinking_particle_P_CaCO3,sinking_particle_P_SiO2,
     &  sinking_particle_dust,sinking_particle_P_iron,
     &  diss,gamma,mass,rhoo
        logical landmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
        common /calcation/landmask
      real(kind=8)    tol, ustarold, q10, t10, cff, cff1,cff2,cff3
      real(kind=8)    ithv, iqhv, vmag, sst, ssq, al2, alz, ssu, ssv
      real(kind=8)    delt, delq, evap, lwup, hfsen, hflat
      real(kind=8)    rdn,ustar,qstar, tstar
      real(kind=8)    cph, hol, psixh, psimh, psixhu, psixh2
      real(kind=8)    rhn, ren, cdn, xqq, qsat, psimhu, cp_a
      real(kind=8)    Cd
       real(kind=8)   fwflx, cff7,cff8
       integer(kind=4) it7,it8,niter,  kiter, i,j
       real(kind=8), parameter ::
     &            umin   = 0.5_8,
     &            CtoK   = 273.16_8,
     &            cpvir  = 0.801640362327584_8,
     &            zvir   = 0.607792984078361_8,
     &           stebol  = 5.670400D-08,
     &            latvap = 2.5D6,
     &            irhoa  = 0.819672_8,
     &            rho_a  = 1.22_8,
     &            a1     = 0.000142_8,
     &            a2     = 0.0000764_8,
     &            a3     = 0.0027_8
      if (istr.eq.iwest .and. .not.west_exchng) then
        imin=istr-1
      else
        imin=istr-1
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        imax=iend+1
      else
        imax=iend+1
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        jmin=jstr-1
      else
        jmin=jstr-1
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        jmax=jend+1
      else
        jmax=jend+1
      endif
      it7=3-itsss
      it8=itsss
      cff7=sss_time(it8)-time
      cff8=time-sss_time(it7)
      cff=1._8/(cff7+cff8)
      cff7=cff7*cff
      cff8=cff8*cff
      al2   = log(zref/ztref)
      alz   = log(zref/10)
      do j=jmin,jmax
        do i=imin,imax
          if (rmask(i,j) .gt. 0.5_8 ) then
            ssu  = coeff_ssuv * 0.5_8*(u(i,j,N,nstp)+u(i+1,j,N,nstp))
            ssv  = coeff_ssuv * 0.5_8*(v(i,j,N,nstp)+v(i,j+1,N,nstp))
            vmag = max( umin,sqrt(
     &             (uwnd(i,j)-ssu)**2
     &            +(vwnd(i,j)-ssv)**2  ) )
            sst  = t(i,j,N,nstp,itemp) + CtoK
            ssq  = irhoA*0.98_8*640380._8*exp(-5107.4_8/sst)
            t10  = tair(i,j) + CtoK
            q10  = qair(i,j)
            delt = t10 - sst
            delq = q10 - ssq
            ren  = 0.0346_8
            rdn   = sqrt( a1 + a2*vmag + a3 / vmag )
            ustar = rdn * vmag
            if (delt .gt. 0) then
              tstar=0.018_8*delt
            else
              tstar=0.0327_8*delt
            endif
            qstar = ren*delq
            do niter = 1,5
              cff  = 1._8/ ( 1._8+ zvir*q10 )
              ithv = cff / t10
              iqhv = zvir*cff
              hol  = vonKar*g*zref*(tstar*ithv+qstar*iqhv)/ustar**2
              if (hol .gt. 0._8) then
                hol    = min(10._8, hol)
                psimh  = -5._8 * hol
                psixh  = -5._8 * hol
                psixh2 = 0.2_8*psixh
                rhn    = 0.018_8
              else
                hol    = max(-10._8, hol)
                xqq    = sqrt( sqrt(1._8-16._8*hol) )
                psimh  = psimhu(xqq)
                psixh  = psixhu(xqq)
                xqq    = sqrt( sqrt(1._8-3.2_8*hol) )
                psixh2 = 2._8*log( 0.5_8*(1._8+xqq*xqq) )
                rhn   = 0.0327_8
              endif
              cff = vonKar*vmag/sqrt( (vonKar + rdn*(alz-psimh))**2 )
              rdn = sqrt( a1 + a2*cff + a3 / cff )
              cff  = ( psixh - psixh2 - al2 - alz ) / vonKar
              t10  = tair(i,j) - cff*tstar + CtoK
              q10  = qair(i,j) - cff*qstar
              delt = t10 - sst
              delq = q10 - ssq
              ustar = vonKar * vmag * rdn / (vonKar + rdn*(alz-psimh))
              tstar = vonKar * delt * rhn / (vonKar + rhn*(alz-psixh))
              qstar = vonKar * delq * ren / (vonKar + ren*(alz-psixh))
            enddo
            cp_a       = 1004.64_8*(1._8 + cpvir*ssq)
            hfsen      =  cp_a*rho_a*ustar*tstar
            hflat      =  latvap*rho_a*ustar*qstar
            evap       = -rho_a*ustar*qstar
            lwup       = -stebol*sst**4
            fwflx      =  0.001_8*(evap-rain(i,j))
            dSdt=1._8/(90._8*day2sec)
              cff=1._8/(rho0*vmag)
              sustr(i,j) = cff*rho_a*ustar*ustar*(uwnd(i,j)-ssu)
              svstr(i,j) = cff*rho_a*ustar*ustar*(vwnd(i,j)-ssv)
              cff=1._8/(rho0*Cp)
              stflx(i,j,itemp) = cff*(
     &          radsw(i,j)+hflat+hfsen+radlw(i,j)+lwup
     &                                 )
              srflx(i,j)       = cff*radsw(i,j)
              iptr_srflx_day(i,j) = iptr_srflx_day(i,j) + 1
              if (iptr_srflx_day(i,j) .gt. num_srflx_day) then
                 iptr_srflx_day(i,j) = 1
                 iptr_srflx_day_set(i,j) = 1
              end if
              srflx_day(i,j,iptr_srflx_day(i,j)) = srflx(i,j)
              stflx(i,j,isalt) = fwflx*t(i,j,N,nstp,isalt)
     &              -dSdt*Hz(i,j,N)*( t(i,j,N,nstp,isalt)
     &              -(cff7*sssg(i,j,it7)+cff8*sssg(i,j,it8))
     &                                                          )
            shflx_net(i,j)= stflx(i,j,itemp)/cff
            shflx_lat(i,j)= hflat
            shflx_sen(i,j)= hfsen
            shflx_rad(i,j)= radsw(i,j)
            shflx_wwk(i,j)= sustr(i,j) *  0.5_8
     &                    * (u(i,j,N,nstp)+u(i+1,j,N,nstp))
     &                    + svstr(i,j) *  0.5_8
     &                    * (v(i,j,N,nstp)+v(i,j+1,N,nstp))
            if (i .lt. imax) then
            surf_u(i,j)   = 0.5_8 * (u(i,j,N,nstp)+u(i+1,j,N,nstp))
            else
               surf_u(i,j)   = u(i,j,N,nstp)
            end if
            if (j .lt. jmax) then
            surf_v(i,j)   = 0.5_8 * (v(i,j,N,nstp)+v(i,j+1,N,nstp))
            else
               surf_v(i,j)   = v(i,j,N,nstp)
            end if
            shflx_wwk(i,j)= sustr(i,j) * surf_u(i,j)
     &                    + svstr(i,j) * surf_v(i,j)
            swflx_emp(i,j)=(  stflx(i,j,isalt)
     &                    /   max(t(i,j,N,nstp,isalt),1._8) )*8640000._8
            sustr_blk(i,j)= sustr(i,j)*rho0
            svstr_blk(i,j)= svstr(i,j)*rho0
          else
            stflx(i,j,isalt)=0._8
            stflx(i,j,itemp)=0._8
            srflx(i,j)=0._8
            sustr(i,j) = 0._8
            svstr(i,j) = 0._8
            swflx_emp(i,j)=0._8
            shflx_rad(i,j)=0._8
            shflx_net(i,j)=0._8
            shflx_lat(i,j)=0._8
            shflx_sen(i,j)=0._8
            sustr_blk(i,j)=0._8
            svstr_blk(i,j)=0._8
            shflx_wwk(i,j)=0._8
            surf_u(i,j)=0._8
            surf_v(i,j)=0._8
          endif
        enddo
      enddo
      call det_srflx_dailyavg_tile(imin,imax,jmin,jmax)
      return
      end
      function psimhu (xd)
      real(kind=8) psimhu, xd
      psimhu=1.570796326794896619D0 -2.D0*atan(xd)
     &         +log(0.125D0*(1._8+xd*xd)*(1._8+xd)**2)
      return
      end
      function psixhu (xd)
      real(kind=8) psixhu, xd
      psixhu = 2._8*log( 0.5_8*(1._8+xd*xd) )
      return
      end
