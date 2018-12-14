      subroutine read_inp (ierr)
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
      real(kind=8) theta_s,theta_b, hc, Cs_w(0:N), Cs_r(N)
      common /scoord_vars/ theta_s,theta_b, hc, Cs_w,Cs_r
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
      integer(kind=4), parameter :: indxTime=1, indxZ=2, indxUb=3, 
     &                              indxVb=4
     &                    , indxU=5, indxV=6, indxO=7, indxW=8
     &                    , indxR=9, indxT=10
     &                    , indxS=indxT+1
       integer(kind=4), parameter :: indxPO4=indxT+ntrc_salt+ntrc_pas+1,
     &           indxNO3 =indxPO4+1, indxSio3=indxPO4+2,
     &           indxNH4 =indxPO4+3, indxFe=indxPO4+4,
     &           indxO2 =indxPO4+5, indxDic=indxPO4+6,
     &           indxAlk =indxPO4+7, indxDoc=indxPO4+8,
     &           indxSpc =indxPO4+9, indxSpchl=indxPO4+10,
     &           indxSpcaco3 =indxPO4+11, indxDiatc=indxPO4+12,
     &           indxDiatchl =indxPO4+13, indxZooc=indxPO4+14,
     &           indxSpfe =indxPO4+15, indxDiatsi=indxPO4+16,
     &           indxDiatfe =indxPO4+17, indxDiazc=indxPO4+18,
     &           indxDiazchl =indxPO4+19, indxDiazfe=indxPO4+20,
     &           indxDon =indxPO4+21, indxDofe=indxPO4+22,
     &           indxDop =indxPO4+23
       integer(kind=4), parameter :: indxNO2 = indxDOP + 1
      integer(kind=4), parameter :: indxN2O = indxNO2 + 1
      integer(kind=4), parameter :: indxN2 = indxN2O + 1
       integer(kind=4), parameter :: indxPH_rst = indxDOP+1
     &      + 3
       integer(kind=4), parameter :: indxPCO2_rst = indxPH_rst+1,
     &      indxPCO2air_rst = indxPCO2_rst+1,
     &      indxPARinc_rst = indxPCO2air_rst+1,
     &      indxPAR_rst = indxPARinc_rst+1
      integer(kind=4), parameter :: indxSedOrgC = indxPAR_rst + 1
      integer(kind=4), parameter :: indxSedCaCO3 = indxSedOrgC + 1
      integer(kind=4), parameter :: indxSedSi = indxSedCaCO3 + 1
      integer(kind=4), parameter :: indxSedFe = indxSedSi + 1
      integer(kind=4), parameter :: indxFreqDomSP_sfc = indxPAR_rst + 1
     &     + NT_sed
      integer(kind=4), parameter :: indxFreqDomDIAT_sfc = 
     &                       indxFreqDomSP_sfc + 1
      integer(kind=4), parameter :: indxFreqDomDIAZ_sfc = 
     &                       indxFreqDomDIAT_sfc+1
      integer(kind=4), parameter :: indxFreqDomSP_int = 
     &                      indxFreqDomDIAZ_sfc + 1
      integer(kind=4), parameter :: indxFreqDomDIAT_int = 
     &                       indxFreqDomSP_int + 1
      integer(kind=4), parameter :: indxFreqDomDIAZ_int = 
     &                       indxFreqDomDIAT_int+1
       integer(kind=4), parameter :: indxSedFirst = indxSedOrgC
       integer(kind=4), parameter :: indxAkv=indxT+NT
     &      + 5
     &     + NT_sed
     &      + 6
       integer(kind=4), parameter :: indxAkt=indxAkv+1
     &                    , indxAks=indxAkt+1
     &                    , indxHbls=indxAks+1
     &                    , indxHbbl=indxHbls+1
      integer(kind=4), parameter :: iaux=6
      integer(kind=4), parameter :: max_blk_file=8
      integer(kind=4) max_blk, ncidbulk(max_blk_file), nrst, ncidrst, 
     &                              nrecrst,
     &      nrrec, nrpfrst, ncidclm,nwrt, ncidhis, nrechis,
     &      nrpfhis
     &     , ntdust, ntiron
      common /ncvars/ max_blk, ncidbulk, nrst, ncidrst, nrecrst,
     &      nrrec, nrpfrst, ncidclm,nwrt, ncidhis, nrechis,
     &      nrpfhis
     &     , ntdust, ntiron
      integer(kind=4) hisSustr_blk,   hisSvstr_blk
     &      , hisShflx_net,   hisShflx_rad
     &      , hisSwflx_emp
     &      , hisShflx_lat,   hisShflx_sen
     &      , hisShflx_wwk
     &      , hisSurf_u, hisSurf_v
      common /ncvars/ hisSustr_blk,   hisSvstr_blk
     &      , hisShflx_net,   hisShflx_rad
     &      , hisSwflx_emp
     &      , hisShflx_lat,    hisShflx_sen
     &      , hisShflx_wwk
     &      , hisSurf_u, hisSurf_v
      integer(kind=4) ntsavg,  navg
      common /ncvars/ ntsavg,  navg
      integer(kind=4) rstTime, rstTstep,      rstZ,   rstUb,  rstVb,
     &        hisTime, hisTstep,      hisZ,   hisUb,  hisVb
      common /ncvars/
     &        rstTime, rstTstep,      rstZ,   rstUb,  rstVb,
     &        hisTime, hisTstep,      hisZ,   hisUb,  hisVb
      integer(kind=4) rst_DU_avg2, rst_DV_avg2
      common /ncvars/ rst_DU_avg2, rst_DV_avg2
      integer(kind=4) rstU, rstV, rstT(NT+1), hisO,   hisW,   hisR,
     &        hisU, hisV, hisT(NT+1), hisAkv, hisAkt, hisAks
     &      , rstPH, rstPCO2, rstPCO2air, rstPARinc, rstPAR
     &      , hisPH, hisPCO2, hisPCO2air, hisPARinc, hisPAR
     &      , avgPH, avgPCO2, avgPCO2air, avgPARinc, avgPAR
     &      , rstTstepFA
     &      , rstTsed(NT_sed), hisTsed(NT_sed)
      common /ncvars/
     &        rstU, rstV, rstT,       hisO,   hisW,   hisR,
     &        hisU, hisV, hisT,       hisAkv, hisAkt, hisAks
     &      , rstPH, rstPCO2, rstPCO2air, rstPARinc, rstPAR
     &      , hisPH, hisPCO2, hisPCO2air, hisPARinc, hisPAR
     &      , avgPH, avgPCO2, avgPCO2air, avgPARinc, avgPAR
     &      , rstTstepFA
     &      , rstTsed, hisTsed
      integer(kind=4) rstHbls, hisHbls
      common /ncvars/ rstHbls, hisHbls
      integer(kind=4) rstHbbl, hisHbbl
      common /ncvars/ rstHbbl, hisHbbl
      integer(kind=4) avgSustr_blk,   avgSvstr_blk
     &      , avgShflx_net,   avgShflx_rad
     &      , avgSwflx_emp
     &      , avgShflx_lat,   avgShflx_sen
     &      , avgShflx_wwk
     &      , avgSurf_u, avgSurf_v
      common /ncvars/ avgSustr_blk,   avgSvstr_blk
     &      , avgShflx_net,   avgShflx_rad
     &      , avgSwflx_emp
     &      , avgShflx_lat,   avgShflx_sen
     &      , avgShflx_wwk
     &      , avgSurf_u, avgSurf_v
      integer(kind=4) ncidavg, nrecavg, nrpfavg, avgTime, avgTstep, 
     &                            avgZ, avgUb,
     &     avgVb
      common /ncvars/  ncidavg, nrecavg, nrpfavg,
     &                                   avgTime, avgTstep, avgZ, avgUb,
     &     avgVb
      integer(kind=4) avgU,  avgV,  avgT(NT+1),  avgR,    avgO,    avgW,
     &                                   avgAkv,  avgAkt,  avgAks
     &      , avgTsed(NT_sed)
      common /ncvars/ avgU, avgV, avgT,  avgR,    avgO,    avgW,
     &                                   avgAkv,  avgAkt,  avgAks
     &      , avgTsed
      integer(kind=4) avgFreqDomSP_sfc, avgFreqDomDIAT_sfc,
     &     avgFreqDomDIAZ_sfc, avgFreqDomSP_int,
     &     avgFreqDomDIAT_int, avgFreqDomDIAZ_int
      common /ncvars/ avgFreqDomSP_sfc, avgFreqDomDIAT_sfc,
     &     avgFreqDomDIAZ_sfc, avgFreqDomSP_int,
     &     avgFreqDomDIAT_int, avgFreqDomDIAZ_int
      integer(kind=4) avgHbls
      common /ncvars/ avgHbls
      integer(kind=4) avgHbbl
      common /ncvars/ avgHbbl
      logical ldefhis, wrthis(100+NT)
      common /ncvars/ ldefhis, wrthis
      logical wrtavg(100+NT)
      common /ncvars/ wrtavg
      logical wrt_pfa(NT)
      common /ncvars/ wrt_pfa
      integer(kind=4), parameter :: r2dvar=0, u2dvar=1, v2dvar=2, 
     &                             p2dvar=3,
     &            r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8
      integer(kind=4) xi_rho, xi_u,   eta_rho, eta_v
      common /ncvars/ xi_rho, xi_u,   eta_rho, eta_v
      integer(kind=4), parameter :: max_name_size=64
      character date_str*44, title*80
      character(len=max_name_size) ininame, grdname,
     &                 hisname, rstname
     &    , blkfile(max_blk_file)
      common /cncvars/ date_str, title,  ininame,
     &        grdname, hisname, rstname
     & , blkfile
      character(len=max_name_size) avgname
      common /cncvars/ avgname
      character(len=max_name_size) clm_file
      common /cncvars/ clm_file
      character(len=max_name_size) bry_file
      common /cncvars/ bry_file
      character(len=42)  vname(3,45+NT
     &     + NT_sed
     &     )
      common /cncvars/ vname
      real(kind=8) bry_time(2), bry_cycle
      integer(kind=4) bry_id, bry_time_id, bry_ncycle, bry_rec, itbry, 
     &                               ntbry
      common /bry_indices/ bry_time, bry_cycle,
     &        bry_id, bry_time_id, bry_ncycle, bry_rec, itbry, ntbry
      integer(kind=4) zeta_west_id
      common /bry_indices/ zeta_west_id
      integer(kind=4) ubar_west_id, vbar_west_id
      common /bry_indices/ ubar_west_id, vbar_west_id
      integer(kind=4) u_west_id, v_west_id
      common /bry_indices/ u_west_id, v_west_id
      integer(kind=4) t_west_id(NT)
      common /bry_indices/ t_west_id
      integer(kind=4) zeta_east_id
      common /bry_indices/ zeta_east_id
      integer(kind=4) ubar_east_id, vbar_east_id
      common /bry_indices/ ubar_east_id, vbar_east_id
      integer(kind=4) u_east_id, v_east_id
      common /bry_indices/ u_east_id, v_east_id
      integer(kind=4) t_east_id(NT)
      common /bry_indices/ t_east_id
      integer(kind=4) zeta_south_id
      common /bry_indices/ zeta_south_id
      integer(kind=4) ubar_south_id, vbar_south_id
      common /bry_indices/ ubar_south_id, vbar_south_id
      integer(kind=4) u_south_id, v_south_id
      common /bry_indices/ u_south_id, v_south_id
      integer(kind=4) t_south_id(NT)
      common /bry_indices/ t_south_id
      integer(kind=4) zeta_north_id
      common /bry_indices/ zeta_north_id
      integer(kind=4) ubar_north_id, vbar_north_id
      common /bry_indices/ ubar_north_id, vbar_north_id
      integer(kind=4) u_north_id, v_north_id
      common /bry_indices/ u_north_id, v_north_id
      integer(kind=4) t_north_id(NT)
      common /bry_indices/ t_north_id
      real(kind=8) zeta_west(0:Mm+1), zeta_west_dt(0:Mm+1,2)
      common /bry_west/ zeta_west, zeta_west_dt
      real(kind=8) ubar_west(0:Mm+1), ubar_west_dt(0:Mm+1,2),
     &     vbar_west(0:Mm+1), vbar_west_dt(0:Mm+1,2)
      common /bry_west/ ubar_west, ubar_west_dt,
     &                  vbar_west, vbar_west_dt
      real(kind=8) u_west(0:Mm+1,N), u_west_dt(0:Mm+1,N,2),
     &     v_west(0:Mm+1,N), v_west_dt(0:Mm+1,N,2)
      common /bry_west/ u_west, u_west_dt,
     &                  v_west, v_west_dt
      real(kind=8) t_west(0:Mm+1,N,NT), t_west_dt(0:Mm+1,N,2,NT)
      common /bry_west/ t_west, t_west_dt
      real(kind=8) zeta_east(0:Mm+1), zeta_east_dt(0:Mm+1,2)
      common /bry_east/ zeta_east, zeta_east_dt
      real(kind=8) ubar_east(0:Mm+1), ubar_east_dt(0:Mm+1,2),
     &     vbar_east(0:Mm+1), vbar_east_dt(0:Mm+1,2)
      common /bry_east/ ubar_east, ubar_east_dt,
     &                  vbar_east, vbar_east_dt
      real(kind=8) u_east(0:Mm+1,N), u_east_dt(0:Mm+1,N,2),
     &     v_east(0:Mm+1,N), v_east_dt(0:Mm+1,N,2)
      common /bry_east/ u_east, u_east_dt,
     &                  v_east, v_east_dt
      real(kind=8) t_east(0:Mm+1,N,NT), t_east_dt(0:Mm+1,N,2,NT)
      common /bry_east/ t_east, t_east_dt
      real(kind=8) zeta_south(0:Lm+1), zeta_south_dt(0:Lm+1,2)
      common /bry_south/ zeta_south, zeta_south_dt
      real(kind=8) ubar_south(0:Lm+1), ubar_south_dt(0:Lm+1,2),
     &     vbar_south(0:Lm+1), vbar_south_dt(0:Lm+1,2)
      common /bry_south/ ubar_south, ubar_south_dt,
     &                   vbar_south, vbar_south_dt
      real(kind=8) u_south(0:Lm+1,N), u_south_dt(0:Lm+1,N,2),
     &     v_south(0:Lm+1,N), v_south_dt(0:Lm+1,N,2)
      common /bry_south/ u_south, u_south_dt,
     &                   v_south, v_south_dt
      real(kind=8) t_south(0:Lm+1,N,NT), t_south_dt(0:Lm+1,N,2,NT)
      common /bry_south/ t_south, t_south_dt
      real(kind=8) zeta_north(0:Lm+1), zeta_north_dt(0:Lm+1,2)
      common /bry_north/ zeta_north, zeta_north_dt
      real(kind=8) ubar_north(0:Lm+1), ubar_north_dt(0:Lm+1,2),
     &     vbar_north(0:Lm+1), vbar_north_dt(0:Lm+1,2)
      common /bry_north/ ubar_north, ubar_north_dt,
     &                   vbar_north, vbar_north_dt
      real(kind=8) u_north(0:Lm+1,N), u_north_dt(0:Lm+1,N,2),
     &     v_north(0:Lm+1,N), v_north_dt(0:Lm+1,N,2)
      common /bry_north/ u_north, u_north_dt,
     &                   v_north, v_north_dt
      real(kind=8) t_north(0:Lm+1,N,NT), t_north_dt(0:Lm+1,N,2,NT)
      common /bry_north/ t_north, t_north_dt
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
      logical new_bgc_flux_his
      integer(kind=4) n_bgc_flux_his, nrpf_bgc_flux_his
     &     , ncid_bgc_flux_his, nrec_bgc_flux_his
      common /scalars_bgc/
     &     new_bgc_flux_his
     &     , n_bgc_flux_his, nrpf_bgc_flux_his
     &     , ncid_bgc_flux_his, nrec_bgc_flux_his
      integer(kind=4) ncid_bgc_flux_avg, nrec_bgc_flux_avg
     &     , nrpf_bgc_flux_avg
      common /scalars_bgc_avg/ ncid_bgc_flux_avg, nrec_bgc_flux_avg
     &     , nrpf_bgc_flux_avg
      character(len=80) bgc_flux_avg_name
      common /c_bgcflux_avg/ bgc_flux_avg_name
      integer(kind=4), parameter :: num_bgcflux_2d = 12
     &     + 3
     &     + NT_sed
      integer(kind=4), parameter :: num_bgcflux = 80
     &     + 7
     &     + NT_sed
      integer(kind=4), dimension(num_bgcflux) :: vid_bec_flux_his
      common /c_bgcflux_bec/ vid_bec_flux_his
      character(len=80) bgc_flux_his_name,
     &     vname_bgcflux(3, num_bgcflux)
      common /c_bgcflux/ bgc_flux_his_name, vname_bgcflux
      integer(kind=4) bgc_flux_hisTime, bgc_flux_hisTstep
     &     , bgc_flux_hisZ
      common /ncids_bgc_flux/ bgc_flux_hisTime, bgc_flux_hisTstep
     &     , bgc_flux_hisZ
      logical new_bgc_flux_avg
      integer(kind=4) nts_bgc_flux_avg, n_bgc_flux_avg
      common /scalars_bgc_avg/
     &     new_bgc_flux_avg,
     &     nts_bgc_flux_avg, n_bgc_flux_avg
      real(kind=8) time_bgc_flux_avg
      common /scalars_bgc_avg_real/ time_bgc_flux_avg
      integer(kind=4), dimension(num_bgcflux) :: vid_bec_flux_avg
      integer(kind=4) :: bgc_flux_avgTstep, bgc_flux_avgTime,
     &    bgc_flux_avgZ
      common /c_bgcflux_avg_bec/ vid_bec_flux_avg, bgc_flux_avgTstep,
     &    bgc_flux_avgTime, bgc_flux_avgZ
      integer(kind=4) NT_PFA
      parameter (NT_PFA = NT)
      logical new_phys_flux_his
      integer(kind=4) n_phys_flux_his
      common /nc_phys_flux_his/n_phys_flux_his, new_phys_flux_his
      integer(kind=4) hisHorXAdvFlux(NT_PFA), hisHorYAdvFlux(NT_PFA)
     &     , hisVertAdvFlux(NT_PFA)
     &     , phys_flux_hisTime, phys_flux_hisTstep, phys_flux_hisZ
     &     , nrpf_phys_flux_his
     &     , ncid_phys_flux_his, nrec_phys_flux_his
      common /inc_phys_flux/ hisHorXAdvFlux, hisHorYAdvFlux
     &     , hisVertAdvFlux
     &     , phys_flux_hisTime, phys_flux_hisTstep, phys_flux_hisZ
     &     , nrpf_phys_flux_his
     &     , ncid_phys_flux_his, nrec_phys_flux_his
      integer(kind=4) indxHorXAdvFlux, indxHorYAdvFlux
      integer(kind=4) indxVertAdvFlux
      parameter(indxHorXAdvFlux = 1)
      parameter(indxHorYAdvFlux = indxHorXAdvFlux + NT_PFA)
      parameter(indxVertAdvFlux = indxHorYAdvFlux + NT_PFA)
      character(len=80) phys_flux_his_name
      character(len=70) vname_phys(4,(3
     &     ) * NT_PFA)
      common /cnc_phys_flux/ phys_flux_his_name, vname_phys
      logical pfa_out(NT)
      common /pfa_by_tracer/ pfa_out
      real(kind=8) time_phys_flux_avg
      common /t_phys_flux_avg/time_phys_flux_avg
      real(kind=8) zeta_phys_flux_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /zeta_phys_flux_avg/zeta_phys_flux_avg
      logical new_phys_flux_avg
      integer(kind=4) nts_phys_flux_avg, n_phys_flux_avg
      common /nc_phys_flux_avg/nts_phys_flux_avg, n_phys_flux_avg
     &        , new_phys_flux_avg
      integer(kind=4) avgHorXAdvFlux(NT_PFA), avgHorYAdvFlux(NT_PFA)
     &     , avgVertAdvFlux(NT_PFA)
     &     , phys_flux_avgTime, phys_flux_avgTstep, phys_flux_avgZ
     &     , nrpf_phys_flux_avg
     &     , ncid_phys_flux_avg, nrec_phys_flux_avg
      common /inc_phys_flux_avg/ avgHorXAdvFlux, avgHorYAdvFlux
     &     , avgVertAdvFlux
     &     , phys_flux_avgTime, phys_flux_avgTstep, phys_flux_avgZ
     &     , nrpf_phys_flux_avg
     &     , ncid_phys_flux_avg, nrec_phys_flux_avg
      character(len=80) phys_flux_avg_name
      common /cnc_phys_flux_avg/ phys_flux_avg_name
      include "mpif.h"
      character keyword*32, fname*64
      character(len=3), parameter :: end_signal='end'
      integer(kind=4), parameter :: kwsize=32, testunit=40, input=15
      integer(kind=4) ierr, iargc, is,ie,  kwlen, lstr, lenstr
     &                                            , itrc
      integer(kind=4), parameter :: ione=1
      ierr=0
      call check_scoord_switches (ierr)
      call check_pre_step_switches (ierr)
      call check_step_uv1_switches (ierr)
      call check_step_uv2_switches (ierr)
      call check_step_t_switches (ierr)
      call check_kpp_switches (ierr)
      if (ierr.ne.0) return
      fname='roms.in'
      if (mynode.eq.0 .and. iargc().eq.1) call getarg(ione,fname)
      call MPI_Bcast(fname,64,MPI_BYTE, 0, ocean_grid_comm, ierr)
      wrthis(indxTime)=.false.
      wrtavg(indxTime)=.false.
      call setup_kwds (ierr)
      open (input, file=fname, status='old', form='formatted', err=97)
   1   keyword='                                '
       read(input,'(A)',err=1,end=99) keyword
       if (ichar(keyword(1:1)).eq.33) goto 1
       is=1
   2   if (is.eq.kwsize) then
         goto 1
       elseif (keyword(is:is).eq.' ') then
         is=is+1
         goto 2
       endif
       ie=is
   3   if (keyword(ie:ie).eq.':') then
         keyword(ie:ie)=' '
         goto 4
       elseif (keyword(ie:ie).ne.' ' .and. ie.lt.kwsize) then
         ie=ie+1
         goto 3
       endif
       goto 1
   4   kwlen=ie-is
       if (is.gt.1) keyword(1:kwlen)=keyword(is:is+kwlen-1)
        if (keyword(1:kwlen).eq.'title') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,'(A)',err=95) title
          lstr=lenstr(title)
          if (mynode.eq.0) write(*,'(/1x,A)') title(1:lstr)
        elseif (keyword(1:kwlen).eq.'time_stepping') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) ntimes, dt, ndtfast, ninfo
          if (mynode.eq.0) write(*,
     &  '(5x,A,I10,3x,A/9x,A,F11.4,2x,A/4x,A,I10,3x,A/6x,A,I10,3x,2A)'
     &    ) 'ntimes =',  ntimes, 'total number of 3D timesteps',
     &          'dt =',     dt,  'time step [sec] for 3D equations',
     &     'ndtfast =', ndtfast, 'mode-splitting ratio',
     &       'ninfo =',   ninfo, 'number of steps between runtime ',
     &                                              'diagnostics'
          dtfast=dt/dble(ndtfast)
        elseif (keyword(1:kwlen) .eq. 'S-coord') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) theta_s, theta_b, hc
          if (mynode.eq.0) write(*,'(2(/4x,A,F10.5,3x,A))')
     &    'theta_s =', theta_s, 'vertical S-coordinate surface',
     &    'theta_b =', theta_b, 'and bottom stretching parameters'
          if (hc.lt.1000._8) then
            if (mynode.eq.0) write(*,'(9x,A,F10.5,3x,2A)')
     &                   'hc =', hc, 'critical depth [m]'
          else
            if (mynode.eq.0) write(*,'(9x,A,ES14.5,3x,2A)')
     &                   'hc =', hc, 'critical depth [m]'
          endif
        elseif (keyword(1:kwlen) .eq. 'rho0') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) rho0
          if (mynode.eq.0) write(*,'(7x,A,F10.4,3x,A)')  'rho0 =',
     &            rho0, 'Boussinesq reference density [kg/m^3].'
        elseif (keyword(1:kwlen) .eq. 'lateral_visc') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) visc2
          if (mynode.eq.0) write(*,'(6x,A,ES10.3,3x,A)') 'visc2 =',
     &      visc2, 'horizontal Laplacian kinematic viscosity [m2/s]'
        elseif (keyword(1:kwlen) .eq. 'gamma2') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) gamma2
          if (mynode.eq.0) write(*,'(5x,A,F4.1,5x,A)') 'gamma2 =',gamma2
     &     , 'slipperiness parameter: free-slip = +1, or no-slip = -1.'
        elseif (keyword(1:kwlen) .eq. 'tracer_diff2') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) (tnu2(itrc),itrc=1,NT)
          do itrc=1,NT
            if (mynode.eq.0) write(*,'(ES10.3,A,I2,2A/32x,A,I2,A)')
     &      tnu2(itrc), '  tnu2(',itrc,')     Horizontal Laplacian ',
     &      'diffusion coefficient (m2/s)', 'for tracer ',  itrc, '.'
          enddo
        elseif (keyword(1:kwlen) .eq. 'bottom_drag') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) rdrg, rdrg2
     &                             , Zob
          if (mynode.eq.0) write(*,'(3(/7x,A,ES10.3,3x,A))')
     &    'rdrg =', rdrg,  'linear bottom drag coefficient [m/s]'
     &  , 'rdrg2=', rdrg2, 'quadratic bottom drag coefficient, nondim'
     &  , ' Zob =', Zob,   'bottom roughness height [m]'
        elseif (keyword(1:kwlen) .eq. 'nudg_cof') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) tauM2_in, tauM2_out, attnM2
     &                       , tauM3_in, tauM3_out
     &                       , tauT_in,  tauT_out
          if (tauM2_in.gt.0._8) then
            tauM2_in =1._8/(tauM2_in *day2sec)
          endif
          if (tauM2_out.gt.0._8) then
            tauM2_out=1._8/(tauM2_out*day2sec)
          endif
          if (tauM3_in.gt.0._8) then
            tauM3_in =1._8/(tauM3_in *day2sec)
          endif
          if (tauM3_out.gt.0._8) then
            tauM3_out=1._8/(tauM3_out*day2sec)
          endif
          if (tauT_in.gt.0._8) then
            tauT_in =1._8/(tauT_in *day2sec)
          endif
          if (tauT_out.gt.0._8) then
            tauT_out=1._8/(tauT_out*day2sec)
          endif
          if (mynode.eq.0) write(*,'(6x,A,ES10.3,ES10.3,2A)')
     &          'tauM2 =', tauM2_in, tauM2_out, '(in/out)  Nudging ',
     &                                   'for barotropic mode [s^-1]'
          if (mynode.eq.0) write(*,'(5x,A,ES10.3,3x,2A)') 'attnM2 =',
     &    attnM2,'open boundary pressure-gradient attenuation [s^-1]'
          if (mynode.eq.0) write(*,'(6x,A,ES10.3,ES10.3,2A)')
     &          'tauM3 =', tauM3_in, tauM3_out, '(in/out)  Nudging ',
     &                                   'for baroclinic mode [s^-1]'
          if (mynode.eq.0) write(*,'(7x,A,ES10.3,ES10.3,2A)')
     &             'tauT =', tauT_in, tauT_out, '(in/out)  Nudging ',
     &                                          'for tracers [s^-1]'
        elseif (keyword(1:kwlen) .eq. 'v_sponge') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) v_sponge
          if (mynode.eq.0) write(*,'(3x,A,F10.2,2x,A)') 'v_sponge =',
     &          v_sponge, 'Maximum viscosity in sponge layer [m2/s]'
        elseif (keyword(1:kwlen) .eq. 'grid') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          open(testunit,file=fname(1:lstr), status='old', err=97)
          close(testunit)
          grdname=fname(1:lstr)
          if (mynode.eq.0) write(*,'(1x,2A)') 'grid file: ',
     &                                     grdname(1:lstr)
        elseif (keyword(1:kwlen) .eq. 'initial') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) nrrec
            read(input,'(A)',err=95) fname
            lstr=lenstr(fname)
            call insert_node (fname, lstr, mynode, NNODES, ierr)
            ininame=fname(1:lstr)
            if (mynode.eq.0) write(*,'(1x,A,I3,2x,3A)')
     &       'initial condition :: rec =', nrrec,  'file = ''',
     &                                   ininame(1:lstr), ''''
        elseif (keyword(1:kwlen) .eq. 'bulk_forcing') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          max_blk=0
          lstr=1
          do while (lstr.gt.0 .and. max_blk.lt.max_blk_file)
            read(input,'(A)',err=95) fname
            lstr=lenstr(fname)
            if (lstr.gt.0) then
              call insert_node (fname, lstr, mynode, NNODES, ierr)
              open (testunit,file=fname(1:lstr),status='old',err=97)
              close(testunit)
              max_blk=max_blk+1
              blkfile(max_blk)=fname(1:lstr)
              if (max_blk.eq.1) then
                 if (mynode.eq.0) write(*,'(1x,A,1x,A)')
     &           'forcing data file(s):', blkfile(max_blk)(1:lstr)
              else
                 if (mynode.eq.0) write(*,'(23x,A)')
     &                                    blkfile(max_blk)(1:lstr)
              endif
            endif
          enddo
        elseif (keyword(1:kwlen) .eq. 'boundary') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          bry_file=fname(1:lstr)
          if (mynode.eq.0) write(*,'(1x,2A)') 'boundary forcing file: ',
     &                                              bry_file(1:lstr)
        elseif (keyword(1:kwlen) .eq. 'restart') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) nrst, nrpfrst
          read(input,'(A)',err=95)  fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          rstname=fname(1:lstr)
          if (mynode.eq.0) write(*,'(1x,A,I8, 2x,A,I5, 2x,3A)')
     &          'restart :: nrst =',nrst, 'rec/file =', nrpfrst,
     &                        'file =''', rstname(1:lstr), ''''
        elseif (keyword(1:kwlen) .eq. 'history') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) ldefhis, nwrt, nrpfhis
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          hisname=fname(1:lstr)
          if (mynode.eq.0) write(*,'(/1x,A,L1,2x,A,I5,2x,A,I4,2x,3A)')
     &       'history :: overwrite = ', ldefhis,      'nwrt =', nwrt,
     &       'rec/file =', nrpfhis, 'file = ''', hisname(1:lstr), ''''
        elseif (keyword(1:kwlen) .eq. 'primary_history_fields') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) wrthis(indxZ),  wrthis(indxUb)
     &                                      ,  wrthis(indxVb)
     &                      ,  wrthis(indxU),  wrthis(indxV)
     &                      , (wrthis(indxT+itrc-1), itrc=1,NT)
     &                      , (wrthis(itrc),
     &                         itrc=indxSedFirst,indxSedFirst+NT_sed-1)
          if ( wrthis(indxZ) .or. wrthis(indxUb) .or. wrthis(indxVb)
     &                         .or. wrthis(indxU) .or. wrthis(indxV)
     &       ) wrthis(indxTime)=.true.
          if (mynode.eq.0) write(*,'(/1x,A,5(/8x,A,T16,L1,T20,A))')
     &                   'fields to write into history file: (T/F)'
     &                 , 'zeta',   wrthis(indxZ),    vname(2,indxZ)
     &                 , 'ubar',   wrthis(indxUb),   vname(2,indxUb)
     &                 , 'vbar',   wrthis(indxVb),   vname(2,indxVb)
     &                 , 'u',      wrthis(indxU),    vname(2,indxU)
     &                 , 'v',      wrthis(indxV),    vname(2,indxV)
          do itrc=1,NT
            if (wrthis(indxT+itrc-1)) wrthis(indxTime)=.true.
            if (mynode.eq.0) write(*,'(8x,A,I2,A,T16,L1,T20,A)') 't(',
     &        itrc, ')', wrthis(indxT+itrc-1), vname(2,indxT+itrc-1)
          enddo
          do itrc = 1, NT_sed
            if (wrthis(indxSedFirst+itrc-1)) wrthis(indxTime)=.true.
            if (mynode.eq.0) write(*, '(6x,L1,2x,A,I2,A,I2,A)')
     &           wrthis(indxSedFirst+itrc-1), 'write T_sed(',
     &           itrc,')  Sediment tracer of index ', itrc,'.'
          end do
        elseif (keyword(1:kwlen) .eq. 'auxiliary_history_fields') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) wrthis(indxR), wrthis(indxO)
     &          ,  wrthis(indxW),  wrthis(indxAkv),  wrthis(indxAkt)
     &                                            ,  wrthis(indxAks)
     &                                            ,  wrthis(indxHbls)
     &                                           ,  wrthis(indxHbbl)
          if ( wrthis(indxR) .or. wrthis(indxO) .or. wrthis(indxW)
     &                     .or. wrthis(indxAkv) .or. wrthis(indxAkt)
     &                                          .or. wrthis(indxAks)
     &                                          .or. wrthis(indxHbls)
     &                                          .or. wrthis(indxHbbl)
     &       ) wrthis(indxTime)=.true.
          if (mynode.eq.0) write(*,'(8(/8x,A,T16,L1,T20,A))')
     &                   'rho',    wrthis(indxR),    vname(2,indxR)
     &                 , 'Omega',  wrthis(indxO),    vname(2,indxO)
     &                 , 'W',      wrthis(indxW),    vname(2,indxW)
     &                 , 'Akv',    wrthis(indxAkv),  vname(2,indxAkv)
     &                 , 'Akt',    wrthis(indxAkt),  vname(2,indxAkt)
     &                 , 'Aks',    wrthis(indxAks),  vname(2,indxAks)
     &                 , 'hbls',   wrthis(indxHbls), vname(2,indxHbls)
     &                 , 'hbbl',   wrthis(indxHbbl), vname(2,indxHbbl)
        elseif (keyword(1:kwlen) .eq. 'averages') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) ntsavg, navg, nrpfavg
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          avgname=fname(1:lstr)
          if (mynode.eq.0) write(*,'(/1x,A,I5,2x,A,I5,2x,A,I4,2x,3A)')
     &       'averages :: ntsavg = ', ntsavg,      'navg =', navg,
     &       'rec/file =', nrpfavg, 'file = ''', avgname(1:lstr), ''''
        elseif (keyword(1:kwlen) .eq. 'primary_averages') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) wrtavg(indxZ),  wrtavg(indxUb)
     &                                        ,  wrtavg(indxVb)
     &                      ,  wrtavg(indxU),     wrtavg(indxV)
     &                      , (wrtavg(indxT+itrc-1), itrc=1,NT)
     &                      , (wrtavg(itrc),
     &                         itrc=indxSedFirst,indxSedFirst+NT_sed-1)
          if ( wrtavg(indxZ) .or. wrtavg(indxUb) .or. wrtavg(indxVb)
     &                       .or. wrtavg(indxU)  .or. wrtavg(indxV)
     &       ) wrtavg(indxTime)=.true.
          if (mynode.eq.0) write(*,'(/1x,A,5(/8x,A,T16,L1,T20,A))')
     &                   'fields to compute time averages of: (T/F)'
     &                 , 'zeta',   wrtavg(indxZ),    vname(2,indxZ)
     &                 , 'ubar',   wrtavg(indxUb),   vname(2,indxUb)
     &                 , 'vbar',   wrtavg(indxVb),   vname(2,indxVb)
     &                 , 'u',      wrtavg(indxU),    vname(2,indxU)
     &                 , 'v',      wrtavg(indxV),    vname(2,indxV)
          do itrc=1,NT
            if (wrtavg(indxT+itrc-1)) wrtavg(indxTime)=.true.
            if (mynode.eq.0) write(*,'(8x,A,I2,A,T16,L1,T20,A)') 't(',
     &        itrc, ')', wrtavg(indxT+itrc-1), vname(2,indxT+itrc-1)
          enddo
        elseif (keyword(1:kwlen) .eq. 'auxiliary_averages') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) wrtavg(indxR), wrtavg(indxO)
     &          ,  wrtavg(indxW),  wrtavg(indxAkv),  wrtavg(indxAkt)
     &                                            ,  wrtavg(indxAks)
     &                                            ,  wrtavg(indxHbls)
     &                                            ,  wrtavg(indxHbbl)
          if ( wrtavg(indxR) .or. wrtavg(indxO) .or. wrtavg(indxW)
     &                     .or. wrtavg(indxAkv) .or. wrtavg(indxAkt)
     &                                          .or. wrtavg(indxAks)
     &                                          .or. wrtavg(indxHbls)
     &                                          .or. wrtavg(indxHbbl)
     &       ) wrtavg(indxTime)=.true.
          if (mynode.eq.0) write(*,'(8(/8x,A,T16,L1,T20,A))')
     &                   'rho',    wrtavg(indxR),    vname(2,indxR)
     &                 , 'Omega',  wrtavg(indxO),    vname(2,indxO)
     &                 , 'W',      wrtavg(indxW),    vname(2,indxW)
     &                 , 'Akv',    wrtavg(indxAkv),  vname(2,indxAkv)
     &                 , 'Akt',    wrtavg(indxAkt),  vname(2,indxAkt)
     &                 , 'Aks',    wrtavg(indxAks),  vname(2,indxAks)
     &                 , 'hbls',   wrtavg(indxHbls), vname(2,indxHbls)
     &                 , 'hbbl',   wrtavg(indxHbbl), vname(2,indxHbbl)
        elseif (keyword(1:kwlen).eq.'bgc_flux_histories') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) new_bgc_flux_his,
     &         n_bgc_flux_his, nrpf_bgc_flux_his
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          bgc_flux_his_name=fname(1:lstr)
          if (mynode.eq.0) write(*,
     &           '(I10,2x,A,1x,A/32x,A/,6x,A,2x,A,1x,A,I3)')
     &      n_bgc_flux_his,
     &         'n_bgc_flux_his     Number of timesteps between',
     &     'writing of',
     &         ' biogeochemical fluxes into file.',
     &     'Biogeochemical flux file:',
     &         bgc_flux_his_name(1:lstr),
     &     'rec/file =', nrpf_bgc_flux_his
        elseif (keyword(1:kwlen).eq.'bgc_flux_averages') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) new_bgc_flux_avg,
     &         nts_bgc_flux_avg, n_bgc_flux_avg, nrpf_bgc_flux_avg
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          bgc_flux_avg_name=fname(1:lstr)
          if (mynode.eq.0) write(*,
     &           '(2(I10,2x,A,1x,A/32x,A/),6x,A,2x,A,1x,A,I3)')
     &      nts_bgc_flux_avg,
     &         'nts_bgc_flux_avg   Starting timestep for the',
     &           ' accumulation of output',
     &         ' time-averaged biogeochemical fluxes.',
     &      n_bgc_flux_avg,
     &         'n_bgc_flux_avg     Number of timesteps between',
     &     'writing of',
     &         ' time-averaged biogeochemical fluxes into file.',
     &     'Averaged biogeochemical flux file:',
     &         bgc_flux_avg_name(1:lstr),
     &     'rec/file =', nrpf_bgc_flux_avg
        elseif (keyword(1:kwlen).eq.'phys_flux_histories') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) new_phys_flux_his,
     &         n_phys_flux_his, nrpf_phys_flux_his
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          phys_flux_his_name=fname(1:lstr)
          if (mynode.eq.0) write(*,
     &           '(I10,2x,A,1x,A/32x,A/,6x,A,2x,A,1x,A,I3)')
     &      n_phys_flux_his,
     &         'n_phys_flux_his    Number of timesteps between',
     &     ' writing of',
     &         ' physical fluxes into file.',
     &     'Physical flux file:',
     &         phys_flux_his_name(1:lstr),
     &     'rec/file =', nrpf_phys_flux_his
        elseif (keyword(1:kwlen).eq.'phys_flux_averages') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) new_phys_flux_avg,
     &         nts_phys_flux_avg, n_phys_flux_avg,
     &         nrpf_phys_flux_avg
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          phys_flux_avg_name=fname(1:lstr)
          read(input,*,err=95) (wrt_pfa(itrc), itrc=1,NT)
          if (mynode.eq.0) write(*,
     &           '(2(I10,2x,A,1x,A/32x,A/),6x,A,2x,A,1x,A,I3)')
     &      nts_phys_flux_avg,
     &         'nts_phys_flux_avg  Starting timestep for the',
     &           ' accumulation of output',
     &         ' time-averaged physical fluxes.',
     &      n_phys_flux_avg,
     &         'n_phys_flux_avg    Number of timesteps between',
     &     ' writing of',
     &         ' time-averaged physical fluxes into file.',
     &     'Averaged physical flux file:',
     &         phys_flux_avg_name(1:lstr),
     &     'rec/file =', nrpf_phys_flux_avg
          do itrc=1,NT
            if (mynode.eq.0) write(*, '(6x,L1,2x,A,I2,A,I2,A)')
     &                       wrt_pfa(itrc),
     &            ' write PFA_avg for tracer of index ', itrc,'.'
          enddo
        elseif (keyword(1:kwlen).eq.'bulk_diags_histories') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) new_bulk_diags_his,
     &         n_bulk_diags_his, nrpf_bulk_diags_his
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          bulk_diags_his_name=fname(1:lstr)
          if (mynode.eq.0) write(*,
     &           '(I10,2x,A,1x,A/32x,A/,6x,A,2x,A,1x,A,I3)')
     &      n_bulk_diags_his,
     &         'n_bulk_diags_his    Number of timesteps between',
     &     ' writing of',
     &         ' bulk diags into file.',
     &     'Bulk Diags file:',
     &         bulk_diags_his_name(1:lstr),
     &     'rec/file =', nrpf_bulk_diags_his
          if (mynode.eq.0) write(*,*) 'new_bulk_diags_his',
     &                                new_bulk_diags_his
        elseif (keyword(1:kwlen).eq.'bulk_diags_averages') then
          call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) new_bulk_diags_avg,
     &         nts_bulk_diags_avg, n_bulk_diags_avg,
     &         nrpf_bulk_diags_avg
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          call insert_node (fname, lstr, mynode, NNODES, ierr)
          bulk_diags_avg_name=fname(1:lstr)
          if (mynode.eq.0) write(*,
     &           '(2(I10,2x,A,1x,A/32x,A/),6x,A,2x,A,1x,A,I3)')
     &      nts_bulk_diags_avg,
     &         'nts_bulk_diags_avg  Starting timestep for the',
     &           ' accumulation of output',
     &         ' time-averaged physical fluxes.',
     &      n_bulk_diags_avg,
     &         'n_bulk_diags_avg    Number of timesteps between',
     &     ' writing of',
     &         ' time-averaged physical fluxes into file.',
     &     'Averaged physical flux file:',
     &         bulk_diags_avg_name(1:lstr),
     &     'rec/file =', nrpf_bulk_diags_avg
          if (mynode.eq.0) write(*,*) 'new_bulk_diags_his',
     &                                new_bulk_diags_his
      elseif (keyword(1:kwlen).eq.'pCO2_atm_file') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) pco2_atm_file
        if (mynode.eq.0) write(*,'(A,A)')
     &        'Read atmospheric pCO2 from ', pco2_atm_file
        else
          if (mynode.eq.0) write(*,'(/1x,4A/)') 'WARNING: ',
     &                'Unrecognized keyword ''', keyword(1:kwlen),
     &                                      ''' --> DISREGARDED.'
        endif
       if (keyword(1:kwlen) .eq. end_signal) goto 99
      goto 1
  95  write(*,'(/1x,4A/)') '### ERROR: read_inp :: Cannot read ',
     &                       'entry ''', keyword(1:kwlen), '''.'
      ierr=ierr+1
      goto 99
  97  lstr=lenstr(fname)
      write(*,'(/1x,4A/)') '### ERROR: read_inp :: Cannot find ',
     &                     'input file ''', fname(1:lstr), '''.'
      ierr=ierr+1
  99  close (input)
      if (ierr.eq.0) call check_kwds (ierr)
      if (ierr.ne.0) then
        write(*,'(/1x,2A,I3,1x,A/)') '### ERROR: read_inp :: ',
     &  'A total of', ierr, 'configuration errors discovered.'
       return
      endif
      call MPI_Barrier (ocean_grid_comm, ierr)
      return
      end
      subroutine cancel_kwd (keyword, ierr)
      implicit none
      integer(kind=4), parameter :: max_opt_size=2048
      character(len=max_opt_size) cpps, srcs, kwds
      common /strings/ cpps, srcs, kwds
      character(len=*) keyword
      integer(kind=4) ierr, is,i,ie, lenkw,lenstr
      lenkw=lenstr(keyword)
      is=1
      do while (is.gt.0 .and. is.lt.max_opt_size)
        do while (kwds(is:is).eq.' ' .and. is.lt.max_opt_size)
          is=is+1
        enddo
        ie=is+1
        do while (kwds(ie:ie).ne.' ' .and. ie.lt.max_opt_size)
          ie=ie+1
        enddo
        if (lenkw.eq.ie-is .and. kwds(is:ie-1).eq.keyword) then
          do i=is,ie-1
            kwds(i:i)=' '
          enddo
          is=0
        else
          is=ie
        endif
      enddo
      if (is.ne.0) then
        write(*,'(/A)') '##### ERROR #####'
        write(*,'(2(1x,A,1x,A,1x,A/)/)') 'cancel_kwd:',
     &         'Can not cancel keyword:',  keyword(1:lenkw),
     &         'check input script for possible',
     &         'duplicated keywords.'
        write(*,'(A/)') '#################'
        ierr=ierr+1
      endif
      return
      end
      subroutine check_kwds (ierr)
      implicit none
      integer(kind=4), parameter :: max_opt_size=2048
      character(len=max_opt_size) cpps, srcs, kwds
      common /strings/ cpps, srcs, kwds
      integer(kind=4) ierr, is,ie
      is=1
      do while (is.lt.max_opt_size)
        do while (kwds(is:is).eq.' ' .and. is.lt.max_opt_size)
          is=is+1
        enddo
        if (is.lt.max_opt_size) then
          ie=is+1
          do while (kwds(ie:ie).ne.' ' .and. ie.lt.max_opt_size)
            ie=ie+1
          enddo
          ierr=ierr+1
          write(*,'(/1x,A,1x,A/)') '### ERROR: keyword not found:',
     &                                               kwds(is:ie-1)
          is=ie
        endif
      enddo
      return
      end
