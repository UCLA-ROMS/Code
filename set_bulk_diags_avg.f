      subroutine set_bulk_diags_avg(tile)
      implicit none
      integer(kind=4) tile
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
      integer(kind=4) istr,iend, jstr,jend, i_X,j_E
      integer(kind=4) inode,jnode,  p_W, p_SW, p_S, p_SE, p_E, p_NE, 
     &                                p_N,
     &                                      p_NW,  exc_call_count
      logical west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      common /hidden_mpi_vars/ inode,jnode, p_W, p_SW, p_S, p_SE,
     &                      p_E, p_NE, p_N, p_NW,  exc_call_count,
     &        west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      save /hidden_mpi_vars/
      integer(kind=4) size_X, margin_X, size_E, margin_E
        j_E=tile/NSUB_X
        i_X=tile-j_E*NSUB_X
        if (mod(j_E,2).eq.1) i_X=NSUB_X-1 -i_X
        if (mod(inode,2).gt.0) then
          i_X=NSUB_X-1 -i_X
        endif
        if (mod(jnode,2).gt.0) then
          j_E=NSUB_E-1 -j_E
        endif
        size_X=(ieast-iwest+NSUB_X)/NSUB_X
        margin_X=(NSUB_X*size_X - ieast+iwest-1)/2
        istr=iwest-margin_X + i_X*size_X
        iend=min( istr + size_X-1 ,ieast)
        istr=max(istr,iwest)
        size_E=(jnorth-jsouth +NSUB_E)/NSUB_E
        margin_E=(NSUB_E*size_E -jnorth+jsouth-1)/2
        jstr=jsouth-margin_E + j_E*size_E
        jend=min( jstr + size_E-1 ,jnorth)
        jstr=max(jstr,jsouth)
      call set_bulk_diags_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end
      subroutine set_bulk_diags_avg_tile(Istr,Iend,Jstr,Jend)
      implicit none
      integer(kind=4) Istr,Iend,Jstr,Jend, i, j, k, iflux, iflux_out
      real(kind=8) cff, cff1
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
      real(kind=8) zeta(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real(kind=8) ubar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real(kind=8) vbar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      common /ocean_zeta/zeta /ocean_ubar/ubar /ocean_vbar/vbar
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
      real(kind=8) zeta_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) ubar_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) vbar_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_zeta/zeta_avg /avg_ubar/ubar_avg
     &                          /avg_vbar/vbar_avg
      real(kind=8) u_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) v_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) t_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real(kind=8) rho_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real(kind=8) w_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /avg_u/u_avg /avg_v/v_avg /avg_t/t_avg
     &                /avg_rho/rho_avg /avg_w/w_avg
      real(kind=8) wvl_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /avg_wvl/wvl_avg
      real(kind=8) akv_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real(kind=8) akt_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /avg_akv/akv_avg /avg_akt/akt_avg
      real(kind=8) aks_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /avg_aks/aks_avg
      real(kind=8) hbl_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_hbl/hbl_avg
      real(kind=8) hbbl_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_hbbl/hbbl_avg
      real(kind=8) sustr_blk_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) svstr_blk_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) shflx_net_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) shflx_lat_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) shflx_sen_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) shflx_rad_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) swflx_emp_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) shflx_wwk_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) surf_u_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) surf_v_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /avg_sustr_blk/sustr_blk_avg
      common /avg_svstr_blk/svstr_blk_avg
      common /avg_shflx_net/shflx_net_avg
      common /avg_shflx_lat/shflx_lat_avg
      common /avg_shflx_sen/shflx_sen_avg
      common /avg_shflx_rad/shflx_rad_avg
      common /avg_swflx_emp/swflx_emp_avg
      common /avg_shflx_wwk/shflx_wwk_avg
      common /avg_surf_u/surf_u_avg
      common /avg_surf_v/surf_v_avg
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
      integer(kind=4) istrU, istrR, iendR
      integer(kind=4) jstrV, jstrR, jendR
      if (istr.eq.iwest .and. .not.west_exchng) then
        istrR=istr-1
        istrU=istr+1
      else
        istrR=istr
        istrU=istr
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        iendR=iend+1
      else
        iendR=iend
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        jstrR=jstr-1
        jstrV=jstr+1
      else
        jstrR=jstr
        jstrV=jstr
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        jendR=jend+1
      else
        jendR=jend
      endif
      if (iic.gt.nts_bulk_diags_avg) then
        cff=1._8/dble(n_bulk_diags_avg)
         if (n_bulk_diags_avg .eq. 1 .or.
     &        mod(iic-nts_bulk_diags_avg,n_bulk_diags_avg).eq.1) then
           if ((istr.eq.iwest .and. jstr.eq.jsouth)) then
            if (mynode.eq.0) write(*,'(7x,A,I8,2x,A,I8,2x,A,I5)')
     &               'set_bdiags_avg :: started averaging at iic=', iic,
     &                          'ntsavg=', ntsavg, 'navg=', navg
           endif
           cff1=0._8
        else
          cff1=1._8
        endif
        if ((istr.eq.iwest .and. jstr.eq.jsouth)) 
     &      time_bulk_diags_avg=cff1*time_bulk_diags_avg + cff*time
            do j = jstrR, jendR
               do i = istrR, iendR
                  zeta_bulk_diags_avg(i,j) =
     &                         cff1*zeta_bulk_diags_avg(i,j) +
     &                         cff*zeta(i,j,knew)
                  sustr_blk_avg(i,j)=cff1*sustr_blk_avg(i,j) +
     &                         cff*sustr_blk(i,j)
                  svstr_blk_avg(i,j)=cff1*svstr_blk_avg(i,j) +
     &                         cff*svstr_blk(i,j)
                  shflx_net_avg(i,j)=cff1*shflx_net_avg(i,j) +
     &                         cff*shflx_net(i,j)
                  shflx_lat_avg(i,j)=cff1*shflx_lat_avg(i,j) +
     &                         cff*shflx_lat(i,j)
                  shflx_sen_avg(i,j)=cff1*shflx_sen_avg(i,j) +
     &                         cff*shflx_sen(i,j)
                  shflx_rad_avg(i,j)=cff1*shflx_rad_avg(i,j) +
     &                         cff*shflx_rad(i,j)
                  swflx_emp_avg(i,j)=cff1*swflx_emp_avg(i,j) +
     &                         cff*swflx_emp(i,j)
                  shflx_wwk_avg(i,j)=cff1*shflx_wwk_avg(i,j) +
     &                         cff*shflx_wwk(i,j)
                  surf_u_avg(i,j)=cff1*surf_u_avg(i,j) +
     &                         cff*surf_u(i,j)
                  surf_v_avg(i,j)=cff1*surf_v_avg(i,j) +
     &                         cff*surf_v(i,j)
            enddo
          enddo
      endif
      return
      end
