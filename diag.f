      subroutine diag (tile)
      implicit none
      integer(kind=4) tile,   i,j
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
      integer(kind=4), parameter :: size_XI=7+(Lm+NSUB_X-1)/NSUB_X,
     &                      size_ETA=7+(Mm+NSUB_E-1)/NSUB_E,
     &         sse=size_ETA/(N+1),  ssz=(N+1)/size_ETA,
     &         N2d=size_XI*(sse*size_ETA+ssz*(N+1))/(sse+ssz),
     &         N3d=size_XI*size_ETA*(N+1)
      real(kind=8) A2d(N2d,32)
      real(kind=8) A3d(N3d,6)
      integer(kind=4) iA2d(N2d)
      common /prv_scrch/ A3d, A2d, iA2d
C$OMP THREADPRIVATE(/prv_scrch/)
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
      j=min(iic,2*ninfo)-1
      i=1
      do while(i.lt.j)
        i=2*i
      enddo
      if (mod(iic-1,min(i,ninfo)) .eq. 0) then
        call diag_tile (istr,iend,jstr,jend, A2d(1,1),A2d(1,2),A2d(1,3)
     &               , A2d(1,4), A2d(1,5), A2d(1,6), A2d(1,7), A2d(1,8)
     &                                                                )
      endif
      return
      end
      subroutine diag_tile (istr,iend,jstr,jend, dVol, ke,pe
     &                            , ke2b, ke3bc, kesrf, ub,vb
     &                                                      )
      implicit none
      integer(kind=4) istr,iend,jstr,jend
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
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: dVol, ke, 
     &                                 pe
     &                             , ke2b, ke3bc, kesrf, ub, vb
      include "mpif.h"
      integer(kind=4) size, step, itag, status(MPI_STATUS_SIZE)
      real(kind=8) buff(16)
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
      real(kind=8) rho1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /eos_rho1/ rho1
      real(kind=8) qp1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /eos_qp1/ qp1
      real(kind=8), parameter :: qp2=0.0000172_8
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
      real(kind=8) area, volume, bc_crss
      common /comm_vars/ area, volume, bc_crss
      real(kind=8) global_sum(0:2*NT+1), global_srf_sum(0:NT)
      common /comm_vars_bio/ global_sum, global_srf_sum
      real(kind=8) hmin,hmax, grdmin,grdmax, rx0,rx1, Cg_min,Cg_max, 
     &                               Cu_Cor
      common /comm_vars/ hmin,hmax, grdmin,grdmax, rx0,rx1,
     &                                        Cg_min,Cg_max, Cu_Cor
      real(kind=4) cpu_all(4)
      integer(kind=4) trd_count, tile_count, bc_count, mcheck, 
     &                             first_time
      common /comm_vars/ cpu_all, trd_count,
     &                   tile_count, bc_count, mcheck, first_time
      real(kind=8)  avzeta, avke,   prev_ke, avpe,
     &                   avke2b, avke3bc, avkesrf
      common /diag_vars/ avzeta, avke,    prev_ke,
     &             avpe, avke2b, avke3bc, avkesrf
      real(kind=8) v2d_max
      common /diag_vars/ v2d_max
      real(kind=8) Cu_Adv,  Cu_W
      integer(kind=4) i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
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
      real(kind=8) cff, dA, v2, my_v2d_max
      real(kind=8) my_avzeta, my_ke,   my_pe
     &        , my_ke3bc,  my_ke2b, my_kesrf
      integer(kind=4) i,j,k, nsubs, ierr, is,ie
      real(kind=8) v2bc
      real(kind=8) my_Cu_Adv, my_Cu_W, ciV, cx, cw
      integer(kind=4) my_i_cmax, my_j_cmax, my_k_cmax
      integer(kind=4) isize,itg, js,jsize,jtg
      integer(kind=4), parameter :: max_check_line=128
      character(len=max_check_line) check_line
      character(len=18) tstring
      ierr=0
      my_v2d_max=0._8
      my_Cu_Adv=0._8  ;  my_Cu_W=0._8
      my_i_cmax=0   ;  my_j_cmax=0  ;  my_k_cmax=0
      do j=jstr,jend+1
        do i=istr,iend+1
          ub(i,j)=(Hz(i,j,N)+Hz(i-1,j,N))*u(i,j,N,nstp)
          vb(i,j)=(Hz(i,j,N)+Hz(i,j-1,N))*v(i,j,N,nstp)
        enddo
        do k=N-1,2,-1
          do i=istr,iend+1
            ub(i,j)=ub(i,j)+(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nstp)
            vb(i,j)=vb(i,j)+(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nstp)
          enddo
        enddo
        do i=istr,iend+1
          ub(i,j)=(ub(i,j)+(Hz(i,j,1)+Hz(i-1,j,1))*u(i,j,1,nstp))
     &         /(z_w(i,j,N)+z_w(i-1,j,N)-z_w(i,j,0)-z_w(i-1,j,0))
          vb(i,j)=(vb(i,j)+(Hz(i,j,1)+Hz(i,j-1,1))*v(i,j,1,nstp))
     &         /(z_w(i,j,N)+z_w(i,j-1,N)-z_w(i,j,0)-z_w(i,j-1,0))
        enddo
      enddo
      cff=g/rho0
      do j=jstr,jend
        do i=istr,iend
          v2=0.5_8*(ub(i,j)**2+ub(i+1,j)**2 +vb(i,j)**2+vb(i,j+1)**2)
          my_v2d_max=max(my_v2d_max, v2)
          ke(i,j)=0._8
          pe(i,j)=0.5_8*g*z_w(i,j,N)*z_w(i,j,N)
          ke2b(i,j)=0.5_8*(z_w(i,j,N)-z_w(i,j,0))*v2
          ke3bc(i,j)=0._8
          kesrf(i,j)=0.25_8*( u(i,j,N,nstp)**2 + u(i+1,j,N,nstp)**2
     &                     +v(i,j,N,nstp)**2 + v(i,j+1,N,nstp)**2)
        enddo
        do k=N,1,-1
          do i=istr,iend
            v2=0.5_8*( u(i,j,k,nstp)**2 + u(i+1,j,k,nstp)**2
     &              +v(i,j,k,nstp)**2 + v(i,j+1,k,nstp)**2)
            v2bc=0.5_8*( (u(i  ,j,k,nstp)-ub(i  ,j))**2
     &                +(u(i+1,j,k,nstp)-ub(i+1,j))**2
     &                +(v(i,j  ,k,nstp)-vb(i,j  ))**2
     &                +(v(i,j+1,k,nstp)-vb(i,j+1))**2)
            ciV=dt*rmask(i,j)*pm(i,j)*pn(i,j)/Hz(i,j,k)
            cw=ciV*( max(We(i,j,k)+Wi(i,j,k), 0._8)
     &                      -min(We(i,j,k-1)+Wi(i,j,k-1), 0._8) )
            cx=cw+ciV*( max(FlxU(i+1,j,k), 0._8) -min(FlxU(i,j,k), 0._8)
     &                 +max(FlxV(i,j+1,k), 0._8) -min(FlxV(i,j,k), 
     &                               0._8))
            if (cx .gt. my_Cu_Adv) then
              my_Cu_Adv=cx ; my_Cu_W=cw
              my_i_cmax=i  ; my_j_cmax=j ; my_k_cmax=k
            endif
            ke(i,j)=ke(i,j) + 0.5_8*v2*Hz(i,j,k)
            pe(i,j)=pe(i,j) + cff*Hz(i,j,k)
     &      *(rho1(i,j,k)+qp1(i,j,k)*(z_w(i,j,N)-z_r(i,j,k)))
     &                               *(z_r(i,j,k)-z_w(i,j,0))
            ke3bc(i,j)=ke3bc(i,j) + 0.5_8*v2bc*Hz(i,j,k)
          enddo
        enddo
        do i=istr,iend
          dA=rmask(i,j)/(pm(i,j)*pn(i,j))
          dVol(i,j) = dA * z_w(i,j,N)
          ke(i,j)   = dA * ke(i,j)
          pe(i,j)   = dA * pe(i,j)
          ke2b(i,j) = dA * ke2b(i,j)
          ke3bc(i,j)= dA * ke3bc(i,j)
          kesrf(i,j)= dA * kesrf(i,j)
        enddo
      enddo
      isize=iend-istr
      jsize=jend-jstr
      do while (isize>0 .or. jsize>0)
        if (jsize>0) then
          js=(jsize+1)/2-1
          do j=0,js
            jtg=jstr+j
            do i=istr,istr+isize
               dVol(i,jtg) =  dVol(i,jtg+j) +  dVol(i,jtg+j+1)
                 ke(i,jtg) =    ke(i,jtg+j) +    ke(i,jtg+j+1)
                 pe(i,jtg) =    pe(i,jtg+j) +    pe(i,jtg+j+1)
               ke2b(i,jtg) =  ke2b(i,jtg+j) +  ke2b(i,jtg+j+1)
              ke3bc(i,jtg) = ke3bc(i,jtg+j) + ke3bc(i,jtg+j+1)
              kesrf(i,jtg) = kesrf(i,jtg+j) + kesrf(i,jtg+j+1)
            enddo
          enddo
          if (2*js+1 < jsize) then
            js=js+1
            jtg=jstr+js
            do i=istr,istr+isize
               dVol(i,jtg) =  dVol(i,jtg+js)
                 ke(i,jtg) =    ke(i,jtg+js)
                 pe(i,jtg) =    pe(i,jtg+js)
               ke2b(i,jtg) =  ke2b(i,jtg+js)
              ke3bc(i,jtg) = ke3bc(i,jtg+js)
              kesrf(i,jtg) = kesrf(i,jtg+js)
            enddo
          endif
          jsize=js
        endif
        if (isize>0) then
          is=(isize+1)/2-1
          do j=jstr,jstr+jsize
            do i=0,is
              itg=istr+i
               dVol(itg,j) =  dVol(itg+i,j) +  dVol(itg+i+1,j)
                 ke(itg,j) =    ke(itg+i,j) +    ke(itg+i+1,j)
                 pe(itg,j) =    pe(itg+i,j) +    pe(itg+i+1,j)
               ke2b(itg,j) =  ke2b(itg+i,j) +  ke2b(itg+i+1,j)
              ke3bc(itg,j) = ke3bc(itg+i,j) + ke3bc(itg+i+1,j)
              kesrf(itg,j) = kesrf(itg+i,j) + kesrf(itg+i+1,j)
            enddo
          enddo
          if (2*is+1 < isize) then
            is=is+1
            itg=istr+is
            do j=jstr,jstr+jsize
               dVol(itg,j) =  dVol(itg+is,j)
                 ke(itg,j) =    ke(itg+is,j)
                 pe(itg,j) =    pe(itg+is,j)
               ke2b(itg,j) =  ke2b(itg+is,j)
              ke3bc(itg,j) = ke3bc(itg+is,j)
              kesrf(itg,j) = kesrf(itg+is,j)
            enddo
          endif
          isize=is
        endif
      enddo
      my_avzeta=dVol(istr,jstr)
      my_ke=ke(istr,jstr)
      my_pe=pe(istr,jstr)
      my_ke2b=ke2b(istr,jstr)
      my_ke3bc=ke3bc(istr,jstr)
      my_kesrf=kesrf(istr,jstr)
      if ((iend-istr.eq.ieast-iwest .and.  jend-jstr.eq.jnorth-jsouth)) 
     &                                then
        nsubs=1
      else
        nsubs=NSUB_X*NSUB_E
      endif
C$OMP CRITICAL (diag_cr_rgn)
      if (tile_count.eq.0) then
        avzeta=my_avzeta
        avke=my_ke
        avpe=my_pe
        v2d_max=my_v2d_max
        avke2b=my_ke2b
        avke3bc=my_ke3bc
        avkesrf=my_kesrf
        Cu_Adv=my_Cu_Adv
        Cu_W=my_Cu_W
        i_cx_max=my_i_cmax
        j_cx_max=my_j_cmax
        k_cx_max=my_k_cmax
      else
        avzeta=avzeta+my_avzeta
        avke=avke+my_ke
        avpe=avpe+my_pe
        v2d_max=max(v2d_max, my_v2d_max)
        avke2b=avke2b+my_ke2b
        avke3bc=avke3bc+my_ke3bc
        avkesrf=avkesrf+my_kesrf
        if (my_Cu_Adv.gt.Cu_Adv) then
          Cu_Adv=my_Cu_Adv
          Cu_W=my_Cu_W
          i_cx_max=my_i_cmax
          j_cx_max=my_j_cmax
          k_cx_max=my_k_cmax
        endif
      endif
      tile_count=tile_count+1
      if (tile_count.eq.nsubs) then
        tile_count=0
        i_cx_max=i_cx_max + iSW_corn
        j_cx_max=j_cx_max + jSW_corn
        if (NNODES.gt.1) then
          size=NNODES
          do while (size.gt.1)
            step=(size+1)/2
            if (mynode.ge.step .and. mynode.lt.size) then
              buff(1)=may_day_flag
              buff(2)=avzeta
              buff(3)=avke
              buff(4)=avpe
              buff(5)=avke2b
              buff(6)=avke3bc
              buff(7)=avkesrf
              buff(8)=v2d_max
              buff(9)=Cu_Adv
              buff(10)=Cu_W
              buff(11)=i_cx_max
              buff(12)=j_cx_max
              buff(13)=k_cx_max
              itag=mynode+400
              call MPI_Send (buff, 14, MPI_REAL8, mynode-step,
     &                       itag, ocean_grid_comm,          ierr)
            elseif (mynode .lt. size-step) then
              itag=mynode+400+step
              call MPI_Recv (buff, 14, MPI_REAL8, mynode+step,
     &                       itag, ocean_grid_comm,  status, ierr)
              i=buff(1)+0.0001_8
              if (i.gt.0 .and. may_day_flag.eq.0) may_day_flag=i
              avzeta=avzeta+buff(2)
              avke=avke+buff(3)
              avpe=avpe+buff(4)
              avke2b=avke2b+buff(5)
              avke3bc=avke3bc+buff(6)
              avkesrf=avkesrf+buff(7)
              v2=buff(8)
              v2d_max=max(v2d_max, v2)
              if (buff(9).gt.Cu_Adv) then
                Cu_Adv=buff(9)
                Cu_W=buff(10)
                i_cx_max=buff(11) +0.0001_8
                j_cx_max=buff(12) +0.0001_8
                k_cx_max=buff(13) +0.0001_8
              endif
            endif
            size=step
          enddo
        endif
        if (mynode.eq.0) then
          avke=avke/(volume+avzeta)
          avpe=avpe/(volume+avzeta)
          avke2b=avke2b/(volume+avzeta)
          avke3bc=avke3bc/(volume+avzeta)
          avkesrf=avkesrf/area
          v2d_max=sqrt(v2d_max)
          avzeta=avzeta/area
          if (first_time.eq.0) then
            first_time=1
            write(*,1)   'STEP', 'time[DAYS]',  'KINETIC_ENRG',
     &                'BAROTR_KE', 'MAX_ADV_CFL', 'MAX_VERT_CFL',
     &                                     'i_cx', 'j_cx', 'k_c'
C$   &                                                   , 'trd'
          endif
          write(check_line,2,iostat=ierr)  avke,  avke2b,  Cu_Adv,
     &                        Cu_W, i_cx_max,  j_cx_max,  k_cx_max
C$   &                                                 , proc(2)
  1       format(/1x,A,2x,A,1x,A,5x,A,8x,A,2(5x,A),3x,A,2(1x,A))
  2       format(ES18.11, ES17.10, 2F16.12, 2I7,I5, I3)
          ie=max_check_line
          do while (check_line(ie:ie).eq.' ' .and. ie.gt.0)
            ie=ie-1
          enddo
          is=1
          do while (check_line(is:is).eq.' ' .and. is.lt.ie)
            is=is+1
          enddo
          i=is-1
          do while (i.lt.ie)
            i=i+1
            if (check_line(i:i).eq.'E' .or.check_line(i:i).eq.'e') then
              check_line(i:ie-1)=check_line(i+1:ie)
              check_line(ie:ie)=' '
              ie=ie-1
            elseif ( check_line(i:i).lt.'0' .or.
     &               check_line(i:i).gt.'9' ) then
              if ( check_line(i:i).ne.' ' .and.
     &             check_line(i:i).ne.'+'  .and.
     &             check_line(i:i).ne.'-'  .and.
     &             check_line(i:i).ne.'.') then
                nwrt=1
                if (may_day_flag.eq.0) then
                  write(*,*) '---> setting may day flag'
                  may_day_flag=1
                endif
              endif
            endif
          enddo
          write(tstring,'(F18.8)') tdays
          i=1
          do while (i.lt.18 .and. tstring(i:i).eq.' ')
            i=i+1
          enddo
          write(*,'(I7,2(1x,A))') iic-1, tstring(i:i+8),
     &                                  check_line(is:ie)
          call flush(6)
        endif
        buff(1)=may_day_flag
        buff(2)=nwrt
        call MPI_Bcast(buff,2, MPI_REAL8, 0,ocean_grid_comm,ierr)
        may_day_flag=buff(1) +0.0001_8
        nwrt=buff(2)         +0.0001_8
      endif
C$OMP END CRITICAL (diag_cr_rgn)
      diag_sync=.true.
      return
      end
