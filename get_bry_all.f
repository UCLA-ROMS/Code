      subroutine get_bry_all (ierr)
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
      integer(kind=4) nf_byte
      integer(kind=4) nf_int1
      integer(kind=4) nf_char
      integer(kind=4) nf_short
      integer(kind=4) nf_int2
      integer(kind=4) nf_int
      integer(kind=4) nf_float
      integer(kind=4) nf_real
      integer(kind=4) nf_double
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      integer(kind=4)           nf_fill_byte
      integer(kind=4)           nf_fill_int1
      integer(kind=4)           nf_fill_char
      integer(kind=4)           nf_fill_short
      integer(kind=4)           nf_fill_int2
      integer(kind=4)           nf_fill_int
      real(kind=8)              nf_fill_float
      real(kind=8)              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer(kind=4) nf_nowrite
      integer(kind=4) nf_write
      integer(kind=4) nf_clobber
      integer(kind=4) nf_noclobber
      integer(kind=4) nf_fill
      integer(kind=4) nf_nofill
      integer(kind=4) nf_lock
      integer(kind=4) nf_share
      integer(kind=4) nf_64bit_offset
      integer(kind=4) nf_sizehint_default
      integer(kind=4) nf_align_chunk
      integer(kind=4) nf_format_classic
      integer(kind=4) nf_format_64bit
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      integer(kind=4) nf_unlimited
      parameter (nf_unlimited = 0)
      integer(kind=4) nf_global
      parameter (nf_global = 0)
      integer(kind=4) nf_max_dims
      integer(kind=4) nf_max_attrs
      integer(kind=4) nf_max_vars
      integer(kind=4) nf_max_name
      integer(kind=4) nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer(kind=4) nf_noerr
      integer(kind=4) nf_ebadid
      integer(kind=4) nf_eexist
      integer(kind=4) nf_einval
      integer(kind=4) nf_eperm
      integer(kind=4) nf_enotindefine
      integer(kind=4) nf_eindefine
      integer(kind=4) nf_einvalcoords
      integer(kind=4) nf_emaxdims
      integer(kind=4) nf_enameinuse
      integer(kind=4) nf_enotatt
      integer(kind=4) nf_emaxatts
      integer(kind=4) nf_ebadtype
      integer(kind=4) nf_ebaddim
      integer(kind=4) nf_eunlimpos
      integer(kind=4) nf_emaxvars
      integer(kind=4) nf_enotvar
      integer(kind=4) nf_eglobal
      integer(kind=4) nf_enotnc
      integer(kind=4) nf_ests
      integer(kind=4) nf_emaxname
      integer(kind=4) nf_eunlimit
      integer(kind=4) nf_enorecvars
      integer(kind=4) nf_echar
      integer(kind=4) nf_eedge
      integer(kind=4) nf_estride
      integer(kind=4) nf_ebadname
      integer(kind=4) nf_erange
      integer(kind=4) nf_enomem
      integer(kind=4) nf_evarsize
      integer(kind=4) nf_edimsize
      integer(kind=4) nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer(kind=4)  nf_fatal
      integer(kind=4) nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character(len=80)   nf_inq_libvers
      external       nf_inq_libvers
      character(len=80)   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer(kind=4)         nf_inq_base_pe
      external        nf_inq_base_pe
      integer(kind=4)         nf_set_base_pe
      external        nf_set_base_pe
      integer(kind=4)         nf_create
      external        nf_create
      integer(kind=4)         nf__create
      external        nf__create
      integer(kind=4)         nf__create_mp
      external        nf__create_mp
      integer(kind=4)         nf_open
      external        nf_open
      integer(kind=4)         nf__open
      external        nf__open
      integer(kind=4)         nf__open_mp
      external        nf__open_mp
      integer(kind=4)         nf_set_fill
      external        nf_set_fill
      integer(kind=4)         nf_set_default_format
      external        nf_set_default_format
      integer(kind=4)         nf_redef
      external        nf_redef
      integer(kind=4)         nf_enddef
      external        nf_enddef
      integer(kind=4)         nf__enddef
      external        nf__enddef
      integer(kind=4)         nf_sync
      external        nf_sync
      integer(kind=4)         nf_abort
      external        nf_abort
      integer(kind=4)         nf_close
      external        nf_close
      integer(kind=4)         nf_delete
      external        nf_delete
      integer(kind=4)         nf_inq
      external        nf_inq
      integer(kind=4)         nf_inq_ndims
      external        nf_inq_ndims
      integer(kind=4)         nf_inq_nvars
      external        nf_inq_nvars
      integer(kind=4)         nf_inq_natts
      external        nf_inq_natts
      integer(kind=4)         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer(kind=4)         nf_inq_format
      external        nf_inq_format
      integer(kind=4)         nf_def_dim
      external        nf_def_dim
      integer(kind=4)         nf_inq_dimid
      external        nf_inq_dimid
      integer(kind=4)         nf_inq_dim
      external        nf_inq_dim
      integer(kind=4)         nf_inq_dimname
      external        nf_inq_dimname
      integer(kind=4)         nf_inq_dimlen
      external        nf_inq_dimlen
      integer(kind=4)         nf_rename_dim
      external        nf_rename_dim
      integer(kind=4)         nf_inq_att
      external        nf_inq_att
      integer(kind=4)         nf_inq_attid
      external        nf_inq_attid
      integer(kind=4)         nf_inq_atttype
      external        nf_inq_atttype
      integer(kind=4)         nf_inq_attlen
      external        nf_inq_attlen
      integer(kind=4)         nf_inq_attname
      external        nf_inq_attname
      integer(kind=4)         nf_copy_att
      external        nf_copy_att
      integer(kind=4)         nf_rename_att
      external        nf_rename_att
      integer(kind=4)         nf_del_att
      external        nf_del_att
      integer(kind=4)         nf_put_att_text
      external        nf_put_att_text
      integer(kind=4)         nf_get_att_text
      external        nf_get_att_text
      integer(kind=4)         nf_put_att_int1
      external        nf_put_att_int1
      integer(kind=4)         nf_get_att_int1
      external        nf_get_att_int1
      integer(kind=4)         nf_put_att_int2
      external        nf_put_att_int2
      integer(kind=4)         nf_get_att_int2
      external        nf_get_att_int2
      integer(kind=4)         nf_put_att_int
      external        nf_put_att_int
      integer(kind=4)         nf_get_att_int
      external        nf_get_att_int
      integer(kind=4)         nf_put_att_real
      external        nf_put_att_real
      integer(kind=4)         nf_get_att_real
      external        nf_get_att_real
      integer(kind=4)         nf_put_att_double
      external        nf_put_att_double
      integer(kind=4)         nf_get_att_double
      external        nf_get_att_double
      integer(kind=4)         nf_def_var
      external        nf_def_var
      integer(kind=4)         nf_inq_var
      external        nf_inq_var
      integer(kind=4)         nf_inq_varid
      external        nf_inq_varid
      integer(kind=4)         nf_inq_varname
      external        nf_inq_varname
      integer(kind=4)         nf_inq_vartype
      external        nf_inq_vartype
      integer(kind=4)         nf_inq_varndims
      external        nf_inq_varndims
      integer(kind=4)         nf_inq_vardimid
      external        nf_inq_vardimid
      integer(kind=4)         nf_inq_varnatts
      external        nf_inq_varnatts
      integer(kind=4)         nf_rename_var
      external        nf_rename_var
      integer(kind=4)         nf_copy_var
      external        nf_copy_var
      integer(kind=4)         nf_put_var_text
      external        nf_put_var_text
      integer(kind=4)         nf_get_var_text
      external        nf_get_var_text
      integer(kind=4)         nf_put_var_int1
      external        nf_put_var_int1
      integer(kind=4)         nf_get_var_int1
      external        nf_get_var_int1
      integer(kind=4)         nf_put_var_int2
      external        nf_put_var_int2
      integer(kind=4)         nf_get_var_int2
      external        nf_get_var_int2
      integer(kind=4)         nf_put_var_int
      external        nf_put_var_int
      integer(kind=4)         nf_get_var_int
      external        nf_get_var_int
      integer(kind=4)         nf_put_var_real
      external        nf_put_var_real
      integer(kind=4)         nf_get_var_real
      external        nf_get_var_real
      integer(kind=4)         nf_put_var_double
      external        nf_put_var_double
      integer(kind=4)         nf_get_var_double
      external        nf_get_var_double
      integer(kind=4)         nf_put_var1_text
      external        nf_put_var1_text
      integer(kind=4)         nf_get_var1_text
      external        nf_get_var1_text
      integer(kind=4)         nf_put_var1_int1
      external        nf_put_var1_int1
      integer(kind=4)         nf_get_var1_int1
      external        nf_get_var1_int1
      integer(kind=4)         nf_put_var1_int2
      external        nf_put_var1_int2
      integer(kind=4)         nf_get_var1_int2
      external        nf_get_var1_int2
      integer(kind=4)         nf_put_var1_int
      external        nf_put_var1_int
      integer(kind=4)         nf_get_var1_int
      external        nf_get_var1_int
      integer(kind=4)         nf_put_var1_real
      external        nf_put_var1_real
      integer(kind=4)         nf_get_var1_real
      external        nf_get_var1_real
      integer(kind=4)         nf_put_var1_double
      external        nf_put_var1_double
      integer(kind=4)         nf_get_var1_double
      external        nf_get_var1_double
      integer(kind=4)         nf_put_vara_text
      external        nf_put_vara_text
      integer(kind=4)         nf_get_vara_text
      external        nf_get_vara_text
      integer(kind=4)         nf_put_vara_int1
      external        nf_put_vara_int1
      integer(kind=4)         nf_get_vara_int1
      external        nf_get_vara_int1
      integer(kind=4)         nf_put_vara_int2
      external        nf_put_vara_int2
      integer(kind=4)         nf_get_vara_int2
      external        nf_get_vara_int2
      integer(kind=4)         nf_put_vara_int
      external        nf_put_vara_int
      integer(kind=4)         nf_get_vara_int
      external        nf_get_vara_int
      integer(kind=4)         nf_put_vara_real
      external        nf_put_vara_real
      integer(kind=4)         nf_get_vara_real
      external        nf_get_vara_real
      integer(kind=4)         nf_put_vara_double
      external        nf_put_vara_double
      integer(kind=4)         nf_get_vara_double
      external        nf_get_vara_double
      integer(kind=4)         nf_put_vars_text
      external        nf_put_vars_text
      integer(kind=4)         nf_get_vars_text
      external        nf_get_vars_text
      integer(kind=4)         nf_put_vars_int1
      external        nf_put_vars_int1
      integer(kind=4)         nf_get_vars_int1
      external        nf_get_vars_int1
      integer(kind=4)         nf_put_vars_int2
      external        nf_put_vars_int2
      integer(kind=4)         nf_get_vars_int2
      external        nf_get_vars_int2
      integer(kind=4)         nf_put_vars_int
      external        nf_put_vars_int
      integer(kind=4)         nf_get_vars_int
      external        nf_get_vars_int
      integer(kind=4)         nf_put_vars_real
      external        nf_put_vars_real
      integer(kind=4)         nf_get_vars_real
      external        nf_get_vars_real
      integer(kind=4)         nf_put_vars_double
      external        nf_put_vars_double
      integer(kind=4)         nf_get_vars_double
      external        nf_get_vars_double
      integer(kind=4)         nf_put_varm_text
      external        nf_put_varm_text
      integer(kind=4)         nf_get_varm_text
      external        nf_get_varm_text
      integer(kind=4)         nf_put_varm_int1
      external        nf_put_varm_int1
      integer(kind=4)         nf_get_varm_int1
      external        nf_get_varm_int1
      integer(kind=4)         nf_put_varm_int2
      external        nf_put_varm_int2
      integer(kind=4)         nf_get_varm_int2
      external        nf_get_varm_int2
      integer(kind=4)         nf_put_varm_int
      external        nf_put_varm_int
      integer(kind=4)         nf_get_varm_int
      external        nf_get_varm_int
      integer(kind=4)         nf_put_varm_real
      external        nf_put_varm_real
      integer(kind=4)         nf_get_varm_real
      external        nf_get_varm_real
      integer(kind=4)         nf_put_varm_double
      external        nf_put_varm_double
      integer(kind=4)         nf_get_varm_double
      external        nf_get_varm_double
      integer(kind=4) nccre
      integer(kind=4) ncopn
      integer(kind=4) ncddef
      integer(kind=4) ncdid
      integer(kind=4) ncvdef
      integer(kind=4) ncvid
      integer(kind=4) nctlen
      integer(kind=4) ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer(kind=4) ncrdwr
      integer(kind=4) nccreat
      integer(kind=4) ncexcl
      integer(kind=4) ncindef
      integer(kind=4) ncnsync
      integer(kind=4) nchsync
      integer(kind=4) ncndirty
      integer(kind=4) nchdirty
      integer(kind=4) nclink
      integer(kind=4) ncnowrit
      integer(kind=4) ncwrite
      integer(kind=4) ncclob
      integer(kind=4) ncnoclob
      integer(kind=4) ncglobal
      integer(kind=4) ncfill
      integer(kind=4) ncnofill
      integer(kind=4) maxncop
      integer(kind=4) maxncdim
      integer(kind=4) maxncatt
      integer(kind=4) maxncvar
      integer(kind=4) maxncnam
      integer(kind=4) maxvdims
      integer(kind=4) ncnoerr
      integer(kind=4) ncebadid
      integer(kind=4) ncenfile
      integer(kind=4) nceexist
      integer(kind=4) nceinval
      integer(kind=4) nceperm
      integer(kind=4) ncenotin
      integer(kind=4) nceindef
      integer(kind=4) ncecoord
      integer(kind=4) ncemaxds
      integer(kind=4) ncename
      integer(kind=4) ncenoatt
      integer(kind=4) ncemaxat
      integer(kind=4) ncebadty
      integer(kind=4) ncebadd
      integer(kind=4) ncests
      integer(kind=4) nceunlim
      integer(kind=4) ncemaxvs
      integer(kind=4) ncenotvr
      integer(kind=4) nceglob
      integer(kind=4) ncenotnc
      integer(kind=4) ncfoobar
      integer(kind=4) ncsyserr
      integer(kind=4) ncfatal
      integer(kind=4) ncverbos
      integer(kind=4) ncentool
      integer(kind=4) ncbyte
      integer(kind=4) ncchar
      integer(kind=4) ncshort
      integer(kind=4) nclong
      integer(kind=4) ncfloat
      integer(kind=4) ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer(kind=4) ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer(kind=4) filbyte
      integer(kind=4) filchar
      integer(kind=4) filshort
      integer(kind=4) fillong
      real(kind=8) filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      integer(kind=4) nf_ubyte
      integer(kind=4) nf_ushort
      integer(kind=4) nf_uint
      integer(kind=4) nf_int64
      integer(kind=4) nf_uint64
      integer(kind=4) nf_string
      integer(kind=4) nf_vlen
      integer(kind=4) nf_opaque
      integer(kind=4) nf_enum
      integer(kind=4) nf_compound
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)
      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)
      integer(kind=4)           nf_fill_ubyte
      integer(kind=4)           nf_fill_ushort
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)
      integer(kind=4) nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)
      integer(kind=4) nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)
      integer(kind=4) nf_netcdf4
      parameter (nf_netcdf4 = 4096)
      integer(kind=4) nf_classic_model
      parameter (nf_classic_model = 256)
      integer(kind=4) nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer(kind=4) nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer(kind=4) nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)
      integer(kind=4) nf_endian_native
      parameter (nf_endian_native = 0)
      integer(kind=4) nf_endian_little
      parameter (nf_endian_little = 1)
      integer(kind=4) nf_endian_big
      parameter (nf_endian_big = 2)
      integer(kind=4) nf_chunked
      parameter (nf_chunked = 0)
      integer(kind=4) nf_contiguous
      parameter (nf_contiguous = 1)
      integer(kind=4) nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer(kind=4) nf_fletcher32
      parameter (nf_fletcher32 = 1)
      integer(kind=4) nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer(kind=4) nf_shuffle
      parameter (nf_shuffle = 1)
      integer(kind=4) nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer(kind=4) nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)
      integer(kind=4) nf_mpiio
      parameter (nf_mpiio = 8192)
      integer(kind=4) nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer(kind=4) nf_pnetcdf
      parameter (nf_pnetcdf = 32768)
      integer(kind=4) nf_independent
      parameter (nf_independent = 0)
      integer(kind=4) nf_collective
      parameter (nf_collective = 1)
      integer(kind=4) nf_ehdferr
      parameter (nf_ehdferr = -101)
      integer(kind=4) nf_ecantread
      parameter (nf_ecantread = -102)
      integer(kind=4) nf_ecantwrite
      parameter (nf_ecantwrite = -103)
      integer(kind=4) nf_ecantcreate
      parameter (nf_ecantcreate = -104)
      integer(kind=4) nf_efilemeta
      parameter (nf_efilemeta = -105)
      integer(kind=4) nf_edimmeta
      parameter (nf_edimmeta = -106)
      integer(kind=4) nf_eattmeta
      parameter (nf_eattmeta = -107)
      integer(kind=4) nf_evarmeta
      parameter (nf_evarmeta = -108)
      integer(kind=4) nf_enocompound
      parameter (nf_enocompound = -109)
      integer(kind=4) nf_eattexists
      parameter (nf_eattexists = -110)
      integer(kind=4) nf_enotnc4
      parameter (nf_enotnc4 = -111)
      integer(kind=4) nf_estrictnc3
      parameter (nf_estrictnc3 = -112)
      integer(kind=4) nf_enotnc3
      parameter (nf_enotnc3 = -113)
      integer(kind=4) nf_enopar
      parameter (nf_enopar = -114)
      integer(kind=4) nf_eparinit
      parameter (nf_eparinit = -115)
      integer(kind=4) nf_ebadgrpid
      parameter (nf_ebadgrpid = -116)
      integer(kind=4) nf_ebadtypid
      parameter (nf_ebadtypid = -117)
      integer(kind=4) nf_etypdefined
      parameter (nf_etypdefined = -118)
      integer(kind=4) nf_ebadfield
      parameter (nf_ebadfield = -119)
      integer(kind=4) nf_ebadclass
      parameter (nf_ebadclass = -120)
      integer(kind=4) nf_emaptype
      parameter (nf_emaptype = -121)
      integer(kind=4) nf_elatefill
      parameter (nf_elatefill = -122)
      integer(kind=4) nf_elatedef
      parameter (nf_elatedef = -123)
      integer(kind=4) nf_edimscale
      parameter (nf_edimscale = -124)
      integer(kind=4) nf_enogrp
      parameter (nf_enogrp = -125)
      integer(kind=4) nf_create_par
      external nf_create_par
      integer(kind=4) nf_open_par
      external nf_open_par
      integer(kind=4) nf_var_par_access
      external nf_var_par_access
      integer(kind=4) nf_inq_ncid
      external nf_inq_ncid
      integer(kind=4) nf_inq_grps
      external nf_inq_grps
      integer(kind=4) nf_inq_grpname
      external nf_inq_grpname
      integer(kind=4) nf_inq_grpname_full
      external nf_inq_grpname_full
      integer(kind=4) nf_inq_grpname_len
      external nf_inq_grpname_len
      integer(kind=4) nf_inq_grp_parent
      external nf_inq_grp_parent
      integer(kind=4) nf_inq_grp_ncid
      external nf_inq_grp_ncid
      integer(kind=4) nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid
      integer(kind=4) nf_inq_varids
      external nf_inq_varids
      integer(kind=4) nf_inq_dimids
      external nf_inq_dimids
      integer(kind=4) nf_def_grp
      external nf_def_grp
      integer(kind=4) nf_def_var_deflate
      external nf_def_var_deflate
      integer(kind=4) nf_inq_var_deflate
      external nf_inq_var_deflate
      integer(kind=4) nf_def_var_fletcher32
      external nf_def_var_fletcher32
      integer(kind=4) nf_inq_var_fletcher32
      external nf_inq_var_fletcher32
      integer(kind=4) nf_def_var_chunking
      external nf_def_var_chunking
      integer(kind=4) nf_inq_var_chunking
      external nf_inq_var_chunking
      integer(kind=4) nf_def_var_fill
      external nf_def_var_fill
      integer(kind=4) nf_inq_var_fill
      external nf_inq_var_fill
      integer(kind=4) nf_def_var_endian
      external nf_def_var_endian
      integer(kind=4) nf_inq_var_endian
      external nf_inq_var_endian
      integer(kind=4) nf_inq_typeids
      external nf_inq_typeids
      integer(kind=4) nf_inq_typeid
      external nf_inq_typeid
      integer(kind=4) nf_inq_type
      external nf_inq_type
      integer(kind=4) nf_inq_user_type
      external nf_inq_user_type
      integer(kind=4) nf_def_compound
      external nf_def_compound
      integer(kind=4) nf_insert_compound
      external nf_insert_compound
      integer(kind=4) nf_insert_array_compound
      external nf_insert_array_compound
      integer(kind=4) nf_inq_compound
      external nf_inq_compound
      integer(kind=4) nf_inq_compound_name
      external nf_inq_compound_name
      integer(kind=4) nf_inq_compound_size
      external nf_inq_compound_size
      integer(kind=4) nf_inq_compound_nfields
      external nf_inq_compound_nfields
      integer(kind=4) nf_inq_compound_field
      external nf_inq_compound_field
      integer(kind=4) nf_inq_compound_fieldname
      external nf_inq_compound_fieldname
      integer(kind=4) nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex
      integer(kind=4) nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset
      integer(kind=4) nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype
      integer(kind=4) nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims
      integer(kind=4) nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes
      integer(kind=4) nf_def_vlen
      external nf_def_vlen
      integer(kind=4) nf_inq_vlen
      external nf_inq_vlen
      integer(kind=4) nf_free_vlen
      external nf_free_vlen
      integer(kind=4) nf_def_enum
      external nf_def_enum
      integer(kind=4) nf_insert_enum
      external nf_insert_enum
      integer(kind=4) nf_inq_enum
      external nf_inq_enum
      integer(kind=4) nf_inq_enum_member
      external nf_inq_enum_member
      integer(kind=4) nf_inq_enum_ident
      external nf_inq_enum_ident
      integer(kind=4) nf_def_opaque
      external nf_def_opaque
      integer(kind=4) nf_inq_opaque
      external nf_inq_opaque
      integer(kind=4) nf_put_att
      external nf_put_att
      integer(kind=4) nf_get_att
      external nf_get_att
      integer(kind=4) nf_put_var
      external nf_put_var
      integer(kind=4) nf_put_var1
      external nf_put_var1
      integer(kind=4) nf_put_vara
      external nf_put_vara
      integer(kind=4) nf_put_vars
      external nf_put_vars
      integer(kind=4) nf_get_var
      external nf_get_var
      integer(kind=4) nf_get_var1
      external nf_get_var1
      integer(kind=4) nf_get_vara
      external nf_get_vara
      integer(kind=4) nf_get_vars
      external nf_get_vars
      integer(kind=4) nf_put_var1_int64
      external nf_put_var1_int64
      integer(kind=4) nf_put_vara_int64
      external nf_put_vara_int64
      integer(kind=4) nf_put_vars_int64
      external nf_put_vars_int64
      integer(kind=4) nf_put_varm_int64
      external nf_put_varm_int64
      integer(kind=4) nf_put_var_int64
      external nf_put_var_int64
      integer(kind=4) nf_get_var1_int64
      external nf_get_var1_int64
      integer(kind=4) nf_get_vara_int64
      external nf_get_vara_int64
      integer(kind=4) nf_get_vars_int64
      external nf_get_vars_int64
      integer(kind=4) nf_get_varm_int64
      external nf_get_varm_int64
      integer(kind=4) nf_get_var_int64
      external nf_get_var_int64
      integer(kind=4) nf_get_vlen_element
      external nf_get_vlen_element
      integer(kind=4) nf_put_vlen_element
      external nf_put_vlen_element
      integer(kind=4) nf_set_chunk_cache
      external nf_set_chunk_cache
      integer(kind=4) nf_get_chunk_cache
      external nf_get_chunk_cache
      integer(kind=4) nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer(kind=4) nf_get_var_chunk_cache
      external nf_get_var_chunk_cache
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
      real(kind=8) cff
      integer(kind=4) lstr, lenstr, ierr, ierr_all, itrc
     &                          , nf_read_bry_EW
     &                          , nf_read_bry_NS
      logical need_bry
     &      , west_bry_active
     &      , east_bry_active
     &      , south_bry_active
     &      , north_bry_active
      integer(kind=4) lvar
      character(len=40) vname_bry
      logical bry_set
      need_bry=.false.
      west_bry_active=.not.west_exchng
      need_bry=need_bry .or. west_bry_active
      east_bry_active=.not.east_exchng
      need_bry=need_bry .or. east_bry_active
      south_bry_active=.not.south_exchng
      need_bry=need_bry .or. south_bry_active
      north_bry_active=.not.north_exchng
      need_bry=need_bry .or. north_bry_active
      ierr_all=0
      ierr=nf_noerr
      lstr=lenstr(bry_file)
      if (iic.eq.ntstart) then
        if (bry_id .eq. -1) then
          ierr=nf_open (bry_file(1:lstr), nf_nowrite, bry_id)
          if (ierr.ne.nf_noerr) write(*,'(/1x,4A/12x,A/)')
     &             '### ERROR: get_bry_all :: Cannot open file ''',
     &              bry_file(1:lstr),   '''.',  nf_strerror(ierr)
        endif
        if (ierr.eq.nf_noerr) then
          ierr=nf_inq_varid (bry_id, 'bry_time',  bry_time_id)
          if (ierr.eq.nf_noerr) then
            if (west_bry_active) then
              ierr=nf_inq_varid (bry_id, 'zeta_west', zeta_west_id)
              if (ierr.ne.nf_noerr) write(*,1) 'zeta_west',
     &                                             bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id, 'ubar_west', ubar_west_id)
              if (ierr.ne.nf_noerr) write(*,1) 'ubar_west',
     &                                              bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id, 'vbar_west', vbar_west_id)
              if (ierr.ne.nf_noerr) write(*,1) 'vbar_west',
     &                                              bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,    'u_west', u_west_id)
              if (ierr.ne.nf_noerr)  write(*,1) 'u_west',
     &                                              bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,    'v_west', v_west_id)
              if (ierr.ne.nf_noerr)  write(*,1) 'v_west',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,'temp_west', t_west_id(itemp))
              if (ierr.ne.nf_noerr) write(*,1) 'temp_west',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,'salt_west',t_west_id(isalt))
              if (ierr.ne.nf_noerr) write(*,1) 'salt_west',
     &                                              bry_file(1:lstr)
              ierr_all=ierr_all+ierr
          do itrc = itemp + ntrc_salt + 1, NT
             lvar=lenstr(vname(1,itrc+indxT-itemp))
             vname_bry = vname(1,itrc+indxT-itemp)(1:lvar)//'_west'
             ierr=nf_inq_varid (bry_id, vname_bry,
     &            t_west_id(itrc))
             if (ierr.ne.nf_noerr) then
                if (mynode.eq.0) write(*,11) vname_bry(1:lvar+5),
     &                                             bry_file(1:lstr)
                t_west_id(itrc) = 0
             end if
          end do
            endif
            if (east_bry_active) then
              ierr=nf_inq_varid (bry_id, 'zeta_east', zeta_east_id)
              if (ierr.ne.nf_noerr) write(*,1) 'zeta_east',
     &                                             bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id, 'ubar_east', ubar_east_id)
              if (ierr.ne.nf_noerr) write(*,1) 'ubar_east',
     &                                              bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id, 'vbar_east', vbar_east_id)
              if (ierr.ne.nf_noerr) write(*,1) 'vbar_east',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,    'u_east',   u_east_id)
              if (ierr.ne.nf_noerr)  write(*,1)   'u_east',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,    'v_east',   v_east_id)
              if (ierr.ne.nf_noerr)  write(*,1)   'v_east',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,'temp_east',t_east_id(itemp))
              if (ierr.ne.nf_noerr) write(*,1) 'temp_east',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,'salt_east',t_east_id(isalt))
              if (ierr.ne.nf_noerr) write(*,1) 'salt_east',
     &                                              bry_file(1:lstr)
              ierr_all=ierr_all+ierr
          do itrc = itemp + ntrc_salt + 1, NT
             lvar=lenstr(vname(1,itrc+indxT-itemp))
             vname_bry = vname(1,itrc+indxT-itemp)(1:lvar)//'_east'
             ierr=nf_inq_varid (bry_id, vname_bry,
     &            t_east_id(itrc))
             if (ierr.ne.nf_noerr) then
                if (mynode.eq.0) write(*,11) vname_bry(1:lvar+5),
     &                                             bry_file(1:lstr)
                t_east_id(itrc) = 0
             end if
          end do
            endif
            if (south_bry_active) then
              ierr=nf_inq_varid (bry_id, 'zeta_south', zeta_south_id)
              if (ierr.ne.nf_noerr) write(*,1) 'zeta_south',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id, 'ubar_south', ubar_south_id)
              if (ierr.ne.nf_noerr) write(*,1) 'ubar_south',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id, 'vbar_south', vbar_south_id)
              if (ierr.ne.nf_noerr) write(*,1) 'vbar_south',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,    'u_south',  u_south_id)
              if (ierr.ne.nf_noerr) write(*,1)    'u_south',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,    'v_south',  v_south_id)
              if (ierr.ne.nf_noerr) write(*,1)    'v_south',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,'temp_south',t_south_id(itemp))
              if (ierr.ne.nf_noerr) write(*,1) 'temp_south',
     &                                                 bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,'salt_south',t_south_id(isalt))
              if (ierr.ne.nf_noerr) write(*,1) 'salt_south',
     &                                                bry_file(1:lstr)
              ierr_all=ierr_all+ierr
          do itrc = itemp + ntrc_salt + 1, NT
             lvar=lenstr(vname(1,itrc+indxT-itemp))
             vname_bry = vname(1,itrc+indxT-itemp)(1:lvar)//'_south'
             ierr=nf_inq_varid (bry_id, vname_bry,
     &            t_south_id(itrc))
             if (ierr.ne.nf_noerr) then
                if (mynode.eq.0) write(*,11) vname_bry(1:lvar+6),
     &                                             bry_file(1:lstr)
                t_south_id(itrc) = 0
             end if
          end do
            endif
            if (north_bry_active) then
              ierr=nf_inq_varid (bry_id, 'zeta_north', zeta_north_id)
              if (ierr.ne.nf_noerr) write(*,1) 'zeta_north',
     &                                              bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id, 'ubar_north', ubar_north_id)
              if (ierr.ne.nf_noerr) write(*,1) 'ubar_north',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id, 'vbar_north', vbar_north_id)
              if (ierr.ne.nf_noerr) write(*,1) 'vbar_north',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,    'u_north',  u_north_id)
              if (ierr.ne.nf_noerr) write(*,1)    'u_north',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,    'v_north',  v_north_id)
              if (ierr.ne.nf_noerr) write(*,1)  'v_north',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,'temp_north',t_north_id(itemp))
              if (ierr.ne.nf_noerr) write(*,1) 'temp_north',
     &                                               bry_file(1:lstr)
              ierr_all=ierr_all+ierr
              ierr=nf_inq_varid (bry_id,'salt_north',t_north_id(isalt))
              if (ierr.ne.nf_noerr) write(*,1) 'salt_north',
     &                                              bry_file(1:lstr)
              ierr_all=ierr_all+ierr
          do itrc = itemp + ntrc_salt + 1, NT
             lvar=lenstr(vname(1,itrc+indxT-itemp))
             vname_bry = vname(1,itrc+indxT-itemp)(1:lvar)//'_north'
             ierr=nf_inq_varid (bry_id, vname_bry,
     &            t_north_id(itrc))
             if (ierr.ne.nf_noerr) then
                if (mynode.eq.0) write(*,11) vname_bry(1:lvar+6),
     &                                             bry_file(1:lstr)
                t_north_id(itrc) = 0
             end if
          end do
            endif
            ierr=ierr_all
            if (ierr.eq.nf_noerr) then
              call set_cycle (bry_id, bry_time_id, ntbry,
     &             bry_cycle, bry_ncycle, bry_rec, ierr)
              itbry=1
              bry_time(1)=-1.D+20
              bry_time(2)=-1.D+20
            endif
          else
             write(*,1) 'bry_time', bry_file(1:lstr)
          endif
        endif
      endif
  1   format(1x,'### ERROR: get_bry_all :: Cannot find variable ''',
     &                                A, ''' in file ''', A, '''.')
  11   format(' NOTE in get_all_bry: cannot find variable ''',
     &                                A, ''' in file ''', A, '''.')
      do while (bry_time(itbry).lt.time+0.5_8*dt .and.ierr.eq.nf_noerr)
        call advance_cycle (bry_cycle,ntbry,bry_ncycle,bry_rec,ierr)
        if (ierr.eq.nf_noerr) then
          ierr=nf_get_var1_double (bry_id, bry_time_id, bry_rec, cff)
          if (ierr.eq.nf_noerr) then
            itbry=min(3-itbry,ntbry)
            bry_time(itbry)=cff*day2sec + bry_cycle*bry_ncycle
            if (need_bry) then
              if (west_bry_active) then
                ierr=nf_read_bry_EW (zeta_west_dt(0,itbry), bry_id,
     &                               zeta_west_id, bry_rec, r2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_EW (ubar_west_dt(0,itbry), bry_id,
     &                               ubar_west_id, bry_rec, u2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_EW (vbar_west_dt(0,itbry), bry_id,
     &                               vbar_west_id, bry_rec, v2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_EW (u_west_dt(0,1,itbry), bry_id,
     &                               u_west_id,  bry_rec,  u3dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_EW (v_west_dt(0,1,itbry), bry_id,
     &                               v_west_id,  bry_rec,  v3dvar)
                ierr_all=ierr_all+ierr
                do itrc=1,NT
                   bry_set = .false.
                   lvar=lenstr(vname(1,indxT+itrc-1))
                   vname_bry=vname(1,indxT+itrc-1)(1:lvar)//'_west'
                   if (t_west_id(itrc) .gt. 0) then
                      ierr=nf_read_bry_EW (t_west_dt(0,1,itbry,itrc),
     &                     bry_id, t_west_id(itrc),  bry_rec, r3dvar)
                      ierr_all=ierr_all+ierr
                      if (ierr.ne.nf_noerr)  then
                         if (itrc .le. 1 + ntrc_salt) then
                            write(*,1)  vname_bry(1:lvar+5)
                            ierr_all=ierr_all+ierr
                         end if
                      else
                         bry_set = .true.
                      end if
                   end if
               if (.not. bry_set .and. itrc .ge. iPO4 .and.
     &              itrc .le. iPO4 + ntrc_bio - 1) then
                  if ( (itrc .ge. iPO4 .and. itrc .le. iSIO3)
     &                 .or. (itrc .ge. iFE .and. itrc .le. iALK)
     &                 ) then
                     write(*,1)  vname_bry(1:lvar+5)
                     ierr_all=ierr_all+ierr
                  else
                     ierr = 0
                     if (mynode.eq.0) write(*,3) vname_bry(1:lvar+5)
                     if (itrc .eq. iNH4 .or. itrc .eq. iDOC .or.
     &                    itrc .eq. iDON .or. itrc. eq. iDOFE
     &                    .or. itrc .eq. iDOP
     &                    .or. itrc .eq. iNO2 .or. itrc .eq. iN2O
     &                    .or. itrc .eq. iN2
     &                    ) then
                        t_west_dt(:,:,itbry,itrc) = 0.0_8
                     else if (itrc .eq. iSPC) then
                        t_west_dt(:,:,itbry,itrc) = 0.10_8
                     else if (itrc .eq. iSPCHL) then
                        t_west_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_west_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iSPFE) then
                        t_west_dt(:,:,itbry,itrc) = 4.D-5 *
     &                     t_west_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iSPCACO3) then
                        t_west_dt(:,:,itbry,itrc) = 0.03_8 *
     &                     t_west_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iDIATC) then
                        t_west_dt(:,:,itbry,itrc) = 0.10_8
                     else if (itrc .eq. iDIATCHL) then
                        t_west_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_west_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIATSI) then
                        t_west_dt(:,:,itbry,itrc) = 0.2_8 *
     &                       t_west_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIATFE) then
                        t_west_dt(:,:,itbry,itrc) = 3.D-5 *
     &                       t_west_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIAZC) then
                        t_west_dt(:,:,itbry,itrc) = 0.01_8
                     else if (itrc .eq. iDIAZCHL) then
                        t_west_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_west_dt(:,:,itbry,iDIAZC)
                     else if (itrc .eq. iDIAZFE) then
                        t_west_dt(:,:,itbry,itrc) = 3.D-5 *
     &                       t_west_dt(:,:,itbry,iDIAZC)
                     else if (itrc .eq. iZOOC) then
                        t_west_dt(:,:,itbry,itrc) = 0.05_8
                     end if
                  end if
               end if
                enddo
              endif
              if (east_bry_active) then
                ierr=nf_read_bry_EW (zeta_east_dt(0,itbry), bry_id,
     &                             zeta_east_id, bry_rec, r2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_EW (ubar_east_dt(0,itbry), bry_id,
     &                               ubar_east_id, bry_rec, u2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_EW (vbar_east_dt(0,itbry), bry_id,
     &                               vbar_east_id, bry_rec, v2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_EW (u_east_dt(0,1,itbry), bry_id,
     &                               u_east_id,  bry_rec,  u3dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_EW (v_east_dt(0,1,itbry), bry_id,
     &                               v_east_id,  bry_rec,  v3dvar)
                ierr_all=ierr_all+ierr
                do itrc=1,NT
                   bry_set = .false.
                   lvar=lenstr(vname(1,indxT+itrc-1))
                   vname_bry = vname(1,indxT+itrc-1)(1:lvar)//'_east'
                   if (t_east_id(itrc) .gt. 0) then
                      ierr=nf_read_bry_EW (t_east_dt(0,1,itbry,itrc),
     &                     bry_id, t_east_id(itrc),  bry_rec, r3dvar)
                      ierr_all=ierr_all+ierr
                      if (ierr.ne.nf_noerr)  then
                         if (itrc .le. 1 + ntrc_salt) then
                            write(*,1)  vname_bry(1:lvar+5)
                            ierr_all=ierr_all+ierr
                         end if
                      else
                         bry_set = .true.
                      end if
                   end if
               if (.not. bry_set .and. itrc .ge. iPO4 .and.
     &              itrc .le. iPO4 + ntrc_bio - 1) then
                  if ( (itrc .ge. iPO4 .and. itrc .le. iSIO3)
     &                 .or. (itrc .ge. iFE .and. itrc .le. iALK)
     &                 ) then
                     write(*,1)  vname_bry(1:lvar+5)
                     ierr_all=ierr_all+ierr
                  else
                     ierr = 0
                     if (mynode.eq.0) write(*,3) vname_bry(1:lvar+5)
                     if (itrc .eq. iNH4 .or. itrc .eq. iDOC .or.
     &                    itrc .eq. iDON .or. itrc. eq. iDOFE
     &                    .or. itrc .eq. iDOP
     &                    .or. itrc .eq. iNO2 .or. itrc .eq. iN2O
     &                    .or. itrc .eq. iN2
     &               ) then
                        t_east_dt(:,:,itbry,itrc) = 0.0_8
                     else if (itrc .eq. iSPC) then
                        t_east_dt(:,:,itbry,itrc) = 0.10_8
                     else if (itrc .eq. iSPCHL) then
                        t_east_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_east_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iSPFE) then
                        t_east_dt(:,:,itbry,itrc) = 4.D-5 *
     &                     t_east_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iSPCACO3) then
                        t_east_dt(:,:,itbry,itrc) = 0.03_8 *
     &                     t_east_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iDIATC) then
                        t_east_dt(:,:,itbry,itrc) = 0.10_8
                     else if (itrc .eq. iDIATCHL) then
                        t_east_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_east_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIATSI) then
                        t_east_dt(:,:,itbry,itrc) = 0.2_8 *
     &                       t_east_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIATFE) then
                        t_east_dt(:,:,itbry,itrc) = 3.D-5 *
     &                       t_east_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIAZC) then
                        t_east_dt(:,:,itbry,itrc) = 0.01_8
                     else if (itrc .eq. iDIAZCHL) then
                        t_east_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_east_dt(:,:,itbry,iDIAZC)
                     else if (itrc .eq. iDIAZFE) then
                        t_east_dt(:,:,itbry,itrc) = 3.D-5 *
     &                       t_east_dt(:,:,itbry,iDIAZC)
                     else if (itrc .eq. iZOOC) then
                        t_east_dt(:,:,itbry,itrc) = 0.05_8
                     end if
                  end if
               end if
                enddo
              endif
              if (south_bry_active) then
                ierr=nf_read_bry_NS (zeta_south_dt(0,itbry), bry_id,
     &                              zeta_south_id, bry_rec,  r2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_NS (ubar_south_dt(0,itbry), bry_id,
     &                               ubar_south_id, bry_rec, u2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_NS (vbar_south_dt(0,itbry), bry_id,
     &                               vbar_south_id, bry_rec, v2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_NS (u_south_dt(0,1,itbry), bry_id,
     &                               u_south_id,  bry_rec,  u3dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_NS (v_south_dt(0,1,itbry), bry_id,
     &                               v_south_id,  bry_rec,  v3dvar)
                ierr_all=ierr_all+ierr
                do itrc=1,NT
                   bry_set = .false.
                   lvar=lenstr(vname(1,indxT+itrc-1))
                   vname_bry = vname(1,indxT+itrc-1)(1:lvar)//'_south'
                   if (t_south_id(itrc) .gt. 0) then
                      ierr=nf_read_bry_NS (t_south_dt(0,1,itbry,itrc),
     &                     bry_id, t_south_id(itrc),  bry_rec, r3dvar)
                      ierr_all=ierr_all+ierr
                      if (ierr.ne.nf_noerr)  then
                         if (itrc .le. 1 + ntrc_salt) then
                            write(*,1)  vname_bry(1:lvar+6)
                            ierr_all=ierr_all+ierr
                         end if
                      else
                         bry_set = .true.
                      end if
                   end if
               if (.not. bry_set .and. itrc .ge. iPO4 .and.
     &              itrc .le. iPO4 + ntrc_bio - 1) then
                  if ( (itrc .ge. iPO4 .and. itrc .le. iSIO3)
     &                 .or. (itrc .ge. iFE .and. itrc .le. iALK)
     &                 ) then
                     write(*,1)  vname_bry(1:lvar+6)
                     ierr_all=ierr_all+ierr
                  else
                     ierr = 0
                     if (mynode.eq.0) write(*,3) vname_bry(1:lvar+6)
                     if (itrc .eq. iNH4 .or. itrc .eq. iDOC .or.
     &                    itrc .eq. iDON .or. itrc. eq. iDOFE
     &                    .or. itrc .eq. iDOP
     &                    .or. itrc .eq. iNO2 .or. itrc .eq. iN2O
     &                    .or. itrc .eq. iN2
     &                    ) then
                        t_south_dt(:,:,itbry,itrc) = 0.0_8
                     else if (itrc .eq. iSPC) then
                        t_south_dt(:,:,itbry,itrc) = 0.10_8
                     else if (itrc .eq. iSPCHL) then
                        t_south_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_south_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iSPFE) then
                        t_south_dt(:,:,itbry,itrc) = 4.D-5 *
     &                     t_south_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iSPCACO3) then
                        t_south_dt(:,:,itbry,itrc) = 0.03_8 *
     &                     t_south_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iDIATC) then
                        t_south_dt(:,:,itbry,itrc) = 0.10_8
                     else if (itrc .eq. iDIATCHL) then
                        t_south_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_south_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIATSI) then
                        t_south_dt(:,:,itbry,itrc) = 0.2_8 *
     &                       t_south_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIATFE) then
                        t_south_dt(:,:,itbry,itrc) = 3.D-5 *
     &                       t_south_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIAZC) then
                        t_south_dt(:,:,itbry,itrc) = 0.01_8
                     else if (itrc .eq. iDIAZCHL) then
                        t_south_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_south_dt(:,:,itbry,iDIAZC)
                     else if (itrc .eq. iDIAZFE) then
                        t_south_dt(:,:,itbry,itrc) = 3.D-5 *
     &                       t_south_dt(:,:,itbry,iDIAZC)
                     else if (itrc .eq. iZOOC) then
                        t_south_dt(:,:,itbry,itrc) = 0.05_8
                     end if
                  end if
               end if
                enddo
              endif
              if (north_bry_active) then
                ierr=nf_read_bry_NS (zeta_north_dt(0,itbry), bry_id,
     &                               zeta_north_id, bry_rec, r2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_NS (ubar_north_dt(0,itbry), bry_id,
     &                               ubar_north_id, bry_rec, u2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_NS (vbar_north_dt(0,itbry), bry_id,
     &                               vbar_north_id, bry_rec, v2dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_NS (u_north_dt(0,1,itbry), bry_id,
     &                               u_north_id,  bry_rec,  u3dvar)
                ierr_all=ierr_all+ierr
                ierr=nf_read_bry_NS (v_north_dt(0,1,itbry), bry_id,
     &                               v_north_id,  bry_rec,  v3dvar)
                ierr_all=ierr_all+ierr
                do itrc=1,NT
                   bry_set = .false.
                   lvar=lenstr(vname(1,indxT+itrc-1))
                   vname_bry = vname(1,indxT+itrc-1)(1:lvar)//'_north'
                   if (t_north_id(itrc) .gt. 0) then
                      ierr=nf_read_bry_NS (t_north_dt(0,1,itbry,itrc),
     &                     bry_id, t_north_id(itrc),  bry_rec, r3dvar)
                      ierr_all=ierr_all+ierr
                      if (ierr.ne.nf_noerr)  then
                         if (itrc .le. 1 + ntrc_salt) then
                            write(*,1)  vname_bry(1:lvar+6)
                            ierr_all=ierr_all+ierr
                         end if
                      else
                         bry_set = .true.
                      end if
                   end if
               if (.not. bry_set .and. itrc .ge. iPO4 .and.
     &              itrc .le. iPO4 + ntrc_bio - 1) then
                  if ( (itrc .ge. iPO4 .and. itrc .le. iSIO3)
     &                 .or. (itrc .ge. iFE .and. itrc .le. iALK)
     &                 ) then
                     write(*,1)  vname_bry(1:lvar+6)
                     ierr_all=ierr_all+ierr
                  else
                     ierr = 0
                     if (mynode.eq.0) write(*,3) vname_bry(1:lvar+6)
                     if (itrc .eq. iNH4 .or. itrc .eq. iDOC .or.
     &                    itrc .eq. iDON .or. itrc. eq. iDOFE
     &                    .or. itrc .eq. iDOP
     &                    .or. itrc .eq. iNO2 .or. itrc .eq. iN2O
     &                    .or. itrc .eq. iN2
     &                    ) then
                        t_north_dt(:,:,itbry,itrc) = 0.0_8
                     else if (itrc .eq. iSPC) then
                        t_north_dt(:,:,itbry,itrc) = 0.10_8
                     else if (itrc .eq. iSPCHL) then
                        t_north_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_north_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iSPFE) then
                        t_north_dt(:,:,itbry,itrc) = 4.D-5 *
     &                     t_north_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iSPCACO3) then
                        t_north_dt(:,:,itbry,itrc) = 0.03_8 *
     &                     t_north_dt(:,:,itbry,iSPC)
                     else if (itrc .eq. iDIATC) then
                        t_north_dt(:,:,itbry,itrc) = 0.10_8
                     else if (itrc .eq. iDIATCHL) then
                        t_north_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_north_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIATSI) then
                        t_north_dt(:,:,itbry,itrc) = 0.2_8 *
     &                       t_north_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIATFE) then
                        t_north_dt(:,:,itbry,itrc) = 3.D-5 *
     &                       t_north_dt(:,:,itbry,iDIATC)
                     else if (itrc .eq. iDIAZC) then
                        t_north_dt(:,:,itbry,itrc) = 0.01_8
                     else if (itrc .eq. iDIAZCHL) then
                        t_north_dt(:,:,itbry,itrc) = 0.1_8 *
     &                       t_north_dt(:,:,itbry,iDIAZC)
                     else if (itrc .eq. iDIAZFE) then
                        t_north_dt(:,:,itbry,itrc) = 3.D-5 *
     &                       t_north_dt(:,:,itbry,iDIAZC)
                     else if (itrc .eq. iZOOC) then
                        t_north_dt(:,:,itbry,itrc) = 0.05_8
                     end if
                  end if
               end if
                enddo
              endif
              ierr=ierr_all
              if (ierr.eq.0) then
                if (mynode.eq.0) then
                  write(*,'(3x,A,3x,A,F12.4,2(1x,A,I4))')
     &           'get_bry_all :: read boundary for all fields',
     &           'bry_time =', cff, 'rec =', bry_rec
                endif
              else
                write(*,'(1x,2A,I4,1x,3A,1x,A,I4)') '### ERROR: ',
     &             'get_bry_all :: Cannot read rec =',  bry_rec,
     &             'of ''',  bry_file(1:lstr),  '''.'
              endif
            endif
          else
            write(*,'(8x,2A,I4)')  '### ERROR: get_bry_all :: ',
     &             'Cannot read variable ''bry_time'''
          endif
        else
           write(*,'(/1x,A,I4,1x,A,I4/7x,3A/7x,2(A,G12.4)/)')
     &   '### ERROR: get_bry_all :: requested time record ', bry_rec,
     &   'exeeds the last record',   ntbry,   'available in file ''',
     &    bry_file(1:lstr),  '''',  'tdays = ', tdays,
     &   '  but the last available  bry_time =',
     &                                   bry_time(itbry)*sec2day
        endif
      enddo
 3    format('  NOTE in get_all_bry: cannot read variable ''',A,'''/',
     &     ' - using default values')
      return
      end
      subroutine set_bry_all_tile (istr,iend,jstr,jend, ierr)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, ierr,
     &        imin,imax,jmin,jmax, it1,it2
     &                           , k, itrc
     &                           , j
     &                           , i
      logical need_bry
      real(kind=8) cff, cff1,cff2
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
      need_bry=.false.
      imin=istr
      if (istr.eq.iwest .and. .not.west_exchng) then
        need_bry=.true.
        imin=istr-1
      endif
      imax=iend
      if (iend.eq.ieast .and. .not.east_exchng) then
        need_bry=.true.
        imax=iend+1
      endif
      jmin=jstr
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        need_bry=.true.
        jmin=jstr-1
      endif
      jmax=jend
      if (jend.eq.jnorth .and. .not.north_exchng) then
        need_bry=.true.
        jmax=jend+1
      endif
      it1=3-itbry
      it2=itbry
      cff=time+0.5_8*dt
      cff1=bry_time(it2)-cff
      cff2=cff-bry_time(it1)
      if (cff1.lt.dt) synchro_flag=.true.
      if (need_bry) then
        cff=1._8/(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff
        if (istr.eq.iwest .and. .not.west_exchng) then
          do j=jmin,jmax
            zeta_west(j)=cff1*zeta_west_dt(j,it1)
     &                  +cff2*zeta_west_dt(j,it2)
          enddo
          do j=jmin,jmax
            ubar_west(j)=cff1*ubar_west_dt(j,it1)
     &                  +cff2*ubar_west_dt(j,it2)
            vbar_west(j)=cff1*vbar_west_dt(j,it1)
     &                  +cff2*vbar_west_dt(j,it2)
          enddo
          do k=1,N
            do j=jmin,jmax
              u_west(j,k)=cff1*u_west_dt(j,k,it1)
     &                   +cff2*u_west_dt(j,k,it2)
              v_west(j,k)=cff1*v_west_dt(j,k,it1)
     &                   +cff2*v_west_dt(j,k,it2)
            enddo
          enddo
          do itrc=1,NT
            do k=1,N
              do j=jmin,jmax
                t_west(j,k,itrc)=cff1*t_west_dt(j,k,it1,itrc)
     &                          +cff2*t_west_dt(j,k,it2,itrc)
              enddo
            enddo
          enddo
        endif
        if (iend.eq.ieast .and. .not.east_exchng) then
          do j=jmin,jmax
            zeta_east(j)=cff1*zeta_east_dt(j,it1)
     &                  +cff2*zeta_east_dt(j,it2)
          enddo
          do j=jmin,jmax
            ubar_east(j)=cff1*ubar_east_dt(j,it1)
     &                  +cff2*ubar_east_dt(j,it2)
            vbar_east(j)=cff1*vbar_east_dt(j,it1)
     &                  +cff2*vbar_east_dt(j,it2)
          enddo
          do k=1,N
            do j=jmin,jmax
              u_east(j,k)=cff1*u_east_dt(j,k,it1)
     &                   +cff2*u_east_dt(j,k,it2)
              v_east(j,k)=cff1*v_east_dt(j,k,it1)
     &                   +cff2*v_east_dt(j,k,it2)
            enddo
          enddo
          do itrc=1,NT
            do k=1,N
              do j=jmin,jmax
                t_east(j,k,itrc)=cff1*t_east_dt(j,k,it1,itrc)
     &                          +cff2*t_east_dt(j,k,it2,itrc)
              enddo
            enddo
          enddo
        endif
        if (jstr.eq.jsouth .and. .not.south_exchng) then
          do i=imin,imax
            zeta_south(i)=cff1*zeta_south_dt(i,it1)
     &                   +cff2*zeta_south_dt(i,it2)
          enddo
          do i=imin,imax
            ubar_south(i)=cff1*ubar_south_dt(i,it1)
     &                   +cff2*ubar_south_dt(i,it2)
            vbar_south(i)=cff1*vbar_south_dt(i,it1)
     &                   +cff2*vbar_south_dt(i,it2)
          enddo
          do k=1,N
            do i=imin,imax
              u_south(i,k)=cff1*u_south_dt(i,k,it1)
     &                    +cff2*u_south_dt(i,k,it2)
              v_south(i,k)=cff1*v_south_dt(i,k,it1)
     &                    +cff2*v_south_dt(i,k,it2)
            enddo
          enddo
          do itrc=1,NT
            do k=1,N
              do i=imin,imax
                t_south(i,k,itrc)=cff1*t_south_dt(i,k,it1,itrc)
     &                           +cff2*t_south_dt(i,k,it2,itrc)
              enddo
            enddo
          enddo
        endif
        if (jend.eq.jnorth .and. .not.north_exchng) then
          do i=imin,imax
            zeta_north(i)=cff1*zeta_north_dt(i,it1)
     &                   +cff2*zeta_north_dt(i,it2)
          enddo
          do i=imin,imax
            ubar_north(i)=cff1*ubar_north_dt(i,it1)
     &                   +cff2*ubar_north_dt(i,it2)
            vbar_north(i)=cff1*vbar_north_dt(i,it1)
     &                   +cff2*vbar_north_dt(i,it2)
          enddo
          do k=1,N
            do i=imin,imax
              u_north(i,k)=cff1*u_north_dt(i,k,it1)
     &                    +cff2*u_north_dt(i,k,it2)
              v_north(i,k)=cff1*v_north_dt(i,k,it1)
     &                    +cff2*v_north_dt(i,k,it2)
            enddo
          enddo
          do itrc=1,NT
            do k=1,N
              do i=imin,imax
                t_north(i,k,itrc)=cff1*t_north_dt(i,k,it1,itrc)
     &                           +cff2*t_north_dt(i,k,it2,itrc)
              enddo
            enddo
          enddo
        endif
        if (cff1.lt.0._8 .or. cff2.lt.0._8) then
           write(*,'(/3A/3(1x,A,F16.10)/)')            '### WARNING: ',
     &      'set_bry_all_tile :: current model time is out of bounds ',
     &      'of ''bry_time''.', 'bry_tstart =', bry_time(it1)*sec2day,
     &      'tdays =',  tdays,  'bry_tend =',   bry_time(it2)*sec2day
          ierr=ierr+1
        endif
      endif
      return
      end
