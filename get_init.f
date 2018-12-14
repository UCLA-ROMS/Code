      subroutine get_init (req_rec, tindx)
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
      integer(kind=4) req_rec, tindx, max_rec,  record,   ncid, varid,
     &        ierr,    lstr,  start(2), count(2), ibuff(iaux),
     &        lvar,    latt,  lenstr,   nf_fread, checkdims
     &                                          , itrc
      real(kind=8) time_scale
      integer(kind=4) init_type
      integer(kind=4), parameter :: init_run=1,  rst_run=2
     &                    , apprx_rst=3, exact_rst=4
      real(kind=8) time_bak
      integer(kind=4) vid_DU, vid_DV
      character(len=64) units
      real(kind=8) zeta(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real(kind=8) ubar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real(kind=8) vbar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      common /ocean_zeta/zeta /ocean_ubar/ubar /ocean_vbar/vbar
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
        real(kind=8) c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1._8, c0=0.0_8,c2=2._8,
     &  c1000=1000._8,p5=0.5_8,
     &  spd = 86400.0_8,  dps = c1 / spd ,
     &  t0_kelvin= 273.16_8)
        real(kind=8)  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C
        parameter ( parm_Red_D_C_P  = 117.0_8,
     &  parm_Red_P_C_P  = 117.0_8,
     &  parm_Red_D_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_P_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_D_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_P_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_Fe_C   = 3.0D-6)
        real(kind=8), parameter :: Q         = 0.137_8
        real(kind=8), parameter :: Qinv      = 1.0_8/0.137_8
        real(kind=8), parameter :: Qp        = 0.00855_8
        real(kind=8), parameter :: Qp_diat   = 0.1_8 * Q
        real(kind=8), parameter :: Qdenit    = 0.9811320755_8
        real(kind=8)  parm_Fe_bioavail,   parm_prod_dissolve,
     &  parm_o2_min,     parm_Rain_CaCO3, parm_Rain_SiO2,
     &  parm_kappa_nitrif,  parm_nitrif_par_lim,  parm_POC_flux_ref,
     &  parm_rest_prod_tau,  parm_rest_prod_z_c,  parm_z_umax_0,
     &  parm_diat_umax_0,     parm_z_mort_0,     parm_z_mort2_0,
     &  parm_sd_remin_0,     parm_sp_kNO3,       parm_diat_kNO3,
     &  parm_sp_kNH4,        parm_diat_kNH4,     parm_sp_kFe,
     &  parm_diat_kFe,       parm_diat_kSiO3,    parm_sp_kPO4,
     &  parm_diat_kPO4,      parm_z_grz,        parm_alphaChl,
     &  parm_labile_ratio,   parm_alphaDiaz,    parm_diaz_umax_0,
     &  gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
       common/eco_para/parm_Fe_bioavail, parm_prod_dissolve,
     &   parm_o2_min, parm_Rain_CaCO3, parm_Rain_SiO2,
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref,
     &   parm_rest_prod_tau, parm_rest_prod_z_c, parm_z_umax_0,
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0,
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe,
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4,
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl,
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
        real(kind=8) parm_kappa_nitrif_no2
        real(kind=8) parm_n2o_prod
        real(kind=8) parm_n2_prod
        common /eco_param_ncycle_anoxic/ parm_kappa_nitrif_no2,
     &       parm_n2o_prod,parm_n2_prod
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
      if (may_day_flag.ne.0) return
      init_type=0
      lstr=lenstr(ininame)
      ierr=nf_open(ininame(1:lstr), nf_nowrite, ncid)
      if (ierr .eq. nf_noerr) then
        ierr=checkdims (ncid, ininame, max_rec)
        if (ierr. eq. nf_noerr) then
          if (max_rec.gt.0) then
            if (req_rec.gt.0) then
              if (req_rec.le.max_rec) then
                record=req_rec
              else
                write(*,'(/1x,2A,I4,1x,A/12x,A,I4,1x,3A/)')
     &              '### ERROR: get_init :: requested restart time ',
     &              'record',req_rec, 'exceeds number', 'of records',
     &                       max_rec,  'available in netCDF file ''',
     &                                        ininame(1:lstr), '''.'
              endif
            else
              record=max_rec
            endif
          else
            record=1
          endif
        endif
      else
        write(*,'(/1x,4A/12x,A/)')  '### ERROR: get_init :: Cannot ',
     &                 'open netCDF file ''', ininame(1:lstr), '''.',
     &                                             nf_strerror(ierr)
      endif
      if (ierr .ne. nf_noerr) goto 99
      lvar=lenstr(vname(1,indxTime))
      ierr=nf_inq_varid (ncid, vname(1,indxTime)(1:lvar), varid)
      if (ierr .ne. nf_noerr) then
        ierr=nf_inq_varid (ncid, 'roms_time', varid)
      endif
      if (ierr .ne. nf_noerr) then
        ierr=nf_inq_varid (ncid, 'scrum_time', varid)
      endif
      if (ierr .eq. nf_noerr) then
        ierr=nf_get_var1_double (ncid, varid, record, start_time)
        if (ierr .eq. nf_noerr) then
          ierr=nf_get_att_text(ncid, varid, 'units', units)
          if (ierr .eq. nf_noerr) then
            latt=lenstr(units)
            if (units(1:6).eq.'second') then
              time_scale=1._8
            elseif (units(1:3).eq.'day') then
              time_scale=day2sec
            else
              write(*,'(/1x,4A/12x,3A/)')    '### ERROR: get_init :',
     &                          ': unknown units of for variable ''',
     &                               vname(1,indxTime)(1:lvar), '''',
     &                   'in netCDF file ''', ininame(1:lstr), '''.'
              goto 99
            endif
            start_time=start_time*time_scale
            if (tindx.eq.2) then
              forw_start=0
              if (record.lt.max_rec) then
                time_bak=start_time
                ierr=nf_get_var1_double (ncid, varid, record+1,
     &                                             start_time)
                if (ierr .eq. nf_noerr) then
                  start_time=start_time*time_scale
                 if (abs(start_time-time_bak-dt) .lt. 0.01_8*dt) then
                    forw_start=1
                  else
            if (mynode.eq.0) write(*,'(1x,2A,2I4/10x,4A/10x,A/)')
     &                'WARNING: Exact restart is requested, but ',
     &                'is not possible: records', record,record+1,
     &                'in ''',   ininame(1:lstr),   ''' are not ',
     &                'consecutive time steps ==> proceeding ',
     &                          'with forward initial time step.'
                  endif
                else
                  write(*,2) vname(1,indxTime)(1:lvar), record,
     &                       ininame(1:lstr), nf_strerror(ierr)
                  goto 99
                endif
              else
                if (mynode.eq.0) write(*,'(1x,2A/10x,4A)')
     &           'WARNING: Exact restart is requested, but is not ',
     &           'possible: initial',  'file ''',   ininame(1:lstr),
     &                    ''' does not contain sufficient records.'
              endif
              if (forw_start.ne.1) return
              forw_start=0
            endif
          else
            write(*,'(/1x,4A/12x,3A/12x,A/)')  '### ERROR: get_init ',
     &           ':: Cannot read attribute ''units'' for variable ''',
     &          vname(1,indxTime)(1:lvar),'''', 'from netCDF file ''',
     &                    ininame(1:lstr),   '''.', nf_strerror(ierr)
          endif
        else
          write(*,2) vname(1,indxTime)(1:lvar), record,
     &               ininame(1:lstr), nf_strerror(ierr)
        endif
      else
        write(*,1) vname(1,indxTime)(1:lvar), ininame(1:lstr),
     &                                      nf_strerror(ierr)
      endif
      if (ierr .ne. nf_noerr)  goto 99
      ierr=nf_inq_varid (ncid, 'time_step', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_inq_vardimid (ncid, varid, ibuff)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_dimlen (ncid, ibuff(1), count(1))
          if (ierr .eq. nf_noerr) then
            start(1)=1
            start(2)=record
            count(2)=1
            ibuff(1:iaux)=0
            ierr=nf_get_vara_int (ncid, varid, start, count, ibuff)
            if (ierr .eq. nf_noerr) then
              ntstart=ibuff(1)+1
              nrecrst=ibuff(2)
              nrechis=ibuff(3)
              nrecavg=ibuff(4)
              if (tindx.eq.2 .and. record.lt.max_rec) then
                start(2)=record+1
                ierr=nf_get_vara_int (ncid, varid, start, count, ibuff)
                if (ierr .eq. nf_noerr) then
                  if (ibuff(1) .eq. ntstart) then
                    forw_start=1
                  else
                  if (mynode.eq.0) write(*,'(1x,2A,2I4/10x,4A/10x,A)')
     &             'WARNING: Exact restart is requested, but is not ',
     &               'possible: records',  record, record+1,  'in ''',
     &                ininame(1:lstr), ''' are not consecutive time ',
     &                                    'steps ==> proceeding with',
     &                                   'forward initial time step.'
                    return
                  endif
                else
                  write(*,2) 'time_step', record, ininame(1:lstr),
     &                                            nf_strerror(ierr)
                endif
              elseif (tindx.eq.1) then
                if (forw_start.eq.1) then
                  init_type=exact_rst
                else
                  init_type=apprx_rst
                endif
              endif
            else
              write(*,2)  'time_step',  record,  ininame(1:lstr),
     &                                           nf_strerror(ierr)
            endif
          else
             write(*,'(/1x,3A/16x,A/)')     '### ERROR: get_init :: ',
     &            'Cannot determine size of the first dimension for ',
     &                  'variable ''time_step'':', nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,3A/16x,A/)') '### ERROR: get_init :: Cannot ',
     &      'determine dimensions for netCDF variable ''time_step'':',
     &                                             nf_strerror(ierr)
        endif
        if (ierr .ne. nf_noerr)  goto 99
      else
        init_type=init_run
        ntstart=1
        nrecrst=0
        nrechis=0
        nrecavg=0
      endif
      if (tindx.eq.1 .and. forw_start.eq.0) forw_start=ntstart
      if (tindx.eq.1) then
        if (ntstart.gt.1) then
          if (mynode.eq.0) write(*,'(/4x,A,I5,A,I8/)')
     &         'Adjusting averaging phase ntsavg =',  ntsavg,
     &                          ' ---> ', ntsavg + ntstart-1
          ntsavg=ntsavg + ntstart-1
          if (mynode.eq.0) write(*,'(/4x,A,I5,A,I8/)')
     &         'Adjusting averaging phase nts_bgc_flux_avg =',
     &         nts_bgc_flux_avg,
     &         ' ---> ', nts_bgc_flux_avg + ntstart-1
          nts_bgc_flux_avg=nts_bgc_flux_avg + ntstart-1
        endif
      endif
      ierr=nf_inq_varid (ncid, 'fa_time_step', varid)
      if (ierr .eq. nf_noerr) then
        start(1)=1
        start(2)=record
        count(1)=4
        count(2)=1
        ierr=nf_get_vara_int (ncid, varid, start, count, ibuff)
        if (ierr .eq. nf_noerr) then
          nrec_bgc_flux_his=ibuff(1)
          if (mynode.eq.0) write(*,'(6x,A,I4)')
     &         'GET_INITIAL: bgc_flux_his rec =',
     &         nrec_bgc_flux_his
          nrec_bgc_flux_avg=ibuff(2)
          if (mynode.eq.0) write(*,'(6x,A,I4)')
     &         'GET_INITIAL: bgc_flux_avg rec =',
     &         nrec_bgc_flux_avg
          nrec_phys_flux_his=ibuff(3)
          write(*,'(6x,A,I4)')
     &         'GET_INITIAL: phys_flux_his rec =',
     &         nrec_phys_flux_his
          nrec_phys_flux_avg=ibuff(4)
          write(*,'(6x,A,I4)')
     &         'GET_INITIAL: phys_flux_avg rec =',
     &         nrec_phys_flux_avg
        else
          write(*,'(/1x,2A/)') 'GET_INITIAL WARNING: Cannot ',
     &          'read record indices for flux analysis.'
        endif
      endif
      lvar=lenstr(vname(1,indxZ))
      ierr=nf_inq_varid (ncid, vname(1,indxZ)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (zeta(-1,-1,1), ncid, varid,
     &                                         record, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(*,2) vname(1,indxZ)(1:lvar),   record,
     &               ininame(1:lstr), nf_strerror(ierr)
          goto 99
        endif
      else
        write(*,1) vname(1,indxZ)(1:lvar), ininame(1:lstr),
     &                                   nf_strerror(ierr)
        goto 99
      endif
      lvar=lenstr(vname(1,indxUb))
      ierr=nf_inq_varid (ncid, vname(1,indxUb)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (ubar(-1,-1,1), ncid, varid,
     &                                         record, u2dvar)
        if (ierr .ne. nf_noerr) then
          write(*,2) vname(1,indxUb)(1:lvar),   record,
     &               ininame(1:lstr), nf_strerror(ierr)
          goto 99
        endif
      else
        write(*,1) vname(1,indxUb)(1:lvar), ininame(1:lstr),
     &                                    nf_strerror(ierr)
        goto 99
      endif
      lvar=lenstr(vname(1,indxVb))
      ierr=nf_inq_varid (ncid, vname(1,indxVb)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (vbar(-1,-1,1), ncid, varid,
     &                                         record, v2dvar)
        if (ierr .ne. nf_noerr) then
          write(*,2) vname(1,indxVb)(1:lvar),   record,
     &               ininame(1:lstr), nf_strerror(ierr)
          goto 99
        endif
      else
        write(*,1) vname(1,indxVb)(1:lvar), ininame(1:lstr),
     &                                    nf_strerror(ierr)
        goto 99
      endif
      ierr=nf_inq_varid (ncid, 'DU_avg2', vid_DU)
      if (ierr .eq. nf_noerr) then
        ierr=nf_inq_varid (ncid, 'DV_avg2', vid_DV)
        if (ierr .eq. nf_noerr) then
          ierr=nf_fread (DU_avg2, ncid, vid_DU, record, u2dvar)
          if (ierr .eq. nf_noerr) then
            ierr=nf_fread (DV_avg2, ncid, vid_DV, record, v2dvar)
            if (ierr .ne. nf_noerr) then
              write(*,2) 'DV_avg2',  record, ininame(1:lstr),
     &                                       nf_strerror(ierr)
              goto 99
            endif
          else
            write(*,2)   'DU_avg2',  record, ininame(1:lstr),
     &                                       nf_strerror(ierr)
            goto 99
          endif
        else
          forw_start=ntstart
        endif
      else
        forw_start=ntstart
      endif
      lvar=lenstr(vname(1,indxU))
      ierr=nf_inq_varid (ncid, vname(1,indxU)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (u(-1,-1,1,tindx), ncid, varid,
     &                                         record, u3dvar)
        if (ierr .ne. nf_noerr) then
          write(*,2) vname(1,indxU)(1:lvar),    record,
     &               ininame(1:lstr), nf_strerror(ierr)
          goto 99
        endif
      else
        write(*,1) vname(1,indxU)(1:lvar), ininame(1:lstr),
     &                                   nf_strerror(ierr)
        goto 99
      endif
      lvar=lenstr(vname(1,indxV))
      ierr=nf_inq_varid (ncid, vname(1,indxV)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (v(-1,-1,1,tindx), ncid, varid,
     &                                         record, v3dvar)
        if (ierr .ne. nf_noerr) then
          write(*,2) vname(1,indxV)(1:lvar),    record,
     &               ininame(1:lstr), nf_strerror(ierr)
          goto 99
        endif
      else
        write(*,1) vname(1,indxV)(1:lvar), ininame(1:lstr),
     &                                   nf_strerror(ierr)
        goto 99
      endif
      do itrc=1,NT
        lvar=lenstr(vname(1,indxT+itrc-1))
        ierr=nf_inq_varid (ncid, vname(1,indxT+itrc-1)(1:lvar),
     &                                                    varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_fread (t(-1,-1,1,tindx,itrc), ncid,
     &                                    varid, record, r3dvar)
          if (ierr .ne. nf_noerr) then
            write(*,2) vname(1,indxT+itrc-1)(1:lvar), record,
     &                     ininame(1:lstr), nf_strerror(ierr)
            goto 99
          endif
        else
           print *, 'tracer ', itrc, ' made it here'
           if (itrc .eq. iNH4 .or. itrc .eq. iDOC .or.
     &          itrc .eq. iDON .or. itrc .eq. iDOFE .or.
     &          itrc .eq. iDOP
     &          .or. itrc .eq. iNO2 .or. itrc .eq. iN2O
     &          .or. itrc .eq. iN2
     &          ) then
              t(:,:,:,tindx,itrc) = 0.0_8
              write(*,3) vname(1,indxT+itrc-1)(1:lvar),
     &             ininame(1:lstr),nf_strerror(ierr)
           else
          write(*,1) vname(1,indxT+itrc-1)(1:lvar), ininame(1:lstr),
     &                                            nf_strerror(ierr)
          goto 99
           endif
        endif
      enddo
      do itrc=1,NT_sed
        lvar=lenstr(vname(1,indxSedFirst+itrc-1))
        ierr=nf_inq_varid (ncid, vname(1,indxSedFirst+itrc-1)(1:lvar),
     &       varid)
        if (ierr .eq. nf_noerr) then
          ierr=nf_fread (t_sed(-1,-1,itrc), ncid,  varid,
     &                                               record, r2dvar)
          if (ierr .ne. nf_noerr) then
            write(*,2) vname(1,indxSedFirst+itrc-1)(1:lvar),
     &            record,  ininame(1:lstr)
            goto 99
          endif
        else
          write(*,3) vname(1,indxSedFirst+itrc-1)(1:lvar),
     &                                    ininame(1:lstr)
          write(*,'(a,a)') 'Note: Previous versions of ROMS did ',
     &         'not include sediment biology.'
          write(*,'(a)') 'Start with 0 concentration.'
        endif
      enddo
      lvar=lenstr(vname(1,indxPH_rst))
      ierr=nf_inq_varid (ncid, vname(1,indxPH_rst)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
         ierr=nf_fread (ph_hist(-1,-1), ncid,  varid,
     &                                               record, r2dvar)
         if (ierr .ne. nf_noerr) then
            write(*,2) vname(1,indxPH_rst)(1:lvar), record,
     &           ininame(1:lstr)
            goto 99
         endif
      else
         if (mynode.eq.0) write(*,*)
     &        'pH field not in init/rst file, initialize to 0'
      endif
      lvar=lenstr(vname(1,indxPCO2_rst))
      ierr=nf_inq_varid (ncid, vname(1,indxPCO2_rst)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
         ierr=nf_fread (pCO2sw(-1,-1), ncid,  varid,
     &                                               record, r2dvar)
         if (ierr .ne. nf_noerr) then
            write(*,2) vname(1,indxPCO2_rst)(1:lvar), record,
     &           ininame(1:lstr)
            goto 99
         endif
      else
         if (mynode.eq.0) write(*,*)
     &        'PCO2 not in init/rst file, initialize to 0'
      endif
      lvar=lenstr(vname(1,indxPCO2air_rst))
      ierr=nf_inq_varid (ncid, vname(1,indxPCO2air_rst)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
         ierr=nf_fread (pCO2air(-1,-1), ncid,  varid,
     &                                               record, r2dvar)
         if (ierr .ne. nf_noerr) then
            write(*,2) vname(1,indxPCO2air_rst)(1:lvar), record,
     &           ininame(1:lstr)
            goto 99
         endif
      else
         if (mynode.eq.0) write(*,*)
     &        'PCO2air not in init/rst file, initialize to 0'
      endif
      lvar=lenstr(vname(1,indxPARinc_rst))
      ierr=nf_inq_varid (ncid, vname(1,indxPARinc_rst)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
         ierr=nf_fread (PARinc(-1,-1), ncid,  varid,
     &        record, r2dvar)
         if (ierr .ne. nf_noerr) then
            write(*,2) vname(1,indxPARinc_rst)(1:lvar), record,
     &           ininame(1:lstr)
            goto 99
         else
            PARinc_rst = PARinc
            PARinc_rst_read = 1
         endif
      else
         if (mynode.eq.0) write(*,*)
     &        'PARinc not in init/rst file, initialize to 0'
      endif
      lvar=lenstr(vname(1,indxPAR_rst))
      ierr=nf_inq_varid (ncid, vname(1,indxPAR_rst)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
         ierr=nf_fread (PAR(-1,-1,1), ncid,  varid,
     &                                               record, r3dvar)
         if (ierr .ne. nf_noerr) then
            write(*,2) vname(1,indxPAR_rst)(1:lvar), record,
     &           ininame(1:lstr)
            goto 99
         endif
      else
         if (mynode.eq.0) write(*,*)
     &        'PAR not in init/rst file, initialize to 0'
      endif
      lvar=lenstr(vname(1,indxHbls))
      ierr=nf_inq_varid (ncid, vname(1,indxHbls)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (hbls(-1,-1,tindx), ncid, varid,
     &                                           record, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(*,2) vname(1,indxHbls)(1:lvar), record,
     &               ininame(1:lstr), nf_strerror(ierr)
          goto 99
        endif
      else
        if (mynode.eq.0) write(*,'(1x,6A)')        'WARNING: netCDF ',
     & 'variable ''', vname(1,indxHbls)(1:lvar), ''' not found in ''',
     &           ininame(1:lstr), ''' ==> initialized to zero state.'
      endif
      lvar=lenstr(vname(1,indxHbbl))
      ierr=nf_inq_varid (ncid, vname(1,indxHbbl)(1:lvar), varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (hbbls(-1,-1,tindx), ncid, varid,
     &                                           record, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(*,2) vname(1,indxHbbl)(1:lvar), record,
     &               ininame(1:lstr), nf_strerror(ierr)
          goto 99
        endif
      else
        if (mynode.eq.0) write(*,'(1x,6A)')        'WARNING: netCDF ',
     & 'variable ''', vname(1,indxHbbl)(1:lvar), ''' not found in ''',
     &           ininame(1:lstr), ''' ==> initialized to zero state.'
      endif
      ierr=nf_close(ncid)
      if (mynode.eq.0) then
       if (tindx.eq.1) then
        if (init_type .eq. init_run) then
          write(*,'(6x,2A,F12.4,1x,A,I4)')   'get_init: Read initial ',
     &   'conditions for day =', start_time*sec2day, 'record =',record
        elseif (init_type .eq. exact_rst) then
          write(*,'(6x,A,F12.4,1x,A,I4,A,I8,A,I4,A,I4,A)')
     &   'get_init: Exact restart from day =',  start_time*sec2day,
     &   'rec =', record, '(',ntstart,',', nrecrst,',', nrechis,').'
        elseif (init_type .eq. apprx_rst) then
          write(*,'(6x,A,F12.4,1x,A,I4,A,I8,A,I4,A,I4,A)')
     &        'get_init: Approximate, single-step restart from day =',
     &         start_time*sec2day,  'rec =', record, '(',ntstart, ',',
     &                                     nrecrst, ',', nrechis,').'
        else
          write(*,'(/1x,A/)') '### ERROR: get_init :: Unknown Error.'
          goto 99
        endif
       endif
      endif
      return
  99  may_day_flag=2
      return
  1   format(/1x,'### ERROR: get_init :: Cannot find variable ''', A,
     &                 '''' /12x, 'in netCDF file ''', A, ''''/12x,A)
  2   format(/1x,'### ERROR: get_init :: Cannot read variable ''', A,
     &               '''' /12x, 'from netCDF file ''', A, ''''/12x,A)
  3   format(/1x,'### WARNING: get_init :: Cannot find variable ''',
     &             A,'''' /12x, ' from netCDF file ''',A, '''' /12x,
     &             ' ==> initialized to zero state.')
      end
