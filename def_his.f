      subroutine def_his (total_rec, ierr)
      implicit none
      logical create_new_file
      integer(kind=4) total_rec, ierr, rec, lncn,lvar,lenstr, timedim
     &      , r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2), checkdims
     &      , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4), itrc
      integer(kind=4) my_nf_def_dim
      character(len=60) text
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
      real(kind=4), parameter :: spv_set=1.D+33
      ierr=0
      lncn=lenstr(hisname)
      if (nrpfhis.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfhis))
        call insert_time_index (hisname, lncn, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefhis
      if (ncidhis .ne. -1)  create_new_file=.false.
  10  if (create_new_file) then
        ierr=nf_create (hisname(1:lncn), nf_clobber+nf_64bit_offset,
     &                                                        ncidhis)
        if (ierr .ne. nf_noerr) then
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_his :: ',
     &          'Cannot create netCDF file ''', hisname(1:lncn),
     &                              '''.',   nf_strerror(ierr)
          goto 99
        endif
        if (nrpfhis .eq. 0) total_rec=0
        call put_global_atts (ncidhis, ierr)
        ierr=my_nf_def_dim (ncidhis, 'xi_rho',  xi_rho,  r2dgrd(1))
        ierr=my_nf_def_dim (ncidhis, 'xi_u',    xi_u,    u2dgrd(1))
        ierr=my_nf_def_dim (ncidhis, 'eta_rho', eta_rho, r2dgrd(2))
        ierr=my_nf_def_dim (ncidhis, 'eta_v',   eta_v,   v2dgrd(2))
        ierr=my_nf_def_dim (ncidhis, 's_rho',   N,       r3dgrd(3))
        ierr=my_nf_def_dim (ncidhis, 's_w',     N+1,     w3dgrd(3))
        ierr=my_nf_def_dim (ncidhis, 'time', nf_unlimited, timedim)
        ierr=my_nf_def_dim (ncidhis, 'auxil',   iaux,     auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2)
        r3dgrd(4)=timedim
        u3dgrd(1)=u2dgrd(1)
        u3dgrd(2)=r2dgrd(2)
        u3dgrd(3)=r3dgrd(3)
        u3dgrd(4)=timedim
        v3dgrd(1)=r2dgrd(1)
        v3dgrd(2)=v2dgrd(2)
        v3dgrd(3)=r3dgrd(3)
        v3dgrd(4)=timedim
        w3dgrd(1)=r2dgrd(1)
        w3dgrd(2)=r2dgrd(2)
        w3dgrd(4)=timedim
        if (total_rec.le.1) call def_grid (ncidhis, r2dgrd)
        ierr=nf_def_var (ncidhis, 'time_step', nf_int, 2, auxil, 
     &                             hisTstep)
        ierr=nf_put_att_text (ncidhis, hisTstep, 'long_name', 48,
     &              'time step and record numbers from initialization')
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncidhis, vname(1,indxTime)(1:lvar), nf_double,
     &                                         1, timedim, hisTime)
        text=vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidhis, hisTime, 'long_name', lvar,
     &                                             text(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncidhis, hisTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        if (wrthis(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_def_var (ncidhis, vname(1,indxZ)(1:lvar), nf_float,
     &                                           3, r2dgrd, hisZ)
          text=vname(2,indxZ)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisZ, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxZ))
          ierr=nf_put_att_text (ncidhis, hisZ, 'units',     lvar,
     &                                  vname(3,indxZ)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisZ, '_FillValue', nf_float,
     &                                                 1, spv_set)
          if (ierr.ne.nf_noerr) then
             write(*,*)  'nf_put_att_XXX:', nf_strerror(ierr)
          endif
        endif
        if (wrthis(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_def_var (ncidhis, vname(1,indxUb)(1:lvar), nf_float,
     &                                           3, u2dgrd, hisUb)
          text=vname(2,indxUb)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisUb, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxUb))
          ierr=nf_put_att_text (ncidhis, hisUb, 'units',     lvar,
     &                                  vname(3,indxUb)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisUb, '_FillValue', nf_float,
     &                                                  1, spv_set)
        endif
        if (wrthis(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_def_var (ncidhis, vname(1,indxVb)(1:lvar), nf_float,
     &                                           3, v2dgrd, hisVb)
          text=vname(2,indxVb)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisVb, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxVb))
          ierr=nf_put_att_text (ncidhis, hisVb, 'units',     lvar,
     &                                  vname(3,indxVb)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisVb, '_FillValue', nf_float,
     &                                                  1, spv_set)
        endif
        if (wrthis(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_def_var (ncidhis, vname(1,indxU)(1:lvar), nf_float,
     &                                           4, u3dgrd, hisU)
          text=vname(2,indxU)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisU, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxU))
          ierr=nf_put_att_text (ncidhis, hisU, 'units',     lvar,
     &                                  vname(3,indxU)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisU, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        if (wrthis(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_def_var (ncidhis, vname(1,indxV)(1:lvar), nf_float,
     &                                           4, v3dgrd, hisV)
          text=vname(2,indxV)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisV, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxV))
          ierr=nf_put_att_text (ncidhis, hisV, 'units',     lvar,
     &                                  vname(3,indxV)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisV, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        do itrc=1,NT
          if (wrthis(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_def_var (ncidhis, vname(1,indxT+itrc-1)(1:lvar),
     &                             nf_float, 4, r3dgrd, hisT(itrc))
            text=vname(2,indxT+itrc-1)
            lvar=lenstr(text)
            ierr=nf_put_att_text (ncidhis, hisT(itrc), 'long_name',
     &                                          lvar, text(1:lvar))
            lvar=lenstr(vname(3,indxT+itrc-1))
            ierr=nf_put_att_text (ncidhis, hisT(itrc), 'units', lvar,
     &                               vname(3,indxT+itrc-1)(1:lvar))
            ierr=nf_put_att_real (ncidhis, hisT(itrc), '_FillValue',
     &                                      nf_float, 1, spv_set)
          endif
        enddo
        do itrc=1,NT_sed
          if (wrthis(indxSedFirst+itrc-1)) then
            lvar=lenstr(vname(1,indxSedFirst+itrc-1))
            ierr=nf_def_var (ncidhis,
     &           vname(1,indxSedFirst+itrc-1)(1:lvar),
     &           nf_float, 3, r2dgrd, hisTsed(itrc))
            lvar=lenstr(vname(2,indxSedFirst+itrc-1))
            ierr=nf_put_att_text (ncidhis, hisTsed(itrc), 'long_name',
     &           lvar, vname(2,indxSedFirst+itrc-1)(1:lvar))
            lvar=lenstr(vname(3,indxSedFirst+itrc-1))
            ierr=nf_put_att_text (ncidhis, hisTsed(itrc), 'units', lvar,
     &           vname(3,indxSedFirst+itrc-1)(1:lvar))
            ierr=nf_put_att_real (ncidhis, hisTsed(itrc), '_FillValue',
     &                                      nf_float, 1, spv_set)
          endif
        enddo
        if (wrthis(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_def_var (ncidhis, vname(1,indxR)(1:lvar), nf_float,
     &                                           4, r3dgrd, hisR)
          text=vname(2,indxR)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisR, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxR))
          ierr=nf_put_att_text (ncidhis, hisR, 'units',     lvar,
     &                                  vname(3,indxR)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisR, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        if (wrthis(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_def_var (ncidhis, vname(1,indxO)(1:lvar), nf_float,
     &                                           4, w3dgrd, hisO)
          text=vname(2,indxO)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisO, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxO))
          ierr=nf_put_att_text (ncidhis, hisO, 'units',     lvar,
     &                                  vname(3,indxO)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisO, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        if (wrthis(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_def_var (ncidhis, vname(1,indxW)(1:lvar), nf_float,
     &                                           4, r3dgrd, hisW)
          text=vname(2,indxW)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisW, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxW))
          ierr=nf_put_att_text (ncidhis, hisW, 'units',     lvar,
     &                                  vname(3,indxW)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisW, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        if (wrthis(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_def_var (ncidhis, vname(1,indxAkv)(1:lvar), nf_float,
     &                                           4, w3dgrd, hisAkv)
          text=vname(2,indxAkv)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisAkv, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkv))
          ierr=nf_put_att_text (ncidhis, hisAkv, 'units',     lvar,
     &                                  vname(3,indxAkv)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisAkv, '_FillValue', nf_float,
     &                                                   1, spv_set)
        endif
        if (wrthis(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_def_var (ncidhis, vname(1,indxAkt)(1:lvar), nf_float,
     &                                           4, w3dgrd, hisAkt)
          text=vname(2,indxAkt)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisAkt, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkt))
          ierr=nf_put_att_text (ncidhis, hisAkt, 'units',     lvar,
     &                                  vname(3,indxAkt)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisAkt, '_FillValue', nf_float,
     &                                                   1, spv_set)
        endif
        if (wrthis(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_def_var (ncidhis, vname(1,indxAks)(1:lvar), nf_float,
     &                                           4, w3dgrd, hisAks)
          text=vname(2,indxAks)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisAks, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAks))
          ierr=nf_put_att_text (ncidhis, hisAks, 'units',     lvar,
     &                                  vname(3,indxAks)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisAks, '_FillValue', nf_float,
     &                                                   1, spv_set)
        endif
        if (wrthis(indxHbls)) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_def_var (ncidhis, vname(1,indxHbls)(1:lvar),
     &                           nf_float, 3, r2dgrd, hisHbls)
          text=vname(2,indxHbls)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisHbls, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxHbls))
          ierr=nf_put_att_text (ncidhis, hisHbls, 'units',     lvar,
     &                                 vname(3,indxHbls)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisHbls, '_FillValue',
     &                                       nf_float, 1, spv_set)
        endif
        if (wrthis(indxHbbl)) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_def_var (ncidhis, vname(1,indxHbbl)(1:lvar),
     &                            nf_float, 3, r2dgrd, hisHbbl)
          text=vname(2,indxHbbl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidhis, hisHbbl, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxHbbl))
          ierr=nf_put_att_text (ncidhis, hisHbbl, 'units',     lvar,
     &                                 vname(3,indxHbbl)(1:lvar))
          ierr=nf_put_att_real (ncidhis, hisHbbl, '_FillValue',
     &                                       nf_float, 1, spv_set)
        endif
        lvar=lenstr(vname(1,indxPH_rst))
        ierr=nf_def_var (ncidhis, vname(1,indxPH_rst)(1:lvar), nf_float,
     &                                           3, r2dgrd, hisPH)
        text=vname(2,indxPH_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidhis, hisPH, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPH_rst))
        ierr=nf_put_att_text (ncidhis, hisPH, 'units',     lvar,
     &                                  vname(3,indxPH_rst)(1:lvar))
        ierr=nf_put_att_real (ncidhis, hisPH, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxPCO2_rst))
        ierr=nf_def_var (ncidhis, vname(1,indxPCO2_rst)(1:lvar), 
     &                             nf_float,
     &                                           3, r2dgrd, hisPCO2)
        text=vname(2,indxPCO2_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidhis, hisPCO2, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPCO2_rst))
        ierr=nf_put_att_text (ncidhis, hisPCO2, 'units',     lvar,
     &                                  vname(3,indxPCO2_rst)(1:lvar))
        ierr=nf_put_att_real (ncidhis, hisPCO2, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxPCO2air_rst))
        ierr=nf_def_var (ncidhis, vname(1,indxPCO2air_rst)(1:lvar), 
     &                             nf_float,
     &                                           3, r2dgrd, hisPCO2air)
        text=vname(2,indxPCO2air_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidhis, hisPCO2air, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPCO2air_rst))
        ierr=nf_put_att_text (ncidhis, hisPCO2air, 'units',     lvar,
     &       vname(3,indxPCO2air_rst)(1:lvar))
        ierr=nf_put_att_real (ncidhis, hisPCO2air, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxPARinc_rst))
        ierr=nf_def_var (ncidhis, vname(1,indxPARinc_rst)(1:lvar), 
     &                             nf_float,
     &       3, r2dgrd, hisPARinc)
        text=vname(2,indxPARinc_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidhis, hisPARinc, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPARinc_rst))
        ierr=nf_put_att_text (ncidhis, hisPARinc, 'units',     lvar,
     &                                  vname(3,indxPARinc_rst)(1:lvar))
        ierr=nf_put_att_real (ncidhis, hisPARinc, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxPAR_rst))
        ierr=nf_def_var (ncidhis, vname(1,indxPAR_rst)(1:lvar), 
     &                             nf_float,
     &                                           4, r3dgrd, hisPAR)
        text=vname(2,indxPAR_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidhis, hisPAR, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPAR_rst))
        ierr=nf_put_att_text (ncidhis, hisPAR, 'units',     lvar,
     &                                  vname(3,indxPAR_rst)(1:lvar))
        ierr=nf_put_att_real (ncidhis, hisPAR, '_FillValue',
     &                                  nf_float, 1, spv_set)
        ierr=nf_enddef(ncidhis)
        if (mynode.eq.0) then
         write(*,'(7x,4A,I4)') 'def_his :: created new netCDF file ''',
     &                          hisname(1:lncn), '''.'
        endif
      elseif (ncidhis.eq.-1) then
        ierr=nf_open (hisname(1:lncn), nf_write, ncidhis)
        if (ierr. eq. nf_noerr) then
          if (mynode.eq.0) then
            write(*,'(7x,5A,I4)')    'def_his :: opened existing ',
     &                     'file ''', hisname(1:lncn),  '''.'
          endif
          ierr=checkdims (ncidhis, hisname, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfhis.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfhis))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write(*,
     &                '(/1x,A,I5,1x,3A/21x,2(A,I5),1x,A/21x,A/)')
     &            'WARNING: def_his :: The actual number of records',
     &             rec, 'present in file ''',  hisname(1:lncn),  '''',
     &            'exceeds record', rec+1-ierr, '/', total_rec,
     &            'specified by restart initial conditions.',
     &        'All records beyond this number will be overwritten.'
              rec=rec-ierr
            elseif (nrpfhis.eq.0) then
              total_rec=rec+1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr. ne. nf_noerr) then
          create_new_file=.true.
          goto 10
        endif
        ierr=nf_inq_varid (ncidhis, 'time_step', hisTstep)
        if (ierr .ne. nf_noerr) then
          write(*,1) 'time_step', hisname(1:lncn)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncidhis,vname(1,indxTime)(1:lvar),hisTime)
        if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxTime)(1:lvar), hisname(1:lncn)
          goto 99
        endif
        if (wrthis(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_inq_varid (ncidhis, vname(1,indxZ)(1:lvar), hisZ)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxZ)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_inq_varid (ncidhis, vname(1,indxUb)(1:lvar), hisUb)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxUb)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_inq_varid (ncidhis, vname(1,indxVb)(1:lvar), hisVb)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxVb)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_inq_varid (ncidhis, vname(1,indxU)(1:lvar), hisU)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxU)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_inq_varid (ncidhis, vname(1,indxV)(1:lvar), hisV)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxV)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        do itrc=1,NT
          if (wrthis(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_inq_varid (ncidhis, vname(1,indxT+itrc-1)(1:lvar),
     &                                                  hisT(itrc))
            if (ierr .ne. nf_noerr) then
              write(*,1) vname(1,indxT+itrc-1)(1:lvar), hisname(1:lncn)
              goto 99
            endif
          endif
        enddo
        do itrc=1,NT_sed
          if (wrthis(indxSedFirst+itrc-1)) then
            lvar=lenstr(vname(1,indxSedFirst+itrc-1))
            ierr=nf_inq_varid (ncidhis,
     &           vname(1,indxSedFirst+itrc-1)(1:lvar),hisTsed(itrc))
            if (ierr .ne. nf_noerr) then
              write(*,1) vname(1,indxSedFirst+itrc-1)(1:lvar),
     &                                       hisname(1:lncn)
              goto 99
            endif
          endif
        enddo
        if (wrthis(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_inq_varid (ncidhis, vname(1,indxR)(1:lvar), hisR)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxR)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_inq_varid (ncidhis, vname(1,indxO)(1:lvar), hisO)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxO)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_inq_varid (ncidhis, vname(1,indxW)(1:lvar), hisW)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxW)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_inq_varid (ncidhis, vname(1,indxAkv)(1:lvar), hisAkv)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAkv)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_inq_varid (ncidhis,vname(1,indxAkt)(1:lvar), hisAkt)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAkt)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_inq_varid (ncidhis,vname(1,indxAks)(1:lvar), hisAks)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAks)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxHbls)) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_inq_varid (ncidhis,vname(1,indxHbls)(1:lvar), hisHbls)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxHbls)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        if (wrthis(indxHbbl)) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_inq_varid (ncidhis,vname(1,indxHbbl)(1:lvar), hisHbbl)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxHbbl)(1:lvar), hisname(1:lncn)
            goto 99
          endif
        endif
        lvar=lenstr(vname(1,indxPH_rst))
        ierr=nf_inq_varid (ncidhis, vname(1,indxPH_rst)(1:lvar), hisPH)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPH_rst)(1:lvar), hisname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxPCO2_rst))
        ierr=nf_inq_varid (ncidhis, vname(1,indxPCO2_rst)(1:lvar),
     &       hisPCO2)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPCO2_rst)(1:lvar), hisname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxPCO2air_rst))
        ierr=nf_inq_varid (ncidhis, vname(1,indxPCO2air_rst)(1:lvar),
     &       hisPCO2air)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPCO2air_rst)(1:lvar), hisname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxPARinc_rst))
        ierr=nf_inq_varid (ncidhis, vname(1,indxPARinc_rst)(1:lvar),
     &       hisPARinc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPARinc_rst)(1:lvar), hisname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxPAR_rst))
        ierr=nf_inq_varid (ncidhis, vname(1,indxPAR_rst)(1:lvar), 
     &                              hisPAR)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPAR_rst)(1:lvar), hisname(1:lncn)
           goto 99
        endif
        if (mynode.eq.0) then
          write(*,'(7x,4A,I4,2x,A,I4)') 'def_his :: opened existing ',
     &   'file ''', hisname(1:lncn), ''' from rec =', rec
        endif
      endif
      ierr=nf_set_fill (ncidhis, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(1x,4A,I4)') '### ERROR: def_his :: Cannot switch ',
     &          'to ''nf_nofill'' mode.', nf_strerror(ierr)
     &
      endif
   1  format(/1x,'### ERROR: def_his :: Cannot find variable ''',
     &                           A, ''' in file ''', A, '''.'/)
        if (total_rec.le.1) call wrt_grid (ncidhis, hisname, lncn)
  99  return
      end
      subroutine def_avg (total_rec, ierr)
      implicit none
      logical create_new_file
      integer(kind=4) total_rec, ierr, rec, lncn,lvar,lenstr, timedim
     &      , r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2), checkdims
     &      , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4), itrc
      integer(kind=4) my_nf_def_dim
      character(len=60) text
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
      real(kind=4), parameter :: spv_set=1.D+33
      ierr=0
      lncn=lenstr(avgname)
      if (nrpfavg.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfavg))
        call insert_time_index (avgname, lncn, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefhis
      if (ncidavg .ne. -1)  create_new_file=.false.
  10  if (create_new_file) then
        ierr=nf_create (avgname(1:lncn), nf_clobber+nf_64bit_offset,
     &                                                        ncidavg)
        if (ierr .ne. nf_noerr) then
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_avg :: ',
     &          'Cannot create netCDF file ''', avgname(1:lncn),
     &                              '''.',   nf_strerror(ierr)
          goto 99
        endif
        if (nrpfavg .eq. 0) total_rec=0
        call put_global_atts (ncidavg, ierr)
        ierr=my_nf_def_dim (ncidavg, 'xi_rho',  xi_rho,  r2dgrd(1))
        ierr=my_nf_def_dim (ncidavg, 'xi_u',    xi_u,    u2dgrd(1))
        ierr=my_nf_def_dim (ncidavg, 'eta_rho', eta_rho, r2dgrd(2))
        ierr=my_nf_def_dim (ncidavg, 'eta_v',   eta_v,   v2dgrd(2))
        ierr=my_nf_def_dim (ncidavg, 's_rho',   N,       r3dgrd(3))
        ierr=my_nf_def_dim (ncidavg, 's_w',     N+1,     w3dgrd(3))
        ierr=my_nf_def_dim (ncidavg, 'time', nf_unlimited, timedim)
        ierr=my_nf_def_dim (ncidavg, 'auxil',   iaux,     auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2)
        r3dgrd(4)=timedim
        u3dgrd(1)=u2dgrd(1)
        u3dgrd(2)=r2dgrd(2)
        u3dgrd(3)=r3dgrd(3)
        u3dgrd(4)=timedim
        v3dgrd(1)=r2dgrd(1)
        v3dgrd(2)=v2dgrd(2)
        v3dgrd(3)=r3dgrd(3)
        v3dgrd(4)=timedim
        w3dgrd(1)=r2dgrd(1)
        w3dgrd(2)=r2dgrd(2)
        w3dgrd(4)=timedim
        if (total_rec.le.1) call def_grid (ncidavg, r2dgrd)
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncidavg, vname(1,indxTime)(1:lvar), nf_double,
     &                                         1, timedim, avgTime)
        text='averaged '/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgTime, 'long_name', lvar,
     &                                             text(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncidavg, avgTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        if (wrtavg(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_def_var (ncidavg, vname(1,indxZ)(1:lvar), nf_float,
     &                                           3, r2dgrd, avgZ)
          text='averaged '/ /vname(2,indxZ)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgZ, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxZ))
          ierr=nf_put_att_text (ncidavg, avgZ, 'units',     lvar,
     &                                  vname(3,indxZ)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgZ, '_FillValue', nf_float,
     &                                                 1, spv_set)
          if (ierr.ne.nf_noerr) then
             write(*,*)  'nf_put_att_XXX:', nf_strerror(ierr)
          endif
        endif
        if (wrtavg(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_def_var (ncidavg, vname(1,indxUb)(1:lvar), nf_float,
     &                                           3, u2dgrd, avgUb)
          text='averaged '/ /vname(2,indxUb)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgUb, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxUb))
          ierr=nf_put_att_text (ncidavg, avgUb, 'units',     lvar,
     &                                  vname(3,indxUb)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgUb, '_FillValue', nf_float,
     &                                                  1, spv_set)
        endif
        if (wrtavg(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_def_var (ncidavg, vname(1,indxVb)(1:lvar), nf_float,
     &                                           3, v2dgrd, avgVb)
          text='averaged '/ /vname(2,indxVb)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgVb, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxVb))
          ierr=nf_put_att_text (ncidavg, avgVb, 'units',     lvar,
     &                                  vname(3,indxVb)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgVb, '_FillValue', nf_float,
     &                                                  1, spv_set)
        endif
        if (wrtavg(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_def_var (ncidavg, vname(1,indxU)(1:lvar), nf_float,
     &                                           4, u3dgrd, avgU)
          text='averaged '/ /vname(2,indxU)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgU, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxU))
          ierr=nf_put_att_text (ncidavg, avgU, 'units',     lvar,
     &                                  vname(3,indxU)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgU, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        if (wrtavg(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_def_var (ncidavg, vname(1,indxV)(1:lvar), nf_float,
     &                                           4, v3dgrd, avgV)
          text='averaged '/ /vname(2,indxV)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgV, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxV))
          ierr=nf_put_att_text (ncidavg, avgV, 'units',     lvar,
     &                                  vname(3,indxV)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgV, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        do itrc=1,NT
          if (wrtavg(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_def_var (ncidavg, vname(1,indxT+itrc-1)(1:lvar),
     &                             nf_float, 4, r3dgrd, avgT(itrc))
            text='averaged '/ /vname(2,indxT+itrc-1)
            lvar=lenstr(text)
            ierr=nf_put_att_text (ncidavg, avgT(itrc), 'long_name',
     &                                          lvar, text(1:lvar))
            lvar=lenstr(vname(3,indxT+itrc-1))
            ierr=nf_put_att_text (ncidavg, avgT(itrc), 'units', lvar,
     &                               vname(3,indxT+itrc-1)(1:lvar))
            ierr=nf_put_att_real (ncidavg, avgT(itrc), '_FillValue',
     &                                      nf_float, 1, spv_set)
          endif
        enddo
        do itrc=1,NT_sed
          if (wrtavg(indxSedFirst+itrc-1)) then
            lvar=lenstr(vname(1,indxSedFirst+itrc-1))
            ierr=nf_def_var (ncidavg,
     &           vname(1,indxSedFirst+itrc-1)(1:lvar),
     &           nf_float, 3, r2dgrd, avgTsed(itrc))
            text='averaged '/ /vname(2,indxSedFirst+itrc-1)
            lvar=lenstr(text)
            ierr=nf_put_att_text (ncidavg, avgTsed(itrc), 'long_name',
     &           lvar, text(1:lvar))
            lvar=lenstr(vname(3,indxSedFirst+itrc-1))
            ierr=nf_put_att_text (ncidavg, avgTsed(itrc), 'units', lvar,
     &           vname(3,indxSedFirst+itrc-1)(1:lvar))
            ierr=nf_put_att_real (ncidavg, avgTsed(itrc), '_FillValue',
     &                                      nf_float, 1, spv_set)
          endif
        enddo
        if (wrtavg(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_def_var (ncidavg, vname(1,indxR)(1:lvar), nf_float,
     &                                           4, r3dgrd, avgR)
          text='averaged '/ /vname(2,indxR)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgR, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxR))
          ierr=nf_put_att_text (ncidavg, avgR, 'units',     lvar,
     &                                  vname(3,indxR)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgR, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        if (wrtavg(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_def_var (ncidavg, vname(1,indxO)(1:lvar), nf_float,
     &                                           4, w3dgrd, avgO)
          text='averaged '/ /vname(2,indxO)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgO, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxO))
          ierr=nf_put_att_text (ncidavg, avgO, 'units',     lvar,
     &                                  vname(3,indxO)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgO, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        if (wrtavg(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_def_var (ncidavg, vname(1,indxW)(1:lvar), nf_float,
     &                                           4, r3dgrd, avgW)
          text='averaged '/ /vname(2,indxW)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgW, 'long_name', lvar,
     &                                            text(1:lvar))
          lvar=lenstr(vname(3,indxW))
          ierr=nf_put_att_text (ncidavg, avgW, 'units',     lvar,
     &                                  vname(3,indxW)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgW, '_FillValue', nf_float,
     &                                                 1, spv_set)
        endif
        if (wrtavg(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_def_var (ncidavg, vname(1,indxAkv)(1:lvar), nf_float,
     &                                           4, w3dgrd, avgAkv)
          text='averaged '/ /vname(2,indxAkv)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgAkv, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkv))
          ierr=nf_put_att_text (ncidavg, avgAkv, 'units',     lvar,
     &                                  vname(3,indxAkv)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgAkv, '_FillValue', nf_float,
     &                                                   1, spv_set)
        endif
        if (wrtavg(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_def_var (ncidavg, vname(1,indxAkt)(1:lvar), nf_float,
     &                                           4, w3dgrd, avgAkt)
          text='averaged '/ /vname(2,indxAkt)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgAkt, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkt))
          ierr=nf_put_att_text (ncidavg, avgAkt, 'units',     lvar,
     &                                  vname(3,indxAkt)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgAkt, '_FillValue', nf_float,
     &                                                   1, spv_set)
        endif
        if (wrtavg(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_def_var (ncidavg, vname(1,indxAks)(1:lvar), nf_float,
     &                                           4, w3dgrd, avgAks)
          text='averaged '/ /vname(2,indxAks)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgAks, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAks))
          ierr=nf_put_att_text (ncidavg, avgAks, 'units',     lvar,
     &                                  vname(3,indxAks)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgAks, '_FillValue', nf_float,
     &                                                   1, spv_set)
        endif
        if (wrtavg(indxHbls)) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_def_var (ncidavg, vname(1,indxHbls)(1:lvar),
     &                           nf_float, 3, r2dgrd, avgHbls)
          text='averaged '/ /vname(2,indxHbls)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgHbls, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxHbls))
          ierr=nf_put_att_text (ncidavg, avgHbls, 'units',     lvar,
     &                                 vname(3,indxHbls)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgHbls, '_FillValue',
     &                                       nf_float, 1, spv_set)
        endif
        if (wrtavg(indxHbbl)) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_def_var (ncidavg, vname(1,indxHbbl)(1:lvar),
     &                            nf_float, 3, r2dgrd, avgHbbl)
          text='averaged '/ /vname(2,indxHbbl)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncidavg, avgHbbl, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxHbbl))
          ierr=nf_put_att_text (ncidavg, avgHbbl, 'units',     lvar,
     &                                 vname(3,indxHbbl)(1:lvar))
          ierr=nf_put_att_real (ncidavg, avgHbbl, '_FillValue',
     &                                       nf_float, 1, spv_set)
        endif
        lvar=lenstr(vname(1,indxPH_rst))
        ierr=nf_def_var (ncidavg, vname(1,indxPH_rst)(1:lvar), nf_float,
     &                                           3, r2dgrd, avgPH)
        text='averaged '/ /vname(2,indxPH_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgPH, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPH_rst))
        ierr=nf_put_att_text (ncidavg, avgPH, 'units',     lvar,
     &                                  vname(3,indxPH_rst)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgPH, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxPCO2_rst))
        ierr=nf_def_var (ncidavg, vname(1,indxPCO2_rst)(1:lvar), 
     &                             nf_float,
     &                                           3, r2dgrd, avgPCO2)
        text='averaged '/ /vname(2,indxPCO2_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgPCO2, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPCO2_rst))
        ierr=nf_put_att_text (ncidavg, avgPCO2, 'units',     lvar,
     &                                  vname(3,indxPCO2_rst)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgPCO2, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxPCO2air_rst))
        ierr=nf_def_var (ncidavg, vname(1,indxPCO2air_rst)(1:lvar), 
     &                             nf_float,
     &                                           3, r2dgrd, avgPCO2air)
        text='averaged '/ /vname(2,indxPCO2air_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgPCO2air, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPCO2air_rst))
        ierr=nf_put_att_text (ncidavg, avgPCO2air, 'units',     lvar,
     &       vname(3,indxPCO2air_rst)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgPCO2air, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxPARinc_rst))
        ierr=nf_def_var (ncidavg, vname(1,indxPARinc_rst)(1:lvar), 
     &                             nf_float,
     &       3, r2dgrd, avgPARinc)
        text='averaged '/ /vname(2,indxPARinc_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgPARinc, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPARinc_rst))
        ierr=nf_put_att_text (ncidavg, avgPARinc, 'units',     lvar,
     &                                  vname(3,indxPARinc_rst)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgPARinc, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxPAR_rst))
        ierr=nf_def_var (ncidavg, vname(1,indxPAR_rst)(1:lvar), 
     &                             nf_float,
     &                                           4, r3dgrd, avgPAR)
        text='averaged '/ /vname(2,indxPAR_rst)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgPAR, 'long_name', lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxPAR_rst))
        ierr=nf_put_att_text (ncidavg, avgPAR, 'units',     lvar,
     &                                  vname(3,indxPAR_rst)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgPAR, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxFreqDomSP_sfc))
        ierr=nf_def_var (ncidavg, vname(1,indxFreqDomSP_sfc)(1:lvar),
     &       nf_float, 3, r2dgrd, avgFreqDomSP_sfc)
        text='averaged '/ /vname(2,indxFreqDomSP_sfc)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgFreqDomSP_sfc, 'long_name', 
     &                               lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomSP_sfc))
        ierr=nf_put_att_text (ncidavg, avgFreqDomSP_sfc, 'units',     
     &                               lvar,
     &       vname(3,indxFreqDomSP_sfc)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgFreqDomSP_sfc, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxFreqDomDIAT_sfc))
        ierr=nf_def_var (ncidavg, vname(1,indxFreqDomDIAT_sfc)(1:lvar),
     &       nf_float, 3, r2dgrd, avgFreqDomDIAT_sfc)
        text='averaged '/ /vname(2,indxFreqDomDIAT_sfc)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgFreqDomDIAT_sfc, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomDIAT_sfc))
        ierr=nf_put_att_text (ncidavg, avgFreqDomDIAT_sfc, 'units', 
     &                               lvar,
     &       vname(3,indxFreqDomDIAT_sfc)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgFreqDomDIAT_sfc, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxFreqDomDIAZ_sfc))
        ierr=nf_def_var (ncidavg, vname(1,indxFreqDomDIAZ_sfc)(1:lvar),
     &       nf_float, 3, r2dgrd, avgFreqDomDIAZ_sfc)
        text='averaged '/ /vname(2,indxFreqDomDIAZ_sfc)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgFreqDomDIAZ_sfc, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomDIAZ_sfc))
        ierr=nf_put_att_text (ncidavg, avgFreqDomDIAZ_sfc, 'units', 
     &                               lvar,
     &       vname(3,indxFreqDomDIAZ_sfc)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgFreqDomDIAZ_sfc, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxFreqDomSP_int))
        ierr=nf_def_var (ncidavg, vname(1,indxFreqDomSP_int)(1:lvar),
     &       nf_float, 3, r2dgrd, avgFreqDomSP_int)
        text='averaged '/ /vname(2,indxFreqDomSP_int)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgFreqDomSP_int, 'long_name', 
     &                               lvar,
     &                                              text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomSP_int))
        ierr=nf_put_att_text (ncidavg, avgFreqDomSP_int, 'units',     
     &                               lvar,
     &       vname(3,indxFreqDomSP_int)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgFreqDomSP_int, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxFreqDomDIAT_int))
        ierr=nf_def_var (ncidavg, vname(1,indxFreqDomDIAT_int)(1:lvar),
     &       nf_float, 3, r2dgrd, avgFreqDomDIAT_int)
        text='averaged '/ /vname(2,indxFreqDomDIAT_int)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgFreqDomDIAT_int, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomDIAT_int))
        ierr=nf_put_att_text (ncidavg, avgFreqDomDIAT_int, 'units', 
     &                               lvar,
     &       vname(3,indxFreqDomDIAT_int)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgFreqDomDIAT_int, '_FillValue',
     &                                  nf_float, 1, spv_set)
        lvar=lenstr(vname(1,indxFreqDomDIAZ_int))
        ierr=nf_def_var (ncidavg, vname(1,indxFreqDomDIAZ_int)(1:lvar),
     &       nf_float, 3, r2dgrd, avgFreqDomDIAZ_int)
        text='averaged '/ /vname(2,indxFreqDomDIAZ_int)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncidavg, avgFreqDomDIAZ_int, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(3,indxFreqDomDIAZ_int))
        ierr=nf_put_att_text (ncidavg, avgFreqDomDIAZ_int, 'units', 
     &                               lvar,
     &       vname(3,indxFreqDomDIAZ_int)(1:lvar))
        ierr=nf_put_att_real (ncidavg, avgFreqDomDIAZ_int, '_FillValue',
     &                                  nf_float, 1, spv_set)
        ierr=nf_enddef(ncidavg)
        if (mynode.eq.0) then
         write(*,'(7x,4A,I4)') 'def_avg :: created new netCDF file ''',
     &                          avgname(1:lncn), '''.'
        endif
      elseif (ncidavg.eq.-1) then
        ierr=nf_open (avgname(1:lncn), nf_write, ncidavg)
        if (ierr. eq. nf_noerr) then
          if (mynode.eq.0) then
            write(*,'(7x,5A,I4)')    'def_avg :: opened existing ',
     &                     'file ''', avgname(1:lncn),  '''.'
          endif
          ierr=checkdims (ncidavg, avgname, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfavg.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfavg))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write(*,
     &                '(/1x,A,I5,1x,3A/21x,2(A,I5),1x,A/21x,A/)')
     &            'WARNING: def_avg :: The actual number of records',
     &             rec, 'present in file ''',  avgname(1:lncn),  '''',
     &            'exceeds record', rec+1-ierr, '/', total_rec,
     &            'specified by restart initial conditions.',
     &        'All records beyond this number will be overwritten.'
              rec=rec-ierr
            elseif (nrpfavg.eq.0) then
              total_rec=rec+1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr. ne. nf_noerr) then
          create_new_file=.true.
          goto 10
        endif
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncidavg,vname(1,indxTime)(1:lvar),avgTime)
        if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxTime)(1:lvar), avgname(1:lncn)
          goto 99
        endif
        if (wrtavg(indxZ)) then
          lvar=lenstr(vname(1,indxZ))
          ierr=nf_inq_varid (ncidavg, vname(1,indxZ)(1:lvar), avgZ)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxZ)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxUb)) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf_inq_varid (ncidavg, vname(1,indxUb)(1:lvar), avgUb)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxUb)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxVb)) then
          lvar=lenstr(vname(1,indxVb))
          ierr=nf_inq_varid (ncidavg, vname(1,indxVb)(1:lvar), avgVb)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxVb)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxU)) then
          lvar=lenstr(vname(1,indxU))
          ierr=nf_inq_varid (ncidavg, vname(1,indxU)(1:lvar), avgU)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxU)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxV)) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf_inq_varid (ncidavg, vname(1,indxV)(1:lvar), avgV)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxV)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        do itrc=1,NT
          if (wrtavg(indxT+itrc-1)) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            ierr=nf_inq_varid (ncidavg, vname(1,indxT+itrc-1)(1:lvar),
     &                                                  avgT(itrc))
            if (ierr .ne. nf_noerr) then
              write(*,1) vname(1,indxT+itrc-1)(1:lvar), avgname(1:lncn)
              goto 99
            endif
          endif
        enddo
        do itrc=1,NT_sed
          if (wrtavg(indxSedFirst+itrc-1)) then
            lvar=lenstr(vname(1,indxSedFirst+itrc-1))
            ierr=nf_inq_varid (ncidavg,
     &           vname(1,indxSedFirst+itrc-1)(1:lvar),avgTsed(itrc))
            if (ierr .ne. nf_noerr) then
              write(*,1) vname(1,indxSedFirst+itrc-1)(1:lvar),
     &                                       avgname(1:lncn)
              goto 99
            endif
          endif
        enddo
        if (wrtavg(indxR)) then
          lvar=lenstr(vname(1,indxR))
          ierr=nf_inq_varid (ncidavg, vname(1,indxR)(1:lvar), avgR)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxR)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxO)) then
          lvar=lenstr(vname(1,indxO))
          ierr=nf_inq_varid (ncidavg, vname(1,indxO)(1:lvar), avgO)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxO)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxW)) then
          lvar=lenstr(vname(1,indxW))
          ierr=nf_inq_varid (ncidavg, vname(1,indxW)(1:lvar), avgW)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxW)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxAkv)) then
          lvar=lenstr(vname(1,indxAkv))
          ierr=nf_inq_varid (ncidavg, vname(1,indxAkv)(1:lvar), avgAkv)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAkv)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxAkt)) then
          lvar=lenstr(vname(1,indxAkt))
          ierr=nf_inq_varid (ncidavg,vname(1,indxAkt)(1:lvar), avgAkt)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAkt)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxAks)) then
          lvar=lenstr(vname(1,indxAks))
          ierr=nf_inq_varid (ncidavg,vname(1,indxAks)(1:lvar), avgAks)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAks)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxHbls)) then
          lvar=lenstr(vname(1,indxHbls))
          ierr=nf_inq_varid (ncidavg,vname(1,indxHbls)(1:lvar), avgHbls)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxHbls)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        if (wrtavg(indxHbbl)) then
          lvar=lenstr(vname(1,indxHbbl))
          ierr=nf_inq_varid (ncidavg,vname(1,indxHbbl)(1:lvar), avgHbbl)
          if (ierr .ne. nf_noerr) then
          write(*,1) vname(1,indxHbbl)(1:lvar), avgname(1:lncn)
            goto 99
          endif
        endif
        lvar=lenstr(vname(1,indxPH_rst))
        ierr=nf_inq_varid (ncidavg, vname(1,indxPH_rst)(1:lvar), avgPH)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPH_rst)(1:lvar), avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxPCO2_rst))
        ierr=nf_inq_varid (ncidavg, vname(1,indxPCO2_rst)(1:lvar),
     &       avgPCO2)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPCO2_rst)(1:lvar), avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxPCO2air_rst))
        ierr=nf_inq_varid (ncidavg, vname(1,indxPCO2air_rst)(1:lvar),
     &       avgPCO2air)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPCO2air_rst)(1:lvar), avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxPARinc_rst))
        ierr=nf_inq_varid (ncidavg, vname(1,indxPARinc_rst)(1:lvar),
     &       avgPARinc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPARinc_rst)(1:lvar), avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxPAR_rst))
        ierr=nf_inq_varid (ncidavg, vname(1,indxPAR_rst)(1:lvar), 
     &                              avgPAR)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxPAR_rst)(1:lvar), avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxFreqDomSP_sfc))
        ierr=nf_inq_varid (ncidavg, vname(1,indxFreqDomSP_sfc)(1:lvar),
     &       avgFreqDomSP_sfc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomSP_sfc)(1:lvar), 
     &                          avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxFreqDomDIAT_sfc))
        ierr=nf_inq_varid (ncidavg, 
     &               vname(1,indxFreqDomDIAT_sfc)(1:lvar),
     &       avgFreqDomDIAT_sfc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAT_sfc)(1:lvar),
     &          avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxFreqDomDIAZ_sfc))
        ierr=nf_inq_varid (ncidavg, 
     &               vname(1,indxFreqDomDIAZ_sfc)(1:lvar),
     &       avgFreqDomDIAZ_sfc)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAZ_sfc)(1:lvar),
     &          avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxFreqDomSP_int))
        ierr=nf_inq_varid (ncidavg, vname(1,indxFreqDomSP_int)(1:lvar),
     &       avgFreqDomSP_int)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomSP_int)(1:lvar), 
     &                          avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxFreqDomDIAT_int))
        ierr=nf_inq_varid (ncidavg, 
     &               vname(1,indxFreqDomDIAT_int)(1:lvar),
     &       avgFreqDomDIAT_int)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAT_int)(1:lvar),
     &          avgname(1:lncn)
           goto 99
        endif
        lvar=lenstr(vname(1,indxFreqDomDIAZ_int))
        ierr=nf_inq_varid (ncidavg, 
     &               vname(1,indxFreqDomDIAZ_int)(1:lvar),
     &       avgFreqDomDIAZ_int)
        if (ierr .ne. nf_noerr) then
           write(*,1) vname(1,indxFreqDomDIAZ_int)(1:lvar),
     &          avgname(1:lncn)
           goto 99
        endif
        if (mynode.eq.0) then
          write(*,'(7x,4A,I4,2x,A,I4)') 'def_avg :: opened existing ',
     &   'file ''', avgname(1:lncn), ''' from rec =', rec
        endif
      endif
      ierr=nf_set_fill (ncidavg, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(1x,4A,I4)') '### ERROR: def_avg :: Cannot switch ',
     &          'to ''nf_nofill'' mode.', nf_strerror(ierr)
     &
      endif
   1  format(/1x,'### ERROR: def_avg :: Cannot find variable ''',
     &                           A, ''' in file ''', A, '''.'/)
        if (total_rec.le.1) call wrt_grid (ncidavg, avgname, lncn)
  99  return
      end
