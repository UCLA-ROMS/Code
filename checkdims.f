      integer(kind=4) function checkdims (ncid, fname, max_rec)
      implicit none
      integer(kind=4) ncid, max_rec, ierr, icount, ndims, recdim,
     &                       dsize, id, ldim, lname, lenstr
      character(len=*) fname
      character(len=16) dname
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
      integer(kind=4) check_scoord, chk_Cs_wr
      chk_Cs_wr=0
      lname=lenstr(fname)
      icount=0
      recdim=0
      max_rec=0
      ierr=nf_inq_ndims (ncid, ndims)
      if (ierr .eq. nf_noerr) then
        ierr=nf_inq_unlimdim (ncid, recdim)
        if (ierr .ne. nf_noerr) then
          write(*,'(/1x,4A/)')        'WARNING: No unlimited ',
     &          'dimension found in ''', fname(1:lname), '''.'
        endif
        do id=1,ndims
          ierr=nf_inq_dim (ncid, id, dname, dsize)
          if (ierr .eq. nf_noerr) then
            ldim=lenstr(dname)
            if ((ldim.eq.6 .and. dname(1:ldim).eq.'xi_rho') .or.
     &          (ldim.eq.4 .and. dname(1:ldim).eq.'xi_v'  )) then
              if (dsize.ne.xi_rho) then
                write(*,1) dname(1:ldim), dsize,fname(1:lname),xi_rho
                icount=icount+1
              endif
            elseif ((ldim.eq.4 .and.dname(1:ldim).eq.'xi_u'  ) .or.
     &              (ldim.eq.6 .and. dname(1:ldim).eq.'xi_psi')) then
              if (dsize.ne.xi_u) then
                write(*,1) dname(1:ldim), dsize, fname(1:lname), xi_u
                icount=icount+1
              endif
            elseif ((ldim.eq.7 .and. dname(1:ldim).eq.'eta_rho') .or.
     &              (ldim.eq.5 .and. dname(1:ldim).eq.'eta_u' )) then
              if (dsize.ne.eta_rho) then
                write(*,1) dname(1:ldim), dsize,fname(1:lname),eta_rho
                icount=icount+1
              endif
            elseif ((ldim.eq.5 .and.dname(1:ldim).eq.'eta_v') .or.
     &              (ldim.eq.7 .and.dname(1:ldim).eq.'eta_psi')) then
              if (dsize.ne.eta_v) then
                write(*,1) dname(1:ldim), dsize, fname(1:lname), eta_v
                icount=icount+1
              endif
            elseif (ldim.eq.3 .and. dname(1:ldim).eq.'s_w') then
              chk_Cs_wr=chk_Cs_wr+1
              if (dsize.ne.N+1) then
                write(*,1) dname(1:ldim), dsize, fname(1:lname), N+1
                icount=icount+1
              endif
            elseif (ldim.eq.5 .and. dname(1:ldim).eq.'s_rho') then
              chk_Cs_wr=chk_Cs_wr+2
              if (dsize.ne.N) then
                write(*,1) dname(1:ldim), dsize, fname(1:lname), N
                icount=icount+1
              endif
            elseif (id.eq.recdim) then
              max_rec=dsize
            endif
          else
            write(*,'(/1x,2A,I3/12x,3A/12x,A/)')      '### ERROR: ',
     &            'checkdims :: Cannot get size of dimension #', id,
     &            'from netCDF file ''',    fname(1:lname),   '''.',
     &                                           nf_strerror(ierr)
            icount=icount+1
          endif
        enddo
      else
        write(*,'(/1x,4A/12x,A/)')       '### ERROR: checkdims :: ',
     &          'Cannot get number of dimensions in netCDF file ''',
     &                     fname(1:lname), '''.', nf_strerror(ierr)
        icount=icount+1
      endif
      if (chk_Cs_wr.gt.0) then
        icount=icount+check_scoord (ncid, fname, chk_Cs_wr)
      endif
  1   format(/' ### ERROR: checkdims :: wrong size of dimension ''',
     &                  A, ''' =', i5 / 12x, 'in netCDF file ''', A,
     &                            ''': must be', i5,1x, 'instead.'/)
      checkdims=icount
      return
      end
      integer(kind=4) function check_scoord (ncid, fname, what_to_check)
      implicit none
      integer(kind=4) ncid, what_to_check
      character(len=*) fname
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
      real(kind=8), parameter :: epsil=1.D-7
      real(kind=8) tst_val, Cs_w_try(0:N), Cs_r_try(N)
      integer(kind=4) ierr, k, icount, lfile, lenstr, read_nc1dat
      logical chk_Cs_w, chk_Cs_r
      if (mod(what_to_check,2).eq.1) then
        chk_Cs_w=.true.
      else
        chk_Cs_w=.false.
      endif
      if (what_to_check.gt.1) then
        chk_Cs_r=.true.
      else
        chk_Cs_r=.false.
      endif
      icount=0
      lfile=lenstr(fname)
      ierr=read_nc1dat(ncid, fname, 'hc', 1, tst_val)
      if (ierr.eq.0) then
        if (abs(hc-tst_val).gt.epsil) then
          write(*,'(/1x,2A,F12.5,1x,A,F12.5,1x,3A/)')  '### ERROR: ',
     &                  'Mismatch in ''hc'': should be', hc, 'found',
     &                   tst_val, 'in file ''', fname(1:lfile), '''.'
          icount=icount+1
        endif
      endif
      if (chk_Cs_w) then
        ierr=read_nc1dat(ncid, fname, 'Cs_w', N+1, Cs_w_try)
        if (ierr.eq.0) then
          do k=N,0,-1
            if (abs(Cs_w(k)-Cs_w_try(k)).gt.epsil) then
              ierr=ierr+1
            endif
          enddo
          if (ierr.ne.0) then
            write(*,'(/1x,4A/)')  '### ERROR: Mismatch in ''Cs_w'' ',
     &                    'values in file ''', fname(1:lfile), '''.'
            icount=icount+1
          endif
        endif
      endif
      if (chk_Cs_r) then
        ierr=read_nc1dat(ncid, fname, 'Cs_r', N, Cs_r_try)
        if (ierr.eq.0) then
          do k=N,1,-1
            if (abs(Cs_r(k)-Cs_r_try(k)).gt.epsil) then
              ierr=ierr+1
            endif
          enddo
          if (ierr.ne.0) then
            write(*,'(/1x,4A/)')  '### ERROR: Mismatch in ''Cs_r'' ',
     &                    'values in file ''', fname(1:lfile), '''.'
            icount=icount+1
          endif
        endif
      endif
      if (icount.ne.0) then
        ierr=read_nc1dat(ncid, fname, 'theta_s', 1,  tst_val)
        if (ierr.eq.0) then
          if (abs(theta_s-tst_val).gt.epsil) then
            write(*,'(/1x,A,2(A,F9.4,1x),3A/)')       '### ERROR: ',
     &              'Mismatch in ''theta_s'': should be',   theta_s,
     &        'found', tst_val, 'in file ''', fname(1:lfile), '''.'
          endif
        endif
        ierr=read_nc1dat(ncid, fname, 'theta_b', 1,  tst_val)
        if (ierr.eq.0) then
          if (abs(theta_b-tst_val).gt.epsil) then
            write(*,'(/1x,A,2(A,F9.4,1x),3A/)')       '### ERROR: ',
     &              'Mismatch in ''theta_b'': should be',   theta_b,
     &        'found', tst_val, 'in file ''', fname(1:lfile), '''.'
          endif
        endif
      endif
      check_scoord=icount
      return
      end
      function read_nc1dat(ncid, fname, vname, nlen, value)
      implicit none
      integer(kind=4) read_nc1dat, ncid, nlen
      character(len=*) fname, vname
      real(kind=8) value(nlen)
      integer(kind=4) ndims, dimids(8), varid, vtype, size,
     &               i,j, ierr,  lvar, lfile, lenstr
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
      lvar=lenstr(vname)
      lfile=lenstr(fname)
      ierr=nf_inq_att(ncid, nf_global, vname, vtype, size)
      if (ierr.eq.nf_noerr) then
        if (size.eq.nlen) then
          ierr=nf_get_att_double (ncid,nf_global, vname(1:lvar),value)
          if (ierr.ne.nf_noerr) then
            write(*,'(/1x,6A/12x,A/)')   '### ERROR: read_nc1dat :: ',
     &         'Cannot retrieve global attribute ''',   vname(1:lvar),
     &         ''' from ''', fname(1:lfile), '''.', nf_strerror(ierr)
          endif
        endif
      else
        size=1
        ierr=nf_inq_varid (ncid, vname, varid)
        if (ierr.eq.nf_noerr) then
          ierr=nf_inq_varndims (ncid, varid, ndims)
          if (ierr.eq.nf_noerr) then
            if (ndims.gt.0) then
              ierr=nf_inq_vardimid (ncid, varid, dimids)
              if (ierr.eq.nf_noerr) then
                i=0
                do while (i.lt.ndims .and. ierr.eq.nf_noerr)
                  i=i+1
                  ierr=nf_inq_dimlen (ncid, dimids(i), j)
                  if (ierr.eq.nf_noerr) then
                    size=size*j
                  endif
                enddo
                if (ierr.ne.nf_noerr) then
                  write(*,'(/1x,2A,I3,1x,3A/12x,A/)')  '### ERROR: ',
     &                  'Cannot find size of dimension #', dimids(i),
     &        'in file ''', fname(1:lfile), '''.', nf_strerror(ierr)
                endif
              else
                write(*,'(/1x,6A/12x,A/)')      '### ERROR: Cannot ',
     &               'determine dimension IDs for ''', vname(1:lvar),
     &    ''' in file ''',  fname(1:lfile), '''.', nf_strerror(ierr)
              endif
            endif
          else
            write(*,'(/1x,6A/12x,A/)')          '### ERROR: Cannot ',
     &        'determine number of dimensions for ''', vname(1:lvar),
     &    ''' in file ''',  fname(1:lfile), '''.', nf_strerror(ierr)
          endif
          if (size.eq.nlen .and. ierr.eq.nf_noerr) then
            ierr=nf_get_var_double (ncid, varid, value)
            if (ierr.ne.nf_noerr) then
              write(*,'(/1x,6A/12x,A/)')   '### ERROR: Cannot read ',
     &              'variable ''',  vname(1:lvar), ''' from file ''',
     &                      fname(1:lfile), '''.', nf_strerror(ierr)
            endif
          endif
        else
          write(*,'(/1x,6A/12x,A/)')  '### ERROR: Cannot determine ',
     &           'netCDF ID for ''',  vname(1:lvar), ''' in file ''',
     &                      fname(1:lfile), '''.', nf_strerror(ierr)
        endif
      endif
      if (size.ne.nlen .and. ierr.eq.nf_noerr) then
        write(*,'(/1x,6A/12x,A,I7,1x,A,I7/)')         '### ERROR: ',
     &      'Unexpected size of ''', vname(1:lvar), ''' in file ''',
     &  fname(1:lfile), ''':', 'should be', nlen, 'but found', size
        ierr=+1
      endif
      read_nc1dat=ierr
      return
      end
      integer(kind=4) function my_nf_def_dim (ncid, name, len, dimid)
      implicit none
      integer(kind=4) ncid, len, dimid, ierr
      character(len=*) name
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
      ierr=nf_def_dim (ncid, name, len, dimid)
      if (ierr.ne.nf_noerr) then
        write(*,'(/1x,A,I3,3A,I4,A,I10/12x,A,I3,3A,I3/)')
     &      '#### ERROR: my_nf_def_dim: ncid =', ncid, ' name = ''',
     &       name,  ''' len =', len, ' dimid =', dimid, '  ierr=',
     &       ierr, ' Cause of error: ', nf_strerror(ierr)
      endif
      my_nf_def_dim=ierr
      return
      end
