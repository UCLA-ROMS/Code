      subroutine setup_kwds (ierr)
      implicit none
      integer(kind=4) ierr, is,ie
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
      do is=1,max_opt_size
        kwds(is:is)=' '
      enddo
      is=1
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='title'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='time_stepping'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='S-coord'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='rho0'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='lateral_visc'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='gamma2'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='tracer_diff2'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bottom_drag'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='nudg_cof'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='v_sponge'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 4
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='grid'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='initial'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bulk_forcing'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='boundary'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='restart'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='history'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +22
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='primary_history_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +24
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='auxiliary_history_fields'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +16
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='primary_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='auxiliary_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bgc_flux_histories'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bgc_flux_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='phys_flux_histories'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='phys_flux_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +20
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bulk_diags_histories'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +19
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='bulk_diags_averages'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      kwds(is:ie)='pCO2_atm_file'
      kwds(ie+1:ie+1)=' '
      is=ie+2
      return
  99  if (mynode.eq.0) write(*,'(/1x,2A/25x,A/)')
     &  '### ERROR: setup_kwds :: Insufficient size of string "kwds" ',
     &  'in file "strings.h".',  'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
