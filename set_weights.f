      subroutine set_weights
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
      integer(kind=4) i,j, iter
      real(kind=8) p,q,r, scale
      real(kind=8) sum, shft, cff
      nfast=0
      do i=1,2*ndtfast
        weight(1,i)=0._8
        weight(2,i)=0._8
      enddo
      p=2.0_8 ; q=4.0_8 ; r=0.284_8*(1._8-2.8_8/(dble(ndtfast)))
      r=0.25_8
      scale=(p+1._8)*(p+q+1._8)/((p+2._8)*(p+q+2._8)*dble(ndtfast))
      do iter=1,16
        nfast=0
        do i=1,2*ndtfast
          cff=scale*dble(i)
          weight(1,i)=cff**p - cff**(p+q) - r*cff
          if (weight(1,i).gt.0._8) nfast=i
          if (nfast.gt.0 .and. weight(1,i).lt.0._8) weight(1,i)=0._8
        enddo
        sum=0._8
        shft=0._8
        do i=1,nfast
          sum=sum+weight(1,i)
          shft=shft+weight(1,i)*dble(i)
        enddo
        scale=scale * shft/(sum*dble(ndtfast))
      enddo
      do iter=1,ndtfast
        sum=0._8
        shft=0._8
        do i=1,nfast
          sum=sum+weight(1,i)
          shft=shft+dble(i)*weight(1,i)
        enddo
        shft=shft/sum
        cff=dble(ndtfast)-shft
        if (cff .gt. 1._8) then
          nfast=nfast+1
          do i=nfast,2,-1
            weight(1,i)=weight(1,i-1)
          enddo
          weight(1,1)=0._8
        elseif (cff .gt. 0._8) then
          sum=1._8-cff
          do i=nfast,2,-1
            weight(1,i)=sum*weight(1,i)+cff*weight(1,i-1)
          enddo
          weight(1,1)=sum*weight(1,1)
        elseif (cff .lt. -1._8) then
          nfast=nfast-1
          do i=1,nfast,+1
            weight(1,i)=weight(1,i+1)
          enddo
          weight(1,nfast+1)=0._8
        elseif (cff .lt. 0._8) then
          sum=1._8+cff
          do i=1,nfast-1,+1
            weight(1,i)=sum*weight(1,i)-cff*weight(1,i+1)
          enddo
          weight(1,nfast)=sum*weight(1,nfast)
        endif
      enddo
      do j=1,nfast
        cff=weight(1,j)
        do i=1,j
          weight(2,i)=weight(2,i)+cff
        enddo
      enddo
      sum=0._8
      cff=0._8
      do i=1,nfast
        sum=sum+weight(1,i)
        cff=cff+weight(2,i)
      enddo
      sum=1._8/sum
      cff=1._8/cff
      do i=1,nfast
        weight(1,i)=sum*weight(1,i)
        weight(2,i)=cff*weight(2,i)
      enddo
      if (mynode.eq.0) write(*,'(/1x,A,I3,4x,A,I4,8x,A,2F5.1,F9.4/)')
     &        'Mode splitting: ndtfast =', ndtfast, 'nfast =', nfast
     &                                             ,'p,q,r =', p,q,r
      return
      end
