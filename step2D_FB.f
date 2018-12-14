      subroutine step2d (tile)
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
      call step2D_FB_tile (istr,iend,jstr,jend, A2d(1,1),  A2d(1,2),
     &                    A2d(1, 3), A2d(1, 4), A2d(1, 5), A2d(1, 6),
     &                    A2d(1, 7), A2d(1, 8), A2d(1, 9), A2d(1,10),
     &                               A2d(1,11), A2d(1,12), A2d(1,13))
      return
      end
      subroutine step2D_FB_tile (istr,iend,jstr,jend, zeta_new,Dnew,
     &                           rubar,rvbar, urhs,vrhs,  DUon,DVom,
     &                                        Drhs, UFx,UFe,VFx,VFe)
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
      integer(kind=4) istr,iend,jstr,jend, i,j, kbak, kold
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: zeta_new, 
     &                               Dnew,
     &                         rubar,rvbar,  urhs,vrhs,  DUon,DVom,
     &                                       Drhs, UFx,UFe,VFx,VFe
      real(kind=8) cff, cff0,cff1,cff2,cff3,  DUnew,DVnew
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
      if (iif.eq.1) then
        kbak=kstp
        kold=kstp
        cff1=1._8
        cff2=0._8
        cff3=0._8
      elseif (iif.eq.1+1) then
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1=1._8
        cff2=0._8
        cff3=0._8
      else
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak-1
        if (kold.lt.1) kold=4
        cff1=1.781105_8  ; cff2=-1.06221_8 ;   cff3=0.281105_8
      endif
      do j=jstrV-2,jend+1
        do i=istrU-2,iend+1
          Drhs(i,j)=h(i,j) +cff1*zeta(i,j,kstp) +cff2*zeta(i,j,kbak)
     &                                          +cff3*zeta(i,j,kold)
        enddo
      enddo
      do j=jstr-1,jend+1
        do i=istrU-1,iend+1
          urhs(i,j)=cff1*ubar(i,j,kstp) +cff2*ubar(i,j,kbak)
     &                                         +cff3*ubar(i,j,kold)
          DUon(i,j)=0.5_8*(Drhs(i,j)+Drhs(i-1,j))*dn_u(i,j)*urhs(i,j)
        enddo
      enddo
      do j=jstrV-1,jend+1
        do i=istr-1,iend+1
          vrhs(i,j)=cff1*vbar(i,j,kstp) +cff2*vbar(i,j,kbak)
     &                                         +cff3*vbar(i,j,kold)
          DVom(i,j)=0.5_8*(Drhs(i,j)+Drhs(i,j-1))*dm_v(i,j)*vrhs(i,j)
        enddo
      enddo
      if (iif.eq.1) then
        cff0=0._8
        cff1=1._8
        cff2=0._8
        cff3=0._8
      elseif (iif.eq.1+1) then
        cff0= 1.0833333333333_8
        cff1=-0.1666666666666_8
        cff2= 0.0833333333333_8
        cff3=0._8
      else
         cff0=0.614_8   ; cff1=0.285_8  ; cff2=0.0880_8  ; cff3=0.013_8
      endif
      do j=jstrV-1,jend
        do i=istrU-1,iend
          zeta_new(i,j)=zeta(i,j,kstp) + dtfast*pm(i,j)*pn(i,j)
     &            *(DUon(i,j)-DUon(i+1,j)+DVom(i,j)-DVom(i,j+1))
          zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
          Dnew(i,j)=zeta_new(i,j)+h(i,j)
          UFx(i,j)=cff0*zeta_new(i,j) +cff1*zeta(i,j,kstp)
     &             +cff2*zeta(i,j,kbak) +cff3*zeta(i,j,kold)
          UFe(i,j)=(1._8+rhoS(i,j))*UFx(i,j)
          VFe(i,j)=UFe(i,j)*UFx(i,j)
          VFx(i,j)=UFx(i,j)*(rhoS(i,j)-rhoA(i,j))
        enddo
      enddo
      call zetabc_tile (istr,iend,jstr,jend, zeta_new)
      do j=jstrR,jendR
        do i=istrR,iendR
          zeta(i,j,knew)=zeta_new(i,j)
        enddo
      enddo
        cff1=weight(1,iif)
        cff2=weight(2,iif)
        if (iif.eq.1) then
          do j=jstrR,jendR
            do i=istrR,iendR
              Zt_avg1(i,j)=cff1*zeta(i,j,knew)
              DU_avg1(i,j)=0._8
              DV_avg1(i,j)=0._8
              DU_avg2(i,j)=cff2*DUon(i,j)
              DV_avg2(i,j)=cff2*DVom(i,j)
            enddo
          enddo
        else
          do j=jstrR,jendR
            do i=istrR,iendR
              Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
              DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
              DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
            enddo
          enddo
        endif
      cff=0.5_8*g
      do j=jstr,jend
        do i=istr,iend
          rubar(i,j)=cff*dn_u(i,j)*( (h(i-1,j)+h(i,j))*(UFe(i-1,j)
     &                        -UFe(i,j)) +VFe(i-1,j)-VFe(i,j)
     &              +(h(i-1,j)-h(i,j))*( VFx(i-1,j)+VFx(i,j)
     &                        +0.333333333333_8*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(UFx(i-1,j)-UFx(i,j)) )
     &                                                              )
          rvbar(i,j)=cff*dm_v(i,j)*( (h(i,j-1)+h(i,j))*(UFe(i,j-1)
     &                        -UFe(i,j)) +VFe(i,j-1)-VFe(i,j)
     &              +(h(i,j-1)-h(i,j))*( VFx(i,j-1)+VFx(i,j)
     &                        +0.333333333333_8*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(UFx(i,j-1)-UFx(i,j)) )
     &                                                              )
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          rubar(i,j)=rubar(i,j) - 0.5_8*(r_D(i,j)+r_D(i-1,j))
     &                  *dm_u(i,j)*dn_u(i,j)*ubar(i,j,kstp)
          rvbar(i,j)=rvbar(i,j) - 0.5_8*(r_D(i,j)+r_D(i,j-1))
     &                  *dm_v(i,j)*dn_v(i,j)*vbar(i,j,kstp)
        enddo
      enddo
      if (iif.eq.1) then
        do j=jstr,jend
          do i=istr,iend
            rufrc(i,j)=rufrc(i,j)-rubar(i,j)
            rvfrc(i,j)=rvfrc(i,j)-rvbar(i,j)
          enddo
        enddo
        do j=jstrV-1,jend
          do i=istrU-1,iend
            UFx(i,j)=zeta_new(i,j)-zeta(i,j,kstp)
            UFe(i,j)=(1._8+rhoS(i,j))*UFx(i,j)
            VFe(i,j)=UFe(i,j)*(zeta_new(i,j)+zeta(i,j,kstp))
            VFx(i,j)=UFx(i,j)*(rhoS(i,j)-rhoA(i,j))
          enddo
        enddo
        cff=0.5_8*g
        do j=jstr,jend
          do i=istr,iend
            rubar(i,j)=rubar(i,j) +cff*dn_u(i,j)*( (h(i-1,j)+h(i,j))
     &          *(UFe(i-1,j)-UFe(i,j)) +VFe(i-1,j)-VFe(i,j)
     &              +(h(i-1,j)-h(i,j))*( VFx(i-1,j)+VFx(i,j)
     &                        +0.333333333333_8*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(UFx(i-1,j)-UFx(i,j)) )
     &                                                              )
            rvbar(i,j)=rvbar(i,j) +cff*dm_v(i,j)*( (h(i,j-1)+h(i,j))
     &          *(UFe(i,j-1)-UFe(i,j)) +VFe(i,j-1)-VFe(i,j)
     &              +(h(i,j-1)-h(i,j))*( VFx(i,j-1)+VFx(i,j)
     &                        +0.333333333333_8*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(UFx(i,j-1)-UFx(i,j)) )
     &                                                              )
          enddo
        enddo
      endif
      do j=jstrV-1,jend
        do i=istrU-1,iend
          DUon(i,j)=zeta(i,j,kstp)+h(i,j)
        enddo
      enddo
      cff=0.5_8*dtfast
      cff1=0.5_8*weight(1,iif)
      do j=jstr,jend
        do i=istrU,iend
          DUnew=( (DUon(i,j)+DUon(i-1,j))*ubar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
     &                             *(rubar(i,j)+rufrc(i,j))
     &                                                    )
     &                                         *umask(i,j)
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
          DU_avg1(i,j)=DU_avg1(i,j) +cff1*DUnew*dn_u(i,j)
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          DVnew=( (DUon(i,j)+DUon(i,j-1))*vbar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
     &                             *(rvbar(i,j)+rvfrc(i,j))
     &                                                    )
     &                                         *vmask(i,j)
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
          DV_avg1(i,j)=DV_avg1(i,j) +cff1*DVnew*dm_v(i,j)
        enddo
      enddo
      call    u2dbc_tile (istr,iend,jstr,jend, UFx)
      call    v2dbc_tile (istr,iend,jstr,jend, UFx)
      if (istr.eq.iwest .and. .not.west_exchng) then
        do j=jstr-1,jendR
          Dnew(istr-1,j)=h(istr-1,j)+zeta_new(istr-1,j)
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        do j=jstr-1,jendR
          Dnew(iend+1,j)=h(iend+1,j)+zeta_new(iend+1,j)
        enddo
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        do i=istr-1,iendR
          Dnew(i,jstr-1)=h(i,jstr-1)+zeta_new(i,jstr-1)
        enddo
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        do i=istr-1,iendR
          Dnew(i,jend+1)=h(i,jend+1)+zeta_new(i,jend+1)
        enddo
      endif
      cff1=0.5_8*weight(1,iif)
      if (istr.eq.iwest .and. .not.west_exchng) then
        do j=jstrR,jendR
          DU_avg1(istrU-1,j)=DU_avg1(istrU-1,j)+cff1*(Dnew(istrU-1,j)
     &         +Dnew(istrU-2,j))*ubar(istrU-1,j,knew)*dn_u(istrU-1,j)
        enddo
        do j=jstrV,jend
          DV_avg1(istr-1,j)=DV_avg1(istr-1,j) +cff1*(Dnew(istr-1,j)
     &       +Dnew(istr-1,j-1) )*vbar(istr-1,j,knew)*dm_v(istr-1,j)
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        do j=jstrR,jendR
          DU_avg1(iend+1,j)=DU_avg1(iend+1,j) +cff1*( Dnew(iend+1,j)
     &            +Dnew(iend,j) )*ubar(iend+1,j,knew)*dn_u(iend+1,j)
        enddo
        do j=jstrV,jend
          DV_avg1(iend+1,j)=DV_avg1(iend+1,j) +cff1*( Dnew(iend+1,j)
     &        +Dnew(iend+1,j-1) )*vbar(iend+1,j,knew)*dm_v(iend+1,j)
        enddo
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        do i=istrU,iend
          DU_avg1(i,jstr-1)=DU_avg1(i,jstr-1) +cff1*( Dnew(i,jstr-1)
     &        +Dnew(i-1,jstr-1) )*ubar(i,jstr-1,knew)*dn_u(i,jstr-1)
        enddo
        do i=istrR,iendR
          DV_avg1(i,jstrV-1)=DV_avg1(i,jstrV-1)+cff1*(Dnew(i,jstrV-1)
     &         +Dnew(i,jstrV-2))*vbar(i,jstrV-1,knew)*dm_v(i,jstrV-1)
        enddo
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        do i=istrU,iend
          DU_avg1(i,jend+1)=DU_avg1(i,jend+1) +cff1*( Dnew(i,jend+1)
     &        +Dnew(i-1,jend+1) )*ubar(i,jend+1,knew)*dn_u(i,jend+1)
        enddo
        do i=istrR,iendR
          DV_avg1(i,jend+1)=DV_avg1(i,jend+1) +cff1*( Dnew(i,jend+1)
     &            +Dnew(i,jend) )*vbar(i,jend+1,knew)*dm_v(i,jend+1)
        enddo
      endif
      if (iif.eq.nfast) then
        do j=jstrR,jendR
          do i=istrR,iendR
            zeta(i,j,knew)=Zt_avg1(i,j)
          enddo
        enddo
        call set_depth_tile (istr,iend,jstr,jend, UFx)
      endif
      call exchange2d_3_tile (istr,iend,jstr,jend,
     &                   zeta(-1,-1,knew),
     &                   ubar(-1,-1,knew),
     &                   vbar(-1,-1,knew))
      return
      end
