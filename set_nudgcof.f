      subroutine set_nudgcof (tile)
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
      call set_nudgcof_tile (istr,iend,jstr,jend, A2d(1,1))
      return
      end
      subroutine set_nudgcof_tile (istr,iend,jstr,jend, wrk)
      implicit none
      integer(kind=4) istr,iend,jstr,jend
      real(kind=8) wrk(istr-2:iend+2,jstr-2:jend+2)
      integer(kind=4) i,j, isp , I1,I0,I2,J1,J0,J2, ibnd,itrc
      real(kind=8) cff,cff1,cff2,cff3, frac, Tau
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
      integer(kind=4) istrR,iendR,jstrR,jendR
      if (istr.eq.iwest) then
        if (west_exchng.and.istr.eq.iwest) then
          istrR=istr-2
        else
          istrR=istr-1
        endif
      else
        istrR=istr
      endif
      if (iend.eq.ieast) then
        if (east_exchng.and.iend.eq.ieast) then
          iendR=iend+2
        else
          iendR=iend+1
        endif
      else
        iendR=iend
      endif
      if (jstr.eq.jsouth) then
        if (south_exchng.and.jstr.eq.jsouth) then
          jstrR=jstr-2
        else
          jstrR=jstr-1
        endif
      else
        jstrR=jstr
      endif
      if (jend.eq.jnorth) then
        if (north_exchng.and.jend.eq.jnorth) then
          jendR=jend+2
        else
          jendR=jend+1
        endif
      else
        jendR=jend
      endif
      isp=(LLm+1)/12
      do j=max(-1,jstrR-1),jendR
        do i=max(-1,istrR-1),iendR
          ibnd=isp
          ibnd=min(ibnd,i + iSW_corn)
          ibnd=min(ibnd,j + jSW_corn)
          ibnd=min(ibnd,MMm+1 - j-jSW_corn)
          wrk(i,j)=dble(isp-ibnd)/dble(isp)
        enddo
      enddo
      do j=jstrR,jendR
        do i=istrR,iendR
          visc2_r(i,j)=visc2_r(i,j)+v_sponge*wrk(i,j)
          diff2(i,j,itemp)=diff2(i,j,itemp)+v_sponge*wrk(i,j)
        enddo
      enddo
      do j=jstr,jendR
        do i=istr,iendR
          visc2_p(i,j)=visc2_p(i,j)+0.25_8*v_sponge*( wrk(i,j)
     &                  +wrk(i-1,j)+wrk(i,j-1)+wrk(i-1,j-1))
        enddo
      enddo
      do itrc=2,NT
        do j=jstrR,jendR
          do i=istrR,iendR
            diff2(i,j,itrc)=diff2(i,j,itemp)
          enddo
        enddo
      enddo
        call exchange2d_tile (istr,iend,jstr,jend, visc2_p)
      return
      end
