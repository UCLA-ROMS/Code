      subroutine u2dbc (tile)
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
      call u2dbc_tile (istr,iend,jstr,jend, A2d(1,1))
      return
      end
      subroutine u2dbc_tile (istr,iend,jstr,jend, grad)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, i,j
      real(kind=8) grad(istr-2:iend+2,jstr-2:jend+2), cx,cy, cff,
     &         dft,dfx,dfy, tau,tau_in,tau_out, zx,hx
      real(kind=8), parameter :: eps=1.D-33
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
      tau_in=dtfast*tauM2_in
      tau_out=dtfast*tauM2_out
      if (istr.eq.iwest .and. .not.west_exchng) then
        do j=jstr,jend
          cff=0.5_8*(h(istr-1,j)+h(istr,j))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5_8*(pm(istr-1,j)+pm(istr,j))
          zx=(0.5_8+cx)*zeta(istr,j,kstp)+(0.5_8-cx)*zeta(istr-1,j,kstp)
          if (cx .gt. 0.292893218813452_8) then
            zx=zx + ( zeta(istr,j,knew) +cx*zeta(istr-1,j,kstp)
     &                               -(1._8+cx)*zeta(istr,j,kstp)
     &                           )*(1._8-0.292893218813452_8/cx)**2
          endif
          ubar(istr,j,knew)=0.5_8*( (1._8-cx)*ubar(istr,j,kstp)
     &                               +cx*ubar(istr+1,j,kstp)
     &                   +ubar_west(j) -hx*(zx-zeta_west(j))
     &                                                    )
          ubar(istr,j,knew)=ubar(istr,j,knew)*umask(istr,j)
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        do j=jstr,jend
          cff=0.5_8*(h(iend,j)+h(iend+1,j))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5_8*(pm(iend,j)+pm(iend+1,j))
          zx=(0.5_8+cx)*zeta(iend,j,kstp)+(0.5_8-cx)*zeta(iend+1,j,kstp)
          if (cx .gt. 0.292893218813452_8) then
            zx=zx + ( zeta(iend,j,knew) +cx*zeta(iend+1,j,kstp)
     &                               -(1._8+cx)*zeta(iend,j,kstp)
     &                           )*(1._8-0.292893218813452_8/cx)**2
          endif
          ubar(iend+1,j,knew)=0.5_8*( (1._8-cx)*ubar(iend+1,j,kstp)
     &                                     +cx*ubar(iend,j,kstp)
     &                       +ubar_east(j) +hx*(zx-zeta_east(j))
     &                                                        )
          ubar(iend+1,j,knew)=ubar(iend+1,j,knew)*umask(iend+1,j)
        enddo
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        do i=istrU-1,iend
          grad(i,jstr-1)=ubar(i+1,jstr-1,kstp)-ubar(i,jstr-1,kstp)
          grad(i,jstr  )=ubar(i+1,jstr  ,kstp)-ubar(i,jstr  ,kstp)
        enddo
        do i=istrU,iend
          cx=-0.125_8*dtfast*(vbar(i,jstr,kstp)+vbar(i-1,jstr,kstp))
     &      *(pn(i,jstr-1)+pn(i-1,jstr-1)+pn(i,jstr)+pn(i-1,jstr))
          cy= 0.125_8*dtfast*(ubar(i,jstr-1,kstp)+ubar(i,jstr,kstp))
     &      *(pm(i,jstr-1)+pm(i-1,jstr-1)+pm(i,jstr)+pm(i-1,jstr))
          if (cx.gt.0._8) then
            tau=0._8
          else
            tau=-cx
            cx=0._8
          endif
          ubar(i,jstr-1,knew)=(1._8-cx)*( ubar(i,jstr-1,kstp)
     &                          -max(cy,0._8)*grad(i-1,jstr-1)
     &                          -min(cy,0._8)*grad(i  ,jstr-1)
     &                                                     )
     &                       +cx*(        ubar(i,jstr,kstp)
     &                            -max(cy,0._8)*grad(i-1,jstr)
     &                            -min(cy,0._8)*grad(i  ,jstr)
     &                                                     )
          ubar(i,jstr-1,knew)=(1._8-tau)*ubar(i,jstr-1,knew)
     &                                  +tau*ubar_south(i)
          ubar(i,jstr-1,knew)=ubar(i,jstr-1,knew)*umask(i,jstr-1)
        enddo
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        do i=istrU-1,iend
          grad(i,jend  )=ubar(i+1,jend  ,kstp)-ubar(i,jend,kstp  )
          grad(i,jend+1)=ubar(i+1,jend+1,kstp)-ubar(i,jend+1,kstp)
        enddo
        do i=istrU,iend
          cx=0.125_8*dtfast*(vbar(i,jend+1,kstp)+vbar(i-1,jend+1,kstp))
     &         *(pn(i,jend)+pn(i-1,jend)+pn(i,jend+1)+pn(i-1,jend+1))
          cy=0.125_8*dtfast*(ubar(i,jend,kstp)+ubar(i,jend+1,kstp))
     &         *(pm(i,jend)+pm(i-1,jend)+pm(i,jend+1)+pm(i-1,jend+1))
          if (cx.gt.0._8) then
            tau=0._8
          else
            tau=-cx
            cx=0._8
          endif
          ubar(i,jend+1,knew)=(1._8-cx)*( ubar(i,jend+1,kstp)
     &                          -max(cy,0._8)*grad(i-1,jend+1)
     &                          -min(cy,0._8)*grad(i  ,jend+1)
     &                                                     )
     &                       +cx*(        ubar(i,jend,kstp)
     &                            -max(cy,0._8)*grad(i-1,jend)
     &                            -min(cy,0._8)*grad(i  ,jend)
     &                                                     )
          ubar(i,jend+1,knew)=(1._8-tau)*ubar(i,jend+1,knew)
     &                                  +tau*ubar_north(i)
          ubar(i,jend+1,knew)=ubar(i,jend+1,knew)*umask(i,jend+1)
        enddo
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        ubar(istr,jstr-1,knew)=0.5_8*( ubar(istr+1,jstr-1,knew)
     &                                  +ubar(istr,jstr,knew))
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        ubar(iend+1,jstr-1,knew)=0.5_8*( ubar(iend,jstr-1,knew)
     &                                +ubar(iend+1,jstr,knew))
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        ubar(istr,jend+1,knew)=0.5_8*( ubar(istr+1,jend+1,knew)
     &                                  +ubar(istr,jend,knew))
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        ubar(iend+1,jend+1,knew)=0.5_8*( ubar(iend,jend+1,knew)
     &                                +ubar(iend+1,jend,knew))
      endif
      return
      end
