      subroutine v3dbc_tile (istr,iend,jstr,jend, grad)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, i,j,k
      real(kind=8) grad(istr-2:iend+2,jstr-2:jend+2), cx,cy, cff,
     &         dtfwd, dft,dfx,dfy, tau,tau_in,tau_out
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
      if (nnew.eq.3) then
        dtfwd=0.5_8*dt
      else
        dtfwd=dt
      endif
      tau_in=dtfwd*tauM3_in
      tau_out=dtfwd*tauM3_out
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        do k=1,N
          do i=istr,iend+1
            grad(i,jstr  )=(v(i,jstr  ,k,nstp)-v(i-1,jstr  ,k,nstp))
     &                                                *pmask(i,jstr)
            grad(i,jstr+1)=(v(i,jstr+1,k,nstp)-v(i-1,jstr+1,k,nstp))
     &                                              *pmask(i,jstr+1)
          enddo
          do i=istr,iend
            dft=v(i,jstr+1,k,nstp)-v(i,jstr+1,k,nnew)
            dfx=v(i,jstr+1,k,nnew)-v(i,jstr+2,k,nnew)
            if (dfx*dft .lt. 0._8) then
              dft=0._8
              tau=tau_in
            else
              tau=tau_out
            endif
            if (dft*(grad(i,jstr+1)+grad(i+1,jstr+1)) .gt. 0._8) then
              dfy=grad(i,jstr+1)
            else
              dfy=grad(i+1,jstr+1)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            v(i,jstr,k,nnew)=( cff*v(i,jstr,k,nstp)
     &                        +cx*v(i,jstr+1,k,nnew)
     &                    -max(cy,0._8)*grad(i  ,jstr)
     &                    -min(cy,0._8)*grad(i+1,jstr)
     &                                   )/(cff+cx)
            v(i,jstr,k,nnew)=(1._8-tau)*v(i,jstr,k,nnew)
     &                               +tau*v_south(i,k)
            v(i,jstr,k,nnew)=v(i,jstr,k,nnew)*vmask(i,jstr)
          enddo
        enddo
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        do k=1,N
          do i=istr,iend+1
            grad(i,jend  )=(v(i,jend  ,k,nstp)-v(i-1,jend  ,k,nstp))
     &                                                *pmask(i,jend)
            grad(i,jend+1)=(v(i,jend+1,k,nstp)-v(i-1,jend+1,k,nstp))
     &                                              *pmask(i,jend+1)
          enddo
          do i=istr,iend
            dft=v(i,jend,k,nstp)-v(i,jend  ,k,nnew)
            dfx=v(i,jend,k,nnew)-v(i,jend-1,k,nnew)
            if (dfx*dft .lt. 0._8) then
              dft=0._8
              tau=tau_in
            else
              tau=tau_out
            endif
            if (dft*(grad(i,jend)+grad(i+1,jend)) .gt. 0._8) then
              dfy=grad(i,jend)
            else
              dfy=grad(i+1,jend)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            v(i,jend+1,k,nnew)=( cff*v(i,jend+1,k,nstp)
     &                              +cx*v(i,jend,k,nnew)
     &                      -max(cy,0._8)*grad(i  ,jend+1)
     &                      -min(cy,0._8)*grad(i+1,jend+1)
     &                                      )/(cff+cx)
            v(i,jend+1,k,nnew)=(1._8-tau)*v(i,jend+1,k,nnew)
     &                                   +tau*v_north(i,k)
            v(i,jend+1,k,nnew)=v(i,jend+1,k,nnew)*vmask(i,jend+1)
          enddo
        enddo
      endif
      if (istr.eq.iwest .and. .not.west_exchng) then
        do k=1,N
          do j=jstrV-1,jend
            grad(istr-1,j)=v(istr-1,j+1,k,nstp)-v(istr-1,j,k,nstp)
            grad(istr  ,j)=v(istr  ,j+1,k,nstp)-v(istr  ,j,k,nstp)
          enddo
          do j=jstrV,jend
          cx=-0.125_8*dtfwd*(u(istr,j,k,nrhs)+u(istr,j-1,k,nrhs))
     &                            *( pm(istr-1,j)+pm(istr-1,j-1)
     &                                +pm(istr,j)+pm(istr,j-1) )
          cy= 0.125_8*dtfwd*(v(istr-1,j,k,nrhs)+v(istr,j,k,nrhs))
     &                            *( pn(istr-1,j)+pn(istr-1,j-1)
     &                                +pn(istr,j)+pn(istr,j-1) )
          if (cx.gt.0._8) then
            tau=0._8
          else
            tau=-cx
            cx=0._8
          endif
          v(istr-1,j,k,nnew)=(1._8-cx)*( v(istr-1,j,k,nstp)
     &                          -max(cy,0._8)*grad(istr-1,j-1)
     &                          -min(cy,0._8)*grad(istr-1,j  )
     &                                                     )
     &                       +cx*(        v(istr,j,k,nstp)
     &                            -max(cy,0._8)*grad(istr,j-1)
     &                            -min(cy,0._8)*grad(istr,j  )
     &                                                     )
            v(istr-1,j,k,nnew)=(1._8-tau)*v(istr-1,j,k,nnew)
     &                                    +tau*v_west(j,k)
            v(istr-1,j,k,nnew)=v(istr-1,j,k,nnew)*vmask(istr-1,j)
          enddo
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        do k=1,N
          do j=jstrV-1,jend
            grad(iend  ,j)=v(iend  ,j+1,k,nstp)-v(iend  ,j,k,nstp)
            grad(iend+1,j)=v(iend+1,j+1,k,nstp)-v(iend+1,j,k,nstp)
          enddo
          do j=jstrV,jend
          cx=0.125_8*dtfwd*(u(iend+1,j,k,nrhs)+u(iend+1,j-1,k,nrhs))
     &                               *( pm(iend+1,j)+pm(iend+1,j-1)
     &                                   +pm(iend,j)+pm(iend,j-1) )
          cy=0.125_8*dtfwd*(v(iend,j,k,nrhs)+v(iend+1,j,k,nrhs))
     &                               *( pn(iend+1,j)+pn(iend+1,j-1)
     &                                   +pn(iend,j)+pn(iend,j-1) )
          if (cx.gt.0._8) then
            tau=0._8
          else
            tau=-cx
            cx=0._8
          endif
          v(iend+1,j,k,nnew)=(1._8-cx)*(   v(iend+1,j,k,nstp)
     &                          -max(cy,0._8)*grad(iend+1,j-1)
     &                          -min(cy,0._8)*grad(iend+1,j  )
     &                                                     )
     &                       +cx*(         v(iend,j,k,nstp)
     &                            -max(cy,0._8)*grad(iend,j-1)
     &                            -min(cy,0._8)*grad(iend,j  )
     &                                                     )
            v(iend+1,j,k,nnew)=(1._8-tau)*v(iend+1,j,k,nnew)
     &                                    +tau*v_east(j,k)
            v(iend+1,j,k,nnew)=v(iend+1,j,k,nnew)*vmask(iend+1,j)
          enddo
        enddo
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        do k=1,N
          v(istr-1,jstr,k,nnew)=0.5_8*( v(istr-1,jstr+1,k,nnew)
     &                               +v(istr  ,jstr  ,k,nnew))
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jstr.eq.jsouth 
     &                   .and. .not.south_exchng) then
        do k=1,N
          v(iend+1,jstr,k,nnew)=0.5_8*( v(iend+1,jstr+1,k,nnew)
     &                               +v(iend  ,jstr  ,k,nnew))
        enddo
      endif
      if (istr.eq.iwest .and. .not.west_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        do k=1,N
          v(istr-1,jend+1,k,nnew)=0.5_8*( v(istr-1,jend,k,nnew)
     &                                 +v(istr,jend+1,k,nnew))
        enddo
      endif
      if (iend.eq.ieast .and. .not.east_exchng .and. jend.eq.jnorth 
     &                   .and. .not.north_exchng) then
        do k=1,N
          v(iend+1,jend+1,k,nnew)=0.5_8*( v(iend+1,jend,k,nnew)
     &                                 +v(iend,jend+1,k,nnew))
        enddo
      endif
      return
      end
