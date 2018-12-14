      subroutine prsgrd (tile)
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
      call prsgrd32AC1_tile (istr,iend,jstr,jend, A3d(1,1),A3d(1,2),
     &                                                      A3d(1,3),
     &                                                      A3d(1,4),
     &                        A2d(1,1), A2d(1,2),
     &                        A2d(1,1), A2d(1,2), A2d(1,3), A2d(1,4))
      return
      end
      subroutine prsgrd32AC1_tile (istr,iend,jstr,jend, ru,rv, P,
     &                                                        rho,
     &                                       dR,dZ, FC,dZx,rx,dRx)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, i,j,k, imin,imax,jmin,jmax
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
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2,N) :: ru,rv, P
     &                                                    , rho
      real(kind=8) dpth
      real(kind=8), dimension(istr-2:iend+2,0:N) :: dR,dZ
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: 
     &                           FC,dZx,rx,dRx
      real(kind=8) grho, HalfGRho, cff, cfr
      real(kind=8), parameter :: OneFifth=0.2_8, OneTwelfth=1._8/12._8, 
     &                             epsil=0._8
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
      real(kind=8) rho1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /eos_rho1/ rho1
      real(kind=8) qp1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /eos_qp1/ qp1
      real(kind=8), parameter :: qp2=0.0000172_8
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
      if (istr.eq.iwest .and. .not.west_exchng) then
        imin=istrU
      else
        imin=istrU-1
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        imax=iend
      else
        imax=iend+1
      endif
        if (jstr.eq.jsouth .and. .not.south_exchng) then
          jmin=jstrV
        else
          jmin=jstrV-1
        endif
        if (jend.eq.jnorth .and. .not.north_exchng) then
          jmax=jend
        else
          jmax=jend+1
        endif
      grho=g/rho0
      HalfGRho=0.5_8*grho
      do j=jstrV-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
            dpth=          -0.5_8*(z_r(i,j,k+1)+z_r(i,j,k))
            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k)
     &              +(qp1(i,j,k+1)-qp1(i,j,k))
     &                     *dpth*(1._8-qp2*dpth)
          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1
          do i=istrU-1,iend
            cff=2._8*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))
            cfr=2._8*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0._8
            endif
            dpth=          -z_r(i,j,k)
            dR(i,k)=dR(i,k)  -qp1(i,j,k)*dZ(i,k)*(1._8-2._8*qp2*dpth)
            rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1._8-qp2*dpth)
          enddo
        enddo
        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + grho*( rho(i,j,N)
     &       +0.5_8*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     &          /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     &                                     *(z_r(i,j,k+1)-z_r(i,j,k))
     &     -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     &                              -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )
     &                -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     &                              -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     &                                                             ))
          enddo
        enddo
      enddo
      do k=N,1,-1
        do j=jstr,jend
          do i=imin,imax
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
     &                              *umask(i,j)
            dpth=-0.5_8*(z_r(i,j,k)+z_r(i-1,j,k))
            rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k)
     &                +(qp1(i,j,k)-qp1(i-1,j,k))
     &                     *dpth*(1._8-qp2*dpth) )
     &                              *umask(i,j)
          enddo
        enddo
        if (istr.eq.iwest .and. .not.west_exchng) then
          do j=jstr,jend
            FC(imin-1,j)=FC(imin,j)
            rx(imin-1,j)=rx(imin,j)
          enddo
        endif
        if (iend.eq.ieast .and. .not.east_exchng) then
          do j=jstr,jend
            FC(imax+1,j)=FC(imax,j)
            rx(imax+1,j)=rx(imax,j)
          enddo
        endif
        do j=jstr,jend
          do i=istrU-1,iend
            cff=2._8*FC(i,j)*FC(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0._8
            endif
            cfr=2._8*rx(i,j)*rx(i+1,j)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0._8
            endif
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &                      *(1._8+2._8*qp2*z_r(i,j,k))
          enddo
          do i=istrU,iend
            ru(i,j,k)=0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &                              P(i-1,j,k)-P(i,j,k)-HalfGRho*(
     &            (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))
     &   -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )
     &              -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     &                                                            )))
          enddo
        enddo
        do j=jmin,jmax
          do i=istr,iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
     &                              *vmask(i,j)
            dpth=-0.5_8*(z_r(i,j,k)+z_r(i,j-1,k))
            rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k)
     &                +(qp1(i,j,k)-qp1(i,j-1,k))
     &                     *dpth*(1._8-qp2*dpth) )
     &                              *vmask(i,j)
          enddo
        enddo
        if (jstr.eq.jsouth .and. .not.south_exchng) then
          do i=istr,iend
            FC(i,jmin-1)=FC(i,jmin)
            rx(i,jmin-1)=rx(i,jmin)
          enddo
        endif
        if (jend.eq.jnorth .and. .not.north_exchng) then
          do i=istr,iend
            FC(i,jmax+1)=FC(i,jmax)
            rx(i,jmax+1)=rx(i,jmax)
          enddo
        endif
        do j=jstrV-1,jend
          do i=istr,iend
            cff=2._8*FC(i,j)*FC(i,j+1)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0._8
            endif
            cfr=2._8*rx(i,j)*rx(i,j+1)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0._8
            endif
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &                         *(1._8+2._8*qp2*z_r(i,j,k))
          enddo
          if (j.ge.jstrV) then
            do i=istr,iend
              rv(i,j,k)=0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &                             P(i,j-1,k)-P(i,j,k) -HalfGRho*(
     &            (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))
     &   -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )
     &              -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     &                                                            )))
            enddo
          endif
        enddo
      enddo
      return
      end
