      subroutine setup_grid2 (tile)
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
      call setup_grid2_tile (istr,iend,jstr,jend, A2d(1,1),A2d(1,3))
      return
      end
      subroutine setup_grid2_tile (istr,iend,jstr,jend, dA,dV)
      implicit none
      integer(kind=4) istr,iend,jstr,jend, i,j, nsub
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: dA,dV
      real(kind=8)  my_area,  my_volume,  my_crss
      real(kind=8) my_hmax, my_grdmax, my_Cu_max, cff,
     &     my_hmin, my_grdmin, my_Cu_min, my_Cu_Cor
      integer(kind=4) is,isize,itg, js,jsize,jtg
     &                    , jnc
     &                    , inc
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
      real(kind=8) area, volume, bc_crss
      common /comm_vars/ area, volume, bc_crss
      real(kind=8) global_sum(0:2*NT+1), global_srf_sum(0:NT)
      common /comm_vars_bio/ global_sum, global_srf_sum
      real(kind=8) hmin,hmax, grdmin,grdmax, rx0,rx1, Cg_min,Cg_max, 
     &                               Cu_Cor
      common /comm_vars/ hmin,hmax, grdmin,grdmax, rx0,rx1,
     &                                        Cg_min,Cg_max, Cu_Cor
      real(kind=4) cpu_all(4)
      integer(kind=4) trd_count, tile_count, bc_count, mcheck, 
     &                             first_time
      common /comm_vars/ cpu_all, trd_count,
     &                   tile_count, bc_count, mcheck, first_time
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
      include "mpif.h"
      integer(kind=4) size, step, itag, status(MPI_STATUS_SIZE), ierr
      real(kind=8) buff(16)
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
      my_hmin=+1.D+20
      my_hmax=-1.D+20
      my_grdmin=+1.D+20
      my_grdmax=-1.D+20
      my_Cu_min=+1.D+20
      my_Cu_max=-1.D+20
      my_Cu_Cor=-1.D+20
      do j=jstrR,jendR
        do i=istrR,iendR
          if (rmask(i,j) .gt. 0.5_8) then
            my_hmin=min(my_hmin, h(i,j))
            my_hmax=max(my_hmax, h(i,j))
            cff=1._8/sqrt(pm(i,j)*pn(i,j))
            my_grdmin=min(my_grdmin, cff)
            my_grdmax=max(my_grdmax, cff)
            cff=dtfast*sqrt(g*h(i,j)*(pm(i,j)*pm(i,j)+pn(i,j)*pn(i,j)))
            my_Cu_min=min(my_Cu_min, cff)
            my_Cu_max=max(my_Cu_max, cff)
            cff=dt*abs(f(i,j))
            my_Cu_Cor=max(my_Cu_Cor, cff)
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          dA(i,j)=rmask(i,j)/(pm(i,j)*pn(i,j))
          dV(i,j)=dA(i,j)*h(i,j)
        enddo
      enddo
      isize=iend-istr  ; jsize=jend-jstr
      do while (isize>0 .or. jsize>0)
        if (jsize>0) then
          js=(jsize+1)/2-1
          do j=0,js
            jtg=jstr+j
            do i=istr,istr+isize
              dA(i,jtg)=dA(i,jtg+j)+dA(i,jtg+j+1)
              dV(i,jtg)=dV(i,jtg+j)+dV(i,jtg+j+1)
            enddo
          enddo
          if (2*js+1 < jsize) then
            js=js+1
            jtg=jstr+js
            do i=istr,istr+isize
              dA(i,jtg)=dA(i,jtg+js)
              dV(i,jtg)=dV(i,jtg+js)
            enddo
          endif
          jsize=js
        endif
        if (isize>0) then
          is=(isize+1)/2-1
          do j=jstr,jstr+jsize
            do i=0,is
              itg=istr+i
              dA(itg,j)=dA(itg+i,j)+dA(itg+i+1,j)
              dV(itg,j)=dV(itg+i,j)+dV(itg+i+1,j)
            enddo
          enddo
          if (2*is+1 < isize) then
            is=is+1
            itg=istr+is
            do j=jstr,jstr+jsize
              dA(itg,j)=dA(itg+is,j)
              dV(itg,j)=dV(itg+is,j)
            enddo
          endif
          isize=is
        endif
      enddo
      my_area  =dA(istr,jstr)
      my_volume=dV(istr,jstr)
      my_crss=0.D0
      if (istr.eq.iwest .and. .not.west_exchng) then
        do j=jstr,jend
          dA(istr,j)=0.5_8*(h(istr-1,j)+h(istr,j))*dn_u(istr,j)
     &                                          *umask(istr,j)
        enddo
        jnc=1
        do while(jstr.le.jend-jnc)
          js=2*jnc
          do j=jstr,jend-jnc,js
            dA(istr,j) = dA(istr,j) + dA(istr,j+jnc)
          enddo
          jnc=js
        enddo
        my_crss=my_crss + dA(istr,jstr)
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        do j=jstr,jend
          dA(iend,j)=0.5_8*(h(iend,j)+h(iend+1,j))*dn_u(iend+1,j)
     &                                         *umask(iend+1,j)
        enddo
        jnc=1
        do while(jstr.le.jend-jnc)
          js=2*jnc
          do j=jstr,jend-jnc,js
            dA(iend,j) = dA(iend,j) + dA(iend,j+jnc)
          enddo
          jnc=js
        enddo
        my_crss=my_crss + dA(iend,jstr)
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        do i=istr,iend
          dA(i,jstr)=0.5_8*(h(i,jstr)+h(i,jstr-1))*dm_v(i,jstr)
     &                                         *vmask(i,jstr)
        enddo
        inc=1
        do while(istr.le.iend-inc)
          is=2*inc
          do i=istr,iend-inc,is
            dA(i,jstr) = dA(i,jstr) + dA(i+inc,jstr)
          enddo
          inc=is
        enddo
        my_crss=my_crss + dA(istr,jstr)
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        do i=istr,iend
          dA(i,jend)=0.5_8*(h(i,jend)+h(i,jend+1))*dm_v(i,jend+1)
     &                                         *vmask(i,jend+1)
        enddo
        inc=1
        do while(istr.le.iend-inc)
          is=2*inc
          do i=istr,iend-inc,is
            dA(i,jend) = dA(i,jend) + dA(i+inc,jend)
          enddo
          inc=is
        enddo
        my_crss=my_crss + dA(istr,jend)
      endif
      if ((iend-istr.eq.ieast-iwest .and.  jend-jstr.eq.jnorth-jsouth)) 
     &                                then
        nsub=1
      else
        nsub=NSUB_X*NSUB_E
      endif
C$OMP CRITICAL (grd2_cr_rgn)
        if (tile_count.eq.0) then
          hmin=my_hmin
          hmax=my_hmax
          grdmin=my_grdmin
          grdmax=my_grdmax
          Cg_min=my_Cu_min
          Cg_max=my_Cu_max
          Cu_Cor=my_Cu_Cor
          area=my_area
          volume=my_volume
          bc_crss=my_crss
        else
          hmin=min(hmin, my_hmin)
          hmax=max(hmax, my_hmax)
          grdmin=min(grdmin, my_grdmin)
          grdmax=max(grdmax, my_grdmax)
          Cg_min=min(Cg_min, my_Cu_min)
          Cg_max=max(Cg_max, my_Cu_max)
          Cu_Cor=max(Cu_Cor, my_Cu_Cor)
          area=area+my_area
          volume=volume+my_volume
          bc_crss=bc_crss+my_crss
        endif
        tile_count=tile_count+1
        if (tile_count.eq.nsub) then
          tile_count=0
          size=NNODES
          do while (size.gt.1)
           step=(size+1)/2
            if (mynode.ge.step .and. mynode.lt.size) then
              buff(1)=hmin
              buff(2)=hmax
              buff(3)=grdmin
              buff(4)=grdmax
              buff(5)=Cg_min
              buff(6)=Cg_max
              buff(7)=Cu_Cor
              buff(8)=area
              buff(9)=volume
              buff(10)=bc_crss
              itag=mynode+300
              call MPI_Send (buff, 10, MPI_REAL8, mynode-step,
     &                       itag, ocean_grid_comm,          ierr)
            elseif (mynode .lt. size-step) then
              itag=mynode+step+300
              call MPI_Recv (buff, 10, MPI_REAL8, mynode+step,
     &                       itag, ocean_grid_comm,  status, ierr)
              cff=buff(1)
              hmin=  min(hmin,   cff)
              cff=buff(2)
              hmax=  max(hmax,   cff)
              cff=buff(3)
              grdmin=min(grdmin, cff)
              cff=buff(4)
              grdmax=max(grdmax, cff)
              cff=buff(5)
              Cg_min=min(Cg_min, cff)
              cff=buff(6)
              Cg_max=max(Cg_max, cff)
              cff= buff(7)
              Cu_Cor=max(Cu_Cor, cff)
              area=area + buff(8)
              volume=volume + buff(9)
              bc_crss=bc_crss + buff(10)
            endif
           size=step
          enddo
          buff(1)=hmin
          buff(2)=hmax
          buff(3)=grdmin
          buff(4)=grdmax
          buff(5)=Cg_min
          buff(6)=Cg_max
          buff(7)=Cu_Cor
          buff(8)=area
          buff(9)=volume
          buff(10)=bc_crss
          call MPI_Bcast(buff, 10, MPI_REAL8, 0,
     &                       ocean_grid_comm, ierr)
          hmin=  buff(1)
          hmax=  buff(2)
          grdmin=buff(3)
          grdmax=buff(4)
          Cg_min=buff(5)
          Cg_max=buff(6)
          Cu_Cor=buff(7)
          area=buff(8)
          volume=buff(9)
          bc_crss=buff(10)
          if (mynode.eq.0) then
            write(*,'(1x,A,F12.6,3x,A,ES14.7,5x,A,ES23.16)')
     &       'hmin =' ,hmin,  'grdmin =', grdmin, 'area =',   area
            write(*,'(1x,A,F12.6,3x,A,ES14.7,3x,A,ES23.16)')
     &       'hmax =' ,hmax,  'grdmax =', grdmax, 'volume =', volume
            write(*,'(43x,A,ES23.16)')        'open_cross =', bc_crss
            write(*,'(1x,A,F10.7,3x,A,F10.7,3x,A,F10.7)')
     &       'Cg_max =',Cg_max, 'Cg_min =',Cg_min, 'Cu_Cor =',Cu_Cor
          endif
        endif
C$OMP END CRITICAL (grd2_cr_rgn)
      return
      end
