      subroutine step3d_uv1 (tile)
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
      call step3d_uv1_tile (istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                                             A2d(1,5),A2d(1,6))
      return
      end
      subroutine step3d_uv1_tile (istr,iend,jstr,jend, ru,rv,
     &                                            FC,WC,CF,DC,
     &                             UFx,UFe,VFx,VFe, wrk1,wrk2)
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
      integer(kind=4) istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2,N) :: ru,rv
      real(kind=8), dimension(istr-2:iend+2,0:N) :: FC,WC,CF,DC
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: 
     &                          UFx,UFe,VFx,VFe,
     &                                                   wrk1,wrk2
      real(kind=8) cff, stress
      real(kind=8), parameter ::  delta=0.125
     &                  , gamma=0.25_8
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
      real(kind=8) sustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) svstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      integer(kind=4) ntrad,nttra,ntprec,ntwnd
      common /blk_nt/ntrad,nttra,ntprec,ntwnd
      integer(kind=4) tra_file_id, rad_file_id, prec_file_id,
     &        wnd_file_id
      common /blk_id/ tra_file_id,rad_file_id,
     &                prec_file_id,wnd_file_id
      real(kind=8) stflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_stflx/stflx
      real(kind=8) ust(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulk_ust/ust
      real(kind=8) tst(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulk_tst/tst
      real(kind=8) qst(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulk_qst/qst
      real(kind=8) prec_scale,tra_scale,wnd_scale,srf_scale
      common /blk_scale/prec_scale,tra_scale,wnd_scale,srf_scale
      real(kind=8) tair(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) qair(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) rain(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) radlw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) radsw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) uwnd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) vwnd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /bulk_tair/tair
     &       /bulk_qair/qair
     &       /bulk_rain/rain
     &       /bulk_radlw/radlw
     &       /bulk_radsw/radsw
     &       /bulk_uwnd/uwnd
     &       /bulk_vwnd/vwnd
      real(kind=8) tairg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) qairg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) raing(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) radlwg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) radswg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) uwndg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real(kind=8) vwndg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulkdat_tairg/tairg
     &       /bulkdat_qairg/qairg
     &       /bulkdat_raing/raing
     &       /bulkdat_radlwg/radlwg
     &       /bulkdat_radswg/radswg
     &       /bulk_uwndg/uwndg
     &       /bulk_vwndg/vwndg
      real(kind=8)    rad_time(2), rad_cycle
      real(kind=8)    tra_time(2), tra_cycle
      real(kind=8)    wnd_time(2), wnd_cycle
      real(kind=8)    prec_time(2),prec_cycle
      real(kind=8)    ztref, zref
      common /bulkdat2/
     &        rad_time, rad_cycle,
     &        tra_time, tra_cycle,
     &        wnd_time, wnd_cycle,
     &        prec_time, prec_cycle
     &        ,ztref, zref
      integer(kind=4) tair_id,uwnd_id,radsw_id,prec_id
      integer(kind=4) qair_id,vwnd_id,radlw_id
      integer(kind=4) itrad,rad_ncycle,rad_rec,rad_tid
      integer(kind=4) ittra,tra_ncycle,tra_rec,tra_tid
      integer(kind=4) itwnd,wnd_ncycle,wnd_rec,wnd_tid
      integer(kind=4) itprec,prec_ncycle,prec_rec,prec_tid
      common /bulkdat1/
     &        tair_id,uwnd_id,radsw_id,prec_id,
     &        qair_id,vwnd_id,radlw_id,
     &        itrad,rad_ncycle,rad_rec,rad_tid,
     &        ittra,tra_ncycle,tra_rec,tra_tid,
     &        itwnd,wnd_ncycle,wnd_rec,wnd_tid,
     &        itprec,prec_ncycle,prec_rec,prec_tid
      real(kind=8) sssg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /sss_dat/sssg
      real(kind=8) sss_cycle, sss_time(2)
      integer(kind=4) sss_ncycle,   sss_rec,  itsss,    ntsss,
     &        sss_file_id,  sss_id,   sss_tid
      common /sssrest_data/ sss_cycle,sss_time,sss_ncycle,
     &  sss_rec,itsss,ntsss,sss_file_id,sss_id, sss_tid
      real(kind=8), parameter :: coeff_ssuv =
     &     0.77_8
      real(kind=8) sustr_blk(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_sustr_blk/sustr_blk
      real(kind=8) svstr_blk(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_svstr_blk/svstr_blk
      real(kind=8) shflx_net(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_net/shflx_net
      real(kind=8) shflx_lat(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_lat/shflx_lat
      real(kind=8) shflx_sen(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_sen/shflx_sen
      real(kind=8) shflx_rad(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_rad/shflx_rad
      real(kind=8) swflx_emp(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_swflx_emp/swflx_emp
      real(kind=8) shflx_wwk(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_wwk/shflx_wwk
      real(kind=8) surf_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_surf_u/surf_u
      real(kind=8) surf_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_surf_v/surf_v
      real(kind=8) time_bulk_diags_his
      common /t_bulk_diags_his/time_bulk_diags_his
      real(kind=8) zeta_bulk_diags_his(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /zeta_bulk_diags_his/zeta_bulk_diags_his
      logical new_bulk_diags_his
      integer(kind=4) n_bulk_diags_his
      common /nc_bulk_diags_his/n_bulk_diags_his,new_bulk_diags_his
      integer(kind=4) nrpf_bulk_diags_his,ncid_bulk_diags_his
     &       , nrec_bulk_diags_his,bulk_diags_hisTstep, bulk_diags_hisZ
     &       , bulk_diags_hisTime
      common /nc_bulk_diags_his/ nrpf_bulk_diags_his
     &     , ncid_bulk_diags_his, nrec_bulk_diags_his
     &     , bulk_diags_hisTstep, bulk_diags_hisZ
     &     , bulk_diags_hisTime
      character(len=80) bulk_diags_his_name
      common /nc_bulk_diags_his/bulk_diags_his_name
      real(kind=8) time_bulk_diags_avg
      common /t_bulk_diags_avg/time_bulk_diags_avg
      real(kind=8) zeta_bulk_diags_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /zeta_bulk_diags_avg/zeta_bulk_diags_avg
      logical new_bulk_diags_avg
      integer(kind=4) nts_bulk_diags_avg, n_bulk_diags_avg
      common /nc_bulk_diags_avg/nts_bulk_diags_avg
     $        ,n_bulk_diags_avg,new_bulk_diags_avg
      integer(kind=4) nrpf_bulk_diags_avg,ncid_bulk_diags_avg
     &       , nrec_bulk_diags_avg,bulk_diags_avgTstep, bulk_diags_avgZ
     &       , bulk_diags_avgTime
      common /nc_bulk_diags_avg/ nrpf_bulk_diags_avg
     &     , ncid_bulk_diags_avg, nrec_bulk_diags_avg
     &     , bulk_diags_avgTstep, bulk_diags_avgZ
     &     , bulk_diags_avgTime
      character(len=80) bulk_diags_avg_name
      common /nc_bulk_diags_avg/bulk_diags_avg_name
         real(kind=8) dust(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
           common /forces_dust/dust
       real(kind=8) dustg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
       common /dustdat_dustg/dustg
        real(kind=8) dustp(2), dust_time(2),dust_cycle, scldqdt
        integer(kind=4) itdust,dust_id,ldustgrd ,dust_ncycle,
     &  dust_rec,dust_tid,dust_file_id,iron_file_id
       common/dustdat/itdust,dust_id,ldustgrd,
     &  dust_ncycle,dust_rec,dust_tid,dust_file_id,iron_file_id
       common/dustdat1/dustp,dust_time,dust_cycle,scldqdt
       real(kind=8) iron(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       common /forces_iron/iron
       real(kind=8) irong(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
       common /irondat_irong/irong
       real(kind=8) ironp(2),iron_time(2),iron_cycle
       integer(kind=4) itiron,iron_id,lirongrd,iron_ncycle,
     &  iron_rec,iron_tid
       common/irondat/ironp,iron_time,iron_cycle
       common/irondat1/itiron,iron_id,lirongrd,
     &  iron_ncycle,iron_rec,iron_tid
      real(kind=8) srflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real(kind=8) srflx_dailyavg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer(kind=4), parameter :: n_srflx_day = 24
      integer(kind=4) :: num_srflx_day
      integer(kind=4), parameter :: MAX_NUM_SRFLX_DAY = 144
      integer(kind=4) iptr_srflx_day(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer(kind=4) iptr_srflx_day_set(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer(kind=4) PARinc_rst_read
      real(kind=8) srflx_day(-1:Lm+2+padd_X,-1:Mm+2+padd_E,
     &                         MAX_NUM_SRFLX_DAY)
      real(kind=8) srflx_frac
      common /forces_srflx/srflx
     &       , srflx_dailyavg
     &     , srflx_day, srflx_frac
      common /i_forces_srflx/ iptr_srflx_day, PARinc_rst_read,
     &     iptr_srflx_day_set,  num_srflx_day
      character(len=100) pco2_atm_file
      common /pco2_atm_file/ pco2_atm_file
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
       do k=1,N
        do j=jstrV-1,jend
          do i=istrU-1,iend
            cff=0.5_8*Hz(i,j,k)*(
     &              fomn(i,j)
     &             +0.5_8*( (v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*dndx(i,j)
     &                   -(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*dmde(i,j))
     &                                                             )
            UFx(i,j)=cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
          enddo
        enddo
        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)+0.5_8*(UFx(i,j)+UFx(i-1,j))
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=rv(i,j,k)-0.5_8*(VFe(i,j)+VFe(i,j-1))
          enddo
        enddo
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
        do j=jstr,jend
          do i=imin,imax
            wrk1(i,j)=u(i-1,j,k,nrhs)-2._8*u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
            wrk2(i,j)=FlxU(i-1,j,k) -2._8*FlxU(i,j,k) +FlxU(i+1,j,k)
          enddo
        enddo
        if (istr.eq.iwest .and. .not.west_exchng) then
          do j=jstr,jend
            wrk1(istrU-1,j) =wrk1(istrU,j)
            wrk2(istrU-1,j)=wrk2(istrU,j)
          enddo
        endif
        if (iend.eq.ieast .and. .not.east_exchng) then
          do j=jstr,jend
            wrk1(iend+1,j) =wrk1(iend,j)
            wrk2(iend+1,j)=wrk2(iend,j)
          enddo
        endif
        do j=jstr,jend
          do i=istrU-1,iend
            cff=FlxU(i,j,k)+FlxU(i+1,j,k)-delta*( wrk2(i  ,j)
     &                                            +wrk2(i+1,j))
            UFx(i,j)=0.25_8*( cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
     &                          -gamma*( max(cff,0._8)*wrk1(i  ,j)
     &                                  +min(cff,0._8)*wrk1(i+1,j)
     &                                                      ))
          enddo
        enddo
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
        do j=jmin,jmax
          do i=istr,iend
            wrk1(i,j)=v(i,j-1,k,nrhs)-2._8*v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
            wrk2(i,j)=FlxV(i,j-1,k) -2._8*FlxV(i,j,k)  +FlxV(i,j+1,k)
          enddo
        enddo
        if (jstr.eq.jsouth .and. .not.south_exchng) then
          do i=istr,iend
            wrk1(i,jstrV-1)=wrk1(i,jstrV)
            wrk2(i,jstrV-1)=wrk2(i,jstrV)
          enddo
        endif
        if (jend.eq.jnorth .and. .not.north_exchng) then
          do i=istr,iend
            wrk1(i,jend+1)=wrk1(i,jend)
            wrk2(i,jend+1)=wrk2(i,jend)
          enddo
        endif
        do j=jstrV-1,jend
          do i=istr,iend
            cff=FlxV(i,j,k)+FlxV(i,j+1,k)-delta*( wrk2(i,j  )
     &                                           +wrk2(i,j+1))
            VFe(i,j)=0.25_8*( cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
     &                          -gamma*( max(cff,0._8)*wrk1(i,j  )
     &                                  +min(cff,0._8)*wrk1(i,j+1)
     &                                                      ))
          enddo
        enddo
        if (jstr.eq.jsouth .and. .not.south_exchng) then
          jmin=jstr
        else
          jmin=jstr-1
        endif
        if (jend.eq.jnorth .and. .not.north_exchng) then
          jmax=jend
        else
          jmax=jend+1
        endif
        do j=jmin,jmax
          do i=istrU,iend
            wrk1(i,j)=u(i,j-1,k,nrhs)-2._8*u(i,j,k,nrhs)+u(i,j+1,k,nrhs)
          enddo
        enddo
        if (jstr.eq.jsouth .and. .not.south_exchng) then
          do i=istrU,iend
            wrk1(i,jstr-1)=wrk1(i,jstr)
          enddo
        endif
        if (jend.eq.jnorth .and. .not.north_exchng) then
          do i=istrU,iend
            wrk1(i,jend+1)=wrk1(i,jend)
          enddo
        endif
        do j=jstr,jend+1
          do i=istrU-1,iend
           wrk2(i,j)=FlxV(i-1,j,k)-2._8*FlxV(i,j,k)+FlxV(i+1,j,k)
          enddo
        enddo
        do j=jstr,jend+1
          do i=istrU,iend
            cff=FlxV(i,j,k)+FlxV(i-1,j,k)-delta*( wrk2(i  ,j)
     &                                           +wrk2(i-1,j))
            UFe(i,j)=0.25_8*( cff*(u(i,j,k,nrhs)+u(i,j-1,k,nrhs))
     &                          -gamma*( max(cff,0._8)*wrk1(i,j-1)
     &                                  +min(cff,0._8)*wrk1(i,j  )
     &                                                      ))
          enddo
        enddo
        if (istr.eq.iwest .and. .not.west_exchng) then
          imin=istr
        else
          imin=istr-1
        endif
        if (iend.eq.ieast .and. .not.east_exchng) then
          imax=iend
        else
          imax=iend+1
        endif
        do j=jstrV,jend
          do i=imin,imax
            wrk1(i,j)=v(i-1,j,k,nrhs)-2._8*v(i,j,k,nrhs)+v(i+1,j,k,nrhs)
          enddo
        enddo
        if (istr.eq.iwest .and. .not.west_exchng) then
          do j=jstrV,jend
            wrk1(istr-1,j)=wrk1(istr,j)
          enddo
        endif
        if (iend.eq.ieast .and. .not.east_exchng) then
          do j=jstrV,jend
            wrk1(iend+1,j)=wrk1(iend,j)
          enddo
        endif
        do j=jstrV-1,jend
          do i=istr,iend+1
           wrk2(i,j)=FlxU(i,j-1,k)-2._8*FlxU(i,j,k)+FlxU(i,j+1,k)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend+1
            cff=FlxU(i,j,k)+FlxU(i,j-1,k)-delta*( wrk2(i,j  )
     &                                           +wrk2(i,j-1))
            VFx(i,j)=0.25_8*( cff*(v(i,j,k,nrhs)+v(i-1,j,k,nrhs))
     &                          -gamma*( max(cff,0._8)*wrk1(i-1,j)
     &                                  +min(cff,0._8)*wrk1(i  ,j)
     &                                                      ))
          enddo
        enddo
        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)-UFx(i,j  )+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i  ,j)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=rv(i,j,k)-VFx(i+1,j)+VFx(i,j  )
     &                         -VFe(i  ,j)+VFe(i,j-1)
          enddo
        enddo
       enddo
      do j=jstr,jend
        do i=istrU,iend
          DC(i,1)=0.5625_8*(Hz(i  ,j,1)+Hz(i-1,j,1))
     &           -0.0625_8*(Hz(i+1,j,1)+Hz(i-2,j,1))
          FC(i,0)=1.5_8*u(i,j,1,nrhs)
          CF(i,1)=0.5_8
        enddo
        do k=1,N-1,+1
          do i=istrU,iend
            DC(i,k+1)=0.5625_8*(Hz(i  ,j,k+1)+Hz(i-1,j,k+1))
     &               -0.0625_8*(Hz(i+1,j,k+1)+Hz(i-2,j,k+1))
            cff=1._8/(2._8*DC(i,k)+DC(i,k+1)*(2._8-CF(i,k)))
            CF(i,k+1)=cff*DC(i,k)
            FC(i,k)=cff*( 3._8*( DC(i,k  )*u(i,j,k+1,nrhs)
     &                        +DC(i,k+1)*u(i,j,k  ,nrhs))
     &                              -DC(i,k+1)*FC(i,k-1))
          enddo
        enddo
        do i=istrU,iend
          FC(i,N)=(3._8*u(i,j,N,nrhs)-FC(i,N-1))/(2._8-CF(i,N))
          DC(i,N)=0._8
        enddo
        do k=N-1,1,-1
          do i=istrU,iend
            FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
            DC(i,k)=FC(i,k) * 0.5_8*( We(i,j,k)+We(i-1,j,k) -0.125_8*(
     &                        (We(i+1,j,k)-We(i  ,j,k))*umask(i+1,j)
     &                       -(We(i-1,j,k)-We(i-2,j,k))*umask(i-1,j)
     &                                                          ))
            ru(i,j,k+1)=ru(i,j,k+1) -DC(i,k+1)+DC(i,k)
          enddo
        enddo
        do i=istrU,iend
          ru(i,j,1)=ru(i,j,1) -DC(i,1)
        enddo
        if (j.ge.jstrV) then
          do i=istr,iend
            DC(i,1)=0.5625_8*(Hz(i  ,j,1)+Hz(i,j-1,1))
     &             -0.0625_8*(Hz(i,j+1,1)+Hz(i,j-2,1))
            FC(i,0)=1.5_8*v(i,j,1,nrhs)
            CF(i,1)=0.5_8
          enddo
          do k=1,N-1,+1
            do i=istr,iend
              DC(i,k+1)=0.5625_8*(Hz(i  ,j,k+1)+Hz(i,j-1,k+1))
     &                 -0.0625_8*(Hz(i,j+1,k+1)+Hz(i,j-2,k+1))
              cff=1._8/(2._8*DC(i,k)+DC(i,k+1)*(2._8-CF(i,k)))
              CF(i,k+1)=cff*DC(i,k)
              FC(i,k)=cff*( 3._8*( DC(i,k  )*v(i,j,k+1,nrhs)
     &                          +DC(i,k+1)*v(i,j,k  ,nrhs))
     &                                -DC(i,k+1)*FC(i,k-1))
            enddo
          enddo
          do i=istr,iend
            FC(i,N)=(3._8*v(i,j,N,nrhs)-FC(i,N-1))/(2._8-CF(i,N))
            DC(i,N)=0._8
          enddo
          do k=N-1,1,-1
            do i=istr,iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
              DC(i,k)=FC(i,k) * 0.5_8*( We(i,j,k)+We(i,j-1,k) -0.125_8*(
     &                         (We(i,j+1,k)-We(i,j  ,k))*vmask(i,j+1)
     &                        -(We(i,j-1,k)-We(i,j-2,k))*vmask(i,j-1)
     &                                                           ))
              rv(i,j,k+1)=rv(i,j,k+1) -DC(i,k+1)+DC(i,k)
            enddo
          enddo
          do i=istr,iend
            rv(i,j,1)=rv(i,j,1) -DC(i,1)
          enddo
        endif
        do i=istrU,iend
          DC(i,0)=dt*0.25_8*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
          FC(i,N-1)= 2._8*dt*(Akv(i,j,N-1)+Akv(i-1,j,N-1))
     &                      /(  Hz(i,j,N  )+Hz(i-1,j,N  )
     &                         +Hz(i,j,N-1)+Hz(i-1,j,N-1))
          WC(i,N-1)= DC(i,0)*0.5_8*(Wi(i,j,N-1)+Wi(i-1,j,N-1))
          cff=1._8/( 0.5_8*(Hz(i,j,N)+Hz(i-1,j,N))
     &                   +FC(i,N-1)-min(WC(i,N-1),0._8) )
          CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0._8) )
          stress = 0.5_8*(sustr(i,j)+sustr(i-1,j))
          DC(i,N)=cff*(u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)+dt*stress)
        enddo
        do k=N-1,2,-1
          do i=istrU,iend
            FC(i,k-1)= 2._8*dt*(Akv(i,j,k-1)+Akv(i-1,j,k-1))
     &                        /(  Hz(i,j,k  )+Hz(i-1,j,k  )
     &                           +Hz(i,j,k-1)+Hz(i-1,j,k-1))
            WC(i,k-1)= DC(i,0)*0.5_8*(Wi(i,j,k-1)+Wi(i-1,j,k-1))
            cff=1._8/( 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
     &                           +FC(i,k-1)-min(WC(i,k-1),0._8)
     &                             +FC(i,k)+max(WC(i,k),0._8)
     &                    -CF(i,k)*(FC(i,k)-min(WC(i,k),0._8))
     &                                                      )
            CF(i,k-1)=cff*(   FC(i,k-1)+max(WC(i,k-1),0._8) )
            DC(i,k)=cff*( u(i,j,k,nnew) +DC(i,0)*ru(i,j,k)
     &                 +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0._8)) )
          enddo
        enddo
        do i=istrU,iend
          stress = 0.5_8*(sustr(i,j)+sustr(i-1,j))
          DC(i,1)=( u(i,j,1,nnew) +DC(i,0)*ru(i,j,1)
     &                      +DC(i,2)*(FC(i,1)-min(WC(i,1),0._8))
     &                          )/( 0.5_8*(Hz(i,j,1)+Hz(i-1,j,1))
     &                            +0.5_8*dt*(r_D(i,j)+r_D(i-1,j))
     &                                +FC(i,1)+max(WC(i,1),0._8)
     &                       -CF(i,1)*(FC(i,1)-min(WC(i,1),0._8))
     &                                                       )
          u(i,j,1,nnew)=DC(i,1) * 0.5_8*(Hz(i,j,1)+Hz(i-1,j,1))
          rufrc(i,j)=ru(i,j,1) +dm_u(i,j)*dn_u(i,j)*( stress
     &                       -0.5_8*(r_D(i-1,j)+r_D(i,j))*DC(i,1) )
        enddo
        do k=2,N,+1
          do i=istrU,iend
            DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)
            u(i,j,k,nnew)=DC(i,k) * 0.5_8*(Hz(i,j,k)+Hz(i-1,j,k))
            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
          enddo
        enddo
        if (j.ge.jstrV) then
          do i=istr,iend
          stress = 0.5_8*(svstr(i,j)+svstr(i,j-1))
            DC(i,0)=dt*0.25_8*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            FC(i,N-1)= 2._8*dt*(Akv(i,j,N-1)+Akv(i,j-1,N-1))
     &                        /(  Hz(i,j,N  )+Hz(i,j-1,N  )
     &                           +Hz(i,j,N-1)+Hz(i,j-1,N-1))
            WC(i,N-1)= DC(i,0)*0.5_8*(Wi(i,j,N-1)+Wi(i,j-1,N-1))
            cff=1._8/( 0.5_8*(Hz(i,j,N)+Hz(i,j-1,N))
     &                     +FC(i,N-1)-min(WC(i,N-1),0._8) )
            CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0._8) )
            DC(i,N)=cff*(v(i,j,N,nnew)+DC(i,0)*rv(i,j,N)+dt*stress)
          enddo
          do k=N-1,2,-1
            do i=istr,iend
              FC(i,k-1)= 2._8*dt*(Akv(i,j,k-1)+Akv(i,j-1,k-1))
     &                          /(  Hz(i,j,k  )+Hz(i,j-1,k  )
     &                             +Hz(i,j,k-1)+Hz(i,j-1,k-1))
              WC(i,k-1)= DC(i,0)*0.5_8*(Wi(i,j,k-1)+Wi(i,j-1,k-1))
              cff=1._8/( 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
     &                              +FC(i,k-1)-min(WC(i,k-1),0._8)
     &                                +FC(i,k)+max(WC(i,k),0._8)
     &                       -CF(i,k)*(FC(i,k)-min(WC(i,k),0._8))
     &                                                        )
              CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0._8) )
              DC(i,k)=cff*( v(i,j,k,nnew) +DC(i,0)*rv(i,j,k)
     &                   +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0._8)) )
            enddo
          enddo
          do i=istr,iend
          stress = 0.5_8*(svstr(i,j)+svstr(i,j-1))
            DC(i,1)=( v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)
     &                       +DC(i,2)*(FC(i,1)-min(WC(i,1),0._8))
     &                            )/( 0.5_8*(Hz(i,j,1)+Hz(i,j-1,1))
     &                              +0.5_8*dt*(r_D(i,j)+r_D(i,j-1))
     &                                  +FC(i,1)+max(WC(i,1),0._8)
     &                         -CF(i,1)*(FC(i,1)-min(WC(i,1),0._8))
     &                                                          )
            v(i,j,1,nnew)=DC(i,1) * 0.5_8*(Hz(i,j,1)+Hz(i,j-1,1))
            rvfrc(i,j)=rv(i,j,1) +dm_v(i,j)*dn_v(i,j)*( stress
     &                         -0.5_8*(r_D(i,j)+r_D(i,j-1))*DC(i,1) )
          enddo
          do k=2,N,+1
            do i=istr,iend
              DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)
              v(i,j,k,nnew)=DC(i,k) * 0.5_8*(Hz(i,j,k)+Hz(i,j-1,k))
              rvfrc(i,j)=rvfrc(i,j) + rv(i,j,k)
            enddo
          enddo
        endif
      enddo
      return
      end
      subroutine check_step_uv1_switches (ierr)
      implicit none
      integer(kind=4) ierr, is,ie, lenstr
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
      ie=lenstr(cpps)
      is=ie+2
      ie=is+13
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='<step3d_uv1.F>'
      is=ie+2
      ie=is+10
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='UPSTREAM_UV'
      is=ie+2
      ie=is+8
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='SPLINE_UV'
      is=ie+2
      ie=is+9
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='NEUMANN_UV'
      return
  99  if (mynode.eq.0) write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      return
      end
