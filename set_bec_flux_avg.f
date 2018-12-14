      subroutine set_bec_flux_avg(tile)
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
      call set_bec_flux_avg_tile(Istr,Iend,Jstr,Jend)
      return
      end
         subroutine set_bec_flux_avg_tile (Istr,Iend,Jstr,Jend)
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
       real(kind=8) tracer(-1:Lm+2+padd_X,N,ntrc_bio,2)
        common /tracers/ tracer
        real(kind=8) tracer_sed(-1:Lm+2+padd_X,NT_sed,2)
        common /tracer_sed/ tracer_sed
        real(kind=8) ifrac(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &    press(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
        common /fic_ap/ifrac,press
        real(kind=8) PH_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   pCO2sw(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   pCO2air(-1:Lm+2+padd_X,-1:Mm+2+padd_E) ,
     &   PAR(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   PARinc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
     &   ,PARinc_rst(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       common /time_averaging/
     &       ph_hist,pCO2sw,pCO2air,
     &       PAR,PARinc
     &  , PARinc_rst
        real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E) :: 
     &                     dom_sp_sfc, dom_diat_sfc,
     &       dom_diaz_sfc, dom_sp_int, dom_diat_int, dom_diaz_int,
     &       spchl_int, diatchl_int, diazchl_int
        common /specdom/ dom_sp_sfc, dom_diat_sfc,
     &       dom_diaz_sfc, dom_sp_int, dom_diat_int, dom_diaz_int,
     &       spchl_int, diatchl_int, diazchl_int
       real(kind=8) WS_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   XKW_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   AP_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   SCHMIDT_O2_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   O2SAT_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   FG_O2_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &    SCHMIDT_CO2_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   CO2STAR_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   DCO2STAR_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &    FG_CO2_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &   IRON_FLUX_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       real(kind=8) fg_n2o_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       real(kind=8) fg_n2_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
        real(kind=8)
     &    PO4_RESTORE_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    NO3_RESTORE_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    SiO3_RESTORE_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   PO4STAR_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    POC_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    POC_PROD_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   POC_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    CaCO3_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    CaCO3_PROD_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    CaCO3_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    SiO2_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    SiO2_PROD_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    SiO2_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    dust_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
          real(kind=8)  dust_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,
     &                                N),
     &    P_iron_FLUX_IN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    P_iron_PROD_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    P_iron_REMIN_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    graze_sp_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    graze_diat_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    graze_tot_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_loss_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_loss_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    zoo_loss_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_agg_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_agg_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
          real(kind=8)  photoC_sp_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    f_ratio_sp_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    photoC_diat_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    f_ratio_diat_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    tot_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    no3_v_sp_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    nh4_v_sp_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    no3_v_diat_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    nh4_v_diat_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DOC_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DOC_remin_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    Fe_scavenge_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_N_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_Fe_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &     sp_PO4_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    sp_light_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_N_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &     diat_Fe_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_PO4_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diat_SiO3_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
         real(kind=8) diat_light_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,
     &                                N),
     &   CaCO3_form_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diaz_Nfix_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    graze_diaz_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diaz_loss_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
        real(kind=8) photoC_diaz_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diaz_P_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    diaz_Fe_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &     diaz_light_lim_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &     Fe_scavenge_rate_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DON_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DON_remin_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DOFe_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
         real(kind=8)  DOFe_remin_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   DOP_prod_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    DOP_remin_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    bSI_form_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    photoFe_diaz_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    photoFe_diat_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    photoFe_sp_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &    nitrif_HIST(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
         real(kind=8) ammox_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &        denitr_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &        n2o_prod_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &        n2_prod_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &        denitr_sed_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       common /time_averaging1/WS_HIST, XKW_HIST,
     &   AP_HIST, SCHMIDT_O2_HIST, O2SAT_HIST, FG_O2_HIST,
     &    SCHMIDT_CO2_HIST, CO2STAR_HIST,
     &    DCO2STAR_HIST,
     &    FG_CO2_HIST, IRON_FLUX_HIST,
     &    PO4_RESTORE_HIST, NO3_RESTORE_HIST,
     &    SiO3_RESTORE_HIST, PO4STAR_HIST,
     &    POC_FLUX_IN_HIST, POC_PROD_HIST, POC_REMIN_HIST,
     &    CaCO3_FLUX_IN_HIST, CaCO3_PROD_HIST,
     &    CaCO3_REMIN_HIST,  SiO2_FLUX_IN_HIST,
     &    SiO2_PROD_HIST, SiO2_REMIN_HIST, dust_FLUX_IN_HIST,
     &    dust_REMIN_HIST, P_iron_FLUX_IN_HIST,
     &    P_iron_PROD_HIST, P_iron_REMIN_HIST,
     &    graze_sp_HIST, graze_diat_HIST, graze_tot_HIST,
     &    sp_loss_HIST, diat_loss_HIST, zoo_loss_HIST,
     &    sp_agg_HIST, diat_agg_HIST,
     &    photoC_sp_HIST, f_ratio_sp_hist,
     &    photoC_diat_HIST, f_ratio_diat_hist, tot_prod_HIST,
     &    no3_v_sp_hist, nh4_v_sp_hist,
     &    no3_v_diat_hist, nh4_v_diat_hist,
     &    DOC_prod_HIST, DOC_remin_HIST, Fe_scavenge_HIST
     &        , fg_n2o_hist
     &        , fg_n2_hist
     &        , denitr_sed_hist
       common /time_averaging2/
     &    sp_N_lim_HIST, sp_Fe_lim_HIST, sp_PO4_lim_HIST,
     &    sp_light_lim_HIST, diat_N_lim_HIST, diat_Fe_lim_HIST,
     &    diat_PO4_lim_HIST, diat_SiO3_lim_HIST,
     &    diat_light_lim_HIST, CaCO3_form_HIST,
     &    diaz_Nfix_HIST, graze_diaz_HIST, diaz_loss_HIST,
     &     photoC_diaz_HIST, diaz_P_lim_HIST,
     &    diaz_Fe_lim_HIST, diaz_light_lim_HIST,
     &     Fe_scavenge_rate_HIST, DON_prod_HIST,
     &    DON_remin_HIST, DOFe_prod_HIST,
     &    DOFe_remin_HIST, DOP_prod_HIST,
     &    DOP_remin_HIST, bSI_form_HIST,
     &    photoFe_diaz_HIST, photoFe_diat_HIST, photoFe_sp_HIST,
     &    nitrif_HIST
     &    , ammox_hist, denitr_hist, n2o_prod_hist, n2_prod_hist
       real(kind=8) bot_flux_poc_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &      bot_flux_caco3_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &      bot_flux_si_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &      bot_flux_fe_hist(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       common /sed_flux/ bot_flux_poc_hist, bot_flux_caco3_hist,
     &      bot_flux_si_hist,bot_flux_fe_hist
       integer(kind=4)  po4_ind , no3_ind,sio3_ind, nh4_ind,fe_ind,
     & o2_ind, dic_ind,alk_ind,doc_ind,spC_ind,spChl_ind,
     & spCaCO3_ind,diatC_ind,diatChl_ind,zooC_ind,spFe_ind,
     &  diatSi_ind,diatFe_ind,diazC_ind,diazChl_ind, diazFe_ind,
     &  don_ind,dofe_ind,dop_ind
       parameter (po4_ind=1 , no3_ind=2,sio3_ind=3, nh4_ind=4,
     &  fe_ind=5,o2_ind=6, dic_ind=7,alk_ind=8,doc_ind=9,
     &  spC_ind=10,spChl_ind=11, spCaCO3_ind=12,diatC_ind=13,
     &  diatChl_ind=14,zooC_ind=15,spFe_ind=16,
     &  diatSi_ind=17,diatFe_ind=18,diazC_ind=19,
     &  diazChl_ind=20, diazFe_ind=21,
     &  don_ind=22,dofe_ind=23,dop_ind=24)
       integer(kind=4), parameter :: no2_ind = dop_ind + 1
       integer(kind=4), parameter :: n2o_ind = no2_ind + 1
       integer(kind=4), parameter :: n2_ind = n2o_ind + 1
       integer(kind=4), parameter :: sed_poc_ind = 1
       integer(kind=4), parameter :: sed_caco3_ind = 2
       integer(kind=4), parameter :: sed_si_ind = 3
       logical lsource_sink,lflux_gas_o2, lflux_gas_co2,
     &  liron_flux,ldust_flux
        common /ecoflag/lsource_sink,lflux_gas_o2,lflux_gas_co2,
     &   liron_flux,ldust_flux
       logical lrest_po4,lrest_no3,lrest_sio3
       real(kind=8) po4_clim(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   no3_clim(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N),
     &   sio3_clim(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
           real(kind=8) nutr_rest_time_inv(N)
        common /restore_flag/lrest_po4,lrest_no3,lrest_sio3
        common /restore_clim/po4_clim,
     &      no3_clim,sio3_clim,nutr_rest_time_inv
      real(kind=8) t_sed(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT_sed)
      common /ocean_t_sed/t_sed
        real(kind=8) sinking_particle_POC(6,-1:Lm+2+padd_X,
     &                          -1:Mm+2+padd_E),
     & sinking_particle_P_CaCO3(6,-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     & sinking_particle_P_sio2(6,-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     & sinking_particle_dust(6,-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     & sinking_particle_P_iron(6,-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       real(kind=8) diss(5),gamma(5),mass(5),rhoo(5)
      common /sinking_part/sinking_particle_POC,
     &  sinking_particle_P_CaCO3,sinking_particle_P_SiO2,
     &  sinking_particle_dust,sinking_particle_P_iron,
     &  diss,gamma,mass,rhoo
        logical landmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
        common /calcation/landmask
       real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E)::
     &     PH_AVG, pCO2_AVG, pCO2air_AVG, PARinc_avg
        real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)::
     &    PAR_avg
        common /time_avg/
     &    PH_AVG, pCO2_AVG, pCO2air_AVG, PARinc_avg,
     &    PAR_avg
        real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E) ::
     &       dom_sp_sfc_avg, dom_diat_sfc_avg, dom_diaz_sfc_avg,
     &       dom_sp_int_avg, dom_diat_int_avg, dom_diaz_int_avg
        common /specdom_avg/
     &       dom_sp_sfc_avg, dom_diat_sfc_avg, dom_diaz_sfc_avg,
     &       dom_sp_int_avg, dom_diat_int_avg, dom_diaz_int_avg
      real(kind=8) t_sed_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT_sed)
      common /ocean_t_sed_avg/t_sed_avg
       real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E):: WS_AVG, 
     &                              XKW_AVG,
     &     AP_AVG, SCHMIDT_O2_AVG, O2SAT_AVG, FG_O2_AVG,
     &     SCHMIDT_CO2_AVG, CO2STAR_AVG, DCO2STAR_AVG,
     &     FG_CO2_AVG, IRON_FLUX_AVG,
     &     PARinc_flux_avg, zeta_bgc_flux_avg
     &     , fg_n2o_avg, fg_n2_avg
        real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)::
     &    PO4_RESTORE_AVG, NO3_RESTORE_AVG,
     &    SiO3_RESTORE_AVG, PAR_flux_avg, PO4STAR_AVG,
     &    POC_FLUX_IN_AVG, POC_PROD_AVG, POC_REMIN_AVG,
     &    CaCO3_FLUX_IN_AVG, CaCO3_PROD_AVG,
     &    CaCO3_REMIN_AVG,  SiO2_FLUX_IN_AVG,
     &    SiO2_PROD_AVG, SiO2_REMIN_AVG, dust_FLUX_IN_AVG,
     &    dust_REMIN_AVG, P_iron_FLUX_IN_AVG,
     &    P_iron_PROD_AVG, P_iron_REMIN_AVG,
     &    graze_sp_AVG, graze_diat_AVG, graze_tot_AVG,
     &    sp_loss_AVG, diat_loss_AVG, zoo_loss_AVG,
     &    sp_agg_AVG, diat_agg_AVG,
     &    photoC_sp_AVG, f_ratio_sp_avg,
     &    photoC_diat_AVG, f_ratio_diat_avg, tot_prod_AVG,
     &    no3_v_sp_avg, nh4_v_sp_avg,
     &    no3_v_diat_avg, nh4_v_diat_avg,
     &    DOC_prod_AVG, DOC_remin_AVG, Fe_scavenge_AVG,
     &    sp_N_lim_AVG, sp_Fe_lim_AVG, sp_PO4_lim_AVG,
     &    sp_light_lim_AVG, diat_N_lim_AVG, diat_Fe_lim_AVG,
     &    diat_PO4_lim_AVG, diat_SiO3_lim_AVG
        real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)::
     &    diat_light_lim_AVG, CaCO3_form_AVG,
     &    diaz_Nfix_AVG, graze_diaz_AVG, diaz_loss_AVG,
     &     photoC_diaz_AVG, diaz_P_lim_AVG,
     &    diaz_Fe_lim_AVG, diaz_light_lim_AVG,
     &     Fe_scavenge_rate_AVG, DON_prod_AVG,
     &    DON_remin_AVG, DOFe_prod_AVG,
     &    DOFe_remin_AVG, DOP_prod_AVG,
     &    DOP_remin_AVG, bSI_form_AVG,
     &    photoFe_diaz_AVG, photoFe_diat_AVG,
     &    photoFe_sp_AVG,nitrif_AVG
         real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N) ::  
     &                             ammox_avg,
     &        denitr_avg, n2o_prod_avg, n2_prod_avg
         real(kind=8),dimension(-1:Lm+2+padd_X,-1:Mm+2+padd_E) :: 
     &                           denitr_sed_avg
       common /time_avg1/WS_AVG, XKW_AVG,
     &    AP_AVG, SCHMIDT_O2_AVG, O2SAT_AVG, FG_O2_AVG,
     &    SCHMIDT_CO2_AVG, CO2STAR_AVG, DCO2STAR_AVG,
     &    FG_CO2_AVG, IRON_FLUX_AVG,
     %    PARinc_flux_avg, zeta_bgc_flux_avg,
     &    PO4_RESTORE_AVG, NO3_RESTORE_AVG,
     &    SiO3_RESTORE_AVG, PAR_flux_avg, PO4STAR_AVG,
     &    POC_FLUX_IN_AVG, POC_PROD_AVG, POC_REMIN_AVG,
     &    CaCO3_FLUX_IN_AVG, CaCO3_PROD_AVG,
     &    CaCO3_REMIN_AVG,  SiO2_FLUX_IN_AVG,
     &    SiO2_PROD_AVG, SiO2_REMIN_AVG, dust_FLUX_IN_AVG,
     &    dust_REMIN_AVG, P_iron_FLUX_IN_AVG,
     &    P_iron_PROD_AVG, P_iron_REMIN_AVG,
     &    graze_sp_AVG, graze_diat_AVG, graze_tot_AVG,
     &    sp_loss_AVG, diat_loss_AVG, zoo_loss_AVG
     &     , fg_n2o_avg, fg_n2_avg
       common /time_avg2/
     &    sp_agg_AVG, diat_agg_AVG,
     &    photoC_sp_AVG, f_ratio_sp_avg,
     &    photoC_diat_AVG, f_ratio_diat_avg, tot_prod_AVG,
     &    no3_v_sp_avg, nh4_v_sp_avg,
     &    no3_v_diat_avg, nh4_v_diat_avg,
     &    DOC_prod_AVG, DOC_remin_AVG, Fe_scavenge_AVG,
     &    sp_N_lim_AVG, sp_Fe_lim_AVG, sp_PO4_lim_AVG,
     &    sp_light_lim_AVG, diat_N_lim_AVG, diat_Fe_lim_AVG,
     &    diat_PO4_lim_AVG, diat_SiO3_lim_AVG,
     &    diat_light_lim_AVG, CaCO3_form_AVG,
     &    diaz_Nfix_AVG, graze_diaz_AVG, diaz_loss_AVG,
     &     photoC_diaz_AVG, diaz_P_lim_AVG,
     &    diaz_Fe_lim_AVG, diaz_light_lim_AVG,
     &     Fe_scavenge_rate_AVG, DON_prod_AVG,
     &    DON_remin_AVG, DOFe_prod_AVG,
     &    DOFe_remin_AVG, DOP_prod_AVG,
     &    DOP_remin_AVG, bSI_form_AVG,
     &    photoFe_diaz_AVG, photoFe_diat_AVG,
     &    photoFe_sp_AVG, nitrif_AVG
     &      , ammox_avg, denitr_avg, n2o_prod_avg, n2_prod_avg,
     &      denitr_sed_avg
       real(kind=8) bot_flux_poc_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &      bot_flux_caco3_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &      bot_flux_si_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E),
     &      bot_flux_fe_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
       common /sed_flux/ bot_flux_poc_avg, bot_flux_caco3_avg,
     &      bot_flux_si_avg,bot_flux_fe_avg
      logical new_bgc_flux_his
      integer(kind=4) n_bgc_flux_his, nrpf_bgc_flux_his
     &     , ncid_bgc_flux_his, nrec_bgc_flux_his
      common /scalars_bgc/
     &     new_bgc_flux_his
     &     , n_bgc_flux_his, nrpf_bgc_flux_his
     &     , ncid_bgc_flux_his, nrec_bgc_flux_his
      integer(kind=4) ncid_bgc_flux_avg, nrec_bgc_flux_avg
     &     , nrpf_bgc_flux_avg
      common /scalars_bgc_avg/ ncid_bgc_flux_avg, nrec_bgc_flux_avg
     &     , nrpf_bgc_flux_avg
      character(len=80) bgc_flux_avg_name
      common /c_bgcflux_avg/ bgc_flux_avg_name
      integer(kind=4), parameter :: num_bgcflux_2d = 12
     &     + 3
     &     + NT_sed
      integer(kind=4), parameter :: num_bgcflux = 80
     &     + 7
     &     + NT_sed
      integer(kind=4), dimension(num_bgcflux) :: vid_bec_flux_his
      common /c_bgcflux_bec/ vid_bec_flux_his
      character(len=80) bgc_flux_his_name,
     &     vname_bgcflux(3, num_bgcflux)
      common /c_bgcflux/ bgc_flux_his_name, vname_bgcflux
      integer(kind=4) bgc_flux_hisTime, bgc_flux_hisTstep
     &     , bgc_flux_hisZ
      common /ncids_bgc_flux/ bgc_flux_hisTime, bgc_flux_hisTstep
     &     , bgc_flux_hisZ
      logical new_bgc_flux_avg
      integer(kind=4) nts_bgc_flux_avg, n_bgc_flux_avg
      common /scalars_bgc_avg/
     &     new_bgc_flux_avg,
     &     nts_bgc_flux_avg, n_bgc_flux_avg
      real(kind=8) time_bgc_flux_avg
      common /scalars_bgc_avg_real/ time_bgc_flux_avg
      integer(kind=4), dimension(num_bgcflux) :: vid_bec_flux_avg
      integer(kind=4) :: bgc_flux_avgTstep, bgc_flux_avgTime,
     &    bgc_flux_avgZ
      common /c_bgcflux_avg_bec/ vid_bec_flux_avg, bgc_flux_avgTstep,
     &    bgc_flux_avgTime, bgc_flux_avgZ
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
      real(kind=8) zeta(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real(kind=8) ubar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real(kind=8) vbar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      common /ocean_zeta/zeta /ocean_ubar/ubar /ocean_vbar/vbar
         integer(kind=4) istr,jstr,iend,jend
         real(kind=8) cff, cff1
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
         if (n_bgc_flux_avg < 1) return
         if (iic.gt.nts_bgc_flux_avg) then
            if (n_bgc_flux_avg .eq. 1 .or.
     &           mod(iic-nts_bgc_flux_avg,n_bgc_flux_avg).eq.1) then
               cff =1.0_8
               cff1=0.0_8
               if ((istr.eq.iwest .and. jstr.eq.jsouth)) then
                  time_bgc_flux_avg=time
                  if (mynode.eq.0) write(*,'(6x,A,I11,2(X,A,I7))')
     &                 'SET_BGC_FLUX_AVG -- Started averaging at iic=',
     &                 iic,' : nts_bgc_flux_avg=',nts_bgc_flux_avg,
     &                 'n_bgc_flux_avg=',n_bgc_flux_avg
               endif
            elseif (mod(iic-nts_bgc_flux_avg,n_bgc_flux_avg).gt.1) then
               cff =1.0_8
               cff1=1.0_8
               if ((istr.eq.iwest .and. jstr.eq.jsouth)) 
     &              time_bgc_flux_avg=time_bgc_flux_avg+time
            elseif (mod(iic-nts_bgc_flux_avg,n_bgc_flux_avg).eq.0) then
               cff=1._8/dble(n_bgc_flux_avg)
               cff1=1.0_8
               if ((istr.eq.iwest .and. jstr.eq.jsouth)) then
                  time_bgc_flux_avg=cff*(time_bgc_flux_avg+time)
                  if (mynode.eq.0) write(*,'(6x,A,I11,2(X,A,I7))')
     &                 'SET_BGC_FLUX_AVG -- Finished averaging at iic=',
     &                 iic,' : nts_bgc_flux_avg=',nts_bgc_flux_avg,
     &                 'n_bgc_flux_avg=',n_bgc_flux_avg
               endif
            endif
            zeta_bgc_flux_avg(istrR:iendR,jstrR:jendR) =
     &        cff * ( cff1*zeta_bgc_flux_avg(istrR:iendR,jstrR:jendR) +
     &        zeta(istrR:iendR,jstrR:jendR,knew) )
            WS_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*WS_AVG(istrR:iendR,jstrR:jendR)
     &           + WS_hist(istrR:iendR,jstrR:jendR) )
            XKW_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*XKW_AVG(istrR:iendR,jstrR:jendR)
     &           + XKW_hist(istrR:iendR,jstrR:jendR) )
            ap_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*ap_AVG(istrR:iendR,jstrR:jendR)
     &           + ap_hist(istrR:iendR,jstrR:jendR) )
            SCHMIDT_O2_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*SCHMIDT_O2_AVG(istrR:iendR,jstrR:jendR)
     &           + SCHMIDT_O2_hist(istrR:iendR,jstrR:jendR) )
            O2SAT_AVG(istrR:iendR,jstrR:jendR)=
     &           cff * ( cff1*O2SAT_AVG(istrR:iendR,jstrR:jendR)
     &           +O2SAT_hist(istrR:iendR,jstrR:jendR) )
            FG_O2_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*FG_O2_AVG(istrR:iendR,jstrR:jendR)
     &           +FG_O2_hist(istrR:iendR,jstrR:jendR) )
            FG_N2O_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*FG_N2O_avg(istrR:iendR,jstrR:jendR)
     &           + FG_N2O_hist(istrR:iendR,jstrR:jendR) )
            FG_N2_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*FG_N2_avg(istrr:iendR,jstrR:jendR)
     &           + FG_N2_hist(istrR:iendR,jstrR:jendR) )
            SCHMIDT_CO2_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*SCHMIDT_CO2_AVG(istrR:iendR,jstrR:jendR)
     &           + SCHMIDT_CO2_hist(istrR:iendR,jstrR:jendR) )
            CO2STAR_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*CO2STAR_AVG(istrR:iendR,jstrR:jendR)
     &           + CO2STAR_hist(istrR:iendR,jstrR:jendR) )
            DCO2STAR_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*DCO2STAR_AVG(istrR:iendR,jstrR:jendR)
     &           + DCO2STAR_hist(istrR:iendR,jstrR:jendR) )
            FG_CO2_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*FG_CO2_AVG(istrR:iendR,jstrR:jendR)
     &           + FG_CO2_hist(istrR:iendR,jstrR:jendR) )
            IRON_FLUX_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*IRON_FLUX_AVG(istrR:iendR,jstrR:jendR)
     &           + IRON_FLUX_hist(istrR:iendR,jstrR:jendR) )
            PARinc_flux_AVG(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*PARinc_flux_AVG(istrR:iendR,jstrR:jendR)
     &           + PARinc(istrR:iendR,jstrR:jendR) )
            bot_flux_poc_avg(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*bot_flux_poc_avg(istrR:iendR,jstrR:jendR)
     &           + bot_flux_poc_hist(istrR:iendR,jstrR:jendR) )
            bot_flux_caco3_avg(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*bot_flux_caco3_avg(istrR:iendR,
     &                            jstrR:jendR)
     &           + bot_flux_caco3_hist(istrR:iendR,jstrR:jendR) )
            bot_flux_si_avg(istrR:iendR,jstrR:jendR) =
     &           cff * ( cff1*bot_flux_si_avg(istrR:iendR,jstrR:jendR)
     &           + bot_flux_si_hist(istrR:iendR,jstrR:jendR) )
            PO4_RESTORE_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff * ( cff1*PO4_RESTORE_AVG(istrR:iendR,jstrR:jendR,:)
     &           + PO4_RESTORE_hist(istrR:iendR,jstrR:jendR,:) )
            NO3_RESTORE_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff * ( cff1*NO3_RESTORE_AVG(istrR:iendR,jstrR:jendR,:)
     &           + NO3_RESTORE_hist(istrR:iendR,jstrR:jendR,:) )
            SiO3_RESTORE_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff * (cff1*SiO3_RESTORE_AVG(istrR:iendR,jstrR:jendR,:)
     &           + SiO3_RESTORE_hist(istrR:iendR,jstrR:jendR,:) )
            PAR_flux_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff * (cff1*PAR_flux_AVG(istrR:iendR,jstrR:jendR,:)
     &           + PAR(istrR:iendR,jstrR:jendR,:) )
            PO4STAR_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff * (cff1*PO4STAR_AVG(istrR:iendR,jstrR:jendR,:)
     &           + PO4STAR_hist(istrR:iendR,jstrR:jendR,:) )
            POC_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff * (cff1*POC_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:)
     &           + POC_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            POC_PROD_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff * (cff1*POC_PROD_AVG(istrR:iendR,jstrR:jendR,:)
     &           + POC_PROD_hist(istrR:iendR,jstrR:jendR,:) )
            POC_REMIN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff * (cff1*POC_REMIN_AVG(istrR:iendR,jstrR:jendR,:)
     &           + POC_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            CaCO3_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff *(cff1*CaCO3_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:)
     &           + CaCO3_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            CaCO3_PROD_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff *(cff1*CaCO3_PROD_AVG(istrR:iendR,jstrR:jendR,:)
     &           + CaCO3_PROD_hist(istrR:iendR,jstrR:jendR,:) )
            CaCO3_REMIN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff *(cff1*CaCO3_REMIN_AVG(istrR:iendR,jstrR:jendR,:) +
     &           CaCO3_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            SiO2_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff *(cff1*SiO2_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:)+
     &           SiO2_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            SiO2_PROD_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff *(cff1*SiO2_PROD_AVG(istrR:iendR,jstrR:jendR,:)
     &           + SiO2_PROD_hist(istrR:iendR,jstrR:jendR,:) )
            SiO2_REMIN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff *(cff1*SiO2_REMIN_AVG(istrR:iendR,jstrR:jendR,:)
     &           + SiO2_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            dust_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff *(cff1*dust_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:)
     &           + dust_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            dust_REMIN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff *(cff1*dust_REMIN_AVG(istrR:iendR,jstrR:jendR,:)
     &           +dust_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            P_iron_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*P_iron_FLUX_IN_AVG(istrR:iendR,jstrR:jendR,:)
     &          +P_iron_FLUX_IN_hist(istrR:iendR,jstrR:jendR,:) )
            P_iron_PROD_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*P_iron_PROD_AVG(istrR:iendR,jstrR:jendR,:)
     &             + P_iron_PROD_hist(istrR:iendR,jstrR:jendR,:) )
            P_iron_REMIN_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*P_iron_REMIN_AVG(istrR:iendR,jstrR:jendR,:)
     &           + P_iron_REMIN_hist(istrR:iendR,jstrR:jendR,:) )
            graze_sp_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*graze_sp_AVG(istrR:iendR,jstrR:jendR,:)
     &           + graze_sp_hist(istrR:iendR,jstrR:jendR,:) )
            graze_diat_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*graze_diat_AVG(istrR:iendR,jstrR:jendR,:)
     &           + graze_diat_hist(istrR:iendR,jstrR:jendR,:) )
            graze_tot_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*graze_tot_AVG(istrR:iendR,jstrR:jendR,:)
     &           + graze_tot_hist(istrR:iendR,jstrR:jendR,:) )
            sp_loss_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*sp_loss_AVG(istrR:iendR,jstrR:jendR,:)
     &           +sp_loss_hist(istrR:iendR,jstrR:jendR,:) )
            diat_loss_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diat_loss_AVG(istrR:iendR,jstrR:jendR,:)
     &           +diat_loss_hist(istrR:iendR,jstrR:jendR,:) )
            zoo_loss_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*zoo_loss_AVG(istrR:iendR,jstrR:jendR,:)
     &           + zoo_loss_hist(istrR:iendR,jstrR:jendR,:) )
            sp_agg_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*sp_agg_AVG(istrR:iendR,jstrR:jendR,:)
     &           +sp_agg_hist(istrR:iendR,jstrR:jendR,:) )
            diat_agg_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diat_agg_AVG(istrR:iendR,jstrR:jendR,:)
     &           +diat_agg_hist(istrR:iendR,jstrR:jendR,:) )
            photoC_sp_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*photoC_sp_AVG(istrR:iendR,jstrR:jendR,:)
     &           + photoC_sp_hist(istrR:iendR,jstrR:jendR,:) )
            f_ratio_sp_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*f_ratio_sp_AVG(istrR:iendR,jstrR:jendR,:)
     &           + f_ratio_sp_hist(istrR:iendR,jstrR:jendR,:) )
            photoC_diat_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*photoC_diat_AVG(istrR:iendR,jstrR:jendR,:)
     &           + photoC_diat_hist(istrR:iendR,jstrR:jendR,:) )
            f_ratio_diat_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*f_ratio_diat_AVG(istrR:iendR,jstrR:jendR,:)
     &           + f_ratio_diat_hist(istrR:iendR,jstrR:jendR,:) )
            tot_prod_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*tot_prod_AVG(istrR:iendR,jstrR:jendR,:)
     &           +tot_prod_hist(istrR:iendR,jstrR:jendR,:) )
            no3_v_sp_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*no3_v_sp_AVG(istrR:iendR,jstrR:jendR,:)
     &           +no3_v_sp_hist(istrR:iendR,jstrR:jendR,:) )
            nh4_v_sp_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*nh4_v_sp_AVG(istrR:iendR,jstrR:jendR,:)
     &           +nh4_v_sp_hist(istrR:iendR,jstrR:jendR,:) )
            no3_v_diat_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*no3_v_diat_AVG(istrR:iendR,jstrR:jendR,:)
     &           +no3_v_diat_hist(istrR:iendR,jstrR:jendR,:) )
            nh4_v_diat_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*nh4_v_diat_AVG(istrR:iendR,jstrR:jendR,:)
     &           +nh4_v_diat_hist(istrR:iendR,jstrR:jendR,:) )
            DOC_prod_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*DOC_prod_AVG(istrR:iendR,jstrR:jendR,:)
     &           + DOC_prod_hist(istrR:iendR,jstrR:jendR,:) )
            DOC_remin_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*DOC_remin_AVG(istrR:iendR,jstrR:jendR,:)
     &           +DOC_remin_hist(istrR:iendR,jstrR:jendR,:) )
            Fe_scavenge_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*Fe_scavenge_AVG(istrR:iendR,jstrR:jendR,:)
     &           + Fe_scavenge_hist(istrR:iendR,jstrR:jendR,:) )
            sp_N_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*sp_N_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           + sp_N_lim_hist(istrR:iendR,jstrR:jendR,:) )
            sp_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*sp_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           + sp_Fe_lim_hist(istrR:iendR,jstrR:jendR,:) )
            sp_PO4_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*sp_PO4_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           + sp_PO4_lim_hist(istrR:iendR,jstrR:jendR,:) )
            sp_light_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*sp_light_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           + sp_light_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_N_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diat_N_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           +diat_N_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diat_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           +diat_Fe_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_PO4_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diat_PO4_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           +diat_PO4_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_SiO3_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diat_SiO3_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           + diat_SiO3_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diat_light_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diat_light_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           + diat_light_lim_hist(istrR:iendR,jstrR:jendR,:) )
            CaCO3_form_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*CaCO3_form_AVG(istrR:iendR,jstrR:jendR,:)
     &           +CaCO3_form_hist(istrR:iendR,jstrR:jendR,:) )
            diaz_Nfix_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diaz_Nfix_AVG(istrR:iendR,jstrR:jendR,:)
     &           +diaz_Nfix_hist(istrR:iendR,jstrR:jendR,:) )
            graze_diaz_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*graze_diaz_AVG(istrR:iendR,jstrR:jendR,:)
     &           + graze_diaz_hist(istrR:iendR,jstrR:jendR,:) )
            diaz_loss_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diaz_loss_AVG(istrR:iendR,jstrR:jendR,:)
     &           + diaz_loss_hist(istrR:iendR,jstrR:jendR,:) )
            photoC_diaz_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*photoC_diaz_AVG(istrR:iendR,jstrR:jendR,:)
     &           + photoC_diaz_hist(istrR:iendR,jstrR:jendR,:) )
            diaz_P_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diaz_P_lim _AVG(istrR:iendR,jstrR:jendR,:)
     &           + diaz_P_lim_hist(istrR:iendR,jstrR:jendR,:) )
            diaz_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diaz_Fe_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           +diaz_Fe_lim _hist(istrR:iendR,jstrR:jendR,:) )
            diaz_light_lim_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*diaz_light_lim_AVG(istrR:iendR,jstrR:jendR,:)
     &           + diaz_light_lim_hist(istrR:iendR,jstrR:jendR,:) )
            Fe_scavenge_rate_AVG(istrR:iendR,jstrR:jendR,:) =
     &         cff*(cff1*Fe_scavenge_rate_AVG(istrR:iendR,jstrR:jendR,:)
     &           + Fe_scavenge_rate_hist(istrR:iendR,jstrR:jendR,:) )
            DON_prod_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*DON_prod_AVG(istrR:iendR,jstrR:jendR,:) +
     &           DON_prod_hist(istrR:iendR,jstrR:jendR,:) )
            DON_remin_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*DON_remin_AVG(istrR:iendR,jstrR:jendR,:)
     &           + DON_remin_hist(istrR:iendR,jstrR:jendR,:) )
            DOFe_prod_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*DOFe_prod_AVG(istrR:iendR,jstrR:jendR,:)
     &           + DOFe_prod_hist(istrR:iendR,jstrR:jendR,:) )
            DOFe_remin_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*DOFe_remin_AVG(istrR:iendR,jstrR:jendR,:)
     &           + DOFe_remin_hist(istrR:iendR,jstrR:jendR,:) )
            DOP_prod_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*DOP_prod_AVG(istrR:iendR,jstrR:jendR,:)
     &           + DOP_prod_hist(istrR:iendR,jstrR:jendR,:) )
            DOP_remin_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*DOP_remin_AVG(istrR:iendR,jstrR:jendR,:)
     &           + DOP_remin_hist(istrR:iendR,jstrR:jendR,:) )
            bSI_form_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*bSI_form_AVG(istrR:iendR,jstrR:jendR,:)
     &           + bSI_form_hist(istrR:iendR,jstrR:jendR,:) )
            photoFe_diaz_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*photoFe_diaz_AVG(istrR:iendR,jstrR:jendR,:)
     &           + photoFe_diaz_hist(istrR:iendR,jstrR:jendR,:) )
            photoFe_diat_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*photoFe_diat_AVG(istrR:iendR,jstrR:jendR,:) +
     &           photoFe_diat_hist(istrR:iendR,jstrR:jendR,:) )
            photoFe_sp_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*photoFe_sp_AVG(istrR:iendR,jstrR:jendR,:)
     &           + photoFe_sp_hist(istrR:iendR,jstrR:jendR,:) )
            nitrif_AVG(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*nitrif_AVG(istrR:iendR,jstrR:jendR,:)
     &           + nitrif_hist(istrR:iendR,jstrR:jendR,:) )
            ammox_avg(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*ammox_avg(istrR:iendR,jstrR:jendR,:)
     &           + ammox_hist(istrR:iendR,jstrR:jendR,:) )
            denitr_avg(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*denitr_avg(istrR:iendR,jstrR:jendR,:)
     &           + denitr_hist(istrR:iendR,jstrR:jendR,:) )
            n2o_prod_avg(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*n2o_prod_avg(istrR:iendR,jstrR:jendR,:)
     &           + n2o_prod_hist(istrR:iendR,jstrR:jendR,:) )
            n2_prod_avg(istrR:iendR,jstrR:jendR,:) =
     &           cff*(cff1*n2_prod_avg(istrR:iendR,jstrR:jendR,:)
     &           + n2_prod_hist(istrR:iendR,jstrR:jendR,:) )
            denitr_sed_avg(istrR:iendR,jstrR:jendR) =
     &           cff*(cff1*denitr_sed_avg(istrR:iendR,jstrR:jendR)
     &           + denitr_sed_hist(istrR:iendR,jstrR:jendR) )
         end if
         return
         end
