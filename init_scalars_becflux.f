       subroutine init_scalars_becflux()
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
       integer::ind
       ncid_bgc_flux_his=-1
       nrec_bgc_flux_his=0
       ncid_bgc_flux_avg=-1
       nrec_bgc_flux_avg=0
       ind=1
       vname_bgcflux(1,ind)='WS'
       vname_bgcflux(2,ind)='Wind speed'
       vname_bgcflux(3,ind)='m/s '
       ind=ind+1
       vname_bgcflux(1,ind)='XKW'
       vname_bgcflux(2,ind)='XKW_AVG'
       vname_bgcflux(3,ind)='m/s '
       ind=ind+1
       vname_bgcflux(1,ind)='ATM_PRESS'
       vname_bgcflux(2,ind)='Atmospheric pressure'
       vname_bgcflux(3,ind)='atm '
       ind=ind+1
       vname_bgcflux(1,ind)='SCHMIDT_O2'
       vname_bgcflux(2,ind)='Schmidt number for O2'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='O2SAT'
       vname_bgcflux(2,ind)='O2SAT_AVG'
       vname_bgcflux(3,ind)='mmol/m3 '
       ind=ind+1
       vname_bgcflux(1,ind)='FG_O2'
       vname_bgcflux(2,ind)='Air-sea flux of O2'
       vname_bgcflux(3,ind)='mmol/m2/s'
       ind=ind+1
       vname_bgcflux(1,ind)='FG_N2O'
       vname_bgcflux(2,ind)='Air-sea flux of excess N2O'
       vname_bgcflux(3,ind)='mmol/m2/s'
       ind=ind+1
       vname_bgcflux(1,ind)='FG_N2'
       vname_bgcflux(2,ind)='Air-sea flux of excess N2'
       vname_bgcflux(3,ind)='mmol/m2/s'
       ind=ind+1
       vname_bgcflux(1,ind)='Sed_denitr'
       vname_bgcflux(2,ind)='Sediment denitrification'
       vname_bgcflux(3,ind)='mmol/m3/s'
       ind=ind+1
       vname_bgcflux(1,ind)='SCHMIDT_CO2'
       vname_bgcflux(2,ind)='Schmidt number for CO2'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='CO2STAR'
       vname_bgcflux(2,ind)='CO2STAR_AVG'
       vname_bgcflux(3,ind)='mmol/m3 '
       ind=ind+1
       vname_bgcflux(1,ind)='DCO2STAR'
       vname_bgcflux(2,ind)='DCO2STAR_AVG'
       vname_bgcflux(3,ind)='mmol/m3'
       ind=ind+1
       vname_bgcflux(1,ind)='FG_CO2'
       vname_bgcflux(2,ind)='Air-sea flux of CO2'
       vname_bgcflux(3,ind)='mmol/m2/s '
       ind=ind+1
       vname_bgcflux(1,ind)='IRON_FLUX'
       vname_bgcflux(2,ind)='Iron surface flux'
       vname_bgcflux(3,ind)='mmol/m2/s'
       ind=ind+1
       vname_bgcflux(1,ind)='PARinc'
       vname_bgcflux(2,ind)=
     &      'Inc. Photosynthetically available radiation'
       vname_bgcflux(3,ind)='W/m2'
       ind=ind+1
       vname_bgcflux(1,ind)='Sed_Flux_POC'
       vname_bgcflux(2,ind)=
     &      'Flux of POC into sediment'
       vname_bgcflux(3,ind)='mmol/m2'
       ind=ind+1
       vname_bgcflux(1,ind)='Sed_Flux_CaCO3'
       vname_bgcflux(2,ind)=
     &      'Flux of CaCO3 into sediment'
       vname_bgcflux(3,ind)='mmol/m2'
       ind=ind+1
       vname_bgcflux(1,ind)='Sed_Flux_Si'
       vname_bgcflux(2,ind)=
     &      'Flux of silicate into sediment'
       vname_bgcflux(3,ind)='mmol/m2'
       ind=ind+1
       vname_bgcflux(1,ind)='PO4_RESTORE'
       vname_bgcflux(2,ind)='PO4 restoring flux'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='NO3_RESTORE'
       vname_bgcflux(2,ind)='NO3 restoring flux'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='SIO3_RESTORE'
       vname_bgcflux(2,ind)='SiO3 restoring flux'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='PAR'
       vname_bgcflux(2,ind)='Photosynthetically available radiation'
       vname_bgcflux(3,ind)='W/m2'
       ind=ind+1
       vname_bgcflux(1,ind)='PO4STAR'
       vname_bgcflux(2,ind)='PO4STAR_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='POC_FLUX_IN'
       vname_bgcflux(2,ind)='POC_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
       ind=ind+1
       vname_bgcflux(1,ind)='POC_PROD'
       vname_bgcflux(2,ind)='POC_PROD_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='POC_REMIN'
       vname_bgcflux(2,ind)='POC remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s'
       ind=ind+1
       vname_bgcflux(1,ind)='CACO3_FLUX_IN'
       vname_bgcflux(2,ind)='CACO3_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
       ind=ind+1
       vname_bgcflux(1,ind)='CACO3_PROD'
       vname_bgcflux(2,ind)='CaCO3 production'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='CACO3_REMIN'
       vname_bgcflux(2,ind)='CaCO3 remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='SIO2_FLUX_IN'
       vname_bgcflux(2,ind)='SIO2_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
       ind=ind+1
       vname_bgcflux(1,ind)='SIO2_PROD'
       vname_bgcflux(2,ind)='SiO2 production'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='SIO2_REMIN'
       vname_bgcflux(2,ind)='SiO2 remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s'
       ind=ind+1
       vname_bgcflux(1,ind)='DUST_FLUX_IN'
       vname_bgcflux(2,ind)='DUST_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
       ind=ind+1
       vname_bgcflux(1,ind)='DUST_REMIN'
       vname_bgcflux(2,ind)='Dust remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='P_IRON_FLUX_IN'
       vname_bgcflux(2,ind)='P_IRON_FLUX_IN_AVG'
       vname_bgcflux(3,ind)='mmol/m2/s '
       ind=ind+1
       vname_bgcflux(1,ind)='P_IRON_PROD'
       vname_bgcflux(2,ind)='P_IRON_PROD_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='P_IRON_REMIN'
       vname_bgcflux(2,ind)='P_IRON remineralization'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='GRAZE_SP'
       vname_bgcflux(2,ind)='GRAZE_SP_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='GRAZE_DIAT'
       vname_bgcflux(2,ind)='GRAZE_DIAT_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s'
       ind=ind+1
       vname_bgcflux(1,ind)='GRAZE_TOT'
       vname_bgcflux(2,ind)='GRAZE_TOT_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s'
       ind=ind+1
       vname_bgcflux(1,ind)='SP_LOSS'
       vname_bgcflux(2,ind)='SP_LOSS_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_LOSS'
       vname_bgcflux(2,ind)='DIAT_LOSS_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='ZOO_LOSS'
       vname_bgcflux(2,ind)='ZOO_LOSS_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='SP_AGG'
       vname_bgcflux(2,ind)='SP_AGG_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_AGG'
       vname_bgcflux(2,ind)='DIAT_AGG_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOC_SP'
       vname_bgcflux(2,ind)='PHOTOC_SP_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='f_ratio_SP'
       vname_bgcflux(2,ind)='f-ratio for small phytoplankton'
       vname_bgcflux(3,ind)='-'
       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOC_DIAT'
       vname_bgcflux(2,ind)='PHOTOC_DIAT_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='f_ratio_DIAT'
       vname_bgcflux(2,ind)='f-ratio for diatoms'
       vname_bgcflux(3,ind)='-'
       ind=ind+1
       vname_bgcflux(1,ind)='TOT_PROD'
       vname_bgcflux(2,ind)='TOT_PROD_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='no3_v_sp'
       vname_bgcflux(2,ind)='NO3 uptake by SP'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='nh4_v_sp'
       vname_bgcflux(2,ind)='NH4 uptake by SP'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='no3_v_diat'
       vname_bgcflux(2,ind)='NO3 uptake by diatoms'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='nh4_v_diat'
       vname_bgcflux(2,ind)='NH4 uptake by diatoms'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='DOC_PROD'
       vname_bgcflux(2,ind)='DOC_PROD_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='DOC_REMIN'
       vname_bgcflux(2,ind)='DOC_REMIN_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='FE_SCAVENGE'
       vname_bgcflux(2,ind)='FE_SCAVENGE_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='SP_N_LIM'
       vname_bgcflux(2,ind)='SP_N_LIM_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='SP_FE_LIM'
       vname_bgcflux(2,ind)='SP_FE_LIM_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='SP_PO4_LIM'
       vname_bgcflux(2,ind)='SP_PO4_LIM_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='SP_LIGHT_LIM'
       vname_bgcflux(2,ind)='SP_LIGHT_LIM_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_N_LIM'
       vname_bgcflux(2,ind)='DIAT_N_LIM_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_FE_LIM'
       vname_bgcflux(2,ind)='DIAT_FE_LIM_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_PO4_LIM'
       vname_bgcflux(2,ind)='DIAT_PO4_LIM_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_SIO3_LIM'
       vname_bgcflux(2,ind)='DIAT_SIO3_LIM_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAT_LIGHT_LIM'
       vname_bgcflux(2,ind)='DIAT_LIGHT_LIM_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='CACO3_FORM'
       vname_bgcflux(2,ind)='CACO3_FORM_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAZ_NFIX'
       vname_bgcflux(2,ind)='DIAZ_NFIX_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s   '
       ind=ind+1
       vname_bgcflux(1,ind)='GRAZE_DIAZ'
       vname_bgcflux(2,ind)='GRAZE_DIAZ_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAZ_LOSS_AGG'
       vname_bgcflux(2,ind)='DIAZ_LOSS_AGG_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOC_DIAZ'
       vname_bgcflux(2,ind)='PHOTOC_DIAZ_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAZ_P_LIM'
       vname_bgcflux(2,ind)='Diazotroph P limitation'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAZ_FE_LIM'
       vname_bgcflux(2,ind)='Diazotroph Fe limitation'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='DIAZ_LIGHT_LIM'
       vname_bgcflux(2,ind)='Diazotroph light limitation'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='FE_SCAVENGE_RATE'
       vname_bgcflux(2,ind)='Iron scavenging rate'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='DON_PROD'
       vname_bgcflux(2,ind)='DON_PROD_AVG'
       vname_bgcflux(3,ind)=' '
       ind=ind+1
       vname_bgcflux(1,ind)='DON_REMIN'
       vname_bgcflux(2,ind)='DON_REMIN_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='DOFE_PROD'
       vname_bgcflux(2,ind)='DOFE_PROD_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='DOFE_REMIN'
       vname_bgcflux(2,ind)='DOFE_REMIN_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='DOP_PROD'
       vname_bgcflux(2,ind)='DOP_PROD_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='DOP_REMIN'
       vname_bgcflux(2,ind)='DOP_REMIN_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='BSI_FORM'
       vname_bgcflux(2,ind)='BSI_FORM_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s '
       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOFE_DIAZ'
       vname_bgcflux(2,ind)='PHOTOFE_DIAZ_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOFE_DIAT'
       vname_bgcflux(2,ind)='PHOTOFE_DIAT_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='PHOTOFE_SP'
       vname_bgcflux(2,ind)='PHOTOFE_SP_AVG'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='NITRIF'
       vname_bgcflux(2,ind)='Nitrification'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='NH4_oxid'
       vname_bgcflux(2,ind)='NH4 oxidation'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='Denitrif'
       vname_bgcflux(2,ind)='Denitrification'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='N2O_prod'
       vname_bgcflux(2,ind)='N2O production'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       ind=ind+1
       vname_bgcflux(1,ind)='N2_prod'
       vname_bgcflux(2,ind)='N2 production'
       vname_bgcflux(3,ind)='mmol/m3/s  '
       end subroutine init_scalars_becflux
