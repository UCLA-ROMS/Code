      subroutine ecosys_tile(Istr,Iend,Jstr,Jend)
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
        real(kind=8) c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1._8, c0=0.0_8,c2=2._8,
     &  c1000=1000._8,p5=0.5_8,
     &  spd = 86400.0_8,  dps = c1 / spd ,
     &  t0_kelvin= 273.16_8)
        real(kind=8)  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C
        parameter ( parm_Red_D_C_P  = 117.0_8,
     &  parm_Red_P_C_P  = 117.0_8,
     &  parm_Red_D_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_P_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_D_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_P_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_Fe_C   = 3.0D-6)
        real(kind=8), parameter :: Q         = 0.137_8
        real(kind=8), parameter :: Qinv      = 1.0_8/0.137_8
        real(kind=8), parameter :: Qp        = 0.00855_8
        real(kind=8), parameter :: Qp_diat   = 0.1_8 * Q
        real(kind=8), parameter :: Qdenit    = 0.9811320755_8
        real(kind=8)  parm_Fe_bioavail,   parm_prod_dissolve,
     &  parm_o2_min,     parm_Rain_CaCO3, parm_Rain_SiO2,
     &  parm_kappa_nitrif,  parm_nitrif_par_lim,  parm_POC_flux_ref,
     &  parm_rest_prod_tau,  parm_rest_prod_z_c,  parm_z_umax_0,
     &  parm_diat_umax_0,     parm_z_mort_0,     parm_z_mort2_0,
     &  parm_sd_remin_0,     parm_sp_kNO3,       parm_diat_kNO3,
     &  parm_sp_kNH4,        parm_diat_kNH4,     parm_sp_kFe,
     &  parm_diat_kFe,       parm_diat_kSiO3,    parm_sp_kPO4,
     &  parm_diat_kPO4,      parm_z_grz,        parm_alphaChl,
     &  parm_labile_ratio,   parm_alphaDiaz,    parm_diaz_umax_0,
     &  gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
       common/eco_para/parm_Fe_bioavail, parm_prod_dissolve,
     &   parm_o2_min, parm_Rain_CaCO3, parm_Rain_SiO2,
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref,
     &   parm_rest_prod_tau, parm_rest_prod_z_c, parm_z_umax_0,
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0,
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe,
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4,
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl,
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
        real(kind=8) parm_kappa_nitrif_no2
        real(kind=8) parm_n2o_prod
        real(kind=8) parm_n2_prod
        common /eco_param_ncycle_anoxic/ parm_kappa_nitrif_no2,
     &       parm_n2o_prod,parm_n2_prod
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
      integer::i,j,k,m,nn,ctime
      integer(kind=4) Istr,Iend,Jstr,Jend
      real(kind=8),dimension(istr:iend)::
     &     SSTT,
     &     SSSS,
     &     SHF_QSW,
     &     QA_dust_def,
     &     PAR_out
      real(kind=8),dimension(istr:iend,ntrc_bio) :: STF
      real(kind=8),dimension(istr:iend,N)::temp
      real(kind=8),dimension(istr:iend)::
     &     FICE_USED,
     &     WS_USED,
     &     XKW,
     &     AP_USED,
     &     XKW_ICE,
     &     SCHMIDT_USED,
     &     PV,
     &     O2SAT_1atm,
     &     O2SAT_USED,
     &     FLUX
      real(kind=8),dimension(istr:iend)::
     &     XCO2,
     &     PHLO,
     &     PHHI,
     &     PH_NEW,
     &     CO2STAR_ROW,
     &     DCO2STAR_ROW,
     &     pCO2SURF_ROW,
     &     DpCO2_ROW
      real(kind=8), PARAMETER ::
     &     a = 8.6D-7,
     &     phlo_init = 5.0_8,
     &     phhi_init = 9.0_8,
     &     del_ph = 0.25_8
      real(kind=8) get_atm_pco2
      XCO2 = get_atm_pco2(tdays)
      ctime = 1
      do j=jstr,jend
         do m=1,ntrc_bio
            do k=1,N
               tracer(istr:iend,k,m,ctime)=
     &              t(istr:iend,j,k,nnew,1+ntrc_salt+ntrc_pas+m)
            enddo
         enddo
         do m = 1,NT_sed
            tracer_sed(istr:iend,m,ctime) = t_sed(istr:iend,j,m)
         end do
         sstt(istr:iend)=t(istr:iend,j,N,nnew,1)
         ssss(istr:iend)=t(istr:iend,j,N,nnew,2)
         stf=0.0_8
         shf_qsw(istr:iend) =
     &        srflx_dailyavg(istr:iend,j)*rho0*Cp
         do k=1,N
            temp(istr:iend,k)=t(istr:iend,j,k,nnew,1)
         enddo
         do i=istr,iend
            if (.NOT. landmask(i,j))  then
               dust(i,j)=0.0_8
               iron(i,j)=0.0_8
            endif
         enddo
         FICE_USED(istr:iend) = IFRAC(istr:iend,j)
         call WS(sustr(istr:iend,j),
     &        svstr(istr:iend,j),landmask(istr:iend,j),
     &        ws_used(istr:iend), istr,iend)
         XKW = a * WS_USED*WS_USED
         AP_USED = PRESS(istr:iend,j)
         WHERE (AP_USED > 1.5_8 .OR. AP_USED < 0.5_8)
            AP_USED = c1
         END WHERE
         WS_HIST(ISTR:IEND,J) = WS_USED
         XKW_HIST(ISTR:IEND,J) = XKW
         AP_HIST(ISTR:IEND,J) = AP_USED
         IF (lflux_gas_o2 .OR. lflux_gas_co2) THEN
            XKW_ICE = XKW
            WHERE (FICE_USED > 0.2_8 .AND. FICE_USED < 0.9999_8)
               XKW_ICE = (c1 - FICE_USED) * XKW_ICE
            END WHERE
            WHERE (FICE_USED >= 0.9999_8)
               XKW_ICE = c0
            END WHERE
         END IF
         IF (lflux_gas_o2) THEN
            call CSCHMIDT_O2(SSTT(istr:iend),
     &           landmask(istr:iend,j),SCHMIDT_USED,
     &           istr,iend)
            SCHMIDT_O2_HIST(istr:iend,j) = SCHMIDT_USED
            call O2SATU(SSTT(istr:iend),
     &           SSSS(istr:iend),landmask(istr:iend,j),
     &           O2SAT_1atm,istr, iend)
            WHERE (LANDMASK(istr:iend,j))
               PV = XKW_ICE * SQRT(660.0_8 / SCHMIDT_USED)
               O2SAT_USED = AP_USED * O2SAT_1atm
               FLUX = PV * (O2SAT_USED -
     &              tracer(istr:iend,N,o2_ind,ctime))
               STF(istr:iend,o2_ind) = FLUX
               tracer(istr:iend,N,o2_ind,ctime)=
     &              tracer(istr:iend,N,o2_ind,ctime)+
     &              stf(istr:iend,o2_ind)
     &              *dt/Hz(istr:iend,j,N)
               O2SAT_HIST(ISTR:IEND,J) = O2SAT_USED
               FG_O2_HIST(ISTR:IEND,J) = FLUX
            ELSEWHERE
               PV = c0
               O2SAT_HIST(ISTR:IEND,J) = c0
               FG_O2_HIST(ISTR:IEND,J) = c0
            END WHERE
         else
            print *, 'warning: no O2 gas exchange!!!'
         END IF
         IF (lflux_gas_o2) THEN
            WHERE (LANDMASK(istr:iend,j))
               flux = -1.0_8 * pv *
     &              tracer(istr:iend,N,n2o_ind,ctime)
               tracer(istr:iend,N,n2o_ind,ctime)=
     &              tracer(istr:iend,N,n2o_ind,ctime)+
     &              flux * dt / Hz(istr:iend,j,N)
               FG_N2O_HIST(ISTR:IEND,J) = FLUX
               flux = -1.0_8 * pv *
     &              tracer(istr:iend,N,n2_ind,ctime)
               tracer(istr:iend,N,n2_ind,ctime)=
     &              tracer(istr:iend,N,n2_ind,ctime)+
     &              flux * dt / Hz(istr:iend,j,N)
               FG_N2_HIST(ISTR:IEND,J) = FLUX
            ELSEWHERE
               FG_N2O_HIST(ISTR:IEND,J) = c0
               FG_N2_HIST(ISTR:IEND,J) = c0
            END WHERE
         else
            print *, 'warning: no N2O and N2 gas exchange!!!'
         END IF
         IF (lflux_gas_co2) THEN
            call CSCHMIDT_CO2(SSTT(istr:iend),
     &           landmask(istr:iend,j),SCHMIDT_USED,
     &           istr,iend)
            SCHMIDT_CO2_HIST(ISTR:IEND,J) = SCHMIDT_USED
            WHERE (LANDMASK(istr:iend,j))
               PV = XKW_ICE * SQRT(660.0_8 / SCHMIDT_USED)
            ELSEWHERE
               PV = c0
            END WHERE
            WHERE (PH_HIST(istr:iend,j) /= c0)
               PHLO = PH_HIST(istr:iend,j) - del_ph
               PHHI = PH_HIST(istr:iend,j) + del_ph
            ELSEWHERE
               PHLO = phlo_init
               PHHI = phhi_init
            END WHERE
            CALL co2calc_row(LANDMASK(istr,j),
     &           SSTT(istr), SSSS(istr),
     &           tracer(istr,N,dic_ind,ctime),
     &           tracer(istr,N,alk_ind,ctime),
     &           tracer(istr,N,po4_ind,ctime),
     &           tracer(istr,N,sio3_ind,ctime),
     &           PHLO, PHHI, PH_NEW, XCO2,
     &           AP_USED(istr), CO2STAR_ROW,
     &           DCO2STAR_ROW, pCO2SURF_ROW,
     &           DpCO2_ROW,istr,iend)
            pH_hist(istr:iend,j) = PH_NEW
            CO2STAR_HIST(istr:iend,j)  = CO2STAR_ROW
            DCO2STAR_HIST(istr:iend,j) = DCO2STAR_ROW
            pCO2sw(istr:iend,j) = pCO2SURF_ROW
            pCO2air(istr:iend,j)  = pCO2SURF_ROW - DpCO2_ROW
            FLUX(istr:iend) = PV(istr:iend) * DCO2STAR_ROW
            STF(istr:iend,dic_ind) =
     &           STF(istr:iend,dic_ind) + FLUX
            tracer(istr:iend,N,dic_ind,ctime)=
     &           tracer(istr:iend,N,dic_ind,ctime)+
     &           stf(istr:iend,dic_ind)
     &           *dt/Hz(istr:iend,j,N)
            FG_CO2_HIST(istr:iend,j) =  FLUX(istr:iend)
         endif
         if (liron_flux) then
            FLUX(istr:iend)=iron(istr:iend,j)
         else
            FLUX = c0
         endif
         FLUX = FLUX * parm_Fe_bioavail
         STF(istr:iend,fe_ind) = FLUX
         tracer(istr:iend,N,fe_ind,ctime)=
     &        tracer(istr:iend,N,fe_ind,ctime)+
     &        stf(istr:iend,fe_ind)
     &        *dt/Hz(istr:iend,j,N)
         IRON_FLUX_HIST(ISTR:IEND,J) = FLUX
            set_interior: do k=N,1,-1
               call ecosys_set_interior(k,temp(istr:iend,k),
     &           SHF_QSW(istr:iend),
     &           PAR_out(istr:iend),
     &           qa_dust_def(istr:iend),istr,iend,j
     &           ,dt,ctime,dust(istr:iend,j)
     &           )
            end do set_interior
         where (tracer(istr:iend,:,o2_ind,ctime) < 0.0_8)
            tracer(istr:iend,:,o2_ind,ctime) = 0.0_8
         end where
         dom_sp_sfc(istr:iend,j) = 0.0_8
         dom_diat_sfc(istr:iend,j) = 0.0_8
         dom_diaz_sfc(istr:iend,j) = 0.0_8
         dom_sp_int(istr:iend,j) = 0.0_8
         dom_diat_int(istr:iend,j) = 0.0_8
         dom_diaz_int(istr:iend,j) = 0.0_8
           spchl_int(i,j) = c0
           diatchl_int(i,j) = c0
           diazchl_int(i,j) = c0
         do i = istr,iend
           if (tracer(i,N,diatchl_ind,ctime) .ge.
     &           tracer(i,N,diazchl_ind,ctime)) then
              if (tracer(i,N,diatchl_ind,ctime) .ge.
     &             tracer(i,N,spchl_ind,ctime) ) then
                 dom_diat_sfc(i,j) = 1.0_8
              else
                 dom_sp_sfc(i,j) = 1.0_8
              endif
           else
              if (tracer(i,N,spchl_ind,ctime) .ge.
     &             tracer(i,N,diazchl_ind,ctime)) then
                 dom_sp_sfc(i,j) = 1.0_8
              else
                 dom_diaz_sfc(i,j) = 1.0_8
              endif
           endif
           do k = 1, N
              spchl_int(i,j) = spchl_int(i,j) +
     &            tracer(i,k,spchl_ind,ctime) * Hz(i,j,k)
              diatchl_int(i,j) = diatchl_int(i,j) +
     &            tracer(i,k,diatchl_ind,ctime) * Hz(i,j,k)
              diazchl_int(i,j) = diazchl_int(i,j) +
     &            tracer(i,k,diazchl_ind,ctime) * Hz(i,j,k)
           end do
           if (diatchl_int(i,j) .ge. diazchl_int(i,j)) then
              if (diatchl_int(i,j) .ge. spchl_int(i,j)) then
                 dom_diat_int(i,j) = 1.0_8
              else
                 dom_sp_int(i,j) = 1.0_8
              endif
           else
              if (spchl_int(i,j) .ge. diazchl_int(i,j)) then
                 dom_sp_int(i,j) = 1.0_8
              else
                 dom_diaz_int(i,j) = 1.0_8
              endif
           endif
        end do
         do k=1,n
            do m=1,ntrc_bio
               t(istr:iend,j,k,nnew,1+ntrc_salt+ntrc_pas+m)=
     &              tracer(istr:iend,k,m,ctime)
            enddo
         enddo
         do m=1,NT_sed
            t_sed(istr:iend,j,m) = tracer_sed(istr:iend,m,ctime)
         enddo
      end do
      return
      end
        subroutine ecosys_set_interior(k,temp,SHF_QSW,
     &  PAR_out,qa_dust_def,istr,iend,j,dt,ctime,dust_flux
     &     )
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
        real(kind=8) c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1._8, c0=0.0_8,c2=2._8,
     &  c1000=1000._8,p5=0.5_8,
     &  spd = 86400.0_8,  dps = c1 / spd ,
     &  t0_kelvin= 273.16_8)
        real(kind=8)  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C
        parameter ( parm_Red_D_C_P  = 117.0_8,
     &  parm_Red_P_C_P  = 117.0_8,
     &  parm_Red_D_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_P_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_D_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_P_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_Fe_C   = 3.0D-6)
        real(kind=8), parameter :: Q         = 0.137_8
        real(kind=8), parameter :: Qinv      = 1.0_8/0.137_8
        real(kind=8), parameter :: Qp        = 0.00855_8
        real(kind=8), parameter :: Qp_diat   = 0.1_8 * Q
        real(kind=8), parameter :: Qdenit    = 0.9811320755_8
        real(kind=8)  parm_Fe_bioavail,   parm_prod_dissolve,
     &  parm_o2_min,     parm_Rain_CaCO3, parm_Rain_SiO2,
     &  parm_kappa_nitrif,  parm_nitrif_par_lim,  parm_POC_flux_ref,
     &  parm_rest_prod_tau,  parm_rest_prod_z_c,  parm_z_umax_0,
     &  parm_diat_umax_0,     parm_z_mort_0,     parm_z_mort2_0,
     &  parm_sd_remin_0,     parm_sp_kNO3,       parm_diat_kNO3,
     &  parm_sp_kNH4,        parm_diat_kNH4,     parm_sp_kFe,
     &  parm_diat_kFe,       parm_diat_kSiO3,    parm_sp_kPO4,
     &  parm_diat_kPO4,      parm_z_grz,        parm_alphaChl,
     &  parm_labile_ratio,   parm_alphaDiaz,    parm_diaz_umax_0,
     &  gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
       common/eco_para/parm_Fe_bioavail, parm_prod_dissolve,
     &   parm_o2_min, parm_Rain_CaCO3, parm_Rain_SiO2,
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref,
     &   parm_rest_prod_tau, parm_rest_prod_z_c, parm_z_umax_0,
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0,
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe,
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4,
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl,
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
        real(kind=8) parm_kappa_nitrif_no2
        real(kind=8) parm_n2o_prod
        real(kind=8) parm_n2_prod
        common /eco_param_ncycle_anoxic/ parm_kappa_nitrif_no2,
     &       parm_n2o_prod,parm_n2_prod
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
        integer(kind=4) ctime,istr,iend,j,i,m
        real(kind=8) dt
        real(kind=8), DIMENSION(istr:iend) :: temp,SHF_QSW,
     &    QA_dust_def ,PAR_out,dust_flux
     &   , fe_flux
       REAL(kind=8), PARAMETER ::
     &    epsC      = 1.00D-8,
     &    epsTinv   = 3.17D-8,
     &    epsnondim = 1.00D-6,
     &    epsN      = 1.D-1
       REAL(kind=8), DIMENSION(istr:iend) ::
     &    PO4_loc,
     &    NO3_loc,
     &    SiO3_loc,
     &    NH4_loc,
     &    Fe_loc,
     &    O2_loc,
     &    DOC_loc,
     &    spC_loc,
     &    spChl_loc,
     &    spCaCO3_loc,
     &    diatC_loc,
     &    diatChl_loc,
     &    zooC_loc,
     &    spFe_loc,
     &    diatSi_loc,
     &    diatFe_loc,
     &    diazC_loc,
     &    diazChl_loc,
     &    diazFe_loc,
     &    DON_loc,
     &    DOFe_loc,
     &    DOP_loc
     &      ,NO2_loc
     &      ,N2O_loc
     &      ,N2_loc
     &      , log10_Fc, log10_den
       REAL(kind=8) ::
     &    PCref,
     &    sp_mort,
     &    sp_mort2,
     &    diat_mort,
     &    diat_mort2,
     &    z_ingest,
     &    z_grz_sqr,
     &    thres_z1,
     &    thres_z2,
     &    PCrefDiaz,
     &    Qp_diaz,
     &    diaz_mort,
     &    diaz_kPO4,
     &    diaz_kFe,
     &    Qfe_zoo
        REAL(kind=8), DIMENSION(istr:iend) ::
     &    PAR_in,
     &    KPARdz,
     &    PAR_lay,
     &    DOC_prod,
     &    DOC_remin,
     &    NITRIF,
     &    RESTORE
     &    , ammox
     &    , denitr
     &    , denitr_sed
     &    , n2o_prod
     &    , n2_prod
        REAL(kind=8), DIMENSION(istr:iend) ::
     &    z_umax,
     &    diat_umax,
     &    z_mort,
     &    C_loss_diaz,
     &    z_mort2,
     &    diaz_umax
       REAL(kind=8), DIMENSION(istr:iend) ::
     &    thetaC_sp,
     &    thetaC_diat,
     &    QCaCO3,
     &    Tfunc,
     &    VNO3_sp,
     &    VNH4_sp,
     &    VNtot_sp,
     &    VFeC_sp,
     &    VPO4_sp,
     &    f_nut,
     &    PCmax,
     &    PCphoto_sp,
     &    photoC_sp,
     &    f_ratio_sp,
     &    NO3_V_sp,
     &    NH4_V_sp,
     &    VNC_sp,
     &    pChl,
     &    photoacc_sp,
     &    CaCO3_prod,
     &    VNO3_diat,
     &    VNH4_diat,
     &    VNtot_diat,
     &    VFeC_diat,
     &    VPO4_diat,
     &    VSiO3_diat,
     &    PCphoto_diat,
     &    photoC_diat,
     &    f_ratio_diat,
     &    NO3_V_diat,
     &    NH4_V_diat,
     &    VNC_diat,
     &    photoacc_diat,
     &    reduceV,
     &    graze_sp,
     &    graze_sp_zoo,
     &    graze_sp_poc,
     &    graze_sp_doc,
     &    graze_sp_dic
        REAL(kind=8), DIMENSION(istr:iend) ::
     &    graze_diat,
     &    graze_diat_zoo,
     &    graze_diat_poc,
     &    graze_diat_doc,
     &    graze_diat_dic,
     &    Pprime,
     &    sp_loss,
     &    sp_loss_poc,
     &    sp_loss_doc,
     &    sp_loss_dic,
     &    sp_agg,
     &    diat_loss,
     &    diat_loss_poc,
     &    diat_loss_doc,
     &    diat_loss_dic,
     &    diat_agg,
     &    f_zoo_detr,
     &    Fe_scavenge,
     &    Zprime,
     &    zoo_loss,
     &    zoo_loss_doc,
     &    zoo_loss_dic,
     &    WORK,
     &    light_lim,
     &    Qsi,
     &    gQsi,
     &    Qfe_sp,
     &    gQfe_sp,
     &    Qfe_diat,
     &    gQfe_diat,
     &    Qfe_diaz,
     &    gQfe_diaz
       REAL(kind=8), DIMENSION(istr:iend) ::
     &    PCphoto_diaz,
     &    photoC_diaz,
     &    Vfec_diaz,
     &    Vpo4_diaz,
     &    photoacc_diaz,
     &    Vnc_diaz,
     &    diaz_Nexcrete,
     &    photoN_diaz,
     &    diaz_Nfix,
     &    thetaC_diaz,
     &    photoFe_diaz,
     &    photoFe_diat,
     &    photoFe_sp,
     &    photoSi_diat,
     &    remaining_diazP,
     &    diaz_agg,
     &    diaz_loss,
     &    diaz_loss_doc,
     &    diaz_loss_dic,
     &    diaz_loss_dop,
     &    diaz_loss_dip,
     &   graze_diaz,
     &   graze_diaz_poc,
     &   graze_diaz_doc,
     &   graze_diaz_dic,
     &   graze_diaz_zoo,
     &   DON_remin,
     &   DOFe_remin,
     &   DOP_remin,
     &   DOM_remin,
     &   Fe_scavenge_rate,
     &   fe_max_scale1,
     &   fe_scavenge_thres2,
     &   CaCO3_temp_thres1,
     &   CaCO3_temp_thres2,
     &   f_prod_sp_CaCO3,
     &   spc_poc_fac
         REAL(kind=8), DIMENSION(istr:iend) ::
     &   DON_prod,
     &   DOFe_prod,
     &   DOP_prod,
     &    C_loss_thres
       REAL(kind=8), DIMENSION(istr:iend) :: Sed_POC, Sed_CaCO3, Sed_Si
       REAL(kind=8), DIMENSION(istr:iend) :: bot_flux_poc, 
     &                          bot_flux_caco3,
     &      bot_flux_si
       REAL(kind=8), DIMENSION(istr:iend) :: remin_sed_poc, 
     &                          remin_sed_caco3,
     &      remin_sed_si
         integer(kind=4) setinterior,curtime,k
          setinterior=1
         IF (.NOT. lsource_sink) RETURN
         curtime=2
       PO4_loc  = MAX(c0, TRACER(istr:iend,k,po4_ind,ctime))
       NO3_loc  = MAX(c0, TRACER(istr:iend,k,no3_ind,ctime))
       SiO3_loc = MAX(c0, TRACER(istr:iend,k,sio3_ind,ctime))
       NH4_loc  = MAX(c0, TRACER(istr:iend,k,nh4_ind,ctime))
       Fe_loc   = MAX(c0, TRACER(istr:iend,k,fe_ind,ctime))
       O2_loc   = MAX(c0, TRACER(istr:iend,k,o2_ind,ctime))
       DOC_loc  = MAX(c0, TRACER(istr:iend,k,doc_ind,ctime))
       spC_loc  = MAX(c0, TRACER(istr:iend,k,spC_ind,ctime))
       spChl_loc= MAX(c0, TRACER(istr:iend,k,spChl_ind,ctime))
       spCaCO3_loc  =
     &           MAX(c0,TRACER(istr:iend,k,spCaCO3_ind,ctime))
       diatC_loc= MAX(c0, TRACER(istr:iend,k,diatC_ind,ctime))
       diatChl_loc  =
     &          MAX(c0, TRACER(istr:iend,k,diatChl_ind,ctime))
       zooC_loc  = MAX(c0, TRACER(istr:iend,k,zooC_ind,ctime))
       spFe_loc  = MAX(c0, TRACER(istr:iend,k,spFe_ind,ctime))
       diatSi_loc =
     &          MAX(c0, TRACER(istr:iend,k,diatSi_ind,ctime))
       diatFe_loc  =
     &         MAX(c0, TRACER(istr:iend,k,diatFe_ind,ctime))
       diazC_loc    =
     &         MAX(c0, TRACER(istr:iend,k,diazC_ind,ctime))
       diazChl_loc  =
     &         MAX(c0, TRACER(istr:iend,k,diazChl_ind,ctime))
       diazFe_loc   =
     &       MAX(c0, TRACER(istr:iend,k,diazFe_ind,ctime))
       DON_loc  = MAX(c0, TRACER(istr:iend,k,don_ind,ctime))
       DOFe_loc = MAX(c0, TRACER(istr:iend,k,dofe_ind,ctime))
       DOP_loc  = MAX(c0, TRACER(istr:iend,k,dop_ind,ctime))
       NO2_loc = MAX(c0, TRACER(istr:iend,k,no2_ind,ctime))
       N2O_loc = MAX(c0, TRACER(istr:iend,k,n2o_ind,ctime))
       N2_loc = MAX(c0, TRACER(istr:iend,k,n2_ind,ctime))
       WHERE (.NOT. LANDMASK(istr:iend,j) )
          PO4_loc      = c0
          NO3_loc      = c0
          SiO3_loc     = c0
          NH4_loc      = c0
          Fe_loc       = c0
          O2_loc       = c0
          DOC_loc      = c0
          spC_loc      = c0
          spChl_loc    = c0
          spCaCO3_loc  = c0
          diatC_loc    = c0
          diatChl_loc  = c0
          zooC_loc     = c0
          spFe_loc     = c0
          diatSi_loc   = c0
          diatFe_loc   = c0
          diazC_loc    = c0
          diazChl_loc  = c0
          diazFe_loc   = c0
          DON_loc      = c0
          DOFe_loc     = c0
          DOP_loc      = c0
          NO2_loc      = c0
          N2O_loc      = c0
          N2_loc       = c0
        END WHERE
        remin_sed_poc = c0
        remin_sed_caco3 = c0
        remin_sed_si = c0
        if (k .eq. 1) then
           Sed_POC = tracer_sed(istr:iend,sed_poc_ind,ctime)
           Sed_CaCO3 = tracer_sed(istr:iend,sed_caco3_ind,ctime)
           Sed_Si = tracer_sed(istr:iend,sed_si_ind,ctime)
        else
           Sed_POC = c0
           Sed_CaCO3 = c0
           Sed_Si = c0
        end if
       WHERE (spC_loc == c0 .OR. spChl_loc == c0 .OR. spFe_loc == c0)
         spC_loc = c0
         spChl_loc = c0
         spCaCO3_loc = c0
         spFe_loc = c0
       END WHERE
       WHERE (diatC_loc == c0 .OR. diatChl_loc == c0 .OR.
     &     diatFe_loc == c0
     &    .OR. diatSi_loc == c0)
          diatC_loc = c0
          diatChl_loc = c0
          diatFe_loc = c0
          diatSi_loc = c0
        END WHERE
       WHERE (diazC_loc == c0 .OR. diazChl_loc == c0 .OR.
     &   diazFe_loc == c0)
         diazC_loc = c0
         diazChl_loc = c0
         diazFe_loc = c0
        END WHERE
         PCref      = 3.0_8 * dps
         sp_mort    = 0.15_8 * dps
         sp_mort2   = 0.0035_8 * dps
         diat_mort  = 0.15_8 * dps
         diat_mort2 = 0.0035_8 * dps
         z_ingest   = 0.3_8
         thres_z1   = 100.0_8
         thres_z2   = 200.0_8
         PCrefDiaz  = 0.4_8 * dps
         diaz_mort  = 0.16_8 * dps
         diaz_kPO4  = 0.005_8
         diaz_kFe   = 0.1D-3
         Qp_diaz   = 0.002735_8
         Qfe_zoo   = 2.5D-6
         thetaC_sp   = spChl_loc / (spC_loc + epsC)
         thetaC_diat = diatChl_loc / (diatC_loc + epsC)
         thetaC_diaz = diazChl_loc / (diazC_loc + epsC)
         Qsi         = diatSi_loc / (diatC_loc + epsC)
         Qfe_diat    = diatFe_loc / (diatC_loc + epsC)
         Qfe_sp      = spFe_loc / (spC_loc + epsC)
        Qfe_diaz     = diazFe_loc / (diazC_loc + epsC)
        WHERE (Qsi > 0.685_8) Qsi = 0.685_8
        gQsi = gQsi_0
        gQfe_diat = 6.0D-6
        gQfe_sp = 6.0D-6
        gQfe_diaz = 42.0D-6
      WHERE (Fe_loc .LT. c2 * parm_diat_kfe)
        gQfe_diat = MAX((gQfe_diat * Fe_loc /(c2 * parm_diat_kfe)),
     &    2.5D-6)
       END WHERE
       WHERE ((Fe_loc .LT. c2 * parm_diat_kfe).AND. fe_loc .gt. c0
     &   .and. (SiO3_loc.GT. (c2 * parm_diat_kSiO3) ) )
         gQsi = MIN(((gQsi*gQsi_coef*c2*parm_diat_kfe/Fe_loc)
     &     - (gQsi_coef-c1)*gQsi_0), 0.685_8)
       END WHERE
       where (Fe_loc == c0)
          gQsi = gQsi_max
       endwhere
       WHERE (Fe_loc .LT. c2 * parm_sp_kfe)
         gQfe_sp = MAX((gQfe_sp*Fe_loc/(c2 * parm_sp_kfe)),
     &     2.5D-6)
        END WHERE
       WHERE (Fe_loc .LT. c2 * diaz_kFe)
           gQfe_diaz = MAX((gQfe_diaz*Fe_loc/(c2 * diaz_kFe)),
     &                  14.0D-6)
         END WHERE
        WHERE (SiO3_loc .LT. (c2 * parm_diat_kSiO3))
         gQsi = MAX((gQsi*SiO3_loc/(c2 * parm_diat_kSiO3)),
     &                  0.0685_8)
        END WHERE
        IF (k == N) THEN
           WHERE (LANDMASK(istr:iend,j))
              PAR_out = MAX(c0, 0.45_8 * SHF_QSW)
           ELSE WHERE
              PAR_out = c0
           END WHERE
           PARinc(istr:iend,j) = PAR_out(istr:iend)
           CALL init_particulate_terms(QA_dust_def,istr,iend,j,
     &          dust_flux)
        END IF
        QCaCO3 = 0.4_8
        where (spCaCO3_loc .gt. c0 .and. spc_loc .gt. c0)
           QCaCO3 = spCaCO3_loc / (spC_loc + epsC)
           WHERE (QCaCO3 > 0.4_8) QCaCO3 = 0.4_8
        else where
              QCaCO3 = 0.01_8
        end where
         PAR_in = PAR_out
         WHERE (.NOT. LANDMASK(ISTR:IEND,J))
     &        PAR_in = c0
        KPARdz = (0.03_8 * (spChl_loc + diatChl_loc +
     &          diazChl_loc) + 0.04_8) * Hz(istr:iend,j,k)
        where (KPARdz .gt. c0)
           PAR_out = PAR_in * EXP(-KPARdz)
           PAR_lay = PAR_in * (c1 - EXP(-KPARdz)) / KPARdz
        elsewhere
           PAR_out = PAR_in
           PAR_lay = PAR_in
        end where
         Tfunc = 2.0_8**(0.1_8 * temp - 3.0_8)
        z_umax = parm_z_umax_0 * Tfunc
        diat_umax = parm_diat_umax_0 * Tfunc
        z_mort2 = parm_z_mort2_0 * Tfunc
        z_mort = parm_z_mort_0 * Tfunc
        diaz_umax = parm_diaz_umax_0 * Tfunc
         DOM_remin= parm_sd_remin_0
        VNO3_sp = (NO3_loc / parm_sp_kNO3) /
     &    (c1 + (NO3_loc / parm_sp_kNO3) + (NH4_loc / parm_sp_kNH4))
        VNH4_sp = (NH4_loc / parm_sp_kNH4) /
     &    (c1 + (NO3_loc / parm_sp_kNO3) + (NH4_loc / parm_sp_kNH4))
        VNtot_sp = VNO3_sp + VNH4_sp
        sp_N_lim_HIST(ISTR:IEND,J,K) = VNtot_sp
           VFeC_sp = Fe_loc / (Fe_loc + parm_sp_kFe)
           VPO4_sp = PO4_loc / (PO4_loc + parm_sp_kPO4)
           sp_Fe_lim_HIST(ISTR:IEND,J,K) = VFeC_sp
           sp_PO4_lim_HIST(ISTR:IEND,J,K) = VPO4_sp
           f_nut = MIN(VNtot_sp, VFeC_sp)
           f_nut = MIN(f_nut, VPO4_sp)
           PCmax = PCref * f_nut * Tfunc
        light_lim = (c1 - EXP((-c1 * parm_alphaChl * thetaC_sp *
     &   PAR_lay) / (PCmax + epsTinv)))
        PCphoto_sp = PCmax * light_lim
        sp_light_lim_HIST(ISTR:IEND,J,K) = light_lim
        photoC_sp = PCphoto_sp * spC_loc
       WHERE (VNtot_sp > c0)
          NO3_V_sp = (VNO3_sp / VNtot_sp) * photoC_sp * Q
          NH4_V_sp = (VNH4_sp / VNtot_sp) * photoC_sp * Q
          VNC_sp = PCphoto_sp * Q
          f_ratio_sp = VNO3_sp / VNtot_sp
       ELSEWHERE
          NO3_V_sp = c0
          NH4_V_sp = c0
          VNC_sp = c0
          f_ratio_sp = c0
       END WHERE
       where (dt*no3_v_sp > no3_loc - 0.001_8)
          no3_v_sp = c0
       end where
       where (dt*nh4_v_sp > nh4_loc - 0.001_8)
          nh4_v_sp = c0
       end where
       photoC_sp = (no3_v_sp + nh4_v_sp)
     &      * Qinv
       where (spC_loc .gt. c0)
          PCphoto_sp = photoC_sp / spC_loc
       else where
          PCphoto_sp = c0
       end where
       VNC_sp = PCphoto_sp * Q
        photoFe_sp = photoC_sp * gQfe_sp
        photoFe_sp_HIST(ISTR:IEND,J,K) = photoFe_sp
         WORK = parm_alphaChl * thetaC_sp * PAR_lay
       WHERE (WORK > c0)
        pChl = 2.5_8 * PCphoto_sp / WORK
        photoacc_sp = (pChl * VNC_sp / thetaC_sp) * spChl_loc
       ELSEWHERE
        photoacc_sp = c0
       END WHERE
        f_prod_sp_CaCO3 = 0.026_8
        CaCO3_prod = f_prod_sp_CaCO3 * photoC_sp
        CaCO3_prod = CaCO3_prod * f_nut
        CaCO3_temp_thres1 = 1.0_8
        CaCO3_temp_thres2 = -2.0_8
       WHERE (temp < CaCO3_temp_thres1)
     &      CaCO3_prod = CaCO3_prod *
     &       MAX((temp - CaCO3_temp_thres2), c0) /
     &       (CaCO3_temp_thres1 - CaCO3_temp_thres2)
        WHERE (spC_loc > 3.0_8)
     &      CaCO3_prod = MIN((CaCO3_prod*spC_loc
     &                     /3.0_8),(0.4_8*photoC_sp))
        CaCO3_form_HIST(ISTR:IEND,J,K) = CaCO3_prod
        VNO3_diat = (NO3_loc / parm_diat_kNO3) /
     &    (c1 + (NO3_loc / parm_diat_kNO3) + (NH4_loc / parm_diat_kNH4))
        VNH4_diat = (NH4_loc / parm_diat_kNH4) /
     &    (c1 + (NO3_loc / parm_diat_kNO3) + (NH4_loc / parm_diat_kNH4))
        VNtot_diat = VNO3_diat + VNH4_diat
        diat_N_lim_HIST(ISTR:IEND,J,K) = VNtot_diat
          VFeC_diat = Fe_loc / (Fe_loc + parm_diat_kFe)
          VPO4_diat = PO4_loc / (PO4_loc + parm_diat_kPO4)
          VSiO3_diat = SiO3_loc / (SiO3_loc + parm_diat_kSiO3)
          diat_Fe_lim_HIST(ISTR:IEND,J,K) = VFeC_diat
          diat_PO4_lim_HIST(ISTR:IEND,J,K) = VPO4_diat
          diat_SiO3_lim_HIST(ISTR:IEND,J,K) = VSiO3_diat
          f_nut = MIN(VNtot_diat, VFeC_diat)
          f_nut = MIN(f_nut, VSiO3_diat)
          f_nut = MIN(f_nut, VPO4_diat)
         PCmax = PCref * f_nut * Tfunc
        light_lim = (c1 - EXP((-c1 * parm_alphaChl * thetaC_diat
     &    * PAR_lay) /
     &              (PCmax + epsTinv)))
         PCphoto_diat = PCmax * light_lim
         diat_light_lim_HIST(ISTR:IEND,J,K) = light_lim
         photoC_diat = PCphoto_diat * diatC_loc
        WHERE (VNtot_diat > c0)
            NO3_V_diat = (VNO3_diat / VNtot_diat) * photoC_diat * Q
            NH4_V_diat = (VNH4_diat / VNtot_diat) * photoC_diat * Q
            VNC_diat = PCphoto_diat * Q
            f_ratio_diat = VNO3_diat / VNtot_diat
        ELSEWHERE
           NO3_V_diat = c0
           NH4_V_diat = c0
           VNC_diat = c0
           f_ratio_diat = c0
         END WHERE
       where (dt*(no3_v_sp+no3_v_diat) > no3_loc - 0.001_8)
          no3_v_diat = c0
       end where
       where (dt*(nh4_v_sp+nh4_v_diat) > nh4_loc - 0.001_8)
          nh4_v_diat = c0
       end where
       photoC_diat = (no3_v_diat + nh4_v_diat)
     &      * Qinv
       where (diatC_loc .gt. c0)
          PCphoto_diat = photoC_diat / diatC_loc
       else where
          PCphoto_diat = c0
       end where
       VNC_diat = PCphoto_diat * Q
         photoFe_diat = photoC_diat * gQfe_diat
         photoSi_diat = photoC_diat * gQsi
          photoFe_diat_HIST(ISTR:IEND,J,K) = photoFe_diat
          bSi_form_HIST(ISTR:IEND,J,K) = photoSi_diat
         WORK = parm_alphaChl * thetaC_diat * PAR_lay
       WHERE (WORK > c0)
           pChl = 4.0_8 * PCphoto_diat / WORK
          photoacc_diat = (pChl * VNC_diat / thetaC_diat) * diatChl_loc
        ELSEWHERE
       photoacc_diat = c0
       END WHERE
        Vfec_diaz = Fe_loc/(Fe_loc + diaz_kFe)
        Vpo4_diaz = PO4_loc / (PO4_loc + diaz_kPO4)
        diaz_Fe_lim_HIST(ISTR:IEND,J,K) = Vfec_diaz
        diaz_P_lim_HIST(ISTR:IEND,J,K) = Vpo4_diaz
        f_nut = MIN(Vpo4_diaz, Vfec_diaz)
        PCmax = PCrefDiaz * f_nut * Tfunc
        light_lim = (c1 - EXP((-c1 * parm_alphaDiaz * thetaC_diaz
     &    * PAR_lay) / (PCmax + epsTinv)))
          PCphoto_diaz = PCmax * light_lim
          diaz_light_lim_HIST(ISTR:IEND,J,K) = light_lim
         photoC_diaz = PCphoto_diaz * diazC_loc
         diaz_Nfix = photoC_diaz * Q * 1.42857142857143_8
         diaz_Nexcrete = diaz_Nfix * 0.3_8
         photoN_diaz   = diaz_Nfix - diaz_Nexcrete
         diaz_Nfix_HIST(ISTR:IEND,J,K) = diaz_Nfix
         Vnc_diaz = PCphoto_diaz * Q
         photoFe_diaz = photoC_diaz * gQfe_diaz
         photoFe_diaz_HIST(ISTR:IEND,J,K) = photoFe_diaz
         WORK = parm_alphaDiaz * thetaC_diaz * PAR_lay
        WHERE (WORK > c0)
            pChl = 3.4_8 * PCphoto_diaz / WORK
          photoacc_diaz = (pChl * Vnc_diaz / thetaC_diaz) * diazChl_loc
        ELSEWHERE
           photoacc_diaz = c0
         END WHERE
         C_loss_thres = 0.001_8
         WHERE (-z_r(istr:iend,j,k) > thres_z1)
            WHERE (-z_r(istr:iend,j,k) < thres_z2)
              C_loss_thres = C_loss_thres *
     &          (thres_z2 + z_r(istr:iend,j,k)) /
     &               (thres_z2 - thres_z1)
             ELSE WHERE
              C_loss_thres = c0
             END WHERE
          END WHERE
             Pprime = MAX(spC_loc - C_loss_thres, c0)
           sp_loss = sp_mort * Pprime
          sp_agg = MIN((0.75_8 * dps) * Pprime,
     &       sp_mort2 * Pprime * Pprime)
         reduceV = Pprime * Pprime
         z_grz_sqr = parm_z_grz * parm_z_grz
          graze_sp = z_umax * zooC_loc * (reduceV /
     &            (reduceV + z_grz_sqr))
        graze_sp_zoo = z_ingest * graze_sp
        spc_poc_fac = 0.22_8
        graze_sp_poc = graze_sp * MAX((0.4_8 * QCaCO3),
     &          MIN((0.18_8 * Pprime),spc_poc_fac))
         graze_sp_doc = 0.34_8 * graze_sp - graze_sp_poc
         graze_sp_dic = 0.36_8 * graze_sp
          sp_loss_poc = QCaCO3 * sp_loss
        sp_loss_doc = (c1 - parm_labile_ratio) *
     &                 (sp_loss - sp_loss_poc)
         sp_loss_dic = parm_labile_ratio * (sp_loss - sp_loss_poc)
          C_loss_thres = 0.02_8
         WHERE (-z_r(istr:iend,j,k) > thres_z1)
            WHERE (-z_r(istr:iend,j,k) < thres_z2)
              C_loss_thres = C_loss_thres *
     &            (thres_z2 + z_r(istr:iend,j,k)) /
     &               (thres_z2 - thres_z1)
             ELSE WHERE
              C_loss_thres = c0
           END WHERE
           END WHERE
          Pprime = MAX(diatC_loc - C_loss_thres, c0)
             diat_loss = diat_mort * Pprime
            diat_agg = MIN((0.75_8 * dps) * Pprime,
     &               diat_mort2 * Pprime * Pprime)
           diat_agg = MAX((0.01_8 * dps) * Pprime, diat_agg)
          reduceV = Pprime * Pprime
          graze_diat = diat_umax *zooC_loc *
     &       (reduceV / (reduceV + z_grz_sqr * 0.81_8))
          graze_diat_zoo = z_ingest * graze_diat
         graze_diat_poc = 0.26_8 * graze_diat
         graze_diat_doc = 0.13_8 * graze_diat
         graze_diat_dic = 0.31_8 * graze_diat
         diat_loss_poc = 0.05_8 * diat_loss
         diat_loss_doc = (c1 - parm_labile_ratio) * 0.95_8 * diat_loss
         diat_loss_dic = parm_labile_ratio * 0.95_8 * diat_loss
            C_loss_diaz = 0.01_8
         WHERE (temp .LT. 15.0_8)
     &           C_loss_diaz = 0.001_8
         WHERE (-z_r(istr:iend,j,k) > thres_z1)
            WHERE (-z_r(istr:iend,j,k) < thres_z2)
             C_loss_diaz = C_loss_diaz *
     &             (thres_z2 + z_r(istr:iend,j,k))/
     &             (thres_z2-thres_z1)
             ELSE WHERE
              C_loss_diaz = c0
           END WHERE
           END WHERE
         Pprime = MAX(diazC_loc - C_loss_diaz, c0)
         diaz_loss = diaz_mort * Pprime
         diaz_agg = 0.0_8
         reduceV = Pprime * Pprime
         graze_diaz = diaz_umax * zooC_loc *
     &        (reduceV / (reduceV + z_grz_sqr))
         graze_diaz_zoo = 0.21_8 * graze_diaz
         graze_diaz_poc = 0.0_8
         graze_diaz_doc = 0.24_8 * graze_diaz
         graze_diaz_dic = 0.55_8 * graze_diaz
        diaz_loss_doc = (c1 - parm_labile_ratio) * diaz_loss
        diaz_loss_dic = parm_labile_ratio * diaz_loss
        remaining_diazP =  ((graze_diaz + diaz_loss + diaz_agg)
     &                                                  * Qp_diaz)-
     &               ((graze_diaz_poc+graze_diaz_zoo) * Qp)
        diaz_loss_dop = (c1 - parm_labile_ratio) * remaining_diazP
        diaz_loss_dip = parm_labile_ratio * remaining_diazP
        f_zoo_detr = (0.1333_8 * (graze_diat + epsC * epsTinv) +
     &    0.0333_8 * (graze_sp + epsC * epsTinv) +
     &    0.0_8 * (graze_diaz + epsC * epsTinv)) /
     &    (graze_diat + graze_sp + graze_diaz + 3.0_8 * epsC * epsTinv)
         C_loss_thres = 0.03_8
         WHERE (-z_r(istr:iend,j,k) > thres_z1)
            WHERE (-z_r(istr:iend,j,k) < thres_z2)
              C_loss_thres = C_loss_thres *
     &            ((-z_r(istr:iend,j,k)-thres_z1)/
     &                 (thres_z2-thres_z1))
             ELSE WHERE
              C_loss_thres = c0
           END WHERE
           END WHERE
        Zprime = MAX(zooC_loc - C_loss_thres, c0)
        zoo_loss = z_mort2 * Zprime * Zprime + z_mort * Zprime
        zoo_loss_doc = (c1 - parm_labile_ratio) * (c1 - f_zoo_detr)
     &                     * zoo_loss
         zoo_loss_dic = parm_labile_ratio * (c1 - f_zoo_detr)
     &                       * zoo_loss
         DOC_prod = sp_loss_doc + graze_sp_doc + zoo_loss_doc
     &       + diat_loss_doc
     &       + graze_diat_doc + diaz_loss_doc + graze_diaz_doc
         DON_prod = (DOC_prod * Q)
     &        + diaz_Nexcrete
         DOP_prod = (sp_loss_doc + graze_sp_doc + zoo_loss_doc
     &       + diat_loss_doc
     &          + graze_diat_doc) * Qp + diaz_loss_dop
         DOFe_prod = (zoo_loss_doc * Qfe_zoo)
     &  + (Qfe_sp * (graze_sp_doc + sp_loss_doc))
     &  + (Qfe_diat * (graze_diat_doc + diat_loss_doc))
     &  + (Qfe_diaz * (graze_diaz_doc + diaz_loss_doc))
        DOC_remin = DOC_loc * DOM_remin
        DON_remin = DON_loc * DOM_remin
        DOFe_remin = DOFe_loc * DOM_remin
        DOP_remin = DOP_loc * DOM_remin
        sinking_particle_POC(3,istr:iend,j)  =
     &      sp_agg + graze_sp_poc
     &      + sp_loss_poc + f_zoo_detr * zoo_loss +
     &   diat_loss_poc + diat_agg + graze_diat_poc + graze_diaz_poc
     &    + diaz_agg
         sinking_particle_P_CaCO3(3,istr:iend,j)=
     &        (0.67_8 * graze_sp + sp_loss + sp_agg) * QCaCO3
        sinking_particle_P_SiO2(3,istr:iend,j) =
     &      (0.5_8 * graze_diat
     &      + diat_agg + 0.05_8 * diat_loss) * Qsi
        sinking_particle_dust(3,istr:iend,j) = c0
         Fe_scavenge_rate = 0.12_8
         fe_max_scale1 = 3.0_8
         Fe_scavenge_rate = Fe_scavenge_rate
     &    *MIN(((sinking_particle_POC(4,istr:iend,j) +
     &     sinking_particle_POC(5,istr:iend,j)+
     &       ((sinking_particle_dust(4,istr:iend,j)
     &    + sinking_particle_dust(5,istr:iend,j)) * 8.33D4) )
     &     / parm_POC_flux_ref), fe_max_scale1)
         WHERE (Fe_loc > 0.6D-3)
     &      Fe_scavenge_rate = Fe_scavenge_rate + (Fe_loc - 0.6D-3)
     &          * (6.0_8 / (1.4D-3))
         fe_scavenge_thres2 = 0.5D-3
        WHERE (Fe_loc < fe_scavenge_thres2)
     &        Fe_scavenge_rate = Fe_scavenge_rate *
     &        (Fe_loc / fe_scavenge_thres2)
          Fe_scavenge = 3.1709792D-8 * Fe_loc * Fe_scavenge_rate
         sinking_particle_P_iron(3,istr:iend,j) =
     &          ((sp_agg + graze_sp_poc
     &       + sp_loss_poc) * Qfe_sp)
     &      + (zoo_loss * f_zoo_detr * Qfe_zoo)
     &    + ((diat_agg + graze_diat_poc + diat_loss_poc) * Qfe_diat)
     &  + ((graze_diaz_poc + diaz_agg) * Qfe_diaz) + (0.1_8 * 
     &                            Fe_scavenge)
        CALL compute_particulate_terms(k, QA_dust_def, temp,istr,
     &               iend,j, bot_flux_poc
     &        , bot_flux_caco3, bot_flux_si
     &        )
        if (k .eq. 1) then
           remin_sed_poc = Sed_POC * t_remin_sed_poc
           remin_sed_caco3 = Sed_CaCO3 * t_remin_sed_caco3
           remin_sed_si = Sed_Si * t_remin_sed_si
        end if
         IF (lrest_no3) THEN
           RESTORE = (NO3_CLIM(istr:iend,j,k) - NO3_loc)
     &          * nutr_rest_time_inv(k)
         ELSE
           RESTORE = c0
         END IF
         NO3_RESTORE_HIST(ISTR:IEND,J,K) = RESTORE
         nitrif = c0
         denitr = c0
         n2o_prod = c0
         n2_prod = c0
         ammox = c0
         denitr_sed = c0
         where (nh4_loc .gt. c0 .and. PAR_out < parm_nitrif_par_lim)
            ammox = parm_kappa_nitrif * NH4_loc
         end where
         where (O2_loc > parm_o2_min)
            where (no2_loc > c0 .and. PAR_out < parm_nitrif_par_lim)
               nitrif = parm_kappa_nitrif_no2 * NO2_loc
            end where
         else where
            where (no3_loc > dt*Sinking_Particle_Poc(6,istr:iend,j)
     &           * Qdenit + epsN)
               denitr = Sinking_Particle_Poc(6,istr:iend,j) * Qdenit
               where (no2_loc > dt*(parm_n2o_prod*denitr + nitrif)
     &              + epsN)
                  n2o_prod = parm_n2o_prod * denitr
               end where
            end where
            where (n2o_loc > dt*parm_n2_prod * N2O_loc + epsN)
               n2_prod = parm_n2_prod * N2O_loc
            end where
         end where
         if (k .eq. 1) then
            where (bot_flux_poc .gt. 0.0_8)
               log10_Fc = log10(8640.0_8 * bot_flux_poc )
               log10_den = -0.9543_8 + 0.7662_8 * log10_Fc -
     &              0.235_8 * log10_Fc * log10_Fc
               denitr_sed = 9.25926D-5 * (10.0_8 ** log10_Den) /
     &              Hz(istr:iend,j,1)
               where (no3_loc < dt*denitr_sed + epsN)
                  denitr_sed = c0
               end where
            else where
               denitr_sed = c0
            end where
         end if
         TRACER(istr:iend,K,no2_ind,curtime) = ammox - nitrif
     &        - n2o_prod + denitr
         TRACER(istr:iend,K,n2o_ind,curtime) = 0.5_8 * n2o_prod -n2_prod
         TRACER(istr:iend,K,n2_ind,curtime) = n2_prod + 0.5_8 * 
     &                             denitr_sed
         NITRIF_HIST(ISTR:IEND,J,K) = NITRIF
         TRACER(istr:iend,K,no3_ind,curtime) = RESTORE + NITRIF-
     &                  (NO3_V_diat + NO3_V_sp
     &        + denitr + denitr_sed
     &        )
         TRACER(istr:iend,K,nh4_ind,curtime) =
     &    - ( NH4_V_diat + NH4_V_sp
     &    ) + Q * ( zoo_loss_dic +
     &    sp_loss_dic +
     &     graze_sp_dic +
     &     diat_loss_dic +
     &    graze_diat_dic +
     &    Sinking_Particle_Poc(6,istr:iend,j) +
     &    diaz_loss_dic +
     &    graze_diaz_dic ) +
     &    DON_remin
     &        - ammox
     &       + q*remin_sed_poc / Hz(istr:iend,j,1)
       TRACER(istr:iend,K,fe_ind,curtime) =
     &    sinking_particle_P_iron(6,istr:iend,j)
     &    - Fe_scavenge +
     &    (Qfe_zoo * zoo_loss_dic) + DOFe_remin - photoFe_diaz
     &     + (Qfe_sp * (sp_loss_dic + graze_sp_dic))
     &     + (Qfe_diat * (diat_loss_dic + graze_diat_dic))
     &     + (Qfe_diaz * (diaz_loss_dic + graze_diaz_dic))
     &     - photoFe_sp - photoFe_diat
     &    + graze_diaz_zoo *(Qfe_diaz-Qfe_zoo)
     &     + graze_diat_zoo *(Qfe_diat-Qfe_zoo)
     &     + graze_sp_zoo * (Qfe_sp-Qfe_zoo)
       if (k .eq. 1) then
          where (TRACER(istr:iend,k,o2_ind,ctime) .gt. 0.0_8)
             fe_flux = 10.0_8**(2.5_8 -
     &            0.0165_8*TRACER(istr:iend,k,o2_ind,ctime))
     &            *0.001_8/86400.0_8/Hz(istr:iend,j,k)
          elsewhere (landmask(istr:iend,j))
             fe_flux = 10.0_8**(2.5_8 - 0.0165_8)
     &            *0.001_8/86400.0_8/Hz(istr:iend,j,k)
          elsewhere
             fe_flux = 0.0_8
          end where
          TRACER(istr:iend,K,fe_ind,curtime) =
     &         TRACER(istr:iend,K,fe_ind,curtime) + fe_flux
     &
       end if
        IF (lrest_sio3) THEN
          RESTORE = (SiO3_CLIM(istr:iend,j,k) - SiO3_loc)
     &             * nutr_rest_time_inv(k)
        ELSE
           RESTORE = c0
        END IF
        SiO3_RESTORE_HIST(ISTR:IEND,J,K) = RESTORE
        TRACER(istr:iend,K,sio3_ind,curtime) = RESTORE +
     &     sinking_particle_P_SiO2(6,istr:iend,j) +
     &      Qsi * (0.5_8 * graze_diat + 0.95_8 * diat_loss)
     &      - photoSi_diat
     &       + remin_sed_si / Hz(istr:iend,j,1)
        IF (lrest_po4) THEN
            RESTORE = (PO4_CLIM(istr:iend,j,k) - PO4_loc)
     &            * nutr_rest_time_inv(k)
        ELSE
           RESTORE = c0
        END IF
        PO4_RESTORE_HIST(ISTR:IEND,J,K) = RESTORE
        TRACER(istr:iend,K,po4_ind,curtime) = RESTORE + (Qp * (
     &    Sinking_Particle_Poc(6,istr:iend,j) +
     &   zoo_loss_dic + sp_loss_dic + graze_sp_dic + diat_loss_dic +
     &    graze_diat_dic - photoC_sp - photoC_diat))
     &    + DOP_remin + diaz_loss_dip - (photoC_diaz * Qp_diaz)
     &       + Qp*remin_sed_poc / Hz(istr:iend,j,1)
        TRACER(istr:iend,K,spC_ind,curtime) = photoC_sp -
     &       graze_sp - sp_loss - sp_agg
        TRACER(istr:iend,K,spChl_ind,curtime) = photoacc_sp -
     &       thetaC_sp * (graze_sp + sp_loss + sp_agg)
       TRACER(istr:iend,K,spCaCO3_ind,curtime) = CaCO3_prod -
     &    (graze_sp + sp_loss + sp_agg) * QCaCO3
       TRACER(istr:iend,K,diatC_ind,curtime) =
     &   photoC_diat - graze_diat -
     &    diat_loss - diat_agg
       TRACER(istr:iend,K,diatChl_ind,curtime) =
     &    photoacc_diat -
     &    thetaC_diat * (graze_diat + diat_loss + diat_agg)
       TRACER(istr:iend,K,zooC_ind,curtime) = graze_sp_zoo
     &  + graze_diat_zoo
     &  + graze_diaz_zoo - zoo_loss
       TRACER(istr:iend,K,doc_ind,curtime) =
     &                 DOC_prod - DOC_remin
       TRACER(istr:iend,K,don_ind,curtime) =
     &                 DON_prod - DON_remin
       TRACER(istr:iend,K,dop_ind,curtime) =
     &                 DOP_prod - DOP_remin
       TRACER(istr:iend,K,dofe_ind,curtime) =
     &                 DOFe_prod - DOFe_remin
       TRACER(istr:iend,K,spFe_ind,curtime) =  photoFe_sp
     &    - (Qfe_sp * (graze_sp+sp_loss+sp_agg))
       TRACER(istr:iend,K,diatFe_ind,curtime) =  photoFe_diat
     &  - (Qfe_diat * (graze_diat+diat_loss+diat_agg))
       TRACER(istr:iend,K,diatSi_ind,curtime) =  photoSi_diat
     &  - (Qsi * (graze_diat+diat_loss+diat_agg))
       TRACER(istr:iend,K,diazC_ind,curtime) =  photoC_diaz
     &      - graze_diaz - diaz_loss - diaz_agg
       TRACER(istr:iend,K,diazChl_ind,curtime) = photoacc_diaz
     &    - thetaC_diaz * (graze_diaz + diaz_loss + diaz_agg)
        TRACER(istr:iend,K,diazFe_ind,curtime) =  photoFe_diaz
     &          - (Qfe_diaz * (graze_diaz + diaz_loss + diaz_agg))
        TRACER(istr:iend,K,dic_ind,curtime) = DOC_remin +
     &       Sinking_Particle_Poc(6,istr:iend,j) +
     &       sinking_particle_P_CaCO3(6,istr:iend,j) +
     &    0.33_8 * graze_sp * QCaCO3 + zoo_loss_dic + sp_loss_dic +
     &    graze_sp_dic + diat_loss_dic + graze_diat_dic -
     &    photoC_sp - photoC_diat
     &   - CaCO3_prod + graze_diaz_dic + diaz_loss_dic - photoC_diaz
     &       + (remin_sed_poc + remin_sed_caco3) / Hz(istr:iend,j,1)
       TRACER(istr:iend,K,alk_ind,curtime) =
     &     -TRACER(istr:iend,K,no3_ind,curtime) +
     &       TRACER(istr:iend,K,nh4_ind,curtime) +
     &    c2 * (
     &       sinking_particle_P_CaCO3(6,istr:iend,j)+
     &    0.33_8 * graze_sp * QCaCO3 - CaCO3_prod)
     &       + c2 * remin_sed_caco3 / Hz(istr:iend,j,1)
        TRACER(istr:iend,K,o2_ind,curtime) =
     &      (photoC_sp + photoC_diat + photoC_diaz)
     &          / parm_Red_D_C_O2
       WHERE (O2_loc > parm_o2_min)
          TRACER(istr:iend,K,o2_ind,curtime) =
     &      TRACER(istr:iend,K,o2_ind,curtime)
     &     +((
     &         - Sinking_Particle_Poc(6,istr:iend,j)
     &       - DOC_remin - zoo_loss_dic - sp_loss_dic - graze_sp_dic-
     &       diat_loss_dic - graze_diat_dic - graze_diaz_dic
     &       - diaz_loss_dic) / parm_Red_P_C_O2)
     &         - 1.5_8 * ammox - 0.5_8 * nitrif
     &       - remin_sed_poc / (parm_Red_P_C_O2 * Hz(istr:iend,j,1))
       END WHERE
       if (k .eq. 1) then
          tracer_sed(istr:iend,sed_poc_ind,curtime) =
     &         bot_flux_poc
     &         - remin_sed_poc
          tracer_sed(istr:iend,sed_caco3_ind,curtime) =
     &         bot_flux_caco3
     &         - remin_sed_caco3
          tracer_sed(istr:iend,sed_si_ind,curtime) =
     &         bot_flux_si
     &         - remin_sed_si
       end if
        PAR(ISTR:IEND,J,K)                  = PAR_lay
        graze_sp_HIST(ISTR:IEND,J,K)        = graze_sp
        graze_diat_HIST(ISTR:IEND,J,K)      = graze_diat
        graze_diaz_HIST(ISTR:IEND,J,K)   = graze_diaz
        graze_tot_HIST(ISTR:IEND,J,K)       = graze_sp
     &                        + graze_diat + graze_diaz
        sp_loss_HIST(ISTR:IEND,J,K)         = sp_loss
        diat_loss_HIST(ISTR:IEND,J,K)       = diat_loss
        diaz_loss_HIST(ISTR:IEND,J,K)       = diaz_loss+diaz_agg
        zoo_loss_HIST(ISTR:IEND,J,K)        = zoo_loss
        sp_agg_HIST(ISTR:IEND,J,K)          = sp_agg
        diat_agg_HIST(ISTR:IEND,J,K)        = diat_agg
        photoC_sp_HIST(ISTR:IEND,J,K)       = photoC_sp
        f_ratio_sp_hist(istr:iend,j,k)      = f_ratio_sp
        photoC_diat_HIST(ISTR:IEND,J,K)     = photoC_diat
        f_ratio_diat_hist(istr:iend,j,k)    = f_ratio_diat
        photoC_diaz_HIST(ISTR:IEND,J,K)     = photoC_diaz
        tot_prod_HIST(ISTR:IEND,J,K)        = photoC_sp +
     &                    photoC_diat + photoC_diaz
        no3_v_sp_hist(ISTR:IEND,J,K)        = no3_v_sp
        nh4_v_sp_hist(ISTR:IEND,J,K)        = nh4_v_sp
        no3_v_diat_hist(ISTR:IEND,J,K)      = no3_v_diat
        nh4_v_diat_hist(ISTR:IEND,J,K)      = nh4_v_diat
        DOC_prod_HIST(ISTR:IEND,J,K)        = DOC_prod
        DOC_remin_HIST(ISTR:IEND,J,K)       = DOC_remin
        DON_prod_HIST(ISTR:IEND,J,K)        = DON_prod
        DON_remin_HIST(ISTR:IEND,J,K)       = DON_remin
        DOP_prod_HIST(ISTR:IEND,J,K)        = DOP_prod
        DOP_remin_HIST(ISTR:IEND,J,K)       = DOP_remin
        DOFe_prod_HIST(ISTR:IEND,J,K)       = DOFe_prod
        DOFe_remin_HIST(ISTR:IEND,J,K)      = DOFe_remin
        Fe_scavenge_HIST(ISTR:IEND,J,K)     = Fe_scavenge
        Fe_scavenge_rate_HIST(ISTR:IEND,J,K) = Fe_scavenge_rate
        ammox_HIST(ISTR:IEND,J,K)          = ammox
        denitr_hist(istr:iend,j,k)         = denitr
        n2o_prod_hist(istr:iend,j,k)       = n2o_prod
        n2_prod_hist(istr:iend,j,k)        = n2_prod
        if (k .eq. 1)
     &       denitr_sed_hist(istr:iend,j)  = denitr_sed
        if (k .eq. 1) then
           bot_flux_poc_hist(istr:iend,j)     = bot_flux_poc
           bot_flux_caco3_hist(istr:iend,j)   = bot_flux_caco3
           bot_flux_si_hist(istr:iend,j)      = bot_flux_si
        end if
        TRACER(istr:iend,k,:,curtime)=
     &      TRACER(istr:iend,k,:,curtime)*dt
     &      +TRACER(istr:iend,k,:,ctime)
        do m = 1, ntrc_bio
           do i = istr,iend
              if (tracer(i,k,m,ctime) .ge. 0.0_8 .and.
     &             tracer(i,k,m,curtime) .lt. -0.001_8) then
                 print *,'i,j,k,m,before/after:',i,j,k,m,
     &                tracer(i,k,m,ctime),
     &                tracer(i,k,m,curtime)
              end if
           end do
        end do
        tracer(istr:iend,k,:,ctime)=
     &       tracer(istr:iend,k,:,curtime)
        if (k .eq. 1) then
           tracer_sed(istr:iend,:,ctime)=
     &          tracer_sed(istr:iend,:,curtime)*dt
     &          +tracer_sed(istr:iend,:,ctime)
        end if
         end subroutine ecosys_set_interior
         SUBROUTINE init_particulate_terms(QA_dust_def,istr,
     &      iend,j,net_dust_in)
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
        real(kind=8) c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1._8, c0=0.0_8,c2=2._8,
     &  c1000=1000._8,p5=0.5_8,
     &  spd = 86400.0_8,  dps = c1 / spd ,
     &  t0_kelvin= 273.16_8)
        real(kind=8)  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C
        parameter ( parm_Red_D_C_P  = 117.0_8,
     &  parm_Red_P_C_P  = 117.0_8,
     &  parm_Red_D_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_P_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_D_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_P_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_Fe_C   = 3.0D-6)
        real(kind=8), parameter :: Q         = 0.137_8
        real(kind=8), parameter :: Qinv      = 1.0_8/0.137_8
        real(kind=8), parameter :: Qp        = 0.00855_8
        real(kind=8), parameter :: Qp_diat   = 0.1_8 * Q
        real(kind=8), parameter :: Qdenit    = 0.9811320755_8
        real(kind=8)  parm_Fe_bioavail,   parm_prod_dissolve,
     &  parm_o2_min,     parm_Rain_CaCO3, parm_Rain_SiO2,
     &  parm_kappa_nitrif,  parm_nitrif_par_lim,  parm_POC_flux_ref,
     &  parm_rest_prod_tau,  parm_rest_prod_z_c,  parm_z_umax_0,
     &  parm_diat_umax_0,     parm_z_mort_0,     parm_z_mort2_0,
     &  parm_sd_remin_0,     parm_sp_kNO3,       parm_diat_kNO3,
     &  parm_sp_kNH4,        parm_diat_kNH4,     parm_sp_kFe,
     &  parm_diat_kFe,       parm_diat_kSiO3,    parm_sp_kPO4,
     &  parm_diat_kPO4,      parm_z_grz,        parm_alphaChl,
     &  parm_labile_ratio,   parm_alphaDiaz,    parm_diaz_umax_0,
     &  gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
       common/eco_para/parm_Fe_bioavail, parm_prod_dissolve,
     &   parm_o2_min, parm_Rain_CaCO3, parm_Rain_SiO2,
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref,
     &   parm_rest_prod_tau, parm_rest_prod_z_c, parm_z_umax_0,
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0,
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe,
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4,
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl,
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
        real(kind=8) parm_kappa_nitrif_no2
        real(kind=8) parm_n2o_prod
        real(kind=8) parm_n2_prod
        common /eco_param_ncycle_anoxic/ parm_kappa_nitrif_no2,
     &       parm_n2o_prod,parm_n2_prod
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
         integer(kind=4) istr,iend,j
            REAL(kind=8), DIMENSION(istr:iend), INTENT(OUT) ::
     &          QA_dust_def
             REAL(kind=8), DIMENSION(istr:iend) ::
     &         net_dust_in
             diss(1)      = 130.0_8
             gamma(1)     = 0.4_8
             mass(1)      = 12.01_8
             rhoo(1)      = c1
             diss(2)  = 600.0_8
             gamma(2) = 0.55_8
             mass(2)  = 100.09_8
             rhoo(2)   = 0.07_8 * mass(2) / mass(1)
           diss(3)   = 210.0_8
           gamma(3)  = 0.37_8
           mass(3)   = 60.08_8
            rhoo(3)    = 0.035_8 * mass(3) / mass(1)
           diss(4)     = 600.0_8
           gamma(4)     = 0.97_8
          mass(4)      = 1.0D6
          rhoo(4)      = 0.07_8 * mass(4)  / mass(1)
            diss(5)   = 600.0_8
            gamma(5)  = c0
            mass(5)   = c0
            rhoo(5)   = c0
         sinking_particle_P_CaCO3(4,istr:iend,j) = c0
         sinking_particle_P_CaCO3(5,istr:iend,j) = c0
        sinking_particle_P_SiO2(4,istr:iend,j) = c0
        sinking_particle_P_SiO2(5,istr:iend,j) = c0
      if (ldust_flux) then
       net_dust_in = net_dust_in * (c1 - parm_fe_bioavail)
       sinking_particle_dust(4,istr:iend,j) =
     &     (c1 - gamma(4)) * net_dust_in
       sinking_particle_dust(5,istr:iend,j) =
     &     gamma(4) * net_dust_in
      ELSE
            sinking_particle_dust(4,istr:iend,j) = c0
            sinking_particle_dust(5,istr:iend,j) = c0
      END IF
            sinking_particle_P_iron(4,istr:iend,j) = c0
            sinking_particle_P_iron(5,istr:iend,j) = c0
           sinking_particle_POC(4,istr:iend,j) = c0
           sinking_particle_POC(5,istr:iend,j) = c0
           QA_dust_def = rhoo(4) *
     &    (sinking_particle_dust(4,istr:iend,j) +
     &           sinking_particle_dust(5,istr:iend,j))
           END SUBROUTINE init_particulate_terms
             SUBROUTINE compute_particulate_terms(k, QA_dust_def,
     &                 temp,istr,iend,j, bot_flux_poc
     &        , bot_flux_caco3, bot_flux_si
     &        )
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
        real(kind=8) c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1._8, c0=0.0_8,c2=2._8,
     &  c1000=1000._8,p5=0.5_8,
     &  spd = 86400.0_8,  dps = c1 / spd ,
     &  t0_kelvin= 273.16_8)
        real(kind=8)  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C
        parameter ( parm_Red_D_C_P  = 117.0_8,
     &  parm_Red_P_C_P  = 117.0_8,
     &  parm_Red_D_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_P_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_D_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_P_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_Fe_C   = 3.0D-6)
        real(kind=8), parameter :: Q         = 0.137_8
        real(kind=8), parameter :: Qinv      = 1.0_8/0.137_8
        real(kind=8), parameter :: Qp        = 0.00855_8
        real(kind=8), parameter :: Qp_diat   = 0.1_8 * Q
        real(kind=8), parameter :: Qdenit    = 0.9811320755_8
        real(kind=8)  parm_Fe_bioavail,   parm_prod_dissolve,
     &  parm_o2_min,     parm_Rain_CaCO3, parm_Rain_SiO2,
     &  parm_kappa_nitrif,  parm_nitrif_par_lim,  parm_POC_flux_ref,
     &  parm_rest_prod_tau,  parm_rest_prod_z_c,  parm_z_umax_0,
     &  parm_diat_umax_0,     parm_z_mort_0,     parm_z_mort2_0,
     &  parm_sd_remin_0,     parm_sp_kNO3,       parm_diat_kNO3,
     &  parm_sp_kNH4,        parm_diat_kNH4,     parm_sp_kFe,
     &  parm_diat_kFe,       parm_diat_kSiO3,    parm_sp_kPO4,
     &  parm_diat_kPO4,      parm_z_grz,        parm_alphaChl,
     &  parm_labile_ratio,   parm_alphaDiaz,    parm_diaz_umax_0,
     &  gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
       common/eco_para/parm_Fe_bioavail, parm_prod_dissolve,
     &   parm_o2_min, parm_Rain_CaCO3, parm_Rain_SiO2,
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref,
     &   parm_rest_prod_tau, parm_rest_prod_z_c, parm_z_umax_0,
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0,
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe,
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4,
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl,
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
        real(kind=8) parm_kappa_nitrif_no2
        real(kind=8) parm_n2o_prod
        real(kind=8) parm_n2_prod
        common /eco_param_ncycle_anoxic/ parm_kappa_nitrif_no2,
     &       parm_n2o_prod,parm_n2_prod
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
           INTEGER(kind=4), INTENT(IN) :: k,istr,iend,j
          REAL(kind=8), DIMENSION(istr:iend), INTENT(INOUT) ::
     &       QA_dust_def
          REAL(kind=8), DIMENSION(istr:iend), INTENT(IN) ::
     &        temp
          REAL::
     &    decay_CaCO3,
     &    decay_dust,
     &    decay_POC_E,
     &    decay_SiO2,
     &    decay_Hard,
     &    POC_prod_avail,
     &    new_QA_dust_def
       REAL(kind=8), DIMENSION(istr:iend) :: bot_flux_poc
       REAL(kind=8), DIMENSION(istr:iend) :: bot_flux_caco3,bot_flux_si
       REAL(kind=8) :: Tfunc_POC, Tfunc_SiO2
       INTEGER:: i
       bot_flux_poc = c0
       bot_flux_caco3 = c0
       bot_flux_si = c0
         Sinking_Particle_P_Caco3(1,istr:iend,j) =
     &         Sinking_Particle_P_Caco3(4,istr:iend,j)
         Sinking_Particle_P_Caco3(2,istr:iend,j) =
     &         Sinking_Particle_P_Caco3(5,istr:iend,j)
          Sinking_Particle_P_Sio2(1,istr:iend,j) =
     &           Sinking_Particle_P_Sio2(4,istr:iend,j)
          Sinking_Particle_P_Sio2(2,istr:iend,j) =
     &           Sinking_Particle_P_Sio2(5,istr:iend,j)
          sinking_particle_dust(1,istr:iend,j) =
     &             sinking_particle_dust(4,istr:iend,j)
           sinking_particle_dust(2,istr:iend,j) =
     &              sinking_particle_dust(5,istr:iend,j)
        SINKING_PARTICLE_POC(1,istr:iend,j) =
     &        SINKING_PARTICLE_POC(4,istr:iend,j)
        SINKING_PARTICLE_POC(2,istr:iend,j) =
     &          SINKING_PARTICLE_POC(5,istr:iend,j)
          sinking_particle_p_iron(1,istr:iend,j) =
     &              sinking_particle_p_iron(4,istr:iend,j)
           sinking_particle_p_iron(2,istr:iend,j) =
     &               sinking_particle_p_iron(5,istr:iend,j)
         DO i = istr,iend
          IF (LANDMASK(i,j) ) THEN
             Tfunc_POC = 1.12_8 ** (0.1_8 * temp(i) - 3.0_8)
             Tfunc_SiO2 = 4.0_8 ** (0.1_8 * temp(i) - 3.0_8)
             decay_POC_E = EXP(-Hz(i,j,k) / diss(1) * Tfunc_POC)
             decay_SiO2  = EXP(-Hz(i,j,k) / diss(3) * Tfunc_SiO2)
             decay_CaCO3 = EXP(-Hz(i,j,k) / Diss(2))
             decay_dust  = EXP(-Hz(i,j,k) / diss(4))
             decay_Hard  = EXP(-Hz(i,j,k) * 2.5D-5)
              Sinking_Particle_P_Caco3(4,i,j) =
     &              Sinking_Particle_P_Caco3(1,i,j) * decay_CaCO3 +
     &              Sinking_Particle_P_Caco3(3,i,j) *
     &             ((c1 - Gamma(2)) * (c1 - decay_CaCO3) * Diss(2))
              Sinking_Particle_P_Caco3(5,i,j) =
     &               Sinking_Particle_P_Caco3(2,i,j) * decay_Hard +
     &               Sinking_Particle_P_Caco3(3,i,j) *
     &               (Gamma(2) * Hz(i,j,k))
               Sinking_Particle_P_Sio2(4,i,j) =
     &               Sinking_Particle_P_Sio2(1,i,j) * decay_SiO2 +
     &               Sinking_Particle_P_Sio2(3,i,j) *
     &               ((c1 - Gamma(3)) * (c1 - decay_SiO2)
     &               * (diss(3) / Tfunc_SiO2))
                Sinking_Particle_P_Sio2(5,i,j) =
     &               Sinking_Particle_P_Sio2(2,i,j) * decay_Hard +
     &               Sinking_Particle_P_Sio2(3,i,j) *
     &              (Gamma(3) * Hz(i,j,k))
                sinking_particle_dust(4,i,j) =
     &                  sinking_particle_dust(1,i,j) * decay_dust
                 sinking_particle_dust(5,i,j) =
     &                  sinking_particle_dust(2,i,j) * decay_Hard
               POC_prod_avail = Sinking_Particle_POC(3,i,j) -
     &             rhoo(2) * Sinking_Particle_P_Caco3(3,i,j) -
     &             rhoo(3) * Sinking_Particle_P_Sio2(3,i,j)
              IF (POC_prod_avail < c0) THEN
               print *,"subroutine compute_particulate_terms:mass ",
     &           " ratio of ballast production exceeds POC production"
               print *, 'POC_prod_avail: ', POC_prod_avail
             END IF
             IF (QA_dust_def(i) > 0) THEN
                new_QA_dust_def = QA_dust_def(i) *
     &               (sinking_particle_dust(4,i,j) +
     &                 sinking_particle_dust(5,i,j)) /
     &                (sinking_particle_dust(1,i,j) +
     &                 sinking_particle_dust(2,i,j))
             ELSE
                new_QA_dust_def = c0
             END IF
             IF (new_QA_dust_def > c0) THEN
                new_QA_dust_def = new_QA_dust_def -
     &                    POC_prod_avail * Hz(i,j,k)
                IF (new_QA_dust_def < c0) THEN
                   POC_prod_avail = -new_QA_dust_def / Hz(i,j,k)
                   new_QA_dust_def = c0
                ELSE
                   POC_prod_avail = c0
                END IF
             END IF
               QA_dust_def(i) = new_QA_dust_def
             IF (SINKING_PARTICLE_POC(2,i,j) == c0 .AND.
     &                  SINKING_PARTICLE_POC(3,i,j) == c0) THEN
                SINKING_PARTICLE_POC(5,i,j) = c0
             ELSE
                SINKING_PARTICLE_POC(5,i,j) = rhoo(2) *
     &                (Sinking_Particle_P_Caco3(4,i,j) +
     &                Sinking_Particle_P_Caco3(5,i,j)) + rhoo(3) *
     &                (Sinking_Particle_P_Sio2(4,i,j) +
     &                Sinking_Particle_P_Sio2(5,i,j)) +
     &                rhoo(4) * (sinking_particle_dust(4,i,j) +
     &                sinking_particle_dust(5,i,j)) -new_QA_dust_def
                 SINKING_PARTICLE_POC(5,i,j) =
     &                   MAX(SINKING_PARTICLE_POC(5,i,j), c0)
             END IF
             SINKING_PARTICLE_POC(4,i,j) = SINKING_PARTICLE_POC(1,i,j)
     &                * decay_POC_E + POC_prod_avail
     &               *((c1 - decay_POC_E) * (diss(1) / Tfunc_POC))
             Sinking_Particle_P_Caco3(6,i,j) =
     &              Sinking_Particle_P_Caco3(3,i,j) +
     &             ((Sinking_Particle_P_Caco3(1,i,j) -
     &             Sinking_Particle_P_Caco3(4,i,j)) +
     &             (Sinking_Particle_P_Caco3(2,i,j) -
     &             Sinking_Particle_P_Caco3(5,i,j))) / Hz(i,j,k)
             Sinking_Particle_P_Sio2(6,i,j) =
     &              Sinking_Particle_P_Sio2(3,i,j) +
     &             ((Sinking_Particle_P_Sio2(1,i,j) -
     &             Sinking_Particle_P_Sio2(4,i,j)) +
     &             (Sinking_Particle_P_Sio2(2,i,j) -
     &             Sinking_Particle_P_Sio2(5,i,j))) / Hz(i,j,k)
             SINKING_PARTICLE_POC(6,i,j) =
     &              SINKING_PARTICLE_POC(3,i,j) +
     &             ((SINKING_PARTICLE_POC(1,i,j) -
     &             SINKING_PARTICLE_POC(4,i,j)) +
     &             (SINKING_PARTICLE_POC(2,i,j) -
     &             SINKING_PARTICLE_POC(5,i,j))) / Hz(i,j,k)
             sinking_particle_dust(6,i,j) =
     &            ((sinking_particle_dust(1,i,j) -
     &            sinking_particle_dust(4,i,j)) +
     &             (sinking_particle_dust(2,i,j) -
     &             sinking_particle_dust(5,i,j))) / Hz(i,j,k)
             IF (SINKING_PARTICLE_POC(1,i,j) +
     &                 SINKING_PARTICLE_POC(2,i,j) == c0) THEN
                sinking_particle_p_iron(6,i,j) =
     &               (SINKING_PARTICLE_POC(6,i,j) * parm_Red_Fe_C)
             ELSE
                sinking_particle_p_iron(6,i,j) =
     &               (SINKING_PARTICLE_POC(6,i,j) *
     &                (sinking_particle_p_iron(1,i,j) +
     &                sinking_particle_p_iron(2,i,j)) /
     &                (SINKING_PARTICLE_POC(1,i,j) +
     &                 SINKING_PARTICLE_POC(2,i,j)))
             END IF
             sinking_particle_p_iron(4,i,j) =
     &            sinking_particle_p_iron(1,i,j) + Hz(i,j,k) *
     &            ((c1 - gamma(5)) * sinking_particle_p_iron(3,i,j)
     &           - sinking_particle_p_iron(6,i,j))
             IF (sinking_particle_p_iron(4,i,j) < c0) THEN
                sinking_particle_p_iron(4,i,j) = c0
                sinking_particle_p_iron(6,i,j) =
     &              sinking_particle_p_iron(1,i,j) / Hz(i,j,k) +
     &               (c1 - gamma(5)) * sinking_particle_p_iron(3,i,j)
             END IF
             sinking_particle_p_iron(6,i,j) =
     &           sinking_particle_p_iron(6,i,j) +
     &           sinking_particle_dust(6,i,j) * 626.712_8
              sinking_particle_p_iron(5,i,j) =
     &               sinking_particle_p_iron(2,i,j)
          ELSE
             Sinking_Particle_P_Caco3(4,i,j) = c0
             Sinking_Particle_P_Caco3(5,i,j) = c0
             Sinking_Particle_P_Caco3(6,i,j) = c0
             Sinking_Particle_P_Sio2(4,i,j) = c0
             Sinking_Particle_P_Sio2(5,i,j) = c0
             Sinking_Particle_P_Sio2(6,i,j) = c0
             sinking_particle_dust(4,i,j) = c0
             sinking_particle_dust(5,i,j) = c0
             sinking_particle_dust(6,i,j) = c0
             SINKING_PARTICLE_POC(4,i,j) = c0
             SINKING_PARTICLE_POC(5,i,j) = c0
             SINKING_PARTICLE_POC(6,i,j) = c0
             sinking_particle_p_iron(4,i,j) = c0
             sinking_particle_p_iron(5,i,j) = c0
             sinking_particle_p_iron(6,i,j) = c0
          END IF
          IF (LANDMASK(i,j) .AND. k == 1) THEN
             bot_flux_caco3(i) = SINKING_PARTICLE_P_CaCO3(4,i,j) +
     &             SINKING_PARTICLE_P_CaCO3(5,i,j)
             bot_flux_si(i) = SINKING_PARTICLE_P_SiO2(4,i,j) +
     &             SINKING_PARTICLE_P_SiO2(5,i,j)
             sinking_particle_dust(6,i,j) =
     &             sinking_particle_dust(6,i,j) +
     &             (sinking_particle_dust(4,i,j) +
     &            sinking_particle_dust(5,i,j)) / Hz(i,j,k)
             sinking_particle_dust(4,i,j) = c0
             sinking_particle_dust(5,i,j) = c0
             bot_flux_poc(i) = SINKING_PARTICLE_POC(4,i,j) +
     &             SINKING_PARTICLE_POC(5,i,j)
             sinking_particle_p_iron(6,i,j) =
     &          sinking_particle_p_iron(6,i,j) +
     &             (sinking_particle_p_iron(4,i,j) +
     &         sinking_particle_p_iron(5,i,j)) / Hz(i,j,k)
             sinking_particle_p_iron(4,i,j) = c0
             sinking_particle_p_iron(5,i,j) = c0
           END IF
         END DO
         POC_FLUX_IN_hist(istr:iend,j,k) =
     &               SINKING_PARTICLE_POC(1,istr:iend,j) +
     &                SINKING_PARTICLE_POC(2,istr:iend,j)
         POC_PROD_hist(istr:iend,j,k)  =
     &                SINKING_PARTICLE_POC(3,istr:iend,j)
         POC_REMIN_hist(istr:iend,j,k)       =
     &                 SINKING_PARTICLE_POC(6,istr:iend,j)
         CaCO3_FLUX_IN_HIST(istr:iend,j,k)  =
     &           Sinking_Particle_P_Caco3(1,istr:iend,j) +
     &           Sinking_Particle_P_Caco3(2,istr:iend,j)
         CaCO3_PROD_hist(istr:iend,j,k)      =
     &           Sinking_Particle_P_Caco3(3,istr:iend,j)
         CaCO3_REMIN_hist(istr:iend,j,k)     =
     &            Sinking_Particle_P_Caco3(6,istr:iend,j)
         SiO2_FLUX_IN_hist(istr:iend,j,k)    =
     &            Sinking_Particle_P_Sio2(1,istr:iend,j) +
     &            Sinking_Particle_P_Sio2(2,istr:iend,j)
         SiO2_PROD_HIST(istr:iend,j,k)      =
     &            Sinking_Particle_P_SiO2(3,istr:iend,j)
         SiO2_REMIN_HIST(istr:iend,j,k)     =
     &             Sinking_Particle_P_SiO2(6,istr:iend,j)
         dust_FLUX_IN_hist(istr:iend,j,k)    =
     &             sinking_particle_dust(1,istr:iend,j) +
     &              sinking_particle_dust(2,istr:iend,j)
         dust_REMIN_hist(istr:iend,j,k)      =
     &               Sinking_Particle_dust(6,istr:iend,j)
         P_iron_FLUX_IN_hist(istr:iend,j,k)  =
     &             sinking_particle_p_iron(1,istr:iend,j) +
     &             sinking_particle_p_iron(2,istr:iend,j)
         P_iron_PROD_HIST(istr:iend,j,k)    =
     &              sinking_particle_p_iron(3,istr:iend,j)
         P_iron_REMIN_hist(istr:iend,j,k)    =
     &             sinking_particle_p_iron(6,istr:iend,j)
          END SUBROUTINE compute_particulate_terms
       subroutine WS(SMFTX, SMFTY,landmask,work,istr,iend)
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
       integer::istr,iend
        REAL(kind=8), DIMENSION(istr:iend) :: WORK
        real(kind=8),parameter::
     &   rho_air   = 1.2_8
        REAL(kind=8), DIMENSION(istr:iend), INTENT(IN) ::
     &    SMFTX,
     &    SMFTY
       REAL(kind=8), PARAMETER ::
     &    coef_1  = 0.0027_8,
     &    coef_2  = 0.000142_8,
     &    coef_3  = 0.0000764_8,
     &    c_d     = 1.7D-3
        REAL(kind=8), DIMENSION(istr:iend) ::
     &    ustar_squared
        LOGICAL,DIMENSION(istr:iend), INTENT(IN) :: landmask
       WHERE (landmask)
         ustar_squared = SQRT(SMFTX**2 + SMFTY**2) * rho0 /
     &         rho_air
         WORK = SQRT(ustar_squared / c_d)
          WORK = WORK - (WORK*(coef_1 + WORK*(coef_2 + WORK*coef_3))
     &         - ustar_squared) /
     &         (coef_1 + WORK*(2*coef_2 + WORK*3*coef_3))
          WORK = WORK - (WORK*(coef_1 + WORK*(coef_2 + WORK*coef_3))
     &          - ustar_squared) /
     &         (coef_1 + WORK*(2*coef_2 + WORK*3*coef_3))
          WORK = WORK - (WORK*(coef_1 + WORK*(coef_2 + WORK*coef_3))
     &           - ustar_squared) /
     &       (coef_1 + WORK*(2*coef_2 + WORK*3*coef_3))
        ELSEWHERE
           WORK = 0.0_8
        END WHERE
        RETURN
        END subroutine WS
         subroutine CSCHMIDT_O2(SSTT,landmask,SCHMIDT_O2,istr,
     &          iend)
        implicit none
        real(kind=8) c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1._8, c0=0.0_8,c2=2._8,
     &  c1000=1000._8,p5=0.5_8,
     &  spd = 86400.0_8,  dps = c1 / spd ,
     &  t0_kelvin= 273.16_8)
        real(kind=8)  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C
        parameter ( parm_Red_D_C_P  = 117.0_8,
     &  parm_Red_P_C_P  = 117.0_8,
     &  parm_Red_D_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_P_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_D_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_P_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_Fe_C   = 3.0D-6)
        real(kind=8), parameter :: Q         = 0.137_8
        real(kind=8), parameter :: Qinv      = 1.0_8/0.137_8
        real(kind=8), parameter :: Qp        = 0.00855_8
        real(kind=8), parameter :: Qp_diat   = 0.1_8 * Q
        real(kind=8), parameter :: Qdenit    = 0.9811320755_8
        real(kind=8)  parm_Fe_bioavail,   parm_prod_dissolve,
     &  parm_o2_min,     parm_Rain_CaCO3, parm_Rain_SiO2,
     &  parm_kappa_nitrif,  parm_nitrif_par_lim,  parm_POC_flux_ref,
     &  parm_rest_prod_tau,  parm_rest_prod_z_c,  parm_z_umax_0,
     &  parm_diat_umax_0,     parm_z_mort_0,     parm_z_mort2_0,
     &  parm_sd_remin_0,     parm_sp_kNO3,       parm_diat_kNO3,
     &  parm_sp_kNH4,        parm_diat_kNH4,     parm_sp_kFe,
     &  parm_diat_kFe,       parm_diat_kSiO3,    parm_sp_kPO4,
     &  parm_diat_kPO4,      parm_z_grz,        parm_alphaChl,
     &  parm_labile_ratio,   parm_alphaDiaz,    parm_diaz_umax_0,
     &  gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
       common/eco_para/parm_Fe_bioavail, parm_prod_dissolve,
     &   parm_o2_min, parm_Rain_CaCO3, parm_Rain_SiO2,
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref,
     &   parm_rest_prod_tau, parm_rest_prod_z_c, parm_z_umax_0,
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0,
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe,
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4,
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl,
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
        real(kind=8) parm_kappa_nitrif_no2
        real(kind=8) parm_n2o_prod
        real(kind=8) parm_n2_prod
        common /eco_param_ncycle_anoxic/ parm_kappa_nitrif_no2,
     &       parm_n2o_prod,parm_n2_prod
       integer::istr,iend
         LOGICAL,DIMENSION(istr:iend) :: landmask
          REAL(kind=8), DIMENSION(istr:iend), INTENT(IN) :: SSTT
          REAL(kind=8), DIMENSION(istr:iend) :: SCHMIDT_O2
          REAL(kind=8), PARAMETER ::
     &      a = 1638.0_8,
     &      b = 81.83_8,
     &      c = 1.483_8,
     &      d = 0.008004_8
         WHERE (LANDMASK)
           SCHMIDT_O2 = a + SSTT* (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
            SCHMIDT_O2 = c0
          END WHERE
          return
          END subroutine CSCHMIDT_O2
          subroutine O2SATU(SSTT, SSSS,landmask, O2SAT,istr,iend)
        implicit none
        real(kind=8) c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1._8, c0=0.0_8,c2=2._8,
     &  c1000=1000._8,p5=0.5_8,
     &  spd = 86400.0_8,  dps = c1 / spd ,
     &  t0_kelvin= 273.16_8)
        real(kind=8)  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C
        parameter ( parm_Red_D_C_P  = 117.0_8,
     &  parm_Red_P_C_P  = 117.0_8,
     &  parm_Red_D_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_P_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_D_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_P_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_Fe_C   = 3.0D-6)
        real(kind=8), parameter :: Q         = 0.137_8
        real(kind=8), parameter :: Qinv      = 1.0_8/0.137_8
        real(kind=8), parameter :: Qp        = 0.00855_8
        real(kind=8), parameter :: Qp_diat   = 0.1_8 * Q
        real(kind=8), parameter :: Qdenit    = 0.9811320755_8
        real(kind=8)  parm_Fe_bioavail,   parm_prod_dissolve,
     &  parm_o2_min,     parm_Rain_CaCO3, parm_Rain_SiO2,
     &  parm_kappa_nitrif,  parm_nitrif_par_lim,  parm_POC_flux_ref,
     &  parm_rest_prod_tau,  parm_rest_prod_z_c,  parm_z_umax_0,
     &  parm_diat_umax_0,     parm_z_mort_0,     parm_z_mort2_0,
     &  parm_sd_remin_0,     parm_sp_kNO3,       parm_diat_kNO3,
     &  parm_sp_kNH4,        parm_diat_kNH4,     parm_sp_kFe,
     &  parm_diat_kFe,       parm_diat_kSiO3,    parm_sp_kPO4,
     &  parm_diat_kPO4,      parm_z_grz,        parm_alphaChl,
     &  parm_labile_ratio,   parm_alphaDiaz,    parm_diaz_umax_0,
     &  gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
       common/eco_para/parm_Fe_bioavail, parm_prod_dissolve,
     &   parm_o2_min, parm_Rain_CaCO3, parm_Rain_SiO2,
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref,
     &   parm_rest_prod_tau, parm_rest_prod_z_c, parm_z_umax_0,
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0,
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe,
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4,
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl,
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
        real(kind=8) parm_kappa_nitrif_no2
        real(kind=8) parm_n2o_prod
        real(kind=8) parm_n2_prod
        common /eco_param_ncycle_anoxic/ parm_kappa_nitrif_no2,
     &       parm_n2o_prod,parm_n2_prod
       integer::istr,iend
         LOGICAL,DIMENSION(istr:iend) :: landmask
       REAL(kind=8), DIMENSION(istr:iend), INTENT(IN) ::
     &    SSTT,
     &    SSSS
        REAL(kind=8), DIMENSION(istr:iend) :: O2SAT
        REAL(kind=8), DIMENSION(istr:iend) :: TS
        REAL(kind=8), PARAMETER ::
     &    a_0 = 2.00907_8,
     &    a_1 = 3.22014_8,
     &    a_2 = 4.05010_8,
     &    a_3 = 4.94457_8,
     &    a_4 = -2.56847D-1,
     &    a_5 = 3.88767_8,
     &    b_0 = -6.24523D-3,
     &    b_1 = -7.37614D-3,
     &    b_2 = -1.03410D-2,
     &    b_3 = -8.17083D-3,
     &    c_0 = -4.88682D-7
        WHERE (LANDMASK)
           TS = LOG( ((T0_Kelvin+25.0_8) - SSTT) / (T0_Kelvin + SSTT) )
         O2SAT = EXP(a_0+TS*(a_1+TS*(a_2+TS*(a_3+TS*(a_4+TS*a_5)))) +
     &         SSSS*( (b_0+TS*(b_1+TS*(b_2+TS*b_3))) + SSSS*c_0 ))
        ELSEWHERE
          O2SAT = c0
       END WHERE
        O2SAT = O2SAT * 44.6596_8
        return
        END subroutine O2SATU
          subroutine CSCHMIDT_CO2(SSTT,landmask,SCHMIDT_CO2,istr,iend)
        implicit none
        real(kind=8) c1, c0, c2,c1000,p5,spd,dps,t0_kelvin
         parameter ( c1=1._8, c0=0.0_8,c2=2._8,
     &  c1000=1000._8,p5=0.5_8,
     &  spd = 86400.0_8,  dps = c1 / spd ,
     &  t0_kelvin= 273.16_8)
        real(kind=8)  parm_Red_D_C_P,parm_Red_P_C_P,parm_Red_D_C_N,
     &  parm_Red_P_C_N,parm_Red_D_C_O2,parm_Red_P_C_O2,parm_Red_Fe_C
        parameter ( parm_Red_D_C_P  = 117.0_8,
     &  parm_Red_P_C_P  = 117.0_8,
     &  parm_Red_D_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_P_C_N  = 117.0_8 / 16.0_8,
     &  parm_Red_D_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_P_C_O2 = 117.0_8 / 150.0_8,
     &  parm_Red_Fe_C   = 3.0D-6)
        real(kind=8), parameter :: Q         = 0.137_8
        real(kind=8), parameter :: Qinv      = 1.0_8/0.137_8
        real(kind=8), parameter :: Qp        = 0.00855_8
        real(kind=8), parameter :: Qp_diat   = 0.1_8 * Q
        real(kind=8), parameter :: Qdenit    = 0.9811320755_8
        real(kind=8)  parm_Fe_bioavail,   parm_prod_dissolve,
     &  parm_o2_min,     parm_Rain_CaCO3, parm_Rain_SiO2,
     &  parm_kappa_nitrif,  parm_nitrif_par_lim,  parm_POC_flux_ref,
     &  parm_rest_prod_tau,  parm_rest_prod_z_c,  parm_z_umax_0,
     &  parm_diat_umax_0,     parm_z_mort_0,     parm_z_mort2_0,
     &  parm_sd_remin_0,     parm_sp_kNO3,       parm_diat_kNO3,
     &  parm_sp_kNH4,        parm_diat_kNH4,     parm_sp_kFe,
     &  parm_diat_kFe,       parm_diat_kSiO3,    parm_sp_kPO4,
     &  parm_diat_kPO4,      parm_z_grz,        parm_alphaChl,
     &  parm_labile_ratio,   parm_alphaDiaz,    parm_diaz_umax_0,
     &  gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
       common/eco_para/parm_Fe_bioavail, parm_prod_dissolve,
     &   parm_o2_min, parm_Rain_CaCO3, parm_Rain_SiO2,
     &   parm_kappa_nitrif, parm_nitrif_par_lim, parm_POC_flux_ref,
     &   parm_rest_prod_tau, parm_rest_prod_z_c, parm_z_umax_0,
     &   parm_diat_umax_0, parm_z_mort_0, parm_z_mort2_0,
     &   parm_sd_remin_0, parm_sp_kNO3, parm_diat_kNO3,
     &   parm_sp_kNH4, parm_diat_kNH4, parm_sp_kFe,
     &   parm_diat_kFe, parm_diat_kSiO3, parm_sp_kPO4,
     &   parm_diat_kPO4, parm_z_grz, parm_alphaChl,
     &   parm_labile_ratio, parm_alphaDiaz, parm_diaz_umax_0,
     &   gQsi_0, gQsi_coef, gQsi_max
     &  , t_remin_sed_poc, t_remin_sed_caco3, t_remin_sed_si
        real(kind=8) parm_kappa_nitrif_no2
        real(kind=8) parm_n2o_prod
        real(kind=8) parm_n2_prod
        common /eco_param_ncycle_anoxic/ parm_kappa_nitrif_no2,
     &       parm_n2o_prod,parm_n2_prod
       integer::istr,iend,jstr,jend
         LOGICAL,DIMENSION(istr:iend) :: landmask
          REAL(kind=8), DIMENSION(istr:iend), INTENT(IN) :: SSTT
          REAL(kind=8), DIMENSION(istr:iend) :: SCHMIDT_CO2
         REAL(kind=8), PARAMETER ::
     &    a = 2073.1_8,
     &    b = 125.62_8,
     &    c = 3.6276_8,
     &    d = 0.043219_8
         WHERE (LANDMASK)
           SCHMIDT_CO2 = a + SSTT * (-b + SSTT * (c + SSTT * (-d)))
         ELSEWHERE
           SCHMIDT_CO2 = c0
         END WHERE
         return
        END subroutine  CSCHMIDT_CO2
       SUBROUTINE co2calc_row(mask, t, s, dic_in, ta_in, pt_in,
     &    sit_in, phlo, phhi, ph, xco2_in, atmpres, co2star,
     &    dco2star, pCO2surf, dpco2,istr,iend)
        implicit none
       real(kind=8),parameter::c0=0.0_8, c1=1.0_8, c10=10.0_8,
     &   c1000=1000.0_8, T0_Kelvin=273.16_8,rho_sw=4.1_8/3.996_8
       integer::istr,iend
       REAL(kind=8), PARAMETER :: xacc = 1e-10
       INTEGER(kind=4), PARAMETER :: maxit = 100
        REAL(kind=8), DIMENSION(istr:iend) ::
     &  k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff,
     &  bt, st, ft, dic, ta, pt, sit
         LOGICAL, DIMENSION(istr:iend), INTENT(IN) :: mask
         REAL(kind=8), DIMENSION(istr:iend), INTENT(IN) ::
     &    t,
     &    s,
     &    dic_in,
     &    ta_in,
     &    pt_in,
     &    sit_in,
     &    phlo,
     &    phhi,
     &    xco2_in,
     &    atmpres
       REAL(kind=8), DIMENSION(istr:iend), INTENT(OUT) ::
     &    ph,
     &    co2star,
     &    dco2star,
     &    pco2surf,
     &    dpco2
       INTEGER(kind=4) :: i
        REAL(kind=8) ::
     &    mass_to_vol,
     &    vol_to_mass,
     &    tk,
     &    is,
     &    scl,
     &    co2starair,
     &    tk100, tk1002, invtk, dlogtk, is2, sqrtis,
     &    s2, sqrts, s15, htotal2
        REAL(kind=8), DIMENSION(istr:iend) ::
     &    xco2,
     &    htotal,
     &    x1, x2
         ph          = c0
         co2star     = c0
         dco2star    = c0
         pCO2surf    = c0
         dpCO2       = c0
         IF (COUNT(mask) == 0) THEN
       RETURN
      END IF
       mass_to_vol = 1e6 * rho_sw
       vol_to_mass = c1 / mass_to_vol
        DO i = istr,iend
          IF (mask(i)) THEN
            dic(i)  = dic_in(i)  * vol_to_mass
            ta(i)   = ta_in(i)   * vol_to_mass
            pt(i)   = pt_in(i)   * vol_to_mass
            sit(i)  = sit_in(i)  * vol_to_mass
           xco2(i) = xco2_in(i) * 1e-6
            tk       = T0_Kelvin + t(i)
            tk100    = tk * 1e-2
            tk1002   = tk100 * tk100
            invtk    = c1 / tk
            dlogtk   = LOG(tk)
            is       = 19.924_8 * s(i) / (c1000 - 1.005_8 * s(i))
            is2      = is * is
            sqrtis   = SQRT(is)
            sqrts    = SQRT(s(i))
            s15      = s(i) ** 1.5_8
            s2       = s(i) ** 2
            scl      = s(i) / 1.80655_8
            ff(i) = EXP(-162.8301_8 + 218.2968_8/tk100 +
     &          90.9241_8*LOG(tk100) -
     &          1.47696_8*tk1002 + s(i)*(.025695_8 - .025225_8*tk100 +
     &          0.0049867_8*tk1002))
            k0(i) = EXP(93.4517_8/tk100 - 60.2409_8 + 
     &                       23.3585_8*LOG(tk100) +
     &          s(i)*(.023517_8 - 0.023656_8 * tk100 + 0.0047036_8 * 
     &                              tk1002))
          k1(i) = 10**(-1*(3670.7_8*invtk - 62.008_8 + 9.7944_8*dlogtk -
     &         0.0118_8*s(i) + 0.000116_8*s2))
          k2(i) = 10**(-1*(1394.7_8*invtk + 4.777_8 -
     &           0.0184_8*s(i) + 0.000118_8*s2))
          kb(i) = EXP((-8966.90_8 - 2890.53_8*sqrts - 77.942_8*s(i) +
     &          1.728_8*s15 - 0.0996_8*s2)*invtk +
     &          (148.0248_8 + 137.1942_8*sqrts + 1.62142_8*s(i)) +
     &          (-24.4344_8 - 25.085_8*sqrts - 0.2474_8*s(i)) *
     &          dlogtk + 0.053105_8*sqrts*tk)
          k1p(i) = EXP(-4576.752_8*invtk + 115.525_8 - 18.453_8 * 
     &                              dlogtk +
     &         (-106.736_8*invtk + 0.69171_8) * sqrts +
     &          (-0.65643_8*invtk - 0.01844_8) * s(i))
          k2p(i) = EXP(-8814.715_8*invtk + 172.0883_8 - 27.927_8 * 
     &                              dlogtk +
     &          (-160.340_8*invtk + 1.3566_8) * sqrts +
     &          (0.37335_8*invtk - 0.05778_8) * s(i))
          k3p(i) = EXP(-3070.75_8*invtk - 18.141_8 +
     &          (17.27039_8*invtk + 2.81197_8) * sqrts +
     &          (-44.99486_8*invtk - 0.09984_8) * s(i))
          ksi(i) = EXP(-8904.2_8*invtk + 117.385_8 - 19.334_8 * dlogtk +
     &         (-458.79_8*invtk + 3.5913_8) * sqrtis +
     &         (188.74_8*invtk - 1.5998_8) * is +
     &          (-12.1652_8*invtk + 0.07871_8) * is2 +
     &          LOG(1.0_8-0.001005_8*s(i)))
          kw(i) = EXP(-13847.26_8*invtk + 148.9652_8 - 23.6521_8 * 
     &                              dlogtk +
     &        (118.67_8*invtk - 5.977_8 + 1.0495_8 * dlogtk) *
     &         sqrts - 0.01615_8 * s(i))
          ks(i) = EXP(-4276.1_8*invtk + 141.328_8 - 23.093_8*dlogtk +
     &          (-13856*invtk + 324.57_8 - 47.986_8*dlogtk) *
     &          sqrtis +
     &          (35474*invtk - 771.54_8 + 114.723_8*dlogtk) * is -
     &          2698*invtk*is**1.5_8 + 1776*invtk*is2 +
     &          LOG(1.0_8 - 0.001005_8*s(i)))
          kf(i) = EXP(1590.2_8*invtk - 12.641_8 + 1.525_8*sqrtis +
     &         LOG(1.0_8 - 0.001005_8*s(i)) +
     &          LOG(1.0_8 + (0.1400_8/96.062_8)*(scl)/ks(i)))
          bt(i) = 0.000232_8 * scl/10.811_8
          st(i) = 0.14_8 * scl/96.062_8
          ft(i) = 0.000067_8 * scl/18.9984_8
          x1(i) = c10 ** (-phhi(i))
          x2(i) = c10 ** (-phlo(i))
        END IF
       END DO
        CALL drtsafe_row(mask, x1, x2, xacc, htotal,istr,
     &    iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
        DO i = istr,iend
          IF (mask(i)) THEN
          htotal2 = htotal(i) ** 2
          co2star(i) = dic(i) * htotal2 /
     &          (htotal2 + k1(i)*htotal(i) + k1(i)*k2(i))
          co2starair = xco2(i) * ff(i) * atmpres(i)
          dco2star(i) = co2starair - co2star(i)
          ph(i) = -LOG10(htotal(i))
          pCO2surf(i) = co2star(i) / ff(i)
          dpCO2(i)    = pCO2surf(i) - xco2(i) * atmpres(i)
          co2star(i)  = co2star(i) * mass_to_vol
          dco2star(i) = dco2star(i) * mass_to_vol
          pCO2surf(i) = pCO2surf(i) * 1e6
          dpCO2(i)    = dpCO2(i) * 1e6
        ELSE
          ph(i)       = c0
          co2star(i)  = c0
          dco2star(i) = c0
          pCO2surf(i) = c0
          dpCO2(i)    = c0
        END IF
       END DO
       END SUBROUTINE co2calc_row
       SUBROUTINE talk_row(mask, x, fn, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
        implicit none
       REAL(kind=8), PARAMETER :: xacc = 1e-10
       INTEGER(kind=4), PARAMETER :: maxit = 100
         integer(kind=4) :: istr,iend
        REAL(kind=8), DIMENSION(istr:iend) ::
     &  k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff,
     &  bt, st, ft, dic, ta, pt, sit
         real(kind=8),parameter::c1=1.0_8, c2=2.0_8, c3=3.0_8
        LOGICAL, DIMENSION(istr:iend), INTENT(IN) :: mask
        REAL(kind=8), DIMENSION(istr:iend), INTENT(IN) :: x
         REAL(kind=8), DIMENSION(istr:iend), INTENT(OUT) :: fn, df
         INTEGER(kind=4) :: i
          REAL(kind=8) ::
     &    x1, x2, x3, k12, k12p, k123p, a, a2, da, b, b2, db, c
      DO i = istr,iend
        IF (mask(i)) THEN
          x1 = x(i)
          x2 = x1 * x1
          x3 = x2 * x1
          k12 = k1(i) * k2(i)
          k12p = k1p(i) * k2p(i)
          k123p = k12p * k3p(i)
          a = x3 + k1p(i) * x2 + k12p * x1 + k123p
          a2 = a * a
          da = c3 * x2 + c2 * k1p(i) * x1 + k12p
          b = x2 + k1(i) * x1 + k12
          b2 = b * b
          db = c2 * x1 + k1(i)
          c = c1 + st(i)/ks(i)
          fn(i) = k1(i) * x1 * dic(i)/b +
     &          c2 * dic(i) * k12/b +
     &          bt(i)/(c1 + x1/kb(i)) +
     &          kw(i)/x1 +
     &          pt(i) * k12p * x1/a +
     &          c2 * pt(i) * k123p/a +
     &          sit(i)/(c1 + x1/ksi(i)) -
     &          x1/c -
     &          st(i)/(c1 + ks(i)/x1/c) -
     &          ft(i)/(c1 + kf(i)/x1) -
     &          pt(i) * x3/a -
     &          ta(i)
          df(i) = ((k1(i)*dic(i)*b) - k1(i)*x1*dic(i)*db)/b2 -
     &          c2 * dic(i) * k12 * db/b2 -
     &          bt(i)/kb(i)/(c1+x1/kb(i)) ** 2 -
     &          kw(i)/x2 +
     &          (pt(i) * k12p * (a - x1 * da))/a2 -
     &          c2 * pt(i) * k123p * da/a2 -
     &          sit(i)/ksi(i)/(c1+x1/ksi(i)) ** 2 -
     &          c1/c +
     &          st(i) * (c1 + ks(i)/x1/c)**(-2) * (ks(i)/c/x2) +
     &          ft(i) * (c1 + kf(i)/x1)**(-2) * kf(i)/x2 -
     &          pt(i) * x2 * (c3 * a - x1 * da)/a2
        END IF
        END DO
      END SUBROUTINE talk_row
        SUBROUTINE drtsafe_row(mask_in, x1, x2, xacc, soln,istr,
     &   iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
        implicit none
       INTEGER(kind=4), PARAMETER :: maxit = 100
       integer::istr,iend
        REAL(kind=8), DIMENSION(istr:iend) ::
     &  k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff,
     &  bt, st, ft, dic, ta, pt, sit
        real(kind=8),parameter::c0=0.0_8, c2=2.0_8
        LOGICAL, DIMENSION(istr:iend), INTENT(IN) :: mask_in
        REAL(kind=8), DIMENSION(istr:iend), INTENT(IN) :: x1, x2
        REAL(kind=8), INTENT(IN) :: xacc
       REAL(kind=8), DIMENSION(istr:iend), INTENT(OUT) :: soln
       LOGICAL :: leave_bracket, dx_decrease
       LOGICAL, DIMENSION(istr:iend) :: mask
       INTEGER(kind=4) ::  i, it
       REAL(kind=8) :: temp
       REAL(kind=8), DIMENSION(istr:iend) :: xlo, xhi, flo, fhi, f,
     &     df, dxold, dx
       mask = mask_in
       CALL talk_row(mask, x1, flo, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
       CALL talk_row(mask, x2, fhi, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
       DO i = istr,iend
         IF (mask(i)) THEN
           IF (flo(i) .LT. c0) THEN
              xlo(i) = x1(i)
              xhi(i) = x2(i)
           ELSE
             xlo(i) = x2(i)
             xhi(i) = x1(i)
             temp = flo(i)
             flo(i) = fhi(i)
             fhi(i) = temp
          END IF
          soln(i) = 0.5_8 * (xlo(i) + xhi(i))
          dxold(i) = ABS(xlo(i) - xhi(i))
          dx(i) = dxold(i)
         END IF
        END DO
        CALL talk_row(mask, soln, f, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
       DO it = 1,maxit
        DO i = istr,iend
           IF (mask(i)) THEN
             leave_bracket = ((soln(i)-xhi(i))*df(i)-f(i)) *
     &             ((soln(i)-xlo(i))*df(i)-f(i)) .GE. 0
             dx_decrease = ABS(c2 * f(i)) .LE. ABS(dxold(i) * df(i))
             IF (leave_bracket .OR. .NOT. dx_decrease) THEN
                dxold(i) = dx(i)
                dx(i) = 0.5_8 * (xhi(i) - xlo(i))
                soln(i) = xlo(i) + dx(i)
                IF (xlo(i) .EQ. soln(i)) mask(i) = .FALSE.
             ELSE
                dxold(i) = dx(i)
                dx(i) = -f(i) / df(i)
                temp = soln(i)
                soln(i) = soln(i) + dx(i)
                IF (temp .EQ. soln(i)) mask(i) = .FALSE.
             END IF
             IF (ABS(dx(i)) .LT. xacc) mask(i) = .FALSE.
          END IF
        END DO
        IF (.NOT. ANY(mask)) RETURN
        CALL talk_row(mask, soln, f, df,istr,iend, k0, k1, k2,
     &   kw, kb, ks, kf, k1p, k2p, k3p, ksi,
     &    ff, bt, st, ft, dic, ta, pt, sit)
        DO i = istr,iend
          IF (mask(i)) THEN
             IF (f(i) .LT. c0) THEN
                xlo(i) = soln(i)
                flo(i) = f(i)
             ELSE
                xhi(i) = soln(i)
                fhi(i) = f(i)
             END IF
          END IF
        END DO
       END DO
       END SUBROUTINE drtsafe_row
