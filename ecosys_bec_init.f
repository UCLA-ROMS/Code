      subroutine ecosys_init()
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
      integer::i,j,k,m,nn
      sinking_particle_POC=0.0_8
      sinking_particle_P_CaCO3=0.0_8
      sinking_particle_P_SiO2=0.0_8
      sinking_particle_dust=0.0_8
      sinking_particle_P_iron =0.0_8
       parm_Fe_bioavail    = 0.02_8
       parm_prod_dissolve  = 0.67_8
       parm_o2_min         = 1.0_8
       parm_Rain_CaCO3     = 0.07_8
       parm_Rain_SiO2      = 0.03_8
       parm_kappa_nitrif_no2 = 0.33333_8 * dps
       parm_n2o_prod       = 0.95_8
       parm_n2_prod        = 0.33333_8 * dps
       parm_kappa_nitrif   = 0.06_8 * dps
       parm_nitrif_par_lim = 4.0_8
       parm_POC_flux_ref   = 2.D-3
       parm_rest_prod_tau  = 30.0_8 * dps
       parm_rest_prod_z_c  = 75._8
       parm_z_umax_0       = 2.5_8 * dps
       parm_diat_umax_0    = 1.95_8 * dps
       parm_z_mort_0       = 0.08_8 * dps
       parm_z_mort2_0      = 0.42_8 * dps
       parm_sd_remin_0     = 0.01_8 * dps
       parm_sp_kNO3        = 0.5_8
       parm_diat_kNO3      = 2.5_8
       parm_sp_kNH4        = 0.01_8
       parm_diat_kNH4      = 0.1_8
       parm_sp_kFe         = 0.035D-3
       parm_diat_kFe       = 0.08D-3
       parm_diat_kSiO3     = 1.0_8
       parm_sp_kPO4        = 0.01_8
       parm_diat_kPO4      = 0.1_8
       parm_z_grz          = 1.05_8
       parm_alphaChl       = 0.3_8 * dps
       parm_labile_ratio   = 0.70_8
       parm_alphaDiaz      = 0.036_8 * dps
       parm_diaz_umax_0    = 1.2_8 * dps
       t_remin_sed_poc   = 0.003_8 * dps
       t_remin_sed_caco3 = 0.003_8 * dps
       t_remin_sed_si    = 0.003_8 * dps
        gQsi_0        = 0.137_8
        gQsi_coef     = 2.5_8
        gQsi_max      = 0.685_8
        lrest_po4      = .FALSE.
        lrest_no3      = .FALSE.
        lrest_sio3     = .FALSE.
        lsource_sink    = .TRUE.
        lflux_gas_o2    = .TRUE.
        lflux_gas_co2   = .TRUE.
        liron_flux      =.TRUE.
        ldust_flux      =.TRUE.
         where(rmask ==1)
              landmask=.TRUE.
            elsewhere
               landmask=.FALSE.
            endwhere
         ifrac=0.0_8
         press=1._8
      IF (lrest_po4 .OR. lrest_no3 .OR. lrest_sio3) THEN
      END IF
       pH_hist              = c0
       pCO2sw               = c0
       pCO2air              = c0
       PAR                  = c0
       PARinc               = c0
       WS_HIST              = c0
       XKW_HIST             = c0
       AP_HIST              = c0
       SCHMIDT_O2_HIST      = c0
       O2SAT_HIST           = c0
       FG_O2_HIST           = c0
       SCHMIDT_CO2_HIST     = c0
       fg_n2o_hist          = c0
       fg_n2_hist           = c0
       CO2STAR_HIST         = c0
       DCO2STAR_HIST        = c0
       FG_CO2_HIST          = c0
       IRON_FLUX_HIST       = c0
       PO4_RESTORE_HIST     = c0
       NO3_RESTORE_HIST     = c0
       SiO3_RESTORE_HIST    = c0
       PO4STAR_HIST         = c0
       POC_FLUX_IN_HIST     = c0
       POC_PROD_HIST        = c0
       POC_REMIN_HIST       = c0
       CaCO3_FLUX_IN_HIST   = c0
       CaCO3_PROD_HIST      = c0
       CaCO3_REMIN_HIST     = c0
       SiO2_FLUX_IN_HIST    = c0
       SiO2_PROD_HIST       = c0
       SiO2_REMIN_HIST      = c0
       dust_FLUX_IN_HIST    = c0
       dust_REMIN_HIST      = c0
       P_iron_FLUX_IN_HIST  = c0
       P_iron_PROD_HIST     = c0
       P_iron_REMIN_HIST    = c0
       graze_sp_HIST        = c0
       graze_diat_HIST      = c0
       graze_tot_HIST       = c0
       sp_loss_HIST         = c0
       diat_loss_HIST       = c0
       zoo_loss_HIST        = c0
       sp_agg_HIST          = c0
       diat_agg_HIST        = c0
       photoC_sp_HIST       = c0
       photoC_diat_HIST     = c0
       tot_prod_HIST        = c0
       no3_v_sp_HIST        = c0
       nh4_v_sp_HIST        = c0
       no3_v_diat_HIST      = c0
       nh4_v_diat_HIST      = c0
       DOC_prod_HIST        = c0
       DOC_remin_HIST       = c0
       Fe_scavenge_HIST     = c0
       sp_N_lim_HIST        = c0
       sp_Fe_lim_HIST       = c0
       sp_PO4_lim_HIST      = c0
       sp_light_lim_HIST    = c0
       diat_N_lim_HIST      = c0
       diat_Fe_lim_HIST     = c0
       diat_PO4_lim_HIST    = c0
       diat_SiO3_lim_HIST   = c0
       diat_light_lim_HIST  = c0
       CaCO3_form_HIST      = c0
       diaz_Nfix_HIST       = c0
       graze_diaz_HIST      = c0
       diaz_loss_HIST       = c0
       photoC_diaz_HIST     = c0
       diaz_P_lim_HIST      = c0
       diaz_Fe_lim_HIST     = c0
       diaz_light_lim_HIST  = c0
       Fe_scavenge_rate_HIST = c0
       DON_prod_HIST     = c0
       DON_remin_HIST    = c0
       DOFe_prod_HIST    = c0
       DOFe_remin_HIST   = c0
       DOP_prod_HIST     = c0
       DOP_remin_HIST    = c0
       bSi_form_HIST     = c0
       photoFe_diat_HIST    = c0
       photoFe_diaz_HIST    = c0
       photoFe_sp_HIST      = c0
       nitrif_hist       = c0
       ammox_hist = c0
       denitr_hist = c0
       n2o_prod_hist = c0
       n2_prod_hist = c0
       pH_avg              = c0
       pCO2_avg          = c0
       pCO2air_avg         = c0
       PAR_avg             = c0
       PARinc_avg          = c0
       WS_AVG              = c0
       XKW_AVG             = c0
       AP_AVG              = c0
       SCHMIDT_O2_AVG      = c0
       O2SAT_AVG           = c0
       FG_O2_AVG           = c0
       SCHMIDT_CO2_AVG     = c0
       fg_n2o_avg          = c0
       fg_n2_avg           = c0
       CO2STAR_AVG         = c0
       DCO2STAR_AVG        = c0
       FG_CO2_AVG          = c0
       IRON_FLUX_AVG       = c0
       PO4_RESTORE_AVG     = c0
       NO3_RESTORE_AVG     = c0
       SiO3_RESTORE_AVG    = c0
       PO4STAR_AVG         = c0
       POC_FLUX_IN_AVG     = c0
       POC_PROD_AVG        = c0
       POC_REMIN_AVG       = c0
       CaCO3_FLUX_IN_AVG   = c0
       CaCO3_PROD_AVG      = c0
       CaCO3_REMIN_AVG     = c0
       SiO2_FLUX_IN_AVG    = c0
       SiO2_PROD_AVG       = c0
       SiO2_REMIN_AVG      = c0
       dust_FLUX_IN_AVG    = c0
       dust_REMIN_AVG      = c0
       P_iron_FLUX_IN_AVG  = c0
       P_iron_PROD_AVG     = c0
       P_iron_REMIN_AVG    = c0
       graze_sp_AVG        = c0
       graze_diat_AVG      = c0
       graze_tot_AVG       = c0
       sp_loss_AVG         = c0
       diat_loss_AVG       = c0
       zoo_loss_AVG        = c0
       sp_agg_AVG          = c0
       diat_agg_AVG        = c0
       photoC_sp_AVG       = c0
       photoC_diat_AVG     = c0
       tot_prod_AVG        = c0
       no3_v_sp_AVG        = c0
       nh4_v_sp_AVG        = c0
       no3_v_diat_AVG      = c0
       nh4_v_diat_AVG      = c0
       DOC_prod_AVG        = c0
       DOC_remin_AVG       = c0
       Fe_scavenge_AVG     = c0
       sp_N_lim_AVG        = c0
       sp_Fe_lim_AVG       = c0
       sp_PO4_lim_AVG      = c0
       sp_light_lim_AVG    = c0
       diat_N_lim_AVG      = c0
       diat_Fe_lim_AVG     = c0
       diat_PO4_lim_AVG    = c0
       diat_SiO3_lim_AVG   = c0
       diat_light_lim_AVG  = c0
       CaCO3_form_AVG      = c0
       diaz_Nfix_AVG       = c0
       graze_diaz_AVG      = c0
       diaz_loss_AVG       = c0
       photoC_diaz_AVG     = c0
       diaz_P_lim_AVG      = c0
       diaz_Fe_lim_AVG     = c0
       diaz_light_lim_AVG  = c0
       Fe_scavenge_rate_AVG = c0
       DON_prod_AVG      = c0
       DON_remin_AVG     = c0
       DOFe_prod_AVG     = c0
       DOFe_remin_AVG    = c0
       DOP_prod_AVG      = c0
       DOP_remin_AVG     = c0
       bSi_form_AVG      = c0
       photoFe_diat_AVG    = c0
       photoFe_diaz_AVG    = c0
       photoFe_sp_AVG      = c0
       nitrif_avg       = c0
       ammox_avg = c0
       denitr_avg = c0
       n2o_prod_avg = c0
       n2_prod_avg = c0
       call init_scalars_becflux()
       return
       end
