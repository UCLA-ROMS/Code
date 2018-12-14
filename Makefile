# makefile for ROMS model 

.SUFFIXES: .o .f
.f.o:
	$(CFT) -c $(FFLAGS) $*.f -o $*.o

CFT = mpif90 -pc80 -align dcommons -auto -stack_temps
LDR = $(CFT)

LARGE_MEM_FLAG = -mcmodel=medium -shared-intel
FFLAGS = -O1 $(LARGE_MEM_FLAG) -mp1 -g -module . 
LCDF = -lnetcdf

SRCS = main.f		step2D_FB.f	read_inp.f\
	set_weights.f	set_scoord.f	init_scalars.f	init_arrays.f\
	setup_grid1.f	setup_grid2.f\
	set_nudgcof.f\
	prsgrd32AC1.f	pre_step3d4S.f	step3d_uv1.f	step3d_uv2.f\
	step3d_t_ISO.f	set_depth.f	omega.f\
	visc3d_S.f	t3dmix_GP.f\
	zetabc.f	u2dbc_im.f	v2dbc_im.f\
	u3dbc_im.f	v3dbc_im.f	t3dbc_im.f	exchange.f\
	rho_eos.f	ab_ratio.f	alfabeta.f\
	lmd_vmix.f      lmd_kpp.f	lmd_swr_frac.f\
        diag.f		timers.f	wvlcty.f	grid_stiffness.f\
        lenstr.f	setup_kwds.f\
        get_date.f	ext_copy_prv2shr.f\
	mpi_setup.f	mpi_exchange8WA.f\
	init_scalars_bec.f ecosys_bec.f	ecosys_bec_init.f \
	det_srflx_dailyavg.f\
	init_scalars_becflux.f set_bec_flux_avg.f\
	init_arrays_physflux.f init_scalars_physflux.f  set_phys_flux_avg.f\
	bio_diag.f\
	nf_fread.f	nf_read_bry.f	set_cycle.f	checkdims.f\
	insert_node.f	closecdf.f	put_global_atts.f\
	get_grid.f	get_init.f	wrt_grid.f\
	def_rst.f	wrt_rst.f	def_his.f	wrt_his.f\
	set_avg.f	wrt_avg.f\
	get_forces.f	get_bry_all.f\
	get_sss.f\
	def_bec_flux.f 	get_dust.f\
	get_iron.f	get_atm_pco2.f\
	wrt_bec_flux_avg.f wrt_bec_flux_his.f\
	def_phys_flux.f wrt_phys_flux_avg.f  wrt_phys_flux_his.f\
	mod_phys_flux.f\
	bulk_flux.f       get_bulk_rad.f     get_bulk_wnd.f\
	get_bulk_prec.f   get_bulk_tra.f     wrt_bulk_diags_avg.f\
	wrt_bulk_diags_his.f def_bulk_diags.f set_bulk_diags_avg.f  

OBJS = $(SRCS:.f=.o)

roms: phys_flux.mod $(OBJS)
	$(LDR) $(FFLAGS) -o roms $(OBJS) $(LCDF) $(LMPI)

phys_flux.mod: mod_phys_flux.o

clean:
	/bin/rm -rf *.o roms *~

