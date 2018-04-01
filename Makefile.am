lib_LTLIBRARIES = libpiaacmcsimul.la

libpiaacmcsimul_la_SOURCES = PIAACMCsimul.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_initpiaacmcconf.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_loadsavepiaacmcconf.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_init.c 
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_free.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_mkFPM_zonemap.c

libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec_compute_image.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec_optimize_lyot_stop_position.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec_optimize_fpmtransmission.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec_computePSF_no_fpm.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec_optimize_PIAA_shapes.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec_optimize_lyot_stops_shapes_positions.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec_multizone_fpm_calib.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec_optimize_fpm_zones.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_exec_optimize_PIAA_shapes_fpmtransm.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_measure_transm_curve.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_eval_poly_design.c

libpiaacmcsimul_la_SOURCES += PIAACMCsimul_mkFocalPlaneMask.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_rings2sectors.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_makePIAAshapes.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_mkPIAAMshapes_from_RadSag.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_run.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_CA2propCubeInt.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_computePSF.c
libpiaacmcsimul_la_SOURCES += PIAACMC_FPM_process.c
libpiaacmcsimul_la_SOURCES += PIAACMC_FPMresp_resample.c
libpiaacmcsimul_la_SOURCES += PIAACMC_FPMresp_rmzones.c
libpiaacmcsimul_la_SOURCES += PIAACMC_f_evalmask.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_achromFPMsol_eval_zonezderivative.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_achromFPMsol_eval.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_optimizeLyotStop.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_mkLyotMask.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_mkSimpleLyotStop.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_init_geomPIAA_rad.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_load2DRadialApodization.c
libpiaacmcsimul_la_SOURCES += PIAACMCsimul_geomProp.c

libpiaacmcsimul_la_SOURCES += PIAACMCsimul.h

EXTRA_DIST = scripts

AM_CPPFLAGS = -I@abs_top_srcdir@/src -fopenmp
