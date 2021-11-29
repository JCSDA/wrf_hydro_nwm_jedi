ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_forecast.x
#                        SOURCES Forecast.cc
#                        LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_hofx.x
                        SOURCES HofX.cc
                        LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_letkf.x
                       SOURCES LETKF.cc
                       LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_error_covariance_training.x
                        SOURCES ErrorCovarianceTraining.cc
                        LIBS    wrf_hydro_nwm_jedi
                                saber )

ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_dirac.x
                        SOURCES Dirac.cc
                        LIBS    wrf_hydro_nwm_jedi
                                saber )

ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_var.x
                        SOURCES Var.cc
                        LIBS    wrf_hydro_nwm_jedi
                                saber )
