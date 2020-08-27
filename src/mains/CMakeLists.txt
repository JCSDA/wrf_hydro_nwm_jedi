ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_hofxnomodel.x
                        SOURCES HofXNoModel.cc
                        LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_make_obs.x
#                         SOURCES MakeObs.cc
#                         LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_forecast.x
#                        SOURCES Forecast.cc
#                        LIBS    wrf_hydro_nwm_jedi )

#ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_hofx.x
#                       SOURCES HofX.cc
#                       LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_letkf.x
#                        SOURCES LETKF.cc
#                        LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_staticbinit.x
#                         SOURCES StaticBInit.cc
#                         LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_executable( TARGET  wrf_hydro_nwm_jedi_var.x
                        SOURCES Var.cc
                        LIBS    wrf_hydro_nwm_jedi )
