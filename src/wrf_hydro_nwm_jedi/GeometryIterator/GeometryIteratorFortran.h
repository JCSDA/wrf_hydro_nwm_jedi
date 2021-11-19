/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_GEOMETRYITERATOR_GEOMETRYITERATORFORTRAN_H_
#define WRF_HYDRO_NWM_JEDI_GEOMETRYITERATOR_GEOMETRYITERATORFORTRAN_H_

#include "wrf_hydro_nwm_jedi/Utilities/interface.h"

namespace wrf_hydro_nwm_jedi {


  extern "C" {
    void wrf_hydro_nwm_jedi_geom_iter_setup_f90(F90iter &, const F90geom &,
                                  const int &, const int &);
    void wrf_hydro_nwm_jedi_geom_iter_clone_f90(F90iter &, const F90iter &);
    void wrf_hydro_nwm_jedi_geom_iter_delete_f90(F90iter &);
    void wrf_hydro_nwm_jedi_geom_iter_equals_f90(const F90iter &, const F90iter&, int &);
    void wrf_hydro_nwm_jedi_geom_iter_current_f90(const F90iter &, double &, double &);
  }
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_GEOMETRYITERATOR_GEOMETRYITERATORFORTRAN_H_
