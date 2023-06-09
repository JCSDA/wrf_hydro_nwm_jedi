/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_GEOMETRY_GEOMETRYFORTRAN_H_
#define WRF_HYDRO_NWM_JEDI_GEOMETRY_GEOMETRYFORTRAN_H_

#include "eckit/config/Configuration.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/Variables.h"

namespace wrf_hydro_nwm_jedi {

  typedef int F90geom;

  extern "C" {
    void wrf_hydro_nwm_jedi_geometry_setup_f90(F90geom &,
                                    const eckit::Configuration * const *);

    void wrf_hydro_nwm_jedi_geometry_set_lonlat_f90(const F90geom &, atlas::field::FieldSetImpl *,
                                     const bool &);
    void wrf_hydro_nwm_jedi_geometry_set_functionspace_pointer_f90(const F90geom &,
                                                    atlas::functionspace::FunctionSpaceImpl *,
                                                    atlas::functionspace::FunctionSpaceImpl *);
    void wrf_hydro_nwm_jedi_geometry_fill_extra_fields_f90(const F90geom &,
                                                    atlas::field::FieldSetImpl *);

    void wrf_hydro_nwm_jedi_geometry_clone_f90(F90geom &, const F90geom &);
    void wrf_hydro_nwm_jedi_geometry_delete_f90(F90geom &);
    void wrf_hydro_nwm_jedi_geometry_info_f90(F90geom, float *, float *, int *, int *, int *);
    void wrf_hydro_nwm_jedi_geometry_get_nn_f90(F90geom, float, float, int, int);
    void wrf_hydro_nwm_jedi_geoval_levels_f90(const F90geom &,
                                   const oops::Variables &,
                                   const std::size_t &,
                                   std::size_t &);
  }
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_GEOMETRY_GEOMETRYFORTRAN_H_
