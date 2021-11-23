/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_INCREMENT_INCREMENTFORTRAN_H_
#define WRF_HYDRO_NWM_JEDI_INCREMENT_INCREMENTFORTRAN_H_

#include <string>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

namespace wrf_hydro_nwm_jedi {

  typedef int F90inc;
  typedef int F90state;

  extern "C" {
    void wrf_hydro_nwm_jedi_increment_create_f90(
        F90inc &, const F90geom &, const oops::Variables &);
    void wrf_hydro_nwm_jedi_increment_read_file_f90(
        const F90geom &, const F90inc &,
        const eckit::Configuration * const *,
        const util::DateTime * const *);
    void wrf_hydro_nwm_jedi_increment_write_file_f90(
        const F90geom &,
        const F90inc &,
        const eckit::Configuration * const *,
        const util::DateTime * const *);
    void wrf_hydro_nwm_jedi_increment_diff_incr_f90(
        F90state &, const F90state &, const F90state &);
    void wrf_hydro_nwm_jedi_increment_create_from_other_f90(F90inc &, const F90inc &);
    void wrf_hydro_nwm_jedi_increment_sub_f90(const F90inc &, const F90inc &);
    void wrf_hydro_nwm_jedi_increment_add_f90(const F90inc &, const F90inc &);
    void wrf_hydro_nwm_jedi_increment_mul_f90(const F90inc &, const float &);

    void wrf_hydro_nwm_jedi_increment_dot_prod_f90(const F90inc &, const F90inc &, double &);
    void wrf_hydro_nwm_jedi_increment_random_f90(const F90inc &);

    void wrf_hydro_nwm_jedi_increment_zero_f90(const F90inc &);
    void wrf_hydro_nwm_jedi_increment_ones_f90(const F90inc &);
    void wrf_hydro_nwm_jedi_increment_dirac_f90(const F90inc &,
         const eckit::Configuration * const *);
    double wrf_hydro_nwm_jedi_increment_rms_f90(const F90inc &);
    void wrf_hydro_nwm_jedi_increment_axpy_f90(const F90inc &, const double &, const F90inc &);
    void wrf_hydro_nwm_jedi_increment_accumul_f90(
        const F90inc &, const double &, const F90state &);
    void wrf_hydro_nwm_jedi_increment_schur_f90(const F90inc &, const F90inc &);
    void wrf_hydro_nwm_jedi_increment_print_f90(const F90inc &, char *string);
    /* void wrf_hydro_nwm_jedi_state_get_mean_stddev_f90(
            const F90state &, int nf, float pstat[][1]); */
    void wrf_hydro_nwm_jedi_increment_copy_f90(const F90inc &, const F90inc &);
    /* void wrf_hydro_nwm_jedi_geometry_clone_f90(F90geom &, const F90geom &); */
    /* void wrf_hydro_nwm_jedi_geometry_delete_f90(F90geom &); */
    void wrf_hydro_nwm_jedi_increment_set_atlas_f90(const F90inc &, const F90geom &,
        const oops::Variables &, atlas::field::FieldSetImpl *);
    void wrf_hydro_nwm_jedi_increment_to_atlas_f90(const F90inc &, const F90geom &,
        const oops::Variables &, atlas::field::FieldSetImpl *);
    void wrf_hydro_nwm_jedi_increment_from_atlas_f90(const F90inc &, const oops::Variables &,
        atlas::field::FieldSetImpl *);
    void wrf_hydro_nwm_jedi_increment_getpoint_f90(const F90inc &, const F90iter &, double &,
        const int &);
  }
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_INCREMENT_INCREMENTFORTRAN_H_
