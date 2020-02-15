/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
// #include "wrf_hydro_nwm_jedi/GeometryIterator/GeometryIterator.h"
#include "wrf_hydro_nwm_jedi/Geometry/GeometryFortran.h"

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"

namespace wrf_hydro_nwm_jedi {

// ----------------------------------------------------------------------------

  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm)
    : comm_(comm) {
    const eckit::Configuration * configc = &conf;
    wrf_hydro_nwm_jedi_geometry_setup_f90(keyGeom_, &configc);
  }

// ----------------------------------------------------------------------------

  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_) {
    wrf_hydro_nwm_jedi_geometry_clone_f90(keyGeom_, other.keyGeom_);
  }

// ----------------------------------------------------------------------------

  Geometry::~Geometry() {
    wrf_hydro_nwm_jedi_geometry_delete_f90(keyGeom_);
  }

// ----------------------------------------------------------------------------

  void Geometry::print(std::ostream & os) const {
    util::abor1_cpp("Geometry::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "Geometry: "
       << "(TODO, print diagnostic info about the geometry here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

  // GeometryIterator Geometry::begin() const {
  //   util::abor1_cpp("Geometry::begin() needs to be implemented.",
  //                   __FILE__, __LINE__);
  //   return GeometryIterator();
  // }

// ----------------------------------------------------------------------------

  // GeometryIterator Geometry::end() const {
  //   util::abor1_cpp("Geometry::end() needs to be implemented.",
  //                   __FILE__, __LINE__);
  //   return GeometryIterator();
  // }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
