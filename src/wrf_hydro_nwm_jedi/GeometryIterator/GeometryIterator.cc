/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/GeometryIterator/GeometryIterator.h"
#include "wrf_hydro_nwm_jedi/GeometryIterator/GeometryIteratorFortran.h"

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point2.h"

#include "oops/util/abor1_cpp.h"

namespace wrf_hydro_nwm_jedi {


// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const GeometryIterator& iter) {
  wrf_hydro_nwm_jedi_geom_iter_clone_f90(keyIter_, iter.toFortran());
}

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const Geometry& geom,
                                       const int & iindex, const int & jindex) {
  wrf_hydro_nwm_jedi_geom_iter_setup_f90(keyIter_, geom.toFortran(), iindex, jindex);
}


// -----------------------------------------------------------------------------

GeometryIterator::~GeometryIterator() {
  wrf_hydro_nwm_jedi_geom_iter_delete_f90(keyIter_);
}

// ----------------------------------------------------------------------------

  bool GeometryIterator::operator==(const GeometryIterator & other) const {
    int equals = 0;
    wrf_hydro_nwm_jedi_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
    return (equals == 1);
  }

// ----------------------------------------------------------------------------

  bool GeometryIterator::operator!=(const GeometryIterator & other) const {
    int equals = 0;
    wrf_hydro_nwm_jedi_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
    return (equals == 0);
  }

// ----------------------------------------------------------------------------

  GeometryIterator& GeometryIterator::operator++() {
    util::abor1_cpp("GeometryIterator::operator++() needs to be implemented.",
                    __FILE__, __LINE__);
    return *this;
  }

// ----------------------------------------------------------------------------

  void GeometryIterator::print(std::ostream  & os) const {
  double lat, lon;
  wrf_hydro_nwm_jedi_geom_iter_current_f90(keyIter_, lon, lat);
  os << "GeometryIterator, lat/lon: " << lat << " / " << lon << std::endl;
  }

// ----------------------------------------------------------------------------
  eckit::geometry::Point2 GeometryIterator::operator*() const {
    util::abor1_cpp("GeometryIterator::operator*() needs to be implemented.",
                    __FILE__, __LINE__);
    return eckit::geometry::Point2(0.0, 0.0);
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
