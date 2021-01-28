/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm-jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm-jedi/Geometry/GeometryFortran.h"
#include "wrf_hydro_nwm_jedi/GeometryIterator/GeometryIterator.h"

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point2.h"

#include "oops/util/abor1_cpp.h"

namespace wrf_hydro_nwm_jedi {


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

  bool GeometryIterator::operator!=(const GeometryIterator &) const {
    util::abor1_cpp("GeometryIterator::operator!=() needs to be implemented.",
                    __FILE__, __LINE__);
    return false;
  }

// ----------------------------------------------------------------------------

  GeometryIterator& GeometryIterator::operator++() {
    util::abor1_cpp("GeometryIterator::operator++() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void GeometryIterator::print(std::ostream  & os) const {
    util::abor1_cpp("GeometryIterator::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "(TODO, print diagnostic info about the GeometryIterator here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------
  eckit::geometry::Point2 GeometryIterator::operator*() const {
    util::abor1_cpp("GeometryIterator::operator*() needs to be implemented.",
                    __FILE__, __LINE__);
    return eckit::geometry::Point2(0.0, 0.0);
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
