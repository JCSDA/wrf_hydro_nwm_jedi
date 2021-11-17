/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_GEOMETRYITERATOR_GEOMETRYITERATOR_H_
#define WRF_HYDRO_NWM_JEDI_GEOMETRYITERATOR_GEOMETRYITERATOR_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
  namespace geometry {
    class Point2;
  }
}

namespace wrf_hydro_nwm_jedi {
  class Geometry;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

  // Geometry class
  class GeometryIterator : public util::Printable,
                           private util::ObjectCounter<GeometryIterator> {
   public:
    static const std::string classname() {return "wrf_hydro_nwm_jedi::GeometryIterator";}

    GeometryIterator(const GeometryIterator &);
    explicit GeometryIterator(const Geometry & geom,
                              const int & iindex = 1, const int & jindex = 1);
    ~GeometryIterator();

    bool operator==(const GeometryIterator &) const;
    bool operator!=(const GeometryIterator &) const;
    GeometryIterator& operator++();
    eckit::geometry::Point2 operator*() const;

   private:
    void print(std::ostream &) const;
  };
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_GEOMETRYITERATOR_GEOMETRYITERATOR_H_
