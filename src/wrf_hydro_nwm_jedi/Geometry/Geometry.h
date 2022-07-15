/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_GEOMETRY_GEOMETRY_H_
#define WRF_HYDRO_NWM_JEDI_GEOMETRY_GEOMETRY_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "wrf_hydro_nwm_jedi/Geometry/GeometryFortran.h"
#include "wrf_hydro_nwm_jedi/GeometryIterator/GeometryIterator.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace wrf_hydro_nwm_jedi {
  class GeometryIterator;
}
namespace oops {
  class Variables;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

  /// Geometry class for wrf_hydro_nwm - jedi
  class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry> {
   public:
    static const std::string classname() {return "wrf_hydro_nwm_jedi::Geometry";}

    explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
    Geometry(const Geometry &);
    ~Geometry();

    // These are needed for the GeometryIterator Interface
    GeometryIterator begin() const;
    GeometryIterator end() const;

    const F90geom & toFortran() const {return keyGeom_;}
    const eckit::mpi::Comm & getComm() const {return comm_;}
    std::vector<size_t> variableSizes(const oops::Variables &) const;
    std::vector<double> verticalCoord(std::string &) const {return {};}

    const atlas::FunctionSpace & functionSpace() const {return functionSpaceIncludingHalo_;}
    atlas::FunctionSpace & functionSpace() {return functionSpaceIncludingHalo_;}
    const atlas::FieldSet & extraFields() const {return extraFields_;}
    atlas::FieldSet & extraFields() {return extraFields_;}
    void latlon(std::vector<double> &, std::vector<double> &, const bool) const;

   private:
    void print(std::ostream &) const;
    Geometry & operator=(const Geometry &);
    F90geom keyGeom_;
    const eckit::mpi::Comm & comm_;
    atlas::FunctionSpace functionSpace_;
    atlas::FunctionSpace functionSpaceIncludingHalo_;
    atlas::FieldSet extraFields_;
  };
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_GEOMETRY_GEOMETRY_H_
