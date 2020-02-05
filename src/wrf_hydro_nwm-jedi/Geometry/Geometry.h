/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_JEDI_GEOMETRY_GEOMETRY_H_
#define WRF_HYDRO_JEDI_GEOMETRY_GEOMETRY_H_

#include <ostream>
#include <string>

#include "eckit/mpi/Comm.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace wrf_hydro_jedi {
  class GeometryIterator;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_jedi {

  // Geometry class
  class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry> {
   public:
    static const std::string classname() {return "wrf_hydro_jedi::Geometry";}

    explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
    Geometry(const Geometry &);
    ~Geometry();

    const eckit::mpi::Comm & getComm() const {return comm_;}

    // //These are needed for the GeometryIterator Interface
    // GeometryIterator begin() const;
    // GeometryIterator end() const;

   private:
    void print(std::ostream &) const;

    int keyGeom_;
    const eckit::mpi::Comm & comm_;
  };
}  // namespace wrf_hydro_jedi

#endif  // WRF_HYDRO_JEDI_GEOMETRY_GEOMETRY_H_
