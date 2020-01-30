/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM-JEDI_MODELAUX_MODELAUXINCREMENT_H_
#define WRF_HYDRO_NWM-JEDI_MODELAUX_MODELAUXINCREMENT_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace wrf_hydro_nwm-jedi {
  class Geometry;
}

//-----------------------------------------------------------------------------

namespace wrf_hydro_nwm-jedi {

  // ModelAuxControl class
  class ModelAuxIncrement : public util::Printable,
                            private util::ObjectCounter<ModelAuxIncrement> {
   public:
    static const std::string classname() {return "wrf_hydro_nwm-jedi::ModelAuxIncrement";}

    ModelAuxIncrement(const ModelAuxIncrement &, const eckit::Configuration &);
    ModelAuxIncrement(const ModelAuxIncrement &, const bool);
    ModelAuxIncrement(const Geometry &, const eckit::Configuration &);
    ~ModelAuxIncrement();

    // Linear algebra operators
    void zero();
    ModelAuxIncrement & operator*=(const double);
    ModelAuxIncrement & operator+=(const ModelAuxIncrement &);
    ModelAuxIncrement & operator-=(const ModelAuxIncrement &);
    double norm() const;
    void axpy(const double, const ModelAuxIncrement &);
    double dot_product_with(const ModelAuxIncrement &) const;

   private:
    void print(std::ostream &) const;
  };
}  // namespace wrf_hydro_nwm-jedi
#endif  // WRF_HYDRO_NWM-JEDI_MODELAUX_MODELAUXINCREMENT_H_
