/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM-JEDI_COVARIANCE_COVARIANCE_H_
#define WRF_HYDRO_NWM-JEDI_COVARIANCE_COVARIANCE_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace oops {
  class Variables;
}
namespace wrf_hydro_nwm-jedi {
  class Geometry;
  class Increment;
  class State;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm-jedi {

  // Fields class
  class Covariance : public util::Printable,
                     private util::ObjectCounter<Covariance> {
   public:
    static const std::string classname() {return "wrf_hydro_nwm-jedi::Covariance";}

    Covariance(const Geometry &, const oops::Variables &,
                    const eckit::Configuration &,
                    const State &, const State &);
    ~Covariance();

    void multiply(const Increment &, Increment &) const;
    void inverseMultiply(const Increment &, Increment &) const;
    void randomize(Increment &) const;

   private:
    void print(std::ostream &) const;
  };

}  // namespace wrf_hydro_nwm-jedi

#endif  // WRF_HYDRO_NWM-JEDI_COVARIANCE_COVARIANCE_H_
