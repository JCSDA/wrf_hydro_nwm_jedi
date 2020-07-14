/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_COVARIANCE_H_
#define WRF_HYDRO_NWM_JEDI_COVARIANCE_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace oops {
  class Variables;
}
namespace wrf_hydro_nwm_jedi {
  class Geometry;
  class Increment;
  class State;
  typedef int F90bmat;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

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
    F90bmat keyFtnConfig_;
    boost::shared_ptr<const Geometry> geom_;
    util::DateTime time_;
  };

}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM-JEDI_COVARIANCE_H_
