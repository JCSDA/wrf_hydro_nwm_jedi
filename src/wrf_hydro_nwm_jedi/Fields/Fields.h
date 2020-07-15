/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_FIELDS_FIELDS_H_
#define WRF_HYDRO_NWM_JEDI_FIELDS_FIELDS_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/base/Variables.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace oops {
  class Variables;
}
namespace wrf_hydro_nwm_jedi {
  class Geometry;
  typedef int F90flds;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

  // Fields class
  class Fields : public util::Printable,
                 private util::ObjectCounter<Fields> {
   public:
    static const std::string classname() {return "wrf_hydro_nwm-jedi::Fields";}

    Fields(const Geometry &, const oops::Variables &,
           const eckit::Configuration &);

    Fields(const Geometry & geom, const oops::Variables & vars,
	   const util::DateTime & vt);
    
    ~Fields();

    const util::DateTime & time() const { return time_; }
    util::DateTime & time() { return time_; }

    double norm() const;

    boost::shared_ptr<const Geometry> geometry() const;

    // To be used to access Fields from Fortran, currently not needed
    F90flds & toFortran() {return keyFlds_;}
    const F90flds & toFortran() const {return keyFlds_;}

   private:
    void print(std::ostream &) const;
    boost::shared_ptr<const Geometry> geom_;
    oops::Variables vars_;
    F90flds keyFlds_;
    util::DateTime time_;
  };

}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_FIELDS_FIELDS_H_
