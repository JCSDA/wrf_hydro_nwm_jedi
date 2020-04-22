/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_STATE_STATE_H_
#define WRF_HYDRO_NWM_JEDI_STATE_STATE_H_

#include <memory>
#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/base/Variables.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
//#include "fv3jedi/Increment/Increment.h"
//#include "wrf_hydro_nwm_jedi/State/State.interface.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace oops {
  class Variables;
}
namespace ufo {
  class GeoVaLs;
  class Locations;
}
namespace wrf_hydro_nwm_jedi {
  class Fields;
  class GetValuesTraj;
  class Geometry;
  class Increment;
  typedef int F90state;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

  // State class
  class State : public util::Printable,
                private util::ObjectCounter<State> {
   public:
    static const std::string classname() {return "wrf_hydro_nwm-jedi::State";}

    // constructor, destructor
    State(const Geometry &, const oops::Variables &,
          const eckit::Configuration &);
    State(const Geometry &, const State &);
    State(const State &);
    State & operator=(const State &);
    ~State();

    // interpolate state to observations locations
    void getValues(const ufo::Locations &,
                   const oops::Variables &,
                   ufo::GeoVaLs &) const;
    void getValues(const ufo::Locations &,
                   const oops::Variables &,
                   ufo::GeoVaLs &,
                   GetValuesTraj &) const;

    // interactions with increment
    State & operator+=(const Increment &);

    const util::DateTime & validTime() const;
    util::DateTime & validTime();
    double norm() const;
    void write(const eckit::Configuration &) const;
    void zero();
    void accumul(const double &, const State &);
    void read(const eckit::Configuration &);

    boost::shared_ptr<const Geometry> geometry() const;// {return fields_->geometry();}
    
    /* F90state & toFortran() {return keyState_;} */
    const F90state & toFortran() const {return keyState_;}
   private:
    void print(std::ostream &) const;
    F90state keyState_;
    boost::scoped_ptr<Fields> fields_;
    boost::shared_ptr<const Geometry> geom_;
    oops::Variables vars_;
    util::DateTime time_;
  };
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_STATE_STATE_H_
