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
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/mpi/Comm.h"

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/base/WriteParametersBase.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "wrf_hydro_nwm_jedi/Fields/Fields.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace oops {
  class Variables;
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
    State(const Geometry &,
          const oops::Variables &,
          const util::DateTime &);
    State(const Geometry &,
          const eckit::Configuration &);
    State(const Geometry &, const State &);
    State(const State &);
    State & operator=(const State &);
    ~State();

    // interactions with increment
    State & operator+=(const Increment &);

    double norm() const;
    void zero();
    void zero(const util::DateTime &);
    void ones();
    void accumul(const double &, const State &);

    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &);
    // void write(const eckit::Configuration &) const;

    std::shared_ptr<const Geometry> geometry() const {return fields_->geometry();}

    const oops::Variables & variables() const {return vars_;}

    // time()
    const util::DateTime & time() const {return time_;}
    util::DateTime & time() {return time_;}

    // Needed by PseudoModel
    void updateTime(const util::Duration & dt) {time_ += dt;}

    /// Serialize and deserialize
    size_t serialSize() const {return 0;}
    void serialize(std::vector<double> &) const {}
    void deserialize(const std::vector<double> &, size_t &) {}

    // validTime()
    const util::DateTime & validTime() const { return time_; }
    util::DateTime & validTime() { return time_; }

    // Accessors to the ATLAS fieldset
    void toFieldSet(atlas::FieldSet &) const;
    void fromFieldSet(const atlas::FieldSet &);

    /* F90state & toFortran() {return keyState_;} */
    const F90state & toFortran() const {return keyState_;}

   private:
     void print(std::ostream &) const;
     F90state keyState_;
     oops::Variables vars_;
     boost::scoped_ptr<Fields> fields_;
     util::DateTime time_;
  };

}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_STATE_STATE_H_
