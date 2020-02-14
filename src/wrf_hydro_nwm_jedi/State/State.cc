/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm-jedi/Fields/Fields.h"
#include "wrf_hydro_nwm-jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm-jedi/State/State.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace wrf_hydro_nwm-jedi {

// ----------------------------------------------------------------------------

  State::State(const Geometry & geom, const oops::Variables & vars,
               const eckit::Configuration & conf)
    : fields_(new Fields(geom, vars, conf)) {
    util::abor1_cpp("State::State() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  State::State(const Geometry &, const State &) {
    util::abor1_cpp("State::State() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  State::State(const State & other)
    : fields_(new Fields(*other.fields_)) {
    util::abor1_cpp("State::State() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  State::~State() {
    util::abor1_cpp("State::~State() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void State::getValues(const ufo::Locations & locs,
                        const oops::Variables & vars,
                        ufo::GeoVaLs & geovals) const {
    util::abor1_cpp("State::getValues() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void State::getValues(const ufo::Locations & locs,
                 const oops::Variables & vars,
                 ufo::GeoVaLs & geovals,
                 GetValuesTraj & traj) const {
    util::abor1_cpp("State::getValues() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void State::print(std::ostream & os) const {
    os << *fields_;
  }

// ----------------------------------------------------------------------------

  void State::write(const eckit::Configuration & conf) const {
    util::abor1_cpp("State::write() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  State & State::operator=(const State & rhs) {
    util::abor1_cpp("State::operator= needs to be implemented.",
                    __FILE__, __LINE__);
    return *this;
  }

// ----------------------------------------------------------------------------

  State & State::operator+=(const Increment & dx)
  {
    util::abor1_cpp("State::operator+=(Increment) needs to be implemented.",
                    __FILE__, __LINE__);
    return *this;
  }

// ----------------------------------------------------------------------------

  void State::zero() {
    util::abor1_cpp("State::zero() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void State::accumul(const double &, const State &) {
    util::abor1_cpp("State::accumul() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  const util::DateTime & State::validTime() const {
    return fields_->time();
  }

// ----------------------------------------------------------------------------

  util::DateTime & State::validTime() {
    return fields_->time();
  }

// ----------------------------------------------------------------------------

  double State::norm() const {
    return fields_->norm();
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm-jedi
