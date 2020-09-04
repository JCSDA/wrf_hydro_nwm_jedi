/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm_jedi/Fields/Fields.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/State/State.h"
#include "wrf_hydro_nwm_jedi/State/StateFortran.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace wrf_hydro_nwm_jedi {


  State::State(const Geometry & geom,
               const eckit::Configuration & conf)
    : vars_(conf, "state variables"),
      fields_(new Fields(geom, vars_)),
      time_(util::DateTime()) {
    oops::Log::trace() << "State::State 1 create from file." << std::endl;
    oops::Variables vars(vars_);
    wrf_hydro_nwm_jedi_state_create_f90(
       keyState_,
       fields_->geometry()->toFortran(),
       vars);
    util::DateTime * dtp = &time_;
    this->read_state_from_file(conf);
    oops::Log::trace() << "State::State 1 create from file done." << std::endl;
  }


  State::State(const Geometry & geom,
	       const oops::Variables & vars,
               const util::DateTime & time)
    : vars_(vars),
      fields_(new Fields(geom, vars_)),
      time_(time) {
    oops::Log::trace() << "State::State 2 create from file." << std::endl;
    oops::Log::trace() << "State::State 2 time_: " << time_ << std::endl;
    wrf_hydro_nwm_jedi_state_create_f90(
	keyState_,
	fields_->geometry()->toFortran(),
	vars);
    oops::Log::trace() << "State::State 2 create from file DONE." << std::endl;
  }


  State::State(const Geometry & geom,
               const State & other)
    : vars_ (other.vars_),
      fields_(new Fields(geom, other.vars_)),
      time_(other.time_) {
    wrf_hydro_nwm_jedi_state_create_f90(keyState_, geom.toFortran(), other.vars_);
    oops::Log::trace() << "State::State created from existing state." << std::endl;
  }


  State::State(const State & other)
    : fields_(new Fields(*other.fields_)) {
    std::cout << "State::State 3 create_from_other from State " << std::endl;
    wrf_hydro_nwm_jedi_state_create_from_other_f90(keyState_, other.keyState_);
    wrf_hydro_nwm_jedi_state_copy_f90(keyState_, other.keyState_);
    time_ = other.time_;
  }

  State::~State() { wrf_hydro_nwm_jedi_state_delete_f90(keyState_); }


  void State::getValues(const ufo::Locations & locs,
                        const oops::Variables & vars,
                        ufo::GeoVaLs & geovals) const {
    // util::abor1_cpp("State::getValues() needs to be implemented.",
    //                 __FILE__, __LINE__);
  }


  void State::getValues(const ufo::Locations & locs,
                 const oops::Variables & vars,
                 ufo::GeoVaLs & geovals,
                 GetValuesTraj & traj) const {
    // util::abor1_cpp("State::getValues() needs to be implemented.",
    //                 __FILE__, __LINE__);
  }


  void State::print(std::ostream & os) const {
    char *string = new char[8192];
    wrf_hydro_nwm_jedi_state_print_f90(keyState_, string);
    os << string;
    delete[] string;
    // os << *fields_;
    // int const nf = 1;
    // float pstat[3][nf];
    // wrf_hydro_nwm_jedi_state_get_mean_stddev_f90(keyState_,nf,pstat);
    // os << std::endl;
    // os << "Mean SNEQV: " << pstat[0][0] << std::endl;
    // os << "Std.dev SNEQV: " << pstat[1][0] << std::endl;
    // os << "RMS SNEQV: " << pstat[2][0] << std::endl;
  }


  // void State::write(const eckit::Configuration & config) const {
  //   const util::DateTime * dtp = &time_;
  //   const eckit::Configuration * conf = &config;
  //   wrf_hydro_nwm_jedi_state_write_file_f90(fields_->geometry()->toFortran(), keyState_, &conf, &dtp);
  void State::write(const eckit::Configuration & conf) const {
    util::abor1_cpp("State::write() needs to be implemented.",
                    __FILE__, __LINE__);
  }
  

  void State::read_state_from_file(const eckit::Configuration & config) {
    oops::Log::trace() << "State read starting" << std::endl;
    const eckit::Configuration * conf = &config;
    util::DateTime * dtp = &time_;
    // oops::Log::trace() << "before time_: " << time_ << std::endl;
    wrf_hydro_nwm_jedi_state_read_file_f90(
        fields_->geometry()->toFortran(),
        keyState_,
	&conf,
        &dtp);
    this->print(std::cout);
    // oops::Log::trace() << "after dtp: " << *dtp << std::endl;
    time_ = *dtp;
    oops::Log::trace() << "after time_: " << time_ << std::endl;
    oops::Log::trace() << "State read done" << std::endl;
  }


  State & State::operator=(const State & rhs) {
    wrf_hydro_nwm_jedi_state_copy_f90(keyState_, rhs.keyState_);
    time_ = rhs.time_;
    return *this;
  }


  State & State::operator+=(const Increment & dx)
  {
    // oops::Log::trace() << "State add increment starting" << std::endl;
    // ASSERT(this->validTime() == dx.validTime());
    // wrf_hydro_nwm_jedi_state_add_incr_f90(geom_->toFortran(), keyState_, dx.toFortran());
    // oops::Log::trace() << "State add increment done" << std::endl;
    util::abor1_cpp("State::operator+=(Increment) needs to be implemented.",
                    __FILE__, __LINE__);
    return *this;
  }


  void State::zero() {
    util::abor1_cpp("State::zero() needs to be implemented.",
                    __FILE__, __LINE__);
  }


  void State::accumul(const double &, const State &) {
    util::abor1_cpp("State::accumul() needs to be implemented.",
                    __FILE__, __LINE__);
  }


  // in State.h
  // const util::DateTime & State::validTime() const { return time_; }
  // util::DateTime & State::validTime() { return time_; }


  double State::norm() const {
    double norm = 0.0;
    norm = wrf_hydro_nwm_jedi_state_rms_f90(toFortran());
    return norm;
  }


  boost::shared_ptr<const Geometry> State::geometry() const {
    return fields_->geometry();
  }


}  // namespace wrf_hydro_nwm_jedi
