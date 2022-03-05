/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field.h"

#include "wrf_hydro_nwm_jedi/Fields/Fields.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"
#include "wrf_hydro_nwm_jedi/State/State.h"
#include "wrf_hydro_nwm_jedi/State/StateFortran.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

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
    this->read(conf);
    oops::Log::trace() << "State::State 1 create from file END." << std::endl;
  }


  State::State(const Geometry & geom,
               const oops::Variables & vars,
               const util::DateTime & time)
    : vars_(vars),
      fields_(new Fields(geom, vars_)),
      time_(time) {
    oops::Log::trace() << "State::State 2 create from file." << std::endl;
    wrf_hydro_nwm_jedi_state_create_f90(
        keyState_,
        fields_->geometry()->toFortran(),
        vars);
    oops::Log::trace() << "State::State 2 create from file END." << std::endl;
  }

  State::State(const Geometry & geom,
               const State & other)
    : vars_(other.vars_),
      fields_(new Fields(geom, other.vars_)),
      time_(other.time_) {
    oops::Log::trace() <<
      "State::State 3 create from existing geom and state." << std::endl;
    wrf_hydro_nwm_jedi_state_create_f90(keyState_, geom.toFortran(), other.vars_);
    wrf_hydro_nwm_jedi_state_create_from_other_f90(keyState_, other.keyState_);
    wrf_hydro_nwm_jedi_state_copy_f90(keyState_, other.keyState_);
    time_ = other.time_;
    oops::Log::trace() <<
      "State::State 3 create from existing geom and state END." << std::endl;
  }

  State::State(const State & other)
    : vars_(other.vars_),
      fields_(new Fields(*other.fields_)) {
    oops::Log::trace() <<
      "State::State 4 create_from_other from State " << std::endl;
    wrf_hydro_nwm_jedi_state_create_from_other_f90(keyState_, other.keyState_);
    wrf_hydro_nwm_jedi_state_copy_f90(keyState_, other.keyState_);
    time_ = other.time_;
    oops::Log::trace() <<
      "State::State 4 create_from_other from State End" << std::endl;
  }

  State::~State() { wrf_hydro_nwm_jedi_state_delete_f90(keyState_); }

  void State::print(std::ostream & os) const {
    char *string = new char[8192];
    wrf_hydro_nwm_jedi_state_print_f90(keyState_, string);
    os << std::endl << "Print State (C++) ------------------------ ";
    os << string;
    os << "End Print State (C++) -------------------- " << std::endl  << std::endl;
    delete[] string;
  }

  void State::write(const eckit::Configuration & config) {  // const {
    oops::Log::trace() << "State::State write start" << std::endl;
    const eckit::Configuration * conf = &config;
    util::DateTime * dtp = &time_;
    wrf_hydro_nwm_jedi_state_write_file_f90(
        fields_->geometry()->toFortran(),
        keyState_,
        &conf,
        &dtp);
    // this->print(std::cout);
    time_ = *dtp;
    oops::Log::trace() << "State write done" << std::endl;
  }


  void State::read(const eckit::Configuration & config) {
    oops::Log::trace() << "State::State read start" << std::endl;
    const eckit::Configuration * conf = &config;
    util::DateTime * dtp = &time_;
    wrf_hydro_nwm_jedi_state_read_file_f90(
        fields_->geometry()->toFortran(),
        keyState_,
        &conf,
        &dtp);
    // this->print(std::cout);
    time_ = *dtp;
    oops::Log::trace() << "State read done" << std::endl;
  }


  State & State::operator=(const State & rhs) {
    wrf_hydro_nwm_jedi_state_copy_f90(keyState_, rhs.keyState_);
    time_ = rhs.time_;
    return *this;
  }


  State & State::operator+=(const Increment & dx)
  {
    oops::Log::trace() << "State add increment starting" << std::endl;
    ASSERT(this->validTime() == dx.validTime());
    wrf_hydro_nwm_jedi_state_add_incr_f90(keyState_, dx.toFortran());
    oops::Log::trace() << "State add increment done" << std::endl;
    return *this;
  }


  void State::zero() {
    wrf_hydro_nwm_jedi_state_zero_f90(keyState_);
  }

  void State::zero(const util::DateTime & time) {
    wrf_hydro_nwm_jedi_state_zero_f90(keyState_);
    time_ = time;
  }

  void State::ones() {
    wrf_hydro_nwm_jedi_state_ones_f90(keyState_);
  }

  void State::accumul(const double & zz, const State & xx) {
    oops::Log::trace() << "State accumul starting" << std::endl;
    wrf_hydro_nwm_jedi_state_axpy_f90(keyState_, zz, xx.keyState_);
    oops::Log::trace() << "State accumul done" << std::endl;
  }

  double State::norm() const {
    double norm = 0.0;
    norm = wrf_hydro_nwm_jedi_state_rms_f90(toFortran());
    return norm;
  }

  void State::getFieldSet(const oops::Variables & vars, atlas::FieldSet & fset) const {
    const bool include_halo = true;
    wrf_hydro_nwm_jedi_state_set_atlas_f90(keyState_, fields_->geometry()->toFortran(),
                                            vars, fset.get(), include_halo);
    wrf_hydro_nwm_jedi_state_to_atlas_f90(keyState_, fields_->geometry()->toFortran(),
                                            vars, fset.get(), include_halo);
  }

}  // namespace wrf_hydro_nwm_jedi

