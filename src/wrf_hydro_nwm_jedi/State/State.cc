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

// ----------------------------------------------------------------------------

  State::State(const Geometry & geom, const oops::Variables & vars,
               const eckit::Configuration & conf)
    : fields_(new Fields(geom, vars, conf)) {

    oops::Log::trace() << "State::State create from analytical or"
      " from file." << std::endl;

    oops::Variables lvars;
    if (conf.has("variables")) {
      oops::Variables lvars(conf);
      this->vars_ = lvars;
      std::cout << "From conf"<<std::endl;
    } else {
      this->vars_ = vars;
      std::cout << "From argument"<<std::endl;
    }

    wrf_hydro_nwm_jedi_state_create_f90(keyState_,
					fields_->geometry()->toFortran(), vars_);
    
    // Analytical or read from file
    // if (conf.has("analytic_init")) {
    //   this->analytic_init(conf, geom);
    // } else {
    this->read(conf);
      //}
    
    oops::Log::trace() << "State::State create from analytical or"
      " from file done." << std::endl;
    // util::abor1_cpp("State::State() needs to be implemented.",
    //                 __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  State::State(const Geometry & geom,
	       const oops::Variables & vars,
	       const util::DateTime & vt)
    : fields_(new Fields(geom, vars, vt)){

    // util::abor1_cpp("State::State(const Geometry &,const Variables &, const DateTime & ) needs to be implemented.", __FILE__, __LINE__);
    std::cout<< "Creating state as increment"<<std::endl;
    wrf_hydro_nwm_jedi_state_create_f90(keyState_, geom.toFortran(), vars);
    vars_ = vars;
    time_ = vt;
  }

  // ----------------------------------------------------------------------------

  State::State(const Geometry & resol, const State & other) :
    fields_(new Fields(resol, other.vars_, other.time_)){
    // util::abor1_cpp("State::State(const Geometry &, const State &) needs to be implemented.", __FILE__, __LINE__);
    wrf_hydro_nwm_jedi_state_create_f90(keyState_, resol.toFortran(), other.vars_);
    vars_ = other.vars_;
    time_ = other.time_;
    oops::Log::trace() << "State::State created by interpolation."
                     << std::endl;
  }

// ----------------------------------------------------------------------------

  State::State(const State & other)
    : fields_(new Fields(*other.fields_)) {

    // wrf_hydro_nwm_jedi_state_create_f90(keyState_, other.fields_->geometry()->toFortran(), other.vars_);

    wrf_hydro_nwm_jedi_state_create_from_other_f90(keyState_, other.keyState_);
    wrf_hydro_nwm_jedi_state_copy_f90(keyState_, other.keyState_);
    time_ = fields_->time();
  }

// ----------------------------------------------------------------------------

  State::~State() {
    
    wrf_hydro_nwm_jedi_state_delete_f90(keyState_);
    
  }

// ----------------------------------------------------------------------------

  void State::getValues(const ufo::Locations & locs,
                        const oops::Variables & vars,
                        ufo::GeoVaLs & geovals) const {
    // util::abor1_cpp("State::getValues() needs to be implemented.",
    //                 __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void State::getValues(const ufo::Locations & locs,
                 const oops::Variables & vars,
                 ufo::GeoVaLs & geovals,
                 GetValuesTraj & traj) const {
    // util::abor1_cpp("State::getValues() needs to be implemented.",
    //                 __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void State::print(std::ostream & os) const {
    char *string = new char[8192];
    wrf_hydro_nwm_jedi_state_print_f90(keyState_,string);
    os << string;
    delete[] string;
  }

// ----------------------------------------------------------------------------

  void State::write(const eckit::Configuration & config) const {
    const util::DateTime * dtp = &time_;
    const eckit::Configuration * conf = &config;
    wrf_hydro_nwm_jedi_state_write_file_f90(fields_->geometry()->toFortran(), keyState_, &conf, &dtp);
  }
  
// ----------------------------------------------------------------------------

  void State::read(const eckit::Configuration & config) {
    oops::Log::trace() << "State read starting" << std::endl;
    const eckit::Configuration * conf = &config;
    util::DateTime * dtp = &time_;
    wrf_hydro_nwm_jedi_state_read_file_f90(fields_->geometry()->toFortran(), keyState_, &conf, &dtp);
    oops::Log::trace() << "State read done" << std::endl;
  }

// ----------------------------------------------------------------------------

  State & State::operator=(const State & rhs) {
    wrf_hydro_nwm_jedi_state_copy_f90(keyState_, rhs.keyState_);
    time_ = rhs.time_;
    return *this;
  }

// ----------------------------------------------------------------------------

  State & State::operator+=(const Increment & dx)
  {
    oops::Log::trace() << "State add increment starting" << std::endl;
    ASSERT(this->validTime() == dx.validTime());
    wrf_hydro_nwm_jedi_state_add_incr_f90(geom_->toFortran(), keyState_, dx.toFortran());
    oops::Log::trace() << "State add increment done" << std::endl;
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
    std::cout << "Norm requested from State"<<std::endl;
    return fields_->norm();
  }

  boost::shared_ptr<const Geometry> State::geometry() const {
    return fields_->geometry();
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
