/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "wrf_hydro_nwm_jedi/LinearVariableChange/LinearVariableChange.h"

#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"
#include "wrf_hydro_nwm_jedi/State/State.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace wrf_hydro_nwm_jedi {

// -----------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom, const Parameters_ & params)
  : geom_(new Geometry(geom)), params_(params), linearVariableChange_() {}

// -------------------------------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::setTrajectory(const State & xbg, const State & xfg) {
  oops::Log::trace() << "LinearVariableChange::setTrajectory starting" << std::endl;
  // Create the variable change
  linearVariableChange_.reset(LinearVariableChangeFactory::create(xbg, xfg, *geom_,
             params_.linearVariableChangeParameters.value()));
  oops::Log::trace() << "LinearVariableChange::setTrajectory done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiply(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiply starting" << std::endl;

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiply(dx, dxout);

  // Copy data from temporary state
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiply done" << dx << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyInverse(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyInverse starting" << std::endl;

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyInverse(dx, dxout);

  // Copy data from temporary state
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiplyInverse done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyAD(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyAD starting" << std::endl;

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyAD(dx, dxout);

  // Copy data from temporary state
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyInverseAD(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyInverseAD starting" << std::endl;

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyInverseAD(dx, dxout);

  // Copy data from temporary state
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::multiplyInverseAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::print(std::ostream & os) const {
  os << "WRF-Hydro variable change";
}

// -----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
