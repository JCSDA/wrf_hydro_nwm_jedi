/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "wrf_hydro_nwm-jedi/GetValuesTraj/GetValuesTraj.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace wrf_hydro_nwm-jedi {

// ----------------------------------------------------------------------------

  GetValuesTraj::GetValuesTraj() {
    util::abor1_cpp("GetValuesTraj::GetValuesTraj() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  GetValuesTraj::~GetValuesTraj() {
    util::abor1_cpp("GetValuesTraj::~GetValuesTraj() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void GetValuesTraj::print(std::ostream & os) const {
    util::abor1_cpp("GetValuesTraj::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "(TODO, print diagnostic info about the GetValuesTraj here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm-jedi
