/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM-JEDI_GETVALUESTRAJ_GETVALUESTRAJ_H_
#define WRF_HYDRO_NWM-JEDI_GETVALUESTRAJ_GETVALUESTRAJ_H_

#include <ostream>

#include "oops/util/Printable.h"

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm-jedi {

  // GetValuesTraj class
  class GetValuesTraj : public util::Printable {
   public:
    GetValuesTraj();
    ~GetValuesTraj();

   private:
    void print(std::ostream & os) const;
  };

}  // namespace wrf_hydro_nwm-jedi

#endif  // WRF_HYDRO_NWM-JEDI_GETVALUESTRAJ_GETVALUESTRAJ_H_
