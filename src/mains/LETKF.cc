/*
 * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "wrf_hydro_nwm-jedi/Traits.h"
#include "oops/runs/LETKF.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::LETKF<wrf_hydro_nwm-jedi::Traits> letkf;
  run.execute(letkf);
  return 0;
}
