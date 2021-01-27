/*
 * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Run.h"
#include "oops/runs/StaticBInit.h"
#include "wrf_hydro_nwm_jedi/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::StaticBInit<wrf_hydro_nwm_jedi::Traits> bmat;
  return run.execute(bmat);
}
