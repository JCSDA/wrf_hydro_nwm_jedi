/*
 * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm_jedi/Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/Increment.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::Increment<wrf_hydro_nwm_jedi::Traits> tests;
  return run.execute(tests);
}

