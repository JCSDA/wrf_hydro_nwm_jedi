/*
 * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm_jedi/Traits.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsFilterFactory<wrf_hydro_nwm_jedi::Traits>();
  oops::Variational<wrf_hydro_nwm_jedi::Traits> var;
  run.execute(var);
  return 0;
}
