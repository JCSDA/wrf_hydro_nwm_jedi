/*
 * (C) Copyright 2019-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Run.h"
#include "saber/oops/ErrorCovarianceTraining.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "wrf_hydro_nwm_jedi/Traits.h"


int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateCovarFactory<wrf_hydro_nwm_jedi::Traits>();
  saber::ErrorCovarianceTraining<wrf_hydro_nwm_jedi::Traits> dir;
  return run.execute(dir);
}
