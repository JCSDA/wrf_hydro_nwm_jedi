/*
  * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "wrf_hydro_nwm_jedi/Traits.h"
#include "oops/runs/Run.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateVariableChangeFactory.h"
#include "test/interface/LinearModel.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateCovarFactory<wrf_hydro_nwm_jedi::Traits>();
  saber::instantiateVariableChangeFactory<wrf_hydro_nwm_jedi::Traits>();
  test::LinearModel<wrf_hydro_nwm_jedi::Traits> tests;
  return run.execute(tests);
}

