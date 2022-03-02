/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
// #include "wrf_hydro_nwm_jedi/GeometryIterator/GeometryIterator.h"
#include "wrf_hydro_nwm_jedi/Geometry/GeometryFortran.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

namespace wrf_hydro_nwm_jedi {

  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm)
    : comm_(comm) {
    const eckit::Configuration * configc = &conf;
    wrf_hydro_nwm_jedi_geometry_setup_f90(keyGeom_, &configc);

    // Set ATLAS lon/lat field
    atlasFieldSet_.reset(new atlas::FieldSet());
    wrf_hydro_nwm_jedi_geometry_set_atlas_lonlat_f90(keyGeom_, atlasFieldSet_->get());
    atlas::Field atlasField = atlasFieldSet_->field("lonlat");

    // Create ATLAS function space
    atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(atlasField));

    // Set ATLAS function space pointer in Fortran
    wrf_hydro_nwm_jedi_geometry_set_atlas_functionspace_pointer_f90(keyGeom_,
        atlasFunctionSpace_->get());

    // Fill ATLAS fieldset
    atlasFieldSet_.reset(new atlas::FieldSet());
    wrf_hydro_nwm_jedi_geometry_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());
  }


  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_) {
    wrf_hydro_nwm_jedi_geometry_clone_f90(keyGeom_, other.keyGeom_);

    // Copy ATLAS function space
    atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(
                              other.atlasFunctionSpace_->lonlat()));

    // Set ATLAS function space pointer in Fortran
    wrf_hydro_nwm_jedi_geometry_set_atlas_functionspace_pointer_f90(keyGeom_,
        atlasFunctionSpace_.get()->get());

    // Copy ATLAS fieldset
    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
      atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
      atlasFieldSet_->add(atlasField);
    }
  }

  Geometry::~Geometry() {
    wrf_hydro_nwm_jedi_geometry_delete_f90(keyGeom_);
  }


  GeometryIterator Geometry::begin() const {
//    util::abor1_cpp("Geometry::begin() needs to be implemented.",
//                    __FILE__, __LINE__);
    return GeometryIterator(*this, 1, 1);
  }


  GeometryIterator Geometry::end() const {
    float dx, dy;
    int npx, npy, npz;

    wrf_hydro_nwm_jedi_geometry_info_f90(keyGeom_, &dx, &dy, &npx, &npy, &npz);
    return GeometryIterator(*this, -1, -1);
//    return GeometryIterator(*this, npx, npy);
  }

// -----------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables &
                                                    vars) const {
    std::vector<size_t> varSizes(vars.size());
    wrf_hydro_nwm_jedi_geoval_levels_f90(keyGeom_, vars, vars.size(),  varSizes[0]);

  return varSizes;
}
// -----------------------------------------------------------------------------
void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
                      const bool halo) const {
  const atlas::functionspace::PointCloud * fspace;

  fspace = atlasFunctionSpace_.get();

  const auto lonlat = atlas::array::make_view<double, 2>(fspace->lonlat());
  const size_t npts = fspace->size();
  lats.resize(npts);
  lons.resize(npts);
  for (size_t jj = 0; jj < npts; ++jj) {
    lats[jj] = lonlat(jj, 1);
    lons[jj] = lonlat(jj, 0);
    if (lons[jj] < 0.0) lons[jj] += 360.0;
  }
}
// -----------------------------------------------------------------------------
  void Geometry::print(std::ostream & os) const {
    float dx, dy;
    int npx, npy, npz;
    // util::abor1_cpp("Geometry::print() needs to be implemented.",
    //                 __FILE__, __LINE__);
    os << "Geometry: "
       << "(TODO, print diagnostic info about the geometry here)"
       << std::endl;
    wrf_hydro_nwm_jedi_geometry_info_f90(keyGeom_, &dx, &dy, &npx, &npy, &npz);
    os << "dx = " << dx << ", dy = " << dy << std::endl;
  }


}  // namespace wrf_hydro_nwm_jedi
