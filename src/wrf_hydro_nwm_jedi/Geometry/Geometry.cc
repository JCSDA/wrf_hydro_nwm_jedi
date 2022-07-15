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

    // Set lon/lat field
    atlas::FieldSet fs;
    const bool include_halo = true;  // include halo so we can set up both function spaces
    wrf_hydro_nwm_jedi_geometry_set_lonlat_f90(keyGeom_, fs.get(), include_halo);

    // Create function space
    const atlas::Field lonlatField = fs.field("lonlat");
    functionSpace_ = atlas::functionspace::PointCloud(lonlatField);

    ASSERT(include_halo);

    // Create function space with halo
    const atlas::Field lonlatFieldInclHalo = fs.field("lonlat_including_halo");
    functionSpaceIncludingHalo_ = atlas::functionspace::PointCloud(lonlatFieldInclHalo);

    // Set function space pointer in Fortran
    wrf_hydro_nwm_jedi_geometry_set_functionspace_pointer_f90(keyGeom_,
                                                   functionSpace_.get(),
                                                   functionSpaceIncludingHalo_.get());

    // Fill extra fields
    extraFields_ = atlas::FieldSet();
    wrf_hydro_nwm_jedi_geometry_fill_extra_fields_f90(keyGeom_, extraFields_.get());
  }

  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_) {
    wrf_hydro_nwm_jedi_geometry_clone_f90(keyGeom_, other.keyGeom_);

    // Copy ATLAS function space
    functionSpace_ = atlas::functionspace::PointCloud(other.functionSpace_.lonlat());
    functionSpaceIncludingHalo_ = atlas::functionspace::PointCloud(
                                         other.functionSpaceIncludingHalo_.lonlat());

    // Set ATLAS function space pointer in Fortran
    wrf_hydro_nwm_jedi_geometry_set_functionspace_pointer_f90(keyGeom_,
        functionSpace_.get(), functionSpaceIncludingHalo_.get());

    // Copy ATLAS fieldset
    extraFields_ = atlas::FieldSet();
    for (auto & field : other.extraFields_) {
      extraFields_->add(field);
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
const atlas::FunctionSpace * fspace;
fspace = &functionSpace_;

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
