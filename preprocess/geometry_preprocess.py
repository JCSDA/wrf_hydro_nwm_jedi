#!/usr/bin/env python

# ./geometry_preprocess.py --help
# ./geometry_preprocess.py \
#     --domain_dir /Users/james/jedi/data/domains/v2.1/taylor_park_v2.1 \
#     --config nwm_long_range_snow

# ipython --pdb -c "%run geometry_preprocess.py     --domain_dir /Users/james/jedi/data/domains/v2.1/taylor_park_v2.1     --config nwm_long_range_snow"

# ./geometry_preprocess.py \
#     --domain_dir /Users/jamesmcc/Downloads/croton_NY \
#     --config nwm_ana \
#     --hrldas_patch_filename hrldas_namelist_patches.json \
#     --hydro_patch_filename = 'hrldas_namelist_patches.json \


import pathlib
import sys
import warnings
import wrfhydropy
import xarray as xr


# This relies on wrfhydropy, make arguments to simply specify the individual files by hand.
def geometry_preprocess(
        domain_dir: pathlib.Path,
        config: str,
        hrldas_patch_filename: str = 'hrldas_namelist_patches.json',
        hydro_patch_filename: str = 'hydro_namelist_patches.json',
        out_dir: pathlib.Path = None):

    hrldas_patch_file = domain_dir / hrldas_patch_filename
    hydro_patch_file = domain_dir / hydro_patch_filename
    hrldas_patches = wrfhydropy.namelist.JSONNamelist(hrldas_patch_file)
    hydro_patches = wrfhydropy.namelist.JSONNamelist(hydro_patch_file)

    # Low-resolution grid / HRLDAS / wrfinput file
    # Consider this part of the geometry mandatory, even if
    # running a channel-only simulation.
    hrldas_namelist = hrldas_patches.get_config(config)
    wrfinput_file = pathlib.Path(
        hrldas_namelist['noahlsm_offline']['hrldas_setup_file'])
    if not wrfinput_file.exists():
        wrfinput_file = domain_dir / wrfinput_file
    if not wrfinput_file.exists():
        raise FileExistsError('wrfinput file not found: ' + wrfinput_file)
    wrfinput_dir = wrfinput_file.parent
    wrfinput_ds = xr.open_dataset(wrfinput_file)
    # variables
    wrfinput_vars = [vv for vv in wrfinput_ds.variables]
    wrfinput_keep_vars = ['XLAT', 'XLONG', 'ZS']
    wrfinput_drop_list = list(
        set(wrfinput_vars).difference(wrfinput_keep_vars))
    lsm_geom = (
        wrfinput_ds.drop_vars(wrfinput_drop_list).squeeze('Time'))
    # attributes
    lsm_keep_attrs = ['DX', 'DY']
    lsm_attrs = {
        key: val for key, val in wrfinput_ds.attrs.items()
        if key in lsm_keep_attrs}
    # Input attrs
    lsm_attrs['domain_dir'] = str(domain_dir)
    lsm_attrs['config'] = config
    lsm_attrs['hrldas_patch_filename'] = str(hrldas_patch_filename)
    lsm_attrs['hydro_patch_filename'] = str(hydro_patch_filename)
    # Actual LSM attributes
    lsm_attrs['lsm_dx'] = lsm_attrs.pop('DX')
    lsm_attrs['lsm_dy'] = lsm_attrs.pop('DY')
    lsm_attrs['lsm_lat_dim_name'] = 'south_north'
    lsm_attrs['lsm_lon_dim_name'] = 'west_east'
    lsm_attrs['lsm_z_dim_name'] = 'soil_layers_stag'
    lsm_attrs['lsm_lat_name'] = 'XLAT'
    lsm_attrs['lsm_lon_name'] = 'XLON'
    lsm_attrs['lsm_z_name'] = 'ZS'
    lsm_attrs['lsm_src_file'] = str(wrfinput_file.absolute())
    lsm_attrs['lsm_src_md5'] = (
        wrfhydropy.core.ioutils.md5(wrfinput_file))
    lsm_geom.attrs = lsm_attrs
    # This is the start of the geometry that is output
    geom = lsm_geom
    

    # Streamflow "grid"  / Routelink file
    # This is not mandatory, throws a warning if missing
    hydro_namelist = hydro_patches.get_config(config)
    routelink_file = pathlib.Path(
        hydro_namelist['hydro_nlist']['route_link_f'])
    if not routelink_file.exists():
        routelink_file = domain_dir / routelink_file
    if not wrfinput_file.exists():
        warnings.warn('routelink file not found: ' + routelink_file)
    else:
        routelink_dir = routelink_file.parent
        routelink_ds = xr.open_dataset(routelink_file)
        # variables
        routelink_vars = [vv for vv in routelink_ds.variables]
        routelink_keep_vars = ['lat', 'lon']
        routelink_drop_vars = (
            set(routelink_vars).difference(set(routelink_keep_vars)))
        routelink_geom = (
            routelink_ds
            .drop_vars(list(routelink_drop_vars))
            .reset_coords(['lat', 'lon']))  # fixed.
        routelink_keep_attrs = []
        routelink_attrs = {
            key: val for key, val in routelink_ds.attrs.items()
            if key in routelink_keep_attrs}
        routelink_attrs['stream_dim_name'] = 'feature_id'
        routelink_attrs['stream_lat_name'] = 'lat'
        routelink_attrs['stream_lon_name'] = 'lon'
        routelink_attrs['stream_src_file'] = str(routelink_file.absolute())
        routelink_attrs['stream_src_md5'] = (
            wrfhydropy.core.ioutils.md5(routelink_file))
        routelink_geom.attrs = routelink_attrs
        # Merge this geometry into the output geom
        _ = geom.attrs.update(routelink_geom.attrs)
        geom_atts = geom.attrs
        geom = xr.merge([geom, routelink_geom])
        geom.attrs = geom_atts

    # Addtional geoms to be added:
    # * High resolution grid: optional
    # * Lakes: optional
    # * Buckets: optional

    # Write out the geometry file.
    if out_dir is None:
        if ('routelink_dir' in locals() and
            routelink_dir == wrfinput_dir):
            out_dir = wrfinput_dir

    out_file = out_dir / ('geometry_' + config + '.nc')
    print("Writing geometry file for config '" +
          config + "' to:\n" + str(out_file))
    geom.to_netcdf(out_file)
    return(0)


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--domain_dir',
        metavar='domain dir',
        type=str,
        help='The domain dir containing the model geometry files'
    )
    parser.add_argument(
        '--config',
        metavar='configuration name to be taken from json files',
        type=str,
        help='The directory to write collected CHANOBS for the full ensemble'
    )
    parser.add_argument(
        '--hrldas_patch_filename',
        metavar='hrldas_patch_filename',
        type=str,
        required=False,
        default='hrldas_namelist_patches.json',
        help='the name of the hrldas patch json file in domain_dir'
    )
    parser.add_argument(
        '--hydro_patch_filename',
        metavar='hydro_patch_filename',
        type=str,
        required=False,
        default='hydro_namelist_patches.json',
        help='the name of the hydro patch json file in domain_dir'
    )
    args = parser.parse_args()
    domain_dir = pathlib.Path(args.domain_dir)
    return domain_dir, args.config, args.hrldas_patch_filename, args.hydro_patch_filename


if __name__ == "__main__":
    domain_dir, config, hrldas_patch_filename, hydro_patch_filename = parse_arguments()
    return_value = geometry_preprocess(
        domain_dir, config, hrldas_patch_filename, hydro_patch_filename)
    sys.exit(return_value)
