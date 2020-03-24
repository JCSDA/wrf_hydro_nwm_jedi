#!/usr/bin/env python

# ./geometry_preprocess.py --help
# ./geometry_preprocess.py \
#     --domain_dir /Users/jamesmcc/Downloads/croton_NY \
#     --config nwm_ana

# ./geometry_preprocess.py \
#     --domain_dir /Users/jamesmcc/Downloads/croton_NY \
#     --config nwm_ana \
#     --hrldas_patch_filename hrldas_namelist_patches.json \
#     --hydro_patch_filename = 'hrldas_namelist_patches.json \


import pathlib
import sys
import wrfhydropy
import xarray as xr


# This relies on wrfhydropy, make arguments to simply specify the individual files by hand.
def geometry_preprocess(
    domain_dir: pathlib.Path,
    config: str,
    hrldas_patch_filename: str = 'hrldas_namelist_patches.json',
    hydro_patch_filename: str = 'hydro_namelist_patches.json'):

    hrldas_patch_file = domain_dir / hrldas_patch_filename
    hydro_patch_file = domain_dir / hydro_patch_filename

    hrldas_patches = wrfhydropy.namelist.JSONNamelist(hrldas_patch_file)
    hydro_patches = wrfhydropy.namelist.JSONNamelist(hydro_patch_file)

    # HRLDAS
    hrldas_namelist = hrldas_patches.get_config(config)
    wrfinput_file = pathlib.Path(
        hrldas_namelist['noahlsm_offline']['hrldas_setup_file'])
    if not wrfinput_file.exists():
        wrfinput_file = domain_dir / wrfinput_file
    if not wrfinput_file.exists():
        raise FileExistsError('wrfinput file not found: ' + wrfinput_file)

    wrfinput_ds = xr.open_dataset(wrfinput_file)

    hrldas_drop_list = [
        'CANWAT', 'DZS', 'HGT', 'ISLTYP', 'IVGTYP',
        'LAI', 'MAPFAC_MX', 'MAPFAC_MY', 'SEAICE', 'SHDMAX',
        'SHDMIN', 'SMOIS', 'SNOW', 'TMN', 'TSK',
         'TSLB', 'XLAND']
    hrldas_geom = wrfinput_ds.drop_vars(hrldas_drop_list).squeeze('Time')

    hrldas_keep_attrs = ['DX', 'DY']
    hrldas_attrs = {
        key: val for key,val in wrfinput_ds.attrs.items() if key in hrldas_keep_attrs }

    hrldas_attrs['lsm_lat_dim_name'] = 'south_north'
    hrldas_attrs['lsm_lon_dim_name'] = 'west_east'
    hrldas_attrs['lsm_z_dim_name'] = 'soil_layers_stag'

    hrldas_attrs['lsm_lat_name'] = 'XLAT'
    hrldas_attrs['lsm_lon_name'] = 'XLON'
    hrldas_attrs['lsm_z_name'] = 'ZS'

    hrldas_geom.attrs = hrldas_attrs

    wrfinput_dir = wrfinput_file.parent
    out_file = wrfinput_dir / ('geometry_' + config + '.nc')
    print("Writing geometry file for config '" + config + "' to:\n" + str(out_file))
    hrldas_geom.to_netcdf(out_file)
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
