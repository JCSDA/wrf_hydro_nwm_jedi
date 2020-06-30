import numpy as np
import pathlib
import xarray as xr

# machine_path = pathlib.Path('/home/vagrant')
machine_path = pathlib.Path('/Users/james/')

repo_path = (machine_path /
             'jedi/bundle-nwm/wrf_hydro_nwm_jedi')

ioda_path = (machine_path /
             'jedi/bundle-ioda/ioda-converters')

out_path = repo_path / 'test/Data/ioda_obs/3d_var_test'

ioda_file = ioda_path / 'test/testoutput/nwm_osse_ioda.nc'
ioda_ds = xr.open_dataset(ioda_file)

# make the obs 5mm more than the modeled values
ioda_ds['swe@ObsValue'] = (
    ioda_ds['swe@ObsValue'] + 5.00)

# Set a variety of obs error variance levels
for pct_oev in [5, 10, 20, 40, 80]:
    ioda_ds['swe@ObsError'] = (
        ioda_ds['swe@ObsValue'] * pct_oev/100)
    out_file = out_path / (
        'nwm_osse_ioda_pct_oev_' +
        str(pct_oev).zfill(2) + '.nc')
    ioda_ds.to_netcdf(out_file)


## Provide a solution -------


def kalman_eqn(x_b, B, H, R, y):
    # Calculate the Kalman equation
    # Solve:
    # x_a =
    #   x_b + B*Ht * 1/(H*B*Ht+R) * (y-H*x_b)
    #   x_b + T0   * T1           * T2
    T0 = np.matmul(B, H.transpose())
    T1 = (
        np.linalg.inv(
            np.matmul(
                np.matmul(H, B),
                H.transpose())
            + R))
    T2 = y - np.matmul(H, x_b)
    x_a = (
        x_b +
        np.matmul(
            np.matmul(T0, T1), T2))
    return(x_a)


# This is the "modeled values" before adding the 5mm
# to the obs
x_b = (xr
       .open_dataset(ioda_file)['swe@ObsValue']
       .values)
B = np.diag(x_b * .2)
H = np.diag(np.ones(len(x_b)))

ans_dict = {}
for pct_oev in [5, 10, 20, 40, 80]:
    obs_file = out_path / (
        'nwm_osse_ioda_pct_oev_' +
        str(pct_oev).zfill(2) + '.nc')
    obs_ds = (xr.open_dataset(obs_file))
    # The obs dont actually change... but oh well.
    y = obs_ds['swe@ObsValue'].values
    R = np.diag(obs_ds['swe@ObsError'].values)
    x_a = kalman_eqn(x_b, B, H, R, y)
    ans_dict[str(pct_oev).zfill(2)] = x_a

# quick look. could put these into the files as "answers"
for key in ans_dict.keys():
    print(ans_dict[key][0:5])
