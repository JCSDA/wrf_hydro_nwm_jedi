Preprocessing
=============


## Geometry 

```
(whp) jedi@vagrant[514]:~/jedi/bundle-nwm/wrf_hydro_nwm_jedi/preprocess> ./geometry_preprocess.py     --domain_dir ~/jedi/domains/private/taylor_park_v2_1     --config nwm_long_range_snow     --hrldas_patch_filename hrldas_namelist_patches.json     --hydro_patch_filename hydro_namelist_patches.json 
Writing geometry file for config 'nwm_long_range_snow' to:
/home/vagrant/jedi/domains/private/taylor_park_v2_1/NWM/DOMAIN/geometry_nwm_long_range_snow.nc

(whp) jedi@vagrant[515]:~/jedi/bundle-nwm/wrf_hydro_nwm_jedi/preprocess> ncdump -h /home/vagrant/jedi/domains/private/taylor_park_v2_1/NWM/DOMAIN/geometry_nwm_long_range_snow.nc
netcdf geometry_nwm_long_range_snow {
dimensions:
	south_north = 24 ;
	west_east = 30 ;
	soil_layers_stag = 4 ;
	feature_id = 189 ;
variables:
	float HGT(south_north, west_east) ;
		HGT:_FillValue = NaNf ;
		HGT:FieldType = 104 ;
		HGT:MemoryOrder = "XY " ;
		HGT:units = "meters MSL" ;
		HGT:description = "Topography height" ;
		HGT:stagger = "M" ;
		HGT:sr_x = 1 ;
		HGT:sr_y = 1 ;
	float XLAT(south_north, west_east) ;
		XLAT:_FillValue = NaNf ;
		XLAT:FieldType = 104 ;
		XLAT:MemoryOrder = "XY " ;
		XLAT:units = "degrees latitude" ;
		XLAT:description = "Latitude on mass grid" ;
		XLAT:stagger = "M" ;
		XLAT:sr_x = 1 ;
		XLAT:sr_y = 1 ;
	float XLONG(south_north, west_east) ;
		XLONG:_FillValue = NaNf ;
		XLONG:FieldType = 104 ;
		XLONG:MemoryOrder = "XY " ;
		XLONG:units = "degrees longitude" ;
		XLONG:description = "Longitude on mass grid" ;
		XLONG:stagger = "M" ;
		XLONG:sr_x = 1 ;
		XLONG:sr_y = 1 ;
	float ZS(soil_layers_stag) ;
		ZS:_FillValue = -1.e+36f ;
		ZS:units = "m" ;
	float Length(feature_id) ;
		Length:_FillValue = NaNf ;
		Length:long_name = "Stream length (m)" ;
	float lat(feature_id) ;
		lat:_FillValue = NaNf ;
		lat:long_name = "latitude of the segment midpoint" ;
		lat:units = "degrees_north" ;
		lat:standard_name = "latitude" ;
	float lon(feature_id) ;
		lon:_FillValue = NaNf ;
		lon:long_name = "longitude of the segment midpoint" ;
		lon:units = "degrees_east" ;
		lon:standard_name = "longitude" ;

// global attributes:
		:domain_dir = "/home/vagrant/jedi/domains/private/taylor_park_v2_1" ;
		:config = "nwm_long_range_snow" ;
		:hrldas_patch_filename = "hrldas_namelist_patches.json" ;
		:hydro_patch_filename = "hydro_namelist_patches.json" ;
		:lsm_dx = 1000.f ;
		:lsm_dy = 1000.f ;
		:lsm_xdim_name = "west_east" ;
		:lsm_ydim_name = "south_north" ;
		:lsm_zdim_name = "soil_layers_stag" ;
		:lsm_lat_name = "XLAT" ;
		:lsm_lon_name = "XLONG" ;
		:lsm_z_name = "ZS" ;
		:lsm_sfc_elev_name = "HGT" ;
		:lsm_src_file = "/home/vagrant/jedi/domains/private/taylor_park_v2_1/NWM/DOMAIN/wrfinput.nc" ;
		:lsm_src_md5 = "849001936d4bd1f45f488e6618d45658" ;
		:stream_dx_name = "Length" ;
		:stream_xdim_name = "feature_id" ;
		:stream_lat_name = "lat" ;
		:stream_lon_name = "lon" ;
		:stream_src_file = "/home/vagrant/jedi/domains/private/taylor_park_v2_1/NWM/DOMAIN/RouteLink.nc" ;
		:stream_src_md5 = "b43e1fa8443092671bf1338556152f04" ;
}
```
