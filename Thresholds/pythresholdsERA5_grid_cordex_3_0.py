##################################################################################################
#  This code is written to detect Gonu Tropical Cyclone according to :
#
#  Camargo, S. J., & Zebiak, S. E. (2002). Improving the Detection and Tracking of Tropical
#  Cyclones in Atmospheric General Circulation Models. Weather and   Forecasting, 17(6), 1152–1162.
#  https://doi.org/10.1175/1520-0434(2002)017<1152:ITDATO>2.0.CO;2
#
#  Written by :Ahmed Homoudi
#  Modified for Python: Ana Aguilar
#  Septiembre 2020
##################################################################################################

######################################### Load Libraries #########################################
print("Loading packages...")
from datetime import datetime
t1 = to = datetime.now()
import sys
import os
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
import xarray as xr
import dask
import pythreshold_support as pth
import fortran_modules.msl_calc as mslc
import fortran_modules.vorticity as vorc

def get_data(filepath):
    '''Reads and organize dataset from netCDF file.
    Returns an ndarray, dimensions and coordinates of the dataset.'''
    with xr.open_mfdataset(filepath, decode_cf=True, decode_times=True, combine='by_coords', parallel=True) as data:
        return data.transpose('longitude', 'latitude', 'time').to_array().values.squeeze(), \
               data.dims, data.coords, data.attrs


def haversine(s_lat, s_lng, e_lat, e_lng):
    '''Calculate distance using the Haversine formula.'''
    # approximate radius of earth in meters
    R = 6371000.0

    s_lat = s_lat*np.pi/180.0
    s_lng = np.deg2rad(s_lng)
    e_lat = np.deg2rad(e_lat)
    e_lng = np.deg2rad(e_lng)

    d = np.sin((e_lat - s_lat)/2)**2 + np.cos(s_lat)*np.cos(e_lat) * np.sin((e_lng - s_lng)/2)**2

    return 2 * R * np.arctan2(np.sqrt(d), np.sqrt(1 - d))

def calc_msl(t850_file, zg850_file):
    '''MSL calculation.'''
    ta850_array, _, _, _ = get_data(t850_file)

    zg850_array, dims, coords, attrs = get_data(zg850_file)

    nlon, nlat, ntime = ta850_array.shape

    msl_array = mslc.p850_to_msl(lon=nlon,
                                 lat=nlat,
                                 time=ntime,
                                 zg_array=zg850_array,
                                 ta_array=ta850_array)

    return msl_array.swapaxes(0, 2)


def calc_vort(ua850_file, va850_file):
    '''Vorticity calculation.'''
    ua850_array, dims, coords, attrs = get_data(ua850_file)

    with xr.open_dataset(ua850_file, decode_cf=True, decode_times=True) as data:
        lons = data.longitude.to_masked_array()
        lats = data.latitude.to_masked_array()
        time = data.time.to_masked_array()

    va850_array, _, _, _ = get_data(va850_file)

    nlon, nlat, ntime = ua850_array.shape

    deltax = haversine(lats, lons[50], lats, lons[51])
    deltay = haversine(lats[35], lons, lats[36], lons)
    deltay = np.average(deltay)

    vorticity = vorc.relative_vorticity(m=nlon,
                                        n=nlat,
                                        o=ntime,
                                        u=ua850_array,
                                        v=va850_array,
                                        deltax=deltax,
                                        deltay=deltay)

    return vorticity.swapaxes(0, 2)


def read_netcdfs(files, dim, dim1=None, values=None, tranf_func=None):
    paths = sorted(glob(files))
    open_kwargs = dict(decode_cf=True, decode_times=True, combine='by_coords', parallel=True)
    open_tasks = [dask.delayed(xr.open_mfdataset)(p, **open_kwargs) for p in paths]
    if tranf_func is None:
        def tranf_func(ds): return ds
    tasks = [dask.delayed(tranf_func)(task) for task in open_tasks]
    datasets = dask.compute(tasks)
    combined = xr.concat(datasets[0], dim)
    if values is not None:
        combined[dim] = (dim, values)
        if dim1 is not None: combined = combined.transpose(dim1, dim, ...)
    return combined

def process_data(year, src=None):
    print(f"================================= PROCESSING YEAR: {year} ====================================")
    print("<< Start processing >>>")
    ########################################### Load Data ############################################
    print("Loading data...")
    to = datetime.now()
    # Set working directory
    WRKSPC = os.path.dirname(os.path.abspath(sys.argv[0])) if src is None else src
    assert os.path.exists(WRKSPC), "Directory {} do not exists.".format(WRKSPC)
    os.chdir(WRKSPC)
    DATADIR = os.path.join(WRKSPC, f'Historical/{year}/projected')  ##change for RCP
    assert os.path.exists(DATADIR), "Directory {} do not exists.".format(DATADIR)
    RESDIR = os.path.join(DATADIR, 'results')
    if not os.path.exists(RESDIR):
        os.mkdir(RESDIR)

    # Set handler for netCDF files
    # 1. East wind, north wind, temperature and vorticity
    filepath = [glob(os.path.join(DATADIR, f"{tag}*_6hr_{year}*_proj.nc"))[0]
                    for tag in ['ta300', 'ta500', 'ta700', 'ta850']]
    temperature = xr.open_mfdataset(filepath, combine='by_coords', compat='override', parallel=True)
    filepath = [glob(os.path.join(DATADIR, f"{tag}*_6hr_{year}*_proj.nc"))[0]
                    for tag in ['ua300', 'ua850', 'va300', 'va850']]
    velocities = xr.open_mfdataset(filepath, combine='by_coords', compat='override', parallel=True)
    # 2. Surface pressure, sea surface temperature and 10 east and north wind
    #nc_forcing = xr.open_mfdataset(os.path.join(DATADIR, 'ERA5_sp_sst_u10_v10_pr.nc'), parallel=True)
    filepath = [glob(os.path.join(DATADIR, f"{tag}*_6hr_{year}*_proj.nc"))[0]
                    for tag in ['uas', 'vas']]
    velocities_as = xr.open_mfdataset(filepath, combine='by_coords', parallel=True)
    filepath = [glob(os.path.join(DATADIR, f"{tag}*_6hr_{year}*_proj.nc"))[0]
                    for tag in ['pr',]]
    precipitance = xr.open_mfdataset(filepath, combine='by_coords', parallel=True)

    #to= datetime.now()
    # Extract latitudes, longitudes, atmospheric levels and length of time
    lat = temperature.latitude
    lon = temperature.longitude
    tframe = temperature.time

    # Extract values from forcing data netCDF files
    # Climatology data

    #t = temperature.t
    data_vars = list(temperature.data_vars)
    data = np.array([temperature[vari] for vari in data_vars]).swapaxes(0,1)
    t = xr.DataArray(data=data,
                     coords={
                         'longitude': temperature.coords.get('longitude'),
                         'latitude': temperature.coords.get('latitude'),
                         'level': np.array([300, 500, 700, 850], dtype=np.int32),
                         'time': temperature.coords.get('time'),
                     },
                     dims=['time', 'level', 'latitude', 'longitude'])

    #ua = velocities.u
    data_vars = list(velocities.data_vars)
    data = np.array([velocities[vari] for vari in data_vars if vari.startswith('ua')]).swapaxes(0,1)
    ua = xr.DataArray(data=data,
                     coords={
                         'longitude': velocities.coords.get('longitude'),
                         'latitude': velocities.coords.get('latitude'),
                         'level': np.array([300, 850], dtype=np.int32),
                         'time': velocities.coords.get('time'),
                     },
                     dims=['time', 'level', 'latitude', 'longitude'])

    #va = velocities.v
    data = np.array([velocities[vari] for vari in data_vars if vari.startswith('va')]).swapaxes(0,1)
    va = xr.DataArray(data=data,
                     coords={
                         'longitude': velocities.coords.get('longitude'),
                         'latitude': velocities.coords.get('latitude'),
                         'level': np.array([300, 850], dtype=np.int32),
                         'time': velocities.coords.get('time'),
                     },
                     dims=['time', 'level', 'latitude', 'longitude'])

    filepath = sorted(glob(os.path.join(DATADIR, f"*[u,v]a850*_6hr_{year}*_proj.nc")))
    vo = calc_vort(*filepath)
    vo = xr.DataArray(data=vo,
                     coords={
                         'longitude': velocities.coords.get('longitude'),
                         'latitude': velocities.coords.get('latitude'),
                         'time': velocities.coords.get('time'),
                     },
                     dims=['time', 'latitude', 'longitude'])

    # Forcing data
    u10 = velocities_as.uas
    v10 = velocities_as.vas
    pr = precipitance.pr * 3.6  # (mm/hr)/(Kg/m2/s)

    filepath = sorted(glob(os.path.join(DATADIR, f"*[t,z]?850*_6hr_{year}*_proj.nc")))
    msl = calc_msl(*filepath)
    msl = xr.DataArray(data=msl,
                     coords={
                         'longitude': velocities.coords.get('longitude'),
                         'latitude': velocities.coords.get('latitude'),
                         'time': velocities.coords.get('time'),
                     },
                     dims=['time', 'latitude', 'longitude'])
    sp = msl

    # Extract values from temperatures in time range from to to tf
    # dimensions -> [time, level, latitude, longitude]
    print("Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    ######################################### Initialization ##########################################
    print("Initializing...", end=' ')
    to = datetime.now()
    # Setting limits for loops, ... etc.
    lonmax = len(lon)
    latmax = len(lat)

    # Extract 850hPa vorticity
    #vo850 = vo.where(vo.level == 850, drop=True).squeeze()
    vo850 = vo.squeeze()  # values for level 850 hPa

    # 10m wind speed
    w10 = np.sqrt(u10**2 + v10**2)

    # Wind speed on pressure levels
    wa = np.sqrt(ua**2 + va**2)

    # Extract wind at 300hPa
    wa300 = wa.where(wa.level == 300, drop=True).squeeze()

    # Extract wind at 850hPa
    wa850 = wa.where(wa.level == 850, drop=True).squeeze()

    # Compute anomaly
    t_357 = t[:,:3,:,:]
    print("Possible TCs: {}".format(len(vo850)))
    print("Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    ###################################### Detection Algorithmn #######################################
    # ===========================================STEP1=================================================
    # The 850-hPa relative vorticity exceeds the vorticity threshold
    # =================================================================================================
    print("Detection Algorithmn...")
    print("   Step 1: Vorticity threshold exceedance...", end=' ', flush=True)
    to = datetime.now()
    threshold_vo850 = 5.8e-04   # [s**-1]
    F1 = pd.DataFrame(np.transpose(np.array(np.where(vo850 > threshold_vo850))), columns=['TIME_index','LAT_index','LON_index'])

    # Ensure every cell chosen has a box of 7x7 around it
    F1 = F1[(F1.LAT_index>2) & (F1.LAT_index<(latmax-3)) & \
       (F1.LON_index>2) & (F1.LON_index<(lonmax-3))].reset_index(drop=True)

    print("   Possible TCs: {}".format(len(F1)))
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ===========================================STEP2=================================================
    # The maximum surface wind speed in a centered 7X7 box exceeds the wind speed threshold
    # =================================================================================================
    print("   Step 2: Wind speed threshold exceedance...", end=' ', flush=True)
    to = datetime.now()
    threshold_sws = 8.970555         #[m*s**-1]
    F1 = pth.thres_step(w10.values, F1, threshold_sws, 'max')
    nop = np.argwhere(F1['sub'].values == 0)
    F2 = F1.drop(nop.flatten()).reset_index(drop=True)
    print("   Possible TCs: {}".format(len(F2)))
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ===========================================STEP3=================================================
    # The sea level pressure is the minimum in a centered 7X7 box
    # ===========================================STEP3=================================================
    print("   Step 3: Centered minimum of sea level pressure...", end=' ', flush=True)
    to = datetime.now()
    F2 = pth.center_step(msl.values, F2)
    nop = np.argwhere(F2['sub'].values == 0)
    F3 = F2.drop(nop.flatten()).reset_index(drop=True)
    F3.insert(3, 'LAT', -999.0)
    F3.insert(4, 'LON', -999.0)
    filler = np.zeros((2, len(F3)), dtype=np.float64)
    for i in range(len(F3)):
         filler[0, i] = float(lon[F3['LON_index'][i]])
         filler[1, i] = float(lat[F3['LAT_index'][i]])
    F3['LON'] = filler[0, :]
    F3['LAT'] = filler[1, :]
    del filler
    print("   Possible TCs: {}".format(len(F3)))
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ===========================================STEP4=================================================
    # The temperature anomaly averaged over the 7X7 box and three pressure levels
    # (300, 500, and 700 hPa) exceeds the temperature anomaly threshold.
    # ===========================================STEP4=================================================
    print("   Step 4: Temperature anomaly averaged and three pressure levels threshold exceedance...", end=' ', flush=True)
    to = datetime.now()
    threshold_T = 1.091           #[Â°C]
    F3 = pth.thres_step_g(t_357.values, F3, threshold_T)
    nop = np.argwhere(F3['sub'].values == 0)
    F4 = F3.drop(nop.flatten()).reset_index(drop=True)
    print("   Possible TCs: {}".format(len(F4)))
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ===========================================STEP5=================================================
    # The local temperature anomaly averaged over the 7X7 box is positive at
    # all three pressure levels (300, 500, and 700 hPa).
    # ===========================================STEP5=================================================
    print("   Step 5: Local temperature anomaly averaged is positive at three pressure levels...", end=' ', flush=True)
    to = datetime.now()
    F4 = pth.levels_positive_step_g(t_357.values, F4)
    nop = np.argwhere(F4['sub'].values == 0)
    F5 = F4.drop(nop.flatten()).reset_index(drop=True)
    print("   Possible TCs: {}".format(len(F5)))
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ===========================================STEP6=================================================
    # The local temperature anomaly, averaged over the 7X7 box, at 300 hPa
    # is greater than at 850 hPa.
    # ===========================================STEP6=================================================
    print("   Step 6: Local temperature anomaly averaged is greater between two pressure levels...", end=' ', flush=True)
    to = datetime.now()
    F5 = pth.levels_compare_step_g(t.values, F5)
    nop = np.argwhere(F5['sub'].values == 0)
    F6 = F5.drop(nop.flatten()).reset_index(drop=True)
    print("   Possible TCs: {}".format(len(F6)))
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ===========================================STEP7=================================================
    # The mean speed averaged over a 7X7 grid box is larger at 850 hPa than at 300 hPa.
    # ===========================================STEP7=================================================
    print("   Step 7: Mean speed averaged is larger between two pressure levels...", end=' ', flush=True)
    to = datetime.now()
    F6 = pth.levels_compare_step(wa850, wa300, F6)
    nop = np.argwhere(F6['sub'].values == 0)
    F7 = results = F6.drop(nop.flatten()).reset_index(drop=True)
    print("   Possible TCs: {}".format(len(F7)))
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ========================================Max Wind Speed===========================================
    # Maximum sustained wind speed over 7x7 grid box
    # ========================================Max Wind Speed===========================================
    print("   Maximum sustained wind speed over 7x7 grid box...", end=' ', flush=True)
    to = datetime.now()
    wmax_grd = pth.max_value_grd(w10, F7)
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ================================Maximum and Sum Precipiction=====================================
    # Maximum precipiction over 7x7 grid box
    # ================================Maximum and Sum Precipiction=====================================
    print("   Maximum precipiction over 7x7 grid box...", end=' ', flush=True)
    to = datetime.now()
    prmax_grd = pth.max_value_grd(pr, F7)
    prsum_grd = pth.max_value_grd(pr, F7, 'sum')
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ===========================================Analysis==============================================
    # Checking effectiveness of filters
    # ===========================================Analysis==============================================
    print("Plotting for checking effectiveness of filters...")
    to = datetime.now()
    #nop = [eval('len(F{})'.format(i+1)) for i in range(7)]
    nop = []
    for i in range(7):
        nop.append(eval('len(F{})'.format(i+1)))

    data = {
        'No.P': nop,
        'Name': [i+1 for i in range(7)],
        'Step': ['S{}'.format(i+1) for i in range(7)]
    }
    analysis = pd.DataFrame(data, columns=['Name', 'No.P', 'Step'])

    results = pth.fill_results(results, tframe.values, msl.values, sp.values)
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # =====================================Data Table================================================
    # Creating data table for selected TCs
    # =====================================Data Table================================================
    print("Creating data table for selected TCs ...")
    to = datetime.now()
    data_table = pth.fill_table(results, w10.values, vo850.values, t.values, wmax_grd, prmax_grd, prsum_grd)
    filename = os.path.join(RESDIR, f"result_table_{year}.csv")
    data_table.to_csv(filename, index=False)
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    # ==================================Plotting-Chart==============================================
    # Plotting for check up , eye pressure ,....etc
    # ==================================Plotting-Charts==============================================
    print("Plotting for check up, eye pressure, etc...")
    to = datetime.now()
    ax = pth.plot_bars_filters(analysis)
    filename = os.path.join(RESDIR, f'bars_filters_{year}.png')
    ax.figure.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0.5)
    plt.close(ax.figure)
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    to = datetime.now()
    ax = pth.plot_points_pressure(results)
    filename = os.path.join(RESDIR, f'points_pressure_{year}.png')
    ax.figure.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0.5)
    plt.close(ax.figure)
    print("   Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))

    #plt.show()

    print("Total time elapsed: {} sec.".format((datetime.now()-t1).total_seconds()))
    print("<<============= FINISHED ==============>>")

## Ends load libraries and defining functions
print("Time elapsed: {} sec.".format((datetime.now()-to).total_seconds()))


if __name__ == '__main__':
    src = None ##change - location
    years = [1990,1991,1992] ##change - years
    for year in years:
        process_data(year, src)
