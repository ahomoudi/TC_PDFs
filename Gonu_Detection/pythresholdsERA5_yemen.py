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
import sys
import os
from calendar import month_name as cal_month_name
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
import xarray as xr
import pythreshold_support as pth
import geopandas as gpd
from shapely.geometry import Point, Polygon

print("<< Start processing >>>")
########################################### Load Data ############################################
print("Loading data...")
# Set working directory
WRKSPC = '/mnt/Shared/000Documentos/001_Master_HSE/FRM_cert/Scripts/Gonu/Fiverr'
#WRKSPC = os.path.dirname(os.path.abspath(sys.argv[0])) if len(sys.argv) == 1 else sys.argv[1]
#assert os.path.exists(WRKSPC), "Directory {} do not exists.".format(WRKSPC)
os.chdir(WRKSPC)
DATADIR = os.path.join(WRKSPC, 'data')

# Set handler for netCDF files
# 1. East wind, north wind, temperature, vorticity and humidity
nc_climatology = xr.open_mfdataset(os.path.join(DATADIR, 'ERA5_u_v_t_vo_q.nc'), parallel=True)
# 2. Surface pressure, sea surface temperature and 10 east and north wind
nc_forcing = xr.open_mfdataset(os.path.join(DATADIR, 'ERA5_sp_sst_u10_v10_pr.nc'), parallel=True)
# 3. Temperatures at May to July, each 6 hours average for 30 years (01-05-2007 to 31-07-2007)
nc_temperatures = xr.open_mfdataset(os.path.join(DATADIR, 'ERA5_t_6hrs30yrs_MAY_JULY.nc'), parallel=True)

# Extract latitudes, longitudes, atmospheric levels and length of time
lat = nc_climatology.latitude
lon = nc_climatology.longitude
lev = nc_climatology.level
tframe = nc_climatology.time

# Extract values from forcing data netCDF files
# Climatology data
q = nc_climatology.q
t = nc_climatology.t
ua = nc_climatology.u
va = nc_climatology.v
vo = nc_climatology.vo

# Forcing data
u10 = nc_forcing.u10
v10 = nc_forcing.v10
msl = nc_forcing.msl
sst = nc_forcing.sst
sp = nc_forcing.sp
pr = nc_forcing.tp

# Extract values from temperatures in time range from to to tf
# dimensions -> [time, level, latitude, longitude]
to = 124  # Corresponds to date 01-06-2007
tf = to + tframe.shape[0]  # Initial date plus time window of input data (climatology + forcing)
t30yrs = nc_temperatures.t[to:tf,:,:,:]

######################################### Initialization ##########################################
print("Initializing...")
# Setting limits for loops, ... etc.
lonmax = len(lon)
latmax = len(lat)

assert t30yrs.time.shape == t.time.shape, "Length of timeseries between climatology and temperature data are not equals."
t30yrs.coords['time'] = t.coords['time']


# Extract 850hPa vorticity
vo850 = vo.where(vo.level == 850, drop=True).squeeze()

# 10m wind speed
w10 = np.sqrt(u10**2 + v10**2)

# Wind speed on pressure levels
wa = np.sqrt(ua**2 + va**2)

# Extract wind at 300hPa
wa300 = wa.where(wa.level == 300, drop=True).squeeze()

# Extract wind at 850hPa
wa850 = wa.where(wa.level == 850, drop=True).squeeze()

# Compute anomaly
t_anomaly = t - t30yrs
t_anomaly357 = t_anomaly[:,:3,:,:]

###################################### Detection Algorithmn #######################################
# ===========================================STEP1=================================================
# The 850-hPa relative vorticity exceeds the vorticity threshold
# =================================================================================================
print("Detection Algorithmn...")
print("   Step 1: Vorticity threshold exceedance...")
threshold_vo850 = 0.0001816065   # [s**-1]
F1 = pd.DataFrame(np.transpose(np.array(np.where(vo850 > threshold_vo850))), columns=['TIME_index','LAT_index','LON_index'])

# Ensure every cell chosen has a box of 7x7 around it
nop = np.argwhere(F1['LON_index'].values <= 3)
F1 = F1.drop(nop.flatten()).reset_index(drop=True)

nop = np.argwhere(F1['LON_index'].values >= lonmax-5)
F1 = F1.drop(nop.flatten()).reset_index(drop=True)

nop = np.argwhere(F1['LAT_index'].values <= 3)
F1 = F1.drop(nop.flatten()).reset_index(drop=True)

nop = np.argwhere(F1['LAT_index'].values >= latmax-5)
F1 = F1.drop(nop.flatten()).reset_index(drop=True)

# ===========================================STEP2=================================================
# The maximum surface wind speed in a centered 7X7 box exceeds the wind speed threshold
# =================================================================================================
print("   Step 2: Wind speed threshold exceedance...")
threshold_sws = 6.670151         #[m*s**-1]
F1 = pth.thres_step(w10, F1, threshold_sws, 'max')
nop = np.argwhere(F1['sub'].values == 0)
F2 = F1.drop(nop.flatten()).reset_index(drop=True)

# ===========================================STEP3=================================================
# The sea level pressure is the minimum in a centered 7X7 box
# ===========================================STEP3=================================================
print("   Step 3: Centered minimum of sea level pressure...")
F2 = pth.center_step(msl, F2)
nop = np.argwhere(F2['sub'].values == 0)
F3 = F2.drop(nop.flatten()).reset_index(drop=True)
F3.insert(3, 'LAT', -999.0)
F3.insert(4, 'LON', -999.0)
filler = np.zeros((2, len(F3)), dtype=np.int0)
for i in range(len(F3)):
     filler[0, i] = float(lon[F3['LON_index'][i]])
     filler[1, i] = float(lat[F3['LAT_index'][i]])
F3['LON'] = filler[0, :]
F3['LAT'] = filler[1, :]
del filler
# ===========================================STEP4=================================================
# The temperature anomaly averaged over the 7X7 box and three pressure levels
# (300, 500, and 700 hPa) exceeds the temperature anomaly threshold.
# ===========================================STEP4=================================================
print("   Step 4: Temperature anomaly averaged and three pressure levels threshold exceedance...")
threshold_T = 0.3292001            #[Â°C]
F3 = pth.thres_step(t_anomaly357, F3, threshold_T, 'mean')
nop = np.argwhere(F3['sub'].values == 0)
F4 = F3.drop(nop.flatten()).reset_index(drop=True)

# ===========================================STEP5=================================================
# The local temperature anomaly averaged over the 7X7 box is positive at
# all three pressure levels (300, 500, and 700 hPa).
# ===========================================STEP5=================================================
print("   Step 5: Local temperature anomaly averaged is positive at three pressure levels...")
F4 = pth.levels_positive_step(t_anomaly357, F4)
nop = np.argwhere(F4['sub'].values == 0)
F5 = F4.drop(nop.flatten()).reset_index(drop=True)

# ===========================================STEP6=================================================
# The local temperature anomaly, averaged over the 7X7 box, at 300 hPa
# is greater than at 850 hPa.
# ===========================================STEP6=================================================
print("   Step 6: Local temperature anomaly averaged is greater between two pressure levels...")
F5 = pth.levels_compare_step(t_anomaly357, t_anomaly, F5)
nop = np.argwhere(F5['sub'].values == 0)
F6 = F5.drop(nop.flatten()).reset_index(drop=True)

# ===========================================STEP7=================================================
# The mean speed averaged over a 7X7 grid box is larger at 850 hPa than at 300 hPa.
# ===========================================STEP7=================================================
print("   Step 7: Mean speed averaged is larger between two pressure levels...")
F6 = pth.levels_compare_step(wa850, wa300, F6)
nop = np.argwhere(F6['sub'].values == 0)
F7 = results = F6.drop(nop.flatten()).reset_index(drop=True)

# ===========================================Analysis==============================================
# Checking effectiveness of filters
# ===========================================Analysis==============================================
print("Plotting for checking effectiveness of filters...")
data = {
    'No.P': [eval('len(F{})'.format(i+1)) for i in range(7)],
    'Name': [i+1 for i in range(7)],
    'Step': ['S{}'.format(i+1) for i in range(7)]
}
analysis = pd.DataFrame(data, columns=['Name', 'No.P', 'Step'])

results = pth.fill_results(results, tframe, msl, sp)

# ==================================Plotting-Chart==============================================
# Plotting fo check up , eye pressure ,....etc
# ==================================Plotting-Charts==============================================
pth.plot_bars_filters(analysis)
pth.plot_poins_pressure(results)

#########################################################################
# read files from IBTrACS csv file
#   Rountine to read and plot IBTrACS data  of the Arabian Sea subbasain of
#   the North indian Ocean basin for the period 1990-2019
#########################################################################
print("Reading csv files with IBTrACS data...")
nabasin = pd.read_csv('data/ibtracs.NI.list.v04r00.csv', header=[0,1], low_memory=False)
nabasin.columns = [hdr[0] for hdr in nabasin.columns.copy()]

# extract month for dataset NA.basin
mom = list(
    map(
        lambda s: cal_month_name[int(s.split()[0][5:7])],
        nabasin['ISO_TIME']
    )
)
nabasin.insert(len(nabasin.columns), 'Month', mom)

# select storms GONU 2007 in AS SubBsasin
substorm = nabasin[nabasin['NAME'] == 'GONU']

# add and ID with name and season
mom = list(
    map(
        lambda d: "{}.{}".format(d[0],d[1]),
        zip(substorm['NAME'], substorm['SEASON'])
    )
)
substorm.insert(len(substorm.columns), 'ID', mom)

#substorm.dtypes

# ================================================================================
# Plotting of different features of GONU2007
# ================================================================================
print("Plotting features of cyclons...")
shapefile = "continent shapefile/continent.shp"
#pth.plot_shape(shapefile)

world = gpd.read_file(shapefile)
#world.geometry.total_bounds # Get spatial extent
#world.crs

colmap = ["USA_LON", "USA_LAT", "USA_WIND"]

substorm[colmap] = substorm[colmap].apply(pd.to_numeric) 
#print(substorm.dtypes)

geom = [Point(xy) for xy in zip(substorm['LON'],substorm['LAT'])]
geom_usa = [Point(xy) for xy in zip(substorm['USA_LON'],substorm['USA_LAT'])]
geom_detect = [Point(xy) for xy in zip(results['LON'],results['LAT'])]
crs = {'init': 'epsg:4326'}

tracks = gpd.GeoDataFrame(substorm['WMO_WIND'], geometry = geom, crs = crs)
tracks_usa = gpd.GeoDataFrame(substorm['USA_WIND'], geometry = geom_usa, crs = crs)
tracks_detect = gpd.GeoDataFrame(results, geometry = geom_detect, crs = crs)

# Create x and y min and max objects to use in the plot boundaries
boundaries = np.array([30,0,80,35])
xlim = ([boundaries[0],  boundaries[2]])
ylim = ([boundaries[1],  boundaries[3]])

fig, ax = plt.subplots(figsize=(15,7))
ax.set_xlim(xlim)
ax.set_ylim(ylim)

world.plot(alpha = .5, ax = ax)
tracks.plot(color='purple', ax=ax, alpha=.1)
tracks_usa.plot(color='green', ax=ax, alpha=.1)
tracks_detect.plot(color='red', ax=ax, alpha=.4, marker = '*')

ax.set(title='TC detection')

plt.show()


print("<<============= FINISHED ==============>>")
