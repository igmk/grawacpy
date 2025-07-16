import xarray as xr
import glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import pyPamtra
from scipy import interpolate as interp
import json
import os

'''
get gas attenuation profile by forward simulating dropsonde profiles from PAMTRA
'''

def getDescriptor():
    
    descriptorFile = np.array([
      #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as'
      # 'beta_as' # 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 
      #'d_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
      ('cwc_q', -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0, 3, 1, 
        'mono', -99.0, -99.0, -99.0, -99.0, 2e-05, -99.0, 'mie-sphere',
        'khvorostyanov01_drops', -99.0)], 
      dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), 
             ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), 
             ('alpha_as', '<f8'), ('beta_as', '<f8'), ('moment_in', '<i8'), 
             ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), 
             ('p_2', '<f8'), ('p_3', '<f8'), ('p_4', '<f8'), ('d_1', '<f8'), 
             ('d_2', '<f8'), ('scat_name', 'S15'), ('vel_size_mod', 'S30'), 
             ('canting', '<f8')])
    
    return descriptorFile 



def interpolate_sounding_to_radarheight(sondexr, gheight):
    '''
    interpolate sounding temperature, relative humidity and pressure to given height array
    input:
    - sondexr: xr dataset of sounding
    - gheight: height array [m] that sounding variables will be interpolated onto
    output:
    - sondeinter: dictionnary with interpolated tdry, rh, pres on hgt
    '''
    sondeinter = {}
    for key in ['temp', 'rh', 'pres']:
        fctinter = interp.interp1d(sondexr['z_gps'], sondexr[key])
        sondeinter[key] = fctinter(gheight)
    sondeinter['hgt'] = gheight
    return sondeinter


def interpolate_era5_to_radarheight(eraxr, gheight):
    '''
    interpolate era5 level temperature, relative humidity and pressure to given height array
    input:
    - eraxr: xr dataset of era5
    - gheight: height array [m] that sounding variables will be interpolated onto
    output:
    - erainter: dictionnary with interpolated t, rh, pres on hgt level array
    '''
    
    #extract height array era:
    #use levels:
    erahgt_lev = eraxr.hgt_lev.values[0,:,0,0]
    erahgt = eraxr.hgt.values[0,:,0,0]
    
    #dictionnary for interpolated variables:
    erainter = {}
    
    gheightlay = (gheight[:-1] + gheight[1:])/2.
    
    #interpolate level variables and swap dimensions in era variables as ERA starts at TOA
    for key in ['t_lev', 'p_lev']:
        #print(key)
        #print(erahgt_lev.shape, eraxr[key].shape)
        fctinter = interp.interp1d(erahgt_lev[::-1], eraxr[key][0,::-1,0,0])
        erainter[key] = fctinter(gheight)
    erainter['hgt'] = gheight
    
    #interpolate rh (which is layer) to new layer array:
    fctinter = interp.interp1d(erahgt[::-1], eraxr['rh'][0,::-1,0,0])
    erainter['rh'] = fctinter(gheightlay)
    
    #get units matched with sounding units:
    erainter['p_lev'] /= 100 #in hPa
    erainter['t_lev'] -= 273.15 #in deg C
    
    
    return erainter



def create_dataset(inputprofs, att_atmo, hgtlay, freqs, attrs, ncout, write = False):
    
    #launchtime = dt.datetime.strptime(attrs['sonde_launch_time'],'%Y%m%d%H%M')
    launchtime = int(attrs['sonde_launch_time'])
    
    coords = {
        'freq':(['freq'], freqs, {'units':'GHz', 'long_name':'simulated radar frequency'}),
        'height':(['height'], hgtlay, {'units':'m', 'long_name':'layer height above msl'}),
        'launchtime':(['launchtime'], [launchtime], {'units':'UTC', 'long_name':'sonde launch time UTC'})}
    
    data_vars = {
        'Att_atmo':(['launchtime','height','freq'], np.expand_dims(att_atmo,axis=0), {'units':'dB', 'long_name': 'atmospheric attenuation profile'}),
        'relh':(['height'], inputprofs['relhum'], {'units':'%', 'long_name':'layer relative humidity'}),
        'p':(['height'], inputprofs['press']/100, {'units':'hPa', 'long_name':'layer pressure'}),
        'temp':(['height'], inputprofs['temp'], {'units':'K', 'long_name':'layer temperature'})
   }
    
    
    
    attrs['creation_time'] = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
    #attrs = attrs
    
    xrds = xr.Dataset(data_vars, coords, attrs)
    
    if write==True:
        #check if ncoutput already available; if yes, delete it.
        if os.path.isfile(ncout):
            print('found old file. deleting old file.')
            os.remove(ncout)
        
        xrds.to_netcdf(ncout, 'w', format='NETCDF4')
    else:
        print('no output saved.')
    return xrds

#load information file:
with open("./attrs_paths_iop4h2o_new.json") as json_file:
    info = json.load(json_file)


#load sounding profiles ===============================
sondefiles = sorted(glob.glob(info['paths']['attenuation']['sondes']+'*.nc'))

#by taking level2 regulargrid soundings, they have no nans in interpolation area.

xrrs = xr.open_dataset(sondefiles[0])
nsondes = xrrs.launchtime.shape[0]
#nsondes = len(sondefiles)
print('found %i sondefiles'%nsondes)

'''
#get ERa5 files:
erafiles = sorted(glob.glob(info['paths']['attenuation']['era']+'/era5_%s*%s*.nc'%(info['global']['mission'], info['global']['RF'])))
nera = len(erafiles)
print('found %i erafiles'%nera)
print('checking that we have era5 profiles for each dropsonde...')
assert nera == nsondes
'''

#load one gband radar file for height grid:
#gpath = '/work/schnitts/gband/hamag/joint_dataset/'
l1 = xr.open_dataset(info['paths']['level1']['output']+'IOP4H2O_draco_v1.0_level1_20250206.nc')

#set up pamtra simulation: ============================
radarfreq = [35.5, 94.1, 167.32, 174.7]

pam = pyPamtra.pyPamtra()
descriptorFile = getDescriptor()
for hyd in descriptorFile:
    pam.df.addHydrometeor(hyd)

#some general settings:
pam.nmlSet['passive'] = False
pam.nmlSet['active'] = True
pam.nmlSet["radar_attenuation"] = 'bottom-up'

pamData = {}

#sondeswithnans = []
#run pamtra simulation for each sounding; save gas attenuation output in netcdf file

xrdsall = []
for n in range(nsondes):
    #for n in range(100):
    
    launchtimestring = str(xrrs.launchtime.values[n])
    
    #vertical grid that attenuation and soundings will be interpolated onto (arbitrary)
    arbheight = np.arange(10,30000,10)
    
    #read sounding:
    #sondexr = xr.open_dataset(sondefiles[n])
    minheight = np.nanmin(xrrs.height.values)
    maxheight = np.nanmax(xrrs.height.values)
    print('modelling sounding %s, min height: %.2f, max height: %.2f'%(launchtimestring, minheight, maxheight))
    
    #interpolate missing values
    
    #and extrapolate missing values:
    if np.any(np.isnan(xrrs.rh.values)) or np.any(np.isnan(xrrs.temp.values)):
        print('found %i nans in rh; and %i nans in temp'%(np.sum(np.isnan(xrrs.rh.values)), np.sum(np.isnan(xrrs.temp.values))))
        xrrs = xrrs.interpolate_na(dim='height',method='linear',fill_value='extrapolate')
    
    
    #prepare pamData dictionnary:
    # define levels:
    pamData['hgt_lev'] = xrrs['height'].values
    pamData['press_lev']= xrrs['pres'].values[n,:] #in Pa
    pamData['temp_lev']= xrrs['temp'].values[n,:] #in K
    relhum_lev = xrrs['rh'].values[n,:]
    
    #define layers (as average between levels measurements)
    pamData['temp'] = (pamData['temp_lev'][1:] + pamData['temp_lev'][:-1])/2.   #in K
    pamData['press'] = (pamData['press_lev'][1:] + pamData['press_lev'][:-1])/2.   #in K
    pamData['relhum'] = (relhum_lev[1:] + relhum_lev[:-1])/2.   #in %
    
    '''
    pamData['temp'] = (sondeinter['temp'].values[validinds][1:] + sondeinter['temp'].values[validinds][:-1])/2.   #in K
    pamData['press'] = (sondeinter['pres'].values[validinds][1:] + sondeinter['pres'].values[validinds][:-1])/2.  #in Pa
    pamData['relhum'] = (sondeinter['rh'].values[validinds][1:] + sondeinter['rh'].values[validinds][:-1])/2.   #%
    '''
    #validindslay = ~np.isnan(pamData['press'].values) & ~np.isnan(['temp'].values) & ~np.isnan(sondeinter['rh'].values)
    
    pam.createProfile(**pamData)
    
    pam.runParallelPamtra(radarfreq, pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=8)
    
    #pamtra output: only up to height where all measurements are available.
    
    #save to netcdf:
    '''attrs = {'author':'Sabrina Schnitt', 
             'title':'bottom-up gas attenuation profile calculated from HAMAG dropsonde profiles complemented by ERA5 data; calculated on arbitrary height grid', 
             'sonde_launch_time':'%s'%(np.datetime_as_string(sondexr.launch_time.values, unit='s'))}
    '''
    attrs = info['attrs']['attenuation']
    #attrs['sonde_launch_time'] = '%s'%np.datetime_as_string(sondexr.launch_time.values, unit='s')
    attrs['sonde_launch_time'] = '%s'%launchtimestring[:] # yyyymmddhhmm
    
    #artificially expand columns again to dimension of arbeight, ie fill up with nans in the attenuation to make sure that output has same height dimension for every sonde:
    '''
    arbheightlay = (arbheight[1:] + arbheight[:-1])/2.
    nz = len(arbheightlay)  #number of arbheight layers
    nf = len(radarfreq)
    #create empty atm attenuation array with dimensions that are always the same
    att_atmo = np.ones((nz, nf))
    
    #find layers in arbheightlay for which attenuation was calculated
    validlays = (arbheightlay >= pam.p['hgt'][0,0,0]) & (arbheightlay <= pam.p['hgt'][0,0,-1])
    valididx = np.where(validlays==True)[0]
    nanidx = np.where(validlays==False)[0]
    
    #save attenuation at correct layer indices; nan if attenuation was not calculated for that layer
    att_atmo[valididx,:] = pam.r['Att_atmo'][0,0,:,:]
    att_atmo[nanidx,:] = np.nan
    
    '''
    #also store input relh, press, temp profiles on arbheightlay dimension:
    nz = pam.p['nlyrs'][0][0]
    inputprofs = {}
    for k in ['temp','relhum','press']:
        inputprofs[k] = np.ones(nz)
        inputprofs[k] = pam.p[k][0,0,:]
        #inputprofs[k][nanidx] = np.nan
    att_atmo = pam.r['Att_atmo'][0,0,:,:]
    
    #ncoutfile = info['paths']['attenuation']['output'] + info['global']['mission'] +'_%s_attenuation.nc'%(launchtimestring[:])
    ncoutfile = info['paths']['attenuation']['output']  + info['global']['mission'] +'_%s_attenuation.nc'%(launchtimestring[:])
    print('saving output to %s'%ncoutfile)
    
    xrds = create_dataset(inputprofs, att_atmo, pam.p['hgt'][0,0,:], radarfreq, attrs, ncoutfile, write=True)
    #xrdsall.append(xrds)
    
    #concatenate along time dimension, save output
    
    
    
    
