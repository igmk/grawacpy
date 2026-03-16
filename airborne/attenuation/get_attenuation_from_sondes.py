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
    for key in ['tdry', 'rh', 'pres']:
        fctinter = interp.interp1d(sondexr['gpsalt'], sondexr[key])
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



def save_output_to_nc(pam, freqs, attrs, ncout, write = False):
    
    
    data_vars = {
        'Att_atmo':(['height','freq'], pam.r['Att_atmo'][0,0,:,:], {'units':'dB', 'long_name': 'atmospheric attenuation profile'}),
        'relh':(['height'], pam.p['relhum'][0,0,:], {'units':'', 'long_name':'layer relative humidity'}),
        'p':(['height'], pam.p['press'][0,0,:]/100, {'units':'hPa', 'long_name':'layer pressure'}),
        'temp':(['height'], pam.p['temp'][0,0,:], {'units':'K', 'long_name':'layer temperature'})
   }
    
    
    coords = {
        'freq':(['freq'], freqs, {'units':'GHz', 'long_name':'simulated radar frequency'}),
        'height':(['height'], pam.p['hgt'][0,0,:], {'units':'m', 'long_name':'layer height above msl'}),}
    
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


def run(info):

    print('calculating two-way attenuation from dropsonde profiles %s %s'%(info['global']['mission'], info['global']['RF']))    

    #load sounding profiles ===============================
    sondefiles = sorted(glob.glob(info['paths']['thermo']+'%s_%s%s%s_*/'%(info['global']['mission'], info['global']['yyyy'], info['global']['mm'], info['global']['dd'])+'*_??????QC.nc'))
    nsondes = len(sondefiles)
    print('found %i sondefiles'%nsondes)

    #get ERa5 files:
    erafiles = sorted(glob.glob(info['paths']['era']+'/era5_%s*%s*.nc'%(info['global']['mission'], info['global']['RF'])))
    nera = len(erafiles)
    print('found %i erafiles'%nera)
    print('checking that we have era5 profiles for each dropsonde...')
    assert nera == nsondes

    #load one gband radar file for height grid:
    #gpath = '/work/schnitts/gband/hamag/joint_dataset/'
    #l1 = xr.open_dataset(info['paths']['level1']['output'])

    #set up pamtra simulation: ============================
    radarfreq = [94.1, 167.32, 174.7]

    pam = pyPamtra.pyPamtra()
    descriptorFile = getDescriptor()
    for hyd in descriptorFile:
        pam.df.addHydrometeor(hyd)

    #some general settings:
    pam.nmlSet['passive'] = False
    pam.nmlSet['active'] = True
    pam.nmlSet["radar_attenuation"] = 'bottom-up'

    pamData = {}

    #run pamtra simulation for each sounding; save gas attenuation output in netcdf file
    for n in range(nsondes):
        
        print('modelling sounding %i'%n)
        
        #vertical grid that attenuation and soundings will be interpolated onto (arbitrary)
        arbheight = np.arange(20,6000,20) #this one is the one that will be used by pamtra
        arbheightatt = arbheight.copy() #this one is the output one filled with nans later if necessary
        #read sounding:
        sondexr = xr.open_dataset(sondefiles[n])
        
        adjustflag = 0 #flag indicating whether arbheight was adjusted for pamtra simulation
        #try interpolating sounding to output attenuation grid: 
        try:
            sondeinter = interpolate_sounding_to_radarheight(sondexr, arbheight)
        except ValueError: #if running into an interpolation error: adjust arbheight
            adjustflag = 1
            minheight = np.nanmin(sondexr.gpsalt.values)
            maxheight = np.nanmax(sondexr.gpsalt.values)

            if minheight > np.nanmin(arbheight): 
                print('DS min height is higher than attenuation grid. skipping.')
                continue
            #if DS minheight is larger than min of new grid: set new grid startpoint to first index above DS min height
            #if minheight > np.nanmin(arbheight):
            #    newind = np.where(minheight > arbheight)[0] #find the first index in arbheight thats larger than minheight
            #    arbheight = arbheight[newind:]
            #if maxheight < np.nanmax(arbheight):
            #    newind = np.where(maxheight < arbheight)[-1] #last index of DS height thats smaller than arbheight
            #    arbheight = arbheight[:newind]
            
            #sondeinter = interpolate_sounding_to_radarheight(sondexr, arbheight)
            
        #still complement layer RH in sondeinter as ERA does not have rh_lev as variable:
        sondeinter['rh_lay'] = (sondeinter['rh'][1:] + sondeinter['rh'][:-1])/2.  #average two levels to one layer
        
        #read matching ERA5 profile:
        eraxr = xr.open_dataset(erafiles[n])
        #interpolate ERA5 profile to arbheight height grid:
        erainter = interpolate_era5_to_radarheight(eraxr, arbheight)
        
        #varPairs = [['p', 'press'], ['p_lev', 'press_lev'], ['t', 'temp'], ['t_lev', 'temp_lev'], ['rh', 'relhum'], ['hgt_lev', 'hgt_lev']]
        varpairs = [['p_lev', 'pres'], ['t_lev', 'tdry'], ['rh', 'rh_lay']]
        
        #find last valid measurement in each pres, tdry; take ERA5 level measurements to complement the array:
        for erakey, sondekey in varpairs:
            
            lastind = np.where(np.isnan(sondeinter[sondekey]))[0][0]
            print(erakey, sondekey, lastind)
            if lastind == 0:            
                lastind = np.where((np.isnan(sondeinter[sondekey])) & (sondeinter['hgt'] > 1000))[0][0]
                #print(lastind)
                1/0
            
            sondeinter[sondekey][lastind:] = erainter[erakey][lastind:]
        
        #prepare pamData dictionnary:
        # define levels:
        pamData['hgt_lev'] = sondeinter['hgt']
        pamData['press_lev']= sondeinter['pres']*100 #in Pa
        pamData['temp_lev']= sondeinter['tdry'] + 273.15 #in K
        
        #define layers (as average between sounding measurements)
        pamData['temp'] = (sondeinter['tdry'][1:] + sondeinter['tdry'][:-1])/2. + 273.15  #in K
        pamData['press'] = (sondeinter['pres'][1:] + sondeinter['pres'][:-1])/2. * 100.  #in Pa
        pamData['relhum'] = sondeinter['rh_lay']    # in %
        
        #create pamtra profile:
        pam.createProfile(**pamData)
        
        #run pamtra:
        pam.runParallelPamtra(radarfreq, pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers=8)
        
        #prepare output and write to netcdf:
        attrs = info['attrs']['attenuation']
        attrs['sonde_launch_time'] = '%s'%np.datetime_as_string(sondexr.launch_time.values, unit='s')
        print('add sondetimestamp as coordinate to file')
        1/0
        ncoutfile = info['paths']['attenuation'] + '%s/%s/%s/'%(info['global']['yyyy'], info['global']['mm'], info['global']['dd']) + info['global']['mission'] + '_'+ info['global']['RF'] + '_sonde0%i_attenuation.nc'%(n+1)
        print('saving output to %s'%ncoutfile)
        
        xrds = save_output_to_nc(pam, radarfreq, attrs, ncoutfile, write=True)
    
    return



    
