import numpy as np
import xarray as xr
import datetime as dt

from level3.watervapor.src import wvfcts as wvfct
import importlib as imp
imp.reload(wvfct)

import pandas as pd

def run_retrieval_ground(radar, thermo, lut, info, write=False):
    '''
    this function runs the retrieval introduced by Roy et al 2018, applied in Millan et al 2024
    input:
    - radar: e.g. level 2 xarray dataset; needs to include Ze1, Ze2, DAR with time, height dimension (ie could also be level-0 input...)
    - thermo: thermo xarray dataset, needs to include profiles of T, p with time, height dimension (e.g. radiosonde dataset, era-5, etc)
    - lut: xarray dataset with look-up table
    
    '''
    
    print('WV retrieval...preparing soundings')
    #rename thermo time dimension and convert launchtime to datetime:
    nt = thermo.launchtime.shape[0]
    thermo = thermo.assign_coords( launchtime=np.asarray([dt.datetime.strptime(str(thermo.launchtime.values[t]), '%Y%m%d%H%M') for t in range(nt)],dtype=object))
    thermo = thermo.rename({'launchtime':'time'})
    
    #interpolate thermo onto radar height grid
    thermoradar = thermo.interp(height=radar.height, time=radar.time)
    
    
    print('WV retrieval...calculating kappa')
    #get differential kappa from LUT based on thermo set of temperature and pressure, output unit: kg m**-2
    diffkappa = wvfct.get_kappa_from_lut(thermoradar.temp, thermoradar.pres, lut)
    
    print('preparing Ze and DAR with quality control')
    #set Ze and Ze2 and DAR to nan where respective SNR is below threshold
    radar['GZe'].values[radar['SNRG'].values < float(info['watervaporsettings']['SNRthresh'])] = np.nan
    radar['G2Ze'].values[radar['SNRG'].values < float(info['watervaporsettings']['SNRthresh'])] = np.nan #for now, set 174.7 channel to nan when 167 is attenuated
    radar['DAR'].values[radar['SNRG'].values < float(info['watervaporsettings']['SNRthresh'])] = np.nan
    #radar['G2Ze'][radar['SNRG2'] < info['watervaporsettings']['SNRthresh']] = np.nan #once part of matlab processing
    #radar['DAR'][(radar['SNR'] < info['watervaporsettings']['SNRthresh']) | (radar['SNRG2'] < info['watervaporsettings']['SNRthresh'])] = np.nan

    #loop through different R, calculate output profile and save output
    for R in np.asarray([info['watervaporsettings']['R']],dtype=int):
        print('WV retrieval...calculating gamma and retrieving with R=%i'%R)
        #get differential gamma calculated based on reflectivities, output unit: 1/m; also get number of averaged volumes per height level
        diffgamma, nvolumes, deltaZerR = wvfct.get_diffgamma(R, radar, outall=False)
        
        #calculate profile based on differential gamma, and differential kappa
        wvprof = np.divide(diffgamma, diffkappa)*1000.  #diffgamma units: 1/m ; diffkappa: m**2 kg**-1 ; wvprof: g m**(-3)
        
        #calculate rhov uncertainty according to Battaglia and Kollias 2019 Eq 7; Roy et al 2018 Eqn 13
        deltarhov = 1/(2 * R * diffkappa) * deltaZerR * 1000. #convert to g m-3
        

        
        #prepare variables to be stored in dataset:
        exportvars = [diffkappa, diffgamma, nvolumes, wvprof, R, radar.time.values, radar.height.values, deltarhov]
        varnames = ['diffkappa', 'diffgamma', 'nranges', 'rhov', 'R', 'time', 'height', 'deltarhov']
        tods = {}
        for vv,var in enumerate(exportvars):
            tods[varnames[vv]] = var

        #create dataset:
        attrs = radar.attrs.copy()
        wvxrds = wvfct.create_wv_dataset(tods, attrs)

        #filter obtained profiles by how many range gates were averaged together; ie set all profile output to nan where nranges is smaller than minnranges given in info json
        wvxrds.rhov.values[:,np.where(wvxrds.nranges < int(info['watervaporsettings']['minNranges']))] = np.nan


        #interpolate in height and time and smooth:
        print('resampling on %s time interval by taking mean of window...'%info['watervaporsettings']['tavg'])

        wvxrdssmooth = wvxrds.resample(time='%s'%info['watervaporsettings']['tavg'],closed='right').mean() #put timestamp at end of averaging period to be consistent with radar time stamping (rpg puts timestamp at end of one measurement)
        
        print('smoothing output vertically on R=%i by using rolling mean'%R)
        wvsmoothed = wvxrdssmooth.rhov.rolling( height= int( info['watervaporsettings']['Nrolling']), min_periods = int(info['watervaporsettings']['minNranges'])).mean()
        wvxrdssmooth.rhov.values = wvsmoothed.values
        
        #print('interpolating smoothing output to regular grid with R')
        #Rgrid = np.arange(0, 11000, R)
        #wvxrdssmooth = wvxrdssmooth.interp( height = Rgrid)
        1/0

        #store to netcdf
        if write==True:
            
            #output of original profiles:
            ncfile = info['global']['mission'] +'_%s_level3b_watervapor_%s%s%s_R%i_tavg%s_original.nc'%(info['global']['version'], info['yyyy'], info['mm'], info['dd'], R, info['watervaporsettings']['tavg'])
            wvxrds.to_netcdf(info['paths']['output'] + ncfile)
            
            #output of smoothed profiles:
            ncfile = info['global']['mission'] +'_%s_level3b_watervapor_%s%s%s_R%i_tavg%s.nc'%(info['global']['version'], info['yyyy'], info['mm'], info['dd'], R, info['watervaporsettings']['tavg'])
            wvxrdssmooth.to_netcdf(info['paths']['output'] + ncfile)
            
        
    return wvxrds, wvxrdssmooth

