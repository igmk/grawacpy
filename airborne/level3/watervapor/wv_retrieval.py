import numpy as np
import xarray as xr
import datetime as dt
import glob
import importlib as imp
from level3.watervapor.src import wvfcts as wvfct
imp.reload(wvfct)

import pandas as pd

def run_retrieval_air(radar, thermo, lut, info, write=False):
    '''
    this function runs the retrieval introduced by Roy et al 2018, applied in Millan et al 2024
    input:
    - radar: e.g. level 2 xarray dataset; needs to include Ze1, Ze2, DAR with time, height dimension (ie could also be level-0 input...)
    - thermo: here: only path to dropsondes. #####thermo xarray dataset, needs to include profiles of T, p with time, height dimension (e.g. radiosonde dataset, era-5, etc)
    - lut: xarray dataset with look-up table
    
    '''
    
    print('WV retrieval...preparing soundings')
    '''this is if thermo is all sounding profiles in one xarra ydataset with height as coordinate.
    #rename thermo time dimension and convert launchtime to datetime:
    nt = thermo.launchtime.shape[0]
    thermo = thermo.assign_coords( launchtime=np.asarray([dt.datetime.strptime(str(thermo.launchtime.values[t]), '%Y%m%d%H%M') for t in range(nt)],dtype=object))
    thermo = thermo.rename({'launchtime':'time'})
    
    #interpolate thermo onto radar height grid
    thermoradar = thermo.interp(height=radar.height, time=radar.time)
    '''
    
    sondefiles = np.array(sorted(glob.glob(thermo+'*_??????QC.nc')))
    #skip broken ones for RF05:
    #sondefiles = sondefiles[np.array([1,2,3])] #skip numbers 0 and 4
    
    #new altitudes:
    newalt = np.arange(10,2800, 5)

    allds = wvfct.get_all_soundings(sondefiles, newalt)

    #interpolate thermo onto radar height grid
    thermoradar = allds.interp(height=radar.height, time=radar.time)
    
    print('WV retrieval...calculating kappa')    
    
    #get differential kappa from LUT based on thermo set of temperature and pressure, output unit: kg m**-2
    ## important assumption here: T and p profile does not vary much between slant and nadir path.
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
        print('WV retrieval...calculating gamma along range dimension and retrieving with R=%i'%R)
        
        #get differential gamma calculated based on reflectivities, output unit: 1/m; also get number of averaged volumes per height level
        diffgamma, nvolumes, deltaZerR = wvfct.get_diffgamma(R, radar, outall=False)
        
        #calculate profile based on differential gamma, and differential kappa
        wvprof = np.divide(diffgamma, diffkappa)*1000.  #diffgamma units: 1/m ; diffkappa: m**2 kg**-1 ; wvprof: g m**(-3)
        
        #calculate rhov uncertainty according to Battaglia and Kollias 2019 Eq 7; Roy et al 2018 Eqn 13
        deltarhov = 1/(2 * R * diffkappa) * deltaZerR * 1000. #convert to g m-3
        
        #prepare variables to be stored in dataset:
        tavg = int(info['watervaporsettings']['tavg'].split('s')[0])
        exportvars = [diffkappa, diffgamma, nvolumes, wvprof, R, radar.time.values, radar.range.values, deltarhov, tavg, radar.height.values]
        varnames = ['diffkappa', 'diffgamma', 'nranges', 'rhov', 'R', 'time', 'range', 'deltarhov','tavg', 'height']
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
        wvsmoothed = wvxrdssmooth.rhov.rolling( range = int( info['watervaporsettings']['Nrolling']), min_periods = int(info['watervaporsettings']['minNranges'])).mean()
        wvxrdssmooth.rhov.values = wvsmoothed.values

        #project profiles on nadir:
        slantfactor = np.cos(np.radians(radar.incangle.interp(time=wvxrdssmooth.time)))
        slantmatrix = np.broadcast_to(slantfactor, (wvxrdssmooth.rhov.shape[1], len(slantfactor))).T
        wvprofnadir = slantmatrix * wvxrdssmooth.rhov

        #print('interpolating smoothing output to regular grid with R')
        #Rgrid = np.arange(0, 11000, R)
        #wvxrdssmooth = wvxrdssmooth.interp( height = Rgrid)
        
        ## now retrieve IWV
        kwargs = {}
        kwargs['groundsigma'] = 'max'        # max means: maximum Ze within specified altitude bins is used; mean means: mean is calculated
        #specify min and max alt used for calculating Ze used for retrieving iwv from
        kwargs['minalt'] = -300    
        kwargs['maxalt'] = 300  
        kwargs['option'] = 'sounding'
        
        zres = np.unique(np.diff(radar.range))[-1]
        freq1 = 167.3
        freq2 = 174.7
        iwv, diffkappaiwv = wvfct.retrieve_iwv(radar.GZe, radar.G2Ze, freq1, freq2, zres, radar.incangle, radar.height, radar.range, allds, lut, kwargs, wvxrds)

        wvxrds = wvxrds.assign(iwv=iwv)
        wvxrdssmooth = wvxrdssmooth.assign(iwv=iwv.interp(time=wvxrdssmooth.time), wvprofnadir=wvprofnadir)

        #store to netcdf
        if write==True:
            
            #output of original profiles:
            ncfile = info['global']['mission'] +'_%s_level3b_watervapor_%s%s%s_R%i_tavg%s_original.nc'%(info['global']['version'], info['yyyy'], info['mm'], info['dd'], R, info['watervaporsettings']['tavg'])
            wvxrds.to_netcdf(info['paths']['output'] + ncfile)
            
            #output of smoothed profiles:
            ncfile = info['global']['mission'] +'_%s_level3b_watervapor_%s%s%s_R%i_tavg%s.nc'%(info['global']['version'], info['yyyy'], info['mm'], info['dd'], R, info['watervaporsettings']['tavg'])
            wvxrdssmooth.to_netcdf(info['paths']['output'] + ncfile)
            
        
    return wvxrds, wvxrdssmooth
