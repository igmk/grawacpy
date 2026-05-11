import numpy as np
import datetime as dt
import xarray as xr
from scipy import interpolate as interp
from src import helpfunctions as helpfct

import importlib as imp
from src import thermo_conversion as tc
imp.reload(tc)

import metpy.calc as mpy 

def val2idx(c, value):
    
    nc = len(c)
    idx = np.round(np.interp(value, c, np.arange(nc)).clip(0,nc-1))
    
    return np.array(idx,dtype='int')


def create_wv_dataset(dsvars, attrs):
    '''
    create xarray dataset from input variables, return xr dataset
    '''
    coords = {'time':(['time'], dsvars['time']),
              'height':(['height'], dsvars['height'], {'units':'m', 'long_name':'Radar range gate height asl (radar range + instrument_altitude'}),
              'R':(['R'], [dsvars['R']], {'units':'m', 'long_name':'retrieval vertical resolution'})
             }
    
    #datavars to include: ze, ldr, startidx, freq, 
    data_vars = { 
    'rhov' : (['time','height'], dsvars['rhov'], {'units': 'g m-3', 'long_name':'absolute humidity profile retrieved from DAR measurements'}),
    'diffkappa' : (['time','height'], dsvars['diffkappa'], {'units': 'mm6 m-3', 'long_name':'radar reflectivity Ze'}),
    'diffgamma' : (['time','height'], dsvars['diffgamma'], {'units': 'm s-1', 'long_name':'radar mean Doppler Velocity'}),
    'nranges' : (['height'], dsvars['nranges'], {'units': '-', 'long_name':'number of range bins that fit into R spacing'})
        }
    
    attrs['title'] = 'Water vapor retrieval output GRaWAC'
    attrs['creation_time'] = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
    
    xrds = xr.Dataset(data_vars, coords, attrs)
    
    return xrds

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
        try:
            fctinter = interp.interp1d(sondexr['gpsalt'], sondexr[key])
            sondeinter[key] = fctinter(gheight)
        except ValueError:
            nansinkey = np.isnan(sondexr[key]) #get all nans in the array
            1/0
            
    sondeinter['hgt'] = gheight
    return sondeinter

def get_all_soundings(sondefiles, newalt):
    '''
    load and read all dropsondefiles, interpolate to regular grid newalt
    input: 
    - sondefiles: filenames including full path for dropsonde data
    - newalt: altitude vector that dropsonde profiles will be interpolated onto
    output:
    - xrds: xr dataset with sonde profiles
    '''
    
    
    nsondes = len(sondefiles)
    nz = len(newalt)
    
    #load each sounding profile; interpolate to regular grid of 10m; save
    allt = np.ones((nsondes, nz))
    allp = np.ones((nsondes, nz))
    allrhov = np.ones((nsondes, nz))
    dstimes = []
    dsiwv = []

    #from each sounding: get T and p profile, interpolate to common grid, copy to 2d array
    for n in range(nsondes):
        #print(n)
        dsdata = xr.open_dataset(sondefiles[n])
        
        #skip broken ones for RF05:
        #if n == 0 or n==4:
        #    dstimes.append(dsdata.launch_time.values)
        #    continue# or n==3 or n==4 or n==5: continue
        
        print('.....TODO: code workaround for interpolation routine: adjust interpolation boundaries such that valid interpolation possible, then fill rest up with nans')
        
        dsinter = interpolate_sounding_to_radarheight(dsdata, newalt)
        allt[n, :] = dsinter['tdry'] + 273.15
        allp[n, :] = dsinter['pres'] 
        allrhov[n, :] = tc.rh_to_qabs(dsinter['rh'], dsinter['tdry']+273.15)*1000.
        dstimes.append(dsdata.launch_time.values)  #also extract launch time as datetime object
        
        #get sounding iwv:
        pw = mpy.precipitable_water(dsdata.pres, dsdata.dp).magnitude
        dsiwv.append(pw)

    dstimes = np.array(dstimes)

    #now put all T and pressure profiles from the different soundings into a xr dataset:
    data_vars = {'temp':(['time', 'height'], allt), 'pres':(['time', 'height'], allp), 'rhov':(['time', 'height'], allrhov), 'iwv':(['time'], dsiwv)}
    coords = {'time':(['time'], dstimes), 'height':(['height'], newalt)}
    xrds = xr.Dataset(data_vars, coords)

    return xrds



def get_kappa_from_lut(temp, press, lut):
    '''
    for a set of temperature and pressure, return differential kappa 175-167
    input:
    - temp: temperature value in K
    - press: pressure value in hPa
    - lut: xr dataset read from look up table
    output:
    - diffkappa (175-167)
    '''
    
    nt = temp.shape[0]
    nz = temp.shape[1]

    kappa1 = np.ones((nt, nz))
    kappa2 = np.ones((nt, nz))
    for t in range(nt):
        #check if the whole column press or temperature entries are nan: if so, put nans in kappa matrices
        if np.all(np.isnan(press[t,:])) or np.all(np.isnan(temp[t,:])):
            kappa1[t,:] = np.nan
            kappa2[t,:] = np.nan
        else:
            pidx = val2idx(lut.press.values, press[t,:])
            tidx = val2idx(lut.temp.values, temp[t,:])
            
            #replace -99999 that val2idx finds with -1 so the assigning works out for finding a value from the LUT
            pidx[pidx < -999] = -1
            tidx[tidx < -999] = -1
            #find LUT value for pidx, tidx; and replace kappa1 kappa2 with invalid press or temp with nan
            kappa1[t,:] = lut.kappa.values[0, tidx, pidx]
            kappa1[t,(pidx==-1)|(tidx==-1)] = np.nan
            
            kappa2[t,:] = lut.kappa.values[1, tidx, pidx]
            kappa2[t,(pidx==-1)|(tidx==-1)] = np.nan
            
    
    #return kappa1, kappa2
    return kappa2-kappa1



def get_diffgamma(R, inradar, outall=False):
    '''
    get differential gamma by looping in the vertical. for each range bin, find the index (resiterated) which is R away starting from index i. 
    note that resiterated-i (number of ranges in R) varies due to different chirp resolution. then, calculate gamma for each i.
    input: 
    - R: gliding window of resolution  in meters
    - radar: xarray dataset with ze1, ze2, dar input on time, height dimensions.
    output:
    -diffgamma: differential gamma calculated by gamma2-gamma1
    
    '''
    
    radar = inradar.copy()
    
    gamma1 = [] #gamma is diff coefficient per each range gate and each time stamp
    gamma2 = []
    nvolumes = [] #number of radar volumes averaged to obtain R spacing
    
    nz = radar.height.shape[0]
    
    
    #check whether ze input is in dBZ; if yes, convert to linear units for dar calculations.
    for var in ['GZe', 'G2Ze']:
        if radar[var].units == 'dBZ':
            radar[var].values = helpfct.get_zlin(radar[var].values)
            radar[var].attrs['units'] = 'mm6 m-3'
    
    
    #convention: i is r1; resiterated is r2
    #now loop through height dimension and find where next radar range is within R:
    for i in range(nz):
        
        resiterated = i
        
        while radar.height.values[resiterated] - radar.height.values[i] < R and resiterated < nz-1:
            resiterated+=1
        
        #count how many ranges are within R distance: (later used for quality flagging):
        nvolumes.append((resiterated-1 - i)-1) #count number of layers, not levels; also is number of bottom levels here though
        
        #calculate gamma coefficient
        gamma1.append(1/(2*R) * np.log( radar.GZe.values[:,i]/radar.GZe.values[:,resiterated-1]))
        gamma2.append(1/(2*R) * np.log( radar.G2Ze.values[:,i]/radar.G2Ze.values[:,resiterated-1]))
        
    
    gamma1 = np.asarray(gamma1).T
    gamma2 = np.asarray(gamma2).T
    
    diffgamma = gamma2 - gamma1
    
    #return swapped time and height dimensions
    return diffgamma, nvolumes

