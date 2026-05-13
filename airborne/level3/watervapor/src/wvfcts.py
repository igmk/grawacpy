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
              'range':(['range'], dsvars['range'], {'units':'m', 'long_name':'Radar range distance from aircraft'}),
              'R':(['R'], [dsvars['R']], {'units':'m', 'long_name':'retrieval vertical resolution'}),
              'tavg':(['tavg'], [dsvars['tavg']], {'units':'s', 'long_name':'retrieval averaging time'})
             }
    
    #datavars to include: ze, ldr, startidx, freq, 
    data_vars = { 
    'rhov' : (['time','range'], dsvars['rhov'], {'units': 'g m-3', 'long_name':'absolute humidity profile retrieved from DAR measurements along slanted path'}),
    'diffkappa' : (['time','range'], dsvars['diffkappa'], {'units': 'mm6 m-3', 'long_name':'radar reflectivity Ze'}),
    'diffgamma' : (['time','range'], dsvars['diffgamma'], {'units': 'm s-1', 'long_name':'radar mean Doppler Velocity'}),
    'deltarhov' : (['time','range'], dsvars['deltarhov'], {'units': 'g m-3', 'long_name':'rhov uncertainty'}),
    'nranges' : (['range'], dsvars['nranges'], {'units': '-', 'long_name':'number of range bins that fit into R spacing'}),
    'height' : (['time','range'], dsvars['height'], {'units': 'm asl', 'long_name':'Altitude asl for each radar range bin.'}),
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
    deltaZerR = [] #variance term (second term) for rhov uncertainty calculation
    nz = radar.range.shape[0]
    
    
    #check whether ze input is in dBZ; if yes, convert to linear units for dar calculations.
    for var in ['GZe', 'G2Ze']:
        if radar[var].units == 'dBZ':
            radar[var].values = helpfct.get_zlin(radar[var].values)
            radar[var].attrs['units'] = 'mm6 m-3'
    
    
    #convention: i is r1; resiterated is r2
    #now loop through height dimension and find where next radar range is within R:
    for i in range(nz):
        
        resiterated = i
        
        while radar.range.values[resiterated] - radar.range.values[i] < R and resiterated < nz-1:
            resiterated+=1
        
        #count how many ranges are within R distance: (later used for quality flagging):
        nvolumes.append((resiterated-1 - i)-1) #count number of layers, not levels; also is number of bottom levels here though
        
        #calculate gamma coefficient
        gamma1.append(1/(2*R) * np.log( radar.GZe.values[:,i]/radar.GZe.values[:,resiterated-1]))
        gamma2.append(1/(2*R) * np.log( radar.G2Ze.values[:,i]/radar.G2Ze.values[:,resiterated-1]))
        
        #this one is the second term of Roy et al 2018 Eq (13)
        #convert deltaZe (in dB) to linear mm6/m3 by dividing with 4.343 factor
        sigmarhovterm2 =  np.sqrt( 
            ((radar.deltaZe.values[:,i]/4.343) / radar.GZe.values[:,i] ) **2 + #deltaZe at r1 for 167
            ((radar.deltaZe.values[:,i]/4.343) / radar.G2Ze.values[:,i] ) **2 + #assume same deltaZe for 174 at r1
            ((radar.deltaZe.values[:,resiterated-1]/4.343) / radar.GZe.values[:,resiterated-1] ) **2 + #deltaZe at r2 for 167
            ((radar.deltaZe.values[:,resiterated-1]/4.343) / radar.G2Ze.values[:,resiterated-1] ) **2 ) #assume same deltaZe for 174 at r2

        deltaZerR.append(sigmarhovterm2)
    
    gamma1 = np.asarray(gamma1).T
    gamma2 = np.asarray(gamma2).T
    deltaZerR = np.asarray(deltaZerR).T

    diffgamma = gamma2 - gamma1
    
    #return swapped time and height dimensions
    return diffgamma, nvolumes, deltaZerR


def calc_sigma0(Zmax, freq, theta, dr, Kw, Fbf):
    '''
    calculate normalized calibrated radar cross section sigma0 of ground return
    input: 
    - zmax: logarithmic Ze used to derive sigma0 from...could be maximum; or mean of all ground bins? [dBZe mm6/m3]
    - freq: frequency in GHz of Ze measurement [GHz]
    - theta: inclination angle off nadir [deg]
    - dr: range resolution of ground return
    - Kw: dielectric constant
    - Fbf: 
    '''
    #convert Ze to linear units, ie mm6/m3
    zlin = helpfct.get_zlin(Zmax)
    
    #calculate inclination angle in radians:
    thetarad = np.radians(theta)
    
    #calculate wavelength: [m]
    c = 3e08
    wvl = c/(1e09*freq)
    print('radar wvl: %.5f m'%wvl)
    
    #factor 1e-18 is there to calculate Zemax in m6/m3 (comes in mm6/m3)
    sigma0 = 1e-18 * zlin * Kw**2 * np.pi**5 * np.cos(thetarad) * dr / (Fbf * wvl**4) 
    
    sigma0dB = helpfct.get_ZdBZ(sigma0)
    
    return sigma0dB



def retrieve_iwv(Ze1, Ze2, freq1, freq2, zres, incangle, radaralt, radarrange, xrds, lut, kwargs, wvprof):
    '''
    retrieve column amount between radar and specific range (ie ground return (air); (or cloud top (air) ?); or cloud base (ground-based) )
    
    version 2.0:
    - changed compared to v1.0 above as diffkappa is available from wvprof input
    
    input:
    - Ze1
    - Ze2
    - freq1
    - freq2
    - zres
    - incangle
    - radaralt
    - radarrange
    - xrds
    - lut
    - kwargs
    - wvprof
    
    '''
    
    if kwargs['groundsigma'] == 'max':
        Zein1 = np.nanmax( Ze1.where((radaralt > kwargs['minalt']) & (radaralt < kwargs['maxalt'])),axis=1)
        Zein2 = np.nanmax( Ze2.where((radaralt > kwargs['minalt']) & (radaralt < kwargs['maxalt'])),axis=1)
    
    
    Kw = 0.86   #frequency and surface dependent!!!
    Fbf = 1
    
    #get sigma0 in dB for each channel: 
    sigma0f1 = calc_sigma0(Zein1, freq1, incangle, zres, Kw, Fbf)
    sigma0f2 = calc_sigma0(Zein2, freq2, incangle, zres, Kw, Fbf)
    
    #calculate iwv
    #diffkappa = 0.07 #v1.0
    #assume diffkappa as a constant factor. (31.10.25). 
    # #checked in differential kappa from LUT that temporal STD of diffkappa throughout RF05 is around 1.5 - 4.0 *10^-4. 
    # Roy et al 2021 Eqn (9) and text below state that Delta Kappa can be taken as constant (independent of humidity profile) if kappa doesn't vary much between aircraft and ground.
    # at time of dropsonde launch 11:16, delta kappa @100m is 0.1075; @ 3500m: 0.101
    # thus, just assume diffkappa as 0.1 for RF05 for now.
    
    #diffkappa = 0.1 #v2.0; 

    diffkappa = np.nanmean(wvprof.diffkappa.values)

    print('found mean diffkappa for flight = %.2f'%diffkappa)
    diffkappastd = np.nanstd(wvprof.diffkappa.values, axis=0) #temporal standard deviation throughout flight
    print('max temporal STD of diffkappa throughout entire flight = %.5f'%np.nanmax(diffkappastd))
    diffkappameanvertical = np.nanmean(wvprof.diffkappa.values[:,-25] - wvprof.diffkappa.values[:,25] ) #mean vertical difference between aircraft and ground level of diffkappa
    diffkappastdvertical = np.nanstd(wvprof.diffkappa.values[:,-25] - wvprof.diffkappa.values[:,25] ) #vertical STD
    print('mean vertical diffkappa difference = %.5f, and its STD = %.5f'%(diffkappameanvertical, diffkappastdvertical))
    
    #iwv = calc_iwv(incangle, diffkappa, sigma0f1, sigma0f2)
    
    thetarad = np.radians(incangle)
    sigmaf1lin = helpfct.get_zlin(sigma0f1)
    sigmaf2lin = helpfct.get_zlin(sigma0f2)
    
    iwv = np.cos(thetarad)/(2*diffkappa) * np.log(sigmaf1lin/sigmaf2lin)
    
    return iwv, diffkappa
    



    return iwv