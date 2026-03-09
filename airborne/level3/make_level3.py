import xarray as xr
import numpy as np
import datetime as dt
import json

from level3.src import level3_functions as fct
from src import helpfunctions as helpfct

def run(level2, att, info, write=False):
    '''
    create level3a dataset, including attenuation correction.
    input:
    - level2: level2 xarray dataset (dual-frequency matched)
    - att: xarray dataset with attenuation information
    - info: xarray dataset with attribtues
    - write: optional, if True: output stored to netcdf
    output:
    - level3a: attenuation corrected.
    '''
    
    level3a = level2.copy()
    
    #convert Ze and DFR to dB and change unit attribute in xarray:
    for zvars in ['WZe', 'GZe', 'G2Ze']:
        level3a[zvars].values = helpfct.get_ZdBZ(level2[zvars].values)
        level3a[zvars].attrs['units'] = 'dBZ'
    
    #level3a['DAR'].values = helpfct.get_ZdBZ(level2['DAR'].values) #corrected on 3/3/26: DAR is in dB already.
    #level3a['DAR'] = level2['DAR'].values
    #level3a['DAR'].attrs['units'] = 'dB'
    
    #also convert DFR:
    level3a.DFR.values = level2.WZe.values - level2.GZe.values
    level3a.DFR.attrs['units'] = 'dB'
    
    #convert time (level1 is in seconds since 1/1/1970) to datetime timestamp:
    nt = level3a.time.shape[0]
    level3a = level3a.assign_coords(time = np.array([dt.datetime.utcfromtimestamp(level3a.time.values[t]) for t in range(nt)], dtype=object))
    
    #interpolate att onto radar time and height:
    print('interpolating attenuation to radar time and height...')
    attinter = att.Att_atmo.interp(time=level3a.time, height=level3a.height)
    
    #calculate attenuation corrected variables. (and rename the old ones first)
    
    for level3akey, attkey in [['W',94], ['G',167], ['G2',174]]:
        level3a['%sZe_attcorr'%level3akey] = level3a['%sZe'%level3akey] + attinter.sel(freq=int(attkey), method='nearest')
    
    level3a = level3a.assign(DFR_attcorr = level3a.WZe_attcorr - level3a.GZe_attcorr)
    level3a.DFR_attcorr.attrs['units'] = 'dB'
    level3a.DFR_attcorr.attrs['long_name'] = 'attenuation corrected DFR W-G'
    
    #change some attributes:
    level3a.attrs['title'] = 'Joint dataset of Grawac and Wband, including attenuation correction (Level3a)'
    level3a.attrs['creation_time'] = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
    level3a.attrs['product-id'] = 'Level-3a'
    
    if write==True:
        ncfile = info['global']['mission'] +'_%s_level3a_draco_%s%s%s.nc'%(info['global']['version'], info['yyyy'], info['mm'], info['dd'])
        level3a.to_netcdf(info['paths']['output'] + ncfile)
    
    return level3a




def attenuation_correction(attfiles, geometry, slant, info, write=False):
    '''
    gas attenuation function. provides 
    attfiles: files with calculated gas attenuation
    geometry: valid option: td (top-down); bu (bottom-up)
    slant: slant factor; can be set to False if nadir/zenith looking
    
    info: info json dataset
    write: optional, if set True: output stored to netcdf
    '''
    
    #open attenuation profiles
    attrs = xr.open_mfdataset(attfiles) #dataset with attenuation profiles for each sounding+freq, on arb height
    
    #cumulate depending on geometry:
    attcum = attrs.copy()
    
    #convert launchtime coordinate to datetime:
    nt = attcum.launchtime.shape[0]
    attcum = attcum.assign_coords( launchtime=np.asarray([dt.datetime.strptime(str(attrs.launchtime.values[t]), '%Y%m%d%H%M') for t in range(nt)],dtype=object))
    
    #and rename the variable to time:
    attcum = attcum.rename({'launchtime':'time'})
    
    if geometry == 'bu':
        attcum.Att_atmo.values = 2*attrs.Att_atmo.cumsum(dim='height').values
    elif geometry == 'td':
        print('top-down attenuation correction...')
        #attcum.Att_atmo.values = 2*attrs.Att_atmo.cumsum(dim='height').values
        
        slantfactor = np.cos(np.radians(l1.incangle))  #define it as a factor, thus slantpath = radarhgt/cos(angle)
        radarhgt = l1.alt - l1.range*slantfactor
        attinter = xratt.interp(time=l1timedt, height=radarhgt) * 1/slantfactor

        1/0
    else:
        print('geometry not specified. returning.')
        return
    
    if slant != False:
        print('applying slant factor.')
    
    if write == True:
        attcum.to_netcdf(info['paths']['output']+'attenuation.nc')
    
    return attcum
