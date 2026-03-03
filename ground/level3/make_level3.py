import xarray as xr
import numpy as np
import datetime as dt
import json

from level3.src import level3_functions as fct
from src import helpfunctions as helpfct

def run(level1, att, info, write=False):
    '''
    create level2 dataset, including attenuation correction.
    input:
    - level1: level1 xarray dataset
    - att: xarray dataset with attenuation information
    - info: xarray dataset with attribtues
    - write: optional, if True: output stored to netcdf:
    '''
    
    level2 = level1.copy()
    
    #convert Ze and DFR to dB and change unit attribute in xarray:
    for zvars in ['WZe', 'GZe', 'G2Ze']:
        level2[zvars].values = helpfct.get_ZdBZ(level1[zvars].values)
        level2[zvars].attrs['units'] = 'dBZ'
    
    level2['DAR'].values = helpfct.get_ZdBZ(level1['DAR'].values)
    level2['DAR'].attrs['units'] = 'dB'
    
    #also convert DFR:
    level2.DFR.values = level2.WZe.values - level2.GZe.values
    level2.DFR.attrs['units'] = 'dB'
    
    #convert time (level1 is in seconds since 1/1/1970) to datetime timestamp:
    nt = level2.time.shape[0]
    level2 = level2.assign_coords(time = np.array([dt.datetime.utcfromtimestamp(level2.time.values[t]) for t in range(nt)], dtype=object))
    
    #quality filtering
    print('quality filtering yet to be done. ie: clutter filter; above SNR for DAR signal')
    
    #13.7. comment attenuation correction just for laptop usage right now; make sure to remove!
    #attenuation correction
    
    #interpolate att onto radar time and height:
    
    attinter = att.Att_atmo.interp(time=level2.time, height=level2.height)
    
    #calculate attenuation corrected variables. (and rename the old ones first)
    
    for level2key, attkey in [['W',94], ['G',167], ['G2',174]]:
        level2['%sZe_attcorr'%level2key] = level2['%sZe'%level2key] + attinter.sel(freq=int(attkey), method='nearest')
    
    level2 = level2.assign(DFR_attcorr = level2.WZe_attcorr - level2.GZe_attcorr)
    level2.DFR_attcorr.attrs['units'] = 'dB'
    level2.DFR_attcorr.attrs['long_name'] = 'attenuation corrected DFR W-G'
    
    
    
    #change some attributes:
    level2.attrs['title'] = 'Joint dataset of Grawac and Wband (Level2)'
    level2.attrs['creation_time'] = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
    level2.attrs['product-id'] = 'Level-2'
    
    if write==True:
        ncfile = info['global']['mission'] +'_%s_level2_draco_%s%s%s.nc'%(info['global']['version'], info['yyyy'], info['mm'], info['dd'])
        level2.to_netcdf(info['paths']['output'] + ncfile)
    
    return level2




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
        print('check calculation.')
        attcum.Att_atmo.values = 2*attrs.Att_atmo.cumsum(dim='height').values
        1/0
    else:
        print('geometry not specified. returning.')
        return
    
    if slant != False:
        print('applying slant factor.')
    
    if write == True:
        attcum.to_netcdf('./attenuation.nc')
    
    return attcum
