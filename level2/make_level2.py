import xarray as xr
import numpy as np
import datetime as dt
import json

from level2.src import level2_functions as fct
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
    
    '''13.7. comment attenuation correction just for laptop usage right now; make sure to remove!
    #attenuation correction
    
    #interpolate att onto radar time and height:
    
    attinter = att.Att_atmo.interp(time=level2.time, height=level2.height)
    
    #calculate attenuation corrected variables. (and rename the old ones first)
    
    for level2key, attkey in [['W',94], ['G',167], ['G2',174]]:
        level2['%sZe_attcorr'%level2key] = level2['%sZe'%level2key] + attinter.Att_atmo.sel(freq=int(attkey))
    
    level2 = level2.assign(DFR_attcorr = level2.WZe_attcorr - level2.GZe_attcorr, {'units':'dB', 'long_name':'attenuation corrected DFR W-G'})
    
    '''
    
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









def ja():


    #==== load level 1 file:
    #datapath = '/work/schnitts/gband/hamag/joint_dataset/'
    #datafile = datapath + 'RF06.nc'





    



    #==== load attenuation profiles
    #load attenuation profiles. should all be on arbitrary model height that reaches further than flight altitude.

    #attpath = '/work/schnitts/gband/hamag/joint_dataset/attenuation/'
    attpath = info['paths']['level2']['att']
    #load all attenuation files available for the entire campaign: this needs some re-thinking for day-to-day processing!
    attfiles = sorted(glob.glob(attpath + '%s_*_attenuation.nc'%(info['global']['mission']))) #files with gas attenuation calculated for each sonde profile
    print(attfiles)
    nsondes = len(attfiles)
    print(nsondes)
    #open each attenuation file; copy att_atmo into new 3d array
    for n in range(nsondes):
        attdata = xr.open_dataset(attfiles[n])
        if n == 0:
            nfreq = attdata.freq.shape[0] #attenuation is calculated for three frequencies
            nz = attdata.height.shape[0]  #number of levels in attenuation file
            allatt =  np.ones((nsondes, nz, nfreq))
            rsrh = rstemp = rspress = rsrhov= np.ones((nsondes, nz))
            rsiwv = np.ones(nsondes)
            dstimes = []

        allatt[n, :, :] = attdata.Att_atmo.values
        #rsrh[n,:] = attdata.relh.values*100. #in %
        #rstemp[n,:] = attdata.temp.values #in K
        #rspress[n,:] = attdata.p.values #in hPa
        #rsrhov[n,:] = wv.rh_to_qabs(attdata.relh.values, attdata.temp.values)*1000.
        #rsdp = mpy.dewpoint_from_relative_humidity(attdata.temp.values*units.K, attdata.relh.values*units.percent)
        #rsiwv[n] = mpy.precipitable_water(attdata.p.values*units.hPa, rsdp).magnitude
        dstimes.append(dt.datetime.strptime(attdata.attrs['sonde_launch_time'], '%Y%m%d%H%M'))  #also extract launch time as datetime object

    dstimes = np.array(dstimes)

    #now put all attenuation profiles from the different files into a xr dataset:
    data_vars = {
        'att_atmo':(['time', 'height', 'freq'], allatt), 'rs_rh':(['time', 'height'], rsrh), 'rs_temp':(['time', 'height'], rstemp), 'rs_press':(['time', 'height'], rspress), 'rs_rhov':(['time', 'height'], rsrhov), 'rs_iwv':(['time'], rsiwv)}

    coords = {'time':(['time'], dstimes), 'height':(['height'], attdata.height.values), 'freq':(['freq'], np.array([94, 167, 174]))}
    xratt = xr.Dataset(data_vars, coords)

    #interpolate the xr dataset to gband level1 times; and to radar range gate altitudes; apply slant factor
    radarhgt = l1.height
    attinter = xratt.interp(time=l1timedt, height=radarhgt)

    #now, attinter has dimensions:  ntime, nrange, nfreq

    #calculate 2-way top-down:
    attintertd = 2*attinter.cumsum(dim='height')
    #somehow lost a dimension range here; re-assign:
    attintertd = attintertd.assign_coords({'height':attinter.height})

    #============= apply attenuation correction to profiles
    #Ze_corrected = Ze_measured + attenuation.  thinking of Ze in dBZ

    W_Ze_attcorr = WZe + attintertd.att_atmo.sel(freq=94)
    G_Ze_attcorr = GZe + attintertd.att_atmo.sel(freq=167)
    G_Ze2_attcorr = GZe2 + attintertd.att_atmo.sel(freq=174)

    DFR_attcorr = W_Ze_attcorr - G_Ze_attcorr

    #==== clutter and quality masks
    print('...TODO: clutter and quality masks including flagging for rotation etc...')



    #==== create dataset; water vapor output will be added later and all will be stored later


    #outputfile = '/work/schnitts/gband/hamag/joint_dataset/HAMAG_draco_lv2_RF06.nc'
    outputfile = info['paths']['level2']['output']

    l2 = make_dataset(l1timedt, l1.height, WZe, GZe, GZe2, Wvd, Gvd, DDV, DFR, DAR, W_Ze_attcorr, G_Ze_attcorr, G_Ze2_attcorr, DFR_attcorr, attintertd,  info['attrs']['level2'], outputfile, write = False)
    return
