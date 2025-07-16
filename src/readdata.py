import xarray as xr
import numpy as np
import datetime as dt
from src import helpfunctions as srcfct

def make_dataset(dsvars, attrs, outfile, withdar =True, write=False):
    
    #create xarray dataset with input data: ===============================
    #dimensions: time, height, nchirps (for the indices)
    
    coords = {'time':(['time'], dsvars['time'], {'units':'UTC', 'long_name':'Timestamp logged at end of chirp sequence of GRaWAC'}),
              'height':(['height'], dsvars['rheight'], {'units':'m', 'long_name':'Radar range gate height asl (radar range + instrument_altitude'}),
              'nchirp':(['nchirp'], np.arange(dsvars['nchirps']))
             }
    
    #datavars to include: ze, ldr, startidx, freq, 
    data_vars = { 
    'Ze' : (['time','height'], dsvars['Ze'], {'units': 'mm6 m-3', 'long_name':'radar reflectivity Ze'}),
    'vd' : (['time','height'], dsvars['vd'], {'units': 'm s-1', 'long_name':'radar mean Doppler Velocity'}),
    'chirpseq_startix' : (['nchirp'], dsvars['chirpstartidx'], {'units': '-', 'long_name':'starting index in range array for where each chirp starts'}),
    'navg' : (['nchirp'], dsvars['navg'], {'units': '-', 'long_name':'number of chirps averaged.'}),
    'maxvel' : (['nchirp'], dsvars['maxvel'], {'units': 'm/s', 'long_name':'Nyquist velocity each chirp sequence.'}),
    'instrument_altitude': (['nchirp'], np.broadcast_to(dsvars['instrumentalt'], dsvars['nchirps']), {'units': 'm', 'long_name':'altitude msl of sensor. variable was added to range to obtain height variable. '}),
    'freq':(['nchirp'], np.broadcast_to(dsvars['freq'], dsvars['nchirps']), {'units':'GHz', 'long_name':'Radar transmitted frequency' })
        }
    
    if withdar == True:
        data_vars['DAR'] = (['time','height'], dsvars['LDR'], {'units': 'dB', 'long_name':'Differential Absorption Radar Dual-Frequency Ratio'})
        data_vars['Ze2'] = (['time','height'], dsvars['Ze2'], {'units': 'mm6 m-3', 'long_name':'GRaWAC Ze at 174.7GHz calculated from Ze167 and LDR'})
    
    
    attrs['creation_time'] = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
    
    
    xrds = xr.Dataset(data_vars, coords, attrs)
    if write==True:
        xrds.to_netcdf(outfile)
    
    return xrds


def read_rpg_lv1(infiles, info, instrumentname, instrumentaltinput=True, withdar=True, write=False):
    '''
    read specified filename(s) of lv1 data set; combine in xr dataset and return. 
    if withdar == True: also read LDR variable.
    '''
    
    #load RPG files:
    try:
        indata = xr.open_mfdataset(infiles)
    except ValueError:
        indata = []
        for infile in infiles:
            infilesingle = xr.open_dataset(infile)
            indata.append(infilesingle)
            infilesingle.close()
            
        indata = xr.concat(indata,dim='Time')
        print('need to remove double time stamps???')
    
    nchirps = indata.Chirp.shape[0]
    
    # edit data to prep for output ==================================
    
    #calculate real times:
    #get timestamp with millisecond precision:
    timesec = indata.Time.values + 0.001*indata.Timems.values
    #calculate total seconds between 1/1/2001 and 1/1/1970
    rpgoffset = (dt.datetime(2001,1,1) - dt.datetime(1970,1,1)).total_seconds()
    #calculate grawactime in seconds since 1/1/1970
    time = timesec+rpgoffset

    
    #concatenate along height dimensions
    Ze = np.concatenate([indata['C%iZE'%i].values for i in np.arange(1,nchirps+1)], axis=1)
    vd = np.concatenate([indata['C%iMeanVel'%i].values for i in np.arange(1,nchirps+1)], axis=1)
    rrange = np.concatenate([indata['C%iRange'%i].values for i in np.arange(1,nchirps+1)])
    
    #set -999 to nan:
    Ze[Ze == -999.] = np.nan
    vd[vd == -999.] = np.nan
    
    if withdar == True:
        #get LDR and calculate Ze2:
        LDR = np.concatenate([indata['C%iLDR'%i].values for i in np.arange(1,nchirps+1)], axis=1)
        #LDR missing values are -100 and -999; set them all to -999.
        LDR[(LDR == -100.) | (LDR == -999)] = np.nan

        #added factor -1 as RPG logs DAR the wrong way around (174-167 instead of 167-174). LDR is in dB
        LDR *= -1.
        #now set missing values back to -999:
        #LDR[LDR>900] = -999.
        
        #calculate Ze 175 channel from Ze@167 and LDR in linear units mm6/m3
        LDRlin = 10**(LDR/10.)
        Ze2 = Ze / LDRlin

    
    #add instrument altitude to range to obtain height as dimension
    if 'instrument_altitude' in list(indata.keys()):
        instrumentalt = np.unique(indata.instrument_altitude.values)[0]
    else:
        instrumentalt = instrumentaltinput
        
    rheight = rrange + instrumentalt
    
    #add index where each chirp sequence starts in overall concatenated height array:
    chirpstartidx = [0]
    for n in range(nchirps-1):
        chirpstartidx.append(chirpstartidx[n] + len(indata['C%iRange'%(n+1)].values))
    
    #frequency:
    freq = np.unique(indata.Freq.values)[0] #returns a number
    
    #average number:
    navg = np.unique(indata.AvgNum.values, axis=0)[0] #returns array with nchirps entries
    
    #max velocity:
    maxvel = np.unique(indata.MaxVel.values, axis=0)[0] #returns array with nchirps entries
    
    #create dataset:
    exportvars = [Ze, vd, rheight, freq, navg, maxvel, time, nchirps, instrumentalt, chirpstartidx]
    varnames = ['Ze', 'vd', 'rheight', 'freq', 'navg', 'maxvel', 'time', 'nchirps', 'instrumentalt', 'chirpstartidx']
    if withdar == True:
        exportvars.append(LDR)
        varnames.append('LDR')
        exportvars.append(Ze2)
        varnames.append('Ze2')
    
    tods = {}
    for vv,var in enumerate(exportvars):
        tods[varnames[vv]] = var
    
    attrs = {'title':'radar data read from RPG LV1.NC output'}
    outfile = info['paths']['output'] + info['global']['mission'] +'_%s_level0_%s_%s%s%s.nc'%(info['global']['version'], instrumentname, info['yyyy'], info['mm'], info['dd'])
    
    xrds = make_dataset(tods, attrs, outfile, withdar=withdar, write=write)
    
    indata.close()
    
    return xrds





def read_compactfiles(indatafiles, inhkfiles, info, instrumentname, withdar = False, write=False ):
    '''
    read compact nc files, return same dataset as above
    '''
    
    try:
        indata = xr.open_mfdataset(indatafiles)
    except ValueError:
        indata = []
        for infile in indatafiles:
            indatasingle = xr.open_dataset(infile)
            indata.append(indatasingle)
            indatasingle.close()
        indata = xr.concat(indata,dim='Time')
        print('need to remove double time stamps???')    
    
    #also housekeeping files:
    inhkdata = xr.open_mfdataset(inhkfiles)
    
    nchirps = len(indata.chirp_sequence.values)
    
    #convert wband time from numpy64datetime to seconds since 1/1/1970:
    unix_epoch = np.datetime64(0,"s")
    one_second = np.timedelta64(1,"s")
    time = (indata.time.values - unix_epoch) / one_second
    print('TODO:check in compact processing whether the time dimension already includes the microseconds stored in the sampletms variable')
    
    #get variables; and set -999 to nan:
    Ze = srcfct.get_zlin(indata.ze.values)
    Ze[Ze == -999.] = np.nan
    vd = indata.vm.values
    vd[vd == -999.] = np.nan
    rheight = indata.height.values #in compact files, stored height coordinate already includes instrumentaltitude
    freq = np.unique(indata.frequency.values)[0]
    
    navg = np.unique(inhkdata.num_avg_chirps.values, axis=0)[0]
    #sampletms = np.unique(inhkdata.sample_tms.values, axis=0)[0]
    maxvel = np.unique(indata.nyquist_velocity.values, axis=0)[0] #returns array with nchirps entries
    
    instrumentalt = np.unique(indata.instrument_altitude)[0]
    chirpstartidx = np.unique(indata.chirpseq_startix.values, axis=0)[0]
    
    if withdar == True:
        #get LDR and calculate Ze2:
        LDR = np.concatenate([indata['C%iLDR'%i].values for i in np.arange(1,nchirps+1)], axis=1)
        #LDR missing values are -100 and -999; set them all to -999.
        LDR[(LDR == -100.) | (LDR == -999)] = np.nan

        #added factor -1 as RPG logs DAR the wrong way around (174-167 instead of 167-174). LDR is in dB
        LDR *= -1.
        #now set missing values back to -999:
        #LDR[LDR>900] = -999.
        
        #calculate Ze 175 channel from Ze@167 and LDR in linear units mm6/m3
        LDRlin = 10**(LDR/10.)
        Ze2 = Ze / LDRlin
    
    
    #create dataset:
    exportvars = [Ze, vd, rheight, freq, navg, maxvel, time, nchirps, instrumentalt, chirpstartidx]
    varnames = ['Ze', 'vd', 'rheight', 'freq', 'navg', 'maxvel', 'time', 'nchirps', 'instrumentalt', 'chirpstartidx']
    if withdar == True:
        exportvars.append(LDR)
        varnames.append('LDR')
        exportvars.append(Ze2)
        varnames.append('Ze2')
    
    tods = {}
    for vv,var in enumerate(exportvars):
        tods[varnames[vv]] = var
    
    attrs = {'title':'radar data read from compact files output produced by matlab processing'}
    outfile = info['paths']['output'] + info['global']['mission'] +'_%s_level0_%s_%s%s%s.nc'%(info['global']['version'], instrumentname, info['yyyy'], info['mm'], info['dd'])
    
    xrds = make_dataset(tods, attrs, outfile, withdar=withdar, write=write )
    
    indata.close()
    inhkdata.close()

    return xrds
