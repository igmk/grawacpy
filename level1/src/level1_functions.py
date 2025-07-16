import numpy as np
import datetime as dt
import xarray as xr
from scipy.spatial import cKDTree

def calculate_chirp_repetition_time(nyqvel, freq):
    '''
    calculate chirp repetition frequency and time given nyquist velocity (ie max velocity set in each chirp sequence) and frequency
    input:
    - nyqvel: Nyquist velocity: max velocity limiting spectrum [m/s]
    - freq: radar frequency [Hz]
    output:
    - chirpreptime: chirp repetition time per chirp [s]
    '''
    c = 299702547 #m/s speed of light
    nyqfreq = nyqvel * 2 * freq / c
    chirpreptime = 1/(2*nyqfreq)
    
    return chirpreptime





def get_seqtime(intime, chirpduration, nchirps):
    '''
    calculate time stamp for each chirp sequence; ie timestamp in the middle of each chirp. assumes that logged timestamp is at the end of all chirp sequences; and that chirpduration is how long each chirp takes
    input: 
    - intime: grawac timestamp in float precision (Time + timems*0.001) [s]
    - chirpduration: nt x nchirps, chirp duration time for each chirp [s]
    output:
    - seqtime: array, dimension nt x nchirp; this time corresponds to timestamp in middle of each chirp sequence
    '''
    nt = len(intime)
    seqtime = np.ones((nt, nchirps))
    
    if nchirps == 2:
        seqtime[:,1] = np.mean([intime, intime - chirpduration[1]], axis=0)
        seqtime[:,0] = np.mean([intime - chirpduration[1], intime - np.sum(chirpduration[:])], axis=0)
    
    elif nchirps == 3:
        seqtime[:,2] = np.mean([intime, intime-chirpduration[-1]], axis=0)  #thats the last one before total timestamp is logged
        seqtime[:,1] = np.mean([intime-chirpduration[-1], intime-(np.sum(chirpduration[1:-1]))], axis=0) #thats the middle one
        seqtime[:,0] = np.mean([intime-(np.sum(chirpduration[1:-1])), intime-np.sum(chirpduration)], axis=0) #first one
    
    elif nchirps == 4:
        seqtime[:,3] = np.mean([intime, intime-chirpduration[-1]], axis=0)  #thats the last one before total timestamp is logged
        seqtime[:,2] = np.mean([intime-chirpduration[-1], intime-(np.sum(chirpduration[2:-1]))], axis=0) #thats the middle one
        seqtime[:,1] = np.mean([intime-np.sum(chirpduration[2:-1]), intime-(np.sum(chirpduration[1:-1]))], axis=0) #thats the middle one
        seqtime[:,0] = np.mean([intime-(np.sum(chirpduration[1:-1])), intime-np.sum(chirpduration)], axis=0) #first one
    
    return seqtime
        



#============= functions


def range_matching(wdata, gdata, modes):
    '''
    matches W and Gband in range
    input:
    - wdata level0
    - gdata level0
    - matchmode: array with nchirp entries specifying 'average' (average all bins within one G-band bin) or 'nearest' (find nearest neighbor in terms of center range)
    '''
    
    '''
    #loop through chirp sequences:
    for n in range(nchirps):
        print('chirp %i'%n)
        print(gchirpidxstart, gchirpidxend)
        wchirpidxstart = int(np.unique(wdata.chirpseq_startix.values[:,n])[0])
        gchirpidxstart = gchirpidxstart
        gchirpidxend = gchirpidxstart + len(gdata['C%iRange'%(n+1)].values)
        print(gchirpidxstart, gchirpidxend)
        
        if n < nchirps-1: #for all but last chirp limit the range with chirpidxend
            wchirpidxend = int(np.unique(wdata.chirpseq_startix.values[:,n+1])[0])
        
            W_Ze_rc_std[:,gchirpidxstart:gchirpidxend],W_Ze_rc[:,gchirpidxstart:gchirpidxend], W_vd_rc_std[:,gchirpidxstart:gchirpidxend],         W_vd_rc[:,gchirpidxstart:gchirpidxend] = get_closest_range_W_G(wdata['height'].values[wchirpidxstart:wchirpidxend], gdata['C%iRange'%(n+1)].values + instrumentaltitude,        wzelin[:,wchirpidxstart:wchirpidxend], wdata['vm'].values[:,wchirpidxstart:wchirpidxend], 
            mode=mode[n])
        
        else:
            #hinten ende offen
            W_Ze_rc_std[:,gchirpidxstart:], W_Ze_rc[:,gchirpidxstart:],         W_vd_rc_std[:,gchirpidxstart:], W_vd_rc[:,gchirpidxstart:] =  get_closest_range_W_G( wdata['height'].values[wchirpidxstart:], gdata['C%iRange'%(n+1)].values + instrumentaltitude, wzelin[:,wchirpidxstart:], wdata['vm'].values[:,wchirpidxstart:], mode=mode[n])
        
        #at the end of each loop, store the end digit as new start index in gchirpidx
        gchirpidxstart = gchirpidxend

    print('both W and G Ze have matched range resolution.')
    '''
    
    
    nchirps = wdata.nchirp.shape[0]
    #get start and end indices in this way: [wstartidx:wendidx] (ie end index is not called)
    wstartidx = [int(x) for x in wdata.chirpseq_startix.values]
    gstartidx = [int(x) for x in gdata.chirpseq_startix.values]
    
    wendidx = [int(wdata.chirpseq_startix.values[i+1]) for i in range(nchirps-1)] 
    wendidx.append((wdata.Ze.shape[1]-1))
    
    gendidx = [int(gdata.chirpseq_startix.values[i+1]) for i in range(nchirps-1)]
    gendidx.append((gdata.Ze.shape[1]-1))
    
    # loop through chirp settings to give option to match differently per chirp (more options depending on all different kinds of chirp settings...); pass only each chirp sequence into the matching function and its respective mode
    dicsmatch = []
    for n in range(nchirps):
        print('range matching chirp %i in mode %s'%(n+1, modes[n]))
        if n < nchirps-1:
            dicsmatch.append(do_range_matching(wdata.height[wstartidx[n]:wendidx[n]], gdata.height[gstartidx[n]:gendidx[n]], wdata.Ze[:, wstartidx[n]:wendidx[n]], wdata.vd[:, wstartidx[n]:wendidx[n]], modes[n]))
        else:
            dicsmatch.append(do_range_matching(wdata.height[wstartidx[n]:], gdata.height[gstartidx[n]:], wdata.Ze[:, wstartidx[n]:], wdata.vd[:, wstartidx[n]:], modes[n]))
    
    #concatenate all chirps:
    WZematched = np.concatenate([dicsmatch[i]['WZe'] for i in np.arange(nchirps)], axis=1)
    Wvdmatched = np.concatenate([dicsmatch[i]['Wvd'] for i in np.arange(nchirps)], axis=1)
    
    
    #output: xarray dataset with coordinates Wtime, Gtime, l1height (means here G-height)
    
    coords = {'gtime':(['gtime'], gdata.time.values, {'units':'UTC', 'long_name':'Grawac time'}),
              'wtime':(['wtime'], wdata.time.values, {'units':'UTC', 'long_name':'Wband time'}),
              'l1height':(['l1height'], gdata.height.values, {'units':'m', 'long_name':'Level-1 height, same as level0 gband height'}),
              'modes':(['modes'], modes, {'long_name':'mode chosen for each chirp sequence for range matching'})
             }
    
    data_vars = { 
    'WZe' : (['wtime','l1height'], WZematched, {'units': 'mm6 m-3', 'long_name':'Wband radar reflectivity Ze matched to Gband range'}),
    'GZe' : (['gtime','l1height'], gdata.Ze.values, {'units': 'mm6 m-3', 'long_name':'Gband Ze'}),
    'Wvd' : (['wtime','l1height'], Wvdmatched, {'units': 'ms-1', 'long_name':'Wband mean Doppler velocity vd matched to Gband range'}),
    'Gvd' : (['gtime','l1height'], gdata.vd.values, {'units': 'ms-1', 'long_name':'Gband vd'}),
    }
    
    
    attrs = {}
    attrs['title'] = 'Radar moments range matched'
    attrs['creation_time'] = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
        
    l1xrds = xr.Dataset(data_vars, coords, attrs)
    
    return l1xrds




def do_range_matching(rangew, rangeg, WZe, Wvd, mode):
    '''
    calculate Wband Ze and vd matched to G-band range resolution. Two modes are available. "nearest" finds for each G-band range bin the closest W-band range bin (vertical resolution should be similar). "average" finds all W-band bins lying within one G-band bin and calculates the mean Wband Ze/vd for each G.band range bin.
    input:
    - rangew: range array of wband [m] [nzw]
    - rangeg: range array of grawac [m] [nzg]
    - Ze: reflectivity array of W-band [nt, nzg]
    - vd: mean doppler vel array of W-.band [nt, nzg]
    
    - wdata: xarray level0 wdata
    - gdata: xarray level0 gdata
    
    output:
    - Zew_rc_std: if mode==nearest: vertical distance between range center of g-band and nearest neighbor in wband; if mode == average: standard deviation of Zew bins averaged for each g-band range bin. dimensions are on wband time as range-match is before time match.
    - Zew_rc: wband reflectivity range-matched to g-band vertical resolution
    - vdw_rc_std: if mode==nearest: all nan (not applicable); if mode==average: standard deviation of vd bins averaged for each g-band range bin. dimensions are on wband time as range-match is before time match.
    - vdw_rc: wband vdmean range-matched to g-band vertical resolution.
    '''
    #dimensions:
    nzg = rangeg.shape[0]
    nt = WZe.shape[0]
    
    #calculate average vertical range resolution of this chirp
    wres = np.mean(np.diff(rangew.values))
    gres = np.mean(np.diff(rangeg.values))
    print('W-band resolution: %.2f, G-band resolution: %.2f'%(wres, gres))
    
    #define outinds: here is the number of wbins in each g-bin
    outinds = np.ones(nzg, dtype=int)
    
    #make empty arrays for new W-band Ze and vd that are on G-band resolution:
    Zew_rc = np.ones((nt, nzg)) #contains best estimate W-band Ze range-matched to G-band; in case of nearest (comparable range resolution): closest neighbor. in case of average (W-band resolution a lot finer than g-band), it is the average of all W-band bins within each G-band bin
    Zew_rc_std = np.ones((nt, nzg))  #in case of average: STD of all W-band bins within each Gband bin; in case of nearest: distance between nearest neighboring bins
    
    #for vd:
    vdw_rc = np.ones((nt, nzg))
    vdw_rc_std = np.ones((nt, nzg))
    
    #start calculating W-band range-matched depending on mode chosen: (good for when resolution is similar)
    #print('mode chosen: %s'%mode)
    #mode nearest: finds nearest neighboring range bin in W-band for each G-band bin
    if mode == 'nearest':
        
        diff = np.ones(nzg)
        #loop through g-band range array; for each find the closest and store difference
        for i in range(nzg):
            diffrange = np.abs(rangeg[i].values - rangew.values)
            #print(diffrange)
            outinds[i] = np.argmin(diffrange)
            diff[i] = rangeg[i] - rangew[outinds[i]]
        
        #get Ze of nearest neighboring range in W
        Zew_rc = WZe.values[:,outinds]
        vdw_rc = Wvd.values[:,outinds]
        Zew_rc_std = np.broadcast_to(diff,(nt, nzg))
        vdw_rc_std[:] = np.nan
    
    #if mode 'average' is chosen: that is all finer resolved w-band bins in one g-band bin are averaged together to one bin. (good for when wband a lot finer than G-band)
    elif mode == 'average':
        #get edges of w and g bins as range given is the center of the range gate
        wup = rangew + wres/2.
        wdown = rangew - wres/2.
        gup = rangeg + gres/2.
        gdown = rangeg - gres/2.
        
        #for each range at each timepoint:
        print('averaging.......')
        for g in range(nzg):
            #find indices in w-band range that lie within boundaries of G-band Ze:
            indsw = (wup <= (gup[g] + gres/5)) & (wdown > (gdown[g]-gres/5)) #because of the &, the bin has to fully lie within the g-band volume. 
            #indsw = (wup <= gup[g]) | (wdown > gdown[g]) #with | the wband volume doesnt need to be fully within the gband volume...
            #print(g, np.sum(indsw), rangew[indsw], gdown[g], gup[g], rangew[indsw]-wres/2, rangew[indsw]+wres/2)
            
            #if multiple bins fit into each g-band bin: for each time in G-Ze, take the associated w-band measurements and average
            #if only one bin fits in to g-band bin: take nearest neighbor.                
            outinds[g] = np.sum(indsw)
            
            
            #edit 07.07.25: get rid of time loop here as vertical ranging does not change with time!
            #for t in range(nt):
            #    if np.sum(indsw) > 1:
            #        Zew_rc[t, g] = np.nanmean(WZe.values[t, indsw])
            #        Zew_rc_std[t, g] = np.nanstd(WZe.values[t, indsw])
            #        
            #        vdw_rc[t, g] = np.nanmean(Wvd.values[t, indsw])
            #        vdw_rc_std[t, g] = np.nanstd(Wvd.values[t, indsw])
            #    elif np.sum(indsw) <= 1:
            #        diffrange = np.abs(rangeg[g] - rangew)
            #        outindsnearest = np.argmin(diffrange)
            #        Zew_rc[t, g] = WZe.values[t, outindsnearest]
            #        vdw_rc[t, g] = Wvd.values[t, outindsnearest]
            if np.sum(indsw) > 1:
                Zew_rc[:, g] = np.nanmean(WZe.values[:, indsw],axis=1)
                Zew_rc_std[:, g] = np.nanstd(WZe.values[:, indsw], axis=1)
                
                vdw_rc[:, g] = np.nanmean(Wvd.values[:, indsw], axis=1)
                vdw_rc_std[:, g] = np.nanstd(Wvd.values[:, indsw], axis=1)
            elif np.sum(indsw) <= 1:
                
                diffrange = np.abs(rangeg[g].values - rangew.values)
                outindsnearest = np.argmin(diffrange)
                Zew_rc[:, g] = WZe.values[:, outindsnearest]
                vdw_rc[:, g] = Wvd.values[:, outindsnearest]
                
    #if mode option is not set properly:
    else:
        print('did not understand your mode. options: nearest or average')
    
    outdic = {'WZe':Zew_rc, 'Wvd':vdw_rc, 'WZestd':Zew_rc_std, 'Wvdstd':vdw_rc_std}
    return outdic



def time_matching(WGrangematch, wdata, gdata):
    '''
    match wband radar data temporally with gband data; ie find closest wband measurement to each gband measurement. output dimension: gband height, gband time.
    input:
    - WGrangematch: xarray dataset with moments matched in range (output of range_match()) on their original time resolution
    - gdata: level0 data
    - wdata: level0 data
    '''
    nchirps = wdata.nchirp.shape[0]
    nt = WGrangematch.gtime.shape[0]
    nz = WGrangematch.l1height.shape[0]
    
    #calculate a timestampe giving out the time in the middle of each chirp sequence. output dimension: (nt x nchirps)
    WSeqTime = get_seqtime(wdata.time.values, wdata.chirpduration.values, nchirps)
    GSeqTime = get_seqtime(gdata.time.values, gdata.chirpduration.values, nchirps)

    #now loop through chirps to find closest wband measurement for each gband timestamp
    timeshift = []
    WZe = []
    Wvd = []
    #starting index:
    #get start and end indices in this way: [wstartidx:wendidx] (ie end index is not called)
    gstartidx = [int(x) for x in gdata.chirpseq_startix.values]
    
    gendidx = [int(gdata.chirpseq_startix.values[i+1]) for i in range(nchirps-1)]
    gendidx.append((gdata.Ze.shape[1]-1))
    
    for n in range(nchirps):
        
        #now find the nearest in W-band time to each g-band time using ckdtree:
        chirpbtree = cKDTree(WSeqTime[:,n].reshape(-1,1))
        chirpdist, chirpidx = chirpbtree.query(GSeqTime[:,n].reshape(-1,1))
        print('chirp %i: %i range gates were identified and shifted'%(n+1,np.sum(np.diff(chirpidx) > 1)))
        
        timeshift.append(chirpdist)
        if n < nchirps-1:
            WZe.append(WGrangematch.WZe.values[chirpidx, gstartidx[n]:gendidx[n]])
            Wvd.append(WGrangematch.Wvd.values[chirpidx, gstartidx[n]:gendidx[n]])
        else:
            WZe.append(WGrangematch.WZe.values[chirpidx, gstartidx[n]:])
            Wvd.append(WGrangematch.Wvd.values[chirpidx, gstartidx[n]:])
    WZematched = np.concatenate(WZe, axis=1)
    Wvdmatched = np.concatenate(Wvd, axis=1)
    
    #now also calculate DFR and DDV
    DFR = WZematched / gdata.Ze.values
    DDV = Wvdmatched / gdata.vd.values
    
    #map timeshift to time, height grid:
    timeshiftoutall = []
    for n in range(nchirps):
        if n < nchirps-1:
            newdim = gendidx[n] - gstartidx[n]
        else:
            newdim = gendidx[n]+1 - gstartidx[n]
        ja = np.expand_dims(timeshift[n],axis=1)
        timeshiftoutall.append(np.broadcast_to(ja, (nt, newdim)))
    
    timeshiftout = np.concatenate(timeshiftoutall, axis=1)
    
    #output: xarray dataset with coordinates Wtime, Gtime, l1height (means here G-height)
    coords = {'l1time':(['l1time'], gdata.time.values, {'units':'UTC', 'long_name':'Grawac time'}),
              'l1height':(['l1height'], gdata.height.values, {'units':'m', 'long_name':'Level-1 height, same as level0 gband height'})
              #'nchirp':(['nchirp'], np.arange(nchirps))
             }
    
    data_vars = { 
    'WZe' : (['l1time','l1height'], WZematched, {'units': 'mm6 m-3', 'long_name':'Wband radar reflectivity Ze matched to Gband range and time'}),
    'GZe' : (['l1time','l1height'], gdata.Ze.values, {'units': 'mm6 m-3', 'long_name':'Gband Ze'}),
    'Wvd' : (['l1time','l1height'], Wvdmatched, {'units': 'ms-1', 'long_name':'Wband mean Doppler velocity vd matched to Gband range and time'}),
    'Gvd' : (['l1time','l1height'], gdata.vd.values, {'units': 'ms-1', 'long_name':'Gband vd'}),
    'timeshift': (['l1time','l1height'], timeshiftout, {'units': 's', 'long_name':'Time shift of Wband radar compared to Gband sequence times'}),
    'DFR' : (['l1time','l1height'], DFR, {'units': '-', 'long_name':'Dual-Frequency Ratio DFR Wband/Gband'}),
    'DDV' : (['l1time','l1height'], DDV, {'units': '-', 'long_name':'Dual-Doppler Velocity DDV Wband/Gband'})
    }
    
    
    attrs = {}
    attrs['title'] = 'Radar moments range and time matched with calculated differential moments.'
    attrs['creation_time'] = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
        
    l1xrds = xr.Dataset(data_vars, coords, attrs)
    
    return l1xrds




def create_dataset(wgtimematch, wdata, gdata, info, write=False):
    '''
    input:
    - wgtimematch: xr dataset with time-range matched Zes, DFR, DDV
    - gdata: xr dataset with level0 gband
    - wdata: xr dataset with level0 wband
    - write: True/False; if True then netcdf output is written
    output:
    - level1: level1 dataset
    '''
    
    coords = {'time':(['time'], wgtimematch.l1time.values, {'units':'seconds since 1/1/1970', 'long_name':'time'}),
              'height':(['height'], wgtimematch.l1height.values, {'units':'m', 'long_name':'height asl (range+instrumentaltitude)'}),
              'nchirp':(['nchirp'], wdata.nchirp.values, {'long_name':'chirp sequences'})
             }
    
    data_vars = { 
    'WZe' : (['time','height'], wgtimematch.WZe.values, {'units': 'mm6 m-3', 'long_name':'Ze at 94GHz'}),
    'GZe' : (['time','height'], wgtimematch.GZe.values, {'units': 'mm6 m-3', 'long_name':'Ze at 167GHz'}),
    'G2Ze' : (['time','height'], gdata.Ze2.values, {'units': 'mm6 m-3', 'long_name':'Ze at 175GHz'}),
    'DAR' : (['time','height'], gdata.DAR.values, {'units': '-', 'long_name':'Differential Absorption Radar Dual-Frequency Ratio 167-175'}),
    'DFR' : (['time','height'], wgtimematch.DFR.values, {'units': '-', 'long_name':'Dual-Frequency Ratio DFR W-G'}),
    'Wvd' : (['time','height'], wgtimematch.Wvd.values, {'units': 'ms-1', 'long_name':'mean Doppler velocity vd at 94GHz'}),
    'Gvd' : (['time','height'], wgtimematch.Gvd.values, {'units': 'ms-1', 'long_name':'mean Doppler velocity vd at 167GHz'}),
    'DDV' : (['time','height'], wgtimematch.DDV.values, {'units': '-', 'long_name':'Dual-Doppler Velocity DDV W-G'}),
    'timeshift': (['time','height'], wgtimematch.timeshift.values, {'units': 's', 'long_name':'Time shift of Wband radar compared to Gband sequence times'}),
    'Gchirpidx': (['nchirp'], np.asarray(gdata.chirpseq_startix.values, dtype=int), {'units': '-', 'long_name':'Chirp sequence start index in height dimension'})
    }
    
    
    attrs = {
        "title": "Joint dataset of GRaWAC and JOYRAD94 (Level1)",
        "doi": "",
        "campaign_id":info['global']['mission'],
        "platform": info['global']['platform'],
        "instrument": "GRaWAC, JOYRAD94",
        "institution": "University of Cologne, Institute for Geophysics and Meteorology",
        "product_id":"Level-1",
        "version":info['global']['version'],
        "references": "AMT",
        "author":"Sabrina Schnitt",
        "author_email":"s.schnitt(at)uni-koeln.de",
        'creation_time':dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
        }
    
    l1 = xr.Dataset(data_vars, coords, attrs)
    
    if write==True:
        ncfile = info['global']['mission'] +'_%s_level1_draco_%s%s%s.nc'%(info['global']['version'], info['yyyy'], info['mm'], info['dd'])
        l1.to_netcdf(info['paths']['output'] + ncfile)
    return l1




