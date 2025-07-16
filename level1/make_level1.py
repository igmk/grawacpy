import numpy as np
import xarray as xr

from level1.src import level1_functions as fct


def run(gdata, wdata, info, write=False):
    '''
    input:
    - grawac: xarray dataset with original gband dataset, dimensions: height, time, chirps
    - wband: xarray dataset with original wdata dataset, dimensions: height, time, chirps
    - info: xarray dataset with info json file imported
    - write: optional; level1 dataset stored if write==True.
    
    '''
    
    nchirps = wdata.nchirp.shape[0]
    
    print('calculating chirp characteristics...')
    
    #calculate chirp repetition time and chirp duration for each input; add them to the dataset. (nchirp dimension)
    ChirpRepTime_W = fct.calculate_chirp_repetition_time(wdata.maxvel.values, wdata.freq.values*1e09)
    ChirpRepTime_G = fct.calculate_chirp_repetition_time(gdata.maxvel.values, gdata.freq.values*1e09)
    
    # calculate chirp duration
    ChirpDuration_W = ChirpRepTime_W * wdata.navg.values
    ChirpDuration_G = ChirpRepTime_G * gdata.navg.values
    
    #calculate chirp sequence time stamps:
    #calculate a timestampe giving out the time in the middle of each chirp sequence. output dimension: (nt x nchirps)
    WSeqTime = fct.get_seqtime(wdata.time, ChirpDuration_W, nchirps)
    GSeqTime = fct.get_seqtime(gdata.time, ChirpDuration_G, nchirps)
    
    
    wdata = wdata.assign(chirpreptime = (['nchirp'], ChirpRepTime_W), chirpduration = ([ 'nchirp'], ChirpDuration_W), chirpseqtime = (['time', 'nchirp'], WSeqTime))
    gdata = gdata.assign(chirpreptime = (['nchirp'], ChirpRepTime_G), chirpduration = ([ 'nchirp'], ChirpDuration_G), chirpseqtime = (['time', 'nchirp'], GSeqTime))
    
    
    
    #now do range matching ================
    ### todo: get the modes from level before
    #modes = ['average', 'average', 'average', 'nearest']
    print('range matching...')
    wgrangematch = fct.range_matching(wdata, gdata, info['level1settings']['mode'])
    
    #now do time matching and also calculate DFR and DDV ===================
    wgtimematch = fct.time_matching(wgrangematch, wdata, gdata)
    
    #create level1 xarray dataset:
    level1 = fct.create_dataset(wgtimematch, wdata, gdata, info, write=write)
    
    
    #flagging: all DFR with timeshift larger than 1 second is nan:
    level1.DFR.values[level1.timeshift.values > 1] = np.nan
    
    return level1

def quicklooks(l1, info, write=False):
    
    
    return
