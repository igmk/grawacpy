import numpy as np

def nearest(items, pivot):
    mintime = min(items, key=lambda x: abs(x - pivot)) 
    minarg = np.where(items == mintime)[0][0]
    return mintime, minarg



def get_zlin(ZdBZ):
    '''
    convert input Ze in dBZ to linear Ze
    in: Ze / dBZ
    out: Ze / mm^6/m3
    '''
    return 10**(ZdBZ/10.)



def get_ZdBZ(zlin):
    '''
    convert input Ze [mm^6/m3] to dB
    in: Ze / mm^6/m3
    out: Ze / dBZ
    '''
    return 10*np.log10(zlin)


def val2idx(c, value):
    
    nc = len(c)
    idx = np.round(np.interp(value, c, np.arange(nc)).clip(0,nc-1))
    
    return np.array(idx,dtype='int')



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
