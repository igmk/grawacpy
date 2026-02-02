def make_dataset(time, rrange, wze, gze, gze2, wvd, gvd, ddv, dfr, dar, wzecorr, gzecorr, gze2corr, dfrcorr, att, attrs, ncoutput, write=False):
    '''
    write all variables into xarray dataset and if set True, export as netcdf
    input:
    
    output:
    - xrds: xr dataset
    '''
    
    #create new data_vars:
    data_vars = { 
        'W_Ze_meas' : (['time','height'], wze, {'units': 'dBZ', 'long_name':'Measured Mirac radar reflectivity matched to GRaWAC range and time; affected by attenuation'}),
        'G_Ze_meas' : (['time','height'], gze, {'units': 'dBZ', 'long_name':'Measured GRaWAC Ze at 167.3GHz; affected by attenuation'}),
        'G_Ze2_meas' : (['time','height'], gze2, {'units': 'dBZ', 'long_name':'Measured GRaWAC Ze at 174.7GHz; affected by attenuation'}),
        'DFR_meas' : (['time','height'], dfr, {'units': 'dB', 'long_name':'Measured W-G Dual-Frequency Ratio; ZeW-ZeG; affected by attenuation'}),
        'DAR' : (['time','height'], dar, {'units': 'dB', 'long_name':'Differential Absorption Radar Dual-Frequency Ratio'}),
        'W_Ze'    : (['time','height'], wzecorr.values, {'units': 'dBZ', 'long_name':'Attenuation-corrected Mirac radar reflectivity matched to GRaWAC range and time; only gas attenuation is corrected'}),
        'G_Ze'    : (['time','height'], gzecorr.values, {'units': 'dBZ', 'long_name':'Attenuation-corrected GRaWAC radar reflectivity at 167.3GHz; only gas attenuation is corrected'}),
        'G_Ze2'    : (['time','height'], gze2corr.values, {'units': 'dBZ', 'long_name':'Attenuation-corrected GRaWAC radar reflectivity at 174.7GHz; only gas attenuation is corrected'}),
        'DFR'    : (['time','height'], dfrcorr.values, {'units': 'dB', 'long_name':'Attenuation-corrected dual-frequency W-G ratio (ZeW-ZeG); only gas attenuation is corrected'}),
        'DDV'    : (['time','height'], ddv, {'units': '-', 'long_name':'Differential Doppler Velocity'}),
        'W_vd' :(['time','height'], wvd, {'units': 'ms-1', 'long_name':'Wband mean doppler velocity'}),
        'G_vd' :(['time','height'], gvd, {'units': 'ms-1', 'long_name':'Gband mean doppler velocity'})    }
    
    coords = {'time':(['time'], time),
              'height':(['height'], rrange.values, {'units': 'm', 'long_name':'Radar range height msl'}),
              'freq':(['freq'], att.freq.values)
             }
    
    attrs['creation_time'] = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
    attrs = attrs
    
    xrds = xr.Dataset(data_vars, coords, attrs)
    
    #variables that already are xr datasets: att; lon; lat; 
    xrds = xrds.assign(att_atmo=att.att_atmo)
    #xrds = xrds.assign(W_Ze = wzecorr)
    
    #save to netcdf:
    print('....TODO: compress output netcdf file')
    if write:
        print('write to %s'%ncoutput)
        #check if ncoutput already available; if yes, delete it.
        if os.path.isfile(ncoutput):
            print('found old file. deleting old file.')
            os.remove(ncoutput)
        
        xrds.to_netcdf(ncoutput, 'w', format="NETCDF4")
    
    return xrds


