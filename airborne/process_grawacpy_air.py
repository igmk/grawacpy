'''
processing routine to run airborne processor
'''

import glob
import sys
import xarray as xr
import numpy as np
#import all submodules
from level2 import make_level2 as l2
from level3 import make_level3 as l3a
from level3.watervapor import wv_retrieval as wv
#from quicklooks import quicklooks_wv_evaluation_grawac_rs as wvql
from src import readdata as importfct
import json

# ===================================================== 
#if code is called with sysargv: replace info['global']['yyyymmdd'] with input timestring.
if len(sys.argv) > 1:
    # dates to be processed:
    yyyy = sys.argv[1] #year
    mm = sys.argv[2] #month
    dd = sys.argv[3] #day
    if int(dd) < 10:
        dd=str('0%s'%dd)
    #config file name:
    configname = sys.argv[4]
    #opening json file for paths:
    with open("./%s.json"%configname) as json_file:
        info = json.load(json_file)
    
    #water vapor retrieval settings:
    rr = sys.argv[5] #R
    tt = sys.argv[6] #tavg
    
else: #specify through info json file:
    with open("./config_hamag.json") as json_file:
        info = json.load(json_file)
    yyyy = info['global']['yyyy']
    mm = info['global']['mm']
    dd = info['global']['dd']
    rr = info['watervaporsettings']['R']
    tt = info['watervaporsettings']['tavg']
    
info['yyyy'] = yyyy
info['mm'] = mm
info['dd'] = dd

yy = yyyy[2:]

#specify water vapor retrieval settings by reading config file variables:
info['watervaporsettings']['R'] = rr
info['watervaporsettings']['tavg'] = '%ss'%tt

# ===================================================== Level - 1: read matlab output
print('processing %s, %s%s%s with water vapor retrieval on R=%sm and tavg=%s'%(info['global']['mission'], yyyy, mm, dd, info['watervaporsettings']['R'], info['watervaporsettings']['tavg']))

print('Level-2....loading G-band')
#gdatafiles = sorted(glob.glob(info['paths']['grawac'] +'%s/%s/%s/'%(yyyy,mm,dd) +'grawac_%s%s%s_*_%s_ZEN.lv1.NC'%(yyyy, mm, dd, info['paths']['grawacchirpprogram'])))
gdatafiles  = sorted(glob.glob(info['paths']['grawac']+'%s/%s/%s/'%(yyyy,mm,dd) + 'grawac_p6_lv1a_%s%s%s*_%s.nc'%(yy, mm, dd, info['paths']['grawacchirpprogram'])))
ghkfiles = sorted(glob.glob(info['paths']['grawac']+'%s/%s/%s/'%(yyyy,mm,dd) + 'grawac_p6_housekeep_lv1a_%s%s%s*_%s.nc'%(yy, mm, dd, info['paths']['grawacchirpprogram'])))
gspfiles = sorted(glob.glob(info['paths']['grawac']+'%s/%s/%s/'%(yyyy,mm,dd) + 'grawac_p6_spectra_lv1a_%s%s%s*_%s.nc'%(yy, mm, dd, info['paths']['grawacchirpprogram'])))
#print(gdatafiles)

#grawacl0data = importfct.read_rpg_lv1(gdatafiles, info, 'GRaWAC', instrumentaltinput = 15, withdar=True, write=True)
grawacl1data = importfct.read_compactfiles(gdatafiles, ghkfiles, gspfiles, info, 'GRaWAC', withdar=True, write=True)


print('Level-2....loading W-band')
wdatafiles  = sorted(glob.glob(info['paths']['wband']+'%s/%s/%s/'%(yyyy,mm,dd) + 'mirac-a_pl6_lv1a_%s%s%s*_%s.nc'%(yy, mm, dd, info['paths']['miracchirpprogram'])))
whkfiles = sorted(glob.glob(info['paths']['wband']+'%s/%s/%s/'%(yyyy,mm,dd) + 'mirac-a_pl6_housekeep_lv1a_%s%s%s*_%s.nc'%(yy, mm, dd, info['paths']['miracchirpprogram'])))
wspfiles = sorted(glob.glob(info['paths']['wband']+'%s/%s/%s/'%(yyyy,mm,dd) + 'mirac-a_pl6_spectra_lv1a_%s%s%s*_%s.nc'%(yy, mm, dd, info['paths']['miracchirpprogram'])))
wbandl1data = importfct.read_compactfiles(wdatafiles, whkfiles, wspfiles, info, 'MIRAC', withdar = False, write=True)

# ===================================================== Level - 2: dual-frequency matching
print('Level-2....writing')
l2data = l2.run(grawacl1data, wbandl1data, info, write=True)

if info['global']['quicklooks'] == True:
    l2.quicklooks(l2data, info, write=True)
print('Level-2....done\n --------------\n')
# ===================================================== Level - 3a: attenuation correction ======
print('Level-3....loading GPS')
# add GPS information to level-2 dataset ==========================
gps = xr.open_dataset(glob.glob(info['paths']['gps']+'%s_*_GPS_INS_%s%s%s_*.nc'%(info['global']['mission'], yyyy, mm, dd))[0])
flightspecs = gps.interp(time=l2data.time)
#add these to l2dataset:
l2data = l2data.assign(lon=flightspecs['lon'], lat=flightspecs['lat'], alt=flightspecs['alt'], pitch=flightspecs['pitch'], roll=flightspecs['roll'])

#now calculate inclination angle and aircraft height
theta = 25. #inclination angle bellypod
incangle = -1*(flightspecs.pitch - theta)

#calculate height [m asl] of each radar range bin:
radarhgt = l2data.alt - l2data.range*np.cos(np.radians(incangle))
radarhgt.attrs['units'] = 'm asl'
radarhgt.attrs['long_name'] = 'Height asl of each radar range bin'

#assign auxiliary gps and position data to l2dataset:
l2data = l2data.assign(height=radarhgt, incangle=incangle)

#attenuation: ================================
# include a check whether attenuation is there already; if not: run the code to produce output in designated file
#if attpath is empty: 
#    run attenuation
#else:
#    attenuation = load attenuation files
print('Level-3....loading attenuation')
attpath = info['paths']['attenuation']
attfiles = sorted(glob.glob(attpath + '%s/%s/%s/%s_*_attenuation.nc'%(yyyy,mm,dd,info['global']['mission']))) #files with gas attenuation calculated for each sonde profile

#check if attenuation profiles are available; if yes, read them and keep going; if no: run attenuation script to produce output
if len(attfiles) < 1:
#if 1 == 1:
    from attenuation import get_attenuation_from_sondes as att
    att.run(info)
    #find files again, and check that they got saved.
    attfiles = sorted(glob.glob(attpath + '%s/%s/%s/%s_*_attenuation.nc'%(yyyy,mm,dd,info['global']['mission']))) #files with gas attenuation calculated for each sonde profile
    assert len(attfiles) > 0
#load all profiles for all campaign duration (for post-processing); this needs to be changed for day-to-day processing!
#print(attfiles)

geometry = 'td' #bu: bottom-up; td: top-down
slant = True #False can be entered for zenith/nadir; otherwise, enter the angle / deg.

att = l3a.attenuation_correction(l2data, attfiles, geometry, l2data, info, write=False)

l3adata = l3a.run(l2data, att, info, write=True)

# ===================================================== Level - 3b: water vapor retrieval
print('Level-3....water vapor retrieval')
#thermo = xr.open_dataset(info['paths']['thermo'])
thermo = info['paths']['thermo'] + '%s_%s%s%s_%s/'%(info['global']['mission'],yyyy,mm,dd, info['global']['RF'] )
lut = xr.open_dataset(info['paths']['lut'])

print('retrieving water vapor...')
wvorig, wvsmooth = wv.run_retrieval_air(l3adata, thermo, lut, info, write=True)  #returns wv as initial output, and as smoothed output

# quicklooks
print('creating quicklooks...')
wvql.make_ql(l3adata, wvorig, wvsmooth, thermo, info, write=True)

print('\n done. \n')
