'''
processing routine to run ground-based processor
'''

import glob
import sys
import xarray as xr

#import all submodules
from level1 import make_level1 as l1
from level2 import make_level2 as l2
from watervapor import wv_retrieval as wv
from quicklooks import quicklooks_wv_evaluation_grawac_rs as wvql
from src import readdata as importfct
import json


print(sys.argv)
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
    yyyy = info['global']['yyyy']
    mm = info['global']['mm']
    dd = info['global']['dd']
    
info['yyyy'] = yyyy
info['mm'] = mm
info['dd'] = dd
#print(rr, tt)
#print(info['watervaporsettings']['R'],info['watervaporsettings']['tavg'] )
info['watervaporsettings']['R'] = rr
info['watervaporsettings']['tavg'] = '%ss'%tt
#print(info['watervaporsettings']['R'],info['watervaporsettings']['tavg'] )
print('processing %s, %s%s%s with water vapor retrieval on R=%sm and tavg=%s'%(info['global']['mission'], yyyy, mm, dd, info['watervaporsettings']['R'], info['watervaporsettings']['tavg']))
1/0
#print(info['yyyy'], info['mm'], info['dd'])

print('loading Gband files....')
gdatafiles = sorted(glob.glob(info['paths']['grawac'] +'%s/%s/%s/'%(yyyy,mm,dd) +'grawac_%s%s%s_*_%s_ZEN.lv1.NC'%(yyyy, mm, dd, info['paths']['grawacchirpprogram'])))
print(gdatafiles)

grawacl0data = importfct.read_rpg_lv1(gdatafiles, info, 'GRaWAC', instrumentaltinput = 15, withdar=True, write=True)

print('loading Wband files....')
wdatafiles  = sorted(glob.glob(info['paths']['mirac']+'%s/%s/%s/'%(yyyy,mm,dd) + 'joyrad94_nya_lv1a_%s%s%s*_%s.nc'%(yyyy, mm, dd, info['paths']['miracchirpprogram'])))
whkfiles = sorted(glob.glob(info['paths']['mirac']+'%s/%s/%s/'%(yyyy,mm,dd) + 'joyrad94_nya_housekeep_lv1a_%s%s%s*_%s.nc'%(yyyy, mm, dd, info['paths']['miracchirpprogram'])))
wbandl0data = importfct.read_compactfiles(wdatafiles, whkfiles, info, 'Mirac', withdar = False, write=True)

#level1 ================================
#todo: add rangematchmodes here
l1data = l1.run(grawacl0data, wbandl0data, info, write=True)

if info['global']['quicklooks'] == True:
    l1.quicklooks(l1data, info, write=True)

#attenuation: ================================
# include a check whether attenuation is there already; if not: run the code to produce output in designated file
#if attpath is empty: 
#    run attenuation
#else:
#    attenuation = load attenuation files

attpath = info['paths']['attenuation']
#load all profiles for all campaign duration (for post-processing); this needs to be changed for day-to-day processing!
attfiles = sorted(glob.glob(attpath + '%s_*_attenuation.nc'%(info['global']['mission']))) #files with gas attenuation calculated for each sonde profile
#print(attfiles)

geometry = 'bu' #bu: bottom-up; td: top-down
slant = False #False can be entered for zenith/nadir; otherwise, enter the angle / deg.

att = l2.attenuation_correction(attfiles, geometry, slant, info, write=False)

#level2  ================================ 
print('getting level2...')
l2data = l2.run(l1data, att, info, write=True)

#water vapor processing  ================================
thermo = xr.open_dataset(info['paths']['thermo'])
lut = xr.open_dataset(info['paths']['lut'])

print('retrieving water vapor...')
wvorig, wvsmooth = wv.run_retrieval_ground(l2data, thermo, lut, info, write=True)  #returns wv as initial output, and as smoothed output

# quicklooks
print('creating quicklooks...')
wvql.make_ql(l2data, wvorig, wvsmooth, thermo, info, write=True)

print('\n done. \n')
