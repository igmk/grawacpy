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
import quicklooks as ql
from src import readdata as importfct
import json

#load retrieval/processing options and input/output paths
##must contain: 
#- inputfiletype
#- perspective (ground/air)
#- inputpath
#- attenuationpath
#- outputpath
#- R (


#opening json file for paths:
with open("./attrs_paths.json") as json_file:
    info = json.load(json_file)

    
#if code is called with sysargv: replace info['global']['yyyymmdd'] with input timestring.
if len(sys.argv) > 1:
    yyyy = sys.argv[1]
    mm = sys.argv[2]
    dd = sys.argv[3]
    if int(dd) < 10:
        dd=str('0%s'%dd)
    
else: #specify through info json file:
    yyyy = info['global']['yyyy']
    mm = info['global']['mm']
    dd = info['global']['dd']
    
print('processing %s, %s%s%s with water vapor retrieval on R=%sm'%(info['global']['mission'], yyyy, mm, dd, info['watervaporsettings']['R']))
info['yyyy'] = yyyy
info['mm'] = mm
info['dd'] = dd
#print(info['yyyy'], info['mm'], info['dd'])
print('loading Gband files....')
gdatafiles = sorted(glob.glob(info['paths']['grawac'] +'grawac_%s%s%s_*_%s_ZEN.lv1.NC'%(yyyy, mm, dd, info['paths']['grawacchirpprogram'])))
#print(gdatafiles)

grawacl0data = importfct.read_rpg_lv1(gdatafiles, info, 'GRaWAC', instrumentaltinput = 15, withdar=True, write=True)

print('loading Wband files....')
wdatafiles  = sorted(glob.glob(info['paths']['mirac'] + 'joyrad94_nya_lv1a_%s%s%s*_%s.nc'%(yyyy, mm, dd, info['paths']['miracchirpprogram'])))
whkfiles = sorted(glob.glob(info['paths']['mirac'] + 'joyrad94_nya_housekeep_lv1a_%s%s%s*_%s.nc'%(yyyy, mm, dd, info['paths']['miracchirpprogram'])))
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

l2data = l2.run(l1data, att, info, write=True)

#water vapor processing  ================================
thermo = xr.open_dataset(info['paths']['thermo'])
lut = xr.open_dataset(info['paths']['lut'])

wvdata = wv.run_retrieval_ground(l2data, thermo, lut, info, write=True)


