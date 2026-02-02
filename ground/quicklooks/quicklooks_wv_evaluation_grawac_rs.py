import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import pandas as pd
import matplotlib.dates as mdates


def nearest(items, pivot):
    mintime = min(items, key=lambda x: abs(x - pivot)) 
    minarg = np.where(items == mintime)[0][0]
    return mintime, minarg


def make_ql(l2, wvorig, wvsmooth, thermo, info, write=False):
    '''
    input: 
    - l2: level2dataset
    - wvorig: water vapor retrieved initial resolution
    - wvsmooth: smoothed wv output
    - thermo: dataset of thermodynamic profiles to be evaluated against, e.g. radiosondes
    - info: paths
    - write: optional, if True: write to png.
    '''
    
    #change rs launchtime to numpydatetime64 coordinates:
    nt = thermo.launchtime.shape[0]
    rslaunchtimestr = np.array(thermo.launchtime.values,dtype='str')
    thermo = thermo.assign_coords( launchtime=np.asarray([dt.datetime.strptime(str(thermo.launchtime.values[t]), '%Y%m%d%H%M') for t in range(nt)],dtype=object))
    
    profilesthisday = [91,92]
    doy = pd.to_datetime(info['yyyy']+info['mm']+info['dd'])
    thermodoy = thermo.sel(launchtime=slice(doy, doy+pd.Timedelta(hours=24)))
    
    #loop through reference profiles:
    for rsind in range(thermodoy.launchtime.shape[0]):
        
        #find closest index to reference profile timestamp in l2, wv datasets
        l2ind = nearest(l2.time.values, thermodoy.launchtime[rsind].values)[1]
        wvind = nearest(wvsmooth.time.values, thermodoy.launchtime[rsind].values)[1]
        wvindo = nearest(wvorig.time.values, thermodoy.launchtime[rsind].values)[1]
        
        #set up figure:
        fig, axcts = plt.subplots(2,1,sharex=True, sharey=True, dpi=100, figsize=(16,9))
        fig.subplots_adjust(left=0.1, bottom=0.55, right=1, top=0.95)

        fig.suptitle('IOP4H2O NYA %s %s:%s UTC'%(rslaunchtimestr[rsind][:8], rslaunchtimestr[rsind][8:10], rslaunchtimestr[rsind][10:12]))

        axze = fig.add_axes([0.1, 0.1, 0.25, 0.37])
        axze.set_xlabel('Ze / dBZ')

        axdar = fig.add_axes([0.4, 0.1, 0.25, 0.37])
        axdar.set_xlabel('DAR / dB')

        axwv = fig.add_axes([0.7, 0.1, 0.25, 0.37])
        axwv.set_xlabel(r'$\rho_v$'+' / $gm^{-3}$' )
        
        #get cloud top height:
        
        for ax in [axcts[0], axcts[1], axze, axdar, axwv]:
            ax.set_ylabel('Height / m')
            #ax.set_ylim(0,8000)

        #now plot data:========================

        #contour plots:
        l2atrs = pd.to_datetime(l2.time.values[l2ind])
        beforers = l2atrs - pd.Timedelta(minutes=5)
        afterrs = l2atrs + pd.Timedelta(minutes=5)

        select = l2.sel(time=slice(beforers, afterrs))
        
        #also limit ylim to cloud top height:
        cth = select.height.values[np.nanmax(np.where(~np.isnan(select.GZe.values))[1])]
        for ax in [axcts[0], axcts[1], axze, axdar, axwv]:
            ax.set_ylim(0,cth+300)
        
        gzect = axcts[0].pcolormesh(select.time, select.height, select.GZe.T, cmap='Spectral', vmin=-55, vmax=20)
        plt.colorbar(gzect, label='Ze 167 / dBZ', ax=axcts[0])

        darct=axcts[1].pcolormesh(select.time, select.height, select.DAR.T, cmap='YlGnBu_r', vmin=0, vmax=10)
        plt.colorbar(darct, label='DAR / dB', ax=axcts[1])

        for i in range(2):
            axcts[i].axvline(l2.time[l2ind].values, linestyle = '--', color='grey')

        # Ze profiles:
        #plot profiles
        axze.plot(l2.GZe[l2ind,:], l2.height,  color='cyan', label='167.3',lw=1)
        axze.plot(l2.G2Ze[l2ind,:], l2.height, color='orange', label='174.7',lw=1)
        axze.set_xlim(-55,20)

        #DAR plot
        axdar.plot(l2.DAR[l2ind, :], l2.height, '.k',lw=1)
        axdar.set_xlim(0,15)
        
        
        #for Ze and DAR plots: plot shading #plot shading of +-std within 30s around launchtime
        tmin = l2atrs - pd.Timedelta(seconds=30)
        tmax = l2atrs + pd.Timedelta(seconds=30)
        zestdg = l2.GZe.sel(time=slice(tmin, tmax)).std(axis=0)
        zestdg2 = l2.G2Ze.sel(time=slice(tmin, tmax)).std(axis=0)
        darstd = l2.DAR.sel(time=slice(tmin, tmax)).std(axis=0)

        axze.fill_betweenx(l2.height.values, l2.GZe[l2ind,:].values - zestdg.values/2., l2.GZe[l2ind,:].values + zestdg.values/2.,color='cyan',alpha=0.4)
        axze.fill_betweenx(l2.height.values, l2.G2Ze[l2ind,:].values - zestdg2.values/2., l2.G2Ze[l2ind,:].values + zestdg2.values/2.,color='orange',alpha=0.4)
        axdar.fill_betweenx(l2.height.values, l2.DAR[l2ind,:].values - darstd.values/2., l2.DAR[l2ind,:].values + darstd.values/2.,color='grey',alpha=0.4)
        
        ## water vapor
        axwv.plot(thermodoy['rhov'][rsind], thermodoy['height'],'-k',label='RS')
        axwv.plot(wvsmooth.rhov[wvind,:], wvsmooth.height, '-r', label='GRaWAC smooth')
        axwv.plot(wvorig.rhov[wvindo,:], wvorig.height, '.', color='grey', ms=1, label='GRaWAC orig')
        axwv.set_xlim(-1,10)
        axwv.text(0.6, 0.7, 'RS IWV = %.2f kgm$^{-2}$'%thermodoy.iwv.values[rsind], transform=axwv.transAxes)
        
        ## legends
        for ax in [axze, axdar, axwv]:
            ax.legend(frameon=False, loc=1)

        myFmt = mdates.DateFormatter('%H:%M')
        axcts[1].xaxis.set_major_formatter(myFmt)
        axcts[1].set_xlabel('Time UTC')
        
        if write==True:
            doylaunchstr = str(thermodoy.launchtime.values[rsind]).split('T')[1][:2] + str(thermodoy.launchtime.values[rsind]).split('T')[1][3:5]
            outfile = info['paths']['quicklooks'] + '%s%s%s_watervapor_profile_grawac_rs%s.png'%(info['yyyy'], info['mm'],info['dd'],doylaunchstr )
            fig.savefig(outfile)
        plt.close()
    return
