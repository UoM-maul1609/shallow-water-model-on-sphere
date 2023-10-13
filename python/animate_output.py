import os
import getpass
import matplotlib
import sys
matplotlib.use('Agg')

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from netCDF4 import Dataset as NetCDFFile

import warnings
warnings.filterwarnings("ignore")


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

username=getpass.getuser()

if not os.path.exists('/tmp/' + username):
    os.mkdir('/tmp/' + username)

#nc=NetCDFFile('/tmp//output.nc')
nc=NetCDFFile('/tmp/' + username + '/output.nc')
lons=nc.variables['phi'][:]
lats=nc.variables['theta'][:]
vort=nc.variables['vort'][:]
h=nc.variables['h'][:]
v=nc.variables['v'][:]
time=nc.variables['time'][:]
iter1=0
for it1 in range(3,len(time)+1,4):
    it=it1-1
    if iter1==0:
        f=plt.figure()
        ax1=f.add_subplot(131)



        # set up orthographic map projection with 
        # perspective of satellite looking down at 50N, 100W.
        # use low resolution coastlines.
        map = Basemap(satellite_height=3000000,projection='nsper',lat_0=90,lon_0=-100,resolution='l')
        # draw the edge of the map projection region (the projection limb)
        # draw lat/lon grid lines every 30 degrees.
        map.drawmeridians(np.arange(0,360,30))
        map.drawparallels(np.arange(-90,90,30))

        (lo,la)=np.meshgrid(lons,lats)
        x,y = map(lo*180./np.pi, la*180./np.pi)
        # contour data over the map.
        #cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
        cs1 = map.pcolor(x,y,h[it,:,:],cmap='jet')
        cbar=map.colorbar(location='bottom')
        
        cbar.ax.set_xticks(cbar.ax.get_xticks())
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        cbar.set_label('h, m')
        ax1.set_title('Height at t=' + '{0:.2f}'.format(time[it]/(86400)) +' days')


        ax2=f.add_subplot(132)


        # set up orthographic map projection with 
        # perspective of satellite looking down at 50N, 100W.
        # use low resolution coastlines.
        map = Basemap(satellite_height=3000000,projection='nsper',lat_0=90,lon_0=-100,resolution='l')
        # draw the edge of the map projection region (the projection limb)
        # draw lat/lon grid lines every 30 degrees.
        map.drawmeridians(np.arange(0,360,30))
        map.drawparallels(np.arange(-90,90,30))

        # contour data over the map.
        #cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
        cs2 = map.pcolor(x,y,vort[it,:,:],cmap='jet')
        cbar=map.colorbar(location='bottom')
        cbar.ax.set_xticks(cbar.ax.get_xticks())
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        cbar.set_label('$\zeta$, s$^{-1}$')
        ax2.set_title('Vorticity')

        ax3=f.add_subplot(133)


        # set up orthographic map projection with 
        # perspective of satellite looking down at 50N, 100W.
        # use low resolution coastlines.
        map = Basemap(satellite_height=3000000,projection='nsper',lat_0=90,lon_0=-100,resolution='l')
        # draw the edge of the map projection region (the projection limb)
        # draw lat/lon grid lines every 30 degrees.
        map.drawmeridians(np.arange(0,360,30))
        map.drawparallels(np.arange(-90,90,30))

        # contour data over the map.
        #cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
        cs3 = map.pcolor(x,y,v[it,:,:],cmap='jet')
        cbar=map.colorbar(location='bottom')
        cbar.set_label('v, m s$^{-1}$')
        ax3.set_title('v')



        #plt.show(block=False)
        #plt.clim((61000,65000))
        # see matplotlib.org/basemap/users/examples.html
    else:
        ax1.set_title('Height at t=' + '{0:.2f}'.format(time[it]/(86400)) +' days')

        cs1.set_array(h[it,0:-1,0:-1].flatten())
        cs2.set_array(vort[it,0:-1,0:-1].flatten())
        cs3.set_array(v[it,0:-1,0:-1].flatten())

        cs1.set_clim(np.min(h[it,0:-1,0:-1].flatten()), np.max(h[it,0:-1,0:-1].flatten())) 
        cs2.set_clim(np.min(vort[it,0:-1,0:-1].flatten()), np.max(vort[it,0:-1,0:-1].flatten()))
        cs3.set_clim(np.min(v[it,0:-1,0:-1].flatten()), np.max(v[it,0:-1,0:-1].flatten()))

        #plt.show(block=False)

#     if not os.path.exists('/tmp/' + username):
#        os.mkdir('/tmp/' + username)
    if iter1==0:
       os.system('rm /tmp/' + username + '/frame*.png')
       os.system('rm /tmp/' + username + '/animation*')

    iter1 += 1
    plt.savefig('/tmp/' +username + '/frame%03d.png' % iter1,format='png', dpi=300) 
#os.system('convert -delay 20 -dispose previous /tmp/'  + username +'/frame*.png /tmp/' +username + '/animation.gif')
os.system('ffmpeg -r 5 -f image2  -i /tmp/' + username +  '/frame%03d.png -vframes 34  -vcodec libx264 -crf 25  -pix_fmt yuv420p  /tmp/' + username + '/animation.mp4')
os.system('ffmpeg -i /tmp/' + username +  '/animation.mp4  /tmp/' + username + '/animation.gif')
os.system('rm /tmp/' + username + '/frame*.png')
os.system('rm /tmp/' + username + '/animation.mp4')



