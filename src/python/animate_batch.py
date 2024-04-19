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
    
os.system('rm /tmp/' + username + '/frame*.png')
os.system('rm /tmp/' + username + '/animation*')
    

u_jet=[50., 100., 150];
c_vis=[0., 0.1, 0.2, 1.0]

fileNames=[['/tmp/' + username + '/output_' + str(i) + '_' + str(j) + '.nc' \
    for j in range(len(c_vis))] for i in range(len(u_jet))]

nc=NetCDFFile(fileNames[0][0])
lons=nc.variables['phi'][:]
lats=nc.variables['theta'][:]
vort=nc.variables['vort'][:]
h=nc.variables['h'][:]
v=nc.variables['v'][:]
time=nc.variables['time'][:]
nc.close()


f=plt.figure(figsize=(10,10))
tot1=len(u_jet)*len(c_vis)
iter1=0
for it1 in range(3,len(time)+1,4):
    it=it1-1
    print(str(it) + ' of ' + str(len(time)))
    if iter1==0:
        ax2=[]
        k=0
        cs=[]
        for i in range(len(u_jet)):
            ax1=[]
            for j in range(len(c_vis)):            
                ax1.append(f.add_subplot(len(u_jet),len(c_vis),k+1))
                k=k+1

                nc=NetCDFFile(fileNames[i][j])
                vort=nc.variables['vort'][:]
                nc.close()
                # set up orthographic map projection with 
                # perspective of satellite looking down at 50N, 100W.
                # use low resolution coastlines.
                map = Basemap(satellite_height=3000000,projection='nsper',\
                    lat_0=90,lon_0=-100,resolution='l')
                # draw the edge of the map projection region (the projection limb)
                # draw lat/lon grid lines every 30 degrees.
                map.drawmeridians(np.arange(0,360,30))
                map.drawparallels(np.arange(-90,90,30))

                (lo,la)=np.meshgrid(lons,lats)
                x,y = map(lo*180./np.pi, la*180./np.pi)
                # contour data over the map.
                #cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
                # contour data over the map.
                #cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
                cs1 = map.pcolor(x,y,vort[it,:,:],cmap='jet')
#                 cbar=map.colorbar(location='bottom')
#                 cbar.ax.set_xticks(cbar.ax.get_xticks())
#                 cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
#                 cbar.set_label('$\zeta$, s$^{-1}$')
                
                if i==0:
                    ax1[-1].set_title('$C_{vis}=' + str(c_vis[j]) + '$')
                if j==0:
                    ax1[-1].set_ylabel('$U_{jet}=' + str(u_jet[i]) + '$ (m s$^{-1}$)')
                
                cs.append(cs1)
    
            ax2 += [ax1]
    else:
        k=0
        for i in range(len(u_jet)):
            for j in range(len(c_vis)):

                nc=NetCDFFile(fileNames[i][j])
                vort=nc.variables['vort'][:]
                nc.close()

                cs[k].set_array(vort[it,0:-1,0:-1].flatten())
                cs[k].set_clim(np.min(vort[it,0:-1,0:-1].flatten()), \
                    np.max(vort[it,0:-1,0:-1].flatten()))
                k=k+1

    f.suptitle('Vorticity at t=' + '{0:.2f}'.format(time[it]/(86400)) +' days',y=0.9)
    plt.subplots_adjust()
    iter1 += 1
    plt.savefig('/tmp/' +username + '/frame%03d.png' % iter1,format='png', dpi=300) 
#os.system('convert -delay 20 -dispose previous /tmp/'  + username +'/frame*.png /tmp/' +username + '/animation.gif')
os.system('ffmpeg -r 5 -f image2  -i /tmp/' + username +  \
    '/frame%03d.png -vframes 34  -vcodec libx264 -crf 25  -pix_fmt yuv420p  /tmp/' + username + '/animation_batch.mp4')
os.system('ffmpeg -i /tmp/' + username +  '/animation_batch.mp4  /tmp/' + username + '/animation_batch.gif')
os.system('rm /tmp/' + username + '/frame*.png')
os.system('rm /tmp/' + username + '/animation_batch.mp4')



