import os
import getpass
import matplotlib
matplotlib.use('Agg')

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from netCDF4 import Dataset as NetCDFFile

from longitude_utils import add_cyclic_longitude, strip_duplicated_endpoint

import warnings
warnings.filterwarnings("ignore")

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

username = getpass.getuser()

if not os.path.exists('/tmp/' + username):
    os.mkdir('/tmp/' + username)

nc = NetCDFFile('/tmp/' + username + '/output.nc')
lons_raw = nc.variables['phi'][:]
lats = nc.variables['theta'][:]
vort_raw = nc.variables['vort'][:]
h_raw = nc.variables['h'][:]
v_raw = nc.variables['v'][:]
time = nc.variables['time'][:]
nc.close()

# Native arrays always contain distinct periodic samples. This also makes the
# script backward-compatible with old files containing both 0 and 2*pi.
lons, h, vort, v = strip_duplicated_endpoint(lons_raw, h_raw, vort_raw, v_raw)

# Plotting gets one temporary cyclic column so the map closes at 360 degrees.
lons_plot, h0 = add_cyclic_longitude(lons, h[0, :, :])
lo, la = np.meshgrid(lons_plot, lats)

iter1 = 0
for it1 in range(3, len(time) + 1, 4):
    it = it1 - 1
    _, h_plot = add_cyclic_longitude(lons, h[it, :, :])
    _, vort_plot = add_cyclic_longitude(lons, vort[it, :, :])
    _, v_plot = add_cyclic_longitude(lons, v[it, :, :])

    if iter1 == 0:
        f = plt.figure()
        ax1 = f.add_subplot(131)

        map1 = Basemap(satellite_height=3000000, projection='nsper',
                       lat_0=90, lon_0=-100, resolution='l')
        map1.drawmeridians(np.arange(0, 360, 30))
        map1.drawparallels(np.arange(-90, 90, 30))
        x1, y1 = map1(lo * 180. / np.pi, la * 180. / np.pi)
        cs1 = map1.pcolor(x1, y1, h_plot, cmap='jet')
        cbar = map1.colorbar(location='bottom')
        cbar.ax.set_xticks(cbar.ax.get_xticks())
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        cbar.set_label('h, m')
        ax1.set_title('Height at t=' + '{0:.2f}'.format(time[it] / 86400) + ' days')

        ax2 = f.add_subplot(132)
        map2 = Basemap(satellite_height=3000000, projection='nsper',
                       lat_0=90, lon_0=-100, resolution='l')
        map2.drawmeridians(np.arange(0, 360, 30))
        map2.drawparallels(np.arange(-90, 90, 30))
        x2, y2 = map2(lo * 180. / np.pi, la * 180. / np.pi)
        cs2 = map2.pcolor(x2, y2, vort_plot, cmap='jet')
        cbar = map2.colorbar(location='bottom')
        cbar.ax.set_xticks(cbar.ax.get_xticks())
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
        cbar.set_label('$\\zeta$, s$^{-1}$')
        ax2.set_title('Vorticity')

        ax3 = f.add_subplot(133)
        map3 = Basemap(satellite_height=3000000, projection='nsper',
                       lat_0=90, lon_0=-100, resolution='l')
        map3.drawmeridians(np.arange(0, 360, 30))
        map3.drawparallels(np.arange(-90, 90, 30))
        x3, y3 = map3(lo * 180. / np.pi, la * 180. / np.pi)
        cs3 = map3.pcolor(x3, y3, v_plot, cmap='jet')
        cbar = map3.colorbar(location='bottom')
        cbar.set_label('v, m s$^{-1}$')
        ax3.set_title('v')
    else:
        ax1.set_title('Height at t=' + '{0:.2f}'.format(time[it] / 86400) + ' days')

        # pcolor with centre-coordinate arrays of the same shape stores the
        # lower-left (M-1)x(N-1) values. The cyclic column means all N native
        # longitudes are retained; only the pre-existing final latitude row is
        # omitted by this plotting convention.
        cs1.set_array(h_plot[:-1, :-1].flatten())
        cs2.set_array(vort_plot[:-1, :-1].flatten())
        cs3.set_array(v_plot[:-1, :-1].flatten())

        cs1.set_clim(np.nanmin(h[it, :, :]), np.nanmax(h[it, :, :]))
        cs2.set_clim(np.nanmin(vort[it, :, :]), np.nanmax(vort[it, :, :]))
        cs3.set_clim(np.nanmin(v[it, :, :]), np.nanmax(v[it, :, :]))

    if iter1 == 0:
        os.system('rm /tmp/' + username + '/frame*.png')
        os.system('rm /tmp/' + username + '/animation*')

    iter1 += 1
    plt.savefig('/tmp/' + username + '/frame%03d.png' % iter1,
                format='png', dpi=300)

os.system('ffmpeg -r 5 -f image2 -i /tmp/' + username +
          '/frame%03d.png -vframes 34 -vcodec libx264 -crf 25 '
          '-pix_fmt yuv420p /tmp/' + username + '/animation.mp4')
os.system('ffmpeg -i /tmp/' + username + '/animation.mp4 /tmp/' +
          username + '/animation.gif')
os.system('rm /tmp/' + username + '/frame*.png')
os.system('rm /tmp/' + username + '/animation.mp4')
