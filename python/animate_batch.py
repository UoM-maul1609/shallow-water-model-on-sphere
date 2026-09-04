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

os.system('rm /tmp/' + username + '/frame*.png')
os.system('rm /tmp/' + username + '/animation*')

u_jet = [50., 100., 150.]
c_vis = [0., 0.1, 0.2, 1.0]

fileNames = [['/tmp/' + username + '/output_' + str(i) + '_' + str(j) + '.nc'
              for j in range(len(c_vis))] for i in range(len(u_jet))]

nc = NetCDFFile(fileNames[0][0])
lons_raw = nc.variables['phi'][:]
lats = nc.variables['theta'][:]
time = nc.variables['time'][:]
nc.close()
lons, = strip_duplicated_endpoint(lons_raw)

# Construct cyclic plotting coordinates once. Data get the same temporary
# cyclic column for each file/frame below.
lons_plot, dummy = add_cyclic_longitude(lons, np.zeros((len(lats), len(lons))))
lo, la = np.meshgrid(lons_plot, lats)

f = plt.figure(figsize=(10, 10))
iter1 = 0
for it1 in range(3, len(time) + 1, 4):
    it = it1 - 1
    print(str(it) + ' of ' + str(len(time)))
    if iter1 == 0:
        ax2 = []
        k = 0
        cs = []
        for i in range(len(u_jet)):
            ax1 = []
            for j in range(len(c_vis)):
                ax1.append(f.add_subplot(len(u_jet), len(c_vis), k + 1))
                k += 1

                nc = NetCDFFile(fileNames[i][j])
                file_lons = nc.variables['phi'][:]
                vort_raw = nc.variables['vort'][it, :, :]
                nc.close()
                file_lons, vort_native = strip_duplicated_endpoint(file_lons, vort_raw)
                _, vort_plot = add_cyclic_longitude(file_lons, vort_native)

                map1 = Basemap(satellite_height=3000000, projection='nsper',
                               lat_0=90, lon_0=-100, resolution='l')
                map1.drawmeridians(np.arange(0, 360, 30))
                map1.drawparallels(np.arange(-90, 90, 30))
                x, y = map1(lo * 180. / np.pi, la * 180. / np.pi)
                cs1 = map1.pcolor(x, y, vort_plot, cmap='jet')

                if i == 0:
                    ax1[-1].set_title('$C_{vis}=' + str(c_vis[j]) + '$')
                if j == 0:
                    ax1[-1].set_ylabel('$U_{jet}=' + str(u_jet[i]) + '$ (m s$^{-1}$)')
                cs.append(cs1)
            ax2 += [ax1]
    else:
        k = 0
        for i in range(len(u_jet)):
            for j in range(len(c_vis)):
                nc = NetCDFFile(fileNames[i][j])
                file_lons = nc.variables['phi'][:]
                vort_raw = nc.variables['vort'][it, :, :]
                nc.close()
                file_lons, vort_native = strip_duplicated_endpoint(file_lons, vort_raw)
                _, vort_plot = add_cyclic_longitude(file_lons, vort_native)

                cs[k].set_array(vort_plot[:-1, :-1].flatten())
                cs[k].set_clim(np.nanmin(vort_native), np.nanmax(vort_native))
                k += 1

    f.suptitle('Vorticity at t=' + '{0:.2f}'.format(time[it] / 86400) + ' days', y=0.9)
    plt.subplots_adjust()
    iter1 += 1
    plt.savefig('/tmp/' + username + '/frame%03d.png' % iter1,
                format='png', dpi=300)

os.system('ffmpeg -r 5 -f image2 -i /tmp/' + username +
          '/frame%03d.png -vframes 34 -vcodec libx264 -crf 25 '
          '-pix_fmt yuv420p /tmp/' + username + '/animation_batch.mp4')
os.system('ffmpeg -i /tmp/' + username + '/animation_batch.mp4 /tmp/' +
          username + '/animation_batch.gif')
os.system('rm /tmp/' + username + '/frame*.png')
os.system('rm /tmp/' + username + '/animation_batch.mp4')
