import os
import getpass
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from netCDF4 import Dataset as NetCDFFile
from scipy.interpolate import griddata

from longitude_utils import add_cyclic_longitude, strip_duplicated_endpoint

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
u_raw = nc.variables['u'][:]
v_raw = nc.variables['v'][:]
time = nc.variables['time'][:]
nc.close()

# Strip a duplicated 360-degree endpoint from old files, if present. Keep
# native distinct samples for interpolation/streamlines.
lons, h, u, v, vort = strip_duplicated_endpoint(
    lons_raw, h_raw, u_raw, v_raw, vort_raw
)

Re = 5.4155760e7
lo_native, la_native = np.meshgrid(lons, np.pi / 2. - lats)
arc_native = la_native * Re
x_native = arc_native * np.cos(lo_native)
y_native = arc_native * np.sin(lo_native)

# Close the height map only for display.
lons_plot, h_plot = add_cyclic_longitude(lons, h[-1, :, :])
lo_plot, la_plot = np.meshgrid(lons_plot, np.pi / 2. - lats)
arc_plot = la_plot * Re
x_plot = arc_plot * np.cos(lo_plot)
y_plot = arc_plot * np.sin(lo_plot)
hmap = plt.pcolor(x_plot, y_plot, h_plot)
plt.axis('square')

xx = np.linspace(-4.5e7, 4.5e7, 100)
yy = np.linspace(-4.5e7, 4.5e7, 100)
xx1, yy1 = np.meshgrid(xx, yy)

u11 = u[-1, :, :] - np.mean(u[-1, :, :], axis=0)
v11 = v[-1, :, :]

# polar-to-Cartesian velocity conversion, using native (nonduplicated) points
u111 = -v11 * np.cos(lo_native) - u11 * np.sin(lo_native)
v111 = -v11 * np.sin(lo_native) + u11 * np.cos(lo_native)

u111 = u111[0::5, 0::5].flatten()
v111 = v111[0::5, 0::5].flatten()

uu = griddata((x_native[0::5, 0::5].flatten(),
               y_native[0::5, 0::5].flatten()), u111, (xx1, yy1))
vv = griddata((x_native[0::5, 0::5].flatten(),
               y_native[0::5, 0::5].flatten()), v111, (xx1, yy1))

plt.streamplot(xx, yy, uu, vv, 4)
plt.xlim((-2e7, 2e7))
plt.ylim((-2e7, 2e7))

cbar = plt.colorbar(hmap)
cbar.set_label('h,m')
plt.title('Height at t=' + '{0:.2f}'.format(time[-1] / 86400) + ' days')

os.system('rm /tmp/' + username + '/frame*.png')
plt.savefig('/tmp/' + username + '/frame.png', format='png', dpi=300)
