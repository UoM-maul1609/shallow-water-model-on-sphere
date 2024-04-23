"""
Module Name: animate_batch.py

Description:
This module contains code to analyse the fourier
wave number of the output.nc netCDF file.
"""

import os
import getpass
import warnings
import matplotlib

import matplotlib.pyplot as plt
import numpy as np

from create_map import create_map_func

# pylint: disable=E0611
# (pylint can't find Dataset in the netCDF4 package for some reason)
from netCDF4 import Dataset as NetCDFFile
from animate_output import nc, lons, lats, vort, time

matplotlib.use("Agg")

warnings.filterwarnings("ignore")

# TODO: FIX AND TEST THIS MODULE pylint: disable=W0511

username = getpass.getuser()

if not os.path.exists("../../output/frames"):
    os.mkdir("../../output/frames")
if not os.path.exists("../../output/animations"):
    os.mkdir("../../output/animations")

u_jet = [50.0, 100.0, 150]
c_vis = [0.0, 0.1, 0.2, 1.0]

fileNames = [
    [
        "../../tests/output_" + str(i) + "_" + str(j) + ".nc"
        for j in range(len(c_vis))
    ]
    for i in range(len(u_jet))
]

f = plt.figure(figsize=(10, 10))
TOT = len(u_jet) * len(c_vis)
ITER = 0
for it1 in range(3, len(time) + 1, 4):
    it = it1 - 1
    print(str(it) + " of " + str(len(time)))
    if ITER == 0:
        ax2 = []
        k = 0
        cs = []
        for i, u_jet_val in enumerate(u_jet):
            ax1 = []
            for j, c_vis_val in enumerate(c_vis):
                ax1.append(f.add_subplot(len(u_jet), len(c_vis), k + 1))
                k = k + 1

                nc = NetCDFFile(fileNames[i][j])
                vort = nc.variables["vort"][:]
                nc.close()
                # set up orthographic map projection
                x, y, basemap = create_map_func(lons, lats)
                # contour data over the map.
                # cs = basemap.contour(x,y,wave+mean,15,linewidths=1.5)
                cs1 = basemap.pcolor(x, y, vort[it, :, :], cmap="jet")
                # cbar=basemap.colorbar(location='bottom')
                # cbar.ax.set_xticks(cbar.ax.get_xticks())
                # cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
                # cbar.set_label('$\zeta$, s$^{-1}$')

                if i == 0:
                    ax1[-1].set_title("$C_{vis}=" + str(c_vis_val) + "$")
                if j == 0:
                    ax1[-1].set_ylabel("$U_{jet}=" + str(u_jet_val) + "$ (m s$^{-1}$)")

                cs.append(cs1)

            ax2 += [ax1]
    else:
        k = 0
        for i in range(len(u_jet)):
            for j in range(len(c_vis)):
                nc = NetCDFFile(fileNames[i][j])
                vort = nc.variables["vort"][:]
                nc.close()

                cs[k].set_array(vort[it, 0:-1, 0:-1].flatten())
                cs[k].set_clim(
                    np.min(vort[it, 0:-1, 0:-1].flatten()),
                    np.max(vort[it, 0:-1, 0:-1].flatten()),
                )
                k = k + 1

    f.suptitle(f"Vorticity at t={time[it]/86400:.2f} days", y=0.9)
    plt.subplots_adjust()
    ITER += 1
    plt.savefig(f"../../output/frames/frame{ITER:03d}.png", format="png", dpi=300)

# os.system('convert -delay 20 -dispose previous /tmp/'  +\
#               username +'/frame*.png /tmp/' +username + '/animation.gif')
os.system(
    "ffmpeg -r 5 -f image2  -i /tmp/"
    + username
    + "/frame%03d.png -vframes 34  -vcodec libx264 -crf 25  -pix_fmt yuv420p  /tmp/"
    + username
    + "/animation_batch.mp4"
)
os.system(
    "ffmpeg -i /tmp/"
    + username
    + "/animation_batch.mp4  /tmp/"
    + username
    + "/animation_batch.gif"
)
