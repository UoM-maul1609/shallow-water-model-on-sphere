"""
Module Name: animate_output.py

Description:
This module contains code to generate frames and an animation
of the height, vorticity, and u and v velocities for the output.nc
netCDF file.
"""

import os
import getpass
import warnings

import matplotlib
import matplotlib.pyplot as plt

from create_map import create_map_func

# pylint: disable=E0611
# (pylint can't find Dataset in the netCDF4 package for some reason)
from netCDF4 import Dataset as NetCDFFile

matplotlib.use("Agg")

warnings.filterwarnings("ignore")

username = getpass.getuser()

nc = NetCDFFile("../../tests/output.nc")
lons = nc.variables["phi"][:]
lats = nc.variables["theta"][:]
vort = nc.variables["vort"][:]
h = nc.variables["h"][:]
v = nc.variables["v"][:]
u = nc.variables["u"][:]
time = nc.variables["time"][:]

if not os.path.exists("../../output/frames"):
    os.mkdir("../../output/frames")
if not os.path.exists("../../output/animations"):
    os.mkdir("../../output/animations")

ITER = 0
for it1 in range(3, len(time) + 1, 4):
    it = it1 - 1
    # if ITER==0:
    f = plt.figure()

    ax1 = f.add_subplot(141)
    # set up orthographic map projection
    x, y, basemap = create_map_func(lons, lats)
    # contour data over the basemap.
    # cs = basemap.contour(x,y,wave+mean,15,linewidths=1.5)
    cs1 = basemap.pcolor(
        x, y, h[it, :, :], cmap="jet", shading="auto", vmin=59000, vmax=64000
    )
    cbar = basemap.colorbar(location="bottom")
    cbar.ax.set_xticks(cbar.ax.get_xticks())
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation="vertical")
    cbar.set_label("h, m")
    ax1.set_title(f"Height at \n t={time[it]/86400:.2f} days")

    ax2 = f.add_subplot(142)
    # set up orthographic map projection
    x, y, basemap = create_map_func(lons, lats)
    # contour data over the basemap.
    # cs = basemap.contour(x,y,wave+mean,15,linewidths=1.5)
    cs2 = basemap.pcolor(
        x, y, vort[it, :, :], cmap="jet", shading="auto", vmin=-0.00008, vmax=0.00008
    )
    cbar = basemap.colorbar(location="bottom")
    cbar.ax.set_xticks(cbar.ax.get_xticks())
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation="vertical")
    cbar.set_label("$\\zeta$, s$^{-1}$")
    ax2.set_title("Vorticity")

    ax3 = f.add_subplot(143)
    # set up orthographic map projection
    x, y, basemap = create_map_func(lons, lats)
    # contour data over the basemap.
    # cs = basemap.contour(x,y,wave+mean,15,linewidths=1.5)
    cs3 = basemap.pcolor(x, y, v[it, :, :], cmap="jet", shading="auto", vmin=-7, vmax=7)
    cbar = basemap.colorbar(location="bottom")
    cbar.set_label("v, m s$^{-1}$")
    ax3.set_title("v")

    ax4 = f.add_subplot(144)
    # set up orthographic map projection
    x, y, basemap = create_map_func(lons, lats)
    # contour data over the basemap.
    # cs = basemap.contour(x,y,wave+mean,15,linewidths=1.5)
    cs4 = basemap.pcolor(
        x, y, u[it, :, :], cmap="jet", shading="auto", vmin=-5, vmax=60
    )
    cbar = basemap.colorbar(location="bottom")
    cbar.set_label("u, m s$^{-1}$")
    ax4.set_title("u")

    ITER += 1
    plt.savefig(f"../../output/frames/frame{ITER:03d}.png", format="png", dpi=300)

os.system(
    "ffmpeg -r 5 -f image2  -i ../../output/frames/frame%03d.png -vframes 34 -vcodec libx264\
          -crf 25 -pix_fmt yuv420p ../../output/animations/animation.mp4"
)
os.system(
    "ffmpeg -i ../../output/animations/animation.mp4  ../../output/animations/animation.gif"
)
