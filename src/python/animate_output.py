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
from matplotlib import ticker

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
    f = plt.figure()

    ax1 = f.add_subplot(141)
    # set up orthographic map projection
    x, y, basemap = create_map_func(lons, lats)
    # contour data over the basemap
    cs1 = basemap.pcolor(
        x, y, h[it, :, :], cmap="jet", shading="auto", vmin=60000, vmax=64000
    )
    cbar = basemap.colorbar(location="bottom")
    cbar.ax.set_xticks(cbar.ax.get_xticks())
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation="horizontal", fontsize=8)
    # Format colorbar tick labels using scientific notation
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 3))  # Adjust the limits as needed
    cbar.ax.xaxis.set_major_formatter(formatter)
    power_text = cbar.ax.xaxis.get_offset_text()
    power_text.set_size(8)
    # Set label and title
    cbar.set_label("h, m", fontsize=8)
    ax1.set_title("Height", fontsize=8)

    ax2 = f.add_subplot(142)
    # set up orthographic map projection
    x, y, basemap = create_map_func(lons, lats)
    # contour data over the basemap.
    cs2 = basemap.pcolor(
        x, y, vort[it, :, :], cmap="jet", shading="auto", vmin=-0.00005, vmax=0.00005
    )
    cbar = basemap.colorbar(location="bottom")
    cbar.ax.set_xticks(cbar.ax.get_xticks())
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation="horizontal", fontsize=8)
    # Format colorbar tick labels using scientific notation
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 3))  # Adjust the limits as needed
    cbar.ax.xaxis.set_major_formatter(formatter)
    power_text = cbar.ax.xaxis.get_offset_text()
    power_text.set_size(8)
    # Set label and title
    cbar.set_label("$\\zeta$, s$^{-1}$", fontsize=8)
    ax2.set_title("Vorticity", fontsize=8)

    ax3 = f.add_subplot(143)
    # set up orthographic map projection
    x, y, basemap = create_map_func(lons, lats)
    # contour data over the basemap.
    cs3 = basemap.pcolor(x, y, v[it, :, :], cmap="jet", shading="auto", vmin=-7, vmax=7)
    cbar = basemap.colorbar(location="bottom")
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation="horizontal", fontsize=8)
    cbar.set_label("v, m s$^{-1}$", fontsize=8)
    ax3.set_title("v", fontsize=8)

    ax4 = f.add_subplot(144)
    # set up orthographic map projection
    x, y, basemap = create_map_func(lons, lats)
    # contour data over the basemap
    cs4 = basemap.pcolor(
        x, y, u[it, :, :], cmap="jet", shading="auto", vmin=-5, vmax=60
    )
    cbar = basemap.colorbar(location="bottom")
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation="horizontal", fontsize=8)
    cbar.set_label("u, m s$^{-1}$", fontsize=8)
    ax4.set_title("u", fontsize=8)

    ITER += 1
    plt.suptitle(f"t={time[it]/86400:.2f} days", fontsize=8, y=0.25)
    plt.savefig(f"../../output/frames/frame{ITER:03d}.png", format="png", dpi=300)

os.system(
    "ffmpeg -r 5 -f image2  -i ../../output/frames/frame%03d.png -vframes 34 -vcodec libx264\
          -crf 25 -pix_fmt yuv420p ../../output/animations/animation.mp4"
)
os.system(
    "ffmpeg -i ../../output/animations/animation.mp4  ../../output/animations/animation.gif"
)
