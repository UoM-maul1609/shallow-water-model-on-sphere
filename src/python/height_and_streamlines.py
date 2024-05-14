"""
Module Name: height_and_stREamlines.py

Description:
This module contains code to analyse the height and
streamlines and produce a frame.png image output at
the final timestep of the data.
"""

import os
import getpass
import matplotlib

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

from scipy.interpolate import griddata
# pylint: disable=E0611
# (pylint can't find Dataset in the netCDF4 package for some reason)
from netCDF4 import Dataset as NetCDFFile

matplotlib.use('Agg')

username=getpass.getuser()

print(os.listdir())

nc=NetCDFFile('../../tests/output.nc')
lons=nc.variables['phi'][:]
lats=nc.variables['theta'][:]
vort=nc.variables['vort'][:]
h=nc.variables['h'][:]
u=nc.variables['u'][:]
v=nc.variables['v'][:]
time=nc.variables['time'][:]

RE=5.4155760e7
(lo,la)=np.meshgrid(lons,np.pi/2.-lats)
arc=la*RE
x=arc*np.cos(lo)
y=arc*np.sin(lo)

hmap=plt.pcolor(x,y,h[-1,:,:])
plt.axis('square')
xx=np.linspace(-4.5e7,4.5e7,100)
yy=np.linspace(-4.5e7,4.5e7,100)

(xx1,yy1)=np.meshgrid(xx,yy)

# uu=griddata((x[0::10,0::10],y[0::10,0::10]),u[-1,0::10,0::10]-np.mean(u[-1,0::10,0::10],axis=0),(xx,yy)) pylint: disable=C0301
# vv=griddata((x[0::10,0::10],y[0::10,0::10]),v[-1,0::10,0::10],(xx,yy))

u11=(u[-1,:,:]-np.mean(u[-1,:,:],axis=0))
v11=(v[-1,:,:])

# pol2cart velocity
# https://math.stackexchange.com/questions/2444965/relationship-between-cartesian-velocity-and-polar-velocity pylint: disable=C0301
u111=-v11*np.cos(lo)-u11*np.sin(lo)
v111=-v11*np.sin(lo)+u11*np.cos(lo)

u111=u111[0::5,0::5].flatten()
v111=v111[0::5,0::5].flatten()

uu=griddata((x[0::5,0::5].flatten(),y[0::5,0::5].flatten()),u111,(xx1,yy1))
vv=griddata((x[0::5,0::5].flatten(),y[0::5,0::5].flatten()),v111,(xx1,yy1))

plt.streamplot(xx,yy,uu,vv,4)
plt.xlim((-2e7,2e7))
plt.ylim((-2e7,2e7))

# Adding the labels
max_h_magnitude = np.max(np.abs(h))
scale = 10 ** np.floor(np.log10(max_h_magnitude))
scaled_label = f'h (10$^{int(np.log10(scale))}$ m)'
cbar=plt.colorbar(hmap)
cbar.set_label(scaled_label)

# Formatter so the colourbar has the correct notation
formatter = FuncFormatter(lambda x, _: f'{x/scale:.2f}')
cbar.ax.yaxis.set_major_formatter(formatter)

# Changing the notation of the x and y axis
plt.ticklabel_format(style='sci', scilimits=(0,0))
plt.gca().ticklabel_format(useMathText=True)

plt.title(f"Height and Velocity Streamlines at t={time[-1]/(86400):.2f} days")

if not os.path.exists('../../output/frames'):
    os.mkdir('../../output/frames')

plt.savefig('../../output/frames/frame.png' ,format='png', dpi=300)
