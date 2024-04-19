import os
import getpass
import matplotlib
import sys
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from netCDF4 import Dataset as NetCDFFile

from scipy.interpolate import griddata

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

Re=5.4155760e7
(lo,la)=np.meshgrid(lons,np.pi/2.-lats)
arc=la*Re
x=arc*np.cos(lo)
y=arc*np.sin(lo)

#plt.ion()
hmap=plt.pcolor(x,y,h[-1,:,:])
plt.axis('square')
xx=np.linspace(-4.5e7,4.5e7,100)
yy=np.linspace(-4.5e7,4.5e7,100)

(xx1,yy1)=np.meshgrid(xx,yy)

#uu=griddata((x[0::10,0::10],y[0::10,0::10]),u[-1,0::10,0::10]-np.mean(u[-1,0::10,0::10],axis=0),(xx,yy))
#vv=griddata((x[0::10,0::10],y[0::10,0::10]),v[-1,0::10,0::10],(xx,yy))

u11=(u[-1,:,:]-np.mean(u[-1,:,:],axis=0))
v11=(v[-1,:,:])

# pol2cart velocity https://math.stackexchange.com/questions/2444965/relationship-between-cartesian-velocity-and-polar-velocity
u111=-v11*np.cos(lo)-u11*np.sin(lo)
v111=-v11*np.sin(lo)+u11*np.cos(lo)

u111=u111[0::5,0::5].flatten()
v111=v111[0::5,0::5].flatten()

uu=griddata((x[0::5,0::5].flatten(),y[0::5,0::5].flatten()),u111,(xx1,yy1))

vv=griddata((x[0::5,0::5].flatten(),y[0::5,0::5].flatten()),v111,(xx1,yy1))

plt.streamplot(xx,yy,uu,vv,4)
plt.xlim((-2e7,2e7))
plt.ylim((-2e7,2e7))

cbar=plt.colorbar(hmap)
cbar.set_label('h,m')

plt.title('Height at t=' + '{0:.2f}'.format(time[-1]/(86400)) +' days')

if not os.path.exists('../../output/frames'):
    os.mkdir('../../output/frames')

plt.savefig('../../output/frames/frame.png' ,format='png', dpi=300) 