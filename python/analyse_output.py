import os
import getpass
import matplotlib
import sys
matplotlib.use('Agg')

import numpy as np

username=getpass.getuser()

if not os.path.exists('/tmp/' + username):
    os.mkdir('/tmp/' + username)


""" do the fourier analysis
"""

u_jet=[5., 10., 20., 30., 40., 50., 60., 70., 80.,\
 90., 100., 125., 150., 175., 200., 250., 300., 350.,];

c_vis=[0.,0.1,0.2,1.0]

plt.ion()
plt.figure()
for j in range(len(c_vis)):
    fileNames=['/tmp/' + username + '/output' + str(i) + '_' + str(j) + '.nc' \
        for i in range(len(u_jet))]
    
    plt.subplot(1,len(c_vis),j)
    
    fourier_wave_number.do_analysis01(fileNames,u_jet[i]);

    # same for each value of c_cvis
    normal_modes_compare(u_jet[i],1);
    plt.title('$C_{vis}$=' + str(c_vis[j]))
    
plt.savefig('/tmp/' + username + '/full_analysis.png')


