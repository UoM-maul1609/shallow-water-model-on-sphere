import os
import getpass
import matplotlib
import sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import fourier_wave_number
import normal_modes_compare

username=getpass.getuser()

if not os.path.exists('/tmp/' + username):
    os.mkdir('/tmp/' + username)


""" do the fourier analysis
"""

u_jet=[50., 70., 100.];

c_vis=[0.2]

plt.ion()
plt.figure()
for j in range(len(c_vis)):
    fileNames=['/tmp/' + username + '/output_' + str(i) + '_' + str(j) + '.nc' \
        for i in range(len(u_jet))]
    
    plt.subplot(1,len(c_vis),j+1)
    
    fourier_wave_number.do_analysis01(fileNames,u_jet);
    
    # same for each value of c_cvis
    normal_modes_compare.normal_modes_compare(u_jet,1);
    plt.title('$C_{vis}$=' + str(c_vis[j]))
    
plt.savefig('/tmp/' + username + '/full_analysis.png')


