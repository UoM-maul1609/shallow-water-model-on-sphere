import os
import getpass
import matplotlib
import sys
matplotlib.use('Agg')

import numpy as np
from netCDF4 import Dataset as NetCDFFile
import scipy as scy

def fourier_wave_number(fileName,u_jet):
    """
        uses findpeaks and fourier analysis to track the current wave number and
        the rotation speed (i.e. using phase information from the fft).
    """
    
    n_files=len(fileName); 

    nc=NetCDFFile(fileName[0]);
    np1=len(nc['time'][:]);
    nc.close()

    wave_number=np.zeros((n_files,np1));
    rotation_rate=np.zeros((n_files,np1-1));
    phase=np.zeros((n_files,np1));
    phaseold=np.zeros((n_files,np1));

    for j in range(n_files):
        nc=NetCDFFile(fileName[j]);
        dt_sec=nc['time'][1]-nc['time'][0];
        r,c = np.shape(nc['vort'][0,:,:]);
    
        Fs=c;
        T=1./Fs;
        L=c;
        t = np.mgrid[0:L]*T;
        f = Fs*np.mgrid[0:(L/2)+1]/L;

        for i in range(np1):
            X=np.mean(nc['v'][i,50:60+1,:],axis=0);
            Y = np.fft.fft(X);
            P2 = np.abs(Y/L);
            P1 = P2[0:int(L/2)+1];
            P1[1:-1] = 2*P1[1:-1];



            locs1=scy.signal.find_peaks(np.concatenate([X,X[0:3]]),prominence=1);
            locs=locs1[0]
            if (len(locs)>0):
                if (locs[0] == locs[-1]-c):
                    locs=np.delete(locs,len(locs)-1)

            ind=len(locs);
        
            phs=np.unwrap(np.angle(Y))*180./np.pi;

            print('wave number is : ' + str(f[ind]) + '; phase: ' + str(phs[ind]));
            phase[j,i]=phs[ind];
            wave_number[j,i]=f[ind];


        nc.close();
        phaseold[j,:]=phase[j,:];
        phase[j,:]=np.unwrap(phaseold[j,:]*np.pi/180.)*180./np.pi/(wave_number[j,:]);
        rotation_rate[j,:]=-np.diff(phase[j,:],n=1,axis=0) /  \
            dt_sec*86400.*365.25/360; # full rotations per year

    return (phase,rotation_rate,wave_number)

if __name__=='__main__':

    fileName=['/tmp/output.nc']
    n_files=len(fileName); 
    
    u_jets=[50.]

    (phase,rotation_rate,wave_number)=fourier_wave_number(fileName,u_jets)

    cmap_lev=64;
    map=plt.get_cmap(lut=cmap_lev,name='ocean')
    u1=np.linspace(u_jets[0],u_jets[-1],cmap_lev);
    int1=scy.interpolate.interp1d(u1,np.mgrid[1:cmap_lev+1],kind='nearest');


    mean_rot=np.zeros((n_files,9));
    std_rot=np.zeros((n_files,9));
    for j in range(n_files):
        for i in range(0,9):
            ind1,=np.where(np.diff(wave_number[j,:])==0)
            ind,=np.where(wave_number[j,ind1]==(i+1));
            ind = ind1[ind]
            if(len(ind)>2):
                mean_rot[j,i]=np.mean(rotation_rate[j,ind[1:-2]]);
#                 print(rotation_rate[j,ind[1:-2]])
                std_rot[j,i]=np.std(rotation_rate[j,ind[1:-2]]);
            else:
                mean_rot[j,i]=np.nan
                std_rot[j,i]=np.nan
                
        h=plt.scatter(np.mgrid[1:10],mean_rot[j,:],s=40,c=u_jets[j]*np.ones(9));
    plt.xlabel('wave number (per full rotation)');
    plt.ylabel('mean rotation rate (rotations per year)');
    h=plt.colorbar()
    h.set_label('Jet speed (m/s)')    

    plt.savefig('/tmp/fourier_wave_number.png' ,format='png', dpi=300) 
