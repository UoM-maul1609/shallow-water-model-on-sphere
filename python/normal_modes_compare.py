import os
import getpass
import matplotlib
import sys
matplotlib.use('Agg')

import numpy as np
import scipy as scy
import scipy.interpolate as sci
import matplotlib.pyplot as plt

username=getpass.getuser()

if not os.path.exists('/tmp/' + username):
    os.mkdir('/tmp/' + username)


def normal_modes_compare(u_jets,flag):
    """
        normal mode analysis of a gaussian jet - see Mak, Atmospheric Dynamics, page 229
    """

    cmap_lev=64;
    map=plt.get_cmap(lut=cmap_lev,name='viridis')
    u1=np.linspace(u_jets[0],u_jets[-1],cmap_lev);
    int1=sci.interp1d(u1,np.mgrid[1:cmap_lev+1],kind='nearest');

    
    L=1000e3;       # measure distance in units of 1000 km
    U=10;           # measure speed in units of 10 m/s
    L_on_U=L/U;
    UL=U*L;



    T=10.55*3600;   # time-period of rotation

    # set up domain.
    ip=200;         # finite difference grid
    jp=100;
    lat_high=85;    # highest latitude
    lat_low=65;     # lowest
    re=5.4155760e7; # radius of saturn in this region (due to squashed spheriod)
    #re=5.8232e7;
    h_jet=1.;      # standard deviation of jet
    lat_jet=77;     # latitude of the jet
    n_ks=36;        # calculate the growth factor for this many k-values



    for i in range(0,len(u_jets)):
        u_jet=u_jets[i];


        # y-grid
        y=np.linspace(re*lat_low*np.pi/180.,re*lat_high*np.pi/180,jp);
        # distance around planet at this latitude
        x_len=2.*np.pi*np.cos((lat_jet)*np.pi/180.)*re;
        x_len2=2.*np.pi*np.cos((lat_jet)*np.pi/180.)*re*4/7;

        # x-grid
        x=np.linspace(0,x_len,ip);

        X,Y=np.meshgrid(x,y);
        dy=y[1]-y[0];
        # the jet:
        a=np.sqrt(2)*h_jet*np.pi/180.*re;
        b=lat_jet*np.pi/180.*re;
        u=u_jet*np.exp(-((y-b)/a)**2);

        # derivative of vorticity wrt y (i.e. d/dy(-du/dy):
        zeta_y=2.*u/a**2.*(1.-2.*(y-b)**2./a**2);

        # beta (df/dy)
        beta1=2.*2.*np.pi/T*np.cos(y/re)/re;
        sigma_max=np.zeros(n_ks);
        sigmas_max=np.zeros((n_ks,3));
        sigmas_max_i=np.zeros((n_ks,3));
        sigmas_min=np.zeros((n_ks,3));
        sigmas_min_i=np.zeros((n_ks,3));
        # now set-up matrix problem
        for n in range(1,n_ks+1):
            k=2*np.pi*n/x_len #*(x_len/x_len2);


            # A matrix:
            # top diag elements + diagonal elements + bottom diag elements
            A=np.diag(1j*u[0:-1]*k/dy**2,1) + \
                np.diag(1j*k*((zeta_y+beta1)-2.*u/dy**2-u*k**2)) + \
                np.diag(1j*u[1:]*k/dy**2,-1) 

            # B matrix:
            # top diag elements + diagonal elements + bottom diag elements
            B=np.diag(-1./dy**2*np.ones(len(u)-1),1) + \
                np.diag((k**2+2./dy**2)*np.ones(len(u))) + \
                np.diag(-1./dy**2*np.ones(len(u)-1),-1); 

            # find eigenvalues:
            E=scy.linalg.eig(A,B,left=True);
#             E,D = scy.linalg.eig(A,B)
            sigma1=E[0];
#             sigma1=np.diag(D)

 
            sigma_max[n-1]=np.max(np.real(sigma1));
            a=np.flip(np.sort(np.real(sigma1)));
            ii=np.flip(np.argsort(np.real(sigma1)));
            sigmas_max[n-1,0:3]=np.real(sigma1[ii[0:3]]);
            sigmas_max_i[n-1,0:3]=np.imag(sigma1[ii[0:3]]);



            a=np.sort(np.real(sigma1));
            ii=np.argsort(np.real(sigma1));
            sigmas_min[n-1,0:3]=np.real(sigma1[ii[0:3]]);
            sigmas_min_i[n-1,0:3]=np.imag(sigma1[ii[0:3]]);
            
            # if growth rate small, set to NaN
            if sigmas_max[n-1,0]<1.e-9:
                sigmas_max_i[n-1,0]=np.nan
                
            if sigmas_max[n-1,1]<1.e-9:
                sigmas_max_i[n-1,1]=np.nan
                
            if sigmas_max[n-1,2]<1.e-9:
                sigmas_max_i[n-1,2]=np.nan

        if flag==1:
            h=[]
            # sigma / k is wave speed in m/s
            # k is 2*pi/lambda and lambda is x_len/n
            # so this is rotations per year
            h.append(plt.plot(np.mgrid[1:n_ks+1], \
                (-sigmas_max_i[:,0]/(2.*np.pi*np.mgrid[1:n_ks+1]/x_len)) \
                *86400.*365.25/x_len));
            h.append(plt.plot(np.mgrid[1:n_ks+1], \
                (-sigmas_max_i[:,1]/(2.*np.pi*np.mgrid[1:n_ks+1]/x_len)) \
                *86400.*365.25/x_len,'--'));
            h.append(plt.plot(np.mgrid[1:n_ks+1], \
                (-sigmas_max_i[:,2]/(2.*np.pi*np.mgrid[1:n_ks+1]/x_len)) \
                *86400.*365.25/x_len,':'));
            plt.ylabel('speed of wave (rotations per year)')
            plt.xlabel('number of peaks')
        elif flag==2:
            h=plt.plot(np.mgrid[1:n_ks+1],sigmas_max[:,0]);
            h2=plt.plot(np.mgrid[1:n_ks+1],sigmas_max[:,1],'--');
            h3=plt.plot(np.mgrid[1:n_ks+1],sigmas_max[:,2],':');
            h4=plt.plot(np.mgrid[1:n_ks+1],np.sum(sigmas_max[:,:],axis=1),lw=3);
            h=[h[0], h2[0], h3[0] ,h4[0]];
            plt.xlabel('number of peaks')
            plt.ylabel('Growth rate')
            plt.legend(['Largest','2nd largest','3rd largest','sum'])
        elif flag==3:
            h=plt.plot(np.mgrid[1:n_ks+1], \
                (-sigmas_max_i[:,0]/(2.*np.pi*np.mgrid[1:n_ks+1]/x_len)));
            h2=plt.plot(np.mgrid[1:n_ks+1], \
                (-sigmas_max_i[:,1]/(2.*np.pi*np.mgrid[1:n_ks+1]/x_len)),'--');
            h=[h, h2];

        if ( len(u_jets) > 1):
            row=int(int1(u_jets[i]))
        else:
            row=1;

        for k in range(len(h)):
            h[k][0].set_color(map(row));

            
        if flag==2:
            plt.legend(['Largest','2nd largest','3rd largest','sum'])

    
if __name__=='__main__':      
    plt.ion()
    plt.figure() 
    u_jets=[50.,70,100.]
#     u_jets=[50.]
    normal_modes_compare(u_jets,1)

    plt.savefig('/tmp/' + username +'/fourier_wave_number.png' ,format='png', dpi=300) 
    