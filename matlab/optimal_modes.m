function normal_modes_compare(u_jets,flag)
% normal mode analysis of a gaussian jet - see Mak, Atmospheric Dynamics


cmap_lev=64;
map=jet(cmap_lev);
u1=linspace(u_jets(1),u_jets(end),cmap_lev);



T=10.55*3600;

% set up domain.
ip=200;
jp=100;
lat_high=85;
lat_low=65;
re=5.4155760e7;
%re=5.8232e7;
h_jet=1.5;
lat_jet=78;
tau=1e5;


for j=1:length(u_jets)
    u_jet=u_jets(j);


    % y-grid
    y=linspace(re*lat_low*pi./180,re*lat_high*pi./180,jp);
    x_len=2.*pi.*cos((lat_jet).*pi./180).*re;

    x=linspace(0,x_len,ip);

    [X,Y]=meshgrid(x,y);
    dy=y(2)-y(1);
    % the jet:
    a=sqrt(2)*h_jet.*pi./180*re;
    b=lat_jet.*pi./180.*re;
    u=u_jet.*exp(-((y-b)./a).^2);

    psi_bar=-cumsum(u.*dy);
    % derivative of vorticity wrt y (i.e. d/dy(-du/dy):
    zeta_y=2.*u./a.^2.*(1-2.*(y-b).^2./a.^2);

    % beta (df/dy)
    beta1=2.*2.*pi./T.*cos(y./re)./re;
    n_ks=1;
    
    
    % now set-up matrix problem
    for n=6
        k=2*pi*n/x_len;

        A=eye(jp,jp);
        B=eye(jp,jp);

        % A matrix:
        A(2:1+jp:jp^2)=1i.*u(1:end-1).*k./dy^2; % top diag elements
        A(1:1+jp:jp^2)=1i.*k.*((zeta_y+beta1)-2.*u./dy.^2-u.*k.^2); % diagonal elements
        A(1+jp:1+jp:jp^2)=1i.*u(2:end).*k./dy^2; % bottom diag elements

        % B matrix:
        B(2:1+jp:jp^2)=-1./dy.^2; % top diag elements
        B(1:1+jp:jp^2)=(k.^2+2./dy.^2); % diagonal elements
        B(1+jp:1+jp:jp^2)=-1./dy^2; % bottom diag elements



        % find eigenvalues:
        [E,D]=eig(A,B);
        sigma1=diag(D);

        
        %--------

        % so now we have everything we need for 
        %phi eigenvalues (P-matrix), and we have lambda-matrix
        % Lambda-matrix
        lambda1=zeros([jp,jp]);
        for jj=1:jp
            lambda1(jj,jj)=exp(1i*sigma1(jj)*tau);
        end
        
        % P-matrix
        P=E;
        
        % now B tau
        Btau=ctranspose(lambda1)*ctranspose(P)*eye(jp,jp)*P*lambda1;
        % B0
        B0=ctranspose(P)*eye(jp,jp)*P;
        
        % find solution for optimal growth:
        [a,lam]=eig(Btau,B0);
%         lam1=diag(lam);
        
        
        % amplification factor
        l=ctranspose(a)*ctranspose(lambda1)*ctranspose(P)* ...
            eye(jp,jp)*P*lambda1*a/(ctranspose(a)*ctranspose(P)*eye(jp,jp)*P*a);

        l1=diag(lam);
        ind=find(max(real(l1))==real(l1));
        
        % Calcualte the stream function for optimal mode:
        psi=zeros(jp,ip);
        for jj=1:jp
            
            for ii=1:ip
                psi(:,ii)=psi(:,ii)+...
                    real(a(ind,jj).*P(:,jj)).*cos(k*x(ii)) - ...
                    imag(a(ind,jj).*P(:,jj)).*sin(k*x(ii));
            end
            
        end
        
        
        % Calcualte the stream function for one normal mode:
        psi1=zeros(jp,ip);
        for jj=1:jp
            
            for ii=1:ip
                psi1(:,ii)=psi(:,ii)+...
                    real(E(:,jj)).*cos(k*x(ii)) - ...
                    imag(E(:,jj)).*sin(k*x(ii));
            end
            
        end
        
        
        
        
    end

end

u=-diff(psi+repmat(psi_bar',[1 ip]),1,1)./dy;
dx=x(2)-x(1);
v=diff(psi+repmat(psi_bar',[1 ip]),1,2)./dx;

subplot(211);
pcolor(x,y./re.*180./pi,real((psi)));shading flat

subplot(212)
pcolor(v);shading flat
