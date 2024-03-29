function normal_modes_compare(u_jets,flag)
% normal mode analysis of a gaussian jet - see Mak, Atmospheric Dynamics


cmap_lev=64;
map=jet(cmap_lev);
u1=linspace(u_jets(1),u_jets(end),cmap_lev);


L=1000e3;
U=10;
L_on_U=L/U;
UL=U*L;


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

    % derivative of vorticity wrt y (i.e. d/dy(-du/dy):
    zeta_y=2.*u./a.^2.*(1-2.*(y-b).^2./a.^2);

    % beta (df/dy)
    beta1=2.*2.*pi./T.*cos(y./re)./re;
    n_ks=36;
    sigma_max=zeros(n_ks,1);
    sigmas_max=zeros(n_ks,3);
    sigmas_max_i=zeros(n_ks,3);
    sigmas_min=zeros(n_ks,3);
    sigmas_min_i=zeros(n_ks,3);
    % now set-up matrix problem
    for n=1:n_ks
    %     n=5;
        k=2*pi*n/x_len;

        A=eye(jp,jp);
        B=eye(jp,jp);

        % A matrix:
        A(2:1+jp:jp^2)=i.*u(1:end-1).*k./dy^2; % top diag elements
        A(1:1+jp:jp^2)=i.*k.*((zeta_y+beta1)-2.*u./dy.^2-u.*k.^2); % diagonal elements
        A(1+jp:1+jp:jp^2)=i.*u(2:end).*k./dy^2; % bottom diag elements

        % B matrix:
        B(2:1+jp:jp^2)=-1./dy.^2; % top diag elements
        B(1:1+jp:jp^2)=(k.^2+2./dy.^2); % diagonal elements
        B(1+jp:1+jp:jp^2)=-1./dy^2; % bottom diag elements



        % find eigenvalues:
        [E,D]=eig(A,B);
        sigma1=diag(D);

    %     sigma1=real(sigma1).*max(E(:,:),[],2)+imag(sigma1);

        sigma_max(n)=max(real(sigma1));
        [a,ii]=sort(real(sigma1),'descend');
        sigmas_max(n,1:3)=a(1:3);
        sigmas_max_i(n,1:3)=imag(sigma1(ii(1:3)));

        [a,ii]=sort(real(sigma1),'ascend');
        sigmas_min(n,1:3)=a(1:3);
        sigmas_min_i(n,1:3)=imag(sigma1(ii(1:3)));
        clear i;

    end
    if flag==1
        h=plot(1:n_ks,...
            (-sigmas_max_i(:,1)./repmat(2.*pi.*[1:n_ks]'./x_len,[1 1])) ...
            .*86400.*365.25./x_len);
    elseif flag==2
        h=plot(1:n_ks,...
            sigmas_max(:,1));
%         h2=plot(1:n_ks,...
%             sum(sigmas_max(:,:),2),'--');
        hold on;
        h2=plot(1:n_ks,...
            sigmas_max(:,2),'--');
        hold on;
        h3=plot(1:n_ks,...
            sigmas_max(:,3),':');
        h=[h, h2, h3];
        hold on;
    elseif flag==3
        h=plot(1:n_ks,...
            (-sigmas_max_i(:,1)./(2.*pi.*[1:n_ks]'./x_len)));
        hold on;
        h2=plot(1:n_ks,...
            (-sigmas_max_i(:,2)./(2.*pi.*[1:n_ks]'./x_len)),'--');
        h=[h, h2];
        hold on;
    end
    
    if ( length(u_jets) > 1)
        row=interp1(u1,1:cmap_lev,u_jets(j),'nearest');
    else
        row=1;
    end
    set(h,'color',map(row,:));
    hold on;
%     ylim([0 30])
%     xlim([0 12])
end
% xlabel('wave numbers per 2\pi')
% ylabel('rotations per year')
% set(gca,'xtick',[0 5 6 8 10 12])
grid on;
%     title(['jet speed: ',num2str(u_jet)])


