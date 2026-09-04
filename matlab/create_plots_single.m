% create some plots / animations
figure('name','jets');
print_image=true;
nc=netcdf(['/tmp/output.nc']);
phi = nc{'phi'}(:);
theta = nc{'theta'}(:);
for n=1:135 % time loop
    m_proj('stereographic','lat',90,'radius',25,'rotangle',45);
    hfield = squeeze(nc{'h'}(n,:,:));
    [phi_plot,h_plot] = cyclic_lon(phi,hfield);
    m_pcolor(phi_plot.*180./pi-180,theta.*180./pi,h_plot);shading flat
    m_grid('fontsize',6,'xticklabels',[],'xtick',[],'ytick',[],'yticklabels',[]);
    title(['time (earth days)=',num2str(nc{'time'}(n)./86400,'%.2f')]);

    if print_image
        if(n==1)
            mkdir /tmp/pics/
        end
        eval(['print -dpng /tmp/pics/output_',num2str(n,'%03d'),'.png']);
    end
end
close(nc);
