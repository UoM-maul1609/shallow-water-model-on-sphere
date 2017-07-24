% create some plots / animations
figure('name','jets');
print_image=true;
nc=netcdf(['/tmp/output.nc']);
for n=1:135 % time loop


    m_proj('stereographic','lat',90,'radius',25,'rotangle',45);
    m_pcolor(nc{'phi'}(:).*180./pi-180,nc{'theta'}(:).*180./pi,nc{'h'}(n,:,:));shading flat
    m_grid('fontsize',6,'xticklabels',[],'xtick',[],'ytick',[],'yticklabels',[]);  
%     h=get(gca,'children');
%     set(h(end),'facealpha',0)
    title(['time (earth days)=',num2str(nc{'time'}(n)./86400,'%.2f')]);
    
    if print_image
        if(n==1)
            mkdir /tmp/pics/
        end
        eval(['print -dpng /tmp/pics/output_', ...
            num2str(n,'%03d'),'.png']);
    end
end
    
close(nc);
