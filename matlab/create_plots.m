% create some plots / animations
u_jet=[5. 10. 20. 30. 40. 50. 60. 70. 80. 90. 100. 125. 150. 175. 200. 250 300 350];
cvis=[0. 0.1 0.2 1.];
% figure('name','jets');
mkdir /tmp/pics/
for ii=1:6 % 6 plots for different jet strengths
    for n=1:135 % time loop
        k=1;
        for j=1:4 % loop over cvis
            for i=1:3 % loop over jet strength
                iii=i-1+(ii-1)*3; % indexes the strength of the jet
                nc=netcdf(['output_',num2str(iii,'%d'),'_',num2str(j-1,'%d'),'.nc']);

                subplot(4,3,k);
                if (ii==1 )
                    m_proj('stereographic','lat',90,'radius',25,'rotangle',45);
                end
%                 m_contourf(nc{'phi'}(:).*180./pi-180,nc{'theta'}(:).*180./pi,nc{'h'}(n,:,:),15);shading flat
                m_pcolor(nc{'phi'}(:).*180./pi-180,nc{'theta'}(:).*180./pi,nc{'h'}(n,:,:));shading flat
                m_grid('fontsize',6,'xticklabels',[],'xtick',[],'ytick',[],'yticklabels',[]);  
                h=get(gca,'children');
                set(h(end),'facealpha',0)
                if(k <= 3 ) 
                    title(['u_{max}=',num2str(u_jet(iii+1)),' m/s']);
                end
                if(mod(k-1,3)==0)
                    ylabel(['c_{vis}=',num2str(cvis(j))]);                    
                end
                if(k==11 && ii==5)
                    xlabel(['time (earth days)=',num2str(nc{'time'}(n)./86400,'%.2f')]);
                end
                close(nc);
                k=k+1;
            end
        end
        eval(['print -dpng /tmp/pics/output_',num2str(ii-1,'%02d'),'_', ...
            num2str(n,'%03d'),'.png']);
    end
end
    