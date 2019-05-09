1;

% Should first call CALCQ0, CALCDT, TRYHC, TRYGOM, TRYFLKEYS, or similar

[yr,mo,dy] = datevec(stn.ndbc_ncep_30a_net_heat_flux.date);
yrmo = datenum(yr,mo,1);
mo30a = grpstats(real(stn.ndbc_ncep_30a_net_heat_flux.data),yrmo);

[yr,mo,dy] = datevec(stn.ndbc_hfbulk_net_heat_flux.date);
hfyrmo = datenum(yr,mo,1);
mohf = grpstats(real(stn.ndbc_hfbulk_net_heat_flux.data),hfyrmo);

[yr,mo,dy] = datevec(stn.ncep_net_heat_flux.date);
ncepyrmo = datenum(yr,mo,1);
moncep = grpstats(real(stn.ncep_net_heat_flux.data),ncepyrmo);

[yr,mo,dy] = datevec(stn.gom_hycom_net_heat_flux.date);
gomyrmo = datenum(yr,mo,1);
mogom = grpstats(real(stn.gom_hycom_net_heat_flux.data),gomyrmo);

figure;
maxigraph;
hold on;
plot(unique(yrmo),mo30a,'m-');
plot(unique(hfyrmo),mohf,'kd');
plot(unique(ncepyrmo),moncep,'cs');
plot(unique(gomyrmo),mogom,'r--');
plot(stn.landy_net_heat_flux.date,stn.landy_net_heat_flux.data,'bo:');
plot(stn.monthly_nocs_net_heat_flux.date,stn.monthly_nocs_net_heat_flux.data,'g^--');
ylabel('Mean Monthly Net Heat Flux [ W / m^2 ]'); 
xlim([datenum(2002,12,1),datenum(2007,2,1)]);
datetick3;
% set(gca,'color',[.95 .95 .95]);
title(['Estimates of Q_0 at ' stn.station_name]);
legend('Gramer & Mariano','Smith (1988) bulk formulae','NCEP NARR 32km reanalysis','GoM 4km HYCOM + NCODA',...
    'Large & Yeager (2009)','NOC Southampton (2009)', 'Location','Best');
% print('-dtiff','../figs/mlrf1-cmp_monthly_fluxes.tiff');
