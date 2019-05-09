1;

% stn = tryhc(stn,false,[],1.20,0.30);

% scatter_hc(stn,'ndbc_sea_t_qc_12_hour_lowpass','netqf',@(x)(1:100:length(x.date)));
% scatter_hc(stn,'ndbc_sea_t_qc_12_hour_lowpass','ndbc_ncep_30a_total_heat_flux_term',@(x)(1:100:length(x.date)));
% help cov
% 0.1*sw_cp(36,25,2)*sw_dens(36,25,2)/3600
% help robustfit
% help bootstrp
% help ksdensity
% m = bootstrp(100,@mean,stn.ndbc_sea_t.data); [fi,xi]=ksdensity(m); figure; plot(xi,fi);
% m = bootstrp(100,@median,stn.ndbc_sea_t.data); [fi,xi]=ksdensity(m); figure; plot(xi,fi);
% m = bootstrp(100,@std,stn.ndbc_sea_t.data); [fi,xi]=ksdensity(m); figure; plot(xi,fi);
% m = bootstrp(100,@var,stn.ndbc_sea_t.data); [fi,xi]=ksdensity(m); figure; plot(xi,fi);
% m=bootstrp(100,@var,stn.ndbc_sea_t.data(ts_boreal_warm(stn.ndbc_sea_t))); [fi,xi]=ksdensity(m); figure; plot(xi,fi);
% m=bootstrp(100,@var,stn.ndbc_sea_t.data(ts_boreal_cool(stn.ndbc_sea_t))); [fi,xi]=ksdensity(m); figure; plot(xi,fi);


Y = lanczosfilter(stn.ndbc_sea_t.data,1,1/20);
stn = verify_variable(stn,'ndbc_sea_t_20_hour_lowpass');
maxigraph(figure); d=stn.ndbc_sea_t.date;
plot(d,stn.ndbc_sea_t.data,'r-',...
     d,Y,'b--',...
     stn.ndbc_sea_t_20_hour_lowpass.date,stn.ndbc_sea_t_20_hour_lowpass.data,'y:');
datetick3; legend('T_S','Lanczos(T_S,20hr)','20hLP(T_S)');
set(gca,'color',[.5 .5 .5]);
