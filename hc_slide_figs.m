1;

figspath = get_thesis_path('../figs');

if ( ~exist('stn','var') || ~isfield(stn,'ndbc_erai_erai_30a_net_flux') )
    stn = optimize_station_heat_budget('smkf1','erai','none','ndbc','tpxo_tide','erai');
end;

%multiplot_station(stn,{'ndbc_air_t','ndbc_sea_t','ndbc_erai_erai_30a_wind_stress_xshore','ndbc_erai_erai_30a_net_flux_1_d_sum'},'Sombrero Key Cold Front 2007',[],{'T_a','T_s','\tau^x^s,\tau^l^s','\Sigma_1_dQ_0'},datenum(2007,[11,12],[21,31]),{[12,32],[12,32],[-.5,+.5],[-15000,15000]},true);
[stn,hls,hxs,hf]=...
    multiplot_station(stn,{'ndbc_air_t','ndbc_sea_t','ndbc_erai_erai_30a_wind_stress_xshore','ndbc_erai_erai_30a_net_flux_term_1_d_sum'},...
    'Sombrero Key Cold Front 2007',[],{'Air Temp [^oC]','Sea Temp [^oC]','Wind Stress [N/m^2]','Heat Flux [K/day]'},...
    datenum(2007,[11,12],[21,31]),{[12,32],[12,32],[-.8,+.8],[-1.5,+1.5]},true,{'k-','k-','k-','k-'});
for hlix=1:numel(hls);
    set(hls{hlix},'LineWidth',2');
end;
axes(hxs(3)); hold on;
plot(stn.ndbc_erai_erai_30a_wind_stress_lshore.date,stn.ndbc_erai_erai_30a_wind_stress_lshore.data,'k--','Color',[.5,.5,.5],'LineWidth',3);
%set(hxs(3),'XTick',[]);
legh=legend('\tau^c^r^o^s^s','\tau^a^l^o^n^g','Location','NorthWest'); set(legh,'FontSize',9);
axes(hxs(4)); hold on;
stn=verify_variable(stn,'ndbc_erai_erai_30a_erai_none_hc_dTdt_1_d_sum');
plot(stn.ndbc_erai_erai_30a_erai_none_hc_dTdt_1_d_sum.date,stn.ndbc_erai_erai_30a_erai_none_hc_dTdt_1_d_sum.data,'k--','Color',[.5,.5,.5],'LineWidth',3);
legh=legend('Air-Sea Q_0/\rhoC_ph','Total Budget','Location','NorthWest'); set(legh,'FontSize',9);

print('-dtiff',fullfile(figspath,'smkf1-cold-front-2007-multiplot.tiff'));


if ( ~exist('looe1','var') || ~isfield(looe1,'adcp_baroclinic_x_40hlp') )
    looe1=get_station_from_station_name('looe1'); looe1=get_looe1_microcat(looe1); looe1=get_looe1_adcp(looe1);
end;
fmg;
[ig,dtix]=min(abs(datenum(2007,12,22)-looe1.adcp_baroclinic_x_40hlp.date));
[cs,ch]=...
    contourf(looe1.adcp_baroclinic_x_40hlp.date(dtix-1440:dtix+1440),...
             looe1.adcp_bin_heights,...
             looe1.adcp_baroclinic_x_40hlp.prof(dtix-1440:dtix+1440,:)',...
             [-0.25:0.025:0.25]);
axis([datenum(2007,[11,12],[21,31]),-2,20.88,-0.25,0.25,-0.25,0.25]);
cbh=colorbar; set(cbh,'FontSize',10);
clabel(cs,ch,'FontSize',9);
xlabel('Date','FontSize',12,'FontWeight','bold');
ylabel('Height above bottom [m]','FontSize',12,'FontWeight','bold');
plot_ts(ts_op(looe1.microcat_seatemp,looe1.adcp_seatemp,'-'));
text(datenum(2007,12,19,6,0,0),1.5,'5m - 20m Sea Temp [K]','FontSize',9,'FontWeight','bold');
text(datenum(2007,12,19,6,0,0),15.2,sprintf('Cross-shore BC Current\n  40hlp [m/s]'),'FontSize',9,'FontWeight','bold');
ylim([-2,20.88]);
plot_ts(stn.ndbc_erai_erai_30a_net_flux_term_1_d_sum,'r');
text(datenum(2007,12,19,6,0,0),-1.2,'Sombrero Q_0/\rhoC_ph [K/day]','FontSize',9,'FontWeight','bold');
datetick3; titlename('Looe Key ADCP Cold Front 2007');
set(datacursormode(gcf),'UpdateFcn',@tyzc_select_cb);

print('-dtiff',fullfile(figspath,'looe1-cold-front-2007.tiff'));

clear figspath hls hxs hf hlix legh ig dtix cbh cs ch ans
