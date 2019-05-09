1;

if ( ~exist('looe1','var') || isempty(looe1) )
  load('../data/looe1_ms.mat');
  if (isfield(looe1,'lon')); looe1 = rmfield(looe1,'lon'); end;
  if (isfield(looe1,'lat')); looe1 = rmfield(looe1,'lat'); end;
  if (isfield(looe1,'depth')); looe1 = rmfield(looe1,'depth'); end;
end;

%q0 = 'ndbc_ncep_30a_net_heat_flux';
%q0 = 'gom_hycom_netqf_heat_flux';
q0 = 'gom_hycom_dt_heat_flux';
dayq0 = [q0 '_24_hour_average'];
looe1 = verify_variable(looe1,dayq0);


if ( ~isfield(looe1,'microcat_seatemp') )
  looe1 = get_looe1_microcat(looe1);
end;

if ( ~isfield(looe1,'adcp_u') )
  looe1 = get_looe1_adcp(looe1,true);
end;

[ix1,ix2]=intersect_dates(looe1.microcat_seatemp.date,looe1.adcp_seatemp.date);
figure; maxigraph; plot(looe1.microcat_seatemp.date(ix1),looe1.microcat_seatemp.data(ix1)-looe1.adcp_seatemp.data(ix2)); datetick3; titlename('\DeltaT vs. ADCP w'); hold on;
plot(looe1.adcp_w_3hlp.date,looe1.adcp_w_3hlp.data*100,'r:');
legend('T_5_m - T_2_0_m [K]','Mean w (3HLP) [cm/s]', 'Location','Best');


btmbin = 2;
sfcbin = 19;

[looe1,hls,hxs]=multiplot_station(looe1,{'ndbc_air_t','microcat_seatemp','adcp_u',q0,'ndbc_wind1_u'},[],[],{'Air Temp','Sea Temps','Offshore Flow','Heat Flux','Winds'},[],{[5,33],[20,33],[-0.5,0.5],[-1000,1000],[-20,20]});
axes(hxs(2)); hold on; plot(looe1.adcp_seatemp.date,looe1.adcp_seatemp.data,'b.');
axes(hxs(3)); hold on; plot(looe1.adcp_u.date,looe1.adcp_u.prof(:,sfcbin),'g.');
axes(hxs(3)); hold on; plot(looe1.adcp_u.date,looe1.adcp_u.prof(:,btmbin),'b.');
axes(hxs(4)); hold on; plot(looe1.(dayq0).date,looe1.(dayq0).data,'b.');
axes(hxs(5)); hold on; plot(looe1.ndbc_wind1_v.date,looe1.ndbc_wind1_v.data,'k.');
xlim([datenum(2007,02,24),datenum(2008,02,20)]); datetick3('x',2,'keeplimits');
% xlim([datenum(2007,12,16),datenum(2007,12,27)]); datetick3('x',2,'keeplimits');
% xlim([datenum(2007,6,29),datenum(2007,7,17)]); datetick3('x',2,'keeplimits');
legend(hxs(2), 'z=5m','z=20m', 'Location','Best','Orientation','horizontal');
legend(hxs(3), 'Mean','Near-surface','Near-bottom', 'Location','Best','Orientation','horizontal');
legend(hxs(4), 'Hourly','24hr Mean', 'Location','Best','Orientation','horizontal');
legend(hxs(5), 'U','V', 'Location','Best','Orientation','horizontal');
clear hls hxs


%%%%%%%%%%%%%%%%%%
% SUMMER EVENTS
%%%%%%%%%%%%%%%%%%

% % figure;
% % evix=2800; dt=60*24; contourf(looe1.adcp_u_3hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_u_3hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Summer: u_3_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% % xlim([datenum(2007,6,29),datenum(2007,7,17)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;
% figure;
% evix=2800; dt=60*24; contourf(looe1.adcp_u_40hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_u_40hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Summer: u_4_0_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% xlim([datenum(2007,6,29),datenum(2007,7,17)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;

% % figure;
% % evix=2800; dt=60*24; contourf(looe1.adcp_v_3hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_v_3hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Summer: v_3_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% % xlim([datenum(2007,6,29),datenum(2007,7,17)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;
% figure;
% evix=2800; dt=60*24; contourf(looe1.adcp_v_40hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_v_40hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Summer: v_4_0_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% xlim([datenum(2007,6,29),datenum(2007,7,17)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;


% figure;
% evix=4177; dt=60*24; contourf(looe1.adcp_u_3hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_u_3hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Summer: u_3_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% xlim([datenum(2007,8,20),datenum(2007,9,9)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;
% % figure;
% % evix=4177; dt=60*24; contourf(looe1.adcp_u_40hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_u_40hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Summer: u_4_0_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% % xlim([datenum(2007,8,20),datenum(2007,9,9)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;

% figure;
% evix=4177; dt=60*24; contourf(looe1.adcp_v_3hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_v_3hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Summer: v_3_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% xlim([datenum(2007,8,20),datenum(2007,9,9)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;
% % figure;
% % evix=4177; dt=60*24; contourf(looe1.adcp_v_40hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_v_40hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Summer: v_4_0_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% % xlim([datenum(2007,8,20),datenum(2007,9,9)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;


%%%%%%%%%%%%%%%%%%
% WINTER EVENTS
%%%%%%%%%%%%%%%%%%

% figure;
% evix=6337; dt=60*24; contourf(looe1.adcp_u_3hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_u_3hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Winter: u_3_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% xlim([datenum(2007,12,16),datenum(2007,12,27)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;
figure;
evix=6337; dt=60*24; contourf(looe1.adcp_u_40hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_u_40hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Winter: u_4_0_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
xlim([datenum(2007,12,16),datenum(2007,12,27)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;

% figure;
% evix=6337; dt=60*24; contourf(looe1.adcp_v_3hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_v_3hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Winter: v_3_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% xlim([datenum(2007,12,16),datenum(2007,12,27)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;
figure;
evix=6337; dt=60*24; contourf(looe1.adcp_v_40hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_v_40hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Winter: v_4_0_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
xlim([datenum(2007,12,16),datenum(2007,12,27)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;

% figure;
% evix=6337; dt=60*24; contourf(looe1.adcp_baroclinic_u_3hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_baroclinic_u_3hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Winter: u^B^C_3_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% xlim([datenum(2007,12,16),datenum(2007,12,27)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;
figure;
evix=6337; dt=60*24; contourf(looe1.adcp_baroclinic_u_40hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_baroclinic_u_40hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Winter: u^B^C_4_0_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
xlim([datenum(2007,12,16),datenum(2007,12,27)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;

% figure;
% evix=6337; dt=60*24; contourf(looe1.adcp_baroclinic_v_3hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_baroclinic_v_3hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Winter: v^B^C_3_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
% xlim([datenum(2007,12,16),datenum(2007,12,27)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;
figure;
evix=6337; dt=60*24; contourf(looe1.adcp_baroclinic_v_40hlp.date(evix-dt:evix+dt),looe1.adcp_bin_heights(1:29),looe1.adcp_baroclinic_v_40hlp.prof(evix-dt:evix+dt,1:29)',[-2 -.2:.02:.2]); caxis([-.3 .3]); maxigraph; ylabel('Height above ADCP [m]'); titlename(['Winter: v^B^C_4_0_h_l_p(',num2str(evix),'\pm',num2str(dt),')']);
xlim([datenum(2007,12,16),datenum(2007,12,27)]); datetick3('x',2,'keeplimits'); colorbar; set_hovmuller_cursor;
