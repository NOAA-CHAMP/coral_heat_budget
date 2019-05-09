1;

c = get_conch_adcp;

stn = get_station_from_station_name('mlrf1');
stn = load_all_ndbc_data(stn);
stn = station_wind_to_wave(stn,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_peakwaveper','ndbc_sigwavehgt','ndbc_peakwavedir','ndbc_peakwave_u','ndbc_peakwave_v');
stn = get_ww3_station(stn);
stn = get_erai_station(stn);
stn = adjust_erai_station_waves(stn);

figspath = get_thesis_path('../figs');
scatter_fit_ts_seasons(stn.ndbc_sigwavehgt,c.adcp_sigwavehgt,[],[],'Bulk Hs adj','Conch',[],[],true);  subplots_set('xlim',[0,3],'ylim',[0,3]);
print('-dtiff',fullfile(figspath,[lower(stn.station_name),'-ndbc-scatter-conch-sigwavehgt.tif']));

scatter_fit_ts_seasons(stn.ww3_sigwavehgt,c.adcp_sigwavehgt,[],[],'WW3 Hs','Conch',[],[],true); subplots_set('xlim',[0,3],'ylim',[0,3]);
print('-dtiff',fullfile(figspath,[lower(stn.station_name),'-ww3-scatter-conch-sigwavehgt.tif']));

scatter_fit_ts_seasons(stn.erai_sigwavehgt,c.adcp_sigwavehgt,[],[],'ERAI Hs','Conch',[],[],true);  subplots_set('xlim',[0,3],'ylim',[0,3]);
print('-dtiff',fullfile(figspath,[lower(stn.station_name),'-erai-scatter-conch-sigwavehgt.tif']));

scatter_fit_ts_seasons(stn.erai_sigwavehgt_adj,c.adcp_sigwavehgt,[],[],'ERAI Hs adj','Conch',[],[],true);  subplots_set('xlim',[0,3],'ylim',[0,3]);
print('-dtiff',fullfile(figspath,[lower(stn.station_name),'-erai-adj-scatter-conch-sigwavehgt.tif']));
