1;

stnm='sanf1';

stn = get_station_from_station_name(stnm);
stn = load_all_ndbc_data(stn);
stn = get_ww3_station(stn);
stn = station_wind_to_wave(stn,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_peakwaveper','ndbc_sigwavehgt','ndbc_peakwavedir','ndbc_peakwave_u','ndbc_peakwave_v');

scatter_fit_ts_seasons(stn.ww3_sigwavehgt,stn.ndbc_sigwavehgt,[],[],'WW3 H_s','NDBC',[],[],true);
scatter_fit_ts_seasons(stn.ww3_peakwaveper,stn.ndbc_peakwaveper,[],[],'WW3 wv_p','NDBC',[],[],true);
scatter_fit_ts_seasons(stn.ww3_peakwave_u,stn.ndbc_peakwave_u,[],[],'WW3 wv_u','NDBC',[],[],true);
scatter_fit_ts_seasons(stn.ww3_peakwave_v,stn.ndbc_peakwave_v,[],[],'WW3 wv_v','NDBC',[],[],true);
