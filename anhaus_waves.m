1;
%%%% SCRIPT to analyze data from Dr. Brian Haus et al. moored wave/currents
%%%% studies off Biscayne Bay, for the Atlantic Ocean Acidification testbed
%%%% Created: Lew.Gramer@noaa.gov, 2011 Sep 22


doPlots = true;

if ( ~exist('stn','var') || isempty(stn) )
    % stn = get_station_from_station_name('aoat_broad_key_2');
    stn = get_station_from_station_name('fwyf1');
    stn = get_ngdc_bathy_station(stn);
    stn = load_all_ndbc_data(stn);
    stn = get_ww3_station(stn);
    stn = station_wind_to_wave(stn,'ndbc_wind1_speed','ndbc_wind1_dir',...
                               'ndbc_peakwaveper','ndbc_sigwavehgt','ndbc_peakwavedir',...
                               'ndbc_peakwave_u','ndbc_peakwave_v');

    stn.ndbc_sigwavehgt_adj.date = stn.ndbc_sigwavehgt.date;
    %stn.ndbc_sigwavehgt_adj.data = ((stn.ndbc_sigwavehgt.data.*0.314) + 0.34);
    stn.ndbc_sigwavehgt_adj.data = ((stn.ndbc_sigwavehgt.data.*0.15) + 0.25);
end;

if ( ~exist('stns','var') || ~isfield(stns,'station_names') )
  stns = load_haus_waves;
  stns.bhwv3 = get_ww3_station(stns.bhwv3);
  stns.bhwv7 = get_ww3_station(stns.bhwv7);
end;

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( ~exist('stns','var') || ~isfield(stns,'station_names') )
  stns.station_names = { ...
      'bhwv1', ...
      'bhwv2', ...
      'bhwv3', ...
      'bhwv4', ...
      'bhwv5', ...
      'bhwv6', ...
      'bhwv7', ...
      'bhwv8', ...
                   };
  stns.station_long_names = { ...
      'Cluster-1 (bottom)', ...
      'Cluster-2 (bottom)', ...
      'Cluster-3 (surface)', ...
      'Cluster-4 (bottom)', ...
      'Offshore-5 (botm-float)', ...
      'Inshore-6 (bottom)', ...
      'South-7 (surface)', ...
      'South-8 (bottom)', ...
                   };
  stns.coords_dms = [ ...
      -80 06 36.3 ;  25 30 02.1 ; ...
      -80 06 22.2 ;  25 29 53.7 ; ...
      -80 06 08.9 ;  25 29 56.7 ; ...
      -80 06 31.2 ;  25 29 55.2 ; ...
      -80 05 53.2 ;  25 23 48.4 ; ...
      -80 08 30.2 ;  25 28 46.0 ; ...
      -80 07 01.2 ;  25 26 08.8 ; ...
      -80 06 48.4 ;  25 28 13.7 ; ...
                    ];
  stns.coords = dms2degrees(stns.coords_dms);
  stns.lons = stns.coords(1:2:end);
  stns.lats = stns.coords(2:2:end);
  stns.depths = [ ...
      10.0 ...
      14.0 ...
      15.0 ...
      8.8 ...
      100 ...
      6.0 ...
      15.0 ...
      9.0 ...
                ];
end;
if ( ~isfield(stns,'bhwv4') )
  load('Hs_SNTK.mat');
  stns.bhwv4.station_name = stns.station_names{4};
  stns.bhwv4.lon = stns.lons(4);
  stns.bhwv4.lat = stns.lats(4);
  stns.bhwv4.depth = stns.depths(4);
  stns.bhwv4.sontek_sigwavehgt.date = datenum(2005,3,19,17,0,0) + ([0:numel(Hs_SNTK)-1]'./24);
  stns.bhwv4.sontek_sigwavehgt.data = Hs_SNTK;
  Hs_SNTK=[]; clear Hs_SNTK;
end;
if ( ~isfield(stns,'bhwv3') )
  stns.bhwv3.station_name = stns.station_names{3};
  stns.bhwv3.lon = stns.lons(3);
  stns.bhwv3.lat = stns.lats(3);
  stns.bhwv3.depth = stns.depths(3);
  stns.bhwv3 = get_triaxys_waves('Summarytas01290.txt');
end;
if ( ~isfield(stns,'bhwv7') )
  stns.bhwv7.station_name = stns.station_names{7};
  stns.bhwv7.lon = stns.lons(7);
  stns.bhwv7.lat = stns.lats(7);
  stns.bhwv7.depth = stns.depths(7);
  stns.bhwv7 = get_triaxys_waves('Summarytab00651.txt');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


if ( doPlots )
    fmg;
    [cs,ch]=contour(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,[-2:-2:-10 -20 -30 -80]);
    clabel(cs,ch);
    alllons=[stns.lons(:);stn.lon]; alllats=[stns.lats(:);stn.lat];
    plot(stns.lons,stns.lats,'bs');
    text(stns.lons,stns.lats,strcat(char(stns.station_names),' ^\rightarrow '),'Color','k','HorizontalAlignment','right');
    plot(stn.lon,stn.lat,'kp','MarkerSize',16); plot(stn.lon,stn.lat,'b.');
    axis([min(alllons)-0.01,max(alllons)+0.01,min(alllats)-0.01,max(alllats)+0.01]);
    daspect([cosd(stn.lat),1,1]);

    clear cs ch alllons alllats
end;

if ( doPlots )
  fmg;
  plot_ts(stns.bhwv3.triaxys_sigwavehgt,stns.bhwv4.sontek_sigwavehgt,stns.bhwv7.triaxys_sigwavehgt);
  legend('Stn3 TA: 15m 25.50^oN','Stn4 Sntk: 8.8m 25.50^oN','Stn7 TA: 15m 25.44^oN', 'Location','best');


  % goodix = 1:length(stns.bhwv3.triaxys_sigwaveper.data);
  goodix = find(stns.bhwv3.triaxys_peakwaveper.data<=15);

  scatter_fit_ts(stns.bhwv3.ww3_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,[],goodix,'HW3 WW3 H_s','HW3 H_s',[],[],true);

  scatter_fit_ts(stn.ww3_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,[],goodix,'WW3 H_s','HW3 H_s',[],[],true);
  scatter_fit_ts(stn.ndbc_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,[],goodix,'NDBC H_s','HW3 H_s',[],[],true);
  scatter_fit_ts(stn.ww3_sigwavehgt,stns.bhwv3.triaxys_avgwavehgt,[],goodix,'WW3 H_s','HW3 avg wv_h',[],[],true);
  scatter_fit_ts(stn.ndbc_sigwavehgt,stns.bhwv3.triaxys_avgwavehgt,[],goodix,'NDBC H_s','HW3 avg wv_h',[],[],true);
  scatter_fit_ts(stn.ww3_sigwavehgt,stns.bhwv3.triaxys_maxwavehgt,[],goodix,'WW3 H_s','HW3 max wv_h',[],[],true);
  scatter_fit_ts(stn.ndbc_sigwavehgt,stns.bhwv3.triaxys_maxwavehgt,[],goodix,'NDBC H_s','HW3 max wv_h',[],[],true);

  scatter_fit_ts(stn.ndbc_sigwavehgt_adj,stns.bhwv3.triaxys_sigwavehgt,[],goodix,'Adj. NDBC H_s','HW3 H_s',[],[],true);
  scatter_fit_ts(stn.ndbc_sigwavehgt_adj,stns.bhwv3.triaxys_avgwavehgt,[],goodix,'Adj. NDBC H_s','HW3 avg wv_h',[],[],true);
  scatter_fit_ts(stn.ndbc_sigwavehgt_adj,stns.bhwv3.triaxys_maxwavehgt,[],goodix,'Adj. NDBC H_s','HW3 max wv_h',[],[],true);


  scatter_fit_ts(stn.ww3_peakwaveper,stns.bhwv3.triaxys_sigwaveper,[],goodix,'WW3 wv_p','HW3 wv_p_S',[],[],true);
  scatter_fit_ts(stn.ndbc_peakwaveper,stns.bhwv3.triaxys_sigwaveper,[],goodix,'NDBC wv_p','HW3 wv_p_S',[],[],true);
  scatter_fit_ts(stn.ww3_peakwaveper,stns.bhwv3.triaxys_avgwaveper,[],goodix,'WW3 wv_p','HW3 avg wv_p',[],[],true);
  scatter_fit_ts(stn.ndbc_peakwaveper,stns.bhwv3.triaxys_avgwaveper,[],goodix,'NDBC wv_p','HW3 avg wv_p',[],[],true);
  scatter_fit_ts(stn.ww3_peakwaveper,stns.bhwv3.triaxys_peakwaveper,[],goodix,'WW3 wv_p','HW3 peak wv_p',[],[],true);
  scatter_fit_ts(stn.ndbc_peakwaveper,stns.bhwv3.triaxys_peakwaveper,[],goodix,'NDBC wv_p','HW3 peak wv_p',[],[],true);

  scatter_fit_ts(stn.ww3_peakwavedir,stns.bhwv3.triaxys_peakwavedir,[],goodix,'WW3 wv_d','HW3 peak wv_d',[],[],true);
  scatter_fit_ts(stn.ndbc_peakwavedir,stns.bhwv3.triaxys_peakwavedir,[],goodix,'NDBC wv_d','HW3 peak wv_d',[],[],true);


  scatter_fit_ts(stn.ndbc_sigwavehgt,stns.bhwv4.sontek_sigwavehgt,[],[],'NDBC H_s','HW4 H_s',[],[],true);
  scatter_fit_ts(stn.ndbc_sigwavehgt_adj,stns.bhwv4.sontek_sigwavehgt,[],[],'Adj. NDBC H_s','HW4 H_s',[],[],true);


  goodix = find(stns.bhwv7.triaxys_sigwaveper.data<=15);

  scatter_fit_ts(stns.bhwv7.ww3_sigwavehgt,stns.bhwv7.triaxys_sigwavehgt,[],goodix,'HW7 WW3 H_s','HW7 H_s',[],[],true);

  scatter_fit_ts(stn.ww3_sigwavehgt,stns.bhwv7.triaxys_sigwavehgt,[],goodix,'WW3 H_s','HW7 H_s',[],[],true);
  scatter_fit_ts(stn.ndbc_sigwavehgt,stns.bhwv7.triaxys_sigwavehgt,[],goodix,'NDBC H_s','HW7 H_s',[],[],true);
  scatter_fit_ts(stn.ww3_sigwavehgt,stns.bhwv7.triaxys_avgwavehgt,[],goodix,'WW3 H_s','HW7 avg wv_h',[],[],true);
  scatter_fit_ts(stn.ndbc_sigwavehgt,stns.bhwv7.triaxys_avgwavehgt,[],goodix,'NDBC H_s','HW7 avg wv_h',[],[],true);
  scatter_fit_ts(stn.ww3_sigwavehgt,stns.bhwv7.triaxys_maxwavehgt,[],goodix,'WW3 H_s','HW7 max wv_h',[],[],true);
  scatter_fit_ts(stn.ndbc_sigwavehgt,stns.bhwv7.triaxys_maxwavehgt,[],goodix,'NDBC H_s','HW7 max wv_h',[],[],true);

  scatter_fit_ts(stn.ww3_peakwaveper,stns.bhwv7.triaxys_sigwaveper,[],goodix,'WW3 wv_p','HW7 wv_p_S',[],[],true);
  scatter_fit_ts(stn.ndbc_peakwaveper,stns.bhwv7.triaxys_sigwaveper,[],goodix,'NDBC wv_p','HW7 wv_p_S',[],[],true);
  scatter_fit_ts(stn.ww3_peakwaveper,stns.bhwv7.triaxys_avgwaveper,[],goodix,'WW3 wv_p','HW7 avg wv_p',[],[],true);
  scatter_fit_ts(stn.ndbc_peakwaveper,stns.bhwv7.triaxys_avgwaveper,[],goodix,'NDBC wv_p','HW7 avg wv_p',[],[],true);
  scatter_fit_ts(stn.ww3_peakwaveper,stns.bhwv7.triaxys_peakwaveper,[],goodix,'WW3 wv_p','HW7 peak wv_p',[],[],true);
  scatter_fit_ts(stn.ndbc_peakwaveper,stns.bhwv7.triaxys_peakwaveper,[],goodix,'NDBC wv_p','HW7 peak wv_p',[],[],true);

  scatter_fit_ts(stn.ww3_peakwavedir,stns.bhwv7.triaxys_peakwavedir,[],goodix,'WW3 wv_d','HW7 peak wv_d',[],[],true);
  scatter_fit_ts(stn.ndbc_peakwavedir,stns.bhwv7.triaxys_peakwavedir,[],goodix,'NDBC wv_d','HW7 peak wv_d',[],[],true);

  scatter_fit_ts(stn.ndbc_sigwavehgt_adj,stns.bhwv7.triaxys_sigwavehgt,[],goodix,'Adj. NDBC H_s','HW7 H_s',[],[],true);
end;

clear doPlots ans
