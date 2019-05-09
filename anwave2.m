1;
%%%% SCRIPT to analyze data from Dr. Brian Haus et al. moored wave/currents
%%%% studies off Biscayne Bay, for the Atlantic Ocean Acidification testbed
%%%% Created: Lew.Gramer@noaa.gov, 2011 Sep 22
%%%% Cleaned up: Lew.Gramer@noaa.gov, 2012 Jan 28


doPlots = true;
doPrint = false;

figspath=get_thesis_path('../figs');

if ( ~exist('fwyf1','var') || ~isfield(fwyf1,'station_names') )
  fwyf1 = get_station_from_station_name('fwyf1');
  fwyf1 = get_ngdc_bathy_station(fwyf1);
  fwyf1 = load_all_ndbc_data(fwyf1);
  fwyf1 = station_wind_to_wave(fwyf1,'ndbc_wind1_speed','ndbc_wind1_dir',...
                               'ndbc_peakwaveper','ndbc_sigwavehgt','ndbc_peakwavedir',...
                               'ndbc_peakwave_u','ndbc_peakwave_v');
  fwyf1.ndbc_sigwavehgt_adj.date = fwyf1.ndbc_sigwavehgt.date;
  fwyf1.ndbc_sigwavehgt_adj.data = ((fwyf1.ndbc_sigwavehgt.data.*0.314) + 0.34);
end;


if ( ~exist('stns','var') || ~isfield(stns,'station_names') )
  stns = load_haus_waves;
  stns.bhwv3 = get_ww3_station(stns.bhwv3);
  stns.bhwv3 = get_erai_station(stns.bhwv3);
  stns.bhwv7 = get_ww3_station(stns.bhwv7);
  stns.bhwv7 = get_erai_station(stns.bhwv7);
end;


if ( doPlots )
  fmg;
  [cs,ch]=contour(fwyf1.ngdc_92m_bathy.lon,fwyf1.ngdc_92m_bathy.lat,fwyf1.ngdc_92m_bathy.field,[-2:-2:-10 -20 -30 -80]);
  clabel(cs,ch);
  alllons=[stns.lons(:);fwyf1.lon]; alllats=[stns.lats(:);fwyf1.lat];
  plot(stns.lons,stns.lats,'bs');
  text(stns.lons,stns.lats,strcat(char(stns.station_names),' ^\rightarrow '),'Color','k','HorizontalAlignment','right');
  plot(fwyf1.lon,fwyf1.lat,'kp','MarkerSize',16); plot(fwyf1.lon,fwyf1.lat,'b.');
  axis([min(alllons)-0.01,max(alllons)+0.01,min(alllats)-0.01,max(alllats)+0.01]);
  daspect([cosd(fwyf1.lat),1,1]);

  clear cs ch alllons alllats


  fmg;
  plot_ts(stns.bhwv3.triaxys_sigwavehgt,stns.bhwv4.sontek_sigwavehgt,stns.bhwv7.triaxys_sigwavehgt);
  legend('Stn3 TA: 15m 25.50^oN','Stn4 Sntk: 8.8m 25.50^oN','Stn7 TA: 15m 25.44^oN', 'Location','best');


  disp('Attempt to parameterize wave attenuation between HW3 and HW4');
  scatter_fit_ts(stns.bhwv3.ww3_sigwavehgt,stns.bhwv4.sontek_sigwavehgt,[],[],'WW3 H_s','HW4 H_s',[],[],true);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv3_triaxys_scatter_bhwv4_sontek_sigwavehgt.tiff'));
  end;


  %% Brian Haus Wave experiment site #3 (bhwv3) - on reef slope

  ix3 = find(0<stns.bhwv3.triaxys_peakwavedir.data & ...
             0.05<stns.bhwv3.triaxys_sigwavehgt.data & ...
             stns.bhwv3.triaxys_peakwaveper.data<=15);
  % ix3 = [];

  scatter_fit_ts(stns.bhwv3.ww3_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,[],ix3,'WW3 H_s','HW3 H_s',[],[],true);
  axis([0,4.7,0,4.7]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv3_ww3_scatter_triaxys_sigwavehgt.tiff'));
  end;
  scatter_fit_ts(stns.bhwv3.ww3_peakwaveper,stns.bhwv3.triaxys_peakwaveper,[],ix3,'WW3 wv_p','HW3 peak wv_p',[],[],true);
  axis([1,15,1,15]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv3_ww3_scatter_triaxys_peakwaveper.tiff'));
  end;
  scatter_fit_ts(stns.bhwv3.ww3_peakwavedir,stns.bhwv3.triaxys_peakwavedir,[],ix3,'WW3 wv_d','HW3 peak wv_d',[],[],true);
  axis([0,360,0,360]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv3_ww3_scatter_triaxys_peakwavedir.tiff'));
  end;

  scatter_fit_ts(fwyf1.ndbc_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,[],ix3,'FWY H_s','HW3 H_s',[],[],true);
  axis([0,4.7,0,4.7]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'fwyf1_scatter_bhwv3_triaxys_sigwavehgt.tiff'));
  end;
  scatter_fit_ts(fwyf1.ndbc_peakwaveper,stns.bhwv3.triaxys_peakwaveper,[],ix3,'FWY wv_p','HW3 peak wv_p',[],[],true);
  axis([1,15,1,15]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv3_ndbc_scatter_triaxys_peakwaveper.tiff'));
  end;
  scatter_fit_ts(fwyf1.ndbc_peakwavedir,stns.bhwv3.triaxys_peakwavedir,[],ix3,'FWY wv_d','HW3 peak wv_d',[],[],true);
  axis([0,360,0,360]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv3_ndbc_scatter_triaxys_peakwavedir.tiff'));
  end;

  scatter_fit_ts(stns.bhwv3.erai_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,[],ix3,'ERAI H_s','HW3 H_s',[],[],true);
  axis([0,4.7,0,4.7]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv3_erai_scatter_triaxys_sigwavehgt.tiff'));
  end;
  scatter_fit_ts(stns.bhwv3.erai_peakwaveper,stns.bhwv3.triaxys_peakwaveper,[],ix3,'ERAI wv_p','HW3 peak wv_p',[],[],true);
  axis([1,15,1,15]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv3_erai_scatter_triaxys_peakwaveper.tiff'));
  end;
  scatter_fit_ts(stns.bhwv3.erai_peakwavedir,stns.bhwv3.triaxys_peakwavedir,[],ix3,'ERAI wv_d','HW3 peak wv_d',[],[],true);
  axis([0,360,0,360]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv3_erai_scatter_triaxys_peakwavedir.tiff'));
  end;


  %% Brian Haus Wave experiment site #4 (bhwv4) - inshore of bhwv3

  [ig,ix4] = intersect_dates(stns.bhwv3.triaxys_sigwavehgt.date(ix3),stns.bhwv4.sontek_sigwavehgt.date);
  % ix4=[];

  scatter_fit_ts(stns.bhwv3.ww3_sigwavehgt,stns.bhwv4.sontek_sigwavehgt,[],ix4,'WW3 H_s','HW4 H_s',[],[],true);
  axis([0,4.7,0,4.7]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv4_ww3_scatter_sontek_sigwavehgt.tiff'));
  end;

  scatter_fit_ts(fwyf1.ndbc_sigwavehgt,stns.bhwv4.sontek_sigwavehgt,[],ix4,'FWY H_s','HW4 H_s',[],[],true);
  axis([0,4.7,0,4.7]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'fwyf1_scatter_bhwv4_sontek_sigwavehgt.tiff'));
  end;

  scatter_fit_ts(stns.bhwv3.erai_sigwavehgt,stns.bhwv4.sontek_sigwavehgt,[],ix4,'ERAI H_s','HW4 H_s',[],[],true);
  axis([0,4.7,0,4.7]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv4_erai_scatter_sontek_sigwavehgt.tiff'));
  end;



  %% Brian Haus Wave experiment site #7 (bhwv7) - on reef slope

  ix7 = find(0<stns.bhwv7.triaxys_peakwavedir.data & ...
             0.05<stns.bhwv7.triaxys_sigwavehgt.data & ...
             stns.bhwv7.triaxys_peakwaveper.data<=15);

  scatter_fit_ts(stns.bhwv7.ww3_sigwavehgt,stns.bhwv7.triaxys_sigwavehgt,[],ix7,'WW3 H_s','HW7 H_s',[],[],true);
  axis([0,4.7,0,4.7]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv7_ww3_scatter_triaxys_sigwavehgt.tiff'));
  end;
  scatter_fit_ts(stns.bhwv7.ww3_peakwaveper,stns.bhwv7.triaxys_peakwaveper,[],ix7,'WW3 wv_p','HW7 peak wv_p',[],[],true);
  axis([1,15,1,15]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv7_ww3_scatter_triaxys_peakwaveper.tiff'));
  end;
  scatter_fit_ts(stns.bhwv7.ww3_peakwavedir,stns.bhwv7.triaxys_peakwavedir,[],ix7,'WW3 wv_d','HW7 peak wv_d',[],[],true);
  axis([0,360,0,360]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv7_ww3_scatter_triaxys_peakwavedir.tiff'));
  end;

  scatter_fit_ts(fwyf1.ndbc_sigwavehgt,stns.bhwv7.triaxys_sigwavehgt,[],ix7,'FWY H_s','HW7 H_s',[],[],true);
  axis([0,4.7,0,4.7]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'fwyf1_scatter_bhwv7_triaxys_sigwavehgt.tiff'));
  end;
  scatter_fit_ts(fwyf1.ndbc_peakwaveper,stns.bhwv7.triaxys_peakwaveper,[],ix7,'FWY wv_p','HW7 peak wv_p',[],[],true);
  axis([1,15,1,15]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv7_ndbc_scatter_triaxys_peakwaveper.tiff'));
  end;
  scatter_fit_ts(fwyf1.ndbc_peakwavedir,stns.bhwv7.triaxys_peakwavedir,[],ix7,'FWY wv_d','HW7 peak wv_d',[],[],true);
  axis([0,360,0,360]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv7_ndbc_scatter_triaxys_peakwavedir.tiff'));
  end;

  scatter_fit_ts(stns.bhwv7.erai_sigwavehgt,stns.bhwv7.triaxys_sigwavehgt,[],ix7,'ERAI H_s','HW7 H_s',[],[],true);
  axis([0,4.7,0,4.7]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv7_erai_scatter_triaxys_sigwavehgt.tiff'));
  end;
  scatter_fit_ts(stns.bhwv7.erai_peakwaveper,stns.bhwv7.triaxys_peakwaveper,[],ix7,'ERAI wv_p','HW7 peak wv_p',[],[],true);
  axis([1,15,1,15]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv7_erai_scatter_triaxys_peakwaveper.tiff'));
  end;
  scatter_fit_ts(stns.bhwv7.erai_peakwavedir,stns.bhwv7.triaxys_peakwavedir,[],ix7,'ERAI wv_d','HW7 peak wv_d',[],[],true);
  axis([0,360,0,360]);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,'bhwv7_erai_scatter_triaxys_peakwavedir.tiff'));
  end;


  cutoff=6.5;
  [ix,ig]=intersect_dates(stns.bhwv3.triaxys_sigwavehgt.date,stns.bhwv3.erai_peakwaveper.date(find(stns.bhwv3.erai_peakwaveper.data<=cutoff)));
  scatter_fit_ts(stns.bhwv3.erai_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,[],ix,['ERAI<=',num2str(cutoff),'s'],'HW3 H_S',[],[],true),
  axis([0,2.5,0,2.5]);
  [ix,ig]=intersect_dates(stns.bhwv3.triaxys_sigwavehgt.date,stns.bhwv3.erai_peakwaveper.date(find(stns.bhwv3.erai_peakwaveper.data>cutoff)));
  scatter_curve_fit_ts(stns.bhwv3.erai_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,'power2',[],ix,['ERAI>',num2str(cutoff),'s'],'HW3 H_S'),
  axis([0,2.5,0,2.5]);

  stns.bhwv3 = adjust_erai_station_waves(stns.bhwv3);
  stns.bhwv7 = adjust_erai_station_waves(stns.bhwv7);

  scatter_fit_ts(stns.bhwv3.erai_sigwavehgt_adj,stns.bhwv3.triaxys_sigwavehgt,[],ix3,'ERAI Adj','HW3 H_S',[],[],true);
  scatter_fit_ts(stns.bhwv3.erai_sigwavehgt_adj,stns.bhwv4.sontek_sigwavehgt,[],ix4,'ERAI Adj','HW4 H_S',[],[],true);
  scatter_fit_ts(stns.bhwv7.erai_sigwavehgt_adj,stns.bhwv7.triaxys_sigwavehgt,[],ix7,'ERAI Adj','HW7 H_S',[],[],true);


  cutoff=8;
  [ix,ig]=intersect_dates(stns.bhwv3.triaxys_sigwavehgt.date,stns.bhwv3.ww3_peakwaveper.date(find(stns.bhwv3.ww3_peakwaveper.data<=cutoff)));
  scatter_fit_ts(stns.bhwv3.ww3_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,[],ix,['WW3<=',num2str(cutoff),'s'],'HW3 H_S',[],[],true),
  axis([0,2.5,0,2.5]);
  [ix,ig]=intersect_dates(stns.bhwv3.triaxys_sigwavehgt.date,stns.bhwv3.ww3_peakwaveper.date(find(stns.bhwv3.ww3_peakwaveper.data>cutoff)));
  scatter_curve_fit_ts(stns.bhwv3.ww3_sigwavehgt,stns.bhwv3.triaxys_sigwavehgt,'power2',[],ix,['WW3>',num2str(cutoff),'s'],'HW3 H_S'),
  axis([0,2.5,0,2.5]);

  stns.bhwv3 = adjust_ww3_station_waves(stns.bhwv3);
  stns.bhwv7 = adjust_ww3_station_waves(stns.bhwv7);

  scatter_fit_ts(stns.bhwv3.ww3_sigwavehgt_adj,stns.bhwv3.triaxys_sigwavehgt,[],ix3,'WW3 Adj','HW3 H_S',[],[],true);
  scatter_fit_ts(stns.bhwv3.ww3_sigwavehgt_adj,stns.bhwv4.sontek_sigwavehgt,[],ix4,'WW3 Adj','HW4 H_S',[],[],true);
  scatter_fit_ts(stns.bhwv7.ww3_sigwavehgt_adj,stns.bhwv7.triaxys_sigwavehgt,[],ix7,'WW3 Adj','HW7 H_S',[],[],true);

end;

clear doPlots doPrint figspath ans
