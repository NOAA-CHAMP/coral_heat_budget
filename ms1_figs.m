function stn = ms1_figs(stnm)

  stn.station_name = upper(stnm);
  [stn.lon,stn.lat,stn.depth]=get_station_coords(stnm);
  stn = load_all_ndbc_data(stn);
  
  stn = verify_variable(stn,'ndbc_wind1_u');
  stn = verify_variable(stn,'ndbc_wind1_v');

  stn.ndbc_wind1_u.data = kts2mps(stn.ndbc_wind1_u.data);
  stn.ndbc_wind1_v.data = kts2mps(stn.ndbc_wind1_v.data);

  stn = station_optimal_isobath_orientation(stn);

  stn = station_reorient_vectors(stn,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');

  %station_monthly_boxplot(stn,fld,ttl,xlbl,ylbl,ylm,doyrs,domos,tol,doPrint)

  fh=figure;
  station_monthly_boxplot(stn,'ndbc_wind1_xshore','none','none','none',[-30,30],1995:2010,[],0.5,true);
  grid on;
  set(fh,'Name',[stn.station_name ' Cross-Shore Wind 1995-2010']);

  fh=figure;
  station_monthly_boxplot(stn,'ndbc_wind1_lshore','none','none','none',[-30,30],1995:2010,[],0.5,true);
  grid on;
  set(fh,'Name',[stn.station_name ' Along-Shore Wind 1995-2010']);

return;
