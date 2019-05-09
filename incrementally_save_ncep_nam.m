1;

more off;

datapath = get_thesis_path('../data');

% for yr = [ 2010 2009 2008 2007 2006 2005 2004 2003 2002 2001 2000 ]
for yr = [ 2009 ]

 switch ( yr )
  case 2009,
   mos = 2;
  case 2010,
   mos = 1:2;
  otherwise,
   mos = 1:12;
 end;

 for mo = mos(:)'

  dat = []; clear dat;

  fname = fullfile(datapath, sprintf('ncep_nam_%04d_%02d.mat', yr, mo));

  disp(['Creating ' fname '...']);

  %DEBUG:    dat = query_ncep_narr_subset([], datenum(yr,mo,1), datenum(yr,mo,2));
  %%%% ??? DEBUG
  dat = query_ncep_nam_subset([], datenum(yr,mo,1), (datenum(yr,(mo+1),1)-1));
  disp(['Saving dat to ' fname]);
  save(fname, 'dat');

  % Subset our world list of stations to those inside our BBOX
  if ( ~exist('stns','var') || isempty(stns) )
    stns = get_all_station_metadata;
    XV = [dat.bbox(3) dat.bbox(4) dat.bbox(4) dat.bbox(3)];
    YV = [dat.bbox(2) dat.bbox(2) dat.bbox(1) dat.bbox(1)];
    goodix = find( inside(stns.lons, stns.lats, XV, YV) );
    clear XV YV;
    stns.codes  = stns.codes(goodix);
    stns.lons   = stns.lons(goodix);
    stns.lats   = stns.lats(goodix);
    stns.depths = stns.depths(goodix);
  end;


  for ix = 1:length(stns.lons)

    stnm = lower(stns.codes{ix});
    disp(stnm);
    [lonix,latix] = gridnbhd_km(dat.lon,dat.lat,stns.lons(ix),stns.lats(ix),0);
    dat.(stnm).lonix = lonix;
    dat.(stnm).latix = latix;

    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_albedo', 'Albedo');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_dewp', 'Dew_point_temperature');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_dlrf', 'Downward_long_wave_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_dsrf', 'Downward_short_wave_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_latent_heat_flux', 'Latent_heat_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_barom', 'Pressure');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_relhumid', 'Relative_humidity');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_sensible_heat_flux', 'Sensible_heat_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_u_star', 'Friction_velocity');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_cloud_cover', 'Total_cloud_cover');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_precip', 'Total_precipitation');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_ulrf', 'Upward_long_wave_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_usrf', 'Upward_short_wave_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_wind_u', 'u_wind_height_above_ground');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'nam_wind_v', 'v_wind_height_above_ground');

  end;

  %%%% ??? DEBUG
  disp(['Saving again to ' fname]);
  save(fname, 'dat');

 end; %for mo

end; %for yr
