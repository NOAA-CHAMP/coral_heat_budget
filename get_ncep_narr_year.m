function [dat,stns] = get_narr_year(yr)
%function [dat,stns] = get_narr_year(yr)
%
% Download and extract one year's worth of NCEP North American Reanalysis
% (NCEP_NARR) meteorological and surface flux data.
%
% Last Saved Time-stamp: <Thu 2010-08-05 12:59:31  Lew.Gramer>

  datapath = get_thesis_path('../data');

  yrst = num2str(yr);

  fname = fullfile(datapath,'ncep',['ncep_narr_' yrst '.mat']);

  if ( exist(fname,'file') )
    disp(['Loading dat from ' fname]);
    load(fname, 'dat');
  else
    %DEBUG:    dat = query_ncep_narr_subset([], datenum(yr,01,01), datenum(yr,1,2));
    %%%% ??? DEBUG
    dat = query_ncep_narr_subset([], datenum(yr,01,01), datenum(yr,12,31));
    disp(['Saving dat to ' fname]);
    save(fname, 'dat');
  end;

  % Subset our world list of stations to those inside our BBOX
  stns = get_all_station_metadata;
  XV = [dat.bbox(3) dat.bbox(4) dat.bbox(4) dat.bbox(3)];
  YV = [dat.bbox(2) dat.bbox(2) dat.bbox(1) dat.bbox(1)];
  goodix = find( inside(stns.lons, stns.lats, XV, YV) );
  clear XV YV;
  stns.codes  = stns.codes(goodix);
  stns.lons   = stns.lons(goodix);
  stns.lats   = stns.lats(goodix);
  stns.depths = stns.depths(goodix);

  for ix = 1:length(stns.lons)

    stnm = lower(stns.codes{ix});
    disp(stnm);
    [lonix,latix] = gridnbhd_km(dat.lon,dat.lat,stns.lons(ix),stns.lats(ix),0);
    dat.(stnm).lonix = lonix;
    dat.(stnm).latix = latix;

    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'albedo', 'Albedo');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'dewp', 'Dew_point_temperature');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'dlrf', 'Downward_longwave_radiation_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'dsrf', 'Downward_shortwave_radiation_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'evap', 'Evaporation');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'latent_heat_flux', 'Latent_heat_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'barom', 'Pressure');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'relhumid', 'Relative_humidity');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'sensible_heat_flux', 'Sensible_heat_flux');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'spechumid', 'Specific_humidity_height_above_ground');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'u_star', 'Surface_friction_velocity');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cloud_cover', 'Total_cloud_cover');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'precip', 'Total_precipitation');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'ulrf', 'Upward_long_wave_radiation_flux_surface');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'usrf', 'Upward_short_wave_radiation_flux_surface');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'wind_u', 'u_wind_height_above_ground');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'wind_v', 'v_wind_height_above_ground');

  end;

  %%%% ??? DEBUG
  disp(['Saving again to ' fname]);
  save(fname, 'dat');

return;


function stn = extract_ncep_field(dat, stn, stfld, ncepfld)
    stn.(['ncep_' stfld]).date = dat.date(:);
    stn.(['ncep_' stfld]).data = squeeze(dat.(ncepfld)(:,stn.lonix,stn.latix));
return;
