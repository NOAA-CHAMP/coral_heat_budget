1;

more off;

datapath = get_thesis_path('../data');

% for yr = [ 2011 2010 2009 2008 2007 2006 2005 2004 2003 2002 2001 2000 1999 1998 1997 1996 1995 1994 1993 1992 1991 1990 1989 1988 1987 ]
for yr = [ 2011 2010 ]

 switch ( yr )
  case 2011,
   mos = 1:2;
  case 2010,
   mos = 11:12;
  otherwise,
   mos = 1:12;
 end;

 for mo = mos(:)'

  dat = []; clear dat;

  fname = fullfile(datapath,'ncep',sprintf('ncep_narr_%04d_%02d.mat', yr, mo));

  if ( exist(fname,'file') )
    disp(['Loading existing ' fname '...']);
    load(fname, 'dat');

    if ( ~isfield(dat,'Temperature_surface') )
      temps_fname = fullfile(datapath,'ncep',sprintf('ncep_narr_temps_%04d_%02d.mat', yr, mo));
      temps = load(temps_fname);
      dat.Temperature_surface = temps.dat.Temperature_surface;
      dat.Temperature_height_above_ground = temps.dat.Temperature_height_above_ground;
      temps = []; clear temps;
    end;

  else
    disp(['Creating ' fname '...']);

    % HACK: Change yr==2010/mo==13 to whatever year/month is currently ONLY PARTIAL at NCEP
    if ( yr==2010 && mo==13 )
      dat = query_ncep_narr_subset([], datenum(yr,mo,1), datenum(yr,mo,13));
    else
      dat = query_ncep_narr_subset([], datenum(yr,mo,1), (datenum(yr,(mo+1),1)-1));
    end;
    disp(['Saving dat to ' fname]);
    save(fname, 'dat');

  end;

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
    %DEBUG:
    disp(stnm);
    [lonix,latix] = gridnbhd_km(dat.lon,dat.lat,stns.lons(ix),stns.lats(ix),0);

    % NOTE: If we add new variables to ALLVARS in QUERY_NCEP_NARR_SUBSET and
    % rerun this script, be sure to change this IF block to simply "if (true)"! 

    if ( ~isfield(dat,stnm) )
      %DEBUG:      disp(['Adding station ' upper(stnm) ' to ' fname]);

      dat.(stnm).lonix = lonix;
      dat.(stnm).latix = latix;

      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'albedo', 'Albedo');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'dewp', 'Dew_point_temperature');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'dlrf', 'Downward_longwave_radiation_flux');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'dsrf', 'Downward_shortwave_radiation_flux');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'evap', 'Evaporation');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'latent_heat_flux', 'Latent_heat_flux');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'apbl_hgt', 'Planetary_boundary_layer_height');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'barom', 'Pressure');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'relhumid', 'Relative_humidity');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'sensible_heat_flux', 'Sensible_heat_flux');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'spechumid', 'Specific_humidity_height_above_ground');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'u_star', 'Surface_friction_velocity');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'sea_t', 'Temperature_surface');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'air_t', 'Temperature_height_above_ground');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cloud_cover', 'Total_cloud_cover');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'precip', 'Total_precipitation');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'ulrf', 'Upward_long_wave_radiation_flux_surface');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'usrf', 'Upward_short_wave_radiation_flux_surface');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'wind_u', 'u_wind_height_above_ground');
      dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'wind_v', 'v_wind_height_above_ground');
    end;

  end;

  %%%% ??? DEBUG
  disp(['Saving again to ' fname]);
  save(fname, 'dat');

 end; %for mo

end; %for yr
