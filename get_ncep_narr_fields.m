function flds = get_ncep_narr_fields(begdt,enddt,minlon,maxlon,minlat,maxlat,varsubset)
%function flds = get_ncep_narr_fields(begdt,enddt,minlon,maxlon,minlat,maxlat,varsubset)

  persistent LONS LATS

  datapath = get_thesis_path('../data');

  ALLVARS = { ...
      'Albedo[0]', ...
      'Dew_point_temperature[0][0]', ...
      'Downward_longwave_radiation_flux[0]', ...
      'Downward_shortwave_radiation_flux[0]', ...
      'Evaporation[0]', ...
      'Latent_heat_flux[0]', ...
      'Pressure[0][0]', ...
      'Relative_humidity[0][0]', ...
      'Sensible_heat_flux[0]', ...
      'Specific_humidity_height_above_ground[0][0]', ...
      'Surface_friction_velocity[0]', ...
      'Total_cloud_cover[0]', ...
      'Total_precipitation[0]', ...
      'Upward_long_wave_radiation_flux_surface[0]', ...
      'Upward_short_wave_radiation_flux_surface[0]', ...
      'u_wind_height_above_ground[0][0]', ...
      'v_wind_height_above_ground[0][0]', ...
            };


  flds.vars = ALLVARS;

  if ( exist('varsubset') && isvector(varsubset) )
    if ( isnumeric(varsubset) )
      flds.vars = flds.vars(varsubset);
    elseif ( iscellstr(varsubset) )
      matchix = [];
      for ix = 1:length(varsubset)
        matchix = [matchix strmatch(lower(varsubset{ix}),lower(flds.vars))];
      end;
      flds.vars = flds.vars(matchix);
    end;
  end;


  flds.date = floor(begdt):(3/24):floor(enddt+1);
  [yrs,mos,dys,hrs,mis,scs] = datevec(flds.date);
  flds.jday = floor(flds.date) - datenum(yrs,1,1) + 1;

  if ( numel(LONS) == 0 || numel(LATS) == 0 )
    LONS = load(fullfile(datapath,'NCEP_NARR_lon.dat'));
    LATS = load(fullfile(datapath,'NCEP_NARR_lat.dat'));
  end;

  [minlonix,minlatix] = gridnbhd_km(LONS,LATS,minlon,minlat,0);
  [maxlonix,maxlatix] = gridnbhd_km(LONS,LATS,maxlon,maxlat,0);

  flds.LONS = LONS(minlatix:maxlatix, minlonix:maxlonix);
  flds.LATS = LATS(minlatix:maxlatix, minlonix:maxlonix);

  % http://nomads.ncdc.noaa.gov/thredds/dodsC/narr/201001/20100119/narr-a_221_20100119_0000_000.grb.dods?Downward_longwave_radiation_flux[0:1:0][lonix][latix],Downward_shortwave_radiation_flux[0:1:0][lonix][latix],Evaporation[0:1:0][lonix][latix],Latent_heat_flux[0:1:0][lonix][latix],Precipitation_rate[0:1:0][lonix][latix],Pressure_nearest_grid_point[0:1:0][lonix][latix],Relative_humidity[0:1:0][0:1:0][lonix][latix],Sensible_heat_flux[0:1:0][lonix][latix],Specific_humidity_height_above_ground[0:1:0][0:1:0][lonix][latix],Total_cloud_cover[0:1:0][lonix][latix],Total_precipitation[0:1:0][lonix][latix],u_wind_height_above_ground[0:1:0][0:1:0][lonix][latix],Upward_long_wave_radiation_flux_surface[0:1:0][lonix][latix],Upward_short_wave_radiation_flux_surface[0:1:0][lonix][latix],v_wind_height_above_ground[0:1:0][0:1:0][lonix][latix]


  for vix = 1:length(flds.vars)
    flds.varnames{vix} = lower(regexprep(flds.vars{vix}, '\[.*\]', ''));
    flds.(lower(flds.varnames{vix})) = repmat(nan, [numel(flds.date) size(LONS)]);
  end;

  for dix = 1:length(flds.date)

    base_dataquery = ...
        sprintf('http://nomads.ncdc.noaa.gov/thredds/dodsC/narr/%04d%02d/%04d%02d%02d/narr-a_221_%04d%02d%02d_%02d00_000.grb.dods', ...
                yrs(dix), mos(dix), yrs(dix), mos(dix), dys(dix),...
                yrs(dix), mos(dix), dys(dix), hrs(dix));

    for vix = 1:length(flds.vars)
      dataquery = sprintf('%s?%s[%g:%g][%g:%g]', base_dataquery, flds.vars{vix}, ...
                          minlonix, maxlonix, minlatix, maxlatix);
      dataquery,
      dat = loaddap(dataquery);
disp('Got it!');
keyboard;
      flds.(lower(flds.varnames{vix}))(dix,:,:) = dat;
    end;

  end;


return;
