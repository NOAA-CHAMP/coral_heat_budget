function flds = get_ncep_narr_daily(begdt,enddt,minlon,maxlon,minlat,maxlat,varsubset)
%function flds = get_ncep_narr_fields(begdt,enddt,minlon,maxlon,minlat,maxlat,varsubset)

  persistent LONS LATS

  datapath = get_thesis_path('../data');


  % NCEP_NARR_DAILIES variables, equivalent NCEP NARR (32km) variables
  ALLVARS = { ...
      'apcp',		'Total_precipitation'				, ...
      'dlwrfsfc',	'Downward_longwave_radiation_flux'		, ...
      'dswrfsfc',	'Downward_shortwave_radiation_flux'		, ...
      'dpt2m',		'Dew_point_temperature[0]'			, ...
      'evpsfc',		'Evaporation'					, ...
      'fricvsfc',	'u_star?'					, ...
      'lhtfl',		'Latent_heat_flux'				, ...
      'rh2m',		'Relative_humidity[0]'				, ...
      'shtfl',		'Sensible_heat_flux'				, ...
      'spfh2m',		'Specific_humidity_height_above_ground[0]'	, ...
      'tcdc',		'Total_cloud_cover'				, ...
      'ulwrfsfc',	'Upward_long_wave_radiation_flux_surface'	, ...
      'uswrfsfc',	'Upward_short_wave_radiation_flux_surface'	, ...
      'ugrd10m',	'u_wind_height_above_ground[0]'			, ...
      'vgrd10m',	'v_wind_height_above_ground[0]'			, ...
      'uflxsfc',	'u_wind_stress_surface?'			, ...
      'vflxsfc',	'v_wind_stress_surface?'			, ...
            };


  flds.vars = ALLVARS(1:2:end);

  if ( exist('varsubset') && isvector(varsubset) )
    if ( isnumeric(varsubset) )
      flds.vars = flds.vars(varsubset);
    elseif ( iscellstr(varsubset) )
      matchix = [];
      for ix = 1:length(varsubset)
        matchix = [matchix strmatch(lower(varsubset{ix}),lower(flds.vars))];
      end;
      flds.vars = flds.vars(matchix);
    elseif ( ischar(varsubset) )
      matchix = strmatch(lower(varsubset),lower(flds.vars));
      flds.vars = flds.vars(matchix);
    else
      error('Optional arg VARSUBSET must be numeric vector, field name or cellstr!');
    end;
  end;


  if ( begdt == enddt )
    flds.date = floor(begdt):(3/24):floor(enddt)+1-eps;
    datestr(flds.date),
  else
    flds.date = floor(begdt):(3/24):ceil(enddt+eps);
  end;
  [yrs,mos,dys,hrs,mis,scs] = datevec(flds.date);
  flds.jday = floor(flds.date) - datenum(yrs,1,1) + 1;

  if ( numel(LONS) == 0 || numel(LATS) == 0 )
    lons = -220.0000:0.3750:-000.6250;
    lats =  +00.0000:0.3750: +89.6250;
    [LONS,LATS] = meshgrid(lons,lats);
  end;

  [minlonix,minlatix] = gridnbhd_km(LONS,LATS,minlon,minlat,0);
  [maxlonix,maxlatix] = gridnbhd_km(LONS,LATS,maxlon,maxlat,0);

  flds.LONS = LONS(minlatix:maxlatix, minlonix:maxlonix);
  flds.LATS = LATS(minlatix:maxlatix, minlonix:maxlonix);

  % http://nomads.ncdc.noaa.gov/dods/NCEP_NARR_DAILY/197912/19791210/narr-a_221_19791210_0000_000?rh2m[0][0:239][0:585]

  for vix = 1:length(flds.vars)
    flds.varnames{vix} = flds.vars{vix};
    flds.(lower(flds.varnames{vix})) = repmat(nan, [numel(flds.date) size(flds.LONS)]);
  end;

  for dix = 1:length(flds.date)

    base_dataquery = ...
        sprintf('http://nomads.ncdc.noaa.gov/dods/NCEP_NARR_DAILY/%04d%02d/%04d%02d%02d/narr-a_221_%04d%02d%02d_0000_000',...
                yrs(dix), mos(dix), yrs(dix), mos(dix), dys(dix),...
                yrs(dix), mos(dix), dys(dix));

    for vix = 1:length(flds.vars)
      dataquery = sprintf('%s?%s[%d][%g:%g][%g:%g]', base_dataquery, flds.vars{vix}, ...
                          round(hrs(dix)/3), minlatix, maxlatix, minlonix, maxlonix);
      %DEBUG:
      dataquery,
      datst = loaddap(dataquery);
      dat = datst.(flds.vars{vix}).(flds.vars{vix});
      dat(dat > 9e20) = nan;
      datst = []; clear datst;
      flds.(lower(flds.varnames{vix}))(dix,:,:) = dat;
    end;

  end;


return;
