function dat = query_ncep_cfsr_subset(yr,mo,bbox,vars)
%function dat = query_ncep_cfsr_subset(yr,mo,bbox,vars)
%
% Extract a spatio-temporal subset of the NCEP Climate Forecast System
% Reanalysis (CFSR) directly from the NOMADS THREDDS site at ncdc.noaa.gov.
% YR and MO determine GRIB2 file to netCDF subset for each variable, e.g.,
% 1987, 01 accesses all variables for January 1987. Optional arg BBOX
% spatially subsets data (DEFAULT: BBOX=[24 27 -83 -79], Florida Reef Tract
% region). Optional arg VARS (cellstr) specifies the variables to extract.
% DEFAULT array of - and the complete list of available - CFSR variables
% (with GRIB2 file base-name, and name of the variable added to returned
% struct DAT, resp., both shown in parentheses afterward) are as follows:
%
%       Momentum_flux_u_component		(wndstrs, wndstrs_u)
%       Momentum_flux_v_component		(wndstrs, wndstrs_v)
%       U-component_of_wind			(wnd10m, wnd10m_u)
%       V-component_of_wind			(wnd10m, wnd10m_v)
%       Upward_Short-Wave_Rad_Flux		(uswsfc)
%       Upward_Long-Wave_Rad_Flux		(ulwsfc)
%       Temperature				(tmpsfc)
%       Temperature				(tmp2m)
%       Sensible_heat_net_flux			(shtfl)
%       Specific_humidity			(q2m)
%       Pressure				(pressfc)
%       Precipitation_rate			(prate)
%       V-component_of_current			(ocnv5)
%       V-component_of_current			(ocnu5)
%       Vertical_velocity_geometric		(ocnvv55)
%       Potential_temperature			(ocnsst)
%       Salinity				(ocnsal5)
%       Geometric_Depth_Below_Sea_Surface	(ocnmld)
%       Latent_heat_net_flux			(lhtfl)
%       Downward_Short-Wave_Rad_Flux		(dswsfc)
%       Downward_Long-Wave_Rad_Flux		(dlwsfc)
%
% EXAMPLE: dat = query_ncep_cfsr_subset(2009,01,[],{'windstr_u','windstr_v'});
%  Returns struct DAT with fields DAT.WINDSTR_U and DAT.WINDSTR_V containing
%  time series of wind stress U- and V components from January 2009 CFSR,
%  subsetted to the region of the Florida reef tract (see above for coords).
%
% CALLS: netCDF-Java Toolbox (matlab-njTbx), via its interface that is like
% old NetCDF Toolbox by Chuck Denham: http://sourceforge.net/apps/trac/njtbx
%
% Last Saved Time-stamp: <Thu 2010-08-05 12:58:17  Lew.Gramer>


error('WAIT! Add air and sea temperature (back) to this code before rerunning!');


  %DEBUG:
  tic,

  if ( ~exist('bbox','var') || isempty(bbox) )
    bbox = [ 24 27 -83 -79 ];
  end;


  % Maximum number of times to retry each URL (i.e., each NOMADS time)
  maxRetries = 3;

%   ALLVARS = { ...
%       'wndstrs_u' ,	'Momentum_flux_u_component' ; ...
%       'wndstrs_v' ,	'Momentum_flux_v_component' ; ...
%       'wnd10m_u' ,	'U-component_of_wind' ; ...
%       'wnd10m_v' ,	'V-component_of_wind' ; ...
%       'uswsfc' ,	'Upward_Short-Wave_Rad_Flux' ; ...
%       'ulwsfc' ,	'Upward_Long-Wave_Rad_Flux' ; ...
%       'tmpsfc' ,	'Temperature' ; ...
%       'tmp2m' ,		'Temperature' ; ...
%       'shtfl' ,		'Sensible_heat_net_flux' ; ...
%       'q2m' ,		'Specific_humidity' ; ...
%       'pressfc' ,	'Pressure' ; ...
%       'prate' ,		'Precipitation_rate' ; ...
%       'ocnv5' ,		'V-component_of_current' ; ...
%       'ocnu5' ,		'V-component_of_current' ; ...
%       'ocnvv55' ,	'Vertical_velocity_geometric' ; ...
%       'ocnsst' ,	'Potential_temperature' ; ...
%       'ocnsal5' ,	'Salinity' ; ...
%       'ocnmld' ,	'Geometric_Depth_Below_Sea_Surface' ; ...
%       'lhtfl' ,		'Latent_heat_net_flux' ; ...
%       'dswsfc' ,	'Downward_Short-Wave_Rad_Flux' ; ...
%       'dlwsfc' ,	'Downward_Long-Wave_Rad_Flux' ; ...
%             };

  ALLVARS = { ...
      'uswsfc' ,	'Upward_Short-Wave_Rad_Flux' ; ...
      'ulwsfc' ,	'Upward_Long-Wave_Rad_Flux' ; ...
      'shtfl' ,		'Sensible_heat_net_flux' ; ...
      'q2m' ,		'Specific_humidity' ; ...
      'prate' ,		'Precipitation_rate' ; ...
      'lhtfl' ,		'Latent_heat_net_flux' ; ...
      'dswsfc' ,	'Downward_Short-Wave_Rad_Flux' ; ...
      'dlwsfc' ,	'Downward_Long-Wave_Rad_Flux' ; ...
            };

  ALLVARS = { ...
      'lhtfl' ,		'Latent_heat_net_flux' ; ...
            };

  if ( ~exist('vars','var') || isempty(vars) )
    vars = ALLVARS(:,2);
    flds = ALLVARS(:,1);

  else
    if ( ischar(vars) )
      vars = { vars };
    end;
    for varix = 1:length(vars)
      ix = find(strcmpi(ALLVARS(:,1),vars{varix}), 1);
      vars{varix} = ALLVARS{ix,2};
      flds{varix} = ALLVARS{ix,1};
    end;
  end;

  fnms = regexprep(flds,'_.*$','');

  dat.bbox = bbox;
  dat.flds = flds;

  dat.date = datenum(yr,mo,1):(1/24):(datenum(yr,(mo+1),1) - (1/(24*60)));

  %DEBUG:    disp(datestr(dat.date([1 end])));

  for vix = 1:length(vars)

    var = vars{vix};
    fnm = fnms{vix};
    fld = flds{vix};

    %http://nomads.ncdc.noaa.gov/thredds/ncss/grid/cfsr1hr/198701/dswsfc.gdas.198701.grb2?var=Downward_Short-Wave_Rad_Flux&temporal=all&spatial=bb&north=27&south=24&west=-83&east=-79&addLatLon=true
    %horizStride=1&

    url = sprintf( 'http://nomads.ncdc.noaa.gov/thredds/ncss/grid/cfsr1hr/%04d%02d/%s.gdas.%04d%02d.grb2?var=%s&temporal=all&spatial=bb&south=%g&north=%g&west=%g&east=%g&addLatLon=true', ...
                   yr, mo, fnm, yr, mo, var, ...
                   bbox(1), bbox(2), bbox(3), bbox(4) ...
                   );
    %DEBUG:    url,
url = '../data/lhtfl.gdas.198701.grb2.nc';

    nc = [];
    tried = 0;
    while ( isempty(nc) && (tried < maxRetries) )
      % This is in here to keep from overheating NCDC's NOMAD servers!
      pause(5);
      nc = mDataset(url);
      tried = tried + 1;
    end;

    if ( isempty(nc) )
      warning('Giving up after %d retries on URL "%s"...', maxRetries, url);

    else

      % First time only(?)
      if ( ~isfield(dat,'sample_url') )
        dat.sample_url = url;
        dat.lon = double( nc{'lon'}(1:end,1:end,1:end,1:end) );
        dat.lon(dat.lon > 180) = dat.lon(dat.lon > 180) - 360;
        dat.lat = double( nc{'lat'}(1:end,1:end,1:end,1:end) );
        % Note: ROWS are LATITUDE, not LONGITUDE! Dyslexics of the world untie...
        dat.nx = length(dat.lat);
        dat.ny = length(dat.lon);
        % Preallocate like a good MATLAB dooby
        for fldix = 1:length(dat.flds)
          dat.(dat.flds{fldix}) = repmat(nan, [numel(dat.date) dat.nx dat.ny]);
        end;
      end;

      % How the heck do you query dimensionality in netCDF-Java??
      x = double( nc{var}(1:end,1:end,1:end,1:end) );

      % Which dimension is time? (CFSR producers added an extra value for,
      % e.g., 01-Feb 00:00 into their "January" files, etc. Nitwits.)
      t = double( nc{'time'}(1:end,1:end) );
      nTimes = numel(t);
      timeDim = find(size(x) == nTimes);
      switch (timeDim)
       case 1, x = double( nc{var}(1:end-1,1:end,1:end,1:end) );
       case 2, x = double( nc{var}(1:end,1:end-1,1:end,1:end) );
       case 3, x = double( nc{var}(1:end,1:end,1:end-1,1:end) );
       case 4, x = double( nc{var}(1:end,1:end,1:end,1:end-1) );
      end;

      if ( numel(x) < (numel(dat.date)*dat.nx*dat.ny) )
        warning('Insufficient data for "%s" ("%s") from "%s"!', var, fld, url);

      else
        y = [];
        switch ( ndims(x) ),
         case 3,
          % We always want the first and nearest the surface
          y = double( squeeze(x(:,:,:)) );
         case 4,
          % We always want the first and nearest the surface
          y = double( squeeze(x(:,:,:,:)) );
         otherwise,
          warning('Unexpected NDIMS %d for %s (%s)!', ndims(x), var, fld);
        end;
        if ( ~isempty(y) )
          dat.(fld)(:,:,:) = y;
        end;
      end; %if numel(x)

      close(nc);

    end; %if isempty(nc) else

  end; %for vix

  %DEBUG:
  toc,

return;
