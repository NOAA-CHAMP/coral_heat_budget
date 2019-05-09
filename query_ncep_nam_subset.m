function dat = query_ncep_narr_subset(bbox,begdt,enddt)
%function dat = query_ncep_narr_subset(bbox,begdt,enddt)
%
% Extract a subset of the NCEP North American Mesoscale Reanalysis (NAM) from
% the NOMADS THREDDS site at ncdc.noaa.gov. Optional arg BBOX will spatially
% subset the data (DEFAULT: Florida Reef Tract region BBOX=[24 27 -83 -79]).
% Optional args BEGDT ENDDT, if present, temporally subset the data.
%
% CALLS: netCDF-Java Toolbox (matlab-njTbx), via its interface that is like
% old NetCDF Toolbox by Chuck Denham: http://sourceforge.net/apps/trac/njtbx
%
% Last Saved Time-stamp: <Thu 2010-08-05 12:57:16  Lew.Gramer>

error('WAIT! Add air and sea temperature to this code before rerunning!');

  %DEBUG:
  tic,

  if ( ~exist('bbox','var') || isempty(bbox) )
    bbox = [ 24 27 -83 -79 ];
  end;

  if ( ~exist('begdt','var') || isempty(begdt) )
    begdt = datenum(1987,02,01);
  end;
  if ( ~exist('enddt','var') || isempty(enddt) )
    %enddt = datenum(2009,12,31);
    %%%% FOR TESTING - just get one day
    enddt = datenum(1987,02,01);
  end;

  % Maximum number of times to retry each URL (i.e., each NOMADS time)
  maxRetries = 3;

  ALLVARS = { ...
      'Albedo', ...
      'Dew_point_temperature', ...
      'Downward_long_wave_flux', ...
      'Downward_short_wave_flux', ...
      'Latent_heat_flux', ...
      'Pressure', ...
      'Relative_humidity', ...
      'Sensible_heat_flux', ...
      'Friction_velocity', ...
      'Total_cloud_cover', ...
      'Total_precipitation', ...
      'Upward_long_wave_flux', ...
      'Upward_short_wave_flux', ...
      'u_wind_height_above_ground', ...
      'v_wind_height_above_ground', ...
            };


  dat.bbox = bbox;
  dat.vars = ALLVARS;

  query_vars = sprintf('&var=%s',dat.vars{:});
  query_vars(1) = [];

  dat.date = floor(begdt):(6/24):(floor(enddt)+1-(1/(24*60)));

  [yrs,mns,dys,hrs] = datevec(dat.date);

  for dix = 1:length(dat.date)

    dt = dat.date(dix);
    %DEBUG:    disp(datestr(dt));

    %http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/198701/19870131/narr-a_221_19870131_2100_000.grb?var=Albedo&var=Surface_friction_velocity&var=Pressure&var=Specific_humidity_height_above_ground&var=Dew_point_temperature&var=Relative_humidity&var=u_wind_height_above_ground&var=v_wind_height_above_ground&var=Downward_longwave_radiation_flux&var=Downward_shortwave_radiation_flux&var=Evaporation&var=Latent_heat_flux&var=Sensible_heat_flux&var=Total_cloud_cover&var=Total_precipitation&var=Upward_long_wave_radiation_flux_surface&var=Upward_short_wave_radiation_flux_surface&latitude=25.01&longitude=-80.38&temporal=all&spatial=bb&north=27&south=24&west=-83&east=-79&east=-79&addLatLon=true

    url = sprintf( 'http://nomads.ncdc.noaa.gov/thredds/ncss/grid/nam/%04d%02d/%04d%02d%02d/nam_218_%04d%02d%02d_%02d00_000.grb?%s&temporal=all&spatial=bb&south=%g&north=%g&west=%g&east=%g&addLatLon=true', ...
                   yrs(dix), mns(dix), yrs(dix), mns(dix), dys(dix), ...
                   yrs(dix), mns(dix), dys(dix), hrs(dix), query_vars, ...
                   bbox(1), bbox(2), bbox(3), bbox(4) ...
                   );
    %DEBUG:    url,


    nc = [];
    tried = 0;
    while ( isempty(nc) && (tried < maxRetries) )
      % This is in here to keep from overheating NCDC's NOMAD servers!
      pause(0.5);
      nc = mDataset(url);
      tried = tried + 1;
    end;

    if ( isempty(nc) )
      warning('Giving up after %d retries on URL "%s"...', maxRetries, url);

    else

      % First time only
      if ( ~isfield(dat,'sample_url') )
        dat.sample_url = url;
        dat.lon = nc{'lon'}(1:end,1:end,1:end,1:end);
        dat.lat = nc{'lat'}(1:end,1:end,1:end,1:end);
        [dat.nx dat.ny] = size(dat.lon);
        % Preallocate like a good MATLAB dooby
        for vix = 1:length(dat.vars)
          dat.(dat.vars{vix}) = repmat(nan, [numel(dat.date) dat.nx dat.ny]);
        end;
      end;

      for vix = 1:length(dat.vars)
        % How the heck do you query dimensionality in netCDF-Java??
        x = nc{dat.vars{vix}}(1:end,1:end,1:end,1:end);
        if ( numel(x) < (dat.nx*dat.ny) )
          warning('Insufficient data for "%s" from "%s"!', dat.vars{vix}, url);

        else
          switch ( ndims(x) ),
           case 2,
            dat.(dat.vars{vix})(dix,:,:) = x;
           case 3,
            % We always want the first and nearest the surface
            dat.(dat.vars{vix})(dix,:,:) = squeeze(x(1,:,:));
           case 4,
            % We always want the first and nearest the surface
            dat.(dat.vars{vix})(dix,:,:) = squeeze(x(1,1,:,:));
          end;
        end; %if numel(x)

      end; %for vix

      close(nc);

    end; %if isempty(nc) else

  end; %for dix

  %DEBUG:
  toc,

return;
