function dat = query_ncep_gfs_subset(bbox,begdt,enddt)


error('FUNCTION CURRENTLY NOT READY TO BE USED!');

error('And WAIT! Add air and sea temperature to this code before rerunning!');


%function dat = query_ncep_gfs_subset(bbox,begdt,enddt)
%
% Extract a subset of the NCEP Global Forecast System (GFS) directly from
% the NOMADS THREDDS site at ncdc.noaa.gov. Optional arg BBOX will spatially
% subset the data (DEFAULT: American Samoa BBOX=[-15 -14 -171 -170]).
% Optional args BEGDT ENDDT, if present, temporally subset the data.
%
% CALLS: netCDF-Java Toolbox (matlab-njTbx), via its interface that is like
% old NetCDF Toolbox by Chuck Denham: http://sourceforge.net/apps/trac/njtbx
%
% Last Saved Time-stamp: <Thu 2010-08-05 12:57:39  Lew.Gramer>

  %DEBUG:
  tic,

  if ( ~exist('bbox','var') || isempty(bbox) )
    bbox = [-15 -14 -171 -170];
  end;

  if ( ~exist('begdt','var') || isempty(begdt) )
    begdt = datenum(2003,09,01);
  end;
  if ( ~exist('enddt','var') || isempty(enddt) )
    %enddt = datenum(2009,12,31);
    %%%% FOR TESTING - just get one day
    enddt = datenum(2003,09,02);
  end;

  ALLVARS = { ...
      'U-component_of_wind_height_above_ground', ...
      'V-component_of_wind_height_above_ground', ...
      'Relative_humidity_height_above_ground', ...
            };


  dat.bbox = bbox;
  dat.vars = ALLVARS;

  query_vars = sprintf('&var=%s',dat.vars{:});
  query_vars(1) = [];

  dat.date = floor(begdt):(3/24):(floor(enddt)+1-(1/(24*60)));

  [yrs,mns,dys,hrs] = datevec(dat.date);

  for dix = 1:length(dat.date)

    dt = dat.date(dix);
    %DEBUG:    disp(datestr(dt));

    %http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/198701/19870131/narr-a_221_19870131_2100_000.grb?var=Albedo&var=Surface_friction_velocity&var=Pressure&var=Specific_humidity_height_above_ground&var=Dew_point_temperature&var=Relative_humidity&var=u_wind_height_above_ground&var=v_wind_height_above_ground&var=Downward_longwave_radiation_flux&var=Downward_shortwave_radiation_flux&var=Evaporation&var=Latent_heat_flux&var=Sensible_heat_flux&var=Total_cloud_cover&var=Total_precipitation&var=Upward_long_wave_radiation_flux_surface&var=Upward_short_wave_radiation_flux_surface&latitude=25.01&longitude=-80.38&temporal=all&spatial=bb&north=27&south=24&west=-83&east=-79&east=-79&addLatLon=true

    %http://nomads.ncdc.noaa.gov/thredds/ncss/grid/gfs4/201002/20100203/gfs_4_20100203_0600_000.grb2?var=U-component_of_wind_height_above_ground&var=V-component_of_wind_height_above_ground&temporal=point&accept=netcdf&time=2010-02-03T06:00:00Z&point=true&latitude=-14&longitude=-171

    url = sprintf( 'http://nomads.ncdc.noaa.gov/thredds/ncss/grid/gfs4/%04d%02d/%04d%02d%02d/gfs_4_%04d%02d%02d_%02d00_000.grb?%s&temporal=point&accept=netcdf&time=   ????    &spatial=bb&south=%g&north=%g&west=%g&east=%g&addLatLon=true', ...
                   yrs(dix), mns(dix), yrs(dix), mns(dix), dys(dix), ...
                   yrs(dix), mns(dix), dys(dix), hrs(dix), query_vars, ...
                   bbox(1), bbox(2), bbox(3), bbox(4) ...
                   );
    %DEBUG:    url,

    nc = mDataset(url);

    % First time only
    if ( dix == 1 )
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
    end;

    close(nc);

  end;

  %DEBUG:
  toc,

return;
