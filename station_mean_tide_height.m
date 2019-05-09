function stn = station_mean_tide_height(stn,mhfld,bfld,hfld,ufld,vfld)
%function stn = station_mean_tide_height(stn,mhfld,bfld,hfld,ufld,vfld)
%
% Calculate the mean tide depth experienced by a watermass moving over
% an M2 tidal ellipse centered on the coordinates of station struct STN.
% Time series is returned in STN.(MHFLD), based on bathymetry STN.(BFLD),
% tide time series STN.(HFLD), and *optional* u and v component time series
% for tidal velocity, STN.(UFLD) and STN.(VFLD). (DEFAULT: 0.06 m/s.) UFLD
% may also be a *numeric scalar* default velocity to use. If BFLD is numeric
% scalar, use that value for the mean seafloor depth.
%
% Last Saved Time-stamp: <Tue 2012-12-04 14:07:28 Eastern Standard Time lew.gramer>


  if ( ~exist('ufld','var') || isempty(ufld) )
    ufld = [];
  end;
  if ( ~exist('vfld','var') || isempty(vfld) )
    vfld = [];
  end;


  %%%
  %% First estimate a mean cross-shore tidal velocity "somehow"

  if ( isscalar(ufld) && isnumeric(ufld) )
    % Caller passed in a guess at the mean velocity
    mean_flood_tide_cross_shore_vel = ufld;

  elseif ( ~isfield(stn,ufld) && ~isfield(stn,vfld) )
    % Good guess for "shelf break" of Florida Keys, per these references:
    %  Shay et al (1998): Effects of low-frequency current variability on
    %                     near-inertial submesoscale vortices
    %  Haus et al (2000): Remote radar measurement of shelf currents off Key
    %                     Largo Florida USA
    %  Haus et al (2004): Southeast Florida Shelf circulation and volume
    %                     exchange, observations of km-scale variability 
    % This estimate also agrees well with 2 weeks of depth-averaged cross
    % shore velocities at SEAKEYS/GLERL/ICON station TNRF1 in July 2010.
    mean_flood_tide_cross_shore_vel = 0.06;

  elseif ( ~isfield(stn,vfld) )
    % Caller passed in only one velocity: U is ALREADY cross-shore component
    mean_flood_tide_cross_shore_vel = iqr(stn.(ufld).data)/2;

  elseif ( ~isfield(stn,ufld) )
    % Caller passed in only one velocity: V is ALREADY cross-shore component
    mean_flood_tide_cross_shore_vel = iqr(stn.(vfld).data)/2;

  else  
    % Calculate cross-shore velocity from tidal U and V velocity estimates

    % Azimuth of lowest local relief relative to True North (0=straight N-S coast)
    orifld = 'isobath_orientation';
    if ( ~isfield(stn,orifld) )
      %DEBUG:
      warning('Trying in-code isobath orientation estimate');
      switch (lower(stn.station_name)),
       case 'fwyf1',	stn.(orifld)=10;
       case 'mlrf1',	stn.(orifld)=65;
       case 'lonf1',	stn.(orifld)=70;
       case 'tnrf1',	stn.(orifld)=65;
       case 'smkf1',	stn.(orifld)=70;
       case 'looe1',	stn.(orifld)=80;  % Or 70 as in GET_LOOE1_ADCP???
       otherwise,	error('No isobath orientation for station %s',stn.station_name);
      end;
    end;

    u = stn.(ufld).data;
    v = stn.(vfld).data;
    xshore = (cosd(stn.(orifld))*u) - (sind(stn.(orifld))*v);
    lshore = (sind(stn.(orifld))*u) + (cosd(stn.(orifld))*v);
    mean_flood_tide_cross_shore_vel = iqr(xshore)/2;
  end;

  %DEBUG:
  disp([mfilename ': Cross-shore Tidal Current ' num2str(mean_flood_tide_cross_shore_vel)]);

  %%%
  %% Now estimate mean bathymetry across one M2 tidal excursion (period/4)

  % Rationale for M2 vs other constituents (e.g., K1) is that we wish to know
  % mean excursion over multiple diurnal cycles, so that semi- and diurnal
  % constituents have small relative coherence. Also the modulation of the
  % *diurnal* warming cycle (the point of this exercise) will tend to be
  % dominated by higher-frequency tidal excursions of 3 hours or less.

  tidal_excursion = mean_flood_tide_cross_shore_vel*(3.1*3600);

  if ( isscalar(bfld) && isnumeric(bfld) )
    mean_bathy_depth = bfld;

  elseif ( ~ischar(bfld) )
    error('BFLD must be scalar numeric, or fieldname string');

  elseif ( ~isfield(stn,bfld) )
    error('Bathymetry field STN.%s not found',bfld);

  else
    % Resolution of the bathymetry in [m]
    bathy_res = min([min(diff(unique(stn.(bfld).lon(:)))),min(diff(unique(stn.(bfld).lat(:))))]).*111e3;

    radx = ceil(tidal_excursion/bathy_res/sqrt(2));

    % %%%% DEBUG???
    % % FOR NOW, SIMPLY TAKE THE MEAN OF THE 9 POINTS SURROUNDING A SITE
    % radx = 1;

    [ig,bathix] = min(sqrt(((stn.(bfld).lon(:)-stn.lon).^2)+((stn.(bfld).lat(:)-stn.lat).^2)));
    [ix,jx] = ind2sub(size(stn.(bfld).field),bathix);

    mean_bathy_depth = -mean(mean(stn.(bfld).field(ix-radx:ix+radx,jx-radx:jx+radx)));
  end;

  %%%
  %% Finally estimate tidal time series based on tide-cycle mean bathymetry

  % Assume tide height already includes our estimate of site depth: subtract
  % from mean before adding, and also add the error in bathymetry estimate of
  % site depth relative to ours, to compensate for overall error in mean.
  % site_bathy_depth = -stn.(bfld).field(ix,jx);
  if ( isfield(stn,hfld) )
    % Nothing to do
  elseif ( isnumeric(hfld) )
    val = hfld;
    hfld = 'fixed_tide_depth';
    stn.(hfld).date = now;
    stn.(hfld).data = val;
  else
    error('HFLD must either be a field name or numeric value');
  end;

  site_bathy_depth = nanmean(stn.(hfld).data);
  stn.(mhfld).date = stn.(hfld).date;
  stn.(mhfld).data = stn.(hfld).data + mean_bathy_depth - site_bathy_depth;

  %DEBUG:  mean_flood_tide_cross_shore_vel, tidal_excursion, mean_bathy_depth, site_bathy_depth,
  %DEBUG:
  disp([mfilename,': Mean Water Depth experienced by water parcel ',num2str(nanmean(stn.(mhfld).data))]);

return;
