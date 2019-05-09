function [udotdelT,dTdx_1,dTdy_1] = calc_udotdelt(u,v,tstr,fldix,interpMethod,lat,lon)
%function [udotdelT,dTdx_1,dTdy_1] = calc_udotdelt(u,v,tstr,fldix,[interpMethod,lat,lon])
%
% Calculate vector dot product between each ordinate of ocean current/wind
% vector U and V and the 2-D field in struct TSTR (which must have fields
% .lon, .lat, .date, and .field). Returns a vector UDOTDELT of dot products
% at central point of TSTR.field, and time series DTDX_1 and DTDY_1 of the X-
% and Y- components of field GRADIENT (v.) for site - all for each date. If
% index or logical vector FLDIX is given, use subset tstr.field(FLDIX,:,:).
% If struct fields TSTR.gradient_x,TSTR.gradient_y already exist, use them.
% Optional INTERPMETHOD (DEFAULT 'midpoint'): if 'midpoint', choose point in
% each gradient matrix nearest the center; otherwise, pass to INTERP_FIELD
% with interpolation coordinates, LAT and LON (which then must be present).
%
% Last Saved Time-stamp: <Thu 2011-05-05 16:07:35  lew.gramer>

  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'midpoint';
  end;

  % Sanity checks
  if ( ~isfield(tstr,'field') || ~isfield(tstr,'lon') || ~isfield(tstr,'lat') )
    error('Third arg should be a struct with .lon, .lat, .field fields!');
  end;
  if ( ~exist('fldix','var') || isempty(fldix) )
    fldix = 1:size(tstr.field,1);
  end;
  fld = tstr.field(fldix,:,:);

  if ( size(u,1) ~= size(v,1) || size(u,1) ~= size(fld,1) )
    error('U, V, and TSTR.field(FLDIX) must all cover equal times (size(:,1))!');
  end;

  midx = round(length(tstr.lon) / 2);
  midy = round(length(tstr.lat) / 2);

  if ( isfield(tstr,'gradient_x') && isfield(tstr,'gradient_y') )
    dTdx = tstr.gradient_x(fldix,:,:);
    dTdy = tstr.gradient_y(fldix,:,:);
  else

    % Inter-gridpoint distances in [m] for gradient calculation
    midlon = repmat(tstr.lon(midy),size(tstr.lat));
    midlat = repmat(tstr.lat(midy),size(tstr.lon));
    dt = tstr.date(fldix);
    dx = [0 ; cumsum(sw_dist(midlat,tstr.lon,'km'))] .* 1e3;
    dy = [0 ; cumsum(sw_dist(tstr.lat,midlon,'km'))] .* 1e3;

    %DEBUG:    disp('Calculating GRADIENT');
    [dTdx,dTdy,dTdt] = gradient(permute(fld,[2 3 1]),dx,dy,dt);
    % [dTdx,dTdy,dTdt] = gradientn(permute(fld,[2 3 1]),5,dx,dy,dt);
    dTdx = permute(dTdx,[3 1 2]);
    dTdy = permute(dTdy,[3 1 2]);

  end;

  if ( strncmpi(interpMethod,'midpoint',3) )
    % We are only interested in the central point
    dTdx_1 = squeeze(dTdx(:,midx,midy));
    dTdy_1 = squeeze(dTdy(:,midx,midy));
  else
    dTdx_1 = interp_field(tstr.lat,tstr.lon,dTdx,lat,lon,interpMethod);
    dTdy_1 = interp_field(tstr.lat,tstr.lon,dTdy,lat,lon,interpMethod);
  end;

  % Field gradient time series at our point of interest
  delT = [dTdx_1 , dTdy_1]';

  % Advective vector at our point of interest
  U = [u , v]';

  udotdelT = dot(delT, U);

return;
