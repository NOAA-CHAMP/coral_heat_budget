function [U, V] = cwave(LON, LAT, U, V, maxV, lon, lat, lambdakm)
%function [U, V] = cwave(LON, LAT, U, V, maxV, lon, lat, lambdakm)
% 
% Generate an internal wave of horizontal wavelength 'lambdakm', on the given
% LON/LAT grid. Generates velocities up to a maximum maxV [m/s], normal to
% the waveguide contained in vectors 'lon' and 'lat' (optionally specifying
% just two end points on the grid, or a single meridional or zonal
% coordinate), and perturbs a mean flow field [U,V] (optionally a pair of
% empty matrices) with those u and v velocities.
% 

  if ( isempty(U) )
    U = repmat(nan, size(LON));
  end;
  if ( isempty(V) )
    V = repmat(nan, size(LON));
  end;
  lambdam = lambdakm * 1e3;
  kappam = 1/lambdam;

  dlon = min(diff(LON(1,:)));
  dlat = min(diff(LAT(:,1)));

  if ( numel(lon) == 1 && numel(lat) == 0 )
    lon(2) = lon;
    lat = [LAT(1) LAT(end)];
  elseif ( numel(lon) == 0 && numel(lat) == 1 )
    lon = [LON(1) LON(end)];
    lat(2) = lat;
  end;

  Y = lon + (i .* lat);
  XI = linspace(1, length(Y), 1000);
  y = interp1(Y, XI, 'pchip');
  lonl = real(y);
  latl = imag(y);
  lon = round(lonl ./ (dlon/2)) .* (dlon/2);
  lat = round(latl ./ (dlat/2)) .* (dlat/2);
  uniqpts = find(diff(lon) | diff(lat));
  lon = lon([1 uniqpts]);
  lat = lat([1 uniqpts]);

  pts = [lon ; lat];
  sct = pts(:,3:end) - pts(:,1:end-2);
  nrm = [sct(2,:) ; -sct(1,:)];


  % Defaults: nautical miles, and degrees CCW from East. WTF, CSIRO??
  [rngkm, ign] = sw_dist(pts(2,:), pts(1,:), 'km');
  pathlenm = [0 (cumsum(rngkm) .* 1e3)];

  Vmult = maxV * sin(2*pi*(pathlenm ./ lambdam));

  for ix = 1:length(lon)
    ptix = find(abs(LON - lon(ix)) <= (dlon/2) ...
                & abs(LAT - lat(ix)) <= (dlat/2));
    [latix, lonix] = ind2sub(size(LON), ptix);
    U(latix, lonix) = U(lonix, latix) + Vmult(ix);
    for subix = 1:size(U,2)
      rngkm = sw_dist([lat(ix) lat(ix)], [lon(ix) LON(latix,subix)], 'km');
      rng = rngkm * 1e3;
      if ( rng <= (lambdam/4) )
        Vl = Vmult(ix) * abs(cos(2*pi*rng/lambdam));
        U(latix, subix) = U(latix, subix) + Vl;
      end;
    end;
  end;

return;
