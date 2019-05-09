function [yix,xix] = query_gom_hycom_indices(lon,lat)
%function [yix,xix] = query_gom_hycom_indices(lon,lat)
%
% Return y- and x-indices (zero-based) for Gulf of Mexico 1/25-degree HYCOM

  datapath = get_thesis_path('../data');

  minlon = -98;
  dlon = 0.0400;
  maxlon = -76.40;
  if ( minlon-dlon > lon || lon > maxlon+dlon )
    error('Ecoforecasts:gomHYCOM:BadLon',...
          'Longitude %f outside value range [%f,%f]',lon,minlon,maxlon);
  end;
  xix = round( (lon - minlon) ./ dlon );
  % Dumb-@ss3d one-based FORTRAN/MATLAB! Shizzle...
  xix = xix + 1;

  % minlat = 18.0916;
  % dlat = 0.0339;  ... dlat = 0.0381;
  % yix = round( (lat - minlat) ./ dlat );
  load(fullfile(datapath, 'gom_hycom_lats.mat'));
  minlat = min(gom_hycom_lats(:));
  dlat = min(abs(diff(gom_hycom_lats(:))));
  maxlat = max(gom_hycom_lats(:));
  if ( minlat-dlat > lat || lat > maxlat+dlat )
    error('Ecoforecasts:gomHYCOM:BadLat',...
          'Latitude %f outside value range [%f,%f]',lat,minlat,maxlat);
  end;
  [ig,yix] = min( abs(gom_hycom_lats - lat) );

return;
