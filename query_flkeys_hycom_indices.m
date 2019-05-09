function [yix,xix] = query_flkeys_hycom_indices(lon,lat)
%function [yix,xix] = query_flkeys_hycom_indices(lon,lat)
%
% Return y- and x-indices (zero-based) for Florida Keys 1/100-degree HYCOM

  datapath = get_thesis_path('../data');

  minlon = -83.36;
  dlon = 0.0100;
  xix = round( (lon - minlon) ./ dlon );
  % Dumb-@ss3d one-based FORTRAN/MATLAB indexing! Shizzle...
  xix = xix + 1;

  % minlat = 22.775372;
  % dlat = ??? 0.0339;  ... dlat = 0.0381;
  % yix = round( (lat - minlat) ./ dlat );
  load(fullfile(datapath, 'flkeys_hycom_lats.mat'));
  [ig,yix] = min( abs(flkeys_hycom_lats - lat) );

return;
