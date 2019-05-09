function [yix,xix] = query_glb_analysis_indices(lon,lat)

  datapath = get_thesis_path('../data');

  lon(lon < 0) = 360 + lon(lon < 0);

  minlon = 74.160034;
  dlon = 0.0800;
  xix = round( (lon - minlon) ./ dlon );

  % minlat = -78.64;
  % dlat = 0.0320;
  % yix = round( (lat - minlat) ./ dlat );
  load(fullfile(datapath, 'glb_analysis_lats.mat'));
  [ig,yix] = min( abs(glb_analysis_lats - lat) );
  yix = yix - 1;

return;
