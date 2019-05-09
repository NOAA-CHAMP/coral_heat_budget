function [yix,xix] = query_fkeys_indices(lon,lat)

  persistent FKEYS_LATS FKEYS_LONS

  warning('This function has been SUPERCEDED by QUERY_FKEYS_HYCOM_INDICES!');

  datapath = get_thesis_path('../data');

  if ( isempty(FKEYS_LATS) )
    load(fullfile(datapath, 'fkeys_coords.mat'));
  end;

  [dx,xix] = min( abs(FKEYS_LONS - lon) );
  [dy,yix] = min( abs(FKEYS_LATS - lat) );

  if ( dx >= 0.02 || dy >= 0.02 )
    error('Lon/Lat %f/%f outside FKEYS domain!',lon,lat);
  end;

  xix = xix - 1;
  yix = yix - 1;

return;
