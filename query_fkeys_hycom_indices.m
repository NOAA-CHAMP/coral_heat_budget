function [yix,xix] = query_fkeys_hycom_indices(lon,lat)
%function [yix,xix] = query_gom_hycom_indices(lon,lat)
%
% Return y- and x-indices (zero-based) for Florida Keys 1/100-degree HYCOM

  global fkeys_hycom_lons fkeys_hycom_lats

  datapath = get_thesis_path('../data');

  if ( isempty(fkeys_hycom_lons) || isempty(fkeys_hycom_lats) )
    load(fullfile(datapath, 'fkeys_hycom_coords.mat'));
  end;

  [dx,xix] = min( abs(fkeys_hycom_lons - lon) );
  [dy,yix] = min( abs(fkeys_hycom_lats - lat) );

  if ( dx >= 0.02 || dy >= 0.02 )
    error('Ecoforecasts:fkeysHYCOM:BadCoords',...
          'Lon/Lat %f/%f outside FKEYS HYCOM domain!',lon,lat);
  end;

return;
