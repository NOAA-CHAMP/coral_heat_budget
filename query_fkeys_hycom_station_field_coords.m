function stn = query_fkeys_hycom_station_field_coords(stn,fldnm)
%function stn = query_fkeys_hycom_station_field_coords(stn,fldnm)
%
% Add lat and lon Nx1 fields for FKEYS 1/100-degree HYCOM to
% STN.(FLDNM). DEFAULT FLDNM: 'fkeys_hycom_seatemp_field'.
%
% Last Saved Time-stamp: <Fri 2011-03-18 14:10:50  lew.gramer>
%
% error('THIS FUNCTION IS NOW DEFUNCT! See GET_FKEYS_HYCOM.m...');

error('THIS FUNCTION IS NOW DEFUNCT! See GET_FKEYS_HYCOM.m...');

  global fkeys_hycom_lons fkeys_hycom_lats

  datapath = get_thesis_path('../data');

  if ( ~exist('fldnm','var') || isempty(fldnm) )
    fldnm = 'fkeys_hycom_seatemp_field';
  end;

  if ( isempty(fkeys_hycom_lons) || isempty(fkeys_hycom_lats) )
    load(fullfile(datapath, 'fkeys_hycom_coords.mat'));
  end;


  lons = fkeys_hycom_lons(:);
  minlon = min(lons(:));
  dlon = min(abs(diff(lons(:))));
  maxlon = max(lons(:));
  if ( minlon-dlon > stn.lon || stn.lon > maxlon+dlon )
    error('Ecoforecasts:fkeysHYCOM:BadLon',...
          'Longitude %f outside value range [%f,%f]',stn.lon,minlon,maxlon);
  end;

  [ig,xix] = min( abs(lons(:) - stn.lon) );


  lats = fkeys_hycom_lats(:);
  minlat = min(lats(:));
  dlat = min(abs(diff(lats(:))));
  maxlat = max(lats(:));
  if ( minlat-dlat > stn.lat || stn.lat > maxlat+dlat )
    error('Ecoforecasts:fkeysHYCOM:BadLat',...
          'Latitude %f outside value range [%f,%f]',stn.lat,minlat,maxlat);
  end;

  [ig,yix] = min( abs(lats(:) - stn.lat) );


  [tsz,ysz,xsz] = size(stn.(fldnm).field);

  begx = xix-floor(xsz/2);
  endx = begx+xsz-1;
  begx = max(begx,1);
  endx = min(endx,length(lons));

  begy = yix-floor(ysz/2);
  endy = begy+ysz-1;
  begy = max(begy,1);
  endy = min(endy,length(lats));

  stn.(fldnm).lon = lons(begx:endx);
  stn.(fldnm).lat = lats(begy:endy);

return;
