function stn = query_gom_hycom_station_field_coords(stn,fldnm)
%function stn = query_gom_hycom_station_field_coords(stn,fldnm)
%
% Add lat and lon Nx1 fields for Gulf of Mexico 1/25-degree HYCOM to
% STN.(FLDNM). DEFAULT FLDNM: 'gom_hycom_seatemp_field'.
%
% Last Saved Time-stamp: <Sun 2011-03-27 17:30:20  lew.gramer>
%
% error('THIS FUNCTION IS NOW DEFUNCT! See GET_GOM_HYCOM.m...');

error('THIS FUNCTION IS NOW DEFUNCT! See GET_GOM_HYCOM.m...');

  datapath = get_thesis_path('../data');

  if ( ~exist('fldnm','var') || isempty(fldnm) )
    fldnm = 'gom_hycom_seatemp_field';
  end;

  minlon = -98.00;
  dlon = 0.0400;
  maxlon = -76.40;
  if ( minlon-dlon > stn.lon || stn.lon > maxlon+dlon )
    error('Ecoforecasts:gomHYCOM:BadLon',...
          'Longitude %f outside value range [%f,%f]',stn.lon,minlon,maxlon);
  end;

  lons = minlon:dlon:maxlon;

  % xix = round( (stn.lon - minlon) ./ dlon ); xix = xix + 1;
  [ig,xix] = min( abs(lons - stn.lon) );


  % minlat = 18.0916;
  % dlat = 0.0339;  ... dlat = 0.0381;
  % yix = round( (lat - minlat) ./ dlat ); yix = yix + 1;
  load(fullfile(datapath, 'gom_hycom_lats.mat'));
  minlat = min(gom_hycom_lats(:));
  dlat = min(abs(diff(gom_hycom_lats(:))));
  maxlat = max(gom_hycom_lats(:));
  if ( minlat-dlat > stn.lat || stn.lat > maxlat+dlat )
    error('Ecoforecasts:gomHYCOM:BadLat',...
          'Latitude %f outside value range [%f,%f]',stn.lat,minlat,maxlat);
  end;

  [ig,yix] = min( abs(gom_hycom_lats - stn.lat) );


  [tsz,ysz,xsz] = size(stn.(fldnm).field);
  % begx = xix-floor(xsz/2);	begx = max(begx,1);
  % endx = xsz-begx+1;		endx = min(endx,xsz);
  begx = xix-floor(xsz/2);
  endx = begx+xsz-1;
  begx = max(begx,1);
  endx = min(endx,length(lons));

  % begy = yix-floor(ysz/2);	begy = max(begy,1);
  % endy = ysz-begy+1;		endy = min(endy,ysz);
  begy = yix-floor(ysz/2);
  endy = begy+ysz-1;
  begy = max(begy,1);
  endy = min(endy,length(gom_hycom_lats));

  stn.(fldnm).lon = lons(begx:endx)';
  stn.(fldnm).lat = gom_hycom_lats(begy:endy)';

return;
