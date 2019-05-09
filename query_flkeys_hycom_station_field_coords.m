function stn = query_flkeys_hycom_station_field_coords(stn,fldnm)
%function stn = query_flkeys_hycom_station_field_coords(stn,fldnm)
%
% Add lat and lon Nx1 fields for Florida Keys 1/100-degree HYCOM to
% STN.(FLDNM), and rename field STN.(FLDNM).data to STN.(FLDNM).field
% instead. DEFAULT FLDNM: 'flkeys_hycom_seatemp_field'.
%
% Last Saved Time-stamp: <Thu 2010-10-07 17:05:45 Eastern Daylight Time gramer>

  datapath = get_thesis_path('../data');

  if ( ~exist('fldnm','var') || isempty(fldnm) )
    fldnm = 'flkeys_hycom_seatemp_field';
  end;

  % Depth[7]
  % 0.0, 5.0, 10.0, 30.0, 50.0, 75.0, 100.0

  minlon = -83.36;
  maxlon = -79.00;
  dlon = 0.0100;
  lons = minlon:dlon:maxlon;
  xix = round( (stn.lon - minlon) ./ dlon );

  % minlat = 18.0916;
  % dlat = 0.0092;  ... dlat = 0.0090;
  % yix = round( (lat - minlat) ./ dlat );
  load(fullfile(datapath, 'flkeys_hycom_lats.mat'));
  [ig,yix] = min( abs(flkeys_hycom_lats - stn.lat) );
  yix = yix - 1;


  [tsz,ysz,xsz] = size(stn.(fldnm).data);
  xrad = floor(xsz/2);
  yrad = floor(ysz/2);

  stn.(fldnm).lon = lons(xix-xrad:xix+xrad)';
  stn.(fldnm).lat = flkeys_hycom_lats(yix-yrad:yix+yrad)';

  stn.(fldnm).field = stn.(fldnm).data;
  stn.(fldnm) = rmfield(stn.(fldnm),'data');

return;
