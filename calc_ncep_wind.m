function stn = calc_ncep_wind(stn,dataset)
%function stn = calc_ncep_wind(stn,dataset)
%
% Calc wind speed (kts.) and direction from NCEP wind U and V
%
% Last Saved Time-stamp: <Sun 2010-03-14 14:34:35 Eastern Daylight Time gramer>

  if ( ~exist('dataset','var') || isempty(dataset) );
    dataset = '';
  elseif ( dataset(end) ~= '_' )
    dataset = [dataset '_'];
  end;

  ufld = ['ncep_' dataset 'wind_u'];
  vfld = ['ncep_' dataset 'wind_v'];
  spdfld = ['ncep_' dataset 'wind_speed'];
  dirfld = ['ncep_' dataset 'wind_dir'];

  if (isfield(stn,spdfld)); stn = rmfield(stn, spdfld); end;
  if (isfield(stn,dirfld)); stn = rmfield(stn, dirfld); end;

  [ix1,ix2] = intersect_dates(stn.(ufld).date, stn.(vfld).date);
  stn.(spdfld).date = stn.(ufld).date(ix1);
  stn.(spdfld).data = uv_to_spd(stn.(ufld).data(ix1), stn.(vfld).data(ix2));
  % Convert wind speed to knots (because we'll convert it back in HFBULKTC!)
  stn.(spdfld).data = stn.(spdfld).data ./ 0.5144444444;

  stn.(dirfld).date = stn.(ufld).date(ix1);
  stn.(dirfld).data = uv_to_dir(stn.(ufld).data(ix1), stn.(vfld).data(ix2));

return;
