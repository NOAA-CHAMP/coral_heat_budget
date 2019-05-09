function stn = estim_ngdc_depth(stn_or_stnm)
%function stn = estim_ngdc_depth(stn_or_stnm)
%
% Find closest gridpoint(s) to location of station named by or contained in
% STN_OR_STNM (string or struct), within the NGDC 3-arcsecond Coastal Relief
% Model; store the median of these values in STN.NGDC_DEPTH for that site.
%
% Last Saved Time-stamp: <Sun 2010-06-13 14:58:06 Eastern Daylight Time gramer>

  if ( ischar(stn_or_stnm) )
    stn.station_name = stn_or_stnm;
  elseif ( isstruct(stn_or_stnm) )
    stn = stn_or_stnm;
  end;
  if ( ~isfield(stn,'station_name') )
    error('No station (5-char) name in specified argument');
  end;

  if ( ~isfield(stn,'lon') )
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
  end;
  if ( ~isfield(stn,'ngdc_92m_bathy') )
    stn = load_ngdc_bathy(stn,true);
  end;

  lons = stn.ngdc_92m_bathy.lon - stn.lon;
  lats = stn.ngdc_92m_bathy.lat - stn.lat;
  % Find all points within one gridpoint of our station
  ix = find(sqrt((lons.^2) + (lats.^2)) < (3.1/3600));
  %DEBUG:
  disp(length(ix));

  stn.ngdc_depth = nanmedian(stn.ngdc_92m_bathy.field(ix));
  stn.ngdc_depth = nanmean(stn.ngdc_92m_bathy.field(ix));

return;
