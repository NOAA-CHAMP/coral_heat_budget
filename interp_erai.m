function [stn,fld] = interp_erai(stn,fld,method)
%function [stn,fld] = interp_erai(stn,fld,method)
%
% Interpolate ERAI gridsquare field in STRUCT FLD (if missing or empty, FLD
% loaded by EXTRACT_ERAI_GRIDSQUARE) to location STN.lon,STN.lat. Optional
% arg METHOD (DEFAULT: 'linear') is passed to INTERP3 (v.)
%
% Last Saved Time-stamp: <Fri 2014-10-17 11:21:17 Eastern Daylight Time gramer>

  if ( ~exist('stn','var') || any(~isfield(stn,{'lon','lat'})) )
    error('First arg must be a STRUCT with .lon and .lat fields');
  end;
  if ( ~exist('fld','var') || isempty(fld) )
    extract_erai_gridsquare;
  end;
  if ( ~exist('method','var') || isempty(method) )
    method = 'linear';
  end;

  stn.erai_wind_u.date = fld.dts;
  stn.erai_wind_u.data = interp3(fld.lats,fld.dts,fld.lons,fld.u,stn.lat,fld.dts,stn.lon,method);

  stn.erai_wind_v.date = fld.dts;
  stn.erai_wind_v.data = interp3(fld.lats,fld.dts,fld.lons,fld.v,stn.lat,fld.dts,stn.lon,method);

  stn.erai_air_t.date = fld.dts;
  stn.erai_air_t.data = interp3(fld.lats,fld.dts,fld.lons,fld.Ta,stn.lat,fld.dts,stn.lon,method);

  stn.erai_dew_t.date = fld.dts;
  stn.erai_dew_t.data = interp3(fld.lats,fld.dts,fld.lons,fld.Td,stn.lat,fld.dts,stn.lon,method);

  stn.erai_albedo.date = fld.dts;
  stn.erai_albedo.data = interp3(fld.lats,fld.dts,fld.lons,fld.A,stn.lat,fld.dts,stn.lon,method);

  stn.erai_sea_t.date = fld.dts;
  stn.erai_sea_t.data = interp3(fld.lats,fld.dts,fld.lons,fld.Ts,stn.lat,fld.dts,stn.lon,method);

  stn.erai_cloud_cover.date = fld.dts;
  stn.erai_cloud_cover.data = interp3(fld.lats,fld.dts,fld.lons,fld.C,stn.lat,fld.dts,stn.lon,method);

return;
