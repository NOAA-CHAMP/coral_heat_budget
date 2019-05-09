function [stn,res] = get_tsg_station(stn_or_stnm,res,rad)
%function [stn,res] = get_tsg_station(stn_or_stnm,res,rad)
%
% Load South Florida Program (SFP) thermosalinograph (TSG) data from Excel
% files and subset for a given station (or station name string) STN_OR_STNM.

  stn = get_station_from_station_name(stn_or_stnm);

  if ( ~exist('res','var') || ~isfield(res,'sst') )
    res = antsg();
  end;
  if ( ~exist('rad','var') || isempty(rad) )
    rad = 15; %[km]
  end;

  if ( ~isfield(stn,'isobath_orientation') )
    stn = station_optimal_isobath_orientation(stn);
  end;


  % Efficiently calculate distance [km] between station and each point in RES
  lon = repmat(stn.lon,[2*numel(res.lon),1]);
  lat = repmat(stn.lat,[2*numel(res.lat),1]);
  lon(1:2:end) = res.lon;
  lat(1:2:end) = res.lat;
  site_d2 = sw_dist(lat,lon,'km');
  site_d = site_d2(1:2:end);
  lon=[]; lat=[]; site_d2=[]; clear lon lat site_d2;

  relhdg = res.hdg - stn.isobath_orientation - 90;
  relhdg(relhdg>90) = 180 - relhdg(relhdg>90);

  % goodix = find( (abs(relhdg)<20) & (site_d<rad) & isfinite(res.crestdsst) );
  % goodix = find( (site_d<rad) & (res.z<-8) );
  % goodix = find( (site_d<rad) & (res.z<-4) & (0.01<=res.dz&res.dz<=0.10) );
  % goodix = find( (site_d<rad) );

  % Pick data points near site, on "nearly cross-shore" headings, near crest,
  % and with reasonable SST changes, lat/lon distances (km/min) and time gaps
  goodix = ...
      find( (site_d<rad) & (abs(relhdg)<20) & (-1>=res.z & res.z>=-100) & ...
            (abs(res.dsst)<5) & (0.05<=res.d & res.d<=0.50) & ([0;diff(res.date)]<15) );


  npts = numel(res.dsst);
  flds = fieldnames(res);
  for fldix=1:length(flds)
    fld = flds{fldix};
    if ( numel(res.(fld)) == npts )
      stn.tsg.(fld) = res.(fld)(goodix);
    else
      stn.tsg.(fld) = res.(fld);
    end;
  end;
  stn.tsg.relhdg = relhdg(goodix);

  % Flip sign of dT/dx for inshore (inbound) transects
  stn.tsg.dsst = sign(sind(stn.tsg.hdg-stn.isobath_orientation)) .* stn.tsg.dsst;

  % Smooth gradients with a 2*N+1 (~AVGPTS/2 km) point average
  avgpts = 1;
  dts=[];
  lon=[];
  lat=[];
  d=[];
  relhdg=[];
  dsst=[];
  for ix = -avgpts:avgpts
    dts    = [dts ,    stn.tsg.date(avgpts+ix+1:end-avgpts+ix)];
    lon    = [lon ,    stn.tsg.lon(avgpts+ix+1:end-avgpts+ix)];
    lat    = [lat ,    stn.tsg.lat(avgpts+ix+1:end-avgpts+ix)];
    d      = [d ,      stn.tsg.d(avgpts+ix+1:end-avgpts+ix)];
    relhdg = [relhdg , stn.tsg.relhdg(avgpts+ix+1:end-avgpts+ix)];
    dsst   = [dsst ,   stn.tsg.dsst(avgpts+ix+1:end-avgpts+ix)];
  end;

  gapix = find(diff(stn.tsg.date)>15);
  dts(gapix)    = nan;
  lon(gapix)    = nan;
  lat(gapix)    = nan;
  d(gapix)      = nan;
  relhdg(gapix) = nan;
  dsst(gapix)   = nan;

  stn.smooth_tsg.date   = mean(dts,2);
  stn.smooth_tsg.lon    = mean(lon,2);
  stn.smooth_tsg.lat    = mean(lat,2);
  stn.smooth_tsg.d      = sum(d,2);
  stn.smooth_tsg.relhdg = mean(relhdg,2);
  stn.smooth_tsg.dsst   = mean(dsst,2);

  % Convert K/km -> K/m
  stn.smooth_tsg_ts.date = stn.smooth_tsg.date;
  stn.smooth_tsg_ts.data = stn.smooth_tsg.dsst ./ 1e3;

  if ( nargout < 2 )
    res = []; clear res;
  end;

return;
