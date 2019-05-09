function [lons,lats,ds,tprof,vprof,topo] = view_hycom_profile(stn,dt,doFKEYS)
%function [lons,lats,ds,tprof,vprof,topo] = view_hycom_profile(stn,dt,doFKEYS)
%
% Extract (via netCDF Java Toolbox call MDATASET - v.) and contour-plot the
% cross-shore transects of sea temperature and along-shore velocity from GoM
% HYCOM 1/25o output, at station STN (struct) or station named STN (char),
% for one date DT (must be a double in the OPeNDAP form YYYYMMDD - *not* a
% MATLAB DATENUM). If DT is during 2008 and DOFKEYS is true, also do FlKeys.
%
% Last Saved Time-stamp: <Fri 2010-10-29 14:04:36 Eastern Daylight Time gramer>

  %DEBUG:
  tic,

  if ( ischar(stn) )
    stnm = lower(stn); clear stn;
    [stn.lon,stn.lat,stn.depth]=get_station_coords(stnm);
    stn.station_name = upper(stnm);
  end;
  if ( ~exist('dt','var') || isempty(dt) )
    dt = 20080522;
  end;
  if ( ~exist('doFKEYS','var') || isempty(doFKEYS) )
    doFKEYS = false;
  end;

  topo = [];
  if ( isfield(stn,'ngdc_92m_bathy') )
    %DEBUG:
    disp('Extracting bathymetry');
    % bth = stn.ngdc_92m_bathy;
    % topo.lon = lons;
    % topo.field = interp2(bth.lon,bth.lat,bth.field',lons,lats,'nearest');
    topo = station_field_transect(stn,'ngdc_92m_bathy',[-16:(0.092):16],90);
  end;


  % minlon = stn.lon - 0.097;
  % maxlon = stn.lon + 0.097;
  minlon = -80.16;
  maxlon = -80.00;

  if ( doFKEYS && 20080101 <= dt && dt < 20090101 )
    [yix,xix] = query_flkeys_hycom_indices(stn.lon,stn.lat);
    xrad = 10;
    nc = mDataset('http://tds.hycom.org/thredds/dodsC/flkeys');
    pause(1);
    ds = cast(nc{'Depth'}(1:end),'double');
    lons = cast(nc{'Longitude'}(xix-xrad:xix+xrad),'double');
    dts = cast(nc{'Date'}(:),'double');
    [ig,dtix] = min(abs(dts - dt));
    pause(1);
    tprof = cast(nc{'temperature'}(dtix,1:end,yix,xix-xrad:1:xix+xrad),'double');
    vprof = cast(nc{'v'}(dtix,1:end,yix,xix-xrad:1:xix+xrad),'double');
    close(nc); clear nc;
    mint = floor(min(tprof(:)));
    figure; maxigraph; hold on;
    contourf(lons,-ds,tprof,[mint:0.5:27]);
    if ( ~isempty(topo) ); plot(topo.lon,topo.field,'k-'); end;
    xlim([minlon maxlon]); ylim([-100 0]); caxis([12 27]); colorbar;
    titlename(['FlKeys HYCOM T ' num2str(dt)]);

    figure; maxigraph; hold on;
    contourf(lons,-ds,vprof);
    if ( ~isempty(topo) ); plot(topo.lon,topo.field,'k-'); end;
    xlim([minlon maxlon]); ylim([-100 0]); caxis([-1 +2.5]); colorbar;
    titlename(['FlKeys HYCOM V ' num2str(dt)]);

    clear xrad ds lons dts dtix tprof vprof;
  end;

  dixes = 1:20;

  [yix,xix] = query_gom_hycom_indices(stn.lon,stn.lat);
  xrad = 3; yrad = 0;
  yr = floor(dt/10000);
  %DEBUG:
  disp('Open dataset...');
  nc = mDataset(['http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1/' num2str(yr)]);
  ds = cast(nc{'Depth'}(dixes),'double');
  lats = cast(nc{'Latitude'}(yix-yrad:yix+yrad),'double');
  lons = cast(nc{'Longitude'}(xix-xrad:xix+xrad),'double');
  dts = cast(nc{'Date'}(:),'double');
  [ig,dtix] = min(abs(dts - dt));
  %DEBUG:
  disp('tprof');
  tprof = cast(nc{'temperature'}(dtix,dixes,yix,xix-xrad:xix+xrad),'double');
  %DEBUG:
  disp('vprof');
  vprof = cast(nc{'v'}(dtix,dixes,yix,xix-xrad:xix+xrad),'double');
  close(nc); clear nc;

  %DEBUG:
  disp('plot');
  mint = floor(min(tprof(:)));
  maxt = ceil(max(tprof(:)));

  figure; maxigraph; hold on;
  contourf(lons,-ds,tprof,[mint:0.5:maxt]);
  if ( ~isempty(topo) ); plot(topo.lon,topo.field,'k-'); end;
  xlim([minlon maxlon]); ylim([-100 0]); caxis([12 27]); colorbar;
  titlename(['Gulf of Mexico HYCOM T ' num2str(dt)]);

  figure; maxigraph; hold on;
  contourf(lons,-ds,vprof);
  if ( ~isempty(topo) ); plot(topo.lon,topo.field,'k-'); end;
  xlim([minlon maxlon]); ylim([-100 0]); caxis([-1 +2.5]); colorbar;
  titlename(['Gulf of Mexico HYCOM V ' num2str(dt)]);

  %DEBUG:
  toc,

return;
