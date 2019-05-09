function view_hycom_tprof(stn,dt)

  if ( ischar(stn) )
    stnm = lower(stn); clear stn;
    [stn.lon,stn.lat,stn.depth]=get_station_coords(stnm);
    stn.station_name = upper(stnm);
  end;
  if ( ~exist('dt','var') || isempty(dt) )
    dt = 20080522;
  end;

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
  figure; maxigraph; contourf(lons,-ds,tprof,[mint:0.5:27]);
  xlim([-80.18 -80.00]); ylim([-100 0]); caxis([24 27]); colorbar;
  title(['FlKeys T ' num2str(dt)]);

  figure; maxigraph; contourf(lons,-ds,vprof);
  xlim([-80.18 -80.00]); ylim([-100 0]); caxis([-1 +2.5]); colorbar;
  title(['FlKeys V ' num2str(dt)]);

  clear xrad ds lons dts dtix tprof;

  [yix,xix] = query_gom_hycom_indices(stn.lon,stn.lat);
  xrad = 2;
  nc = mDataset('http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1/2008');
  ds = cast(nc{'Depth'}(1:10),'double');
  lons = cast(nc{'Longitude'}(xix-xrad:xix+xrad),'double');
  dts = cast(nc{'Date'}(:),'double');
  [ig,dtix] = min(abs(dts - dt));
  tprof = cast(nc{'temperature'}(dtix,1:10,yix,xix-xrad:xix+xrad),'double');
  close(nc); clear nc;
  mint = floor(min(tprof(:)));
  figure; maxigraph; contourf(lons,-ds,tprof,[mint:0.5:27]);
  xlim([-80.18 -80.00]); ylim([-100 0]); caxis([24 27]); colorbar;
  title(['GoM T ' num2str(dt)]);

return;
