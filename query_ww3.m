function stn = query_ww3(stn,begdt,enddt)
%function stn = query_ww3(stn)

  if ( ~isfield(stn,'lon') )
    [stn.lon,stn.lat,stn.depth]=get_station_coords(stn.station_name);
  end;

  % x = read_grib('../data/ww3/multi_1.at_4m.hs.200502.grb2',-1, 'ScreenDiag',0);
  x = read_grib('../data/ww3/wna.hs.200711.grb',-1, 'ScreenDiag',0);

  stn.ww3_latix = round( (x(1).gds.La1 - stn.lat) / x(1).gds.Di ) + 2;
  lon=stn.lon; lon(lon<0)=lon(lon<0)+360;
  stn.ww3_lonix = round( (lon - x(1).gds.Lo1) / x(1).gds.Dj ) + 1;

  for ix = 1:length(x)
    Hs = reshape(x(ix).fltarray, [x(1).gds.Ni x(1).gds.Nj])';
    Hs(Hs > 100) = nan;

    stn.ww3_sigwvhgt.date = ???
    stn.ww3_sigwvhgt.data(end+1:end+ndat) = Hs(stn.wwiii_latix+1,stn.wwiii_lonix);

    % DEBUG:    if (ix==1); dix=3; figure; contourf(Hs(stn.ww3_latix-dix:stn.ww3_latix+dix,stn.ww3_lonix-dix:stn.ww3_lonix+dix)); set(gca,'ydir','reverse'); maxigraph; colorbar; grid on; shading faceted; end;

  end;

return;
