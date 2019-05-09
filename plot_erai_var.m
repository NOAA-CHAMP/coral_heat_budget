function [dat,dts,lon,lat] = plot_erai_var(nc,var)
%function [dat,dts,lon,lat] = plot_erai_var(nc,var)

  atrs = getAttributes(nc);

  dtstr = strrep(strrep(atrs.CoordinateModelRunDate,'T',' '),'Z','');
  basedt = datenum(dtstr,'yyyy-mm-dd HH:MM:SS');

  hrs = cast( nc{'time'}(:,:,:,:), 'double' );
  dts = basedt + (hrs./24.0);
  lon = cast( nc{'lon'}(:,:,:,:), 'double' );
  lat = cast( nc{'lat'}(:,:,:,:), 'double' );
  midx = round(length(lon)/2);
  midy = round(length(lat)/2);
  dat = squeeze(cast( nc{var}(:,midy,midx), 'double' ));

  figure;
  maxigraph;
  plot(dts,dat);
  datetick3;
  titlename([strrep(var,'_','\_') ' from ' dtstr]);

  if ( nargout < 3 ); clear lat lon; end;
  if ( nargout < 2 ); clear dts; end;
  if ( nargout < 1 ); clear dat; end;

return;
