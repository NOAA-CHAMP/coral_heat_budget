function stn = plot_xg(stn,fnm,interpMethod)
%function stn = plot_xg(stn,fnm,interpMethod)
% Plot cross-shore gradient in field STN.(FNM) (DEFAULT fkeys_hycom_seatemp_field).

  if ( ~exist('fnm','var') || isempty(fnm) )
    fnm = 'fkeys_hycom_seatemp_field';
  end;

  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;

  if ( ~isfield(stn.(fnm),'gradient_x') )
    stn = calc_field_terms(stn,fnm);
  end;

  gx = interp_field(stn.(fnm).lat,stn.(fnm).lon,stn.(fnm).gradient_x,stn.lat,stn.lon,interpMethod);
  gy = interp_field(stn.(fnm).lat,stn.(fnm).lon,stn.(fnm).gradient_y,stn.lat,stn.lon,interpMethod);

  ori=52.6;
  [xg,lg] = reorient_vectors(ori,gx,gy);

  % xg = xg*100;
  xg = 12*3600*0.01*xg;

  % fmg; plot(stn.(fnm).date,xg); datetick3;
  % scatter_fit(get_yearday(stn.(fnm).date),xg);
  fmg; plot(1:366,grpstats(xg,get_jday(stn.(fnm).date))); datetick3;

  % titlename([stn.station_name ' ' strrep(fnm,'_','\_') ' over 100m']);
  titlename([stn.station_name ' ' strrep(fnm,'_','\_') ' * 0.01 m/s * 3600s * 12h']);

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
