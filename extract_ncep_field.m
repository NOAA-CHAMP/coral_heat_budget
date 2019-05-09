function stn = extract_ncep_field(dat, stn, stfld, ncepfld)
%function stn = extract_ncep_field(dat, stn, stfld, ncepfld)

  if ( ~isfield(dat,'date') )
    warning('No struct field "date" in first argument!');
  elseif ( ~isfield(dat,ncepfld) )
    warning('No struct field "%s" in first argument!', ncepfld);
  else
    stn.(['ncep_' stfld]).date = dat.date(:);
    % stn.(['ncep_' stfld]).data = squeeze(dat.(ncepfld)(:,stn.lonix,stn.latix));
    stn.(['ncep_' stfld]).data = squeeze(dat.(ncepfld)(:,stn.latix,stn.lonix));
  end;

return;
