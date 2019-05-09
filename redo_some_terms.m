function stn = redo_some_terms(stn)

error('Run-once fix-it code!');

  set_more off;

  if ( ~isfield(stn,'ncep_swhf') )
    [yr,ig,ig] = datevec(stn.ncep_dsrf.date);
    yd = stn.ncep_dsrf.date - datenum(yr,1,1);

    stn.ncep_swhf.date = stn.ncep_dsrf.date;
    stn.ncep_swhf_albedo.date = stn.ncep_dsrf.date;
    stn.ncep_swhf_sorad.date = stn.ncep_dsrf.date;
    stn.ncep_swhf_sunalt.date = stn.ncep_dsrf.date;
    stn.ncep_swhf_trans.date = stn.ncep_dsrf.date;

    %[stn.ncep_swhf.data,stn.ncep_swhf_albedo.data] = swhf(yd,yr,-stn.lon,stn.lat,stn.ncep_dsrf.data);

    [stn.ncep_swhf_sunalt.data, stn.ncep_swhf_sorad.data] = soradna1(yd,yr,-stn.lon,stn.lat);

    stn.ncep_swhf_trans.data = repmat(inf, size(stn.ncep_swhf_sorad.data));
    goodix = find(stn.ncep_swhf_sorad.data > 0);
    stn.ncep_swhf_trans.data(goodix) = stn.ncep_dsrf.data(goodix) ./ stn.ncep_swhf_sorad.data(goodix);

    % HACK!
    badix = find(isfinite(stn.ncep_swhf_trans.data) & stn.ncep_swhf_trans.data > 1);
    stn.ncep_swhf_trans.data(badix) = 1.0;

    stn.ncep_swhf_albedo.data = albedo(stn.ncep_swhf_trans.data,stn.ncep_swhf_sunalt.data);

    stn.ncep_swhf.data = (1 - stn.ncep_swhf_albedo.data) .* stn.ncep_dsrf.data;
    badix = find(isnan(stn.ncep_swhf.data));
    stn.ncep_swhf.data(badix) = 0;
  end;

  if ( ~isfield(stn,'ncep_lwhf') )
    [ix1,ix2,ix3] = intersect_all_dates([], stn.ndbc_sea_t.date, stn.ncep_dlrf.date, ...
                                        stn.ncep_dsrf.date);

    stn.ncep_lwhf.date = stn.ndbc_sea_t.date(ix1);
    stn.ncep_lwhf.data = lwhf(stn.ndbc_sea_t.data(ix1), stn.ncep_dlrf.data(ix2), ...
                              stn.ncep_dsrf.data(ix3));
  end;

  rhfld = 'ndbc_relhumid';
  [stn,Q0_factor] = ...
      station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                        'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lwhf',...
                        'ndbc_ncep_26','ncep_dsrf','ncep_dlrf');
  [stn,Q0_factor] = ...
      station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                        'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lwhf',...
                        'ndbc_ncep_30','ncep_dsrf','ncep_dlrf','ncep_precip');
  [stn,Q0_factor] = ...
      station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                        'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lwhf',...
                        'ndbc_ncep_30a','ncep_dsrf','ncep_dlrf','ncep_precip',...
                        'ndbc_wind1_dir','default','default','default','default','default',true);


  stn.ncep_swhf_term.date = stn.ncep_swhf.date;
  stn.ncep_swhf_term.data = (stn.ncep_swhf.data ./ Q0_factor) .* 3600;
  stn.ncep_srf_term.date = stn.ncep_srf.date;
  stn.ncep_srf_term.data = (stn.ncep_srf.data ./ Q0_factor) .* 3600;
  stn.ncep_lwhf_term.date = stn.ncep_lwhf.date;
  stn.ncep_lwhf_term.data = (stn.ncep_lwhf.data ./ Q0_factor) .* 3600;
  stn.ncep_lrf_term.date = stn.ncep_lrf.date;
  stn.ncep_lrf_term.data = (stn.ncep_lrf.data ./ Q0_factor) .* 3600;

  set_more;

return;
