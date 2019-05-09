function stn = compsrf(stn)
%function stn = compsrf(stn)

  [yr,ig,ig] = datevec(stn.ncep_dsrf.date);
  jd = stn.ncep_dsrf.date - datenum(yr,1,1);
  [qsw,a] = swhf(jd,yr,-stn.lon,stn.lat,stn.ncep_dsrf.data);

  stn.albedo.date = stn.ncep_dsrf.date;
  stn.albedo.data = a;

  q0f = Q0factor(stn.ndbc_sea_t.data,[],2);
  Le = vapor(nanmean(stn.ndbc_sea_t.data));

  stn.ndbc_air_sea_spechumid_term.date = stn.ndbc_air_sea_spechumid.date;
  stn.ndbc_air_sea_spechumid_term.data = stn.ndbc_air_sea_spechumid.data.*Le.*10.*1e-3;

  stn.ndbc_bulk_rh_lrf_term.date = stn.ndbc_bulk_rh_lrf.date;
  stn.ndbc_bulk_rh_lrf_term.data = (stn.ndbc_bulk_rh_lrf.data ./ q0f) .* (60*60);

  stn.qsw.date = stn.ncep_dsrf.date;
  stn.qsw.data = qsw;
  stn.qswterm.date = stn.qsw.date;
  stn.qswterm.data = (qsw ./ q0f) .* (60*60);

  [ix1,ix2,ix3,ix4] = ...
      intersect_all_dates([], ...
                          stn.ndbc_bulk_30_latent_heat_flux.date, ...
                          stn.ndbc_bulk_30_sensible_heat_flux.date, ...
                          stn.ndbc_bulk_rh_lrf.date, ...
                          stn.ncep_dsrf.date);

  stn.q0.date = stn.ndbc_bulk_30_latent_heat_flux.date(ix1);
  stn.q0.data = stn.ndbc_bulk_30_latent_heat_flux.data(ix1) + ...
      stn.ndbc_bulk_30_sensible_heat_flux.data(ix2) + ...
      stn.ndbc_bulk_rh_lrf.data(ix3) + ...
      qsw(ix4);

  stn.q0term.date = stn.q0.date;
  stn.q0term.data = (stn.q0.data ./ q0f) .* (60*60);

  figure;
  plot(stn.nocs_heat_flux_term.date, stn.nocs_heat_flux_term.data, 'b-', ...
       stn.ncep_heat_flux_term.date, stn.ncep_heat_flux_term.data, 'g:', ...
       stn.q0term.date, stn.q0term.data, 'c--');
  maxigraph;
  datetick3;
  colormap(gray);
  legend('NOCS','NCEP','AIR\_SEA SWHF');
  titlename(sprintf('%s - Comparing net shortwave radiative flux estimates'));

return;
