function stn = do30(stn)

error('Run-once testing script!');

  [stn,Q0_factor] = ...
      station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid',...
                        'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
                        'ndbc_bulk_30',...
                        'ncep_dsrf','ndbc_bulk_rh_dlrf','ncep_precip');

  stn.ndbc_bulk_30_heat_flux_sum.date = stn.ndbc_bulk_30_heat_flux_term.date;
  stn.ndbc_bulk_30_heat_flux_sum.data = cumsum(stn.ndbc_bulk_30_heat_flux_term.data);

return;
