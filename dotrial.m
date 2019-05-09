function [R,coh] = dotrial(stn)

  flds = { ...
      'ncep_net_heat_flux', ...
      'ncep_bulk_net_heat_flux', ...
      'ndbc_hfbulk_net_heat_flux', ...
      'ndbc_bulk_rh_net_heat_flux', ...
      'ndbc_bulk_26_net_heat_flux', ...
      'ndbc_bulk_30_net_heat_flux', ...
      'ndbc_bulk_30a_net_heat_flux', ...
      'ndbc_ncep_26_net_heat_flux', ...
      'ndbc_ncep_30_net_heat_flux', ...
      'ndbc_ncep_30a_net_heat_flux', ...
         };

  q0 = Q0factor(stn.ndbc_sea_t.data,[],2);

  for fldix = 1:length(flds)
    fld = flds{fldix};
    [ix1,ix2] = intersect_dates(stn.ndbc_sea_t.date,stn.(fld).date);

    seat.date = stn.ndbc_sea_t.date(ix1);
    seat.data = stn.ndbc_sea_t.data(ix1);
    flux.date = stn.(fld).date(ix2);
    flux.data = stn.(fld).data(ix2);

    badix = find(~isfinite(seat.data) | ~isfinite(flux.data));
    seat.date(badix) = [];
    seat.data(badix) = [];
    flux.date(badix) = [];
    flux.data(badix) = [];

    term = real( (flux.data ./ q0) * (3600) );

    R(fldix) = corr2(seat.data,cumsum(term));
    cxy = mscohere(seat.data,cumsum(term));
    coh(fldix) = nanmean(cxy);
  end;

  figure;
  plot(coh);
  maxigraph;
  title(sprintf('%s - MS coher. T_s_e_a vs. \\Sigma Q_0', stn.station_name));

return;
