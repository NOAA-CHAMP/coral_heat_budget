function stn = tweak_term(stn)
%function stn = tweak_term(stn)
%
% This was a one-time (per station) fix-it function!

error('This was a one-time (per station) fix-it function!');

  if ( isfield(stn,'ncep_heat_flux_term') )
%     stn.ncep_heat_flux_term.data = stn.ncep_heat_flux_term.data .* 3600;
    stn.ncep_heat_flux_sum.date = stn.ncep_heat_flux_term.date;
    stn.ncep_heat_flux_sum.data = cumsum(stn.ncep_heat_flux_term.data);
    if ( isfield(stn,'ncep_heat_flux_term_1_day_sum') )
      stn = rmfield(stn,'ncep_heat_flux_term_1_day_sum');
    end;
    stn = verify_variable(stn,'ncep_heat_flux_term_1_day_sum');
  end;
  if ( isfield(stn,'ncep_bulk_heat_flux_term') )
%     stn.ncep_bulk_heat_flux_term.data = stn.ncep_bulk_heat_flux_term.data .* 3600;
    stn.ncep_bulk_heat_flux_sum.date = stn.ncep_bulk_heat_flux_term.date;
    stn.ncep_bulk_heat_flux_sum.data = cumsum(stn.ncep_bulk_heat_flux_term.data);
    if ( isfield(stn,'ncep_bulk_heat_flux_term_1_day_sum') )
      stn = rmfield(stn,'ncep_bulk_heat_flux_term_1_day_sum');
    end;
    stn = verify_variable(stn,'ncep_bulk_heat_flux_term_1_day_sum');
  end;
  if ( isfield(stn,'ndbc_bulk_heat_flux_term') )
%     stn.ndbc_bulk_heat_flux_term.data = stn.ndbc_bulk_heat_flux_term.data .* 3600;
    stn.ndbc_bulk_heat_flux_sum.date = stn.ndbc_bulk_heat_flux_term.date;
    stn.ndbc_bulk_heat_flux_sum.data = cumsum(stn.ndbc_bulk_heat_flux_term.data);
    if ( isfield(stn,'ndbc_bulk_heat_flux_term_1_day_sum') )
      stn = rmfield(stn,'ndbc_bulk_heat_flux_term_1_day_sum');
    end;
    stn = verify_variable(stn,'ndbc_bulk_heat_flux_term_1_day_sum');
  end;
  if ( isfield(stn,'sat_heat_flux_term') )
%     stn.sat_heat_flux_term.data = stn.sat_heat_flux_term.data .* 3600;
    stn.sat_heat_flux_sum.date = stn.sat_heat_flux_term.date;
    stn.sat_heat_flux_sum.data = cumsum(stn.sat_heat_flux_term.data);
    if ( isfield(stn,'sat_heat_flux_term_1_day_sum') )
      stn = rmfield(stn,'sat_heat_flux_term_1_day_sum');
    end;
    stn = verify_variable(stn,'sat_heat_flux_term_1_day_sum');
  end;

return;
