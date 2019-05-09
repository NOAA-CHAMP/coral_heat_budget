1;

error('THIS SCRIPT IS NO LONGER USEFUL - SHOULD HAVE RUN ONCE ONLY!');

  mlrf1.ncep_bulk_heat_flux_sum.date = mlrf1.ncep_bulk_heat_flux_term.date;
  mlrf1.ncep_bulk_heat_flux_sum.data = cumsum(mlrf1.ncep_bulk_heat_flux_term.data);


  mlrf1.ncep_heat_flux_sum.date = mlrf1.ncep_heat_flux_term.date;
  mlrf1.ncep_heat_flux_sum.data = cumsum(mlrf1.ncep_heat_flux_term.data);


    mlrf1.ndbc_bulk_heat_flux_sum.date = mlrf1.ndbc_bulk_heat_flux_term.date;
    mlrf1.ndbc_bulk_heat_flux_sum.data = cumsum(mlrf1.ndbc_bulk_heat_flux_term.data);



  mlrf1.sat_heat_flux_sum.date = mlrf1.sat_heat_flux_term.date;
  mlrf1.sat_heat_flux_sum.data = cumsum(mlrf1.sat_heat_flux_term.data);


  % Diagnostic fields
  mlrf1.ncep_dsrf_sum.date = mlrf1.ncep_dsrf.date;
  mlrf1.ncep_dsrf_sum.data = cumsum(mlrf1.ncep_dsrf.date);
  mlrf1.licor_surf_dsrf_sum.date = mlrf1.licor_surf_dsrf.date;
  mlrf1.licor_surf_dsrf_sum.data = cumsum(mlrf1.licor_surf_dsrf.data);

  mlrf1.ncep_usrf_sum.date = mlrf1.ncep_usrf.date;
  mlrf1.ncep_usrf_sum.data = cumsum(mlrf1.ncep_usrf.data);
  mlrf1.licor_surf_usrf_sum.date = mlrf1.licor_surf_usrf.date;
  mlrf1.licor_surf_usrf_sum.data = cumsum(mlrf1.licor_surf_usrf.data);

  mlrf1.ncep_dlrf_sum.date = mlrf1.ncep_dlrf.date;
  mlrf1.ncep_dlrf_sum.data = cumsum(mlrf1.ncep_dlrf.data);
  mlrf1.ndbc_bulk_dlrf_sum.date = mlrf1.ndbc_bulk_dlrf.date;
  mlrf1.ndbc_bulk_dlrf_sum.data = cumsum(mlrf1.ndbc_bulk_dlrf.data);

  mlrf1.ncep_ulrf_sum.date = mlrf1.ncep_ulrf.date;
  mlrf1.ncep_ulrf_sum.data = cumsum(mlrf1.ncep_ulrf.data);
  mlrf1.ndbc_bulk_ulrf_sum.date = mlrf1.ndbc_bulk_ulrf.date;
  mlrf1.ndbc_bulk_ulrf_sum.data = cumsum(mlrf1.ndbc_bulk_ulrf.data);
