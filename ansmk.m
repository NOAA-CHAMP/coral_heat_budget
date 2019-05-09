1;


set_more off

stn = calcdt('smkf1');


% Modify to plot individual terms for ndbc_ncep_30a
dy=110; plot_fluxes(stn,2005,dy,dy+40);

multiplot_station(stn,{'ndbc_wind1_dir','ndbc_wind1_u','quasi_eulerian_u','ndbc_wind1_v','quasi_eulerian_v','advected_heat'},[],[],[],[datenum(2005,4,10) datenum(2005,4,21)]);

multiplot_station(stn,{'ndbc_ncep_30a_heat_flux_term','ndbc_ncep_30a_latent_flux_term','ndbc_ncep_30a_sensible_flux_term','ndbc_ncep_30a_shortwave_flux_term','ndbc_ncep_30a_longwave_flux_term','advected_heat'},[],[],[],[datenum(2005,4,10) datenum(2005,4,21)]);

multiplot_station(stn,{'ndbc_ncep_30a_heat_flux_term','ndbc_ncep_30a_latent_flux_term','ndbc_ncep_30a_sensible_flux_term','ndbc_ncep_30a_shortwave_flux_term','ndbc_ncep_30a_longwave_flux_term','advected_heat'},[],[],[],[datenum(2005,4,10) datenum(2005,4,21)],[-.5 .5]);

station_plot_fields(stn,'avhrr_weekly_sst',datenum(2005,4,10):(1/24):datenum(2005,4,17));

set_more
