1;

% ori = 52.25;
% ori = 52.4;
ori = 52.6;

stn = station_reorient_vectors(stn,ori,'fkeys_hycom_u','fkeys_hycom_v','fkeys_hycom_x','fkeys_hycom_l');
multiplot_station(stn,{'fkeys_hycom_u','fkeys_hycom_v','fkeys_hycom_x','fkeys_hycom_l'},[],[],[],[],[-1,1.5]);
titlename([stn.station_name ' ' num2str(ori) '^o']);

stn.fkeys_hycom_seatemp_x.date = stn.fkeys_hycom_seatemp_field.date;
stn.fkeys_hycom_seatemp_x.data = interp_field(stn.fkeys_hycom_seatemp_field.lat,stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.gradient_x,stn.lat,stn.lon);

stn.fkeys_hycom_seatemp_y.date = stn.fkeys_hycom_seatemp_field.date;
stn.fkeys_hycom_seatemp_y.data = interp_field(stn.fkeys_hycom_seatemp_field.lat,stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.gradient_y,stn.lat,stn.lon);

stn = station_reorient_vectors(stn,ori,'fkeys_hycom_seatemp_x','fkeys_hycom_seatemp_y','fkeys_hycom_seatemp_cross','fkeys_hycom_seatemp_along');
multiplot_station(stn,{'fkeys_hycom_seatemp_x','fkeys_hycom_seatemp_y','fkeys_hycom_seatemp_cross','fkeys_hycom_seatemp_along'},[],[],[],[],[-1.5e-3,2e-3]);
titlename([stn.station_name ' ' num2str(ori) '^o']);

stn.x.date=stn.fkeys_hycom_seatemp_cross.date;
stn.x.data=-6*3600.*stn.fkeys_hycom_seatemp_cross.data.*stn.fkeys_hycom_x.data;

stn.l.date=stn.fkeys_hycom_seatemp_along.date;
stn.l.data=-6*3600.*stn.fkeys_hycom_seatemp_along.data.*stn.fkeys_hycom_l.data;

stn.a.date=stn.x.date;
stn.a.data=stn.x.data+stn.l.data;

% fmg; plot_ts(stn.x,stn.l,stn.a);

% fmg; plot(stn.x.date,[cumsum(stn.x.data),cumsum(stn.l.data),cumsum(stn.a.data)]); datetick3;
% legend('\partial_xu\times\partial_xT','\partial_lu\times\partial_lT','u^.\nablaT');
% titlename([stn.station_name ' ' num2str(ori) '^o']);

[ig,ix] = intersect_dates(stn.a.date,stn.benthic_ndbc_erai_30a_net_heat_flux_term.date);
fmg; plot(stn.x.date,[cumsum(stn.x.data),cumsum(stn.l.data),cumsum(stn.a.data)],stn.benthic_ndbc_erai_30a_net_heat_flux_term.date(ix(1):ix(end)),cumsum(stn.benthic_ndbc_erai_30a_net_heat_flux_term.data(ix(1):ix(end))),'y:'); datetick3;
legend('\partial_xu\times\partial_xT','\partial_lu\times\partial_lT','u^.\nablaT','(Q_0+Q_b)/\rhoC_ph');
titlename([stn.station_name ' ' num2str(ori) '^o']);


return;




error('This is a one-time fixit script! Not to be rerun...');

fclose('all'); close all; clear all; dbclear all;
load('../data/smkf1_thesis.mat');
[wz,az,pz,stz] = station_instrument_heights(smkf1.station_name);
Q0_factor = Q0factor(smkf1.ncep_sea_t.data,[],stz);
  smkf1.ncep_heat_flux_term.date = smkf1.ncep_net_heat_flux.date;
  smkf1.ncep_heat_flux_term.data = (smkf1.ncep_net_heat_flux.data ./ Q0_factor) .* (60*60);
  smkf1.ncep_srf_term.date = smkf1.ncep_srf.date;
  smkf1.ncep_srf_term.data = (smkf1.ncep_srf.data ./ Q0_factor) .* (60*60);
  smkf1.ncep_lrf_term.date = smkf1.ncep_lrf.date;
  smkf1.ncep_lrf_term.data = (smkf1.ncep_lrf.data ./ Q0_factor) .* (60*60);
  smkf1.ncep_latent_flux_term.date = smkf1.ncep_latent_heat_flux.date;
  smkf1.ncep_latent_flux_term.data = (smkf1.ncep_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  smkf1.ncep_sensible_flux_term.date = smkf1.ncep_sensible_heat_flux.date;
  smkf1.ncep_sensible_flux_term.data = (smkf1.ncep_sensible_heat_flux.data ./ Q0_factor) .* (60*60);

  smkf1.ncep_bulk_heat_flux_term.date = smkf1.ncep_bulk_net_heat_flux.date;
  smkf1.ncep_bulk_heat_flux_term.data = (smkf1.ncep_bulk_net_heat_flux.data ./ Q0_factor) .* (60*60);
  smkf1.ncep_bulk_latent_flux_term.date = smkf1.ncep_bulk_latent_heat_flux.date;
  smkf1.ncep_bulk_latent_flux_term.data = (smkf1.ncep_bulk_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  smkf1.ncep_bulk_sensible_flux_term.date = smkf1.ncep_bulk_sensible_heat_flux.date;
  smkf1.ncep_bulk_sensible_flux_term.data = (smkf1.ncep_bulk_sensible_heat_flux.data ./ Q0_factor) .* (60*60);
disp('Saving smkf1');
save('../data/smkf1_thesis.mat','smkf1');


fclose('all'); close all; clear all; dbclear all;
load('../data/mlrf1_thesis.mat');
[wz,az,pz,stz] = station_instrument_heights(mlrf1.station_name);
Q0_factor = Q0factor(mlrf1.ncep_sea_t.data,[],stz);
  mlrf1.ncep_heat_flux_term.date = mlrf1.ncep_net_heat_flux.date;
  mlrf1.ncep_heat_flux_term.data = (mlrf1.ncep_net_heat_flux.data ./ Q0_factor) .* (60*60);
  mlrf1.ncep_srf_term.date = mlrf1.ncep_srf.date;
  mlrf1.ncep_srf_term.data = (mlrf1.ncep_srf.data ./ Q0_factor) .* (60*60);
  mlrf1.ncep_lrf_term.date = mlrf1.ncep_lrf.date;
  mlrf1.ncep_lrf_term.data = (mlrf1.ncep_lrf.data ./ Q0_factor) .* (60*60);
  mlrf1.ncep_latent_flux_term.date = mlrf1.ncep_latent_heat_flux.date;
  mlrf1.ncep_latent_flux_term.data = (mlrf1.ncep_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  mlrf1.ncep_sensible_flux_term.date = mlrf1.ncep_sensible_heat_flux.date;
  mlrf1.ncep_sensible_flux_term.data = (mlrf1.ncep_sensible_heat_flux.data ./ Q0_factor) .* (60*60);

  mlrf1.ncep_bulk_heat_flux_term.date = mlrf1.ncep_bulk_net_heat_flux.date;
  mlrf1.ncep_bulk_heat_flux_term.data = (mlrf1.ncep_bulk_net_heat_flux.data ./ Q0_factor) .* (60*60);
  mlrf1.ncep_bulk_latent_flux_term.date = mlrf1.ncep_bulk_latent_heat_flux.date;
  mlrf1.ncep_bulk_latent_flux_term.data = (mlrf1.ncep_bulk_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  mlrf1.ncep_bulk_sensible_flux_term.date = mlrf1.ncep_bulk_sensible_heat_flux.date;
  mlrf1.ncep_bulk_sensible_flux_term.data = (mlrf1.ncep_bulk_sensible_heat_flux.data ./ Q0_factor) .* (60*60);
disp('Saving mlrf1');
save('../data/mlrf1_thesis.mat','mlrf1');


fclose('all'); close all; clear all; dbclear all;
load('../data/fwyf1_thesis.mat');
[wz,az,pz,stz] = station_instrument_heights(fwyf1.station_name);
Q0_factor = Q0factor(fwyf1.ncep_sea_t.data,[],stz);
  fwyf1.ncep_heat_flux_term.date = fwyf1.ncep_net_heat_flux.date;
  fwyf1.ncep_heat_flux_term.data = (fwyf1.ncep_net_heat_flux.data ./ Q0_factor) .* (60*60);
  fwyf1.ncep_srf_term.date = fwyf1.ncep_srf.date;
  fwyf1.ncep_srf_term.data = (fwyf1.ncep_srf.data ./ Q0_factor) .* (60*60);
  fwyf1.ncep_lrf_term.date = fwyf1.ncep_lrf.date;
  fwyf1.ncep_lrf_term.data = (fwyf1.ncep_lrf.data ./ Q0_factor) .* (60*60);
  fwyf1.ncep_latent_flux_term.date = fwyf1.ncep_latent_heat_flux.date;
  fwyf1.ncep_latent_flux_term.data = (fwyf1.ncep_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  fwyf1.ncep_sensible_flux_term.date = fwyf1.ncep_sensible_heat_flux.date;
  fwyf1.ncep_sensible_flux_term.data = (fwyf1.ncep_sensible_heat_flux.data ./ Q0_factor) .* (60*60);

  fwyf1.ncep_bulk_heat_flux_term.date = fwyf1.ncep_bulk_net_heat_flux.date;
  fwyf1.ncep_bulk_heat_flux_term.data = (fwyf1.ncep_bulk_net_heat_flux.data ./ Q0_factor) .* (60*60);
  fwyf1.ncep_bulk_latent_flux_term.date = fwyf1.ncep_bulk_latent_heat_flux.date;
  fwyf1.ncep_bulk_latent_flux_term.data = (fwyf1.ncep_bulk_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  fwyf1.ncep_bulk_sensible_flux_term.date = fwyf1.ncep_bulk_sensible_heat_flux.date;
  fwyf1.ncep_bulk_sensible_flux_term.data = (fwyf1.ncep_bulk_sensible_heat_flux.data ./ Q0_factor) .* (60*60);
disp('Saving fwyf1');
save('../data/fwyf1_thesis.mat','fwyf1');


fclose('all'); close all; clear all; dbclear all;
load('../data/sanf1_thesis.mat');
[wz,az,pz,stz] = station_instrument_heights(sanf1.station_name);
Q0_factor = Q0factor(sanf1.ncep_sea_t.data,[],stz);
  sanf1.ncep_heat_flux_term.date = sanf1.ncep_net_heat_flux.date;
  sanf1.ncep_heat_flux_term.data = (sanf1.ncep_net_heat_flux.data ./ Q0_factor) .* (60*60);
  sanf1.ncep_srf_term.date = sanf1.ncep_srf.date;
  sanf1.ncep_srf_term.data = (sanf1.ncep_srf.data ./ Q0_factor) .* (60*60);
  sanf1.ncep_lrf_term.date = sanf1.ncep_lrf.date;
  sanf1.ncep_lrf_term.data = (sanf1.ncep_lrf.data ./ Q0_factor) .* (60*60);
  sanf1.ncep_latent_flux_term.date = sanf1.ncep_latent_heat_flux.date;
  sanf1.ncep_latent_flux_term.data = (sanf1.ncep_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  sanf1.ncep_sensible_flux_term.date = sanf1.ncep_sensible_heat_flux.date;
  sanf1.ncep_sensible_flux_term.data = (sanf1.ncep_sensible_heat_flux.data ./ Q0_factor) .* (60*60);

  sanf1.ncep_bulk_heat_flux_term.date = sanf1.ncep_bulk_net_heat_flux.date;
  sanf1.ncep_bulk_heat_flux_term.data = (sanf1.ncep_bulk_net_heat_flux.data ./ Q0_factor) .* (60*60);
  sanf1.ncep_bulk_latent_flux_term.date = sanf1.ncep_bulk_latent_heat_flux.date;
  sanf1.ncep_bulk_latent_flux_term.data = (sanf1.ncep_bulk_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  sanf1.ncep_bulk_sensible_flux_term.date = sanf1.ncep_bulk_sensible_heat_flux.date;
  sanf1.ncep_bulk_sensible_flux_term.data = (sanf1.ncep_bulk_sensible_heat_flux.data ./ Q0_factor) .* (60*60);
disp('Saving sanf1');
save('../data/sanf1_thesis.mat','sanf1');


fclose('all'); close all; clear all; dbclear all;
load('../data/lonf1_thesis.mat');
[wz,az,pz,stz] = station_instrument_heights(lonf1.station_name);
Q0_factor = Q0factor(lonf1.ncep_sea_t.data,[],stz);
  lonf1.ncep_heat_flux_term.date = lonf1.ncep_net_heat_flux.date;
  lonf1.ncep_heat_flux_term.data = (lonf1.ncep_net_heat_flux.data ./ Q0_factor) .* (60*60);
  lonf1.ncep_srf_term.date = lonf1.ncep_srf.date;
  lonf1.ncep_srf_term.data = (lonf1.ncep_srf.data ./ Q0_factor) .* (60*60);
  lonf1.ncep_lrf_term.date = lonf1.ncep_lrf.date;
  lonf1.ncep_lrf_term.data = (lonf1.ncep_lrf.data ./ Q0_factor) .* (60*60);
  lonf1.ncep_latent_flux_term.date = lonf1.ncep_latent_heat_flux.date;
  lonf1.ncep_latent_flux_term.data = (lonf1.ncep_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  lonf1.ncep_sensible_flux_term.date = lonf1.ncep_sensible_heat_flux.date;
  lonf1.ncep_sensible_flux_term.data = (lonf1.ncep_sensible_heat_flux.data ./ Q0_factor) .* (60*60);

  lonf1.ncep_bulk_heat_flux_term.date = lonf1.ncep_bulk_net_heat_flux.date;
  lonf1.ncep_bulk_heat_flux_term.data = (lonf1.ncep_bulk_net_heat_flux.data ./ Q0_factor) .* (60*60);
  lonf1.ncep_bulk_latent_flux_term.date = lonf1.ncep_bulk_latent_heat_flux.date;
  lonf1.ncep_bulk_latent_flux_term.data = (lonf1.ncep_bulk_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  lonf1.ncep_bulk_sensible_flux_term.date = lonf1.ncep_bulk_sensible_heat_flux.date;
  lonf1.ncep_bulk_sensible_flux_term.data = (lonf1.ncep_bulk_sensible_heat_flux.data ./ Q0_factor) .* (60*60);
disp('Saving lonf1');
save('../data/lonf1_thesis.mat','lonf1');
