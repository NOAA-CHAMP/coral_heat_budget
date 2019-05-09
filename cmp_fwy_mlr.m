1;

stn.ndbc_air_sea_t = ts_op(stn.ndbc_air_t,stn.ndbc_sea_t,'-');
stn.erai_ndbc_air_sea_spechumid = ts_op(stn.erai_spechumid,stn.ndbc_sea_spechumid,'-');
stn = station_reorient_vectors(stn,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');

fwyf1.ndbc_air_sea_t = ts_op(fwyf1.ndbc_air_t,fwyf1.ndbc_sea_t,'-');
fwyf1.erai_ndbc_air_sea_spechumid = ts_op(fwyf1.erai_spechumid,fwyf1.ndbc_sea_spechumid,'-');
fwyf1 = station_reorient_vectors(fwyf1,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');

for cfld={...
    'ndbc_sea_t',...
    'ndbc_air_t',...
    'ndbc_air_sea_t',...
    'ndbc_wind1_xshore',...
    'ndbc_wind1_lshore',...
    'erai_ndbc_air_sea_spechumid',...
    'ndbc_erai_erai_30a_sensible_flux',...
    'ndbc_erai_erai_30a_latent_flux',...
    'ndbc_erai_erai_30a_net_flux',...
    'erai_avhrr_advected_heat',...
    'avhrr_weekly_diffused_heat_flux',...
    'ndbc_erai_erai_30a_avhrr_qtadv_flux',...
    'b_ndbc_erai_erai_30a_avhrr_dt',...
    'ndbc_erai_erai_30a_avhrr_hc_dTdthc',...
    'ndbc_erai_erai_30a_avhrr_hc_dTdt_flux',...
        };
  fld = cfld{:};
  scatter_fit_ts_seasons(stn.(fld),fwyf1.(fld),[],[],strrep(['MLRF1 ',fld],'_','\_'),'FWYF1',[],[],true),
  linkaxes
  for ix=1:4; spt(2,2,ix); legend('Location','Best'); end;
  pause(0.5);
end;
