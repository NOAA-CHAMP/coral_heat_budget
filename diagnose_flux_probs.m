1;

clims_hft=anclim(smkf1,'ncep_heat_flux_term',@nansum,@sum);
figure; plot(1:366, [clims_sea_t.hrclim ; clims_hft.hrclim]); maxigraph;
figure; plot(1:366, [clims_sea_t.hrclim , clims_hft.hrclim]); maxigraph;
figure; plot(1:(366*24), [clims_sea_t.hrclim , clims_hft.hrclim]); maxigraph;
figure; plot(1:(366*24), [clims_sea_t.hrclim ; clims_hft.hrclim]); maxigraph;
figure; plot(1:(366), [clims_sea_t.dyclim ; clims_hft.dyclim]); maxigraph;
figure; plot(1:(366*24), [clims_sea_t.hrclim ; clims_hft.hrclim]); maxigraph;
figure; plot(1:(366*24), [clims_sea_t.hranom ; clims_hft.hranom]); maxigraph;
figure; plot(1:numel(clims_sea_t.hranom), [clims_sea_t.hranom ; clims_hft.hranom]); maxigraph;
size(clims_sea_t.hranom)
figure; plot(1:(366*24), [clims_sea_t.hranom(1,:) ; clims_hft.hranom(1,:)]); maxigraph;
figure; plot(1:(366), [clims_sea_t.dyanom(1,:) ; clims_hft.dyanom(1,:)]); maxigraph;
clims_hft=anclim(smkf1,'ncep_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyanom(1,:) ; clims_hft.dyanom(1,:)]); maxigraph;
figure; plot(1:(366), [clims_sea_t.dyanom(1,:) ; clims_hft.dyanom(1,:)]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyanom(2,:) ; clims_hft.dyanom(2,:)]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyclim ; clims_hft.dyclim]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ncep_heat_flux_term',@nansum,@sum);
figure; plot(1:(366), [clims_sea_t.dyclim ; clims_hft.dyclim]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ncep_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyanom ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyclim ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366*24), [clims_sea_t.hrclim-clims_sea_t.hrclim(1) ; cumsum(clims_hft.hrclim)]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ncep_bulk_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ncep_bulk_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ncep_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dymean(1,:)-clims_sea_t.dymean(1,1) ; cumsum(clims_hft.dymean(1,:))]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dymean(2,:)-clims_sea_t.dymean(2,1) ; cumsum(clims_hft.dymean(2,:))]); maxigraph; legend('SeaT','Q_0');
clims_sea_t.dymean(:,1)
yr
yr=3; figure; plot(1:(366), [clims_sea_t.dymean(yr,:)-clims_sea_t.dymean(yr,1) ; cumsum(clims_hft.dymean(yr,:))]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ncep_bulk_heat_flux_term');
yr=3; figure; plot(1:(366), [clims_sea_t.dymean(yr,:)-clims_sea_t.dymean(yr,1) ; cumsum(clims_hft.dymean(yr,:))]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ncep_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
ylim([-2 14])
smkf1 = station_heat_flux(smkf1,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid','ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lrf','ndbc');
clims_hft=anclim(smkf1,'ndbc_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
smkf1.ndbc_heat_flux_term
smkf1.ndbc_heat_flux_term.data
smkf1.ndbc_bulk_heat_flux_term.data
smkf1.ncep_heat_flux_term.data
smkf1.ndbc_bulk_heat_flux_term.data
smkf1.ndbc_heat_flux_term.data
scatter_fit_station(smkf1,'ndbc_relhumid','ncep_relhumid');
ylim([0 100]); xlim([0 100]);
help scatter_fit_ts
scatter_fit_ts(smkf1.ndbc_relhumid,smkf1.ncep_relhumid,@(x)(find(x.data < 98)),@(x)(find(x.data < 98)));
scatter_fit_ts(smkf1.ndbc_relhumid,smkf1.ncep_relhumid,@(x)(find(x.data < 93)),@(x)(find(x.data < 93)));
ylim([0 100]); xlim([0 100]);
scatter_fit_ts(smkf1.ndbc_relhumid,smkf1.ncep_relhumid,@(x)(find(x.data < 93)),[]);
xlim([0 100]); ylim([0 100]);
smkf1
smkf1 = rmfield(smkf1,{'ndbc_wind_stress ndbc_sensible_heat_flux ndbc_latent_heat_flux ndbc_net_heat_flux ndbc_heat_flux_term
smkf1 = rmfield(smkf1,{'ndbc_wind_stress','ndbc_sensible_heat_flux','ndbc_latent_heat_flux','ndbc_net_heat_flux','ndbc_heat_flux_term'});
smkf1
smkf1 = station_heat_flux(smkf1,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid','ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lrf','ncep_bulk_rh');
clims_hft=anclim(smkf1,'ncep_bulk_rh_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim).*1e-3]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim).*1e3]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim).*3e3]); maxigraph; legend('SeaT','Q_0');
size(clims_hft.yrs)
clims_hft.yrs
clims_sea_t.yrs'
figure; plot(1:(366), [clims_sea_t.dymean(end,:)-clims_sea_t.dymean(end,1) ; cumsum(clims_hft.dymean(end,:))]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ncep_bulk_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dymean(end,:)-clims_sea_t.dymean(end,1) ; cumsum(clims_hft.dymean(end,:))]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dymean(end-1,:)-clims_sea_t.dymean(end-1,1) ; cumsum(clims_hft.dymean(end-1,:))]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ncep_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dymean(end-1,:)-clims_sea_t.dymean(end-1,1) ; cumsum(clims_hft.dymean(end-1,:))]); maxigraph; legend('SeaT','Q_0');
clims_hft=anclim(smkf1,'ndbc_bulk_heat_flux_term');
clims_hft.yrs
clims_sea_t.yrs'
figure; plot(1:(366), [clims_sea_t.dymean(15,:)-clims_sea_t.dymean(15,1) ; cumsum(clims_hft.dymean(2,:))]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dymean(14,:)-clims_sea_t.dymean(14,1) ; cumsum(clims_hft.dymean(1,:))]); maxigraph; legend('SeaT','Q_0');
figure; plot(1:(366), [clims_sea_t.dymean(16,:)-clims_sea_t.dymean(16,1) ; cumsum(clims_hft.dymean(3,:))]); maxigraph; legend('SeaT','Q_0');
clims_hft.dymean
cumsum([1 2 3 4 5])
cumsum([1 nan 3 4 5])
help nancumsum
help nansum
help nanvar
smkf1
smkf1 = rmfield(smkf1,{'ncep_bulk_rh_wind_stress','ncep_bulk_rh_sensible_heat_flux','ncep_bulk_rh_latent_heat_flux','ncep_bulk_rh_net_heat_flux','ncep_bulk_rh_heat_flux_term'});
smkf1
save('../data/smkf1_thesis.mat');
clims_hft=anclim(smkf1,'ncep_heat_flux_term');
smkf1 = station_heat_flux(smkf1,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid','ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lrf','ncep_bulk_rh');
clims_hft=anclim(smkf1,'ncep_bulk_rh_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
titlename('Sombrero Key daily climatologies: T_s_e_a vs. Q_0 (NCEP Bulk, in situ RH)');
print('-dpng','../figs/smkf1_sea_t_and_ncep_bulk_rh_forcing_climatology.png');
clims_hft=anclim(smkf1,'ncep_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
titlename('Sombrero Key daily climatologies: T_s_e_a vs. Q_0 (pure NCEP)');
print('-dpng','../figs/smkf1_sea_t_and_ncep_forcing_climatology.png');
clims_hft=anclim(smkf1,'ncep_bulk_heat_flux_term');
figure; plot(1:(366), [clims_sea_t.dyclim-clims_sea_t.dyclim(1) ; cumsum(clims_hft.dyclim)]); maxigraph; legend('SeaT','Q_0');
titlename('Sombrero Key daily climatologies: T_s_e_a vs. Q_0 (NCEP Bulk)');
print('-dpng','../figs/smkf1_sea_t_and_ncep_bulk_forcing_climatology.png');
ylim([-2 14])
multiplot_station(smkf1,{'ncep_bulk_rh_srf','ncep_bulk_rh_lrf','ncep_bulk_rh_sensible_heat_flux','ncep_bulk_rh_latent_heat_flux'});
multiplot_station(smkf1,{'ncep_srf','ncep_lrf','ncep_bulk_rh_sensible_heat_flux','ncep_bulk_rh_latent_heat_flux'});
xlim([datenum(2006,1,1) datenum(2006,12,31,23,59,59)])
multiplot_station(smkf1,{'ncep_srf','ncep_lrf','ndbc_bulk_lrf','ncep_bulk_rh_sensible_heat_flux','ncep_bulk_rh_latent_heat_flux'});
xlim([datenum(2001,1,1) datenum(2001,12,31,23,59,59)])
scatter_fit_ts(smkf1.ndbc_bulk_cloud_cover,smkf1.ncep_cloud_cover);
smkf1
scatter_fit_ts(smkf1.ndbc_cloud_cover,smkf1.ncep_cloud_cover);
xlim([datenum(2006,1,1) datenum(2006,12,31,23,59,59)])
datetick3
multiplot_station(smkf1,{'ncep_srf','ncep_lrf','ndbc_bulk_lrf','ncep_bulk_rh_sensible_heat_flux','ncep_bulk_rh_latent_heat_flux'});
multiplot_station(smkf1,{'ncep_srf','ncep_lrf','ncep_sensible_heat_flux','ncep_latent_heat_flux','ncep_bulk_rh_sensible_heat_flux','ncep_bulk_rh_latent_heat_flux'});
xlim([datenum(2006,1,1) datenum(2006,12,31,23,59,59)])
titlename('Sombrero Key comparative heat flux terms');
print('-dpng','../figs/smkf1_comparative_forcing_2006.png');
smkf1 = verify_variable(smkf1,'ndbc_sea_t_40_hour_highpass');
fclose('all'); close all; clear all; dbclear all; pack;
smkf1 = calcq0('smkf1');
cns = anclim(stn,'ncep_sensible_heat_flux');
cns = anclim(smkf1,'ncep_sensible_heat_flux');
cnbs = anclim(smkf1,'ncep_bulk_sensible_heat_flux');
cnbs
cbs
cbs = anclim(smkf1,'ndbc_bulk_sensible_heat_flux');
crs = anclim(smkf1,'ndbc_bulk_rh_sensible_heat_flux');
figure; plot(1:(366), [cumsum(cns.dyclim) ; cumsum(crs.dyclim)]); maxigraph; legend('Q_0 NCEP','Q_0 NCEP Bulk RH');
clear cnbs cbs
figure; plot(1:(366), [cumsum(cns.dyclim) ; cumsum(crs.dyclim)]); maxigraph; legend('Q_0 NCEP','Q_0 NCEP Bulk RH');
cnl = anclim(smkf1,'ncep_latent_heat_flux');
crl
crl = anclim(smkf1,'ndbc_bulk_rh_latent_heat_flux');
figure; plot(1:(366), [cumsum(cnl.dyclim) ; cumsum(crl.dyclim)]); maxigraph; legend('Q_L_H NCEP','Q_L_H NCEP Bulk RH');
Q0f = Q0factor(smkf1.ndbc_sea_t.data,[],2);
Q0f
400/Q0f
ans*3600
figure; plot(1:(366), [cumsum(cnl.dyclim)/Q0f ; cumsum(crl.dyclim)/Q0f]); maxigraph; legend('Q_L_H NCEP','Q_L_H NCEP Bulk RH');
figure; plot(1:(366), [cumsum(cnl.dyclim) ; cumsum(crl.dyclim)]); maxigraph; legend('Q_L_H NCEP','Q_L_H NCEP Bulk RH');
figure; plot(1:(366), [cumsum(cnl.dyclim)/Q0f ; cumsum(crl.dyclim)/Q0f]); maxigraph; legend('Q_L_H NCEP','Q_L_H NCEP Bulk RH');
figure; plot(1:(366), [cumsum(cnl.dyclim)*3600/Q0f ; cumsum(crl.dyclim)*3600/Q0f]); maxigraph; legend('Q_L_H NCEP','Q_L_H NCEP Bulk RH');
figure; plot(1:(366), [cumsum(cns.dyclim)*3600/Q0f ; cumsum(crs.dyclim)*3600/Q0f]); maxigraph; legend('Q_S_H NCEP','Q_S_H NCEP Bulk RH');
ct = anclim(smkf1,'ndbc_sea_t');
cr=anclim(smkf1,'ndbc_bulk_rh_heat_flux_term');
cn=anclim(smkf1,'ncep_heat_flux_term');
figure; plot(1:(366), [ct.dyclim-ct.dyclim(1) ; cumsum(cn.dyclim) ; cumsum(cr.dyclim)]); maxigraph; legend('SeaT','Q_0 NCEP','Q_0 Bulk RH');
crw
crlw = anclim(smkf1,'ndbc_bulk_rh_lrw');
crlw = anclim(smkf1,'ndbc_bulk_rh_lrf');
cnlw = anclim(smkf1,'ncep_lrf');
figure; plot(1:(366), [cumsum(cnlw.dyclim)*3600/Q0f ; cumsum(crlw.dyclim)*3600/Q0f]); maxigraph; legend('Q_L_W NCEP','Q_L_W NCEP Bulk RH');
nsh
nsh = anclim(smkf1,'ncep_spechumid');
bsh = anclim(smkf1,'ndbc_spechumid');
figure; plot(1:(366), [nsh.dyclim ; bsh.dyclim ]); maxigraph; legend('NCEP q_s','Bulk q_s');
x=cumsum(cn.dyclim); y=cumsum(cr.dyclim); figure; plot(1:(366), [ct.dyclim-ct.dyclim(209) ; x-x(209) ; y-y(209) ]); maxigraph; legend('SeaT','Q_0 NCEP','Q_0 Bulk RH');
