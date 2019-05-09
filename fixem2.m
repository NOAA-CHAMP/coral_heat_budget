1;

% ori = 52.25;
% ori = 52.4;
ori = 52.6;

stn = station_reorient_vectors(stn,ori,'gom_hycom_u','gom_hycom_v','gom_hycom_x','gom_hycom_l');
multiplot_station(stn,{'gom_hycom_u','gom_hycom_v','gom_hycom_x','gom_hycom_l'},[],[],[],[],[-1,1.5]);
titlename([stn.station_name ' ' num2str(ori) '^o']);

stn.gom_hycom_seatemp_x.date = stn.gom_hycom_seatemp_field.date;
stn.gom_hycom_seatemp_x.data = interp_field(stn.gom_hycom_seatemp_field.lat,stn.gom_hycom_seatemp_field.lon,stn.gom_hycom_seatemp_field.gradient_x,stn.lat,stn.lon);

stn.gom_hycom_seatemp_y.date = stn.gom_hycom_seatemp_field.date;
stn.gom_hycom_seatemp_y.data = interp_field(stn.gom_hycom_seatemp_field.lat,stn.gom_hycom_seatemp_field.lon,stn.gom_hycom_seatemp_field.gradient_y,stn.lat,stn.lon);

stn = station_reorient_vectors(stn,ori,'gom_hycom_seatemp_x','gom_hycom_seatemp_y','gom_hycom_seatemp_cross','gom_hycom_seatemp_along');
multiplot_station(stn,{'gom_hycom_seatemp_x','gom_hycom_seatemp_y','gom_hycom_seatemp_cross','gom_hycom_seatemp_along'},[],[],[],[],[-1.5e-3,2e-3]);
titlename([stn.station_name ' ' num2str(ori) '^o']);

stn.gx.date=stn.gom_hycom_seatemp_cross.date;
stn.gx.data=-6*3600.*stn.gom_hycom_seatemp_cross.data.*stn.gom_hycom_x.data;

stn.gl.date=stn.gom_hycom_seatemp_along.date;
stn.gl.data=-6*3600.*stn.gom_hycom_seatemp_along.data.*stn.gom_hycom_l.data;

stn.ga.date=stn.gx.date;
stn.ga.data=stn.gx.data+stn.gl.data;

% fmg; plot_ts(stn.gx,stn.gl,stn.ga);

% fmg; plot(stn.gx.date,[cumsum(stn.gx.data),cumsum(stn.gl.data),cumsum(stn.ga.data)]); datetick3;
% legend('\partial_xu\times\partial_xT','\partial_lu\times\partial_lT','u^.\nablaT');
% titlename([stn.station_name ' ' num2str(ori) '^o']);

[gix,ix] = intersect_dates(stn.ga.date,stn.benthic_ndbc_erai_30a_net_heat_flux_term.date);
fmg; plot(stn.gx.date,[cumsum(stn.gx.data),cumsum(stn.gl.data),cumsum(stn.ga.data)],stn.benthic_ndbc_erai_30a_net_heat_flux_term.date(ix(1):ix(end)),cumsum(stn.benthic_ndbc_erai_30a_net_heat_flux_term.data(ix(1):ix(end))),'y:'); datetick3;
legend('\partial_xu\partial_xT','\partial_lu\partial_lT','u^.\nablaT','(Q_0+Q_b)/\rhoC_ph');
titlename([stn.station_name ' ' num2str(ori) '^o']);

fmg; hold on; plot(stn.gx.date,[cumsum(stn.gx.data)],stn.benthic_ndbc_erai_30a_net_heat_flux_term.date(ix(1):ix(end)),cumsum(stn.benthic_ndbc_erai_30a_net_heat_flux_term.data(ix(1):ix(end))),'y:'); datetick3;
stn.ad = ts_op(stn.gx,stn.benthic_ndbc_erai_30a_net_heat_flux_term,'+');
plot(stn.ad.date,cumsum(stn.ad.data),'r');
legend('\partial_xu\partial_xT','(Q_0+Q_b)/\rhoC_ph','\partial_xu\partial_xT + (Q_0+Q_b)/\rhoC_ph');
titlename([stn.station_name ' ' num2str(ori) '^o']);

return;
