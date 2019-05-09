1;


stn = get_station_from_station_name('looe1');
stn = get_looe1_adcp(stn);
% stn = get_looe1_microcat(stn);
stn = station_optimal_isobath_orientation(stn);

stn = get_gom_hycom(stn);
stn = station_reorient_vectors(stn,'isobath_orientation','gom_hycom_u','gom_hycom_v');
stn.hourly_gom_hycom_u = interp_ts(stn.gom_hycom_u);
stn.hourly_gom_hycom_v = interp_ts(stn.gom_hycom_v);
stn = station_reorient_vectors(stn,'isobath_orientation','hourly_gom_hycom_u','hourly_gom_hycom_v');

stn = get_fkeys_hycom(stn);
stn.hourly_fkeys_hycom_u = interp_ts(stn.fkeys_hycom_u);
stn.hourly_fkeys_hycom_v = interp_ts(stn.fkeys_hycom_v);
stn = station_reorient_vectors(stn,'isobath_orientation','hourly_fkeys_hycom_u','hourly_fkeys_hycom_v');
stn = station_reorient_vectors(stn,'isobath_orientation','fkeys_hycom_u','fkeys_hycom_v');

stn = calc_field_terms(stn,'gom_hycom_seatemp_field','gom_hycom_seatemp','linear',stn.lat,stn.lon);
stn = calc_field_terms(stn,'fkeys_hycom_seatemp_field','fkeys_hycom_seatemp','linear',stn.lat,stn.lon);

scatter_fit_ts(stn.adcp_seatemp,stn.gom_hycom_seatemp,[],[],'ADCP T_s','GOM',[],[],true)
scatter_fit_ts(stn.adcp_seatemp,stn.hourly_gom_hycom_seatemp,,[],[],'ADCP T_s','Hourly GOM',[],[],true)
scatter_fit_ts(stn.adcp_seatemp,stn.fkeys_hycom_seatemp,[],[],'ADCP T_s','FKEYS',[],[],true)
scatter_fit_ts(stn.adcp_seatemp,stn.hourly_fkeys_hycom_seatemp,[],[],'ADCP T_s','Hourly FKEYS',[],[],true)


stn = get_avhrr_weekly_field(stn,true);
stn.hourly_avhrr_weekly_sst = interp_ts(stn.avhrr_weekly_sst);
stn.hourly_avhrr_weekly_sst_x = interp_ts(stn.avhrr_weekly_sst_x);
stn.hourly_avhrr_weekly_sst_y = interp_ts(stn.avhrr_weekly_sst_y);

scatter_fit_ts(stn.adcp_seatemp,stn.avhrr_weekly_sst,[],[],'ADCP T_s','AVHRR',[],[],true)
scatter_fit_ts(stn.adcp_seatemp,stn.hourly_avhrr_weekly_sst,[],[],'ADCP T_s','Hourly AVHRR',[],[],true)


stn = verify_variable(stn,'adcp_u_6_hour_lowpass');
fmg; plot_ts(stn.hourly_fkeys_hycom_xshore,stn.adcp_u_6_hour_lowpass,stn.hourly_gom_hycom_xshore); legend('FKEYS','ADCP 6hlp','GOM'); titlename('Looe Key Cross-shore Current [m/s]');
annotline([],0)

% for yr=2005:2008; for mo=1:1:12; axis([datenum(yr,[mo mo],[1 30]),-.2,.2]); datetick3('x',2,'keeplimits'); pause; end; end;


fn=@(x)(x.data(ismember(get_year(x.date),[2005:2008]))); fmg; spt(3,2,1); hist(fn(stn.adcp_u),1000); title('ADCP U'); xlim([-.2,.2]); grid on; spt(3,2,2); hist(fn(stn.adcp_v),1000); title('ADCP V'); xlim([-1,1]); grid on; spt(3,2,3); hist(fn(stn.hourly_gom_hycom_xshore),1000); title('GOM X 1hr'); xlim([-.2,.2]); grid on; spt(3,2,4); hist(fn(stn.hourly_gom_hycom_lshore),1000); title('GOM L 1hr'); xlim([-1,1]); grid on; spt(3,2,5); hist(fn(stn.hourly_fkeys_hycom_xshore),1000); title('FKEYS X 1hr'); xlim([-.2,.2]); grid on; spt(3,2,6); hist(fn(stn.hourly_fkeys_hycom_lshore),1000); title('FKEYS L 1hr'); xlim([-1,1]); grid on;

for seas=1:4; fn=@(x)(x.data(get_season(x.date)==seas & ismember(get_year(x.date),[2005:2008]))); fmg; suptitlename(num2str(seas)); spt(3,2,1); hist(fn(stn.adcp_u),1000); title('ADCP U'); xlim([-.2,.2]); grid on; spt(3,2,2); hist(fn(stn.adcp_v),1000); title('ADCP V'); xlim([-1,1]); grid on; spt(3,2,3); hist(fn(stn.hourly_gom_hycom_xshore),1000); title('GOM X 1hr'); xlim([-.2,.2]); grid on; spt(3,2,4); hist(fn(stn.hourly_gom_hycom_lshore),1000); title('GOM L 1hr'); xlim([-1,1]); grid on; spt(3,2,5); hist(fn(stn.hourly_fkeys_hycom_xshore),1000); title('FKEYS X 1hr'); xlim([-.2,.2]); grid on; spt(3,2,6); hist(fn(stn.hourly_fkeys_hycom_lshore),1000); title('FKEYS L 1hr'); xlim([-1,1]); grid on; end; clear fn seas



stn = station_reorient_vectors(stn,'isobath_orientation','gom_hycom_seatemp_x','gom_hycom_seatemp_y');
stn = station_reorient_vectors(stn,'isobath_orientation','fkeys_hycom_seatemp_x','fkeys_hycom_seatemp_y');
stn = station_reorient_vectors(stn,'isobath_orientation','avhrr_weekly_sst_x','avhrr_weekly_sst_y');
stn = station_reorient_vectors(stn,'isobath_orientation','hourly_avhrr_weekly_sst_x','hourly_avhrr_weekly_sst_y');

stn.hourly_gom_hycom_seatemp_xshore = interp_ts(stn.gom_hycom_seatemp_xshore);
stn.hourly_gom_hycom_seatemp_lshore = interp_ts(stn.gom_hycom_seatemp_lshore);
stn.hourly_fkeys_hycom_seatemp_xshore = interp_ts(stn.fkeys_hycom_seatemp_xshore);
stn.hourly_fkeys_hycom_seatemp_lshore = interp_ts(stn.fkeys_hycom_seatemp_lshore);

stn.agx = ts_op(stn.hourly_gom_hycom_seatemp_xshore,stn.adcp_u,'*');
stn.afx = ts_op(stn.hourly_fkeys_hycom_seatemp_xshore,stn.adcp_u,'*');
stn.ggx = ts_op(stn.hourly_gom_hycom_seatemp_xshore,stn.hourly_gom_hycom_xshore,'*');
stn.ffx = ts_op(stn.hourly_fkeys_hycom_seatemp_xshore,stn.hourly_fkeys_hycom_xshore,'*');

stn.agl = ts_op(stn.hourly_gom_hycom_seatemp_lshore,stn.adcp_v,'*');
stn.afl = ts_op(stn.hourly_fkeys_hycom_seatemp_lshore,stn.adcp_v,'*');
stn.ggl = ts_op(stn.hourly_gom_hycom_seatemp_lshore,stn.hourly_gom_hycom_lshore,'*');
stn.ffl = ts_op(stn.hourly_fkeys_hycom_seatemp_lshore,stn.hourly_fkeys_hycom_lshore,'*');

stn.aax = ts_op(stn.hourly_avhrr_weekly_sst_xshore,stn.adcp_u,'*');
stn.aal = ts_op(stn.hourly_avhrr_weekly_sst_lshore,stn.adcp_v,'*');

fn=@(x)(3600.*x.data(ismember(get_year(x.date),[2005:2008]))); fmg; suptitlename(['LOOE1 Heat Advection [K/hr]']); spt(3,2,1); hist(fn(stn.agx),1000); title('ADCP U*GOM X'); xlim([-0.015,0.015]); grid on; spt(3,2,2); hist(fn(stn.agl),1000); title('ADCP V*GOM L'); xlim([-0.015,0.015]); grid on; spt(3,2,3); hist(fn(stn.ggx),1000); title('GOM X 1hr*GOM X'); xlim([-0.015,0.015]); grid on; spt(3,2,4); hist(fn(stn.ggl),1000); title('GOM L 1hr*GOM L'); xlim([-0.015,0.015]); grid on; spt(3,2,5); hist(fn(stn.ffx),1000); title('FKEYS X 1hr*FKEYS X'); xlim([-0.015,0.015]); grid on; spt(3,2,6); hist(fn(stn.ffl),1000); title('FKEYS L 1hr*FKEYS L'); xlim([-0.015,0.015]); grid on; clear fn seas

for seas=1:4; fn=@(x)(3600.*x.data(get_season(x.date)==seas & ismember(get_year(x.date),[2005:2008]))); fmg; suptitlename(['LOOE1 Heat Advection [K/hr] Season ' num2str(seas)]); spt(3,2,1); hist(fn(stn.agx),1000); title('ADCP U*GOM X'); xlim([-0.015,0.015]); grid on; spt(3,2,2); hist(fn(stn.agl),1000); title('ADCP V*GOM L'); xlim([-0.015,0.015]); grid on; spt(3,2,3); hist(fn(stn.ggx),1000); title('GOM X 1hr*GOM X'); xlim([-0.015,0.015]); grid on; spt(3,2,4); hist(fn(stn.ggl),1000); title('GOM L 1hr*GOM L'); xlim([-0.015,0.015]); grid on; spt(3,2,5); hist(fn(stn.ffx),1000); title('FKEYS X 1hr*FKEYS X'); xlim([-0.015,0.015]); grid on; spt(3,2,6); hist(fn(stn.ffl),1000); title('FKEYS L 1hr*FKEYS L'); xlim([-0.015,0.015]); grid on; end; clear fn seas


fn=@(x)(3600.*x.data(ismember(get_year(x.date),[2005:2008]))); fmg; suptitlename(['LOOE1 Heat Advection [K/hr]']); spt(3,2,1); hist(fn(stn.aax),1000); title('ADCP U*GOM X'); xlim([-0.015,0.015]); grid on; spt(3,2,2); hist(fn(stn.aal),1000); title('ADCP V*GOM L'); xlim([-0.015,0.015]); grid on; spt(3,2,3); hist(fn(stn.ggx),1000); title('GOM X 1hr*GOM X'); xlim([-0.015,0.015]); grid on; spt(3,2,4); hist(fn(stn.ggl),1000); title('GOM L 1hr*GOM L'); xlim([-0.015,0.015]); grid on; spt(3,2,5); hist(fn(stn.ffx),1000); title('FKEYS X 1hr*FKEYS X'); xlim([-0.015,0.015]); grid on; spt(3,2,6); hist(fn(stn.ffl),1000); title('FKEYS L 1hr*FKEYS L'); xlim([-0.015,0.015]); grid on; clear fn seas

for seas=1:4; fn=@(x)(3600.*x.data(get_season(x.date)==seas & ismember(get_year(x.date),[2005:2008]))); fmg; suptitlename(['LOOE1 Heat Advection [K/hr] Season ' num2str(seas)]); spt(3,2,1); hist(fn(stn.aax),1000); title('ADCP U*GOM X'); xlim([-0.015,0.015]); grid on; spt(3,2,2); hist(fn(stn.aal),1000); title('ADCP V*GOM L'); xlim([-0.015,0.015]); grid on; spt(3,2,3); hist(fn(stn.ggx),1000); title('GOM X 1hr*GOM X'); xlim([-0.015,0.015]); grid on; spt(3,2,4); hist(fn(stn.ggl),1000); title('GOM L 1hr*GOM L'); xlim([-0.015,0.015]); grid on; spt(3,2,5); hist(fn(stn.ffx),1000); title('FKEYS X 1hr*FKEYS X'); xlim([-0.015,0.015]); grid on; spt(3,2,6); hist(fn(stn.ffl),1000); title('FKEYS L 1hr*FKEYS L'); xlim([-0.015,0.015]); grid on; end; clear fn seas


fnd=@(x)(3600.*x.data(ismember(get_year(x.date),[2005:2008])));
fnt=@(x)(x.date(ismember(get_year(x.date),[2005:2008])));
fmg; spt(1,3,1); boxplot(fnd(stn.aax),get_season(fnt(stn.aax)),'notch','on'); title('ADCP U*AVHRR X'); ylim([-1,1]*3600*1e-4); spt(1,3,2); boxplot(fnd(stn.ggx),get_season(fnt(stn.ggx)),'notch','on'); title('GOM U*GOM X'); ylim([-1,1]*3600*1e-4); spt(1,3,3); boxplot(fnd(stn.ffx),get_season(fnt(stn.ffx)),'notch','on'); title('FKEYS U*FKEYS X'); ylim([-1,1]*3600*1e-4);


% Meteorology and cross-analysis

smkf1 = get_station_from_station_name('smkf1');
smkf1 = load_all_ndbc_data(smkf1);
smkf1 = load_station_data(smkf1);

fmg; plot_ts(smkf1.ndbc_sea_t,smkf1.ct_shallow_seatemp,smkf1.ndbc_ct_shallow_seatemp);

smkf1 = verify_variable(smkf1,'ndbc_wind1_u');
smkf1 = verify_variable(smkf1,'ndbc_wind1_v');
smkf1 = station_reorient_vectors(smkf1,smkf1.isobath_orientation,'ndbc_wind1_u','ndbc_wind1_v');

smkf1 = get_avhrr_weekly_field(smkf1,true);
smkf1.avhrr_weekly_sst.date = smkf1.avhrr_weekly_sst_field.date;
smkf1.avhrr_weekly_sst.data = interp_field(smkf1.avhrr_weekly_sst_field.lat,smkf1.avhrr_weekly_sst_field.lon,smkf1.avhrr_weekly_sst_field.field,smkf1.lat,smkf1.lon,'linear');
smkf1.hourly_avhrr_weekly_sst = interp_ts(smkf1.avhrr_weekly_sst);
scatter_fit_ts(smkf1.ct_shallow_seatemp,smkf1.ndbc_sea_t,[],[],'CT','SMKF1 T_s',[],[],true)
scatter_fit_ts(smkf1.avhrr_weekly_sst,smkf1.ndbc_sea_t,[],[],'AVHRR','SMKF1 T_s',[],[],true)
scatter_fit_ts(smkf1.hourly_avhrr_weekly_sst,smkf1.ndbc_sea_t,[],[],'Hourly AVHRR','SMKF1 T_s',[],[],true)

scatter_fit_ts(stn.adcp_seatemp,smkf1.ndbc_sea_t,[],[],'LOOE ADCP','SMKF1 T_s',[],[],true)




stn = get_ww3_station(stn);
shi=find(stn.ww3_sigwavehgt.data>=2);
scatter_fit_ts(stn.ww3_sigwavehgt,stn.adcp_sfc_speed,shi,[],'WW3 H_s (>2.5m)','Upper 5m ADCP 40hlp')


stn = verify_variable(stn,'adcp_u_40_hour_lowpass');
stn = verify_variable(stn,'adcp_v_40_hour_lowpass');
scatter_fit_ts(smkf1.ndbc_wind1_xshore,stn.adcp_u_40_hour_lowpass,@(x)(find(x.data<-20)),[],'W xshore (<-20)','ADCP xshore 40hlp')

shi=find(stn.ww3_sigwavehgt.data>=2);
scatter_fit_ts(stn.ww3_sigwavehgt,stn.adcp_sfc_speed_40_hour_lowpass,shi,[],'WW3 Sig wave hgt (>2m)','Near-surface ADCP 40hlp')
scatter_fit_ts(stn.ww3_sigwavehgt,stn.adcp_speed_40_hour_lowpass,shi,[],'WW3 Sig wave hgt (>2m)','ADCP 40hlp')
% print('-dtiff','../figs/looe1-adcp_speed_40_hour_lowpass-scatter-ww3_sigwavehgt.tiff');
