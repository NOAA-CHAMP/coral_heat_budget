1;

disp(stnm);

figspath = get_thesis_path('../figs');

stn = get_station_from_station_name(stnm);
stn = load_all_ndbc_data(stn);
stn = verify_variable(stn,'ndbc_sea_t_1_day_average');
nix = find(stn.ndbc_sea_t_1_day_average);

stn = get_fkeys_hycom(stn);
stn = calc_field_terms(stn,'fkeys_hycom_seatemp_field','fkeys_hycom_seatemp','linear',[],[],[],true);
stn = verify_variable(stn,'fkeys_hycom_seatemp_1_day_average');
kix = find(stn.fkeys_hycom_seatemp_1_day_average);

stn = get_gom_hycom(stn);
stn = calc_field_terms(stn,'gom_hycom_seatemp_field','gom_hycom_seatemp','linear',[],[],[],true);

stn = get_avhrr_weekly_field(stn,true);

scatter_fit_ts(stn.ndbc_sea_t,stn.gom_hycom_seatemp,[],[],'NDBC','GoM',[],[],true), axis([15,33,15,33]);
print('-dtiff',fullfile(figspath,[lower(stn.station_name) '-ndbc_sea_t-scatter-gom.tiff']));
scatter_fit_ts(stn.ndbc_sea_t,stn.fkeys_hycom_seatemp,[],[],'NDBC','FKEYS',[],[],true), axis([15,33,15,33]);
print('-dtiff',fullfile(figspath,[lower(stn.station_name) '-ndbc_sea_t-scatter-fkeys.tiff']));
scatter_fit_ts(stn.ndbc_sea_t,stn.avhrr_weekly_sst,[],[],'NDBC','AVHRR',[],[],true), axis([15,33,15,33]);
print('-dtiff',fullfile(figspath,[lower(stn.station_name) '-ndbc_sea_t-scatter-avhrr.tiff']));

if (0)
scatter_fit_ts(stn.avhrr_weekly_sst_x,stn.fkeys_hycom_seatemp_x,[],[],'AVHRR','FKEYS \partial_xT',[],[],true), axis([-6,10,-6,10]*1e-4);
scatter_fit_ts(stn.avhrr_weekly_sst_y,stn.fkeys_hycom_seatemp_y,[],[],'AVHRR','FKEYS \partial_yT',[],[],true), axis([-6,10,-6,10]*1e-4);

scatter_fit_ts(stn.avhrr_weekly_sst_x,stn.gom_hycom_seatemp_x,[],[],'AVHRR','GoM \partial_xT',[],[],true), axis([-6,10,-6,10]*1e-4);
scatter_fit_ts(stn.avhrr_weekly_sst_y,stn.gom_hycom_seatemp_y,[],[],'AVHRR','GoM \partial_yT',[],[],true), axis([-6,10,-6,10]*1e-4);
end;
