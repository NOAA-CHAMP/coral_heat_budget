1;

lonf1=get_station_from_station_name('lonf1'); lonf1=load_all_ndbc_data(lonf1); lonf1=load_station_data(lonf1);
[lonf1.ts,lonf1.ta,lonf1.u10]=intersect_tses(lonf1.ndbc_sea_t,lonf1.ndbc_air_t,lonf1.ndbc_wind1_speed);
find_date_ranges(lonf1.ts.date,(3.1/24))
lonf1.ts.date=lonf1.ts.date(98983:116959); lonf1.ts.data=lonf1.ts.data(98983:116959);
[lonf1.ts,lonf1.ta,lonf1.u10]=intersect_tses(lonf1.ts,lonf1.ta,lonf1.u10);
find_date_ranges(lonf1.ts.date,(1.1/24))
plot_spec(lonf1,{'u10','ta','ts'},[],[],[],[1e-4,1e7]);

smkf1=get_station_from_station_name('smkf1'); smkf1=load_all_ndbc_data(smkf1); smkf1=load_station_data(smkf1);
[smkf1.ts,smkf1.ta,smkf1.u10]=intersect_tses(smkf1.ndbc_sea_t,smkf1.ndbc_air_t,smkf1.ndbc_wind1_speed);
find_date_ranges(smkf1.ts.date,(3.1/24))
[smkf1.ts,smkf1.ta,smkf1.u10]=intersect_tses(smkf1.ndbc_sea_t,smkf1.ndbc_air_t,smkf1.ndbc_wind1_speed);
%%smkf1.ts.date=smkf1.ts.date(45224:53636); smkf1.ts.data=smkf1.ts.data(45224:53636);
smkf1.ts.date=smkf1.ts.date(138259:148421); smkf1.ts.data=smkf1.ts.data(138259:148421);
[smkf1.ts,smkf1.ta,smkf1.u10]=intersect_tses(smkf1.ts,smkf1.ta,smkf1.u10);
plot_spec(smkf1,{'u10','ta','ts'},[],[],[],[1e-4,1e7]);

looe1=get_station_from_station_name('looe1'); looe1=get_looe1_microcat(looe1); looe1=get_looe1_adcp(looe1);
plot_spec(looe1,{'adcp_baroclinic_btm_x','adcp_btm_x','adcp_btm_l'},[],[],[],[]);
[looe1.bc_btm_x,looe1.btm_l]=intersect_tses(looe1.adcp_baroclinic_btm_x,looe1.adcp_btm_l);
find_date_ranges(looe1.bc_btm_x.date,(3.1/24))
looe1.bc_btm_x.date=looe1.bc_btm_x.date(1959:10508); looe1.bc_btm_x.data=looe1.bc_btm_x.data(1959:10508);
[looe1.bc_btm_x,looe1.btm_l]=intersect_tses(looe1.bc_btm_x,looe1.btm_l);
find_date_ranges(looe1.bc_btm_x.date,(1.1/24))
plot_spec(looe1,{'bc_btm_x','btm_l'},[],[],[],[]);
[looe1.bc_btm_x,looe1.btm_x,looe1.btm_l]=intersect_tses(looe1.bc_btm_x,looe1.btm_x,looe1.btm_l);
plot_spec(looe1,{'bc_btm_x','btm_x','btm_l'},[],[],[],[]);
plot_spec(looe1,{'adcp_baroclinic_sfc_x','adcp_sfc_x','adcp_sfc_l'},[],[],[],[]);

mlrf1=get_station_from_station_name('mlrf1'); mlrf1=load_all_ndbc_data(mlrf1); mlrf1=load_station_data(mlrf1);
plot_spec(mlrf1,{'ndbc_wind1_speed','ndbc_air_t','ndbc_sea_t'},[],[],[],[1e-4,1e7]);

looe1.adcpseatemp_ndbc_erai_erai_30a_wind_stress=load('../data/looe1_adcpseatemp_ndbc_erai_erai_30a_wind_stress.mat');
looe1=verify_variable(looe1,'adcpseatemp_ndbc_erai_erai_30a_wind_stress_1_d_sum');

fmg; plot_ts(looe1.adcpseatemp_ndbc_erai_erai_30a_wind_stress,looe1.adcp_baroclinic_btm_x,looe1.adcp_baroclinic_sfc_x); legend('\tau','BC Btm X','BC Sfc X');
for yr=2006:2008; for dt=datenum(yr,[2:3:12],15); disp(datestr(dt)); figure(26); axis([dt-50,dt+50,-.4,.4]); datetick3; pause; end; end;

dryf1=get_station_from_station_name('dryf1'); dryf1=load_all_ndbc_data(dryf1); dryf1=load_station_data(dryf1);
plot_spec(dryf1,{'ndbc_wind1_speed','ndbc_air_t','ndbc_sea_t'},[],[],[],[1e-4,1e7],'tiff');
find_date_ranges(dryf1.ndbc_sea_t.date,(8.1/24))
[dryf1.ts,dryf1.ta,dryf1.u10]=intersect_tses(dryf1.ndbc_sea_t,dryf1.ndbc_air_t,dryf1.ndbc_wind1_speed);
dryf1.ts.date=dryf1.ts.date(34232:57289); dryf1.ts.data=dryf1.ts.data(34232:57289);
[dryf1.ts,dryf1.ta,dryf1.u10]=intersect_tses(dryf1.ts,dryf1.ta,dryf1.u10);
plot_spec(dryf1,{'u10','ta','ts'},[],[],[],[1e-4,1e7]);
