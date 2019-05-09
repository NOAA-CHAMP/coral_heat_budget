1;
% SCRIPT to calculate values for Appendix II of Lew Gramer's Dissertation
%
% Last Saved Time-stamp: <Mon 2013-07-01 17:56:43 Eastern Daylight Time gramer>

npts = 5;
grid_interp_method = 'linear';
% EXCEPTION: SANF1
%npts = 3;
%grid_interp_method = 'nearest';
% EXCEPTION: DRYF1
%npts = 3;
%grid_interp_method = {'nanmean',2,2};


if ( ~exist('fwyf1','var') )
  fwyf1 = get_station_from_station_name('fwyf1');
  fwyf1 = load_all_ndbc_data(fwyf1);

  fwyf1 = get_erai_station(fwyf1);
  fwyf1 = adjust_erai_station(fwyf1);
  fwyf1.erai_dsrf_daily_dose = par_dose(fwyf1.erai_dsrf);
  fwyf1.erai_dlrf_daily_dose = par_dose(fwyf1.erai_dlrf);
  fwyf1.erai_dsrf_adj_daily_dose = par_dose(fwyf1.erai_dsrf_adj);
  fwyf1.erai_dlrf_adj_daily_dose = par_dose(fwyf1.erai_dlrf_adj);

  fwyf1 = get_ncep_station(fwyf1,'narr');
  fwyf1.ncep_dsrf_daily_dose = par_dose(fwyf1.ncep_dsrf);
  fwyf1.ncep_dlrf_daily_dose = par_dose(fwyf1.ncep_dlrf);

  fwyf1 = get_avhrr_weekly_field(fwyf1,true,grid_interp_method,npts);
end;

if ( ~exist('mlrf1','var') )
  mlrf1 = get_station_from_station_name('mlrf1');
  mlrf1 = station_optimal_isobath_orientation(mlrf1);
  mlrf1 = load_all_ndbc_data(mlrf1);
  mlrf1 = load_station_data(mlrf1);

  mlrf1 = get_erai_station(mlrf1);
  mlrf1 = adjust_erai_station(mlrf1);
  mlrf1 = adjust_erai_station_waves(mlrf1);
  mlrf1.erai_dsrf_daily_dose = par_dose(mlrf1.erai_dsrf);
  mlrf1.erai_dlrf_daily_dose = par_dose(mlrf1.erai_dlrf);
  mlrf1.erai_dsrf_adj_daily_dose = par_dose(mlrf1.erai_dsrf_adj);
  mlrf1.erai_dlrf_adj_daily_dose = par_dose(mlrf1.erai_dlrf_adj);

  mlrf1 = get_fkeys_hycom(mlrf1);
  mlrf1 = calc_field_terms(mlrf1,'fkeys_hycom_seatemp_field','fkeys_hycom_seatemp',mlrf1.fkeys_hycom_interp_method,[],[],npts);
  mlrf1 = station_reorient_vectors(mlrf1,'isobath_orientation','fkeys_hycom_seatemp_x','fkeys_hycom_seatemp_y');

  mlrf1 = get_gom_hycom(mlrf1);
  mlrf1 = calc_field_terms(mlrf1,'gom_hycom_seatemp_field','gom_hycom_seatemp',mlrf1.gom_hycom_interp_method,[],[],npts);
  mlrf1 = station_reorient_vectors(mlrf1,'isobath_orientation','gom_hycom_seatemp_x','gom_hycom_seatemp_y');

  mlrf1 = get_avhrr_weekly_field(mlrf1,true,grid_interp_method,npts);
  mlrf1 = station_reorient_vectors(mlrf1,'isobath_orientation','hourly_avhrr_weekly_sst_x','hourly_avhrr_weekly_sst_y');
end;

if ( ~exist('lonf1','var') )
  lonf1 = get_station_from_station_name('lonf1');
  lonf1 = station_optimal_isobath_orientation(lonf1);
  lonf1 = load_all_ndbc_data(lonf1);
  lonf1 = station_dewp_to_relhumid(lonf1,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
  lonf1 = station_relhumid_to_spechumid(lonf1,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');

  lonf1 = get_erai_station(lonf1);
  lonf1 = adjust_erai_station(lonf1);

  lonf1 = get_avhrr_weekly_field(lonf1,true,grid_interp_method,npts);
  lonf1 = station_reorient_vectors(lonf1,'isobath_orientation','hourly_avhrr_weekly_sst_x','hourly_avhrr_weekly_sst_y');
end;

if ( ~exist('smkf1','var') )
  smkf1 = get_station_from_station_name('smkf1');
  smkf1 = load_all_ndbc_data(smkf1);
  smkf1 = station_dewp_to_relhumid(smkf1,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
  smkf1 = station_relhumid_to_spechumid(smkf1,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
  smkf1 = get_erai_station(smkf1);
  smkf1 = adjust_erai_station(smkf1);
end;

if ( ~exist('looe1','var') )
  looe1 = get_station_from_station_name('looe1');
  looe1 = station_optimal_isobath_orientation(looe1);
  looe1 = get_looe1_adcp(looe1);
  looe1 = get_fkeys_hycom(looe1);
  looe1 = station_reorient_vectors(looe1,'isobath_orientation','fkeys_hycom_u','fkeys_hycom_v');
  looe1 = get_gom_hycom(looe1);
  looe1 = station_reorient_vectors(looe1,'isobath_orientation','gom_hycom_u','gom_hycom_v');
end;

if ( ~exist('rsmas','var') )
  rsmas = read_rsmas_weatherpack_data;
  rsmas.rsmas_dsrf_daily_dose = par_dose(rsmas.rsmas_dsrf);
  rsmas.rsmas_dlrf_daily_dose = par_dose(rsmas.rsmas_dlrf);
end;

if ( ~exist('bhwv3','var') )
  stns = load_haus_waves;
  bhwv3 = stns.bhwv3;
  stns=[]; clear stns;
  bhwv3 = get_erai_station(bhwv3);
  bhwv3 = adjust_erai_station_waves(bhwv3);
  bhwv3 = get_ww3_station(bhwv3);
end;

% Statistical Comparison	Winter (JFM)	Spring (AMJ)	Summer  (JAS)	Autumn (OND)	Overall
% RSMAS vs. ERAI QSWI
% RSMAS vs. NARR QSWI
% MLRF1 vs. adj. QSWI
% RSMAS vs. ERAI QLWI
% RSMAS vs. NARR QLWI
% RSMAS vs. ERAI wh
% RSMAS vs. adjusted wh
% RSMAS vs. WW3 wh
% MLRF1 vs. ERAI Ta
% SMKF1 vs. ERAI qa
% MLRF1 vs. AVHRR Ts
% TSG vs. AVHRR dxTs
% AVHRR dxTs, MLRF1
% AVHRR dyTs, MLRF1
% AVHRR dxTs, LONF1
% AVHRR dyTs, LONF1
% AVHRR vs. GoM dxTs
% AVHRR vs.FKEYS dxTs
% LOOE1 vs. FKEYS u
% LOOE1 vs. FKEYS v
% LOOE1 vs. GoM u
% LOOE1 vs. GoM v

%COMPARE PAR to ERAI and NARR insolation
%% ECMWF Reanalysis-Interim
%adjust_reanalysis=true
%mlrf1_par=phd_ch3_parinsoldose;
%adjust_reanalysis=false
%mlrf1_par=phd_ch3_parinsoldose;
%% North American Regional Reanalysis
%mlrf1_par=phd_ch3_parinsoldose([],'ncep');

if (1)
diary('app_ii.log');
more off;
for cfn={@ts_jfm,@ts_amj,@ts_jas,@ts_ond,[],};
fn=cfn{:};
disp({char(fn)});
app_ii_fit_disp(rsmas.rsmas_dsrf_daily_dose,fwyf1.erai_dsrf_daily_dose,fn,fn,'RSMAS QSWI','ERAI QSWI dose');
app_ii_fit_disp(rsmas.rsmas_dsrf_daily_dose,fwyf1.ncep_dsrf_daily_dose,fn,fn,'RSMAS QSWI','NARR QSWI dose');
app_ii_fit_disp(rsmas.rsmas_dsrf_daily_dose,fwyf1.erai_dsrf_adj_daily_dose,fn,fn,'RSMAS QSWI','ERAI adj QSWI dose');

app_ii_fit_disp(rsmas.rsmas_dlrf_daily_dose,fwyf1.erai_dlrf_daily_dose,fn,fn,'RSMAS QLWI','ERAI QLWI dose');
app_ii_fit_disp(rsmas.rsmas_dlrf_daily_dose,fwyf1.ncep_dlrf_daily_dose,fn,fn,'RSMAS QLWI','NARR QLWI dose');
app_ii_fit_disp(rsmas.rsmas_dlrf_daily_dose,fwyf1.erai_dlrf_adj_daily_dose,fn,fn,'RSMAS QLWI','ERAI adj QLWI dose');

try,
app_ii_fit_disp(bhwv3.triaxys_sigwavehgt,bhwv3.erai_sigwavehgt,fn,fn,'RSMAS wh','ERAI wh');
app_ii_fit_disp(bhwv3.triaxys_sigwavehgt,bhwv3.erai_sigwavehgt_adj,fn,fn,'RSMAS wh','ERAI adj wh');
app_ii_fit_disp(bhwv3.triaxys_sigwavehgt,bhwv3.ww3_sigwavehgt,fn,fn,'RSMAS wh','WW3 wh');
catch,
end;

app_ii_fit_disp(mlrf1.ndbc_air_t,mlrf1.erai_air_t,fn,fn,'MLRF1 Ta','ERAI Ta');
app_ii_fit_disp(smkf1.ndbc_spechumid,smkf1.erai_spechumid,fn,fn,'SMKF1 qa','ERAI qa');
%app_ii_fit_disp(lonf1.ndbc_spechumid,lonf1.erai_spechumid,fn,fn,'LONF1 qa','ERAI qa');
app_ii_fit_disp(mlrf1.ndbc_sea_t,mlrf1.hourly_avhrr_weekly_sst,fn,fn,'MLRF1 Ts','AVHRR Ts');
app_ii_fit_disp(lonf1.ndbc_sea_t,lonf1.hourly_avhrr_weekly_sst,fn,fn,'LONF1 Ts','AVHRR Ts');

%app_ii_fit_disp(tsg.dTdx,mlrf1.hourly_avhrr_weekly_sst,fn,fn,'MLRF1 Ts','AVHRR Ts');
app_ii_fit_disp(mlrf1.hourly_avhrr_weekly_sst_xshore,mlrf1.gom_hycom_seatemp_xshore,fn,fn,'AVHRR dxsTs','GoM dxsTs');
%app_ii_fit_disp(mlrf1.hourly_avhrr_weekly_sst_lshore,mlrf1.gom_hycom_seatemp_lshore,fn,fn,'AVHRR dlsTs','GoM dlsTs');
app_ii_fit_disp(mlrf1.hourly_avhrr_weekly_sst_xshore,mlrf1.fkeys_hycom_seatemp_xshore,fn,fn,'AVHRR dxsTs','FKEYS dxsTs');
%app_ii_fit_disp(mlrf1.hourly_avhrr_weekly_sst_lshore,mlrf1.fkeys_hycom_seatemp_lshore,fn,fn,'AVHRR dlsTs','FKEYS dlsTs');

app_ii_fit_disp(looe1.adcp_x,looe1.gom_hycom_xshore,fn,fn,'LOOE1 uxs','GoM uxs');
app_ii_fit_disp(looe1.adcp_l,looe1.gom_hycom_lshore,fn,fn,'LOOE1 uls','GoM uls');
app_ii_fit_disp(looe1.adcp_x,looe1.fkeys_hycom_xshore,fn,fn,'LOOE1 uxs','FKEYS uxs');
app_ii_fit_disp(looe1.adcp_l,looe1.fkeys_hycom_lshore,fn,fn,'LOOE1 uls','FKEYS uls');

end;

diary off;
more on;
end;


if (0)
diary('app_ii.log');
more off;
%% WHOLE-RECORD
%
%disp({'ERAI QSWI ',nanmedian(mlrf1.erai_dsrf.data),iqr(mlrf1.erai_dsrf.data),});
disp({'ERAI adj QSWI     ',nanmedian(mlrf1.erai_dsrf_adj.data),iqr(mlrf1.erai_dsrf_adj.data),});
disp({'ERAI adj QLWI     ',nanmedian(mlrf1.erai_dlrf_adj.data),iqr(mlrf1.erai_dlrf_adj.data),});
disp({'ERAI adj QSWI dos ',nanmedian(mlrf1.erai_dsrf_adj_daily_dose.data),iqr(mlrf1.erai_dsrf_adj_daily_dose.data),});
disp({'ERAI adj QLWI dos ',nanmedian(mlrf1.erai_dlrf_adj_daily_dose.data),iqr(mlrf1.erai_dlrf_adj_daily_dose.data),});
disp({'ERAI adj wh       ',nanmedian(mlrf1.erai_sigwavehgt_adj.data),iqr(mlrf1.erai_sigwavehgt_adj.data),});
disp({'ERAI Ta           ',nanmedian(mlrf1.erai_air_t.data),iqr(mlrf1.erai_air_t.data),});
disp({'ERAI U10          ',nanmedian(mlrf1.erai_wind_speed.data),iqr(mlrf1.erai_wind_speed.data),});
disp({'SMKF1 ERAI qa     ',nanmedian(smkf1.erai_spechumid.data),iqr(smkf1.erai_spechumid.data),});
disp({'MLRF1 AVHRR Ts    ',nanmedian(mlrf1.hourly_avhrr_weekly_sst.data),iqr(mlrf1.hourly_avhrr_weekly_sst.data),});
disp({'LONF1 AVHRR Ts    ',nanmedian(lonf1.hourly_avhrr_weekly_sst.data),iqr(lonf1.hourly_avhrr_weekly_sst.data),});
disp({'MLRF1 AVHRR dxsTs ',nanmedian(mlrf1.hourly_avhrr_weekly_sst_xshore.data),iqr(mlrf1.hourly_avhrr_weekly_sst_xshore.data),});
disp({'MLRF1 AVHRR dlsTs ',nanmedian(mlrf1.hourly_avhrr_weekly_sst_lshore.data),iqr(mlrf1.hourly_avhrr_weekly_sst_lshore.data),});
disp({'MLRF1 AVHRR D2Ts  ',nanmedian(mlrf1.hourly_avhrr_weekly_sst_l.data),iqr(mlrf1.hourly_avhrr_weekly_sst_l.data),});
disp({'MLRF1 GoM dxsTs   ',nanmedian(mlrf1.gom_hycom_seatemp_xshore.data),iqr(mlrf1.gom_hycom_seatemp_xshore.data),});
disp({'MLRF1 GoM dlsTs   ',nanmedian(mlrf1.gom_hycom_seatemp_lshore.data),iqr(mlrf1.gom_hycom_seatemp_lshore.data),});
disp({'MLRF1 GoM D2Ts    ',nanmedian(mlrf1.gom_hycom_seatemp_l.data),iqr(mlrf1.gom_hycom_seatemp_l.data),});
disp({'MLRF1 FKEYS dxsTs ',nanmedian(mlrf1.fkeys_hycom_seatemp_xshore.data),iqr(mlrf1.fkeys_hycom_seatemp_xshore.data),});
disp({'MLRF1 FKEYS dlsTs ',nanmedian(mlrf1.fkeys_hycom_seatemp_lshore.data),iqr(mlrf1.fkeys_hycom_seatemp_lshore.data),});
disp({'MLRF1 FKEYS D2Ts  ',nanmedian(mlrf1.fkeys_hycom_seatemp_l.data),iqr(mlrf1.fkeys_hycom_seatemp_l.data),});
disp({'LONF1 AVHRR dxsTs ',nanmedian(lonf1.hourly_avhrr_weekly_sst_xshore.data),iqr(lonf1.hourly_avhrr_weekly_sst_xshore.data),});
disp({'LONF1 AVHRR dlsTs ',nanmedian(lonf1.hourly_avhrr_weekly_sst_lshore.data),iqr(lonf1.hourly_avhrr_weekly_sst_lshore.data),});
disp({'LOOE1 ADCP uxs    ',nanmedian(looe1.adcp_x.data),iqr(looe1.adcp_x.data),});
disp({'LOOE1 ADCP uls    ',nanmedian(looe1.adcp_l.data),iqr(looe1.adcp_l.data),});
disp({'LOOE1 GoM uxs     ',nanmedian(looe1.gom_hycom_xshore.data),iqr(looe1.gom_hycom_xshore.data),});
disp({'LOOE1 GoM uls     ',nanmedian(looe1.gom_hycom_lshore.data),iqr(looe1.gom_hycom_lshore.data),});
disp({'LOOE1 FKEYS uxs   ',nanmedian(looe1.fkeys_hycom_xshore.data),iqr(looe1.fkeys_hycom_xshore.data),});
disp({'LOOE1 FKEYS uls   ',nanmedian(looe1.fkeys_hycom_lshore.data),iqr(looe1.fkeys_hycom_lshore.data),});

%% WHOLE-RECORD
%    'ERAI adj QSWI     '    49, 405
%    'ERAI adj QLWI     '    393, 43
%    'ERAI adj wh       '    0.50, 0.45
%    'ERAI Ta           '    25.6, 4.6
%    'ERAI U10          '    8.1, 5.5
%    'SMKF1 ERAI qa     '    0.01, 0.005
%    'MLRF1 AVHRR Ts    '    26.2, 4.1
%    'LONF1 AVHRR Ts    '    25.9, 6.3
%    'MLRF1 AVHRR dxsTs '    3.8x10-5, 1.7x10-4
%    'MLRF1 AVHRR dlsTs '    -1.2x10-5, 0.9x10-4
%    'MLRF1 AVHRR D2Ts  '    1.4x10-8, 2.5x10-7
%    'MLRF1 GoM dxsTs   '    -2.3x10-5, 2.2x10-4
%    'MLRF1 GoM dlsTs   '    -0.8x10-5, 0.4x10-4
%    'MLRF1 GoM D2Ts    '    1.6x10-8, 0.8x10-7
%    'MLRF1 FKEYS dxsTs '    4.3x10-5, 1.5x10-4
%    'MLRF1 FKEYS dlsTs '    -0.1x10-5, 0.2x10-4
%    'MLRF1 FKEYS D2Ts  '    -0.8x10-8, 1.1x10-7
%    'LONF1 AVHRR dxsTs '    -0.07x10-5, 1.4x10-4
%    'LONF1 AVHRR dlsTs '    -1.35x10-5, 1.5x10-4
%    'LOOE1 ADCP uxs    '    0.00, 0.03
%    'LOOE1 ADCP uls    '    -0.01, 0.28
%    'LOOE1 GoM uxs     '    0.00, 0.06
%    'LOOE1 GoM uls     '    -0.04, 0.26
%    'LOOE1 FKEYS uxs   '    -0.01, 0.07
%    'LOOE1 FKEYS uls   '    0.14, 0.47
diary off;
more on;
end;

if (0)
diary('app_ii.log');
more off;
%% WINTER,SPRING,SUMMER,AUTUMN
%
for cfn={@ts_jfm,@ts_amj,@ts_jas,@ts_ond};
fn=cfn{:};
disp(char(fn));
%app_ii_disp('ERAI QSWI ',mlrf1.erai_dsrf,fn);
app_ii_disp('ERAI adj QSWI     ',mlrf1.erai_dsrf_adj,fn);
app_ii_disp('ERAI adj QLWI     ',mlrf1.erai_dlrf_adj,fn);
app_ii_disp('ERAI adj QSWI dos ',mlrf1.erai_dsrf_adj_daily_dose,fn);
app_ii_disp('ERAI adj QLWI dos ',mlrf1.erai_dlrf_adj_daily_dose,fn);
app_ii_disp('ERAI adj wh       ',mlrf1.erai_sigwavehgt_adj,fn);
app_ii_disp('ERAI Ta           ',mlrf1.erai_air_t,fn);
app_ii_disp('ERAI U10          ',mlrf1.erai_wind_speed,fn);
app_ii_disp('SMKF1 ERAI qa     ',smkf1.erai_spechumid,fn);
app_ii_disp('MLRF1 AVHRR Ts    ',mlrf1.hourly_avhrr_weekly_sst,fn);
app_ii_disp('LONF1 AVHRR Ts    ',lonf1.hourly_avhrr_weekly_sst,fn);
app_ii_disp('MLRF1 AVHRR dxsTs ',mlrf1.hourly_avhrr_weekly_sst_xshore,fn);
app_ii_disp('MLRF1 AVHRR dlsTs ',mlrf1.hourly_avhrr_weekly_sst_lshore,fn);
app_ii_disp('MLRF1 AVHRR D2Ts  ',mlrf1.hourly_avhrr_weekly_sst_l,fn);
app_ii_disp('MLRF1 GoM dxsTs   ',mlrf1.gom_hycom_seatemp_xshore,fn);
app_ii_disp('MLRF1 GoM dlsTs   ',mlrf1.gom_hycom_seatemp_lshore,fn);
app_ii_disp('MLRF1 GoM D2Ts    ',mlrf1.gom_hycom_seatemp_l,fn);
app_ii_disp('MLRF1 FKEYS dxsTs ',mlrf1.fkeys_hycom_seatemp_xshore,fn);
app_ii_disp('MLRF1 FKEYS dlsTs ',mlrf1.fkeys_hycom_seatemp_lshore,fn);
app_ii_disp('MLRF1 FKEYS D2Ts  ',mlrf1.fkeys_hycom_seatemp_l,fn);
app_ii_disp('LONF1 AVHRR dxsTs ',lonf1.hourly_avhrr_weekly_sst_xshore,fn);
app_ii_disp('LONF1 AVHRR dlsTs ',lonf1.hourly_avhrr_weekly_sst_lshore,fn);
app_ii_disp('LOOE1 ADCP uxs    ',looe1.adcp_x,fn);
app_ii_disp('LOOE1 ADCP uls    ',looe1.adcp_l,fn);
app_ii_disp('LOOE1 GoM uxs     ',looe1.gom_hycom_xshore,fn);
app_ii_disp('LOOE1 GoM uls     ',looe1.gom_hycom_lshore,fn);
app_ii_disp('LOOE1 FKEYS uxs   ',looe1.fkeys_hycom_xshore,fn);
app_ii_disp('LOOE1 FKEYS uls   ',looe1.fkeys_hycom_lshore,fn);
end;
diary off;
more on;
end;


% TRIMMED
% ts_jfm
%     'ERAI adj QSWI     '    37, 356
%     'ERAI adj QLWI     '    365, 39
%     'ERAI adj wh       '    0.65, 0.46
%     'ERAI Ta           '    22.3, 3.8
%     'ERAI U10          '    9.5, 5.4
%     'SMKF1 ERAI qa     '    0.01, 0.005
%     'MLRF1 AVHRR Ts    '    23.4, 1.3
%     'LONF1 AVHRR Ts    '    22.1, 3.0
%     'MLRF1 AVHRR dxsTs '    9.3x10-5, 1.5x10-4
%     'MLRF1 AVHRR dlsTs '    -1.6x10-5, 0.7x10-4
%     'MLRF1 AVHRR D2Ts  '    -0.2x10-8, 2.2x10-7
%     'MLRF1 GoM dxsTs   '    6.7x10-5, 2.1x10-4
%     'MLRF1 GoM dlsTs   '    -0.7x10-5, 0.4x10-4
%     'MLRF1 GoM D2Ts    '    -0.9x10-8, 0.7x10-7
%     'MLRF1 FKEYS dxsTs '    7.3x10-5, 2.2x10-4
%     'MLRF1 FKEYS dlsTs '    0.0x10-5, 0.3x10-4
%     'MLRF1 FKEYS D2Ts  '    0.1x10-8, 1.5x10-7
%     'LONF1 AVHRR dxsTs '    -0.6x10-5, 1.3x10-4
%     'LONF1 AVHRR dlsTs '    -3.4x10-5, 1.4x10-4
%     'LOOE1 ADCP uxs    '    0.01, 0.03
%     'LOOE1 ADCP uls    '    0.01, 0.28
%     'LOOE1 GoM uxs     '    -0.01, 0.06
%     'LOOE1 GoM uls     '    0.03, 0.36
%     'LOOE1 FKEYS uxs   '    -0.01, 0.08
%     'LOOE1 FKEYS uls   '    0.13, 0.47
%
% ts_amj
%     'ERAI adj QSWI     '    98, 517
%     'ERAI adj QLWI     '    395, 30
%     'ERAI adj wh       '    0.44, 0.40
%     'ERAI Ta           '    26.2, 3.0
%     'ERAI U10          '    7.6, 5.4
%     'SMKF1 ERAI qa     '    0.02, 0.004
%     'MLRF1 AVHRR Ts    '    26.7, 2.6
%     'LONF1 AVHRR Ts    '    27.4, 3.3
%     'MLRF1 AVHRR dxsTs '    -1.2x10-5, 1.2x10-4
%     'MLRF1 AVHRR dlsTs '    -0.8x10-5, 0.8x10-4
%     'MLRF1 AVHRR D2Ts  '    4.2x10-8, 2.1x10-7
%     'MLRF1 GoM dxsTs   '    -10.7x10-5, 1.9x10-4
%     'MLRF1 GoM dlsTs   '    -1.9x10-5, 0.5x10-4
%     'MLRF1 GoM D2Ts    '    4.8x10-8, 0.7x10-7
%     'MLRF1 FKEYS dxsTs '    2.9x10-5, 1.2x10-4
%     'MLRF1 FKEYS dlsTs '    -0.3x10-5, 0.2x10-4
%     'MLRF1 FKEYS D2Ts  '    -0.7x10-8, 1.0x10-7
%     'LONF1 AVHRR dxsTs '    0.5x10-5, 1.4x10-4
%     'LONF1 AVHRR dlsTs '    -0.9x10-5, 1.6x10-4
%     'LOOE1 ADCP uxs    '    0.01, 0.03
%     'LOOE1 ADCP uls    '    0.01, 0.26
%     'LOOE1 GoM uxs     '    0.00, 0.06
%     'LOOE1 GoM uls     '    -0.01, 0.26
%     'LOOE1 FKEYS uxs   '    -0.01, 0.08
%     'LOOE1 FKEYS uls   '    0.14, 0.50
%
% ts_jas
%     'ERAI adj QSWI     '    90, 460
%     'ERAI adj QLWI     '    416, 11
%     'ERAI adj wh       '    0.34, 0.27
%     'ERAI Ta           '    27.9, 2.1
%     'ERAI U10          '    6.3, 4.0
%     'SMKF1 ERAI qa     '    0.02, 0.001
%     'MLRF1 AVHRR Ts    '    29.2, 1.0
%     'LONF1 AVHRR Ts    '    29.8, 1.4
%     'MLRF1 AVHRR dxsTs '    -2.1x10-5, 1.2x10-4
%     'MLRF1 AVHRR dlsTs '    -0.9x10-5, 1.1x10-4
%     'MLRF1 AVHRR D2Ts  '    4.1x10-8, 2.6x10-7
%     'MLRF1 GoM dxsTs   '    -11.3x10-5, 1.3x10-4
%     'MLRF1 GoM dlsTs   '    0.3x10-5, 0.4x10-4
%     'MLRF1 GoM D2Ts    '    5.3x10-8, 0.5x10-7
%     'MLRF1 FKEYS dxsTs '    1.4x10-5, 1.0x10-4
%     'MLRF1 FKEYS dlsTs '    -0.1x10-5, 0.1x10-4
%     'MLRF1 FKEYS D2Ts  '    -1.6x10-8, 0.9x10-7
%     'LONF1 AVHRR dxsTs '    -0.7x10-5, 1.6x10-4
%     'LONF1 AVHRR dlsTs '    -0.6x10-5, 1.5x10-4
%     'LOOE1 ADCP uxs    '    0.01, 0.03
%     'LOOE1 ADCP uls    '    -0.01, 0.28
%     'LOOE1 GoM uxs     '    0.02, 0.06
%     'LOOE1 GoM uls     '    -0.05, 0.20
%     'LOOE1 FKEYS uxs   '    -0.01, 0.07
%     'LOOE1 FKEYS uls   '    0.21, 0.48
%
% ts_ond
%     'ERAI adj QSWI     '    15, 312
%     'ERAI adj QLWI     '    384, 36
%     'ERAI adj wh       '    0.66, 0.48
%     'ERAI Ta           '    24.4, 3.8
%     'ERAI U10          '    9.3, 5.5
%     'SMKF1 ERAI qa     '    0.01, 0.005
%     'MLRF1 AVHRR Ts    '    25.7, 2.3
%     'LONF1 AVHRR Ts    '    24.0, 3.7
%     'MLRF1 AVHRR dxsTs '    11.9x10-5, 1.8x10-4
%     'MLRF1 AVHRR dlsTs '    -1.6x10-5, 1.0x10-4
%     'MLRF1 AVHRR D2Ts  '    -2.5x10-8, 2.8x10-7
%     'MLRF1 GoM dxsTs   '    7.7x10-5, 2.0x10-4
%     'MLRF1 GoM dlsTs   '    -0.9x10-5, 0.4x10-4
%     'MLRF1 GoM D2Ts    '    -2.1x10-8, 0.6x10-7
%     'MLRF1 FKEYS dxsTs '    8.3x10-5, 2.1x10-4
%     'MLRF1 FKEYS dlsTs '    -0.1x10-5, 0.2x10-4
%     'MLRF1 FKEYS D2Ts  '    -0.8x10-8, 1.2x10-7
%     'LONF1 AVHRR dxsTs '    0.1x10-5, 1.4x10-4
%     'LONF1 AVHRR dlsTs '    -0.7x10-5, 1.3x10-4
%     'LOOE1 ADCP uxs    '    0.00, 0.03
%     'LOOE1 ADCP uls    '    -0.05, 0.30
%     'LOOE1 GoM uxs     '    0.00, 0.06
%     'LOOE1 GoM uls     '    -0.08, 0.21
%     'LOOE1 FKEYS uxs   '    -0.01, 0.06
%     'LOOE1 FKEYS uls   '    0.08, 0.42
