1;
%%%% SCRIPT to create spectral analysis figures for JMR MS. from Ch. 2

if ~exist('fwyf1','var'); fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1); end;
if ~exist('lonf1','var'); lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1); end;
if ~exist('mlrf1','var'); mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1); end;
if ~exist('smkf1','var'); smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1); end;
if ~exist('sanf1','var'); sanf1 = get_station_from_station_name('sanf1'); sanf1 = load_all_ndbc_data(sanf1); end;
if ~exist('dryf1','var'); dryf1 = get_station_from_station_name('dryf1'); dryf1 = load_all_ndbc_data(dryf1); end;
if ~exist('looe1','var'); looe1 = get_station_from_station_name('looe1'); looe1 = get_looe1_adcp(looe1); end;
if ~isfield(smkf1,'ndbc_spechumid')
  smkf1 = station_dewp_to_relhumid(smkf1,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
  smkf1 = station_relhumid_to_spechumid(smkf1,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
end;
if ~isfield(smkf1,'ndbc_wind1_speed_10m')
  fwyf1 = station_wind_at_height(fwyf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  lonf1 = station_wind_at_height(lonf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  mlrf1 = station_wind_at_height(mlrf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  smkf1 = station_wind_at_height(smkf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  sanf1 = station_wind_at_height(sanf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  dryf1 = station_wind_at_height(dryf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
end;
if ~isfield(smkf1,'ndbc_tide_m')
  lonf1.ndbc_tide_m = lonf1.ndbc_tide; lonf1.ndbc_tide_m.data = unitsratio('m','ft')*(lonf1.ndbc_tide.data-nanmean(lonf1.ndbc_tide.data));
  smkf1.ndbc_tide_m = smkf1.ndbc_tide; smkf1.ndbc_tide_m.data = unitsratio('m','ft')*(smkf1.ndbc_tide.data-nanmean(smkf1.ndbc_tide.data));
end;
if ~isfield(looe1,'adcp_x_contig')
  looe1.adcp_x_contig = subset_ts(looe1.adcp_x,@(x)(find(datenum(2005,3,25)<=x.date&x.date<=datenum(2009,1,25))));
  looe1.adcp_x_contig = subset_ts(looe1.adcp_x_contig,@ts_isfinite);
  looe1.adcp_l_contig = subset_ts(looe1.adcp_l,@(x)(find(datenum(2005,3,25)<=x.date&x.date<=datenum(2009,1,25))));
  looe1.adcp_l_contig = subset_ts(looe1.adcp_l_contig,@ts_isfinite);
end;


%% MAKE FIGURES

CI = 0.85;

doPrint = true;
%doPrint = false;

%CPDbands = 20;
CPDbands = 5;
CPHbands = 5;

if 0
%ch2_spec(fwyf1.ndbc_wind1_speed_10m,CPDbands,true,CI,'FWYF1 U10'); axis([3/24,3000,0,220]);
ch2_spec(fwyf1.ndbc_wind1_speed_10m,CPDbands,true,CI,'FWYF1 U10'); axis([3/24,3000,0,80]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'fwyf1-pmtm-U10.tif')); end;
ch2_spec(mlrf1.ndbc_air_t,CPDbands,true,CI,'MLRF1 T_a'); axis([3/24,3000,0,150]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'mlrf1-pmtm-Ta.tif')); end;
ch2_spec(smkf1.ndbc_spechumid,CPDbands,true,CI,'SMKF1 q_a'); axis([3/24,3000,0,1e-4]);
text(3.5/24,0.95e-4,'\times10^-^4','FontSize',16);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-qa.tif')); end;
ch2_spec(smkf1.ndbc_barom,CPDbands,true,CI,'SMKF1 P_a'); axis([3/24,3000,0,400]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-Pa.tif')); end;
end;

if 0
ch2_spec(smkf1.ndbc_tide_m,CPDbands,true,CI,'SMKF1 h_t_i_d_e'); axis([3/24,3000,0,27]);
set(gca,'XTickLabel',[' 1  ';' 10 ';' 100';'1000']);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-h.tif')); end;
end;

if 1
ch2_spec(smkf1.ndbc_tide_m,CPDbands,true,CI,'SMKF1 h_t_i_d_e'); axis([3/24,30,0,0.3]);
set(gca,'XTickLabel',[' 1  ';' 10 ']);
set(gca,'YTick',[0.0,0.1,0.2,0.3]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-h-INSET.tif')); end;
end;

if 0
ch2_spec(lonf1.ndbc_sea_t,5,true,CI,'LONF1 T_s'); axis([3/24,3000,0,220]);
set(gca,'XTickLabel',[' 1  ';' 10 ';' 100';'1000']);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'lonf1-pmtm-Ts.tif')); end;

ch2_spec(smkf1.ndbc_sea_t,5,true,CI,'SMKF1 T_s'); axis([3/24,3000,0,220]);
set(gca,'XTickLabel',[' 1  ';' 10 ';' 100';'1000']);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-Ts.tif')); end;

ch2_spec(smkf1.ndbc_sea_t,CPDbands,false,CI,'SMKF1 T_s'); axis([23,31,0,1.8]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-Ts-diurnal.tif')); end;
end;

if 0
ch2_spec(smkf1.ndbc_sea_t,CPDbands,false,CI,'SMKF1 T_s'); axis([26,30,0,0.11]);
set(gcf,'OuterPos',[.2,0,.6,1]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-Ts-inertial-INSET.tif')); end;
end;

if 0
%ch2_spec(looe1.adcp_x_contig,5,true,CI,'LOOE1 u_x_s'); axis([3/24,3000,0,0.06]); % 0.023]);
ch2_spec(looe1.adcp_x_contig,5,true,CI,'LOOE1 u_x_s'); axis([3/24,3000,0,0.035]);
set(gca,'XTickLabel',[' 1  ';' 10 ';' 100';'1000']);
set(gca,'YTick',[0.00,0.01,0.02,0.03]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'looe1-pmtm-xs.tif')); end;
end;

if 0
ch2_spec(looe1.adcp_x_contig,CPDbands,false,CI,'LOOE1 u_x_s'); axis([10,40,0,0.005]);
set(gcf,'OuterPos',[.2,0,.6,1]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'looe1-pmtm-xs-inertial-INSET.tif')); end;
end;

if 0
%ch2_spec(looe1.adcp_l_contig,5,true,CI,'LOOE1 u_l_s'); axis([3/24,3000,0,0.06]);
ch2_spec(looe1.adcp_l_contig,5,true,CI,'LOOE1 u_l_s'); axis([3/24,400,0,0.50]);
set(gca,'XTickLabel',[' 1  ';' 10 ';' 100';'1000']);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'looe1-pmtm-ls.tif')); end;
end;

if 0
ch2_spec(looe1.adcp_l_contig,CPDbands,false,CI,'LOOE1 u_l_s'); axis([10,40,0,0.05]);
set(gcf,'OuterPos',[.2,0,.6,1]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'looe1-pmtm-ls-inertial-INSET.tif')); end;
end;
