1;
%% SCRIPT to analyze annual and seasonal means for SEAKEYS stations
%
% Last Saved Time-stamp: <Wed 2014-01-08 17:13:39 Eastern Standard Time gramer>


if ( ~exist('fwyf1','var') )
  fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1); %fwyf1 = load_station_data(fwyf1);
end;
if ( ~exist('mlrf1','var') )
  mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1); %mlrf1 = load_station_data(mlrf1);
end;
if ( ~exist('lonf1','var') )
  lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1); %lonf1 = load_station_data(lonf1);
end;
if ( ~exist('smkf1','var') )
  smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1); smkf1 = load_station_data(smkf1);
end;
if ( ~exist('sanf1','var') )
  sanf1 = get_station_from_station_name('sanf1'); sanf1 = load_all_ndbc_data(sanf1); %sanf1 = load_station_data(sanf1);
end;
if ( ~exist('dryf1','var') )
  dryf1 = get_station_from_station_name('dryf1'); dryf1 = load_all_ndbc_data(dryf1); %dryf1 = load_station_data(dryf1);
end;

if 0
fwyf1 = station_filter_bad_dates(fwyf1);
mlrf1 = station_filter_bad_dates(mlrf1);
lonf1 = station_filter_bad_dates(lonf1);
smkf1 = station_filter_bad_dates(smkf1);
sanf1 = station_filter_bad_dates(sanf1);
dryf1 = station_filter_bad_dates(dryf1);
end;


% Years with at least 5600 data points or 233 full days of data
%for yr=unique(get_year(f.date))'; disp([yr,numel(find(get_year(f.date)==yr))]); end;
yrs = [1993,1995,1998,1999,2000,2001,2003];


if 0
[f,m,l,k,s,d] = intersect_tses([],subset_ts(fwyf1.ndbc_sea_t,@ts_jfm),mlrf1.ndbc_sea_t,lonf1.ndbc_sea_t,smkf1.ndbc_sea_t,sanf1.ndbc_sea_t,dryf1.ndbc_sea_t);
[P,ANOVATAB,STATS] = anova1([f.data,m.data,l.data,k.data,s.data,d.data],{'FWY','MLR','LON','SMK','SAN','DRY'},'off');
%appendtitlename(' T_s');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename(' for intersecting T_s JFM');

[f,m,l,k,s,d] = intersect_tses([],subset_ts(fwyf1.ndbc_sea_t,@ts_amj),mlrf1.ndbc_sea_t,lonf1.ndbc_sea_t,smkf1.ndbc_sea_t,sanf1.ndbc_sea_t,dryf1.ndbc_sea_t);
[P,ANOVATAB,STATS] = anova1([f.data,m.data,l.data,k.data,s.data,d.data],{'FWY','MLR','LON','SMK','SAN','DRY'},'off');
%appendtitlename(' T_s');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename(' for intersecting T_s AMJ');

[f,m,l,k,s,d] = intersect_tses([],subset_ts(fwyf1.ndbc_sea_t,@ts_jas),mlrf1.ndbc_sea_t,lonf1.ndbc_sea_t,smkf1.ndbc_sea_t,sanf1.ndbc_sea_t,dryf1.ndbc_sea_t);
[P,ANOVATAB,STATS] = anova1([f.data,m.data,l.data,k.data,s.data,d.data],{'FWY','MLR','LON','SMK','SAN','DRY'},'off');
%appendtitlename(' T_s');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename(' for intersecting T_s JAS');

[f,m,l,k,s,d] = intersect_tses([],subset_ts(fwyf1.ndbc_sea_t,@ts_ond),mlrf1.ndbc_sea_t,lonf1.ndbc_sea_t,smkf1.ndbc_sea_t,sanf1.ndbc_sea_t,dryf1.ndbc_sea_t);
[P,ANOVATAB,STATS] = anova1([f.data,m.data,l.data,k.data,s.data,d.data],{'FWY','MLR','LON','SMK','SAN','DRY'},'off');
%appendtitlename(' T_s');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename(' for intersecting T_s OND');
end;


if 1

[f,m,l,k,s,d] = intersect_tses([],fwyf1.ndbc_sea_t,mlrf1.ndbc_sea_t,lonf1.ndbc_sea_t,smkf1.ndbc_sea_t,sanf1.ndbc_sea_t,dryf1.ndbc_sea_t);
[P,ANOVATAB,STATS] = anova1([f.data,m.data,l.data,k.data,s.data,d.data],{'FWY','MLR','LON','SMK','SAN','DRY'},'off');
%appendtitlename(' T_s');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename(' for intersecting T_s');

alldat = [fwyf1.ndbc_sea_t.data;...
          mlrf1.ndbc_sea_t.data;...
          lonf1.ndbc_sea_t.data; ...
          smkf1.ndbc_sea_t.data;...
          sanf1.ndbc_sea_t.data;...
          dryf1.ndbc_sea_t.data];
allgps = [repmat('FWY',[numel(fwyf1.ndbc_sea_t.data),1]);...
          repmat('MLR',[numel(mlrf1.ndbc_sea_t.data),1]);...
          repmat('LON',[numel(lonf1.ndbc_sea_t.data),1]);...
          repmat('SMK',[numel(smkf1.ndbc_sea_t.data),1]);...
          repmat('SAN',[numel(sanf1.ndbc_sea_t.data),1]);...
          repmat('DRY',[numel(dryf1.ndbc_sea_t.data),1])];

[P,ANOVATAB,STATS] = anova1(alldat,allgps,'off');
%appendtitlename(' T_s');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak');
appendtitlename(' for all T_s');

end;


if 0

%%%%%%%%%%%%%%%%%%%%
%% Air Temperature


[f,m,l,k,s,d] = intersect_tses([],fwyf1.ndbc_air_t,mlrf1.ndbc_air_t,lonf1.ndbc_air_t,smkf1.ndbc_air_t,sanf1.ndbc_air_t,dryf1.ndbc_air_t);
[P,ANOVATAB,STATS] = anova1([f.data,m.data,l.data,k.data,s.data,d.data],{'FWY','MLR','LON','SMK','SAN','DRY'},'off');
%appendtitlename(' T_a');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' for intersecting T_a, N=',num2str(numel(f.data))]);

alldat = [fwyf1.ndbc_air_t.data;...
          mlrf1.ndbc_air_t.data;...
          lonf1.ndbc_air_t.data; ...
          smkf1.ndbc_air_t.data;...
          sanf1.ndbc_air_t.data;...
          dryf1.ndbc_air_t.data];
allgps = [repmat('FWY',[numel(fwyf1.ndbc_air_t.data),1]);...
          repmat('MLR',[numel(mlrf1.ndbc_air_t.data),1]);...
          repmat('LON',[numel(lonf1.ndbc_air_t.data),1]);...
          repmat('SMK',[numel(smkf1.ndbc_air_t.data),1]);...
          repmat('SAN',[numel(sanf1.ndbc_air_t.data),1]);...
          repmat('DRY',[numel(dryf1.ndbc_air_t.data),1])];

[P,ANOVATAB,STATS] = anova1(alldat,allgps,'off');
%appendtitlename(' T_a');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak');
appendtitlename(' for all T_a');



%%%%%%%%%%%%%%%%%%%%
%% Wind Speed (10m)

wfld = ['ndbc_wind1_vel'];
u10fld = [wfld,'_10m'];
if ~isfield(fwyf1,u10fld); fwyf1.(wfld) = ts_fun(fwyf1.ndbc_wind1_speed,@kts2mps); fwyf1 = station_wind_at_height(fwyf1,wfld,'ndbc_wind1_dir','ndbc_air_t',[]); end;
if ~isfield(mlrf1,u10fld); mlrf1.(wfld) = ts_fun(mlrf1.ndbc_wind1_speed,@kts2mps); mlrf1 = station_wind_at_height(mlrf1,wfld,'ndbc_wind1_dir','ndbc_air_t',[]); end;
if ~isfield(lonf1,u10fld); lonf1.(wfld) = ts_fun(lonf1.ndbc_wind1_speed,@kts2mps); lonf1 = station_wind_at_height(lonf1,wfld,'ndbc_wind1_dir','ndbc_air_t',[]); end;
if ~isfield(smkf1,u10fld); smkf1.(wfld) = ts_fun(smkf1.ndbc_wind1_speed,@kts2mps); smkf1 = station_wind_at_height(smkf1,wfld,'ndbc_wind1_dir','ndbc_air_t',[]); end;
if ~isfield(sanf1,u10fld); sanf1.(wfld) = ts_fun(sanf1.ndbc_wind1_speed,@kts2mps); sanf1 = station_wind_at_height(sanf1,wfld,'ndbc_wind1_dir','ndbc_air_t',[]); end;
if ~isfield(dryf1,u10fld); dryf1.(wfld) = ts_fun(dryf1.ndbc_wind1_speed,@kts2mps); dryf1 = station_wind_at_height(dryf1,wfld,'ndbc_wind1_dir','ndbc_air_t',[]); end;


[f,m,l,k,s,d] = intersect_tses([],fwyf1.(u10fld),mlrf1.(u10fld),lonf1.(u10fld),smkf1.(u10fld),sanf1.(u10fld),dryf1.(u10fld));
[P,ANOVATAB,STATS] = anova1([f.data,m.data,l.data,k.data,s.data,d.data],{'FWY','MLR','LON','SMK','SAN','DRY'},'off');
%appendtitlename(' T_a');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' for intersecting U10, N=',num2str(numel(f.data))]);

alldat = [fwyf1.(u10fld).data;...
          mlrf1.(u10fld).data;...
          lonf1.(u10fld).data; ...
          smkf1.(u10fld).data;...
          sanf1.(u10fld).data;...
          dryf1.(u10fld).data];
allgps = [repmat('FWY',[numel(fwyf1.(u10fld).data),1]);...
          repmat('MLR',[numel(mlrf1.(u10fld).data),1]);...
          repmat('LON',[numel(lonf1.(u10fld).data),1]);...
          repmat('SMK',[numel(smkf1.(u10fld).data),1]);...
          repmat('SAN',[numel(sanf1.(u10fld).data),1]);...
          repmat('DRY',[numel(dryf1.(u10fld).data),1])];

[P,ANOVATAB,STATS] = anova1(alldat,allgps,'off');
%appendtitlename(' U10');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak');
appendtitlename(' for all U10');


% Wind Vector Components (10m)

ufld = [u10fld,'_u'];
vfld = [u10fld,'_v'];

[f,m,l,k,s,d] = intersect_tses([],fwyf1.(ufld),mlrf1.(ufld),lonf1.(ufld),smkf1.(ufld),sanf1.(ufld),dryf1.(ufld));
[P,ANOVATAB,STATS] = anova1([f.data,m.data,l.data,k.data,s.data,d.data],{'FWY','MLR','LON','SMK','SAN','DRY'},'off');
%appendtitlename(' T_a');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' for intersecting u at 10 m, N=',num2str(numel(f.data))]);

[f,m,l,k,s,d] = intersect_tses([],fwyf1.(vfld),mlrf1.(vfld),lonf1.(vfld),smkf1.(vfld),sanf1.(vfld),dryf1.(vfld));
[P,ANOVATAB,STATS] = anova1([f.data,m.data,l.data,k.data,s.data,d.data],{'FWY','MLR','LON','SMK','SAN','DRY'},'off');
%appendtitlename(' T_a');
fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' for intersecting v at 10 m, N=',num2str(numel(f.data))]);

end;

