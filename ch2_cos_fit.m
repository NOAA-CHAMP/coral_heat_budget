1;
%% SCRIPT to fit single cosine curves to SEAKEYS station annual cycles
%
% Last Saved Time-stamp: <Mon 2013-11-04 15:39:54 Eastern Standard Time gramer>

%fld = 'ndbc_wind1_speed';
%fld = 'ndbc_wind1_xshore';
%fld = 'ndbc_air_t';
fld = 'ndbc_sea_t';

%if ~exist('lkwf1'); lkwf1 = get_station_from_station_name('lkwf1'); lkwf1 = load_all_ndbc_data(lkwf1); end;
if ~exist('fwyf1'); fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1); end;
if ~exist('mlrf1'); mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1); end;
if ~exist('lonf1'); lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1); end;
if ~exist('smkf1'); smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1); end;
if ~exist('sanf1'); sanf1 = get_station_from_station_name('sanf1'); sanf1 = load_all_ndbc_data(sanf1); end;
if ~exist('dryf1'); dryf1 = get_station_from_station_name('dryf1'); dryf1 = load_all_ndbc_data(dryf1); end;

if (0)
%lkwf1 = station_filter_bad_dates(lkwf1);
fwyf1 = station_filter_bad_dates(fwyf1);
mlrf1 = station_filter_bad_dates(mlrf1);
lonf1 = station_filter_bad_dates(lonf1);
smkf1 = station_filter_bad_dates(smkf1);
sanf1 = station_filter_bad_dates(sanf1);
dryf1 = station_filter_bad_dates(dryf1);
end;

% % SUBSET_TS calls could be used to balance out annual cycle
%                                subset_ts(smkf1.(fld),@(x)(x.date<datenum(2005,8,1))),...
%                                subset_ts(dryf1.(fld),@(x)(get_year(x.date)~=2004)) ...

[f,m,l,k,s] = intersect_tses(fwyf1.(fld),...
                             mlrf1.(fld),...
                             lonf1.(fld),...
                             smkf1.(fld),...
                             sanf1.(fld)...
                             );

[ig,d] = intersect_tses(f,dryf1.(fld)); clear ig


if(0)
% Check seasonal cycle and interannual bias
%fmg; hist(get_week(m.date),52); xlim([1,52]); titlename('Data Density by Week');
fmg; hist(get_month(m.date),12); xlim([1,12]); titlename('Data Density by Month');
fmg; hist(get_season(m.date),4); xlim([1, 4]); titlename('Data Density by Season');
fmg; hist(get_year(m.date),numel(unique(get_year(m.date)))); titlename('Data Density by Year');
end;

if(1)
warning('OFF','curvefit:fit:noStartPoint');
%scatter_curve_fit(get_yearday_no_leap(w.date),w.data,fittype('(a1/2)*cos(2*pi*(x-c1)/365)+d1'),'Year-day','LKWF1 T');
%axis([0,366,5,35]);
scatter_curve_fit(get_yearday_no_leap(f.date),f.data,fittype('(a1/2)*cos(2*pi*(x-c1)/365)+d1'),'Year-day','FWYF1 T');
axis([0,366,5,35]);
scatter_curve_fit(get_yearday_no_leap(m.date),m.data,fittype('(a1/2)*cos(2*pi*(x-c1)/365)+d1'),'Year-day','MLRF1 T');
axis([0,366,5,35]);
scatter_curve_fit(get_yearday_no_leap(l.date),l.data,fittype('(a1/2)*cos(2*pi*(x-c1)/365)+d1'),'Year-day','LONF1 T');
axis([0,366,5,35]);
scatter_curve_fit(get_yearday_no_leap(k.date),k.data,fittype('(a1/2)*cos(2*pi*(x-c1)/365)+d1'),'Year-day','SMKF1 T');
axis([0,366,5,35]);
scatter_curve_fit(get_yearday_no_leap(s.date),s.data,fittype('(a1/2)*cos(2*pi*(x-c1)/365)+d1'),'Year-day','SANF1 T');
axis([0,366,5,35]);
scatter_curve_fit(get_yearday_no_leap(d.date),d.data,fittype('(a1/2)*cos(2*pi*(x-c1)/365)+d1'),'Year-day','DRYF1 T');
axis([0,366,5,35]);
warning('ON','curvefit:fit:noStartPoint');
end;
