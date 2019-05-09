1;

if ~exist('smkf1'); smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1); end;
if ~exist('lonf1'); lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1); end;

if (0)
smkf1 = station_filter_bad_dates(smkf1);
lonf1 = station_filter_bad_dates(lonf1);
end;

stormdates = datenum(2005,[6,8,9,9,10,10,10],[2,26,20,21,24,25,26]);
fn = @(x)(~ismember(floor(x.date),stormdates));
smkf1.ndbc_tide = subset_ts(smkf1.ndbc_tide,fn);
smkf1.ndbc_sea_t = subset_ts(smkf1.ndbc_sea_t,fn);
lonf1.ndbc_tide = subset_ts(lonf1.ndbc_tide,fn);
lonf1.ndbc_sea_t = subset_ts(lonf1.ndbc_sea_t,fn);
clear fn;

disp('Tide reporting accuracy [mm]');
disp(min(diff(unique(smkf1.ndbc_tide.data)))*unitsratio('m','ft')*1e3);

% Filter incomplete years (<333 d worth)
minhrs = (333*24);
%minhrs = (300*24);
%minhrs = (270*24);

warning('OFF','run_length_encode:generic');

fn = @ts_isfinite; yroff=-213;
%fn = @ts_isfinite; yroff=0;
%fn = @ts_boreal_warm; minhrs = (150*24); yroff=0;
%fn = @ts_boreal_cool; minhrs = (150*24); yroff=0;

% SMKF1 TIDE data starts 2000-08-01, ends 2013-08-31, thus -213 d
ts = smkf1.ndbc_tide;
ts = subset_ts(ts,fn);
yrs = get_year(ts.date+yroff); h = ts.data.*unitsratio('m','ft'); h = h-nanmean(h);
[nyrs,uyrs] = run_length_encode(yrs); badix = find(ismember(yrs,uyrs(nyrs < minhrs))); yrs(badix) = []; h(badix) = [];
[P,ANOVATAB,STATS] = anova1(h,yrs,'off'); fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' SMKF1 Tide ',num2str(P)]); view(-90,90);
scatter_fit(str2num(char(STATS.gnames)),STATS.means','Year','SMKF1 Tide')

[B,BINT]=regress(h,[repmat(1,size(yrs)),yrs],.05);
{'SMKF1 Tide [mm] ',numel(smkf1.ndbc_tide.data),': ',B(2).*1000,' +/- ',(B(2)-BINT(2,1)).*1000},


% LONF1 TIDE data starts 2000-08-01, ends 2008-01-19
ts = lonf1.ndbc_tide;
ts = subset_ts(ts,fn);
yrs = get_year(ts.date); h = ts.data.*unitsratio('m','ft'); h = h-nanmean(h);
[nyrs,uyrs] = run_length_encode(yrs); badix = find(ismember(yrs,uyrs(nyrs < minhrs))); yrs(badix) = []; h(badix) = [];
[P,ANOVATAB,STATS] = anova1(h,yrs,'off'); fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' LONF1 Tide ',num2str(P)]); view(-90,90);
scatter_fit(str2num(char(STATS.gnames)),STATS.means','Year','LONF1 Tide')

[B,BINT]=regress(h,[repmat(1,size(yrs)),yrs],.05);
{'LONF1 Tide [mm] ',numel(smkf1.ndbc_tide.data),': ',B(2).*1000,' +/- ',(B(2)-BINT(2,1)).*1000},


% % SMKF1 SEA_T data starts ?
% ts = smkf1.ndbc_sea_t;
% ts = subset_ts(ts,fn);
% yrs = get_year(ts.date+yroff); h = ts.data;
% [nyrs,uyrs] = run_length_encode(yrs); badix = find(ismember(yrs,uyrs(nyrs < minhrs))); yrs(badix) = []; h(badix) = [];
% [P,ANOVATAB,STATS] = anova1(h,yrs,'off'); fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' ',upper(smkf1.station_name),' Sea T']); view(-90,90);
% scatter_fit(str2num(char(STATS.gnames)),STATS.means','Year','SMKF1 Sea T')

warning('ON','run_length_encode:generic');



if (1)

%fn = @ts_isfinite;
fn = @(x)(find(x.date<datenum(2008,1,18)));
%fn = @(x)(find(x.date<datenum(2009,1,1)));
%fn = @(x)(find(datenum(2001,1,1)<=x.date & x.date<datenum(2010,1,1)));
disp(['Hourly fit with ',upper(char(fn))]);

ts = smkf1.ndbc_tide;
ix = fn(ts);
yrs = yearfrac(ts.date(ix)); h = ts.data(ix).*unitsratio('m','ft'); h = h-nanmean(h);
[B,Stats,fh,lh]=scatter_fit(yrs,h,'Year','Tide [m]'); titlename(['SMKF1 Tide Height ',char(fn)]);
set(lh(1),'Color',[.5,.5,.5],'LineStyle',':');
set(lh(2),'Color','k','LineWidth',2);
ylim([-1,1]);
%legend hide;
%axis([2001,2010,-0.6,1.6]);

ts = lonf1.ndbc_tide;
ix = fn(ts);
yrs = yearfrac(ts.date(ix)); h = ts.data(ix).*unitsratio('m','ft'); h = h-nanmean(h);
[B,Stats,fh,lh]=scatter_fit(yrs,h,'Year','Tide [m]'); titlename(['LONF1 Tide Height ',char(fn)]);
set(lh(1),'Color',[.5,.5,.5],'LineStyle',':');
set(lh(2),'Color','k','LineWidth',2);
ylim([-1,1]);
%legend hide;
%axis([2001,2010,-0.6,1.6]);

end;





if (0)

[P,ANOVATAB,STATS] = anova1(smkf1.ndbc_tide.data,get_year(smkf1.ndbc_tide.date),'off'); fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' ',upper(smkf1.station_name),' Tide']); view(-90,90);
[P,ANOVATAB,STATS] = anova1(lonf1.ndbc_tide.data,get_year(lonf1.ndbc_tide.date),'off'); fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' ',upper(lonf1.station_name),' Tide']); view(-90,90);

fn = @(x)(ts_jas(x));
disp(fn);
ix = fn(smkf1.ndbc_tide);
[P,ANOVATAB,STATS] = anova1(smkf1.ndbc_tide.data(ix),get_year(smkf1.ndbc_tide.date(ix)),'off'); fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' ',upper(smkf1.station_name),' Tide ',char(fn)]); view(-90,90);
ix = fn(lonf1.ndbc_tide);
[P,ANOVATAB,STATS] = anova1(lonf1.ndbc_tide.data(ix),get_year(lonf1.ndbc_tide.date(ix)),'off'); fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' ',upper(lonf1.station_name),' Tide ',char(fn)]); view(-90,90);

fn = @(x)(ts_ond(x));
disp(fn);
ix = fn(smkf1.ndbc_tide);
[P,ANOVATAB,STATS] = anova1(smkf1.ndbc_tide.data(ix),get_year(smkf1.ndbc_tide.date(ix)),'off'); fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' ',upper(smkf1.station_name),' Tide ',char(fn)]); view(-90,90);
ix = fn(lonf1.ndbc_tide);
[P,ANOVATAB,STATS] = anova1(lonf1.ndbc_tide.data(ix),get_year(lonf1.ndbc_tide.date(ix)),'off'); fmg; multcompare(STATS,'ctype','tukey-kramer dunn-sidak'); appendtitlename([' ',upper(lonf1.station_name),' Tide ',char(fn)]); view(-90,90);

end;


if (0)

[B,BINT]=regress(smkf1.ndbc_tide.data,[repmat(1,size(smkf1.ndbc_tide.date)),yearfrac(smkf1.ndbc_tide.date)],.05);
{'SMKF1 Tide ',numel(smkf1.ndbc_tide.data),': ',B(2).*100,' +/- ',(B(2)-BINT(2,1)).*100},
[B,BINT]=regress(lonf1.ndbc_tide.data,[repmat(1,size(lonf1.ndbc_tide.date)),yearfrac(lonf1.ndbc_tide.date)],.05);
{'LONF1 Tide ',numel(lonf1.ndbc_tide.data),': ',B(2).*100,' +/- ',(B(2)-BINT(2,1)).*100},
[B,BINT]=regress(smkf1.ndbc_sea_t.data,[repmat(1,size(smkf1.ndbc_sea_t.date)),yearfrac(smkf1.ndbc_sea_t.date)],.05);
{'SMKF1 Ts ',numel(smkf1.ndbc_sea_t.data),': ',B(2),' +/- ',(B(2)-BINT(2,1))},
[B,BINT]=regress(lonf1.ndbc_sea_t.data,[repmat(1,size(lonf1.ndbc_sea_t.date)),yearfrac(lonf1.ndbc_sea_t.date)],.05);
{'LONF1 Ts ',numel(lonf1.ndbc_sea_t.data),': ',B(2),' +/- ',(B(2)-BINT(2,1))},


%fn = @(x)(find(x.date<datenum(2009,1,1)));
fn = @(x)(find(datenum(2001,1,1)<=x.date & x.date<datenum(2010,1,1)));
disp(fn);

ix = fn(smkf1.ndbc_tide);
[B,BINT]=regress(smkf1.ndbc_tide.data(ix),[repmat(1,size(smkf1.ndbc_tide.date(ix))),yearfrac(smkf1.ndbc_tide.date(ix))],.05);
{'SMKF1 Tide ',numel(ix),': ',B(2).*100,' +/- ',(B(2)-BINT(2,1)).*100},

ix = fn(lonf1.ndbc_tide);
[B,BINT]=regress(lonf1.ndbc_tide.data(ix),[repmat(1,size(lonf1.ndbc_tide.date(ix))),yearfrac(lonf1.ndbc_tide.date(ix))],.05);
{'LONF1 Tide ',numel(ix),': ',B(2).*100,' +/- ',(B(2)-BINT(2,1)).*100},

ix = fn(smkf1.ndbc_sea_t);
[B,BINT]=regress(smkf1.ndbc_sea_t.data(ix),[repmat(1,size(smkf1.ndbc_sea_t.date(ix))),yearfrac(smkf1.ndbc_sea_t.date(ix))],.05);
{'SMKF1 Ts ',numel(ix),': ',B(2),' +/- ',(B(2)-BINT(2,1))},

ix = fn(lonf1.ndbc_sea_t);
[B,BINT]=regress(lonf1.ndbc_sea_t.data(ix),[repmat(1,size(lonf1.ndbc_sea_t.date(ix))),yearfrac(lonf1.ndbc_sea_t.date(ix))],.05);
{'LONF1 Ts ',numel(ix),': ',B(2),' +/- ',(B(2)-BINT(2,1))},

end;


if (0)

%fn = @(x)(find(x.date<datenum(2009,1,1)));
fn = @(x)(find(datenum(2001,1,1)<=x.date & x.date<datenum(2010,1,1)));
disp(fn);

ix = fn(smkf1.ndbc_tide);
[B,Stats,fh,lh]=scatter_fit(yearfrac(smkf1.ndbc_tide.date(ix)),smkf1.ndbc_tide.data(ix).*12*.0254,'Year','Tide [m]'); titlename(['SMKF1 Tide Height ',char(fn)]);
set(lh(1),'Color',[.5,.5,.5],'LineStyle',':');
set(lh(2),'Color','k','LineWidth',2);
%legend hide;
axis([2001,2010,-0.6,1.6]);

ix = fn(lonf1.ndbc_tide);
[B,Stats,fh,lh]=scatter_fit(yearfrac(lonf1.ndbc_tide.date(ix)),lonf1.ndbc_tide.data(ix).*12*.0254,'Year','Tide [m]'); titlename(['LONF1 Tide Height ',char(fn)]);
set(lh(1),'Color',[.5,.5,.5],'LineStyle',':');
set(lh(2),'Color','k','LineWidth',2);
%legend hide;
axis([2001,2010,-0.6,1.6]);

end;



if (0)
% SMKF1 Tide
% scatter_fit(smkf1.ndbc_tide.date,smkf1.ndbc_tide.data.*12*.0254,'Date','Tide [m]'), datetick3; titlename('SMKF1 Tide Height');
scatter_fit(yearfrac(smkf1.ndbc_tide.date),smkf1.ndbc_tide.data.*12*.0254,'Year','Tide [m]'), titlename('SMKF1 Tide Height');
% % Not-so-perfect mix of dates
% fmg; hist(get_jday(smkf1.ndbc_tide.date),366); titlename('SMKF1 tide seasonal data density'); xlim([1,366]);
% find_date_ranges(smkf1.ndbc_tide.date)
ix = find(smkf1.ndbc_tide.date<datenum(2013,1,1));
[B,Stats,fh,lh]=scatter_fit(yearfrac(smkf1.ndbc_tide.date(ix)),smkf1.ndbc_tide.data(ix).*12*.0254,'Year','Tide [m]'), titlename('SMKF1 Tide Height <2013');

ix = find(smkf1.ndbc_tide.date>=datenum(2001,07,31));
scatter_fit(yearfrac(smkf1.ndbc_tide.date(ix)),smkf1.ndbc_tide.data(ix).*12*.0254,'Year','Tide [m]'), titlename('SMKF1 Tide Height >Jul-2001');


% LONF1 Tide
% scatter_fit(lonf1.ndbc_tide.date,lonf1.ndbc_tide.data.*12*.0254,'Date','Tide [m]'), datetick3; titlename('LONF1 Tide Height');
scatter_fit(yearfrac(lonf1.ndbc_tide.date),lonf1.ndbc_tide.data.*12*.0254,'Year','Tide [m]'), titlename('LONF1 Tide Height');
% % Not-so-perfect mix of dates
% fmg; hist(get_jday(lonf1.ndbc_tide.date),366); titlename('LONF1 tide seasonal data density'); xlim([1,366]);
% find_date_ranges(lonf1.ndbc_tide.date)
% Limiting to whole years
ix = find(lonf1.ndbc_tide.date<datenum(2008,1,1));
scatter_fit(yearfrac(lonf1.ndbc_tide.date(ix)),lonf1.ndbc_tide.data(ix).*12*.0254,'Year','Tide [m]'), titlename('LONF1 Tide Height');
end;

if (0)
% SMKF1 Ts
% scatter_fit(smkf1.ndbc_sea_t.date,smkf1.ndbc_sea_t.data,'Date','T_s [^oC]'), titlename('SMKF1 Sea Temperature'); datetick3;
scatter_fit(yearfrac(smkf1.ndbc_sea_t.date),smkf1.ndbc_sea_t.data,'Year','T_s [^oC]'), titlename('SMKF1 Sea Temperature');
% % Nearly perfect mix of dates
% fmg; hist(get_jday(smkf1.ndbc_sea_t.date),366); titlename('SMKF1 T_s seasonal data density'); xlim([1,366]);


% LONF1 Ts
scatter_fit(lonf1.ndbc_sea_t.date,lonf1.ndbc_sea_t.data,'Date','T_s [^oC]'), titlename('LONF1 Sea Temperature'); datetick3;
scatter_fit(yearfrac(lonf1.ndbc_sea_t.date),lonf1.ndbc_sea_t.data,'Year','T_s [^oC]'), titlename('LONF1 Sea Temperature');
% Not-so-perfect mix of dates
fmg; hist(get_jday(lonf1.ndbc_sea_t.date),366); titlename('LONF1 T_s seasonal data density'); xlim([1,366]);
find_date_ranges(lonf1.ndbc_sea_t.date)
end;
