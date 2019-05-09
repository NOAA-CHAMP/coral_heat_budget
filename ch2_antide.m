1;

if ( ~exist('smkf1','var') )
  smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1);
end;
if ( ~isfield(smkf1,'ndbc_tide_m') )
  smkf1.ndbc_tide_m=smkf1.ndbc_tide; smkf1.ndbc_tide_m.data=unitsratio('m','ft').*(smkf1.ndbc_tide_m.data-nanmean(smkf1.ndbc_tide_m.data));
end;

if ( ~exist('lonf1','var') )
  lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1);
end;
if ( ~isfield(lonf1,'ndbc_tide_m') )
  lonf1.ndbc_tide_m=lonf1.ndbc_tide; lonf1.ndbc_tide_m.data=unitsratio('m','ft').*(lonf1.ndbc_tide_m.data-nanmean(lonf1.ndbc_tide_m.data));
end;

fmg; boxplot_ts(lonf1.ndbc_tide_m,[],'mean',true); ylim([-1.1,1.6]); titlename('LONF1');
fmg; boxplot_ts(smkf1.ndbc_tide_m,[],'mean',true); ylim([-1.1,1.6]); titlename('SMKF1');

%%%% T_TIDE harmonic analysis and analysis of residuals
ch2_tide;


% fmg; boxplot_ts(lonf1.ndbc_tide_m,'week'); ylim([-1.1,1.6]); titlename('LONF1');
% fmg; boxplot_ts(smkf1.ndbc_tide_m,'week'); ylim([-1.1,1.6]); titlename('SMKF1');

% [l,s] = intersect_tses(lonf1.ndbc_tide_m,smkf1.ndbc_tide_m);
% fmg; boxplot_ts(l,[],'mean',true); ylim([-1.1,1.6]); titlename('LONF1 inter');
% fmg; boxplot_ts(s,[],'mean',true); ylim([-1.1,1.6]); titlename('SMKF1 inter');
