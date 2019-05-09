1;

if ( ~exist('lonf1','var') )
  lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1);
end;
if ( ~exist('smkf1','var') )
  smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1);
end;

[lts,sts] = intersect_tses(lonf1.ndbc_sea_t,smkf1.ndbc_sea_t);

% fmg; grpplot_ts(lts,@get_yearday);
% fmg; grpplot_ts(sts,@get_yearday);

[cum,tid] = grp_ts(lts.data,lts.date,@get_yearday,@nanmean,0);
cum(tid>=365)=[]; tid(tid>=365)=[];
fmg; plot(tid,cum); titlename('LONF1 climatological sea temperature');
axis([0,366,18,34]);

[cum,tid] = grp_ts(lts.data,lts.date,@get_jday,@nanmean,0);
cum(tid>=365)=[]; tid(tid>=365)=[];
fmg; plot(tid,cum); titlename('LONF1 year-day sea temperature');
axis([0,366,18,34]);

[cum,tid] = grp_ts(sts.data,sts.date,@get_yearday,@nanmean,0);
cum(tid>=365)=[]; tid(tid>=365)=[];
fmg; plot(tid,cum); titlename('SMKF1 climatological sea temperature');
axis([0,366,18,34]);

[cum,tid] = grp_ts(sts.data,sts.date,@get_jday,@nanmean,0);
cum(tid>=365)=[]; tid(tid>=365)=[];
fmg; plot(tid,cum); titlename('SMKF1 year-day sea temperature');
axis([0,366,18,34]);
