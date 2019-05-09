1;

% AVG = @nanmean;
AVG = @nanmedian;

% stnm = 'fwyf1';
% stnm = 'mlrf1';
% stnm = 'lonf1';
% stnm = 'smkf1';
stnm = 'sanf1';

%FWYF1 mean annual 23.0-29.5 JF- A: peak above 1.2 NDJ, mean 0.9 MAMJJ, low 0.5 SO
%MLRF1 mean annual 23.5-30.0 JF- A: peak above 1.1 JFM, mean 0.9 AMJ, low 0.6 JASON
%LONF1 mean annual 20.0-31.0 J -JA: peak above 1.1 DJF, mean 0.9 AMJ, low 0.6 ASON
%SMKF1 mean annual 23.5-30.5 JF- A: peak above 1.5 DJF, mean 1.1 AMJJA, low 0.6 Oc?
%SANF1 mean annual 22.5-30.0 JF- A: peak above 1.1 JF, mean 0.8 MAMJA, low 0.4 ON


ix=1;
if ( ~exist('stns','var') )
  stns = { get_station_from_station_name(stnm) };
  stns{ix} = load_all_ndbc_data(stns{ix});
end;


% Analyze mean and mean diurnal amplitude of sea temperature by year-day

[cum,tid]=grp_ts(stns{ix}.ndbc_sea_t.data,stns{ix}.ndbc_sea_t.date,@get_yearday,AVG,0);
fmg; plot(tid,cum); axis([0,366,18,34]);
ttl = [upper(stns{ix}.station_name) ' Interannual range night-day'];
titlename(ttl); disp(ttl);

[mxcum,mxtid]=grp_ts(stns{ix}.ndbc_sea_t.data,stns{ix}.ndbc_sea_t.date,@floor,@max,20);
[mdcum,mdtid]=grp_ts(stns{ix}.ndbc_sea_t.data,stns{ix}.ndbc_sea_t.date,@floor,AVG,20);
[mncum,mntid]=grp_ts(stns{ix}.ndbc_sea_t.data,stns{ix}.ndbc_sea_t.date,@floor,@min,20);
[jdfcum,jdftid]=grp_ts(mxcum-mncum,mxtid,@get_jday,AVG,0);
fmg; plot(jdftid,jdfcum); axis([0,366,0.4,2.0]);
datetick3('x',3,'keeplimits');
ttl = [upper(stns{ix}.station_name) ' Diurnal variability'];
titlename(ttl); disp(ttl);

[jmdcum,jmdtid]=grp_ts(mdcum,mdtid,@get_jday,AVG,0);

[jdxcum,jdxtid]=grp_ts(mxcum-mdcum,mxtid,@get_jday,AVG,0);
[jdncum,jdntid]=grp_ts(mdcum-mncum,mxtid,@get_jday,AVG,0);
fmg; plot(jmdtid,jmdcum,jmdtid,jmdcum+jdxcum,jmdtid,jmdcum-jdncum); axis([0,366,18,34]);
ttl = [upper(stns{ix}.station_name) ' Mean diurnal range'];
titlename(ttl); disp(ttl);

[jmxcum,jmxtid]=grp_ts(mxcum,mxtid,@get_jday,AVG,0);
[jmncum,jmntid]=grp_ts(mncum,mntid,@get_jday,AVG,0);
fmg; plot(jmdtid,jmdcum,jmxtid,jmxcum,jmntid,jmncum); axis([0,366,18,34]);
ttl = [upper(stns{ix}.station_name) ' Interannual range min-max'];
titlename(ttl); disp(ttl);
