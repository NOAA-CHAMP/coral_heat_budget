1;

station = []; clear station; clear f ix1 ix2 pct;

load('../data/smkf1_ms.mat');
station = verify_variable(station,'qf_40_hour_lowpass');
station = verify_variable(station,'netqf_40_hour_lowpass');
station = verify_variable(station,'ndbc_sea_t_40_hour_lowpass');
f = findiff_ts(station.ndbc_sea_t_40_hour_lowpass,5);

[ix1,ix2]=intersect_dates(station.qf_40_hour_lowpass.date,f.date);
pct = (abs(station.qf_40_hour_lowpass.data(ix1))./(abs(f.data(ix2))));

length(find(pct < 1)),
figure;
hist(pct(pct<1));
nanmean(pct(pct<1)),

station = []; clear station; clear f ix1 ix2 pct;
