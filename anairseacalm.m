1;

stn = verify_variable(stn,{'ndbc_air_t_1_d_avg','ndbc_sea_t_1_d_avg','ndbc_wind1_speed_1_d_avg'});

nd = ts_op(stn.ndbc_air_t_1_d_avg,stn.ndbc_sea_t_1_d_avg,'-',@(x)(intersect_dates(x.date,stn.daily_oaflux_air_t.date)));
[nix,wix] = intersect_dates(nd.date,stn.ndbc_wind1_speed_1_d_avg.date(stn.ndbc_wind1_speed_1_d_avg.data<7));

nd.date = nd.date(nix);
nd.data = nd.data(nix);

od = ts_op(stn.daily_oaflux_air_t,stn.daily_oaflux_seatemp,'-',@(x)(intersect_dates(x.date,nd.date)));

fmg; grpplot_ts(nd,@get_month,@nanmean,0,'b'); grpplot_ts(od,@get_month,@nanmean,0,'r'); legend('NDBC','OAFlux'); ylabel('K'); titlename([upper(stn.station_name),' Air-Sea Temperature Difference']);
