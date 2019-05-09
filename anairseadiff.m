function stn = anairseadiff(stn_or_stnm)
%function stn = anairseadiff(stn_or_stnm)
%
% Analyze seasonality in air-sea temperature and specific humidity differences

  stn = get_station_from_station_name(stn_or_stnm);
  if ( ~isfield(stn,'ndbc_air_t') )
    stn = load_all_ndbc_data(stn);
  end;
  if ( ~isfield(stn,'daily_oaflux_air_t') )
    stn = station_load_oaflux(stn);
  end;

  stn = verify_variable(stn,{'ndbc_sea_t_1_d_avg','ndbc_air_t_1_d_avg'});
  nd = ts_op(stn.ndbc_air_t_1_d_avg,stn.ndbc_sea_t_1_d_avg,'-',@(x)(intersect_dates(x.date,stn.daily_oaflux_air_t.date)));
  od = ts_op(stn.daily_oaflux_air_t,stn.daily_oaflux_seatemp,'-',@(x)(intersect_dates(x.date,stn.ndbc_air_t_1_d_avg.date)));

  fmg; grpplot_ts(nd,@get_week,@nanmean,0,'b'); grpplot_ts(od,@get_week,@nanmean,0,'r'); legend('NDBC','OAFlux'); ylabel('K'); titlename([upper(stn.station_name),' Air-Sea Temperature Difference']);
  ylim([-5,1]);

  if ( isfield(stn,'ndbc_dew_t') )
    stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
    stn = station_relhumid_to_spechumid(stn,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
    stn = station_relhumid_to_spechumid(stn,'ndbc_sea_t',100,'ndbc_sea_spechumid');
    % Stommel's salinity adjustment 0.98
    stn.ndbc_sea_spechumid.data = stn.ndbc_sea_spechumid.data .* 0.98;
    stn.daily_oaflux_spechumid.data = stn.daily_oaflux_spechumid.data./1e6;
    stn = station_relhumid_to_spechumid(stn,'daily_oaflux_seatemp',100,'daily_oaflux_sea_spechumid');
    % Stommel's salinity adjustment 0.98
    stn.daily_oaflux_sea_spechumid.data = stn.daily_oaflux_sea_spechumid.data .* 0.98;

    stn = verify_variable(stn,{'ndbc_spechumid_1_d_avg','ndbc_sea_spechumid_1_d_avg'});
    nhd = ts_op(stn.ndbc_spechumid_1_d_avg,stn.ndbc_sea_spechumid_1_d_avg,'-',@(x)(intersect_dates(x.date,stn.daily_oaflux_spechumid.date)));
    ohd = ts_op(stn.daily_oaflux_spechumid,stn.daily_oaflux_sea_spechumid,'-',@(x)(intersect_dates(x.date,stn.ndbc_spechumid_1_d_avg.date)));

    fmg; grpplot_ts(nhd,@get_week,@nanmean,0,'b'); grpplot_ts(ohd,@get_week,@nanmean,0,'r'); legend('NDBC','OAFlux'); ylabel('kg/kg'); titlename([upper(stn.station_name),' Air-Sea Humidity Difference']);
    ylim([-14,2]*1e-3);

    fmg; lh=boxplot_ts(nhd,@get_week,'allcol','b'); rh=boxplot_ts(ohd,@get_week,'allcol','r'); legend([lh(1),rh(1)],'NDBC','OAFlux','Location','South'); ylabel('kg/kg'); titlename([upper(stn.station_name),' Air-Sea Humidity Difference']);
    ylim([-14,2]*1e-3);
  end;

return;
