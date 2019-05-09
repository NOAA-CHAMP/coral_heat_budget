1;

  perfun = @floor; minN=24; perstr='daily';
  %perfun = @get_yearweek; minN=7*23; perstr='weekly';

  q0udt = ts_op(stn.b_ndbc_erai_erai_30a_net_flux_term,stn.erai_avhrr_advected_heat,'+');
  full_dt = ts_op(q0udt,stn.avhrr_weekly_diffused_heat,'+');

  [sqt.data,sqt.date] = grp_ts(stn.simple_ndbc_erai_erai_30a_net_flux_term.data,stn.simple_ndbc_erai_erai_30a_net_flux_term.date,perfun,@nansum,minN);
  [bqt.data,bqt.date] = grp_ts(stn.b_ndbc_erai_erai_30a_net_flux_term.data,stn.b_ndbc_erai_erai_30a_net_flux_term.date,perfun,@nansum,minN);
  [dt.data,dt.date] = grp_ts(stn.b_ndbc_erai_erai_30a_avhrr_dt.data,stn.b_ndbc_erai_erai_30a_avhrr_dt.date,perfun,@nansum,minN);
  [full_dt.data,full_dt.date] = grp_ts(full_dt.data,full_dt.date,perfun,@nansum,minN);
  [hchc.data,hchc.date] = grp_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc.data,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc.date,perfun,@nansum,minN);
  [hc.data,hc.date] = grp_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt.data,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt.date,perfun,@nansum,minN);

  [dft.data,dft.date] = grp_ts(stn.ndbc_sea_t_diff.data,stn.ndbc_sea_t_diff.date,perfun,@nansum,minN);
  % dft.date = dft.date(2:end);
  % dft.data = dft.data(1:end-1);

  [sqt,bqt,dt,full_dt,hchc,hc,dft] = intersect_tses(sqt,bqt,dt,full_dt,hchc,hc,dft);

  fmg; plot_ts(sqt,bqt,dt,full_dt,hchc,hc,dft,'k-o','LineWidth',1.5); legend('Simple','Sfc./Btm.','Sfc./Btm./F_qAdv./Dif.','Sfc./Btm./Adv./Dif.','HC','\partial_tT','Implied');
  % xlim(datenum(2001,[7,10],1));
  % xlim(datenum(2007,[7,10],1));
  xlim(datenum(2008,[7,10],1));
  datetick3;
  ylim([-3,3]);
  titlename([upper(stn.station_name),' ',perstr,' temperature change (budget versions)']);
