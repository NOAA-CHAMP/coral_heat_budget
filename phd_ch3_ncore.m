1;

more off;

if ( ~exist('klgf1','var') )
  [klgf1,mrtf1]=an_ncore;

  klgf1.cm_shallow_seatemp_diff.date = klgf1.cm_shallow_seatemp.date(1:end-1);
  klgf1.cm_shallow_seatemp_diff.data = diff(klgf1.cm_shallow_seatemp.data);
  klgf1 = filter_gaps(klgf1,'cm_shallow_seatemp','cm_shallow_seatemp_diff',(1.5/24));
  klgf1.cm_deep_seatemp_diff.date = klgf1.cm_deep_seatemp.date(1:end-1);
  klgf1.cm_deep_seatemp_diff.data = diff(klgf1.cm_deep_seatemp.data);
  klgf1 = filter_gaps(klgf1,'cm_deep_seatemp','cm_deep_seatemp_diff',(1.5/24));

  klgf1.cm_sheer_xshore_10min = ts_op(klgf1.cm_shallow_xshore_10min,klgf1.cm_deep_xshore_10min,'-',[],[],[],(4.5/60/24));
  klgf1.cm_dTdz_10min = ts_op(klgf1.cm_shallow_seatemp_10min,klgf1.cm_deep_seatemp_10min,'-',[],[],[],(4.5/60/24));
  klgf1.cm_dxTdz_10min = ts_op(klgf1.cm_sheer_xshore_10min,klgf1.cm_dTdz_10min,'*',[],[],[],(4.5/60/24));
  klgf1.cm_sheer_xshore = ts_op(klgf1.cm_shallow_xshore,klgf1.cm_deep_xshore,'-');
  klgf1.cm_dTdz = ts_op(klgf1.cm_shallow_seatemp,klgf1.cm_deep_seatemp,'-');
  klgf1.cm_dxTdz = ts_op(klgf1.cm_sheer_xshore,klgf1.cm_dTdz,'*');
  klgf1 = verify_variable(klgf1,'cm_sheer_xshore_1_d_avg');
  klgf1 = verify_variable(klgf1,'cm_sheer_xshore_40_h_lp');
  klgf1 = verify_variable(klgf1,'cm_dxTdz_1_d_avg');
  klgf1 = verify_variable(klgf1,'cm_dxTdz_40_h_lp');

  klgf1.cm_deep_xdTdz = ts_op(klgf1.cm_deep_xshore,klgf1.cm_dTdz,'*');
  klgf1.cm_deep_xdTdz_10min = ts_op(klgf1.cm_deep_xshore_10min,klgf1.cm_dTdz_10min,'*',[],[],[],(4.5/60/24));


  % How much of cross-shore current sheer is explicable by tides?
  klgf1.cm_deep_seapres_dt = ts_dt(klgf1.cm_deep_seapres);
  %scatter_fit_ts(klgf1.cm_deep_seapres_dt,klgf1.cm_sheer_xshore,[],[],'KLGF1 dP_S/dt','\Deltau Cross-shore');


  mrtf1.cm_shallow_seatemp_diff.date = mrtf1.cm_shallow_seatemp.date(1:end-1);
  mrtf1.cm_shallow_seatemp_diff.data = diff(mrtf1.cm_shallow_seatemp.data);
  mrtf1 = filter_gaps(mrtf1,'cm_shallow_seatemp','cm_shallow_seatemp_diff',(1.5/24));
  mrtf1.cm_deep_seatemp_diff.date = mrtf1.cm_deep_seatemp.date(1:end-1);
  mrtf1.cm_deep_seatemp_diff.data = diff(mrtf1.cm_deep_seatemp.data);
  mrtf1 = filter_gaps(mrtf1,'cm_deep_seatemp','cm_deep_seatemp_diff',(1.5/24));

  mrtf1.cm_sheer_xshore_10min = ts_op(mrtf1.cm_shallow_xshore_10min,mrtf1.cm_deep_xshore_10min,'-',[],[],[],(4.5/60/24));
  mrtf1.cm_dTdz_10min = ts_op(mrtf1.cm_shallow_seatemp_10min,mrtf1.cm_deep_seatemp_10min,'-',[],[],[],(4.5/60/24));
  mrtf1.cm_dxTdz_10min = ts_op(mrtf1.cm_sheer_xshore_10min,mrtf1.cm_dTdz_10min,'*',[],[],[],(4.5/60/24));
  mrtf1.cm_sheer_xshore = ts_op(mrtf1.cm_shallow_xshore,mrtf1.cm_deep_xshore,'-');
  mrtf1.cm_dTdz = ts_op(mrtf1.cm_shallow_seatemp,mrtf1.cm_deep_seatemp,'-');
  mrtf1.cm_dxTdz = ts_op(mrtf1.cm_sheer_xshore,mrtf1.cm_dTdz,'*');
  mrtf1 = verify_variable(mrtf1,'cm_sheer_xshore_1_d_avg');
  mrtf1 = verify_variable(mrtf1,'cm_sheer_xshore_40_h_lp');
  mrtf1 = verify_variable(mrtf1,'cm_dxTdz_1_d_avg');
  mrtf1 = verify_variable(mrtf1,'cm_dxTdz_40_h_lp');

  mrtf1.cm_deep_xdTdz = ts_op(mrtf1.cm_deep_xshore,mrtf1.cm_dTdz,'*');
  mrtf1.cm_deep_xdTdz_10min = ts_op(mrtf1.cm_deep_xshore_10min,mrtf1.cm_dTdz_10min,'*',[],[],[],(4.5/60/24));

  % How much of cross-shore current sheer is explicable by tides?
  mrtf1.cm_deep_seapres_dt = ts_dt(mrtf1.cm_deep_seapres);
  %scatter_fit_ts(mrtf1.cm_deep_seapres_dt,mrtf1.cm_sheer_xshore,[],[],'MRTF1 dP_S/dt','\Deltau Cross-shore');

end;

if ( ~exist('stn','var') || ~isfield(stn,'simple_ndbc_erai_erai_30a_net_flux') || ~strcmpi(stn.station_name,'MLRF1') )
  stn=[]; clear stn
  matfname = fullfile(get_thesis_path('../data'),['mlrf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat']);
  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname,'stn');
  else
    stn = optimize_station_heat_budget('mlrf1','erai','avhrr_weekly','ndbc','tpxo_tide','erai');
    disp(['SAVING ',matfname]);
    save(matfname,'stn');
  end;
  clear matfname

  x=load_station_data('mlrf1');
  stn.bic_surf_par = x.bic_surf_par;
  x=[]; clear x;

  x = get_ncep_station('mlrf1','narr');
  stn.ncep_par = x.ncep_par;
  stn.ncep_dsrf = x.ncep_dsrf;
  stn.ncep_srf = x.ncep_srf;
  x=[]; clear x;

  %stn = station_optimal_isobath_orientation(stn);
  stn = station_reorient_vectors(stn,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');

  [stn.ndbc_wind1_xshore_3hlp.date,stn.ndbc_wind1_xshore_3hlp.data] = lanczos_ts(stn.ndbc_wind1_xshore.date,stn.ndbc_wind1_xshore.data);
  [stn.ndbc_wind1_lshore_3hlp.date,stn.ndbc_wind1_lshore_3hlp.data] = lanczos_ts(stn.ndbc_wind1_lshore.date,stn.ndbc_wind1_xshore.data);
end;



stn = verify_variable(stn,'ndbc_wind1_xshore_1_d_var');
stn = verify_variable(stn,'ndbc_wind1_lshore_1_d_var');
stn = verify_variable(stn,'ndbc_wind1_speed_1_d_var');
stn = verify_variable(stn,'simple_ndbc_erai_erai_30a_net_flux_1_d_var');

klgf1 = verify_variable(klgf1,'cm_shallow_seatemp_1_d_var');
klgf1 = verify_variable(klgf1,'cm_deep_seatemp_1_d_var');
klgf1 = verify_variable(klgf1,'cm_deep_xshore_1_d_var');
klgf1 = verify_variable(klgf1,'cm_deep_lshore_1_d_var');
klgf1 = verify_variable(klgf1,'cm_deep_speed_1_d_var');
klgf1 = verify_variable(klgf1,'cm_sheer_xshore_1_d_var');

mrtf1 = verify_variable(mrtf1,'cm_shallow_seatemp_1_d_var');
mrtf1 = verify_variable(mrtf1,'cm_deep_seatemp_1_d_var');
mrtf1 = verify_variable(mrtf1,'cm_deep_xshore_1_d_var');
mrtf1 = verify_variable(mrtf1,'cm_deep_lshore_1_d_var');
mrtf1 = verify_variable(mrtf1,'cm_deep_speed_1_d_var');
mrtf1 = verify_variable(mrtf1,'cm_sheer_xshore_1_d_var');


if (0)

  %fmg; plot_ts(klgf1.cm_shallow_seapres_raw,klgf1.cm_deep_seapres_raw);

  fmg;
  legs={};
  plot_ts(klgf1.cm_sheer_xshore); legs{end+1} = '\deltau/\deltaz';
  plot_ts(klgf1.cm_sheer_xshore_40_h_lp,'m--','LineWidth',1.5); legs{end+1} = '\deltau/\deltaz 40hlp';
  %plot_ts(klgf1.cm_dxTdz_40_h_lp,'k-.','LineWidth',1.5); legs{end+1} = '\delta(u^.T)/\deltaz 40hlp';
  plot_ts(klgf1.cm_dxTdz_1_d_avg,'k-.','LineWidth',1.5); legs{end+1} = '\mu_1_d\delta(u^.T)/\deltaz';
  plot_ts(ts_op(stn.ndbc_wind1_speed,100,'/'),'r'); legs{end+1} = 'Wind/100';
  plot_ts(ts_op(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_avg,1000,'/'),'color',[0,1,0]); legs{end+1} = '\mu_1_dQ_0/1000';
  legend(legs);
  xlim(datenum([2007,2008],[1,12],[1,1])); datetick3;


  scatter_fit_ts(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_avg,klgf1.cm_dxTdz_1_d_avg,[],[],'\mu_1_dQ_0','\mu_1_d\delta(u^.T)/\deltaz');

  fmg;
  plot_ts(klgf1.cm_shallow_seatemp_10min,...
          ts_op(klgf1.cm_shallow_xshore_10min,100,'*'),...
          klgf1.cm_deep_seatemp_10min,...
          ts_op(klgf1.cm_deep_xshore_10min,100,'*'),...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/'),...
          stn.ndbc_wind1_speed...
          );
  legend('T_7','u_7','T_2_2','u_2_2','Q_0/100','Wind');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,7,[7,16])); datetick3;
  titlename('KLGF1 summer warming response');

  fmg;
  plot_ts(klgf1.cm_dTdz,...
          ts_op(klgf1.cm_shallow_xshore_10min,10,'*'),...
          ts_op(klgf1.cm_deep_xshore_10min,10,'*'),...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/'),...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/')...
          );
  legend('T_7-T_2_2','u_7','u_2_2','Q_0/100','\partial_tT/100');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,7,[7,16])); datetick3;
  titlename('KLGF1 summer warming response');

  fmg;
  plot_ts(klgf1.cm_dTdz,...
          ts_op(klgf1.cm_shallow_xshore,10,'*'),...
          ts_op(klgf1.cm_deep_xshore,10,'*'),...
          ts_op(klgf1.cm_sheer_xshore,10,'*'),...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/'),...
          ts_op(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux,100,'/')...
          );
  legend('T_7-T_2_2','u_7','u_2_2','u_7-u_2_2','Q_0/100','\partial_tT/100');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,7,[7,16])); datetick3;
  titlename('KLGF1 summer warming response');

  fmg;
  plot_ts(klgf1.cm_dTdz,...
          ts_op(klgf1.cm_sheer_xshore,10,'*'),...
          ts_op(ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux,'-'),100,'/')...
          );
  legend('T_7-T_2_2','u_7-u_2_2','(Q_0-\partial_tT)/100');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,7,[7,16])); datetick3;
  titlename('KLGF1 summer warming response');

  fmg;
  plot_ts(klgf1.cm_dxTdz,...
          ts_op(ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux,'-'),100,'/')...
          );
  legend('(T_7-T_2_2) \times (u_7-u_2_2)','(Q_0-\partial_tT)/100');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,7,[7,16])); datetick3;
  titlename('KLGF1 summer warming response');

  fmg;
  plot_ts(klgf1.cm_shallow_xshore,'b',...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux_term,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt,'-')...
          );
  legend('u_7','(Q_0/\rhoC_ph)-\partial_tT');
  xlim(datenum(2007,8,[22,32])); datetick3;
  titlename('KLGF1 summer warming response');


  dts=datenum(2007,8,[1,32]);
  scatter_fit_ts(ts_op(stn.simple_ndbc_erai_erai_30a_net_flux_term,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt,'-'),...
                 klgf1.cm_shallow_xshore,[],ts_date_range(klgf1.cm_shallow_xshore,dts),...
                 '(Q_0/\rhoC_ph)-\partial_tT',['u_7 ',datestr(dts(1)),'-',datestr(dts(2))]);


  fmg; hist(klgf1.cm_sheer_xshore.data(ts_boreal_cool(klgf1.cm_sheer_xshore)),1000); xlim([-.5,.5]); titlename('Klg Cool');
  fmg; hist(klgf1.cm_sheer_xshore.data(ts_boreal_warm(klgf1.cm_sheer_xshore)),1000); xlim([-.5,.5]); titlename('Klg Warm');

  fmg; hist(mrtf1.cm_sheer_xshore.data(ts_boreal_cool(mrtf1.cm_sheer_xshore)),1000); xlim([-.5,.5]); titlename('Mrt Cool');
  fmg; hist(mrtf1.cm_sheer_xshore.data(ts_boreal_warm(mrtf1.cm_sheer_xshore)),1000); xlim([-.5,.5]); titlename('Mrt Warm');

  scatter_fit_ts(stn.ndbc_erai_erai_30a_wind_stress_xshore,klgf1.cm_shallow_xshore,ts_date_range(stn.ndbc_erai_erai_30a_wind_stress_xshore,datenum([2007,2008],[11,4],1)),[],'\tau^x^s Cool','u');


  fmg;
  plot_ts(ts_op(klgf1.cm_shallow_seatemp,30,'-'),...
          ts_op(klgf1.cm_shallow_xshore,100,'*'),...
          ts_op(klgf1.cm_deep_seatemp,30,'-'),...
          ts_op(klgf1.cm_deep_xshore,100,'*'),...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/'),...
          stn.ndbc_wind1_speed,...
          ts_op(stn.ndbc_sea_t,30,'-'),'k-'...
          );
  legend('T_7-30','u_7','T_2_2-30','u_2_2','Q_0/100','Wind','T_s-30');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,7,[7,16])); datetick3;
  titlename('KLGF1 summer warming response');


  fmg;
  plot_ts(klgf1.cm_shallow_seatemp_diff,...
          klgf1.cm_shallow_xshore,...
          klgf1.cm_deep_seatemp_diff,...
          klgf1.cm_deep_xshore,...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux_term,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt,'-'),...
          stn.ndbc_erai_erai_30a_wind_stress,...
          stn.ndbc_sea_t_diff,'k-'...
          );
  legend('T_7-30','u_7','T_2_2-30','u_2_2','Q_0/\rhoC_ph','\tau','\deltaT_s/\deltat');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,7,[7,16])); datetick3;
  titlename('KLGF1 summer warming response');


  fmg;
  plot_ts(klgf1.cm_deep_seatemp_1_d_var,...
          ts_op(klgf1.cm_deep_xshore,100,'*'),...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/'),...
          stn.ndbc_wind1_xshore,...
          stn.ndbc_wind1_lshore,...
          stn.ndbc_wind1_speed...
          );
  legend('\sigma_1_dT_2_2','u^x^s_2_2\times100','Q_0/100','U^x^s','U^l^s','|U|');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,[7,9],[1,15])); datetick3;
  titlename('KLGF1 summer warming response at depth');

  fmg;
  plot_ts(klgf1.cm_shallow_seatemp_10min,...
          ts_op(klgf1.cm_shallow_xshore_10min,100,'*'),...
          klgf1.cm_deep_seatemp_10min,...
          ts_op(klgf1.cm_deep_xshore_10min,100,'*'),...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/'),...
          stn.ndbc_wind1_xshore,...
          stn.ndbc_wind1_lshore...
          );
  legend('T_7','u^x^s_7','T_2_2','u^x^s_2_2','Q_0/100','U^x^s','U^l^s');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,7,[7,16])); datetick3;
  titlename('KLGF1 summer warming response');


  fmg;
  plot_ts(ts_op(klgf1.cm_shallow_seatemp_1_d_var,10,'*'),...
          ts_op(klgf1.cm_deep_seatemp_1_d_var,10,'*'),...
          ts_op(klgf1.cm_deep_xshore_1_d_var,1000,'*'),...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_var,10000,'/'),...
          ts_op(stn.ndbc_wind1_xshore_1_d_var,10,'/'),...
          ts_op(stn.ndbc_wind1_speed_1_d_var,10,'/')...
          );
  legend('\sigma_1_dT_7\times10','\sigma_1_dT_2_2\times10','\sigma_1_du^x^s_2_2\times1000','\sigma_1_dQ_0/10000','\sigma_1_dU^x^s/10','\sigma_1_d|U|/10');
  %xlim(datenum(2007,6,[6,15])); datetick3;
  xlim(datenum(2007,[7,9],[1,15])); datetick3;
  titlename('KLGF1 summer warming response at depth');


  fs = @(x)(find(ismember(get_hour(x.date),[7,19])));
  fmg; plot_ts(subset_ts(klgf1.cm_dTdz,fs),subset_ts(klgf1.cm_sheer_xshore,fs),subset_ts(stn.simple_ndbc_erai_erai_30a_net_flux,fs)); legend('dTdz','Sheer','Q0/\rhoC_ph');


  fmg; plot_ts(stn.simple_ndbc_erai_erai_30a_net_flux_dly,stn.ndbc_erai_erai_30a_net_flux_dly,stn.b_ndbc_erai_erai_30a_net_flux_dly,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly,stn.ndbc_sea_t_dly_diff_flux); legend('Simple','Sfc.','Sfc./Btm.','HC','Implied');


  fmg;
  fs = @(x)(find(datenum(2007,7,1)<=x.date&x.date<=datenum(2007,7,22)));
  ts=subset_ts(klgf1.cm_shallow_seatemp,fs); plot(get_hour(ts.date),ts.data,'b.');
  ts=subset_ts(klgf1.cm_deep_seatemp,fs); plot(get_hour(ts.date),ts.data,'r.');
  ts=subset_ts(stn.simple_ndbc_erai_erai_30a_net_flux,fs); plot(get_hour(ts.date),ts.data./100,'m.');
  ts=subset_ts(klgf1.cm_sheer_xshore,fs); plot(get_hour(ts.date),ts.data.*100,'.','Color',[0,.5,0]);


  fmg;
  fs = @(x)(find(datenum(2007,7,1)<=x.date&x.date<=datenum(2007,7,22)));
  ts=subset_ts(klgf1.cm_dxTdz,fs); plot(get_hour(ts.date),ts.data.*10,'b.');
  ts=subset_ts(stn.simple_ndbc_erai_erai_30a_net_flux,fs); plot(get_hour(ts.date)+0.5,ts.data./100,'m.');

end; %if (0)



if (1)
  % % lagdif.date = stn.ndbc_sea_t_dly_diff_flux.date(2:end);
  % % lagdif.data = stn.ndbc_sea_t_dly_diff_flux.data(1:end-1);
  % % fmg; plot_ts(stn.simple_ndbc_erai_erai_30a_net_flux_dly,stn.ndbc_erai_erai_30a_net_flux_dly,stn.b_ndbc_erai_erai_30a_net_flux_dly,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly,lagdif,'k-o','LineWidth',1.5); legend('Simple','Sfc.','Sfc./Btm.','HC','Implied');
  % % % xlim(datenum(2007,[7,10],1));
  % % xlim(datenum(2008,[7,10],1));
  % % datetick3;
  % % titlename([upper(stn.station_name),' daily net heat fluxes (budget versions)']);

  % lagdif.date = stn.ndbc_sea_t_dly_diff.date(2:end);
  % lagdif.data = stn.ndbc_sea_t_dly_diff.data(1:end-1);
  % [sqt,qt,bqt,hc,dt] = intersect_tses(stn.simple_ndbc_erai_erai_30a_net_flux_term_dly,stn.ndbc_erai_erai_30a_net_flux_term_dly,stn.b_ndbc_erai_erai_30a_net_flux_term_dly,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly,lagdif);
  % fmg; plot_ts(sqt,qt,bqt,hc,dt,'k-o','LineWidth',1.5); legend('Simple','Sfc.','Sfc./Btm.','HC','Implied');
  % % xlim(datenum(2007,[7,10],1));
  % xlim(datenum(2008,[7,10],1));
  % datetick3;
  % titlename([upper(stn.station_name),' daily temperature change (budget versions)']);

  % (Corrects a bug NOW FIXED in OPTIMIZE_STATION_HEAT_BUDGET.m - see TS_DT.m in Ecoforecasts toolkit)
  lagdif.date = stn.ndbc_sea_t_dly_diff.date(2:end);
  lagdif.data = stn.ndbc_sea_t_dly_diff.data(1:end-1);

  [sqt,bqt,dt,hc,dft] = intersect_tses(stn.simple_ndbc_erai_erai_30a_net_flux_term_dly,stn.b_ndbc_erai_erai_30a_net_flux_term_dly,stn.b_ndbc_erai_erai_30a_avhrr_dt_dly,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly,lagdif);

  fmg; plot_ts(sqt,bqt,dt,hc,dft,'k-o','LineWidth',1.5); legend('Simple','Sfc./Btm.','Sfc./Btm./Adv./Dif.','HC','Implied');
  % xlim(datenum(2007,[7,10],1));
  xlim(datenum(2008,[7,10],1));
  datetick3;
  titlename([upper(stn.station_name),' daily temperature change (budget versions)']);


  fmg;
  plot_ts(stn.ndbc_sea_t,'k','LineWidth',1.5,...
          klgf1.cm_shallow_seatemp,...
          mrtf1.cm_shallow_seatemp,':','Color',[0,.5,0],'LineWidth',1.5,...
          ts_op(klgf1.cm_sheer_xshore,100,'*'),'r',...
          klgf1.cm_deep_seatemp,'Color',[0,.75,.75],...
          mrtf1.cm_deep_seatemp,'Color',[0,.75,.75],'LineWidth',1.5,...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/'),...
          klgf1.cm_deep_seapres,'Color',[.5,.5,.5],'LineWidth',1.5,...
          stn.ndbc_wind1_xshore,'k','LineWidth',1.5...
          );
  legend('T_3','T_7','Mrt T_7','u^x^s_7-u^x^s_2_2','T_2_2','Mrt T_2_2','Q_0/100','P_2_2','U^x^s', 'Location','SouthEast');
  xlim(datenum(2007,7,[1,40])); datetick3;
  % % First cooling event of the season - 2007 Oct 01
  % xlim(datenum(2007,9,[15,45])); datetick3;
  % % xlim(datenum(2008,7,[1,40])); datetick3;
  titlename('NCORE summer warming response - direct observations');

  fmg;
  plot_ts(ts_op(klgf1.cm_shallow_seatemp_1_d_var,10,'*'),...
          ts_op(klgf1.cm_deep_seatemp_1_d_var,10,'*'),...
          ts_op(klgf1.cm_deep_xshore_1_d_var,1000,'*'),...
          ts_op(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_var,10000,'/'),...
          ts_op(stn.ndbc_wind1_xshore_1_d_var,10,'/'),...
          ts_op(stn.ndbc_wind1_speed_1_d_var,10,'/')...
          );
  legend('\sigma_1_dT_7\times10','\sigma_1_dT_2_2\times10','\sigma_1_du^x^s_2_2\times1000','\sigma_1_dQ_0/10000','\sigma_1_dU^x^s/10','\sigma_1_d|U|/10');
  % %xlim(datenum(2007,6,[6,15])); datetick3;
  % xlim(datenum(2007,[7,9],[1,15])); datetick3;
  xlim(datenum(2008,[7,9],[1,15])); datetick3;
  titlename('KLGF1 summer warming response - daily variability');


  phd_ch3_daily_diagnostic(stn,klgf1);

end; %if (1)


if (0)
  Jul = datenum(2007,7,1);
  %Jul = datenum(2008,7,1);
  Oct = datenum(2007,10,1);
  %Oct = datenum(2008,10,1);

  fs=@(x)(find(get_yearmonth(x.date)==Jul)); mo=datestr(Jul,12);

  x.lat = klgf1.lat;
  %x.ts=subset_ts(klgf1.cm_deep_seapres_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['KLG ',mo,' Deep P']);
  x.ts=subset_ts(klgf1.cm_deep_seatemp_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['KLG ',mo,' Deep T']);
  %x.ts=subset_ts(klgf1.cm_shallow_seapres_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['KLG ',mo,' Shallow P']);
  x.ts=subset_ts(klgf1.cm_shallow_seatemp_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['KLG ',mo,' Shallow T']);

  % x.lat = mrtf1.lat;
  % %x.ts=subset_ts(mrtf1.cm_deep_seapres_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['MRT ',mo,' Deep P']);
  % x.ts=subset_ts(mrtf1.cm_deep_seatemp_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['MRT ',mo,' Deep T']);
  % %x.ts=subset_ts(mrtf1.cm_shallow_seapres_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['MRT ',mo,' Shallow P']);
  % x.ts=subset_ts(mrtf1.cm_shallow_seatemp_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['MRT ',mo,' Shallow T']);


  fs=@(x)(find(get_yearmonth(x.date)==Oct)); mo=datestr(Oct,12); 

  x.lat = klgf1.lat;
  %x.ts=subset_ts(klgf1.cm_deep_seapres_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['KLG ',mo,' Deep P']);
  x.ts=subset_ts(klgf1.cm_deep_seatemp_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['KLG ',mo,' Deep T']);
  %x.ts=subset_ts(klgf1.cm_shallow_seapres_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['KLG ',mo,' Shallow P']);
  x.ts=subset_ts(klgf1.cm_shallow_seatemp_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['KLG ',mo,' Shallow T']);

  % x.lat = mrtf1.lat;
  % %x.ts=subset_ts(mrtf1.cm_deep_seapres_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['MRT ',mo,' Deep P']);
  % x.ts=subset_ts(mrtf1.cm_deep_seatemp_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['MRT ',mo,' Deep T']);
  % %x.ts=subset_ts(mrtf1.cm_shallow_seapres_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['MRT ',mo,' Shallow P']);
  % x.ts=subset_ts(mrtf1.cm_shallow_seatemp_10min,fs); plot_spec(x,'ts',[],[],[],[],[],true); titlename(['MRT ',mo,' Shallow T']);

end; %if (0) %%%% Second block

more on;
