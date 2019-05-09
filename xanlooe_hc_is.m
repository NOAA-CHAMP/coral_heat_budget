1;

%-- 4/5/12 11:02 PM --%
fclose('all'); close all; clear all; clear import; clear classes; clear java; dbclear all; pack;
help trim_looe1_adcp
type trim_looe1_adcp
stn = optimize_station_heat_budget('looe1','erai','avhrr_weekly','ndbc','tpxo_tide','erai','ad_seatemp');
stn = trim_looe1_adcp(stn);
scatter_fit_ts(stn.adcp_btm_x,stn.adseatemp_ndbc_erai_erai_30a_cooling_flux),
scatter_fit_ts(stn.adcp_btm_x,stn.adseatemp_ndbc_erai_erai_30a_cooling_flux_term)
scatter_fit_ts(stn.adcp_btm_x,stn.adseatemp_ndbc_erai_erai_30a_cooling_flux_1_d_avg)
scatter_fit_ts(stn.adcp_v_40hlp,stn.adseatemp_ndbc_erai_erai_30a_cooling_flux_1_d_avg)
scatter_fit_ts(stn.adcp_v_40hlp,stn.adseatemp_ndbc_erai_erai_30a_net_flux_1_d_avg)
scatter_fit_ts(stn.adcp_v_40hlp,stn.adseatemp_ndbc_erai_erai_30a_net_flux)
stn = verify_variable(stn,'adseatemp_ndbc_erai_erai_30a_net_flux_1_d_avg');
scatter_fit_ts(stn.adcp_v_40hlp,stn.adseatemp_ndbc_erai_erai_30a_net_flux_1_d_avg)
25350/60
25350/60/24
25350/60/4
sw_cp(35,27,1000)
pi*(160e3^2)
pi*(80e3^2)
scatter_fit_ts(stn.adcp_u_40hlp,stn.adseatemp_ndbc_erai_erai_30a_net_flux_1_d_avg)
scatter_fit_ts(stn.adcp_u_40_h_decimate,stn.adseatemp_ndbc_erai_erai_30a_net_flux_40_h_decimate)
stn=verify_variable(stn,{'adcp_u_40_h_decimate','adseatemp_ndbc_erai_erai_30a_net_flux_40_h_decimate'});
scatter_fit_ts(stn.adcp_u_40_h_decimate,stn.adseatemp_ndbc_erai_erai_30a_net_flux_40_h_decimate)
stn=verify_variable(stn,{'adcp_u_40_h_decimate','adcp_v_40_h_decimate','adseatemp_ndbc_erai_erai_30a_net_flux_40_h_decimate'});
scatter_fit_ts(stn.adcp_v_40_h_decimate,stn.adseatemp_ndbc_erai_erai_30a_net_flux_40_h_decimate)
scatter_fit_ts(stn.adcp_v_40_h_decimate,stn.adseatemp_ndbc_erai_erai_30a_net_flux_40_h_decimate,[],[],'Barotropic V','40hLP Q_0',[],[],true)
scatter_fit_ts(stn.adcp_v_40_h_decimate,stn.adseatemp_ndbc_erai_erai_30a_net_flux_40_h_decimate,[],[],'Barotropic V','40hLP Q_0',[],'resid',true)
fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date,stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.profile); clabel(cs,ch);
stn = get_looe1_adcp(stn);
fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date,stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.profile); clabel(cs,ch);
stn.adcp_baroclinic_u_40hlp
fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date,stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.prof); clabel(cs,ch);
size(stn.adcp_bin_heights)
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.prof(1:strd:end,:)'); clabel(cs,ch);
datetick3;
colorbar;
stn = verify_variable(stn,'adseatemp_ndbc_erai_erai_30a_net_flux_40_h_decimate');
stn = verify_variable(stn,'adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_decimate');
fmg; plot_ts(stn.adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_decimate);
stn = verify_variable(stn,'adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_sum');
fmg; plot_ts(stn.adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_decimate);
fmg; plot_ts(stn.adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_sum);
plot(stn.adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_sum.date,stn.adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_sum.data+25,'b-');
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_v_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_baroclinic_v_40hlp.prof(1:strd:end,:)'); clabel(cs,ch); colorbar; datetick3;
plot(stn.adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_sum.date,stn.adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_sum.data+25,'b-');
titlename('V vs. Q_0 40hlp');
titlename('U vs. Q_0 40hlp');
plot_ts(stn.adcp_seatemp,'r');
ylim([0,34])
plot(stn.adseatemp_ndbc_erai_erai_30a_wind_stress.date,stn.adseatemp_ndbc_erai_erai_30a_wind_stress.data+30,'k-');
plot_ts(stn.ndbc_wind1_u_72_h_lp,'k',stn.ndbc_wind1_v_72_h_lp,'g');
plot_ts(stn.ndbc_wind1_speed_72_h_lp,'k');
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_btm_x,[],[],'u_H_C','u_B_T_M')
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_btm_x_40_h_lp,[],[],'u_H_C','u_B_T_M')
stn=verify_variable(stn,{'adcp_btm_x_40_h_lp','adcp_btm_l_40_h_lp','adcp_btm_u_40_h_lp','adcp_btm_v_40_h_lp'});
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_btm_x_40_h_lp,[],[],'u_H_C','u_B_T_M')
stn=verify_variable(stn,{'adcp_btm_x_40_h_dec','adcp_btm_l_40_h_dec','adcp_btm_u_40_h_dec','adcp_btm_v_40_h_dec'});
stn=verify_variable(stn,{'adcp_btm_x_40_h_decimate','adcp_btm_l_40_h_decimate','adcp_btm_u_40_h_decimate','adcp_btm_v_40_h_decimate'});
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_btm_x_40_h_lp,[],[],'u_H_C','u_B_T_M')
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_btm_x_40_h_decimate,[],[],'u_H_C','u_B_T_M')
stn.adcp_baroclinic_btm_x
stn=verify_variable(stn,{'adcp_baroclinic_btm_x_40_h_decimate','adcp_baroclinic_btm_l_40_h_decimate','adcp_baroclinic_btm_u_40_h_decimate','adcp_baroclinic_btm_v_40_h_decimate'});
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_baroclinic_btm_x_40_h_decimate,[],[],'u_H_C','u_B_T_M')
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_baroclinic_btm_x,[],[],'u_H_C','u_B_T_M')
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_baroclinic_btm_x_40_h_decimate,[],[],'u_H_C','u_B_T_M')
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_baroclinic_btm_l_40_h_decimate,[],[],'u_H_C','u_B_T_M')
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_baroclinic_btm_l_40_h_decimate,[],[],'u_H_C','v_B_T_M')
nansummary(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u.data)
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_baroclinic_btm_l_40_h_decimate,[],[],'u_H_C','u_B_T_M')
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_baroclinic_btm_l_40_h_decimate,[],[],'u_H_C','v_B_T_M')
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,stn.adcp_baroclinic_btm_x_40_h_decimate,[],[],'u_H_C','u_B_T_M')
stn = get_looe1_microcat(stn);
stn.mc_seatemp = interp_ts(stn.microcat_seatemp);
stn.dTdz = ts_op(stn.mc_seatemp,stn.ad_seatemp,'-');
fmg; plot_ts(stn.dTdz);
fmg; boxplot_ts(stn.dTdz);
stn.dTdz = ts_op(ts_op(stn.mc_seatemp,stn.ad_seatemp,'-'),17,'./');
fmg; boxplot_ts(stn.dTdz);
stn.dTdt_hc = ts_op(stn.adcp_baroclinic_btm_x,stn.dTdz,'.*
stn.dTdt_hc = ts_op(stn.adcp_baroclinic_btm_x,stn.dTdz,'.*');
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,ts_fun(stn.dTdt_hc,@abs),[],[],'HC','In situ'),
stn.dTdt_hc = ts_op(3600,ts_op(stn.adcp_baroclinic_btm_x,stn.dTdz,'.*'),'.*');
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_u,ts_fun(stn.dTdt_hc,@abs),[],[],'HC','In situ'),
scatter_fit_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc,ts_fun(stn.dTdt_hc,@abs),[],[],'HC','In situ'),
fmg; plot_ts(stn.dTdt_hc);
fmg; grpplot_ts(stn.dTdt_hc);
fmg; grpplot_ts(stn.dTdt_hc); grpplot_ts(stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc,[],[],[],'color','r');
stn.dTdt_is=stn.dTdt_hc;
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_baroclinic_btm_x,stn.dTdz,'.*'),'.*');
stn.dTdt_hc = ts_op(3600,stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc,'.*');
fmg; grpplot_ts(stn.dTdt_is); grpplot_ts(stn.dTdt_hc,[],[],[],'color','r');
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_btm_x,stn.dTdz,'.*'),'.*');
fmg; grpplot_ts(stn.dTdt_is); grpplot_ts(stn.dTdt_hc,[],[],[],'color','r');
fmg; grpplot_ts(stn.dTdt_is,'weekly'); grpplot_ts(stn.dTdt_hc,'weekly',[],[],'color','r');
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_btm_x,stn.dTdz,'.*'),'.*'); nansummary(stn.dTdt_is.data),
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_baroclinic_btm_x,stn.dTdz,'.*'),'.*'); nansummary(stn.dTdt_is.data),
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_btm_x,stn.dTdz,'.*'),'.*'); nansummary(stn.dTdt_is.data),
scatter_fit_ts(ts_fun(stn.dTdt_hc,@abs),ts_fun(stn.dTdt_is,@abs),[],[],'HC','In situ'),
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_btm_x_40_h_decimate,stn.dTdz_40_h_decimate,'.*'),'.*'); nansummary(stn.dTdt_is.data),
stn=verify_variable(stn,'dTdz_40_h_decimate');
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_btm_x_40_h_decimate,stn.dTdz_40_h_decimate,'.*'),'.*'); nansummary(stn.dTdt_is.data),
scatter_fit_ts(ts_fun(stn.dTdt_hc,@abs),ts_fun(stn.dTdt_is,@abs),[],[],'HC','In situ'),
stn.dTdt_hc = ts_op(3600,stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc_40_h_decimate,'.*');
stn=verify_variable(stn,'adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc_40_h_decimate');
stn.dTdt_hc = ts_op(3600,stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc_40_h_decimate,'.*');
scatter_fit_ts(ts_fun(stn.dTdt_hc,@abs),ts_fun(stn.dTdt_is,@abs),[],[],'HC','In situ'),
stn=verify_variable(stn,'adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc_1_d_sum');
stn=verify_variable(stn,'dTdt_is_1_d_sum');
stn.dTdt_hc = ts_op(3600,stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc,'.*');
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_btm_x_40_h_decimate,stn.dTdz,'.*'),'.*'); nansummary(stn.dTdt_is.data),
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_btm_x,stn.dTdz,'.*'),'.*'); nansummary(stn.dTdt_is.data),
stn=verify_variable(stn,{'dTdt_is_1_d_sum','dTdt_hc_1_d_sum'});
scatter_fit_ts(ts_fun(stn.dTdt_hc_1_d_sum,@abs),ts_fun(stn.dTdt_is_1_d_sum,@abs),[],[],'HC','In situ'),
stn.dTdt_hc = stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc;
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_btm_x,stn.dTdz,'.*'),'.*'); nansummary(stn.dTdt_is.data),
stn.dTdt_hc = stn.adseatemp_ndbc_erai_erai_30a_erai_avhrr_hc_dTdthc; nansummary(stn.dTdt_hc.data),
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_baroclinic_btm_x,stn.dTdz,'.*'),'.*'); nansummary(stn.dTdt_is.data),
stn=verify_variable(stn,{'adcp_baroclinic_btm_x_40_h_lp','dTdz_40_h_lp'});
stn.dTdt_is = ts_op(3600,ts_op(stn.adcp_baroclinic_btm_x,stn.dTdz,'.*'),'.*'); nansummary(stn.dTdt_is.data),
stn=verify_variable(stn,{'dTdt_is_40_h_lp'});
help verify_variable
stn=verify_variable(stn,{'dTdt_is_1_d_sum','dTdt_hc_1_d_sum'},true);
scatter_fit_ts(ts_fun(stn.dTdt_hc_1_d_sum,@abs),ts_fun(stn.dTdt_is_1_d_sum,@abs),[],[],'HC','In situ'),
[stn.dTdt_is_1d.data,stn.dTdt_is_1d.date]=grp_ts(stn.dTdt_is.data,stn.dTdt_is.date,@floor,@nansum);
[stn.dTdt_hc_1d.data,stn.dTdt_hc_1d.date]=grp_ts(stn.dTdt_hc.data,stn.dTdt_hc.date,@floor,@nansum);
scatter_fit_ts(ts_fun(stn.dTdt_hc_1d,@abs),ts_fun(stn.dTdt_is_1d,@abs),[],[],'HC','In situ'),
[stn.dTdt_hc_5d.data,stn.dTdt_hc_5d.date]=grp_ts(stn.dTdt_hc.data,stn.dTdt_hc.date,@get_yearpentad,@nansum);
[stn.dTdt_is_5d.data,stn.dTdt_is_5d.date]=grp_ts(stn.dTdt_is.data,stn.dTdt_is.date,@get_yearpentad,@nansum);
scatter_fit_ts(ts_fun(stn.dTdt_hc_5d,@abs),ts_fun(stn.dTdt_is_5d,@abs),[],[],'HC','In situ'),
stn.dTdt_isx = ts_op(3600,ts_op(stn.adcp_baroclinic_btm_l,stn.dTdz,'.*'),'.*'); nansummary(stn.dTdt_isx.data),
[stn.dTdt_isx_5d.data,stn.dTdt_isx_5d.date]=grp_ts(stn.dTdt_isx.data,stn.dTdt_isx.date,@get_yearpentad,@nansum);
scatter_fit_ts(ts_fun(stn.dTdt_hc_5d,@abs),ts_fun(stn.dTdt_isx_5d,@abs),[],[],'HC','In situ CONTROL'),
nansummary(stn.adcp_baroclinic_x_3hlp.data)
nansummary(stn.adcp_x_3hlp.data)
nansummary(stn.adcp_x.data),nansummary(stn.adcp_x_3hlp.data),nansummary(stn.adcp_x_40hlp.data)
stn.dTdt_is_40hlp = ts_op(3600,ts_op(stn.adcp_baroclinic_btm_x_40_h_lp,stn.dTdz_40_h_lp,'.*'),'.*'); nansummary(stn.dTdt_is.data),
[stn.dTdt_is_40hlp_5d.data,stn.dTdt_is_40hlp_5d.date]=grp_ts(stn.dTdt_is.data,stn.dTdt_is.date,@get_yearpentad,@nansum);
[stn.dTdt_is_40hlp_5d.data,stn.dTdt_is_40hlp_5d.date]=grp_ts(stn.dTdt_is_40hlp.data,stn.dTdt_is_40hlp.date,@get_yearpentad,@nansum);
scatter_fit_ts(ts_fun(stn.dTdt_hc_5d,@abs),ts_fun(stn.dTdt_is_40hlp_5d,@abs),[],[],'HC','In situ'),
scatter_fit_ts(ts_fun(stn.dTdt_hc_5d,@abs),ts_fun(stn.dTdt_is_40hlp_5d,@abs),[],[],'\Sigma_5_dHC','In situ 40hlp'),
scatter_fit_ts(stn.dTdt_hc_5d,stn.dTdt_is_40hlp_5d,[],[],'\Sigma_5_dHC','In situ 40hlp'),
scatter_fit_ts(stn.dTdt_hc_5d,stn.dTdt_is_5d,[],[],'\Sigma_5_dHC','In situ 40hlp'),
scatter_fit_ts(stn.dTdt_hc_5d,stn.dTdt_is_5d,[],[],'\Sigma_5_dHC','In situ'),
grepstruct(stn,'_u')
fmg; plot_spec(stn,'adcp_baroclinic_btm_v');
fmg; plot_spec(stn,'adcp_baroclinic_btm_x');
fmg; plot_spec(stn,'adcp_baroclinic_sfc_x');
help plot_spec
fmg; plot_spec(stn,'adcp_sfc_x');
plot_spec(stn,'adcp_sfc_x');
plot_spec(stn,'adcp_btm_x');
plot_spec(stn,'adcp_baroclinic_btm_x');
plot_spec(stn,'adcp_x');
plot_spec(stn,'adcp_baroclinic_x');
for fld={'adcp_baroclinic_btm_x','adcp_btm_x','adcp_baroclinic_sfc_x','adcp_sfc_x','adcp_x'}; plot_spec(stn,fld{:}); end;
for fld={'adcp_baroclinic_btm_l','adcp_btm_l','adcp_baroclinic_sfc_l','adcp_sfc_l','adcp_l'}; plot_spec(stn,fld{:}); end;
12*40/24
stn=verify_variable(stn,{'adcp_baroclinic_btm_x_20_h_lp','dTdz_20_h_lp'});
stn.dTdt_is_20hlp = ts_op(3600,ts_op(stn.adcp_baroclinic_btm_x_20_h_lp,stn.dTdz_20_h_lp,'.*'),'.*'); nansummary(stn.dTdt_is_20hlp.data),
fmg; hist(stn.dTdt_is_20hlp.data,1000),
fmg; hist(stn.dTdt_is_20hlp.data*5*24,1000),
[stn.dTdt_is_20hlp_5d.data,stn.dTdt_is_20hlp_5d.date]=grp_ts(stn.dTdt_is_20hlp.data,stn.dTdt_is_20hlp.date,@get_yearpentad,@nansum);
scatter_fit_ts(stn.dTdt_hc_5d,stn.dTdt_is_20hlp_5d,[],[],'\Sigma_5_dHC','In situ 20hlp'),
edit xanlooe_hc_is.m
scatter_fit_ts(stn.adcp_btm_x,stn.adcp_baroclinic_btm_x),
scatter_fit_ts(stn.adcp_btm_x,stn.adcp_baroclinic_btm_x,[],[],'u_B_T_M^x^s','u_B_T_M^x^s-\integral_0^zu^x^sdz'),
scatter_fit_ts(stn.adcp_btm_x,stn.adcp_baroclinic_btm_x,[],[],'u_B_T_M^x^s','u_B_T_M^x^s-\int_0^zu^x^sdz'),
scatter_fit_ts(stn.adcp_btm_x,stn.adcp_baroclinic_btm_x,[],[],'u_B_T_M^x^s','u_B_T_M^x^s - \int_0^zu^x^sdz'),
fmg; plot_ts(stn.adcp_btm_x,stn.adcp_baroclinic_btm_x);
fmg; plot_ts(stn.adcp_btm_x,stn.adcp_baroclinic_btm_x); legend('u_B_T_M^x^s','u_B_T_M^x^s - \int_0^zu^x^sdz');
datetick3;
fmg; plot_ts(stn.adcp_btm_x,stn.adcp_baroclinic_btm_x,stn.adcp_x); legend('u_B_T_M^x^s','u_B_T_M^x^s - \int_0^zu^x^sdz','u^x^s');
datetick3;
fmg; plot_ts(stn.adcp_btm_x,stn.adcp_sfc_x); legend('u_B_T_M^x^s','u_S_F_C^x^s');
datetick3;
fmg; hist(stn.adcp_btm_x.data,1000);
fmg; hist(stn.adcp_btm_x.data,1000); titlename('u_B_T_M^X^S all seasons');
fmg; hist(stn.adcp_btm_x.data(ts_jas(stn.adcp_btm_x)),1000); titlename('u_B_T_M^X^S all seasons');
fmg; hist(stn.adcp_btm_x.data(ts_jas(stn.adcp_btm_x)),1000); titlename('u_B_T_M^X^S JAS');
fmg; hist(stn.adcp_btm_x.data(ts_jfm(stn.adcp_btm_x)),1000); titlename('u_B_T_M^X^S JFM');
stn.adcp_bin_heights
22-5
stn.adcp_bin_heights(22)
stn.adcp_bin_heights(21)
fmg; plot(stn.adcp_x.date,stn.adcp_x.prof([1,2])); datetick3;
stn.adcp_x
fmg; plot(stn.adcp_x.date,stn.adcp_x.prof(:,[1,2])); datetick3;
scatter_fit(stn.adcp_x.prof(:,1),stn.adcp_x.prof(:,2)),
help scatter_fit
scatter_fit(stn.adcp_x.prof(:,1),stn.adcp_x.prof(:,2),'Bin1','Bin2'),
scatter_fit(stn.adcp_x.prof(:,1),stn.adcp_x.prof(:,2),'Bin1','Bin2',[],[],[],true),
stn.hc_is = ts_op( ts_op(stn.adcp_bin2_x,stn.ad_seatemp,'.*') , ts_op(stn.adcp_bin21_x,stn.mc_seatemp,'.*') , '-' );
stn.adcp_bin2_x.date=stn.adcp_x.date; stn.adcp_bin2_x.date=stn.adcp_x.data(:,2);
stn.adcp_bin2_x.date=stn.adcp_x.date; stn.adcp_bin2_x.date=stn.adcp_x.prof(:,2);
stn.adcp_bin21_x.date=stn.adcp_x.date; stn.adcp_bin21_x.date=stn.adcp_x.prof(:,21);
stn.adcp_bin_heights
stn.hc_is = ts_op( ts_op(stn.adcp_bin2_x,stn.ad_seatemp,'.*') , ts_op(stn.adcp_bin21_x,stn.mc_seatemp,'.*') , '-' );
stn.adcp_bin2_x.date=stn.adcp_x.date; stn.adcp_bin2_x.data=stn.adcp_x.prof(:,2);
stn.adcp_bin21_x.date=stn.adcp_x.date; stn.adcp_bin21_x.data=stn.adcp_x.prof(:,21);
stn.hc_is = ts_op( ts_op(stn.adcp_bin2_x,stn.ad_seatemp,'.*') , ts_op(stn.adcp_bin21_x,stn.mc_seatemp,'.*') , '-' );
scatter_fit_ts(stn.dTdt_hc,stn.hc_is,[],[],'HC','u_3_mT_3_m - u_1_7_mT_1_7_m'),
nansummary(stn.dTdt_hc.data)
nansummary(stn.dTdt_hc.data*24)
stn.station_name
stn=[]; ans=[]; clear stn ans; pack
clear all
pack
stn = optimize_station_heat_budget('smkf1','erai','avhrr_weekly','ndbc','tpxo_tide','erai');
looe1 = get_looe1_microcat;
looe1.mc_seatemp = interp_ts(looe1.microcat_seatemp);
numel(getfield(intersect_tses(looe1.mc_seatemp,stn.ndbc_erai_erai_30a_erai_avhrr_hc_dTdt),'data'))
v=intersect_tses(looe1.mc_seatemp,stn.ndbc_erai_erai_30a_erai_avhrr_hc_dTdt); v,
v=intersect_tses(looe1.mc_seatemp,stn.ndbc_erai_erai_30a_erai_avhrr_hc_dTdt); v(1)
v=intersect_tses(looe1.mc_seatemp,stn.ndbc_erai_erai_30a_erai_avhrr_hc_dTdt); v{1}
v=intersect_tses(looe1.mc_seatemp,stn.ndbc_erai_erai_30a_erai_avhrr_hc_dTdt); v{2}
reviewanim([],0,1)
reviewanim([],1,0)
print('-dtiff','../figs/looe1-hc-scatter-baroclinic_btm_x-BEST.tiff');
print('-dtiff','../figs/looe1-hc-scatter-baroclinic_btm_x_40hlp.tiff');
print('-dtiff','../figs/looe1-hc_5d-scatter-baroclinic_btm_x_20hlp.tiff');
print('-dtiff','../figs/looe1-hc_5d-scatter-baroclinic_btm_x_40hlp.tiff');
delete('../figs/looe1-hc-scatter-baroclinic_btm_x_40hlp.tiff');
print('-dtiff','../figs/looe1-scatter-bin2T-minus-bin21T.tiff');
fmg; plot_ts(stn.ndbc_wind1_speed_72_h_lp,'k');
