1;

% Optional full functionality, if you first run script PHD_CH3_NCORE

if ( ~exist('klgf1','var') )
  klgf1=[]; mrtf1=[]; clear klgf1 mrtf1
  [klgf1,mrtf1]=an_ncore;

  klgf1.cm_dTdz = ts_op(klgf1.cm_shallow_seatemp,klgf1.cm_deep_seatemp,'-');
  klgf1.cm_deep_xdTdz = ts_op(klgf1.cm_deep_xshore,klgf1.cm_dTdz,'*');
  klgf1.cm_dTdz_10min = ts_op(klgf1.cm_shallow_seatemp_10min,klgf1.cm_deep_seatemp_10min,'-',[],[],[],(5.5/60/24));
  klgf1.cm_deep_xdTdz_10min = ts_op(klgf1.cm_deep_xshore_10min,klgf1.cm_dTdz_10min,'*',[],[],[],(5.5/60/24));

  mrtf1.cm_dTdz = ts_op(mrtf1.cm_shallow_seatemp,mrtf1.cm_deep_seatemp,'-');
  mrtf1.cm_deep_xdTdz = ts_op(mrtf1.cm_deep_xshore,mrtf1.cm_dTdz,'*');
  mrtf1.cm_dTdz_10min = ts_op(mrtf1.cm_shallow_seatemp_10min,mrtf1.cm_deep_seatemp_10min,'-',[],[],[],(5.5/60/24));
  mrtf1.cm_deep_xdTdz_10min = ts_op(mrtf1.cm_deep_xshore_10min,mrtf1.cm_dTdz_10min,'*',[],[],[],(5.5/60/24));
end;

if ( ~exist('looe1','var') || ~isfield(looe1,'adcp_seatemp') )
  looe1 = []; clear looe1
  looe1 = get_station_from_station_name('looe1');
  looe1 = get_looe1_microcat(looe1);
  looe1 = get_looe1_adcp(looe1);
  %looe1 = rmfield(looe1,grepstruct(looe1,'_(baroclinic|[uvw])'));
  looe1 = rmfield(looe1,grepstruct(looe1,'_([uvw])'));
end;

looe1.dTdz = ts_op(looe1.microcat_seatemp,looe1.adcp_seatemp,'-');

looe1.adcp_btm_xdTdz = ts_op(looe1.dTdz,looe1.adcp_btm_x,'*');
looe1.adcp_baroclinic_btm_xdTdz = ts_op(looe1.dTdz,looe1.adcp_baroclinic_btm_x,'*');
looe1 = verify_variable(looe1,{'dTdz_3_h_lp','adcp_baroclinic_btm_x_3_h_lp'});
looe1.adcp_baroclinic_btm_xdTdz_3_h_lp = ts_op(looe1.dTdz_3_h_lp,looe1.adcp_baroclinic_btm_x_3_h_lp,'*');

[ix1,ix2] = intersect_dates(looe1.dTdz.date,looe1.adcp_baroclinic_x.date);
% looe1.adcp_baroclinic_x1dTdz.date = looe1.dTdz.date(ix1);
% looe1.adcp_baroclinic_x1dTdz.data = looe1.dTdz.data(ix1) .* looe1.adcp_baroclinic_x.prof(ix2,1);
for ix=[1:6,11:16,20:22];
  fld = ['adcp_baroclinic_x',num2str(ix)];
  looe1.([fld]).date = looe1.dTdz.date;
  looe1.([fld]).data = looe1.adcp_baroclinic_x.prof(:,ix);
  looe1.([fld,'dTdz']).date = looe1.dTdz.date(ix1);
  looe1.([fld,'dTdz']).data = looe1.dTdz.data(ix1) .* looe1.adcp_baroclinic_x.prof(ix2,ix);
end;

% Positive sheer means upslope relative movement of near-bottom water
looe1.adcp_mcd_sheer_x = ts_op(looe1.adcp_mcd_x,looe1.adcp_btm_x,'-');
% Negative value indicates cooler water moving inshore
looe1.adcp_mcd_sheer_xdTdz = ts_op(looe1.dTdz,ts_op(0,looe1.adcp_mcd_sheer_x,'-'),'*');

looe1 = verify_variable(looe1,{'adcp_mcd_sheer_x_3_h_lp'});
looe1.adcp_mcd_sheer_xdTdz_3_h_lp = ts_op(looe1.dTdz_3_h_lp,ts_op(0,looe1.adcp_mcd_sheer_x_3_h_lp,'-'),'*');


if (1)
  fmg;
  for hr=0:23;
    spt(4,6,hr+1);
    boxplot_ts(subset_ts(looe1.adcp_mcd_sheer_xdTdz_3_h_lp,@(x)(find(ismember(get_month(x.date),[6:9])&(get_hour(x.date)==hr)))));
    ylim([-.1,.1]); xlabel(num2str(hr)); grid on;
  end;
  suptitlename('LOOE1 mid-level sheer by hour and month (summer)');
  %print('-dtiff',fullfile(get_thesis_path('../figs'),'looe1-adcp_mcd_sheer_xdTdz_3_h_lp-boxplot-hour-month-summer.tif'));
end;

if (0)

  fmg; grpplot_ts(klgf1.cm_dTdz,@floor,@nanmin);
  grpplot_ts(klgf1.cm_dTdz,@floor,@nanmean,'Color','k');
  grpplot_ts(klgf1.cm_dTdz,@floor,@nanmax,'Color','r');
  if ( exist('stn','var') )
    [stn.t.data,stn.t.date] = grp_ts(stn.ndbc_sea_t.data,stn.ndbc_sea_t.date,@floor,@nanmean,20);
    stn.dt.date = stn.t.date(1:end-1);
    stn.dt.data = diff(stn.t.data);
    stn = filter_gaps(stn,'t','dt',1.5);

    grpplot_ts(stn.simple_ndbc_erai_erai_30a_net_flux_term,@floor,@nansum,'Color','m:');
    %plot_ts(stn.simple_ndbc_erai_erai_30a_net_flux_term_dly,'m.');
    plot_ts(stn.dt,'k.');
    legend('Min daily dT/dz','Mean daily dT/dz','Max daily dT/dz','Q0/\rhoC_ph','\Delta_1_dT_s','Location','North')
  else
    legend('Min daily dT/dz','Mean daily dT/dz','Max daily dT/dz','Location','North')
  end;
  xlim(klgf1.cm_dTdz.date([1,end])); datetick3;
  ylim([-2,6]);
  %legend('Min daily dT/dz','Mean daily dT/dz','Max daily dT/dz','Location','North')
  titlename('NCORE Key Largo daily surface-bottom temperature differences');


  fmg; grpplot_ts(klgf1.cm_dTdz_10min,@floor,@nanmin);
  grpplot_ts(klgf1.cm_dTdz_10min,@floor,@nanmean,'Color','k');
  grpplot_ts(klgf1.cm_dTdz_10min,@floor,@nanmax,'Color','r');
  if ( exist('stn','var') )
    grpplot_ts(stn.simple_ndbc_erai_erai_30a_net_flux_term,@floor,@nansum,'Color','m');
  end;
  xlim(klgf1.cm_dTdz.date([1,end])); datetick3;
  ylim([-2,6]);
  legend('Min daily 10min dT/dz','Mean daily dT/dz','Max daily 10min dT/dz','Location','North')
  titlename('NCORE Key Largo daily surface-bottom 10-minute temp. differences');

  % [sh,dp]=intersect_tses(klgf1.cm_shallow_xshore,klgf1.cm_deep_xshore); z=[-7,-22]; t=[sh.date]; u=[sh.data,dp.data]';
  % fmg; contourf(t,z,u); xlim(datenum(2007,7,[7,25])); datetick3; colorbar; caxis([-0.4,0.4]); set_hovmuller_cursor;


  fmg; grpplot_ts(mrtf1.cm_dTdz,@floor,@nanmin);
  grpplot_ts(mrtf1.cm_dTdz,@floor,@nanmean,'Color','k');
  grpplot_ts(mrtf1.cm_dTdz,@floor,@nanmax,'Color','r');
  if ( exist('stn','var') )
    grpplot_ts(stn.simple_ndbc_erai_erai_30a_net_flux_term,@floor,@nansum,'Color','m:');
    plot_ts(stn.simple_ndbc_erai_erai_30a_net_flux_term_dly,'m.');
    plot_ts(stn.ndbc_sea_t_dly_diff,'k.');
    legend('Min daily dT/dz','Mean daily dT/dz','Max daily dT/dz','Q0/\rhoC_ph','\Delta_1_dT_s','Location','North')
  else
    legend('Min daily dT/dz','Mean daily dT/dz','Max daily dT/dz','Location','North')
  end;
  xlim(mrtf1.cm_dTdz.date([1,end])); datetick3;
  ylim([-2,6]);
  %legend('Min daily dT/dz','Mean daily dT/dz','Max daily dT/dz','Location','North')
  titlename('NCORE Marathon daily surface-bottom temperature differences');

  fmg; grpplot_ts(mrtf1.cm_dTdz_10min,@floor,@nanmin);
  grpplot_ts(mrtf1.cm_dTdz_10min,@floor,@nanmean,'Color','k');
  grpplot_ts(mrtf1.cm_dTdz_10min,@floor,@nanmax,'Color','r');
  if ( exist('stn','var') )
    grpplot_ts(stn.simple_ndbc_erai_erai_30a_net_flux_term,@floor,@nansum,'Color','m');
  end;
  xlim(mrtf1.cm_dTdz.date([1,end])); datetick3;
  ylim([-2,6]);
  legend('Min daily 10min dT/dz','Mean daily 10min dT/dz','Max daily 10min dT/dz','Location','North')
  titlename('NCORE Marathon daily surface-bottom 10-minute temp. differences');

  % [sh,dp]=intersect_tses(mrtf1.cm_shallow_xshore,mrtf1.cm_deep_xshore); z=[-7,-22]; t=[sh.date]; u=[sh.data,dp.data]';
  % fmg; contourf(t,z,u); xlim(datenum(2007,7,[7,25])); datetick3; colorbar; caxis([-0.4,0.4]); set_hovmuller_cursor;



  fmg; grpplot_ts(looe1.dTdz,@floor,@nanmin); grpplot_ts(looe1.dTdz,@floor,@nanmean,'k'); grpplot_ts(looe1.dTdz,@floor,@nanmax,'Color','r');
  legend('Min daily dT/dz','Mean daily dT/dz','Max daily dT/dz', 'Location','North')
  titlename('LOOE1 daily surface-bottom temperature differences');

  if ( exist('stn','var') )
    fmg; plot_ts(ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/'),'m',looe1.microcat_seatemp,'b',looe1.adcp_seatemp,'r');
    xlim(datenum([2005,2010],[1,6],1)); datetick3;
  end;

end; %if (0)

clear ans dix dts dy dys fs hc hrs hs h ht ix minN mo perfun t0 ts ylm yr
