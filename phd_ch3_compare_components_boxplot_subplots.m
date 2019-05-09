function phd_ch3_compare_components_boxplot_subplots(stn)
%function phd_ch3_compare_components_boxplot_subplots(stn)
%
% Produce a 2x3 matrix of boxplots comparing components of the total heat budget
%
% Last Saved Time-stamp: <Sun 2013-07-21 16:00:45 Eastern Daylight Time gramer>

  doPrint = false;

  [t,hc,hchc,dt,bq,q,qlh] = intersect_tses([],stn.ndbc_sea_t_dly_diff, ...
                                           stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly, ...
                                           stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc_dly, ...
                                           stn.b_ndbc_erai_erai_30a_avhrr_dt_dly, ...
                                           stn.b_ndbc_erai_erai_30a_net_flux_term_dly, ...
                                           stn.simple_ndbc_erai_erai_30a_net_flux_term_dly, ...
                                           stn.ndbc_erai_erai_30a_latent_flux_term_dly); 

  mn=-2.5;
  mx=+2.5;

  fmg;
  n=2;
  m=3;
  ix=0;

  ix=ix+1;
  subplot(n,m,ix);
  boxplot_ts(t);
  ylim([mn,mx]);
  grid on;
  ylabel([]);
  set(gca,'xticklabel',[]);
  xlabel('\Delta_1_dT_s');

  ix=ix+1;
  subplot(n,m,ix);
  boxplot_ts(hc);
  ylim([mn,mx]);
  grid on;
  ylabel([]);
  set(gca,'xticklabel',[]);
  xlabel('\Sigma_1_d \partial_tT_s');

  ix=ix+1;
  subplot(n,m,ix);
  boxplot_ts(dt);
  ylim([mn,mx]);
  grid on;
  ylabel([]);
  set(gca,'xticklabel',[]);
  xlabel('\Sigma_1_d u\bullet\nablaT_s+K\nabla^2T_s+Q/\rhoC_ph');

  ix=ix+1;
  subplot(n,m,ix);
  boxplot_ts(bq);
  ylim([mn,mx]);
  grid on;
  ylabel([]);
  xlabel('\Sigma_1_d(Q_0(\gamma)+Q_b)/\rhoC_ph');

  ix=ix+1;
  subplot(n,m,ix);
  boxplot_ts(q);
  ylim([mn,mx]);
  grid on;
  ylabel([]);
  xlabel('\Sigma_1_dQ_0/\rhoC_ph');

  ix=ix+1;
  subplot(n,m,ix);
  boxplot_ts(qlh);
  ylim([mn,mx]);
  grid on;
  ylabel([]);
  xlabel('\Sigma_1_dQ_L_H/\rhoC_ph');

  set(gcf,'Name',upper(stn.station_name));

  if ( doPrint )
    print('-dtiff',fullfile(get_thesis_path('../figs'), ...
                            [lower(stn.station_name),'-',mfilename,'.tif']));
  end;

  fmg;
  subplot(1,2,1);
  boxplot_ts(dt);
  ylim([mn,mx]);
  grid on;
  ylabel([]);
  set(gca,'xticklabel',[]);
  xlabel('\Sigma_1_d u\bullet\nablaT_s+K\nabla^2T_s+Q/\rhoC_ph');

  subplot(1,2,2);
  boxplot_ts(hchc);
  ylim([mn,mx]);
  grid on;
  ylabel([]);
  set(gca,'xticklabel',[]);
  xlabel('\Sigma_1_d u_H_C\bullet\nablaT_s(h,\beta,Q_0)');

return;
