1;

%{
  looe1.dTdz = ts_op(looe1.microcat_seatemp,looe1.adcp_seatemp,'-');
  % looe1.adcp_btm_xdTdz = ts_op(looe1.adcp_baroclinic_btm_x,looe1.dTdz,'*');
  % looe1.adcp_baroclinic_btm_xdTdz = ts_op(looe1.adcp_baroclinic_btm_x,looe1.dTdz,'*');

  ix=find(datenum(2007,6,1)<=looe1.adcp_baroclinic_x_3hlp.date&looe1.adcp_baroclinic_x_3hlp.date<datenum(2007,10,1));
  ts.date=looe1.adcp_baroclinic_x_3hlp.date(ix); ts.prof=looe1.adcp_baroclinic_x_3hlp.prof(ix,:);

  fmg;
  contourf(ts.date,looe1.adcp_bin_heights,ts.prof');
  xlim(datenum(2007,[6,8],[15,22])); ylim([-1,22]); caxis([-.5,.5]); datetick3;
  colorbar; set_hovmuller_cursor;
  plot_ts(looe1.dTdz,'b.-');
  plot_ts(ts_op(stn.simple_ndbc_erai_erai_30a_net_flux_term,10,'*'),'r');
  plot_ts(looe1.adcp_seatemp,'k-','linew',2);
  plot_ts(looe1.microcat_seatemp,'m');
  plot(looe1.adcp_x.date,(looe1.adcp_x.prof(:,2)*10)+25,'r');
  xlim(datenum(2007,[7,8],[20,1])); ylim([-1,32]); caxis([-.2,.2]); datetick3;
  legend('BC 3hlp Cross','dT/dz','Q0/\rhoC_ph','T_b_t_m','T_t_o_p','X_b_t_m', 'Location','East');

  scatter_fit_ts(looe1.dTdz,stn.b_ndbc_erai_erai_30a_avhrr_dt_dly,@ts_jas,[],'LOOE1 dT/dz','MLRF1 Q_0+Q_b');
%}

if (0)
  %dts=datenum(2007,[7,9],1); hrs=[15:23];
  %% CONSISTENT RESULT?! Afternoon/evening deep heat transport (1.4m deep) matches horizontal convection
  dts=datenum(2007,[7,8],1); hrs=[17:23];
  fs = @(x)(find(ismember(get_hour(x.date),hrs)&dts(1)<=x.date&x.date<dts(end)));
  ht=klgf1.cm_deep_xdTdz; ix=fs(ht); ht.date=ht.date(ix); ht.data=ht.data(ix)*1.4; nansum(ht.data)
  hc=stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc; ix=fs(hc); hc.date=hc.date(ix); hc.data=hc.data(ix); nansum(hc.data)

  scatter_fit_ts(hc,ht,[],[],'Heat Transport','Horizontal Convection'),
end;

if (0)
  fmg; plot_ts(looe1.adcp_btm_x,looe1.adcp_x); legend('BTM XS','XS'); xlim(datenum(2007,[7,8],[20,1])); datetick3;

  perfun = @floor; minN = 5; ylm=[-1.1,+1.1];
  %perfun = @get_yeartriad; minN = 5*3; ylm=[-2,+2];
  %perfun = @get_yearpentad; minN = 5*4; ylm=[-3,+3];
  %perfun = @get_yearweek; minN = 5*6; ylm=[-4,+4];

  fmg;
  %grpplot_ts(ts_op(0,klgf1.cm_deep_xdTdz,'-'),perfun,@nansum,minN,'LineStyle','none','Marker','.','Color','b');
  %klgf1.cum_heat_transport=grpplot_ts(klgf1.cm_deep_xdTdz_10min,perfun,@nansum,minN,'LineStyle','none','Marker','.','Color','b');
  cum_ht=grpplot_ts(ht,perfun,@nansum,minN,'LineStyle','none','Marker','.','Color','b');
  cum_hc=grpplot_ts(hc,perfun,@nansum,minN,'LineStyle','none','Marker','.','Color','r');
  xlim(ht.date([1,end])); datetick3;
  ylim(ylm);
  legend(['\Sigmau_2_2(T_7 - T_2_2) ',num2str(hrs,'%d,'),'h UT'],['\SigmaQ_vdT_H_C/dx ',strrep(char(perfun),'_','\_')],'Location','NorthWest');
  titlename([upper(stn.station_name),' PM horizontal convection vs. ',upper(klgf1.station_name),' heat transport']);

  scatter_fit_ts(cum_hc,cum_ht);


  [all_ht,all_hc]=intersect_tses(klgf1.cm_deep_xdTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
  fmg;
  spt(2,1,1); boxplot_ts(subset_ts(all_ht,@(x)(find(ismember(get_hour(x.date),hrs)))),@(d)(datestr(d,3))); ylim([-.1,+.1]); xlabel('Heat Transport'); grid on; ylabel('^oC/hour');
  spt(2,1,2); boxplot_ts(subset_ts(all_hc,@(x)(find(ismember(get_hour(x.date),hrs)))),@(d)(datestr(d,3))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
  titlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(klgf1.station_name),' heat transport climatology']);

  [all_ht,all_hc]=intersect_tses(klgf1.cm_deep_xdTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
  fmg;
  spt(2,1,1); boxplot_ts(subset_ts(all_ht,@(x)(find(ismember(get_hour(x.date),hrs)))),@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Heat Transport'); grid on; ylabel('^oC/hour');
  spt(2,1,2); boxplot_ts(subset_ts(all_hc,@(x)(find(ismember(get_hour(x.date),hrs)))),@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
  titlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(klgf1.station_name),' heat transport']);

  % [all_ht,all_hc]=intersect_tses(klgf1.cm_dxTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
  % fmg;
  % spt(2,1,1); boxplot_ts(subset_ts(all_ht,@(x)(find(ismember(get_hour(x.date),hrs))))); ylim([-.1,+.1]); xlabel('Sheered Heat Transport'); grid on;
  % spt(2,1,2); boxplot_ts(subset_ts(all_hc,@(x)(find(ismember(get_hour(x.date),hrs))))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on;

  [all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_btm_xdTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Baroclinic Heat Transport'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7);
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7);
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport']);

  [all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_x1dTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Baroclinic Bottom Heat Transport'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7);
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7);
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport']);

  [all_ht,all_hc]=intersect_tses(looe1.adcp_btm_xdTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Heat Transport'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7);
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7);
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport']);


  allmos={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

  [all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_btm_xdTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); xlabel('Heat Transport'); grid on; ylabel('^oC/hour');
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport climatology']);

  [all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_x1dTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); xlabel('Bottom Heat Transport'); grid on; ylabel('^oC/hour');
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport climatology']);


  looe1 = verify_variable(looe1,'adcp_baroclinic_x1dTdz_1_d_avg');
  stn = verify_variable(stn,'ndbc_erai_erai_30a_avhrr_hc_dTdthc_1_d_sum');
  [all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_x1dTdz_1_d_sum,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc_1_d_sum);
  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); xlabel('Bottom Heat Transport'); grid on; ylabel('^oC/day');
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/day');
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport daily sums']);



  fs = @(x)(find(ismember(get_hour(x.date),[15:23])|~ismember(get_jday(x.date),[105:248])));

  [all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_btm_xdTdz,subset_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc,fs));
  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); xlabel('Baroclinic Heat Transport (summer eve)'); grid on; ylabel('^oC/hour');
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport climatology']);

  [all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_btm_xdTdz,subset_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc,fs));
  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Baroclinic Heat Transport (summer eve)'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7);
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7);
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport']);


  fmg; plot_ts(looe1.adcp_seatemp,looe1.microcat_seatemp); %datetick3('x',12); 
  %plot_ts(ts_op(looe1.adcp_speed,16,'+'),'r');
  plot_ts(ts_op(ts_op(looe1.adcp_x,50,'*'),15,'+'),'r');
  plot_ts(ts_op(stn.ndbc_wind1_xshore,15,'+'),'k');
  % plot_ts(ts_op(ts_op(looe1.adcp_speed,50,'*'),0,'+'),'r');
  % plot_ts(ts_op(stn.ndbc_wind1_speed,0,'+'),'k');
  plot_ts(ts_op(ts_op(stn.simple_ndbc_erai_erai_30a_net_flux,100,'/'),10,'+'),'m');
  plot_ts(stn.ndbc_sea_t,'k','LineWidth',1.5);
  xlim(datenum(2009,[6,7],[25,15])); datetick3; ylim([0,40]); legend('T_B_T_M','T_T_O_P','u*50+15','U+15','Q0/100+10','T_s');
  titlename('LOOE1 rapid summer bottom cooling');


  fs = @(x)(find(ismember(get_hour(x.date),[15:23])));
  [all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_btm_xdTdz,subset_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc,fs));

  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); title('Baroclinic Heat Transport (15-23h)'); grid on; ylabel('^oC/hour');
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,3)),'grouporder',allmos); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport climatology']);

  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,12))); ylim([-.1,+.1]); title('Baroclinic Heat Transport (15-23h)'); grid on; ylabel('^oC/hour');
  hs=findobj(gca,'type','text'); for h=hs(:)'; set(h,'FontSize',6); end;
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,12))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
  hs=findobj(gca,'type','text'); for h=hs(:)'; set(h,'FontSize',6); end;
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport']);

end; %if (0)


if (0)
  allmos={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
  fs = @(x)(find(ismember(get_hour(x.date),[15:23])));
  [all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_btm_xdTdz_3_h_lp,subset_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc,fs));

  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,3)),'grouporder',allmos,'mean',true); ylim([-.1,+.1]); title('3hlp Baroclinic Heat Transport (15-23h)'); grid on; ylabel('^oC/hour');
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,3)),'grouporder',allmos,'mean',true); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport climatology']);

  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,12)),'mean',true); ylim([-.1,+.1]); title('3hlp Baroclinic Heat Transport (15-23h)'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7.5);
  hs=findobj(gca,'type','text'); for h=hs(:)'; set(h,'FontSize',6); end;
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,12)),'mean',true); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7.5);
  hs=findobj(gca,'type','text'); for h=hs(:)'; set(h,'FontSize',6); end;
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport']);
end; %if (0)


if (1)
  allmos={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
  %hrs=[0:23];
  hrs=[15:23];
  fs = @(x)(find(ismember(get_hour(x.date),hrs)));
    % lags = [-6,-5,-4,-3,-2,-1,0,-1,-2,-3,-4,-5]';
    % lag = lags(get_month(looe1.adcp_mcd_sheer_xdTdz_3_h_lp.date));
    %% Minus six to zero
    %lag = -3*cos(2*pi*get_yearday(looe1.adcp_mcd_sheer_xdTdz_3_h_lp.date)/365)-3;
    % Minus six to plus six
    lag = -6*cos(2*pi*get_yearday(looe1.adcp_mcd_sheer_xdTdz_3_h_lp.date)/365);

    ht.date = looe1.adcp_mcd_sheer_xdTdz_3_h_lp.date + (lag/24);
    ht.data = looe1.adcp_mcd_sheer_xdTdz_3_h_lp.data;
    [all_ht,all_hc]=intersect_tses(ht,subset_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc,fs));

    fmg;
    spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,3)),'grouporder',allmos,'mean',true); ylim([-.1,+.1]); title(['3hlp Sheered Heat Transport var. lag']); grid on; ylabel('^oC/hour');
    spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,3)),'grouporder',allmos,'mean',true); ylim([-.1,+.1]); xlabel(['Horizontal Convection']); grid on; ylabel('^oC/hour');
    suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport climatology ',num2str(numel(hrs)),' hrs']);
    print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-hc-vs-looe1-ht-var-lag.tiff']));

  if (0)
  %for lag=0;
  for lag=[-6,0,6];
  %for lag=[0,3];
  %for lag=-6:6;
    ht.date = looe1.adcp_mcd_sheer_xdTdz_3_h_lp.date + (lag/24);
    ht.data = looe1.adcp_mcd_sheer_xdTdz_3_h_lp.data;
    [all_ht,all_hc]=intersect_tses(ht,subset_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc,fs));

    fmg;
    spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,3)),'grouporder',allmos,'mean',true); ylim([-.1,+.1]); title(['3hlp Sheered Heat Transport ',num2str(lag),'hr lag']); grid on; ylabel('^oC/hour');
    spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,3)),'grouporder',allmos,'mean',true); ylim([-.1,+.1]); xlabel(['Horizontal Convection']); grid on; ylabel('^oC/hour');
    suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport climatology ',num2str(numel(hrs)),' hrs']);
    print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-hc-vs-looe1-ht-',num2str(lag),'hr.tiff']));
  end;
  end;

  if (0)
    fmg;
    spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,12)),'mean',true); ylim([-.1,+.1]); title('3hlp Sheered Heat Transport'); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7.5);
    hs=findobj(gca,'type','text'); for h=hs(:)'; set(h,'FontSize',6); end;
    spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,12)),'mean',true); ylim([-.1,+.1]); xlabel(['Horizontal Convection ',num2str(lag)','hr lag']); grid on; ylabel('^oC/hour'); set(gca,'FontSize',7.5);
    hs=findobj(gca,'type','text'); for h=hs(:)'; set(h,'FontSize',6); end;
    suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport ',num2str(numel(hrs)),' hrs']);
  end;

  %hrs=[0:23];
  hrs=[15:23];
  % April and December have very limited data
  fs = @(x)(find(ismember(get_hour(x.date),hrs)&ismember(get_month(x.date),[5:11])));
  %lag=0;
  lag=3;
  ht.date = klgf1.cm_deep_xdTdz.date + (lag/24);
  ht.data = klgf1.cm_deep_xdTdz.data;
  dp.date = klgf1.cm_deep_seapres_dt.date + (lag/24);
  dp.data = klgf1.cm_deep_seapres_dt.data;
  [all_ht,all_hc,all_dp]=intersect_tses(ht,subset_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc,fs),dp);

  fmg;
  spt(2,1,1); boxplot_ts(all_ht,@(d)(datestr(d,3)),'mean',true); ylim([-.1,+.1]); title(['Sheered Heat Transport ',num2str(lag),'hr lag']); grid on; ylabel('^oC/hour');
  spt(2,1,2); boxplot_ts(all_hc,@(d)(datestr(d,3)),'mean',true); ylim([-.1,+.1]); xlabel(['Horizontal Convection']); grid on; ylabel('^oC/hour');
  suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(klgf1.station_name),' heat transport climatology ',num2str(numel(hrs)),' hrs']);
  print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-hc-vs-klgf1-ht-',num2str(lag),'hr.tiff']));

end; %if (1)


clear ans dix dts dy dys fs hc hrs hs h ht ix minN mo perfun t0 ts ylm yr
