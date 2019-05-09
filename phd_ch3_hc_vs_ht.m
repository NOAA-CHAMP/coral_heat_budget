1;

dt = ts_op(stn.b_ndbc_erai_erai_30a_avhrr_dt,stn.ndbc_sea_t_diff,'-');
ht.data=looe1.adcp_mcd_sheer_xdTdz.data; ht.date=looe1.adcp_mcd_sheer_xdTdz.date;
% ht.data=looe1.adcp_mcd_sheer_xdTdz_3_h_lp.data; ht.date=looe1.adcp_mcd_sheer_xdTdz_3_h_lp.date;
% % dt = ts_op(stn.b_ndbc_erai_erai_30a_avhrr_dt_dly,stn.ndbc_sea_t_dly_diff,'-');
% % [ht.data,ht.date] = grp_ts(looe1.adcp_mcd_sheer_xdTdz.data,looe1.adcp_mcd_sheer_xdTdz.date,@floor,@nansum,24);

[dt,ht]=intersect_tses(dt,ht);

allmos={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

fmg;
boxplot_ts(dt,@(x)(datestr(x,3)),'symbol','k.','grouporder',allmos,'mean','x','allcol',[.6,.6,.6]);
hold on;
%% Works in 2010a but not 2007a: 'boxstyle','filled'
%boxplot_ts(ts_op(-1.4,ht,'*'),[],'boxstyle','filled','symbol','w','allcol','k');
boxplot_ts(ts_op(-4,ht,'*'),@(x)(datestr(x,3)),'widths',0.1,'symbol','k.','grouporder',allmos,'allcol','k');
ylim([-.3,.3]);
%ylim([-1,1]);
titlename([upper(stn.station_name),' net heating (wide gray) vs. convective cooling (thin black)']);
%titlename([upper(stn.station_name),' net heating (wide gray) vs. convective cooling 3hlp (thin black)']);

% OLD figure had to be "Saved As" a .tif otherwise boxplot bug garbled the XTickLabels...
print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-hc-boxplot-looe1-ht.tif']));
