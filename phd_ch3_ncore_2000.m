1;

if ( ~exist('stns','var') )
    stns = get_ncore_2000;
end;

shr = ts_op(stns.csh.cm_xshore_3hlp,stns.cdp.cm_xshore_3hlp,'-');
dTdz = ts_op(stns.csh.cm_seatemp_3hlp,stns.cdp.cm_seatemp_3hlp,'-');

udTdz = ts_op(shr,dTdz,'*');

fmg; plot_ts(stns.cdp.cm_xshore_3hlp,stns.csh.cm_xshore_3hlp,stns.cdp.cm_seatemp_3hlp,stns.csh.cm_seatemp_3hlp,stn.ndbc_wind1_speed);
datetick3;
legend('Deep XS','Shal XS','Deep T','Shal T','Wind');

% [anomts,clim,tid] = anomalize_ts(stns.cdp.cm_seatemp_3hlp,@floor);
% udTdz = ts_op(shr,anomts,'*');
% [u,U,V,Hs]=intersect_tses(udTdz,stn.ndbc_wind1_u,stn.ndbc_wind1_v,stn.ww3_sigwavehgt);
% fmg; plot_ts(stns.cdp.cm_xshore_3hlp,stns.csh.cm_xshore_3hlp,udTdz,stn.ndbc_wind1_speed);
% legend('Deep XS','Shal XS','\DeltaT','Wind');
% datetick3;

%[udT,U,V,Hs]=intersect_tses(udTdz,stn.ndbc_wind1_u,stn.ndbc_wind1_v,stn.ww3_sigwavehgt);
%fmg; boxplot_tses({U,V,Hs},@(x)(datestr(x,12))); boxplot_ts(udT,@(x)(datestr(x,12)),'allcol','r'); ylim([-20,20]);


[anom1d,cum1d,tid1d] = anomalize_ts(stns.cdp.cm_seatemp_3hlp,@floor);
clim1d.date = tid1d; clim1d.data = cum1d; 
[anom3d,cum3d,tid3d] = anomalize_ts(stns.cdp.cm_seatemp_3hlp,@get_yeartriad);
clim3d.date = tid3d; clim3d.data = cum3d; 
[anom7d,cum7d,tid7d] = anomalize_ts(stns.cdp.cm_seatemp_3hlp,@get_yearweek);
clim7d.date = tid7d; clim7d.data = cum7d; 

[anom1dsh,cum1dsh,tid1dsh] = anomalize_ts(stns.csh.cm_seatemp_3hlp,@floor);
clim1dsh.date = tid1dsh; clim1dsh.data = cum1dsh; 
[anom3dsh,cum3dsh,tid3dsh] = anomalize_ts(stns.csh.cm_seatemp_3hlp,@get_yeartriad);
clim3dsh.date = tid3dsh; clim3dsh.data = cum3dsh; 
[anom7dsh,cum7dsh,tid7dsh] = anomalize_ts(stns.csh.cm_seatemp_3hlp,@get_yearweek);
clim7dsh.date = tid7dsh; clim7dsh.data = cum7dsh; 

anom1dshdp = ts_op(clim1dsh,stns.cdp.cm_seatemp_3hlp,'-');
anom3dshdp = ts_op(clim3dsh,stns.cdp.cm_seatemp_3hlp,'-');
anom7dshdp = ts_op(clim7dsh,stns.cdp.cm_seatemp_3hlp,'-');

udTdz1d = ts_op(shr,anom1dshdp,'*');
udTdz3d = ts_op(shr,anom3dshdp,'*');
udTdz7d = ts_op(shr,anom7dshdp,'*');

allmos = datestr(datenum(1,1:12,1),3);
fmg; boxplot_ts(udTdz1d,@(x)(datestr(x,3)),'grouporder',allmos); ylim([-20,20]); titlename('\DeltaV/\Deltaz \times anom_1_d T_d');
fmg; boxplot_ts(udTdz3d,@(x)(datestr(x,3)),'grouporder',allmos); ylim([-20,20]); titlename('\DeltaV/\Deltaz \times anom_3_d T_d');
fmg; boxplot_ts(udTdz7d,@(x)(datestr(x,3)),'grouporder',allmos); ylim([-20,20]); titlename('\DeltaV/\Deltaz \times anom_7_d T_d');


clear -regexp (anom|clim|cum|tid)[137]*


stn = get_station_from_station_name('mlrf1'); stn = get_avhrr_weekly_field(stn,true); stn = load_all_ndbc_data(stn); stn = get_ngdc_bathy_station(stn); stn = get_erai_station(stn); stn = adjust_erai_station_waves(stn); stn = get_ww3_station(stn); stn=station_spddir_to_uv(stn,'ndbc_wind1_speed','ndbc_wind1_dir'); stn=station_optimal_isobath_orientation(stn); stn=station_reorient_vectors(stn,stn.isobath_orientation,'ndbc_wind1_u','ndbc_wind1_v');

if ( exist('stn','var') && exist('lonf1','var') )
  plot_ngdc_bathy_station(stn,-[0:2:100]); plot(stns.a.lon,stns.a.lat,'wd'); plot(stns.b.lon,stns.b.lat,'wd'); plot(stns.csh.lon,stns.csh.lat,'wd'); plot(stns.csh1.lon,stns.csh1.lat,'ws'); plot(stns.csh2.lon,stns.csh2.lat,'ws'); plot(stns.csh3.lon,stns.csh3.lat,'ws'); plot(stns.tl1.lon,stns.tl1.lat,'wv'); plot(stns.tl2.lon,stns.tl2.lat,'wv'); plot(stns.tl3.lon,stns.tl3.lat,'wv'); plot(stns.tl4.lon,stns.tl4.lat,'wv'); axis([-80.42,-80.28,24.98,25.12]); daspect([cosd(25),1,1]);

  fmg; plot_ts(stns.a.cm_seatemp_3hlp,stns.b.cm_seatemp_3hlp,stns.tl4.tl_seatemp,'-',stns.tl3.tl_seatemp,'-',stns.tl2.tl_seatemp,'-',stns.tl1.tl_seatemp,'-',stns.csh.cm_seatemp_3hlp,stns.cdp.cm_seatemp_3hlp,'.-','Color',[.5,.5,.5],stn.ndbc_sea_t,'k-','LineW',1.5,lonf1.ndbc_sea_t,'y-','LineW',1.5); legend('A','B','TL4','TL3','TL2','TL1','Csh','Cdp','MLRF1','LONF1');
  xlim(datenum(2001,7,[12,16])); datetick3;
  titlename('Sea temperatures from NCORE 2000-2002');

  %fmg; plot_ts(stns.a.cm_u_3hlp,stns.a.cm_v_3hlp,stns.b.cm_u_3hlp,stns.b.cm_v_3hlp,stns.csh.cm_u_3hlp,stns.csh.cm_v_3hlp,stns.cdp.cm_u_3hlp,stns.cdp.cm_v_3hlp); legend('A u','A v','B u','B v','Csh u','Csh v','Cdp u','Cdp v');

  % fmg; plot_ts(stns.a.cm_xshore_3hlp,stns.a.cm_lshore_3hlp,stns.b.cm_xshore_3hlp,stns.b.cm_lshore_3hlp,stns.csh.cm_xshore_3hlp,stns.csh.cm_lshore_3hlp,stns.cdp.cm_xshore_3hlp,stns.cdp.cm_lshore_3hlp); legend('A x','A l','B x','B l','Csh x','Csh l','Cdp x','Cdp l');
  % xlim(datenum(2001,7,[12,16])); datetick3;

  fmg; plot_ts(stns.a.cm_xshore_3hlp,stns.b.cm_xshore_3hlp,stns.csh.cm_xshore_3hlp,stns.cdp.cm_xshore_3hlp,stn.ndbc_wind1_xshore,'k-','LineW',1.5); legend('A x','B x','Csh x','Cdp x','MLRF1 x');
  ylim([-30,30]);
  xlim(datenum(2001,7,[12,16])); datetick3; ylim([-12,12]);
  titlename('Cross-shore currents from NCORE 2000-2002');
end;
