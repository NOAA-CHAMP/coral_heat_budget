function phd_ch3_daily_diagnostic(stn,ncr,yr,mo,dys)

  if ( nargin<3 )
    yr=2007; mo=7; dys=1:19;
    %yr=2007; mo=6; dys=11:29;
  end;

  ndys = numel(dys);
  nrows = floor(sqrt(ndys));
  ncols = ceil(sqrt(ndys));

  fmg;
  for dix=1:ndys
    dy = dys(dix);
    dt = datenum(yr,mo,dy);
    fs = @(x)(find(dt==floor(x.date)));
    t0 = stn.ndbc_sea_t.data(fs(stn.ndbc_sea_t)); t0=t0(1);

    if ( dix==ndys )
      subplot_tight(nrows,ncols,[dix,dix+1]);
    else
      subplot_tight(floor(sqrt(numel(dys))),ceil(sqrt(numel(dys))),dix);
    end;

    hold on;
    % ts=subset_ts(ncr.cm_deep_xdTdz,fs); plot(get_hour(ts.date),ts.data.*10,'b-');
    % %ts=subset_ts(ncr.cm_dxTdz,fs); plot(get_hour(ts.date),ts.data.*10,'b-');
    % %ts=subset_ts(ncr.cm_deep_seatemp,fs); plot(get_hour(ts.date),ts.data-nanmean(ts.data),'r-');

    ts=subset_ts(ncr.cm_dTdz,fs); plot(get_hour(ts.date),ts.data./3,'c-');
    ts=subset_ts(ncr.cm_deep_xshore,fs); plot(get_hour(ts.date),ts.data.*10,'b-');
%    text(2,-0.9,'Upslope','Color','b');

    % % 0.1K/hour warming at MLRF1 == approx. 1200 W/m^2 net surface heat flux
    % ts=subset_ts(stn.simple_ndbc_erai_erai_30a_net_flux,fs); plot(get_hour(ts.date),ts.data./1000,'m-');
    ts=subset_ts(stn.simple_ndbc_erai_erai_30a_net_flux_term,fs); plot(get_hour(ts.date),cumsum(ts.data),'m-','LineWidth',1.5);
    %ts=subset_ts(stn.b_ndbc_erai_erai_30a_avhrr_dt,fs); plot(get_hour(ts.date),cumsum(ts.data),'r-');
    ts=subset_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt,fs); plot(get_hour(ts.date),cumsum(ts.data),'r-');

    % ts=subset_ts(stn.ndbc_sea_t,fs); plot(get_hour(ts.date),ts.data-nanmean(ts.data),'k-','LineWidth',1.5);
    ts=subset_ts(stn.ndbc_sea_t,fs); plot(get_hour(ts.date),ts.data-t0,'k-','LineWidth',1.5);

    % %ts=subset_ts(stn.ndbc_erai_erai_30a_wind_stress_xshore,fs); plot(get_hour(ts.date),ts.data,'c-');
    % ts=subset_ts(stn.ndbc_erai_erai_30a_wind_stress,fs); plot(get_hour(ts.date),ts.data,'c-');

    ts=subset_ts(ncr.cm_deep_seapres,fs); plot(get_hour(ts.date),ts.data-nanmean(ts.data),'-','Color',[.8,.8,.8],'LineWidth',1.5);

    % ts=subset_ts(stn.bic_surf_par,fs); plot(get_hour(ts.date),ts.data./1500,'-','Color',[.9,.9,0]);
    % % %ts=subset_ts(stn.ncep_par,fs); plot(get_hour(ts.date)+2,ts.data./1500,'-','Color',[0,.5,0]);
    % % %ts=subset_ts(stn.ncep_dsrf,fs); plot(get_hour(ts.date)+2,ts.data./1500,'-','Color',[0,.5,0]);

    %axis([0,23,-10,+10]);
    axis([0,23,-1.2,+1.2]);
    grid on;
    xlabel(num2str(get_dom(dt)));
  end;

  suptitlename([upper(ncr.station_name),' vs. ',upper(stn.station_name),' ',datestr(datenum(yr,mo,dys(1))),' - ',datestr(datenum(yr,mo,dys(end)))]);

  %legend('dT/dz/3','u^B^T^M_x_s*10','\SigmaQ_0/\rhoC_ph','\Sigma\partial_tT','T_s-T_0','P^B^T^M_s','PAR/1500', 'Location','EastOutside');
  legend('dT/dz/3','u^B^T^M_x_s*10','\SigmaQ_0/\rhoC_ph','\Sigma\partial_tT','T_s-T_0','P^B^T^M_s', 'Location','EastOutside');

return;
