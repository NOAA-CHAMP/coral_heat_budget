1;
% Script COMPARE_DAILY_FLUX.M
%
% Compare heat budget results to "daily variability" in sea temperature
%
% Last Saved Time-stamp: <Sun 2016-07-24 16:03:10 Eastern Daylight Time gramer>


if ( ~exist('stnm','var') )
  stnm = 'mlrf1';
end;

if ( ~exist('stn','var') || ~isfield(stn,'station_name') || ~strcmpi(stn.station_name,stnm) )
  stn=[]; clear stn
  matfname = fullfile(get_thesis_path('../data'),[stnm,'-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat']);
  if ( ~exist(matfname,'file') )
    error('Not sure what heat budget to load for %s!',stnm);
  end;
  disp(['Loading ',matfname]);
  load(matfname)
end;

hq0 = stn.ndbc_erai_erai_30a_net_flux_term;
hqb = stn.b_ndbc_erai_erai_30a_net_flux_term;
hth = stn.b_ndbc_erai_erai_30a_avhrr_dt;
hhc = stn.ndbc_erai_erai_30a_avhrr_hc_dTdt;
ht = stn.ndbc_sea_t;

dq0 = stn.ndbc_erai_erai_30a_net_flux_term_dly;
dqb = stn.b_ndbc_erai_erai_30a_net_flux_term_dly;
dth = stn.b_ndbc_erai_erai_30a_avhrr_dt_dly;
dhc = stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly;
dt = stn.ndbc_sea_t_dly_diff;


if 0;
  % Net surface heat flux: R^2 ~x%, beta ~ x
  scatter_fit_ts(dq0,dt,[],[],'Q_0/\rhoC_ph','\Delta\mu_1_dT',[],[],true,0);

  % Net surface heat flux: R^2 ~x%, beta ~ x
  scatter_fit_ts(dqb,dt,[],[],'(Q_0+Q_b)/\rhoC_ph','\Delta\mu_1_dT',[],[],true,0);

  % Total Heat Budget without HC, R^2 ~ x%, beta ~ x ???
  scatter_fit_ts(dth,dt,[],[],'\partial_tT - HC','\Delta\mu_1_dT',[],[],true,0);

  % Total Heat Budget: R^2 ~ x%, beta ~ x
  scatter_fit_ts(dhc,dt,[],[],'\partial_tT','\Delta\mu_1_dT',[],[],true,0);
end;

% Eliminate interannual variability for comparison purposes
[q0,qb,th,hc,t] = intersect_tses(dq0,dqb,dth,dhc,dt);

lm = abs(nanmax([abs(nanmin([t.data;q0.data;hc.data])),nanmax([t.data;q0.data;hc.data])]));

sqrt(sum((q0.data-t.data).^2))
sqrt(sum((qb.data-t.data).^2))
sqrt(sum((th.data-t.data).^2))
sqrt(sum((hc.data-t.data).^2))

if 1;
  % Net surface heat flux: R^2 ~x%, beta ~ x
  scatter_fit_ts(q0,t,[],[],'Q_0/\rhoC_ph','\Delta\mu_1_dT',[],[],true,0); axis([-lm,+lm,-lm,+lm]);
  scatter_fit_ts_seasons(q0,t,[],[],'Q_0/\rhoC_ph','\Delta\mu_1_dT',[],[],true,0,[-lm,+lm,-lm,+lm]);

  % Net vertical heat flux: R^2 ~x%, beta ~ x
  scatter_fit_ts(qb,t,[],[],'(Q_0+Q_b)/\rhoC_ph','\Delta\mu_1_dT',[],[],true,0); axis([-lm,+lm,-lm,+lm]);
  scatter_fit_ts_seasons(qb,t,[],[],'(Q_0+Q_b)/\rhoC_ph','\Delta\mu_1_dT',[],[],true,0); axis([-lm,+lm,-lm,+lm]);

  % Total Heat Budget without HC: R^2 ~ x%, beta ~ x
  scatter_fit_ts(th,t,[],[],'\partial_tT - HC','\Delta\mu_1_dT',[],[],true,0); axis([-lm,+lm,-lm,+lm]);
  scatter_fit_ts_seasons(th,t,[],[],'\partial_tT - HC','\Delta\mu_1_dT',[],[],true,0); axis([-lm,+lm,-lm,+lm]);

  % Total Heat Budget: R^2 ~ x%, beta ~ x
  scatter_fit_ts(hc,t,[],[],'\partial_tT','\Delta\mu_1_dT',[],[],true,0); axis([-lm,+lm,-lm,+lm]);
  scatter_fit_ts_seasons(hc,t,[],[],'\partial_tT','\Delta\mu_1_dT',[],[],true,0,[-lm,+lm,-lm,+lm]);
end;

if 0;
  fmg; boxplot_ts(t); titlename('\Delta\mu_1_dT'); ylim([-lm,+lm]);
  fmg; boxplot_ts(q0); titlename('Q_0/\rhoC_ph'); ylim([-lm,+lm]);
  fmg; boxplot_ts(hc); titlename('\partial_tT'); ylim([-lm,+lm]);
end;

if 0;
  % Total Heat Budget with +/- 1 day lags: R^2 ~ x%, beta ~ x
  scatter_fit_ts(dhc,dt,[],[],'\partial_tT + 1 d','\Delta\mu_1_dT',[],[],true,+1);
  scatter_fit_ts(dhc,dt,[],[],'\partial_tT - 1 d','\Delta\mu_1_dT',[],[],true,-1);

  % Must recalculate the *_dly values to compare one- to 23-hour lags...
  for lagh = 0:6:21;
    [rawts.data,rawts.date] = grp_ts(ht.data,ht.date,@floor,@nanmean,23);
    rawt.date = rawts.date(2:end);
    rawt.data = diff(rawts.data);

    [lagq0,lagth,laghc,lagt] = intersect_tses(dq0,dth,dhc,dt);

    [lagt.data,lagt.date] = grp_ts(rawt.data(lagh+1:end),rawt.date(lagh+1:end)-(lagh/24),@floor,@nansum,24);
    scatter_fit_ts(q0,lagt,[],[],'Q_0/\rhoC_ph',sprintf('%s %g h lag','\Delta\mu_1_dT',lagh),[],[],true,0); axis([-lm,+lm,-lm,+lm]);
    scatter_fit_ts(hc,lagt,[],[],'\partial_tT',sprintf('%s %g h lag','\Delta\mu_1_dT',lagh),[],[],true,0); axis([-lm,+lm,-lm,+lm]);
  end;
end;

if 1;
  [tmx.data,tmx.date,ig,ig,ig,ig,fdt.data,fdt.date] = grp_ts(ht.data,ht.date,@floor,@nanmax,24);
  [tmn.data,tmn.date] = grp_ts(ht.data,ht.date,@floor,@nanmin,24);

  % Daily range (max - min)
  tr.date=tmn.date; tr.data=tmx.data-tmn.data;

  fmg; boxplot_ts(tr); titlename('Daily sea temperature range');
  scatter_fit_ts(q0,tr)
  scatter_fit_ts(qb,tr)
  scatter_fit_ts(th,tr)
  scatter_fit_ts(hc,tr)

  % Daily change (hour 23 - hour 0 of each full day)
  fdd.date = fdt.date(1:24:end);
  fdd.data = fdt.data(24:24:end) - fdt.data(1:24:end);

  fmg; boxplot_ts(fdd); titlename('Daily sea temperature change');
  scatter_fit_ts(q0,fdd); axis([-3,6,-3,6])
  scatter_fit_ts(qb,fdd); axis([-3,6,-3,6])
  scatter_fit_ts(th,fdd); axis([-3,6,-3,6])
  scatter_fit_ts(hc,fdd); axis([-3,6,-3,6])
end;
