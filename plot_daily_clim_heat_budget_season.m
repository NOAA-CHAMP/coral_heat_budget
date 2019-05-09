function plot_daily_clim_heat_budget_season(stn,jds,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,doWarmix,advix,kthix)
%function plot_daily_clim_heat_budget_season(stn,jds,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,doWarmix,advix,kthix)
%
% Plot annual climatology resulting from one of the optimization attempts
% stored in struct STN.optim (v. OPTIMIZE_STATION_HEAT_BUDGET.M), starting on
% MIN(JDS) and ending on MAX(JDS) (JDS whole numbers between 1 and 366).
%
% Last Saved Time-stamp: <Fri 2016-06-17 19:26:36 Eastern Daylight Time gramer>

  if ( ~exist('doWarmix','var') ); doWarmix = 1; end;
  if ( ~exist('advix','var') ); advix = 1; end;
  if ( ~exist('kthix','var') ); kthix = 1; end;

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

  begyr = get_year(stn.optim.sq.date(1));
  endyr = get_year(stn.optim.sq.date(end)+1);

  % doWarm = stn.optim.doWarm(1,doWarmix,1,1);
  % if ( doWarm ); doWarmStr = 'Warm Layer';
  % else;	         doWarmStr = 'No Warm Layer';
  % end;

  % E.g., stn.ndbc_sea_t_dly vs. ndbc_erai_erai_30a_avhrr_hc_dTdt_dly vs. ndbc_erai_erai_30a_avhrr_hc_dTdt_err_1_d_sum
  [t,dt,q,qe] = ...
      intersect_tses(stn.(sfld),stn.(bdTfld),stn.(hcdTdt),stn.([hcdTdt,'_err']));


  fmg;
  climsq = squeeze(stn.optim.climsq(:,doWarmix,advix,kthix));
  climsq_err = stn.optim.climsq_err(end,doWarmix,advix,kthix);
  climsq_minus_err = ts_op(climsq(end),climsq_err,'-');
  climsq_plus_err = ts_op(climsq(end),climsq_err,'+');
  % % climsq_minus_err = ts_op(climsq(end),ts_op(climsq_err,24,'*'),'-');
  % % climsq_plus_err = ts_op(climsq(end),ts_op(climsq_err,24,'*'),'+');
  % % % climsq_minus_err.date = climsq_err.date;
  % % % % climsq_minus_err.data = climsq(end).data - cumsum(climsq_err.data);
  % % % climsq_minus_err.data = climsq(end).data - (climsq_err.data.*24);
  % % % climsq_plus_err.date = climsq_err.date;
  % % % % climsq_plus_err.data = climsq(end).data + cumsum(climsq_err.data);
  % % % % plot_ts(stn.optim.climt,climsq,climsq_minus_err,'k^',climsq_plus_err,'kv');
  % % % climsq_plus_err.data = climsq(end).data + (climsq_err.data.*24);
  % lhs=plot_ts(stn.optim.climt,'k','LineWidth',3,climsq,climsq_minus_err,'k:',climsq_plus_err,'k:');
  lhs=plot_ts(stn.optim.climt,'k','LineWidth',3,climsq,'Color',[.5,.5,.5],'LineWidth',1.5,climsq_minus_err,'k:','LineWidth',1.5,climsq_plus_err,'k:','LineWidth',1.5);
  lhs(end) = [];

  % climssqt = squeeze(stn.optim.climssqt(end,doWarmix,advix,kthix));
  % climsbqt = squeeze(stn.optim.climsbqt(end,doWarmix,advix,kthix));
  % climsdt = squeeze(stn.optim.climsdt(end,doWarmix,advix,kthix));
  % lhs(end+1:end+3) = plot_ts(climssqt,'r--',climsbqt,'co',climsdt,'m-.');

  datetick3('x',3);
  % titlename([upper(stn.station_name) ' Daily Clim: ' doWarmStr ' ' stn.optim.cbdstrs{1}]);
  titlename([upper(stn.station_name) ' Daily Clim: ' stn.optim.cbdstrs{1}]);
  % % legend({'T_s',stn.optim.kdstrs{:},[stn.optim.kdstrs{end},' \pm error']});
  kdstrs = strcat( stn.optim.kdstrs',{' Errs '},...
                   cellstr(num2str(stn.optim.dayrmse(:,doWarmix,advix,kthix),'%.1f')),{', '},...
                   cellstr(num2str(stn.optim.dayerror(:,doWarmix,advix,kthix),'%.1f')),{', '},...
                   cellstr(num2str(stn.optim.climerror(:,doWarmix,advix,kthix),'%.1f')) );
  %%%%DEBUG:
  for kdstrix=1:numel(kdstrs); disp(kdstrs{kdstrix}); end;
  % legend(lhs,{'T_s',kdstrs{:},[stn.optim.kdstrs{end},' \pm error']}, 'Location','South');
  legend(lhs,...
         {'T_s',kdstrs{:},[stn.optim.kdstrs{end},' \pm error'],...
         % 'Q_0/\rhoC_ph','(Q_0(\gamma)+Q_b)/\rhoC_ph','Non-HC',...
         }, 'Location','South');
  xlim(stn.optim.climt.date([1 end]));
  ylim([minmin([stn.optim.climsq.data]),maxmax([stn.optim.climsq.data])]);
  ylim([16,34]);
  % ylim([16,50]);
  ylabel('^oC');
  % % % appendtitlename([' (' strrep(QEPFX,'_','\_') ' Adv=' stn.optim.advfacstrs{advix} ' K_\theta=' stn.optim.kthstrs{kthix} ')' stn.commentstr]);
  % % appendtitlename([' (' strrep(hcdTdt,'_','\_') ' Adv=' stn.optim.advfacstrs{advix} ' K_\theta=' stn.optim.kthstrs{kthix} ')' stn.commentstr]);
  % appendtitlename([' ' strrep(hcdTdt,'_','\_') ': ' stn.commentstr]);
  appendtitlename([' ' stn.commentstr]);
  appendtitlename([' (' num2str(begyr) '-' num2str(endyr) ')']);

return;
