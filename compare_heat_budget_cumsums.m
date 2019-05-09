function compare_heat_budget_cumsums(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names)
%function compare_heat_budget_cumsums(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names)
%
% Plot comparison of simple cumulative sums of daily climatology vs. heat budget
% estimates. This compares both annual amplitude and interannual variability.
%
% USES: dsffld,climq0fld,sq0fld,bq0fld,qtAdvffld,bdTffld,hcdTdtf
%
% Last Saved Time-stamp: <Tue 2012-07-31 16:57:48  lew.gramer>


  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

  err_dt = 60*24;

  [dix,cix,six,bix,aix,kix,hix] = ...
      intersect_all_dates([],stn.(dsffld).date,stn.(climq0fld).date,stn.(sq0fld).date,stn.(bq0fld).date,stn.(qtAdvffld).date,stn.(bdTffld).date,stn.(hcdTdtf).date);

  fmg;
  lh=[];
  lh(end+1)=plot(stn.(dsffld).date(dix(1):end),nancumsum(stn.(dsffld).data(dix(1):end)),'k');
  lh(end+1)=plot(stn.(climq0fld).date(cix(1):end),nancumsum(stn.(climq0fld).data(cix(1):end).*24),'c');

  cs.date = stn.(sq0fld).date(six(1):end);
  cs.data = nancumsum(stn.(sq0fld).data(six(1):end));
  lh(end+1)=plot(cs.date,cs.data,'r-.');
  if ( isfield(stn,[sq0fld,'_err']) )
    [csix,serrix] = intersect_dates(cs.date,stn.([sq0fld,'_err']).date);
    csp = cs.data(csix) + nancumsum(real(stn.([sq0fld,'_err']).data(serrix)));
    plot(stn.([sq0fld,'_err']).date(serrix(1:1:end)),csp(1:1:end),'r-.');
    plot(stn.([sq0fld,'_err']).date(serrix(1:err_dt:end)),csp(1:err_dt:end),'rv','MarkerSize',5);
    csn = cs.data(csix) - nancumsum(real(stn.([sq0fld,'_err']).data(serrix)));
    plot(stn.([sq0fld,'_err']).date(serrix(1:1:end)),csn(1:1:end),'r-.');
    plot(stn.([sq0fld,'_err']).date(serrix(1:err_dt:end)),csn(1:err_dt:end),'r^','MarkerSize',5);
  end;

  bcs.date = stn.(bq0fld).date(bix(1):end);
  bcs.data = nancumsum(stn.(bq0fld).data(bix(1):end));
  lh(end+1)=plot(stn.(bq0fld).date(bix(1):end),nancumsum(stn.(bq0fld).data(bix(1):end)),'m:');
  if ( isfield(stn,[bq0fld,'_err']) )
    [bcsix,qerrix] = intersect_dates(bcs.date,stn.([bq0fld,'_err']).date);
    bcsp = bcs.data(bcsix) + nancumsum(real(stn.([bq0fld,'_err']).data(qerrix)));
    plot(stn.([bq0fld,'_err']).date(qerrix(1:1:end)),bcsp(1:1:end),'m-.');
    plot(stn.([bq0fld,'_err']).date(qerrix(1:err_dt:end)),bcsp(1:err_dt:end),'mv','MarkerSize',5);
    bcsn = bcs.data(bcsix) - nancumsum(real(stn.([bq0fld,'_err']).data(qerrix)));
    plot(stn.([bq0fld,'_err']).date(qerrix(1:1:end)),bcsn(1:1:end),'m-.');
    plot(stn.([bq0fld,'_err']).date(qerrix(1:err_dt:end)),bcsn(1:err_dt:end),'m-.^','MarkerSize',5);
  end;

  lh(end+1)=plot(stn.(qtAdvffld).date(aix(1):end),nancumsum(stn.(qtAdvffld).data(aix(1):end)),'y--');
  lh(end+1)=plot(stn.(bdTffld).date(kix(1):end),nancumsum(stn.(bdTffld).data(kix(1):end)),'y:');
  lh(end+1)=plot(stn.(hcdTdtf).date(hix(1):end),nancumsum(stn.(hcdTdtf).data(hix(1):end)),'b--');
  legend(lh,'Actual \rhoC_ph\partial_tT_s','OAFlux/ISCCP Q_0',...
         'G&M Q_0 \pm\sigmaQ_0','G&M Q_0(\gamma)+Q_b \pm\sigmaQ','G&M Q_0(\gamma)+Q_b + \rhoC_phu^.\nablaT',...
         'G&M Q_0(\gamma)+Q_b + \rhoC_ph(u^.\nablaT + K\nabla^2T)',...
         'G&M \rhoC_ph\partial_tT_s',...
         'Location','Best');
  datetick3('x',2,'keeplimits');
  % yl = ylim;  yl = max([abs(yl(:));3e6]);
  yl = 3e6;
  ylim([-yl,+yl]);    ylabel('W/m^2');
  titlename([upper(stn.station_name) ' simple cumulative sums: Hourly fluxes']);
  if ( isfield(stn,'commentstr') )
    appendtitlename([stn.commentstr]);
  end;

return;
