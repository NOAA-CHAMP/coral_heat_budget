function res = compare_flux_climatologies(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,commentstr)
%function res = compare_flux_climatologies(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,commentstr)
%
% Plot annual or interannual heat-budget term climatologies from Gramer and
% Mariano (2012), OAFlux (Yu and Weller 2007) and ISCCP (Zhang et al. 2004).
%
% Last Saved Time-stamp: <Tue 2013-06-18 15:57:25 Eastern Daylight Time gramer>

  %DEBUG:  tic,

  if ( ~exist('per','var') || isempty(per) )
    per = 'daily';
  end;
  if ( exist('commentstr','var') && ~isempty(commentstr) )
    stn.commentstr = commentstr;
  elseif ( ~isfield(stn,'commentstr') )
    stn.commentstr = '';
  end;

  forceImplied = false;
  %forceImplied = true;

  doOtherTerms = false;
  %doOtherTerms = true;


  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

  % Would prefer MEDIAN, but published climatologies generally use MEAN (and
  % for, e.g., daily insolation, these give very different results!)
  %sumfun = @nanmedian;
  sumfun = @nanmean;

  yrs = get_year(stn.(sfld).date);
  uyrs = unique(yrs);
  nyrs = numel(uyrs);

  % Allowing leap-days might cause problems when comparing cum stats on
  % *different time series*, e.g., a time series which includes data for
  % one or more 29th's of Feb, and one which does not. If caller does not
  % like this behavior, they can simply specify custom CUMFUN and MINN.
  switch ( per ),
   case 'hourly',   n=365*24; minN=1;      mrk=24; N=1;
   case 'daily',    n=365;    minN=23;     mrk=30; N=24;
   case 'pentad',   n=73;     minN=23*3;   mrk= 5; N=24*5;
   case 'weekly',   n=52;     minN=23*4;   mrk= 4; N=24*7;
   case 'monthly',  n=12;     minN=23*16;  mrk= 1; N=24*[31,28,31,30,31,30,31,31,30,31,30,31]';
   case 'seasonal', n=4;      minN=23*46;  mrk= 1; N=24*[31+28+31,30+31+30,31+31+30,31+30+31]';
   % With interannual comparison, everything is different...
   case 'yearly',   n=nyrs;   minN=23*273; mrk= 1; N=1;
   otherwise,       error('Do not know how to do a "%s"-period climatology!',char(per));
  end;

  rawt = subset_ts(stn.(sfld),@(x)(find(ismember(get_year(x.date),unique(get_year(stn.(climsrfld).date))))));
  [t,dtf,q0,qt,sq0,sqt,bq0,bq0t,bdT,bdTf,udT,kd2T,hcdTdt,hc,sr,asr,lr,qlh,qsh,qrh,qbo] = ...
      intersect_tses(rawt,stn.(dsffld),stn.(q0fld),stn.(qtfld),stn.(sqtfld),stn.(bq0tfld),stn.(sq0fld),stn.(bq0fld),stn.(bdTfld),stn.(bdTffld),stn.(fqudTffld),stn.(kd2Tffld),stn.(hcdTdt),stn.(hcdTdthcf),stn.(srfld),stn.(asrfld),stn.(lrfld),stn.(qlhfld),stn.(qshfld),stn.(qrhfld),stn.(qbofld));


  %SEE function [cum,tid] = grp_ts(dat,dts,per_or_cumfun,sumfun,minN)

  [res.t.data,res.t.date,nPerYr,delt,cumfun] = grp_ts(t.data,t.date,per,sumfun,minN);

  [res.dtf.data,res.dtf.date] = grp_ts(dtf.data,dtf.date,per,sumfun,minN);

  [res.q0.data,res.q0.date] = grp_ts(q0.data,q0.date,per,sumfun,minN);

  [res.raw_qt.data,res.raw_qt.date] = grp_ts(qt.data,qt.date,per,sumfun,minN);
  res.qt.date=res.raw_qt.date; res.qt.data=N.*cumsum(res.raw_qt.data);

  [res.sq0.data,res.sq0.date] = grp_ts(sq0.data,sq0.date,per,sumfun,minN);

  [res.raw_sqt.data,res.raw_sqt.date] = grp_ts(sqt.data,sqt.date,per,sumfun,minN);
  res.sqt.date=res.raw_sqt.date; res.sqt.data=N.*cumsum(res.raw_sqt.data);

  [res.bq0.data,res.bq0.date] = grp_ts(bq0.data,bq0.date,per,sumfun,minN);

  [res.raw_bq0t.data,res.raw_bq0t.date] = grp_ts(bq0t.data,bq0t.date,per,sumfun,minN);
  res.bq0t.date=res.raw_bq0t.date; res.bq0t.data=N.*cumsum(res.raw_bq0t.data);

  % % Special handling for huge advection swings???
  % [res.raw_udT.data,res.raw_udT.date] = grp_ts(udT.data,udT.date,per,sumfun,minN);
  % res.udT.date=res.raw_udT.date; res.udT.data=N.*cumsum(res.raw_udT.data);

  % [res.raw_kd2T.data,res.raw_kd2T.date] = grp_ts(kd2T.data,kd2T.date,per,sumfun,minN);
  % res.kd2T.date=res.raw_kd2T.date; res.kd2T.data=N.*cumsum(res.raw_kd2T.data);

  [res.raw_bdT.data,res.raw_bdT.date] = grp_ts(bdT.data,bdT.date,per,sumfun,minN);
  res.bdT.date=res.raw_bdT.date; res.bdT.data=N.*cumsum(res.raw_bdT.data);

  [res.raw_hcdTdt.data,res.raw_hcdTdt.date] = grp_ts(hcdTdt.data,hcdTdt.date,per,sumfun,minN);
  res.hcdTdt.date=res.raw_hcdTdt.date; res.hcdTdt.data=N.*cumsum(res.raw_hcdTdt.data);


  [res.sr.data,res.sr.date] = grp_ts(sr.data,sr.date,per,sumfun,minN);

  [res.asr.data,res.asr.date] = grp_ts(asr.data,asr.date,per,sumfun,minN);

  [res.lr.data,res.lr.date] = grp_ts(lr.data,lr.date,per,sumfun,minN);

  [res.qlh.data,res.qlh.date] = grp_ts(qlh.data,qlh.date,per,sumfun,minN);

  [res.qsh.data,res.qsh.date] = grp_ts(qsh.data,qsh.date,per,sumfun,minN);

  [res.qrh.data,res.qrh.date] = grp_ts(qrh.data,qrh.date,per,sumfun,minN);


  [res.qbo.data,res.qbo.date] = grp_ts(qbo.data,qbo.date,per,sumfun,minN);

  [res.udT.data,res.udT.date] = grp_ts(udT.data,udT.date,per,sumfun,minN);

  [res.kd2T.data,res.kd2T.date] = grp_ts(kd2T.data,kd2T.date,per,sumfun,minN);

  [res.bdTf.data,res.bdTf.date] = grp_ts(bdTf.data,bdTf.date,per,sumfun,minN);

  [res.hc.data,res.hc.date] = grp_ts(hc.data,hc.date,per,sumfun,minN);


  % Pare down all fields to just those "dates" that mutually intersect
  [res.t,res.dtf,res.q0,res.qt,res.sq0,res.sqt,res.bq0,res.bq0t,...
   res.udT,res.kd2T,res.hc,res.bdT,res.bdTf,...
   res.hcdTdt,res.sr,res.asr,res.lr,res.qlh,res.qsh,res.qrh,res.qbo] = ...
      intersect_tses(res.t,res.dtf,res.q0,res.qt,res.sq0,res.sqt,res.bq0,res.bq0t,...
                     res.udT,res.kd2T,res.hc,res.bdT,res.bdTf,...
                     res.hcdTdt,res.sr,res.asr,res.lr,res.qlh,res.qsh,res.qrh,res.qbo);

  res.radif = ts_op(res.sr,res.lr,'+');
  res.aradif = ts_op(res.asr,res.lr,'+');
  res.turif = ts_op(ts_op(res.qlh,res.qsh,'+'),res.qrh,'+');
  res.coolif = ts_op(res.lr,res.turif,'+');

  % Limit climatological comparisons to only valid year/periods based on MINN
  t.data(~ismember(cumfun(t.date),res.t.date))=[];
  t.date(~ismember(cumfun(t.date),res.t.date))=[];

  % Climatology has one value per day!
  raw_climsr = subset_ts(stn.(climsrfld),@(x)(find(ismember(get_year(x.date),get_year(t.date)))));
  [climsr,climlr,climqlh,climqsh] = ...
      intersect_tses(raw_climsr,stn.(climlrfld),stn.(climqlhfld),stn.(climqshfld));

  [res.climsr.data,res.climsr.date] = grp_ts(climsr.data,climsr.date,per,sumfun,1);
  [res.climlr.data,res.climlr.date] = grp_ts(climlr.data,climlr.date,per,sumfun,1);
  [res.climqlh.data,res.climqlh.date] = grp_ts(climqlh.data,climqlh.date,per,sumfun,1);
  [res.climqsh.data,res.climqsh.date] = grp_ts(climqsh.data,climqsh.date,per,sumfun,1);

  % Make sure calculated and climatology vectors are the same length!
  [ig,res.climsr,res.climlr,res.climqlh,res.climqsh] = ...
      intersect_tses(res.t,res.climsr,res.climlr,res.climqlh,res.climqsh);

  res.climq0 = ts_op(ts_op(res.climsr,res.climlr,'+'),ts_op(res.climqlh,res.climqsh,'+'),'+');

  res.climradif = ts_op(res.climsr,res.climlr,'+');
  res.climturif = ts_op(res.climqlh,res.climqsh,'+');


  mrkix = [1:mrk:length(res.t.date)-1,length(res.t.date)];


  fh = fmg;

if (0)
  ax(1)=subplot_tight(3,1,[1]);
  hold on; grid on;
  % plot(res.t.date,[res.t.data,res.t.data(1)+res.sqt.data-res.sqt.data(1),res.t.data(1)+res.bq0t.data-res.bq0t.data(1),res.t.data(1)+res.hcdTdt.data-res.hcdTdt.data(1)]);
  % legend(ax(1),'T ','T(0)+\SigmaQ_0/\rhoC_ph ','T(0)+\Sigma(Q_0(\gamma)+Q_b)/\rhoC_ph ','T(0)+\Sigma\partial_tT_H_C ',...
  %        'Location','SouthEast', 'Orientation','horizontal');
  plot(res.t.date,res.t.data,'k-',...
       res.t.date,res.t.data(1)+res.bq0t.data-res.bq0t.data(1),'b-',...
       res.t.date,res.t.data(1)+res.hcdTdt.data-res.hcdTdt.data(1),'r-');
  plot(res.t.date(mrkix),res.t.data(mrkix),'k.',...
       res.t.date(mrkix),res.t.data(1)+res.bq0t.data(mrkix)-res.bq0t.data(1),'bs',...
       res.t.date(mrkix),res.t.data(1)+res.hcdTdt.data(mrkix)-res.hcdTdt.data(1),'ro');
  plh=plot(1,res.t.data(1),'k.-',...
           1,res.t.data(1)+res.bq0t.data(1)-res.bq0t.data(1),'bs-',...
           1,res.t.data(1)+res.hcdTdt.data(1)-res.hcdTdt.data(1),'ro-');
  legend(plh,'T','T(0)+\Sigma(Q_0(\gamma)+Q_b)/\rhoC_ph','T(0)+\Sigma\partial_tT_H_C',...
         'Location','SouthEast', 'Orientation','horizontal');
  axis([min(res.t.date),max(res.t.date), 5,45]);

  titlename([upper(stn.station_name) ': ' ...
             strrep(bdTfld,'_','\_') ' ' upper(per) ' ' upper(char(sumfun)) ' climatology ' ...
             stn.commentstr]);

  ax(2)=subplot_tight(3,1,[2,3]);
  hold on; grid on;
end;

  % Lines
  % PLOT (v.) does not accept property-value pairs as non-final arguments!
  if ( forceImplied || ~strcmp(per,'daily') )
    % Daily implied actual flux is way too noisy - blots out other plots
    plot(res.dtf.date,res.dtf.data,'k-','LineWidth',2);
  end;
  plot(res.radif.date,res.radif.data,'r-');
  plot(res.turif.date,res.turif.data,'b-');
  plot(res.aradif.date,res.aradif.data,'m-');
  plot(res.qbo.date,res.qbo.data,'k-.','LineWidth',1.5);
  if ( doOtherTerms )
    % plot(res.udT.date,res.udT.data,'k-','LineWidth',1.5,'Color',[.5,.5,.5]);
    % plot(res.kd2T.date,res.kd2T.data,'k:','LineWidth',1.5,'Color',[.5,.5,.5]);
    plot(res.bdTf.date,res.bdTf.data,'k-','LineWidth',1.5,'Color',[.5,.5,.5]);
    plot(res.hc.date,res.hc.data,'k:','LineWidth',1.5,'Color',[.5,.5,.5]);
  end;
  plot(res.climradif.date,res.climradif.data,'r:','LineWidth',1.5);
  plot(res.climturif.date,res.climturif.data,'b:','LineWidth',1.5);

  % Occasional markers (stupid MATLAB)
  if ( forceImplied || ~strcmp(per,'daily') )
    % Daily implied actual flux is way too noisy - blots out other plots
    plot(res.dtf.date(mrkix),res.dtf.data(mrkix),'k.','LineWidth',2);
  end;
  plot(res.radif.date(mrkix),res.radif.data(mrkix),'rs');
  plot(res.turif.date(mrkix),res.turif.data(mrkix),'bs');
  plot(res.aradif.date(mrkix),res.aradif.data(mrkix),'m^');
  plot(res.qbo.date(mrkix),res.qbo.data(mrkix),'k+');
  if ( doOtherTerms )
    % plot(res.udT.date(mrkix),res.udT.data(mrkix),'k+','LineWidth',1.5,'Color',[.5,.5,.5]);
    % plot(res.kd2T.date(mrkix),res.kd2T.data(mrkix),'k+','LineWidth',1.5,'Color',[.5,.5,.5]);
    plot(res.bdTf.date(mrkix),res.bdTf.data(mrkix),'kv','LineWidth',1.5,'Color',[.5,.5,.5]);
    plot(res.hc.date(mrkix),res.hc.data(mrkix),'k^','LineWidth',1.5,'Color',[.5,.5,.5]);
  end;
  plot(res.climradif.date(mrkix),res.climradif.data(mrkix),'ro');
  plot(res.climturif.date(mrkix),res.climturif.data(mrkix),'bo');

  % First marker and line - so Legend displays correctly (stupid MATLAB)
  plh=[];
  legs={};
  if ( forceImplied || ~strcmp(per,'daily') )
    % Daily implied actual flux is way too noisy - blots out other plots
    plh(end+1)=plot(res.dtf.date(1),res.dtf.data(1),'k.-','LineWidth',2);
    legs(end+1) = {'Actual'};
  end;
  plh(end+1)=plot(res.radif.date(1),res.radif.data(1),'rs-');
  legs(end+1) = {'G&M: Q_S_W+Q_L_W'};
  plh(end+1)=plot(res.turif.date(1),res.turif.data(1),'bs-');
  legs(end+1) = {'G&M: Q_L_H+Q_S_H'};
  plh(end+1)=plot(res.aradif.date(1),res.aradif.data(1),'m^-');
  legs(end+1) = {'G&M: \gammaQ_S_W+Q_L_W'};
  plh(end+1)=plot(res.qbo.date(1),res.qbo.data(1),'k-.+');
  legs(end+1) = {'G&M: Q_b'};
  if ( doOtherTerms )
    % plh(end+1)=plot(res.udT.date(1),res.udT.data(1),'k-+','LineWidth',1.5,'Color',[.5,.5,.5]);
    % legs(end+1) = {'G&M: F_qu\bullet\nablaT_s'};
    % plh(end+1)=plot(res.kd2T.date(1),res.kd2T.data(1),'k:+','LineWidth',1.5,'Color',[.5,.5,.5]);
    % legs(end+1) = {'G&M: K_H\nabla^2T_s'};
    plh(end+1)=plot(res.bdTf.date(1),res.bdTf.data(1),'k-v','LineWidth',1.5,'Color',[.5,.5,.5]);
    legs(end+1) = {'G&M: KM-scale'};
    plh(end+1)=plot(res.hc.date(1),res.hc.data(1),'k:^','LineWidth',1.5,'Color',[.5,.5,.5]);
    legs(end+1) = {'G&M: Hor. conv.'};
  end;
  plh(end+1)=plot(res.climradif.date(1),res.climradif.data(1),'ro:');
  legs(end+1) = {'ISCCP: Q_S_W+Q_L_W'};
  plh(end+1)=plot(res.climturif.date(1),res.climturif.data(1),'bo:');
  legs(end+1) = {'OAFlux: Q_L_H+Q_S_H'};

  legh=legend(plh,legs);
  % % set(legh, 'orientation','vertical', 'Location','Best');
  % set(legh,'FontSize',7);
  set(legh,'FontSize',9);
  axis([min(res.t.date),max(res.t.date), -400,500]);

  if ( strcmp(per,'daily') )
    datetick('x',3,'keeplimits');
  end;

  titlename([upper(stn.station_name),': ', ...
             strrep(bdTfld,'_','\_'),' ',upper(per),' ',upper(char(sumfun)),' climatology ', ...
             stn.commentstr,' (',num2str(get_year(t.date(1))),'-',num2str(get_year(t.date(end))),')']);

  if ( nargout < 1 )
    res = []; clear res;
  end;

  %DEBUG:  toc,

return;
