function s = calc_heat_budget_rmse(stn,cf,gp,minN)
%function s = calc_heat_budget_rmse(stn,cf,gp,minN)
%
% Calculate cumulative values for sea temperature (@NANMEAN) and various heat
% budget terms (@NANSUM) by calling GRP_TS (v.) with accumulator function CF
% (DEFAULT: @FLOOR), max gap GP (DEFAULT: 1 d) and minimum period  data count
% MINN (DEFAULT: 24 - for hourly values). Estimate Root Mean Squared Error
% between sea temperature variability and heat budget terms; display results
% on Command Window by calling DUMP_HEAT_BUDGET_RMSE. Returns struct S with
% accumulated fields for basic heat budget terms, sea temperature, and RMSE.
% RMSE is estimated both between change in cumulative (DEFAULT daily mean)
% sea temperature and budget terms, and between climatological (DEFAULT
% year-day mean) sea temperature and cumulative budget terms.
%
% Last Saved Time-stamp: <Sat 2013-04-06 00:03:25 Eastern Daylight Time gramer>

  if ( ~exist('stn','var') || ~isfield(stn,'optim') || ~isfield(stn,'sfld') || ~isfield(stn,stn.sfld) )
    error('First arg STN must contain a heat budget result');
  end;

  % Having saved these fieldnames in OPTIMIZE_STATION_HEAT_BUDGET, we are
  % saved from having to specify all the *PFX args to this function again.
  s.sfld  = stn.sfld;
  s.qfld  = stn.qfld;
  s.bqfld = stn.bqfld;
  s.dtfld = stn.dtfld;
  s.hcfld = stn.hcfld;

  % Save in return struct S anything we may need to reproduce these estimates
  s.station_name = stn.station_name;
  s.optim        = stn.optim;
  s.(s.sfld)     = stn.(s.sfld);
  s.(s.qfld)     = stn.(s.qfld);
  s.(s.bqfld)    = stn.(s.bqfld);
  s.(s.dtfld)    = stn.(s.dtfld);
  s.(s.hcfld)    = stn.(s.hcfld);

  if ( ~exist('cf','var') || isempty(cf) )
    cf = @floor; % Daily
    %cf = @get_yeartriad; % Average every three days
    %cf = @get_yearpentad; % Average every five days
    %cf = @get_yearweek; % Weekly
  end;
  if ( ~exist('gp','var') || strcmpi(gp,'default') )
    switch (lower(char(cf))),
     case 'get_yearhour',	gp=(1.1/24); ccf=@get_yearday;
     case 'floor',		gp=1.1;      ccf=@get_jday_no_leap;
     case 'get_yeartriad',	gp=3.1;      ccf=@get_triad;
     case 'get_yearpentad',	gp=5.1;      ccf=@get_pentad;
     case 'get_yearweek',	gp=7.1;      ccf=@get_week;
     case 'get_yearmonth',	gp=31.1;     ccf=@get_month;
     case 'get_yearseason',	gp=92.1;     ccf=@get_season;
     case 'get_year',		gp=366.1;    ccf=@get_jday_no_leap;
     otherwise,			gp=[];       ccf=@get_jday_no_leap;
    end;
  end;
  if ( ~exist('ccf','var') || isempty(ccf) )
    ccf=@get_jday_no_leap;
  end;
  if ( ~exist('minN','var') || strcmpi(minN,'default') )
    switch (lower(char(cf))),
     case 'get_yearhour',	minN=[]; % We have no guess of the sampling frequency!
     case 'floor',		minN=23;
     case 'get_yeartriad',	minN=(24*1)+(23*2);
     case 'get_yearpentad',	minN=(24*2)+(23*3);
     case 'get_yearweek',	minN=(24*3)+(23*4);
     case 'get_yearmonth',	minN=23*27;
     case 'get_yearseason',	minN=23*87;
     case 'get_year',		minN=23*340;
     otherwise,			minN=[];
    end;
  end;

  s.cumfun     = cf;
  s.climcumfun = ccf;
  s.maxGap     = gp;
  s.minN       = minN;

  disp([char(cf),' & ',char(ccf)]);


  %% Calculate running means and sums (DEFAULT: daily)

  [ts,q,bq,dt,hc] = intersect_tses(s.(s.sfld),s.(s.qfld),s.(s.bqfld),s.(s.dtfld),s.(s.hcfld));

  [s.raw_ts.data,s.raw_ts.date] = grp_ts(ts.data,ts.date,cf,@nanmean,minN);
  if ( isempty(gp) )
    gp = 1.1*median(diff(s.raw_ts.date));
  end;
  if ( isempty(minN) )
    minN = 24*gp;
    s.minN       = minN;
    [s.raw_ts.data,s.raw_ts.date] = grp_ts(ts.data,ts.date,cf,@nanmean,minN);
  end;
  %disp({'DEBUG ',gp,minN});
  s.raw_td.date=s.raw_ts.date(2:end);
  s.raw_td.data=diff(s.raw_ts.data);
  s = filter_gaps(s,'raw_ts','raw_td',gp);

  % [s.raw_q.data,s.raw_q.date] = grp_ts(q.data,q.date,cf,@nansum,minN);
  % [s.raw_bq.data,s.raw_bq.date] = grp_ts(bq.data,bq.date,cf,@nansum,minN);
  % [s.raw_dt.data,s.raw_dt.date] = grp_ts(dt.data,dt.date,cf,@nansum,minN);
  % [s.raw_hc.data,s.raw_hc.date] = grp_ts(hc.data,hc.date,cf,@nansum,minN);

  % Use 24*NANMMEAN rather than NANSUM: makes days with <24 hours useable
  [cum,tid] = grp_ts(q.data,q.date,cf,@nanmean,minN);
  s.raw_q.data = cum*24; s.raw_q.date = tid;
  [cum,tid] = grp_ts(bq.data,bq.date,cf,@nanmean,minN);
  s.raw_bq.data = cum*24; s.raw_bq.date = tid;
  [cum,tid] = grp_ts(dt.data,dt.date,cf,@nanmean,minN);
  s.raw_dt.data = cum*24; s.raw_dt.date = tid;
  [cum,tid] = grp_ts(hc.data,hc.date,cf,@nanmean,minN);
  s.raw_hc.data = cum*24; s.raw_hc.date = tid;

  % [s.ts,s.td,s.q,s.bq,s.dt,s.hc]=intersect_tses(s.raw_ts,s.raw_td,s.raw_q,s.raw_bq,s.raw_dt,s.raw_hc);
  [s.ts,s.q,s.bq,s.dt,s.hc]=intersect_tses(s.raw_ts,s.raw_q,s.raw_bq,s.raw_dt,s.raw_hc);
  [ig,s.td]=intersect_tses(s.raw_ts,s.raw_td);

  s.qse =ts_fun(ts_op(s.q,s.td,'-'),@(x)(x.^2));
  s.bqse=ts_fun(ts_op(s.bq,s.td,'-'),@(x)(x.^2));
  s.dtse=ts_fun(ts_op(s.dt,s.td,'-'),@(x)(x.^2));
  s.hcse=ts_fun(ts_op(s.hc,s.td,'-'),@(x)(x.^2));

  s.dtbqse=ts_op(s.dtse,s.bqse,'-');

  s.N     =numel(s.qse.data);
  s.qrmse =sqrt(sum(s.qse.data)/s.N);
  s.bqrmse=sqrt(sum(s.bqse.data)/s.N);
  s.dtrmse=sqrt(sum(s.dtse.data)/s.N);
  s.hcrmse=sqrt(sum(s.hcse.data)/s.N);


  %% Calculate climatology (default, year-day)

  [s.raw_tsc.data,s.raw_tsc.date] = grp_ts(s.raw_ts.data,s.raw_ts.date,ccf,@nanmean,[]);
  [s.raw_qc.data,s.raw_qc.date] = grp_ts(s.raw_q.data,s.raw_q.date,ccf,@nanmean,[]);
  [s.raw_bqc.data,s.raw_bqc.date] = grp_ts(s.raw_bq.data,s.raw_bq.date,ccf,@nanmean,[]);
  [s.raw_dtc.data,s.raw_dtc.date] = grp_ts(s.raw_dt.data,s.raw_dt.date,ccf,@nanmean,[]);
  [s.raw_hcc.data,s.raw_hcc.date] = grp_ts(s.raw_hc.data,s.raw_hc.date,ccf,@nanmean,[]);

  [s.tsc,s.qcx,s.bqcx,s.dtcx,s.hccx]=intersect_tses(s.raw_tsc,s.raw_qc,s.raw_bqc,s.raw_dtc,s.raw_hcc);
  s.qc.date = s.qcx.date;   s.qc.data = s.tsc.data(1) + cumsum(s.qcx.data) - s.qcx.data(1);
  s.bqc.date = s.bqcx.date; s.bqc.data = s.tsc.data(1) + cumsum(s.bqcx.data) - s.bqcx.data(1);
  s.dtc.date = s.dtcx.date; s.dtc.data = s.tsc.data(1) + cumsum(s.dtcx.data) - s.dtcx.data(1);
  s.hcc.date = s.hccx.date; s.hcc.data = s.tsc.data(1) + cumsum(s.hccx.data) - s.hccx.data(1);

  s.qsec  =ts_fun(ts_op(s.qc,s.tsc,'-'),@(x)(x.^2));
  s.bqsec =ts_fun(ts_op(s.bqc,s.tsc,'-'),@(x)(x.^2));
  s.dtsec =ts_fun(ts_op(s.dtc,s.tsc,'-'),@(x)(x.^2));
  s.hcsec =ts_fun(ts_op(s.hcc,s.tsc,'-'),@(x)(x.^2));

  s.Nc      =numel(s.qsec.data);
  s.qrmsec  =sqrt(sum(s.qsec.data)/s.Nc);
  s.bqrmsec =sqrt(sum(s.bqsec.data)/s.Nc);
  s.dtrmsec =sqrt(sum(s.dtsec.data)/s.Nc);
  s.hcrmsec =sqrt(sum(s.hcsec.data)/s.Nc);


  %% Display results to command window (or CSV file, figures, or other?)

  s = dump_heat_budget_rmse(s);

  pause(0.5);

return;
