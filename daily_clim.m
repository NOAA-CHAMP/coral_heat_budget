function daily_clim(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%function daily_clim(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%
% Plot three-year-at-a-time climatologies of sea temperature and heat budget
%
% Last Saved Time-stamp: <Wed 2012-03-07 12:47:12  lew.gramer>

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;


  [rawt,rawq,rawqe] = intersect_tses(stn.(sfld),stn.(hcdTdt),stn.([hcdTdt,'_err']));

  [cum,tid] = grp_ts(rawt.data,rawt.date,@floor,@nanmean,24);
  dt.data = cum;
  dt.date = tid;
  [cum,tid] = grp_ts(rawq.data,rawq.date,@floor,@nansum,24);
  dq.data = cum;
  dq.date = tid;
  [cum,tid] = grp_ts(rawqe.data,rawqe.date,@floor,@nansum,24);
  dqe.data = cum;
  dqe.date = tid;

  dyr=3;
  for yr=2000:dyr:2011
    yrs=yr:(yr+dyr-1);
    fun = @(x)(find(ismember(get_year(x.date),yrs)));
    yrsstr = [num2str(yrs(1)),'-',num2str(yrs(end))];

    t = subset_ts(dt,fun);
    q = subset_ts(dq,fun);
    qe = subset_ts(dqe,fun);

    if ( isempty(t.date) )
      disp(['Skipping ',num2str(yrs)]);
      continue;
    end;

    % Plot daily climatology
    [cum,tid] = grp_ts(t.data,t.date,'daily',@nanmean,0);
    climt.data = cum;
    climt.date = tid;

    fmg;

    lhs=[]; legs={};
    lhs(end+1) = plot_ts(climt,'k','LineWidth',3);
    legs{end+1} = ['T_s ',yrsstr];

    c = {'b','r','m','c','g'};
    for subyrix=1:length(yrs)
      subyr = yrs(subyrix);
      yrix = find(get_year(t.date)==subyr);
      if ( isempty(yrix) )
      else
        t0 = t.data(yrix(1));
        sq.date = q.date(yrix);
        sq.data = t0 + nancumsum(q.data(yrix)) - q.data(yrix(1));

        sq_minus_err.date = qe.date(yrix);
        sq_minus_err.data = sq.data - (qe.data(yrix).*24);
        sq_plus_err.date = qe.date(yrix);
        sq_plus_err.data = sq.data + (qe.data(yrix).*24);

        lhs(end+1) = plot(get_yearday(sq.date),sq.data,[c{subyrix},'-']);
        legs{end+1} = ['\Sigma\partial_tT \pm error ',num2str(subyr)];

        plot(get_yearday(sq_minus_err.date),sq_minus_err.data,[c{subyrix},':']);
        plot(get_yearday(sq_plus_err.date),sq_plus_err.data,[c{subyrix},':']);
      end; %if isempty(yrix) else
    end; %for subyrix=1:length(yrs)

    datetick3('x',3);
    titlename([upper(stn.station_name),' Daily Clim ',yrsstr]);
    legend(lhs,legs, 'Location','South');
    xlim(climt.date([1 end]));
    ylim([16,34]);
    ylabel('^oC');
    appendtitlename(stn.commentstr);
    appendtitlename([' (' num2str(yrs(1)) '-' num2str(yrs(end)) ')']);

  end; %for yr=2000:dyr:2011

return;
