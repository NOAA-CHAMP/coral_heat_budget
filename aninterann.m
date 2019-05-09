function aninterann

  datapath = get_thesis_path('../data');
  figspath = get_thesis_path('../figs');

  stnms = { 'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1' };
  stnms = { 'smkf1' };
  stnms = { 'lonf1' };

  for stix=1:numel(stnms)
    stnm = lower(stnms{stix});
    stn = load_all_ndbc_data([],stnm);

    switch (stnm)
     case 'smkf1',
      ix = find(datenum(1998,7,14)<stn.ndbc_sea_t.date&stn.ndbc_sea_t.date<datenum(1998,7,18));
      stn.ndbc_sea_t.date(ix) = [];
      stn.ndbc_sea_t.data(ix) = [];
    end;

    % % Only analyze "warm" months (MJJASO)
    % stn.ndbc_sea_t.data = stn.ndbc_sea_t.data(ts_boreal_warm(stn.ndbc_sea_t));
    % stn.ndbc_sea_t.date = stn.ndbc_sea_t.date(ts_boreal_warm(stn.ndbc_sea_t));
    % % Only analyze months AMJJASOND
    % stn.ndbc_sea_t.data = stn.ndbc_sea_t.data(ismember(get_month(stn.ndbc_sea_t.date),[4:12]));
    % stn.ndbc_sea_t.date = stn.ndbc_sea_t.date(ismember(get_month(stn.ndbc_sea_t.date),[4:12]));

    stn.t.date = [];
    stn.t.data = [];

    ts = [];
    legs = {};

    yrs = get_year(stn.ndbc_sea_t.date);

    allyrs = 1987:2011;
    for yrix=1:numel(allyrs)
      yr = allyrs(yrix);
      ix = find(yrs == yr);
      if ( min(yrs) > yr || yr > max(yrs) || numel(ix) < ((365-45)*24) )
      % if ( min(yrs) > yr || yr > max(yrs) || numel(ix) < ((365-135)*24) )
      % if ( min(yrs) > yr || yr > max(yrs) || numel(ix) < ((183-30)*24) )
      % if ( min(yrs) > yr || yr > max(yrs) || numel(ix) < (260*24) )
        % Ignore missing leap years - does not matter
        stn.t.date(end+1:end+8760,1) = datenum(yr,1,1):(1/24):(datenum(yr,1,1)+364.99);
        stn.t.data(end+1:end+8760,1) = repmat(0,[8760,1]);
      else
        stn.t.date(end+1:end+numel(ix),1) = stn.ndbc_sea_t.date(ix);
        stn.t.data(end+1:end+numel(ix),1) = stn.ndbc_sea_t.data(ix);
        ts(end+1).date = get_yearday(stn.ndbc_sea_t.date(ix));
        ts(end).data = stn.ndbc_sea_t.data(ix);
        legs{end+1} = num2str(yr);
      end;
    end;

    if (0)
    station_anova_multcompare(stn,'t',[],'^oC',[1,23],[1,0,0,0]);
    view(270,90);
    xlim([25,28]);
    % xlim([27,31]);
    % xlim([25,31]);
    ylim([1,23]);

    ytks = get(gca,'YTick');
    ytls = get(gca,'YTickLabel');
    set(gca,'YTick',ytks(1:2:end));
    set(gca,'YTickLabel',ytls(1:2:end));
    set(gca,'FontSize',13);

    print('-dpng',fullfile(figspath,[stnm,'-ndbc_sea_t-interannual-multcompare.png']));
    print('-dtiff',fullfile(figspath,[stnm,'-ndbc_sea_t-interannual-multcompare.tiff']));
    end;

    if (0)
    fmg;
    lhs=plot_ts(ts);
    % legend(legs,'Location','NorthEastOutside');
    legend(legs,'Location','South');
    set(lhs,'Marker','none');
    datetick3('x',3,'keeplimits');
    xlim([0,366]);
    h = datacursormode(gcf);
    set(h,'UpdateFcn',@yearday_select_cb);
    end;

    stn.t = stn.ndbc_sea_t;
    stn.t.date(stn.t.data==0) = [];
    stn.t.data(stn.t.data==0) = [];
    % stn.t.data(get_year(stn.t.date)==2007) = [];
    % stn.t.date(get_year(stn.t.date)==2007) = [];
    [B,Stats,fh,lh] = scatter_fit(stn.t.date,stn.t.data);
    datetick3('x',10,'keeplimits');
    ylim([15,35]);
    {numel(Stats.w),Stats.se(2),Stats.t(2),Stats.p(2),},
    set(lh(1),'Marker','.','MarkerSize',1,'LineStyle','none','LineWidth',1,'Color',[.3,.3,.3]);
    set(lh(2),'Marker','none','LineStyle','-','LineWidth',3,'Color','k');

    if (0)
    [cum,tid] = grp_ts(stn.t.data,stn.t.date,@get_year,@nanmean,0);
    [B,Stats] = scatter_fit([1:numel(tid)]',cum);
    ylim([25,28]);
    {numel(Stats.w),Stats.se(2),Stats.t(2),Stats.p(2),},
    end;
    

    stn.ndbc_tide.data(get_year(stn.ndbc_tide.date)>=2011) = [];
    stn.ndbc_tide.date(get_year(stn.ndbc_tide.date)>=2011) = [];
    [B,Stats,fh,lh] = scatter_fit(stn.ndbc_tide.date,stn.ndbc_tide.data);
    datetick3('x',10,'keeplimits');
    ylim([-1.5,+3.5]);
    {numel(Stats.w),Stats.se(2),Stats.t(2),Stats.p(2),},
    set(lh(1),'Marker','.','MarkerSize',1,'LineStyle','none','LineWidth',1,'Color',[.3,.3,.3]);
    set(lh(2),'Marker','none','LineStyle','-','LineWidth',3,'Color','k');

    stn = []; clear stn;
  end;

return;
