function annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,commentstr,mos,begyr,endyr)
%function annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,commentstr,mos,begyr,endyr)
%
% Plot yearly cumulative comparison (one subplot per year of sea temperature)
% for various forms (intermediate results) of the reef ocean heat budget.
% Args STN to SUBSTITUTE_FIELD_NAMES are described in STATION_HEAT_BUDGET (v.).
% COMMENTSTR (DEFAULT: STN.commentstr or '') is a string appended to plot
% title; MOS, BEGYR, ENDYR constrain which months of which years are plotted;
% finally STN.(COMPARISONSFLD) plots a sea temperature other than STN.(SFLD)
% for comparison with the heat budget results. NOTE: Resets on any gap>=7d.
%
% Last Saved Time-stamp: <Tue 2012-07-31 16:34:30  lew.gramer>


  if ( exist('commentstr','var') && ~isempty(commentstr) )
    % Hide the user's comment string away in our STN struct, safe from
    % manipulation by script STATION_HEAT_BUDGET_FIELD_NAMES (below).
    stn.commentstr = commentstr;
  end;
  if ( ~isfield(stn,'commentstr') )
    stn.commentstr = '';
  end;


  station_heat_budget_field_names;

  if ( ~exist('mos','var') || isempty(mos) )
    mos=1:12;
  end;

  if ( ~exist('comparisonsfld','var') || isempty(comparisonsfld) )
    comparisonsfld = sfld;
  end;


  % What to accumulate: air-sea, -unabsorbed, +benthic, +km-scale, +horizontal-convection

  q0f = sqtfld;
  % q0f = qtfld;
  bq0f = bq0tfld;
  % dtf = dTfld;
  dtf = bdTfld;
  qf = hcdTdt;

  if ( ~exist('begyr','var') || isempty(begyr) )
    begyr = max([get_year(stn.(q0f).date(1)),get_year(stn.(bq0f).date(1)),get_year(stn.(dtf).date(1)),get_year(stn.(qf).date(1))]);
  end;
  if ( ~exist('endyr','var') || isempty(endyr) )
    % Last complete year of data (presumably)
    endyr = get_year(now) - 1;
  end;

  nyrs = numel(begyr:endyr);
  switch (nyrs)
   case {1,2},
    nrows = 1;
    ncols = nyrs;
   case {3,4},
    nrows = 2;
    ncols = 2;
   case {5,6},
    nrows = 3;
    ncols = 2;
   case {7,8,9},
    nrows = 3;
    ncols = 3;
   case {10,11,12},
    nrows = 4;
    ncols = 3;
   case {13,14,15,16},
    nrows = 4;
    ncols = 4;
   case {17,18,19,20},
    nrows = 5;
    ncols = 4;
   otherwise,
    nrows = 5;
    ncols = 4;
    begyr = endyr - 19;
  end;


  yrs=begyr:endyr;

  disp(['Comparing ' comparisonsfld ' with ' qf]);

  fmg;
  suptitlename([ upper(stn.station_name) ' ' ...
                 strrep(comparisonsfld,'_','\_') ' vs ' strrep(qf,'_','\_') ...
                 ' ' num2str([min(mos),max(mos)],'%g-') ' ' stn.commentstr]);

  [cix,qix,bix,dix,tix] = intersect_all_dates([],...
                                              stn.(comparisonsfld).date,...
                                              stn.(q0f).date,...
                                              stn.(bq0f).date,...
                                              stn.(dtf).date,...
                                              stn.(qf).date);
  stn.(comparisonsfld).date = stn.(comparisonsfld).date(cix);
  stn.(comparisonsfld).data = stn.(comparisonsfld).data(cix);
  stn.(q0f).date = stn.(q0f).date(qix);
  stn.(q0f).data = stn.(q0f).data(qix);
  stn.(bq0f).date = stn.(bq0f).date(bix);
  stn.(bq0f).data = stn.(bq0f).data(bix);
  stn.(dtf).date = stn.(dtf).date(dix);
  stn.(dtf).data = stn.(dtf).data(dix);
  stn.(qf).date = stn.(qf).date(tix);
  stn.(qf).data = stn.(qf).data(tix);

  %%%% ??? DEBUG
  miny=floor(nanmin(stn.(comparisonsfld).data)-2);
  maxy=ceil(nanmax(stn.(comparisonsfld).data)+2);
  %%%% ??? DEBUG


  for yrix=1:length(yrs)
    yr = yrs(yrix);
    mindt = datenum(yr,mos(1),1);
    maxdt = datenum(yr,mos(end)+1,1);

    dtix = find(mindt<=stn.(comparisonsfld).date & stn.(comparisonsfld).date<=maxdt);
    t.date = stn.(comparisonsfld).date(dtix);
    t.data = stn.(comparisonsfld).data(dtix);

    dtix = find(mindt<=stn.(q0f).date & stn.(q0f).date<=maxdt);
    q0.date = stn.(q0f).date(dtix);
    q0.data = stn.(q0f).data(dtix);

    dtix = find(mindt<=stn.(bq0f).date & stn.(bq0f).date<=maxdt);
    bq0.date = stn.(bq0f).date(dtix);
    bq0.data = stn.(bq0f).data(dtix);

    dtix = find(mindt<=stn.(dtf).date & stn.(dtf).date<=maxdt);
    dt.date = stn.(dtf).date(dtix);
    dt.data = stn.(dtf).data(dtix);

    dtix = find(mindt<=stn.(qf).date & stn.(qf).date<=maxdt);
    q.date = stn.(qf).date(dtix);
    q.data = stn.(qf).data(dtix);

    if ( ~isempty(t.data(isfinite(t.data))) && ~isempty(q0.data(isfinite(q0.data))) && ~isempty(bq0.data(isfinite(bq0.data))) && ~isempty(dt.data(isfinite(dt.data))) && ~isempty(q.data(isfinite(q.data))) )

      subplot_tight(nrows,ncols,yrix);
      hold on;

      % Assume equilibrium is achieved again after any time series gap >=7 days
      gapix = [1 ; (find(diff(q.date) >= 7)+1) ; length(q.date)+1];
      for gapixix = 1:length(gapix)-1
        ixes = gapix(gapixix):(gapix(gapixix+1)-1);

        q0s.date = q0.date(ixes);
        q0s.data = t.data(ixes(1)) + nancumsum(q0.data(ixes)) - q0.data(ixes(1));

        bq0s.date = bq0.date(ixes);
        bq0s.data = t.data(ixes(1)) + nancumsum(bq0.data(ixes)) - bq0.data(ixes(1));

        dts.date = dt.date(ixes);
        dts.data = t.data(ixes(1)) + nancumsum(dt.data(ixes)) - dt.data(ixes(1));

        qs.date = q.date(ixes);
        qs.data = t.data(ixes(1)) + nancumsum(q.data(ixes)) - q.data(ixes(1));

        plot(t.date(ixes),t.data(ixes),'k',q0s.date,q0s.data,'c-.',bq0s.date,bq0s.data,'y--',dts.date,dts.data,'r:',qs.date,qs.data,'b:', 'LineWidth',2);
      end;

      xlim([mindt,maxdt-(1/1000)]);
      if ( length(mos) < 12 )
        datetick('x',2,'keeplimits');
      else
        datetick('x',17,'keeplimits');
      end;
      set_datetick_cursor;
      grid on;

      % miny = nanmin([t.data;q0s.data;bq0s.data;dts.data;qs.data]);
      % maxy = nanmax([t.data;q0s.data;bq0s.data;dts.data;qs.data]);

      if ( nrows < 3 )
        lh = legend('T_s','Q_0/\rhoC_ph','(Q_0(\gamma)+Q_b)/\rhoC_ph','(Q_0(\gamma)+Q_b)/\rhoC_ph+u^.\nablaT+K\nabla^2T','\partial_tT_s',...
               'Location','Best');
        set(lh,'FontSize',7);
        % miny=13; maxy=33;
      else
        % % miny=10; maxy=40;
        % miny=18; maxy=33;
      end;

      %%%%DEBUG???
      miny = min([miny,12]);
      maxy = max([maxy,35]);

      ylim([miny,maxy]);

    end; %if ( ~isempty(t.data(isfinite(t.data))) && ...

  end; %for yrix=1:length(yrs)

  if ( nrows >= 3 )
    subplot_tight(nrows,ncols,1);
    lh = legend('T_s','Q_0/\rhoC_ph','(Q_0(\gamma)+Q_b)/\rhoC_ph','(Q_0(\gamma)+Q_b)/\rhoC_ph+u^.\nablaT+K\nabla^2T','\partial_tT_s',...
           'Location','Best');
    set(lh,'FontSize',7);
  end;


return;
