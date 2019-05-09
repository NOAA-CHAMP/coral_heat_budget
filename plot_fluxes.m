function fh = plot_fluxes(stn,yr,begdy,enddy,relflds,absflds,accflds,fh,legends)
%function fh = plot_fluxes(stn,yr,begdy,enddy,relflds,absflds,accflds,fh,legends)
%
% Plot temperature or similar field(s) RELFLDS vs. multiple estimates of
% cumulative heat flux ACCFLDS, T0+sum(dT/dt), over the time period from
% DATENUM(YR,1,0)+BEGDY to DATENUM(YR,1,0)+ENDDY. Additional non-relative
% fields ABSFLDS{:} can also be plotted. Returns figure handle used for all
% plots (value in FH, or DEFAULT: new Figure). Show a LEGEND using LEGENDS
% (CELLSTR) if specified, otherwise using all fieldnames passed in.
%
% Last Saved Time-stamp: <Sat 2010-12-11 15:20:42 Eastern Standard Time gramer>


  if ( ~exist('fh','var') || isempty(fh) )
    fh = figure;
    hold on;
  elseif ( ishandle(fh) )
    figure(fh);
    hold on;
  else
    error('Arg FH must either be missing, empty, or a valid figure handle!');
  end;

  if ( ~exist('legends','var') || isempty(legends) )
    legends = [];
  end;

  % How many hours between Markers in each plot?
  ndy = enddy - begdy;
  if ( ndy > 10*365 )
    dt = 8*7*24;
  elseif ( ndy > 365 )
    dt = 2*7*24;
  elseif ( ndy > 90 )
    dt = 7*24;
  elseif ( ndy > 30 )
    dt = 2*24;
  elseif ( ndy > 14 )
    dt = 24;
  else
    dt = 12;
  end;


  % Fields to inter-compare

  % Relative (T - T_0) time series evolution
  if ( ~exist('relflds','var') || isempty(relflds) || (~iscell(relflds)&&strcmpi(relflds,'default')) )
    relflds = { ...
        'ndbc_sea_t', ...
              };

        % 'ndbc_sea_t_12_hour_lowpass', ...
        % 'ndbc_sea_t', ...
        % 'ndbc_air_t', ...
  end;

  % Absolute time series
  if ( ~exist('absflds','var') || (~iscell(absflds)&&strcmpi(absflds,'default')) )
    absflds = { ...
        'ndbc_ncep_30a_wind_stress_30_day_maximum', ...
        'gom_hycom_speed_7_day_maximum', ...
              };

        % 'gom_hycom_speed', ...
        % 'gom_hycom_u', ...
        % 'ndbc_sea_t_1_day_deviation_3_day_average', ...
        % 'ndbc_wind1_u_7_day_deviation_sum_ndbc_wind1_v_7_day_deviation', ...

        % 'ndbc_wind1_u_3_day_deviation_sum_ndbc_wind1_v_3_day_deviation', ...
        % 'ndbc_wind1_u_7_day_deviation_sum_ndbc_wind1_v_7_day_deviation', ...
        % 'ndbc_wind1_u', ...
        % 'ndbc_wind1_v', ...
        % 'ndbc_wind1_v_24_hour_lowpass', ...
        % 'ndbc_wind1_v_3_day_maximum', ...
        % 'stokes_drift_v_3_day_maximum', ...
        % 'stokes_drift_speed_3_day_maximum', ...
        % 'ndbc_ncep_30a_wind_stress_30_day_maximum', ...
        % 'ww3_sigwavehgt', ...
        % 'ww3_sigwavehgt_30_day_maximum', ...
        % 'quasi_eulerian_speed_30_day_maximum', ...
        % 'ndbc_ncep_30a_wind_stress_30_day_maximum', ...

        % 'ndbc_wind1_speed_30_day_maximum', ...
        % 'ndbc_wind1_u_3_day_average', ...
        % 'ndbc_wind1_v_3_day_average', ...
        % 'ndbc_ncep_30a_wind_stress', ...
  end;

  % Cumulative (Sum(H)) time series evolution
  if ( ~exist('accflds','var') || isempty(accflds) || (~iscell(accflds)&&strcmpi(accflds,'default')) )
    % accflds = { ...
    %     'ndbc_ncep_26_heat_flux_term', ...
    %     'ndbc_ncep_26_dt', ...
    %     'ndbc_ncep_30_heat_flux_term', ...
    %     'ndbc_ncep_30_dt', ...
    %     'ndbc_ncep_30a_heat_flux_term', ...
    %     'ndbc_ncep_30a_dt', ...
    %     'bic_ncep_26_heat_flux_term', ...
    %     'bic_ncep_26_dt', ...
    %           };
    accflds = { ...
        'gom_hycom_netqf', ...
        'gom_hycom_dt_netqf', ...
        'netqf', ...
              };

        % 'ndbc_ncep_30a_heat_flux_term', ...
        % 'gom_hycom_dt', ...

        % 'gom_hycom_qvf', ...
        % 'gom_hycom_qf', ...
        % 'ndbc_ncep_30a_total_dt', ...
  end;


  % linspec = {'k.-','k.:','k*-','k*:','ko-','ko:','ks-','ks:','k^-','k^:','kp-','kp:','kd-','kd:','kh-','kh:','kv-','kv:'};
  %DEBUG:  linspec = {'k-','k:','k.-','k.:','k*-','k*:','ko-','ko:','ks-','ks:','k^-','k^:','kp-','kp:','kd-','kd:','kh-','kh:','kv-','kv:'};

  % Maximum visibility in PowerPoints
  %%%% linspec = {'k-','r.-','b*-','go-','ms-','c^-','kp-','rd-','bh-','gv-','m-','c.-',};
  linspec = {'k-','r.-','b*-','go-','ms-','c^-','yp-','kd-','rh-','bv-','gx-','m+-','c<-','y>-'};

  % Line handles we want to show legends for
  lh = [];

  % cncspec = {'k-','k:',};
  % mrkspec = {'k-','k:','k.','k.','k*','k*','ko','ko','ks','ks','k^','k^','kp','kp','kd','kd','kh','kh','kv','kv'};


  begdt = datenum(yr,1,1) + begdy - 1;
  enddt = datenum(yr,1,1) + enddy - 1;

  cix = 1;
  legs = {};
  for ix = 1:length(relflds)
    stn = verify_variable(stn,relflds{ix});
    if (~isfield(stn,relflds{ix}))
      warning('No field "%s" could be found or calculated!', relflds{ix});
      continue;
    end;
    fld = stn.(relflds{ix});
    dtix = find(begdt <= fld.date & fld.date < enddt);
    if (isempty(dtix))
      warning('No matching dates in field "%s"!', relflds{ix});
      continue;
    end;
    legs = { legs{:} relflds{ix} };
    dts = fld.date(dtix);
    dat = fld.data(dtix);
    % Make *all* RELFLDS fields relative to the *first* field (e.g., sea temp.)
    if ( ix == 1 )
      firstix = find(isfinite(dat), 1);
      if (isempty(firstix)); firstix=1; end;  % It's all bad??
                                              % val0 = dat(firstix);
      firstdts = dts(firstix:end);
      firstdat = dat(firstix:end);
    end;
% %     %DEBUG:
% %     % plot(dts,real(dat-val0),linspec{mod(cix-1,length(linspec))+1});
% %     plot(dts,real(dat),linspec{mod(cix-1,length(linspec))+1});
%     plot(dts(1:dt:end),real(dat(1:dt:end)),linspec{mod(cix-1,length(linspec))+1},'LineWidth',1.5);
    lh(end+1)=...
    plot(dts(1       ),real(dat(1       )),linspec{mod(cix-1,length(linspec))+1});
    plot(dts(1: 1:end),real(dat(1: 1:end)),linspec{mod(cix-1,length(linspec))+1},'LineWidth',1.5,'Marker','none');
    plot(dts(1:dt:end),real(dat(1:dt:end)),linspec{mod(cix-1,length(linspec))+1},'LineWidth',1.5,'LineStyle','none');
    cix = cix + 1;
  end;
  for ix = 1:length(absflds)
    stn = verify_variable(stn,absflds{ix});
    if (~isfield(stn,absflds{ix}))
      warning('No field "%s" could be found or calculated!', absflds{ix});
      continue;
    end;
    fld = stn.(absflds{ix});
    dtix = find(begdt <= fld.date & fld.date < enddt);
    if (isempty(dtix))
      warning('No matching dates in field "%s"!', absflds{ix});
      continue;
    end;
    legs = { legs{:} absflds{ix} };
    dts = fld.date(dtix);
    dat = fld.data(dtix);
% %     plot(dts,real(dat),linspec{mod(cix-1,length(linspec))+1});
%     plot(dts(1:dt:end),real(dat(1:dt:end)),linspec{mod(cix-1,length(linspec))+1});
    lh(end+1)=...
    plot(dts(1       ),real(dat(1       )),linspec{mod(cix-1,length(linspec))+1});
    plot(dts(1: 1:end),real(dat(1: 1:end)),linspec{mod(cix-1,length(linspec))+1},'Marker','none');
    plot(dts(1:dt:end),real(dat(1:dt:end)),linspec{mod(cix-1,length(linspec))+1},'LineStyle','none');
    %%%% ??? HACK!
    % plot(dts,15+( 10*real(dat) ),linspec{mod(cix-1,length(linspec))+1});
    cix = cix + 1;
  end;
  for ix = 1:length(accflds)
    stn = verify_variable(stn,accflds{ix});
    if (~isfield(stn,accflds{ix}))
      warning('No field "%s" could be found or calculated!', accflds{ix});
      continue;
    end;
    fld = stn.(accflds{ix});
    dtix = find(begdt <= fld.date & fld.date < enddt);
    if (isempty(dtix))
      warning('No matching dates in field "%s"!', accflds{ix});
      continue;
    end;
    legs = { legs{:} accflds{ix} };
    dts = fld.date(dtix);
    dat = fld.data(dtix);
    dts(~isfinite(dat)) = [];
    dat(~isfinite(dat)) = [];
    [rix,aix] = intersect_dates(firstdts,dts);
    dts = dts(aix);
    dat = dat(aix);
    % % Do cumulative sums only between any big gaps in time series
    % bigd = [ 1 ; find(diff(dts) > 14) ; (length(dat)+1) ];
    % Do cumulative sums only between any big gaps in time series date OR data
    bigd = [ 1 ; find(diff(dts) > 14 | abs(diff(dat)) > (3*std(dat))) ; (length(dat)+1) ];
    %DEBUG:    bigd=[1;find(diff(get_year(dts))~=0);(length(dat)+1)]; % New Years reset only
    %DEBUG:
    bigd = [1 ; (length(dat)+1)]; % NO gap resets!
    for segix = 1:(length(bigd)-1)
      dat( bigd(segix):(bigd(segix+1)-1) ) = ...
          firstdat( rix(bigd(segix)) ) + ...
          cumsum( dat( bigd(segix):(bigd(segix+1)-1) ) );
    end;
% %     plot(dts,real(dat),linspec{mod(cix-1,length(linspec))+1});
%     plot(dts(1:dt:end),real(dat(1:dt:end)),linspec{mod(cix-1,length(linspec))+1});
    lh(end+1)=...
    plot(dts(1       ),real(dat(1       )),linspec{mod(cix-1,length(linspec))+1});
    plot(dts(1: 1:end),real(dat(1: 1:end)),linspec{mod(cix-1,length(linspec))+1},'Marker','none');
    plot(dts(1:dt:end),real(dat(1:dt:end)),linspec{mod(cix-1,length(linspec))+1},'LineStyle','none');
    cix = cix + 1;
  end;
  hold off;

  maxigraph;
  set(gca,'position',[0.05 0.05 0.90 0.90]);
  datetick3;

  legs = strrep(legs,'_','\_');
  if ( isempty(legends) )
    legends = legs;
  end;
  % legend(legends, 'Location','Best');
  legend(lh,legends, 'Location','Best');

  % Very light grey
  set(gca,'Color',[.9 .9 .9]);

  titlename([ stn.station_name ': ' datestr(begdt,1) ' - ' datestr(enddt,1) ]);
  grid on;

  xlim(dts([1 end]));
%%%%%%%%%% ??? DEBUG
%%%% ylim([-20 20]);
%%%%%%%%%% ??? DEBUG
ylabel('Sea temperature T vs. estimates of T(0) + \Sigma(\delta_tT) [^oC]');

return;
