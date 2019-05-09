function stn = station_climatology_heat_budget(stn,CLIMPFX)
%function stn = station_climatology_heat_budget(stn,CLIMPFX)

  commentstr = 'CLIM ';

  KMPFX = 'gom_hycom';

  switch (CLIMPFX),
   case 'erai',		WAVEPFX = 'erai';
   case 'oaflux',	WAVEPFX = 'ww3';
   otherwise,		WAVEPFX = 'ww3';
  end;

  %%%
  %% Call SCRIPT to set:
  %% Variable-name prefixes ("PFX") for various input and output datasets; AND,
  %% All station struct fieldnames used to produce heat budget 
  station_heat_budget_field_names;

  % Climatologies may be daily or monthly - expand to hourly time series
  switch (CLIMPFX),
   case 'oaflux',
    for cfld = {climq0fld,climsrfld}
      fld = cfld{:};
      stn.(fld) = ts_daily_mean_to_hourly(stn.(['daily_' fld]));
    end;
  end;

  %%%
  %% Turbulent and Net Fluxes

  stn = station_heat_flux_term(stn,climq0fld,climqtfld,sfld,[],mhfld);


  %%%
  %% Eulerian (km + Stokes) Heat Advection, Km-scale Heat Diffusion

  stn = station_cross_shore_advection(stn,bathorifld,...
                                      qeufld,qevfld,Tfld,kmtfld,...
                                      ['raw_' udTfld],udTfld,...
                                      climqtfld,climqtAdvfld,'nearest');
  stn = station_calc_kdel2t(stn,K_theta,Tfld,...
                            ['raw_' kd2Tfld],kd2Tfld,...
                            climqtAdvfld,climdTfld,'nearest');
  stn = station_heat_flux_term_inverse(stn,climdTffld,...
                                       climdTfld,sfld,[],mhfld);


  %%%
  %% Benthic Heat Exchanges

  stn.(climbq0fld) = stn.(climq0fld);
  stn.(climbq0tfld) = stn.(climqtfld);
  stn.(climbdTfld) = stn.(climdTfld);
  stn.(climbdTffld) = stn.(climdTffld);


  %%%
  %% Horizontal Convection

    bet = stn.(slopefld);
    % stn = station_horizontal_convection(stn,sfld,[],mhfld,tspdfld,bq0fld,bet,hc_opts);

    [tix,hix,tspdix,aix,Wix,q0ix,bq0ix,dTix] = ...
      intersect_all_dates([],stn.(sfld).date,stn.(mhfld).date,stn.(tspdfld).date,stn.(afld).date,stn.(Wfld).date,stn.(climq0fld).date,stn.(climbq0fld).date,stn.(climbdTffld).date);
    t = stn.(sfld).data(tix);
    s = repmat(36,size(stn.(sfld).data(tix)));
    h = stn.(mhfld).data(hix);
    tspd = stn.(tspdfld).data(tspdix);
    at = stn.(afld).data(aix);
    W = stn.(Wfld).data(Wix);
    q0 = stn.(climq0fld).data(q0ix);
    bq0 = stn.(climbq0fld).data(bq0ix);
    dT = stn.(climbdTffld).data(dTix);

    qstr=climbdTffld;
    %%%% ??? DEBUG: Base horizontal convection on surface heating only!
    q=bq0; commentstr = [commentstr ' HC(Q0+Qb) '];
    %%%% ??? DEBUG: Base horizontal convection on total (km-scale) budget
    % q=dT;

    hc_opts.R = (1.00-0.08);

    hc_opts.scaling = 'SS';
    % hc_opts.scaling = 'US';
    % hc_opts.scaling = 'SU';
    % hc_opts.scaling = 'UU';
    hc_opts.maximum_onset_time = 6*3600;
    %%%% ??? DEBUG
    commentstr = [commentstr ' (HC ' hc_opts.scaling ' maxT:9h) '];
    res = horizontal_convection(t,s,h,q,bet,hc_opts);

    dts = stn.(sfld).date(tix);
    flds = fieldnames(res);
    for fldix = 1:length(flds)
      fld = flds{fldix};
      dat = res.(fld);
      res.(fld) = [];
      res.(fld).date = dts;
      res.(fld).data = dat;

      stnfld = [CLIMHCPFX '_' fld];
      stn.(stnfld).date = dts;
      stn.(stnfld).data = dat;
    end;

    %%%% ??? DEBUG
    % Truncate start of all time series for display purposes
    startdt = datenum(2004,1,1);
    % startdt = datenum(2005,1,1);
    % startdt = datenum(2005,1,23);

    startix = find(dts>=startdt,1);

    %%%% DEBUG???
    ix = find(~isfinite(q0(startix:end)),1,'last');
    if ( ~isempty(ix) )
      startix = startix + ix;
    end;

    dts = dts(startix:end);
    t = t(startix:end);
    h = h(startix:end);
    tspd = tspd(startix:end);
    at = at(startix:end);
    W = W(startix:end);
    q0 = q0(startix:end);
    bq0 = bq0(startix:end);
    dT = dT(startix:end);
    fac = stn.(climhcfactor).data(startix:end);
    dTdt.date = stn.(climhcdTdt).date(startix:end);
    dTdt.data = stn.(climhcdTdt).data(startix:end);

    bt = stn.(btfld);
    [ig,startix] = min(abs(bt.date-dts(1)));
    [ig,endix] = min(abs(bt.date-dts(end)));
    bt.date = bt.date(startix:endix);
    bt.data = bt.data(startix:endix);
    %%%% ??? DEBUG

    dsr_lpfld = [climsrfld '_24_hour_sum'];
    stn = verify_variable(stn,dsr_lpfld);
    [ig,dsr_lpix] = intersect_dates(dts,stn.(dsr_lpfld).date);

    % lhf_lpfld = [qlhfld '_24_hour_sum'];
    % stn = verify_variable(stn,lhf_lpfld);
    % [ig,lhf_lpix] = intersect_dates(dts,stn.(lhf_lpfld).date);

    udT_lpfld = [udTfld];
    stn = verify_variable(stn,udT_lpfld);
    [ig,udT_lpix] = intersect_dates(dts,stn.(udT_lpfld).date);

    %%%% ??? DEBUG
    T0 = t(1);
    bigfh=figure, maxigraph; hold on;
    plot(dts,t,dts,at,'k:',bt.date,bt.data,'m--',dts,T0+cumsum(q0.*fac),dts,T0+cumsum(bq0.*fac),dts,T0+cumsum(dT.*fac),dTdt.date,T0+cumsum(dTdt.data));
    % plot(stn.(dsr_lpfld).date(dsr_lpix),(stn.(dsr_lpfld).data(dsr_lpix)./1000),stn.(lhf_lpfld).date(lhf_lpix),-(stn.(lhf_lpfld).data(lhf_lpix)./1000),'y:','Color',[.8,.8,0]);
    % plot(stn.(dsr_lpfld).date(dsr_lpix),(stn.(dsr_lpfld).data(dsr_lpix)./1000),stn.(udT_lpfld).date(udT_lpix),T0+cumsum(stn.(udT_lpfld).data(udT_lpix)),'y:','Color',[.8,.8,0]);
    plot(stn.(dsr_lpfld).date(dsr_lpix),T0+(stn.(dsr_lpfld).data(dsr_lpix).*fac),stn.(udT_lpfld).date(udT_lpix),T0+cumsum(stn.(udT_lpfld).data(udT_lpix)),'y:','Color',[.8,.8,0]);
    % plot(dts,h,dts,tspd,'k:',dts,W,'o','Color',[0,.8,.2]);
    plot(dts,h,dts,tspd,'k:','Color',[0,.8,.2]);
    datetick3('x',2,'keeplimits');
    % legend('T_s','T_a','T_b','T_0+Q_0/\rhoC_ph','T_0+(Q_0+Q_b)/\rhoC_ph','T_0+\partial_tT_k_m',['T_0+\partial_tT ' hc_opts.scaling],'T_0+\Sigma_1_d Q_S_W\times1h/\rho^.C_p^.h','T_0+\Sigmau^.\nablaT_k_m','h_t_i_d_e','SPD_t_i_d_e','W', 'Location','SouthWest'); %'Best');
    legend('T_s','T_a','T_b','T_0+Q_0/\rhoC_ph','T_0+(Q_0+Q_b)/\rhoC_ph','T_0+\partial_tT_k_m',['T_0+\partial_tT ' hc_opts.scaling],'T_0+\Sigma_1_d Q_S_W\times1h/\rho^.C_p^.h','T_0+\Sigmau^.\nablaT_k_m','h_t_i_d_e','SPD_t_i_d_e', 'Location','SouthWest'); %'Best');
    titlename([ commentstr stn.station_name ' ' strrep(qstr,'_','\_') ]);
    disp(commentstr);
    grid on;
    % axis('tight');
    % axis('tight'); ylim([-50 150]);
    % axis([startdt,startdt+365,0,60]);
    % axis([startdt,startdt+93,0,26]);
    % axis([startdt,startdt+5,0,26]);
    axis([startdt,datenum(2008,12,31),-100,100]);
    %%%% ??? DEBUG
    datetick3('x',2,'keeplimits');


    % Plot daily climatology comparisons
    sumfun = @nanmedian;
    cumfun = @get_jday;  n = 366; N = 24;
    x.t = grpstats(stn.(sfld).data,cumfun(stn.(sfld).date),sumfun);
    x.q0 = grpstats(stn.(climq0fld).data,cumfun(stn.(climq0fld).date),sumfun);
    x.qt = N.*cumsum(grpstats(stn.(climqtfld).data,cumfun(stn.(climqtfld).date),sumfun));
    x.bdT = N.*cumsum(grpstats(stn.(climbdTfld).data,cumfun(stn.(climbdTfld).date),sumfun));
    x.hc_dTdt = N.*cumsum(grpstats(stn.(climhcdTdt).data,cumfun(stn.(climhcdTdt).date),sumfun));


return;
