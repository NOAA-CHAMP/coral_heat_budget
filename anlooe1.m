function stn = anlooe1(stn,doHeatBudget)
%function stn = anlooe1(stn,doHeatBudget)
%
% Load in situ data for Looe Key spar (and meteorology data from nearby SMKF1
% SEAKEYS station), and appropriate climatology, reanalysis, and model data,
% then calculate a full heat budget based on LOOE1 MicroCAT sea temperature.
%
% Last Saved Time-stamp: <Tue 2011-04-19 10:18:28  Lew.Gramer>

  set_more off
  %DEBUG:
  tic,

  datapath = get_thesis_path('../data');

  if ( ~exist('stn','var') )
    stn = [];
  end;
  if ( ~exist('doHeatBudget','var') || isempty(doHeatBudget) )
    doHeatBudget = true;
  end;

  if ( ~isfield(stn,'station_name') )
    stn.station_name = 'LOOE1';
  end;
  if ( ~isfield(stn,'lon') )
    [stn.lon,stn.lat,stn.depth] = get_station_coords('looe1');
  end;

  if ( ~isfield(stn,'ndbc_air_t') )
    % Use Sombrero Key for in situ meteorological data
    smkf1 = load_all_ndbc_data([],'smkf1');
    flds = grepstruct(smkf1,'ndbc_');
    for ix=1:length(flds)
      fld = flds{ix};
      stn.(fld) = smkf1.(fld);
    end;
    smkf1 = []; clear smkf1;
  end;

  stn = verify_variable(stn,'ndbc_wind1_u');
  stn = verify_variable(stn,'ndbc_wind1_v');

  if ( ~isfield(stn,'microcat_seatemp') )
    stn = get_looe1_microcat(stn);
  end;

  if ( ~isfield(stn,'adcp_speed') )
    stn = get_looe1_adcp(stn);
  end;

  if ( doHeatBudget )

    %%%
    %% Call SCRIPT to set:
    %% Variable-name prefixes ("PFX") for various input and output datasets; AND,
    %% All station struct fieldnames used to produce heat budget 
    station_heat_budget_field_names;

    commentstr = '';

    % NOTE: Replace SMKF1 sea temperature with Looe MicroCAT data!
    % Don't forget to check sea temperature sensor depth in STATION_INSTRUMENT_HEIGHTS.m!
    stn.(sfld) = stn.microcat_seatemp;
    matfname = fullfile(datapath,[lower(stn.station_name) '_microcat_' TURPFX '_' QEPFX '_heat_budget.mat']);
    startdt = datenum(2007,2,25);

    % % NOTE: Replace SMKF1 sea temperature with Looe ADCP data!
    % % Don't forget to check sea temperature sensor depth in STATION_INSTRUMENT_HEIGHTS.m!
    % stn.(sfld) = stn.adcp_seatemp;
    % matfname = fullfile(datapath,[lower(stn.station_name) '_adcp_' TURPFX '_' QEPFX '_heat_budget.mat']);
    % startdt = datenum(2005,3,26);


    if ( exist(matfname,'file') )

      disp(['Loading from ' matfname]);
      load(matfname,'station');
      flds = fieldnames(station);
      for fldix=1:length(flds)
        fld = flds{fldix};
        stn.(fld) = station.(fld);
      end;
      station = []; clear station;


    else

      % % Get air-sea flux climatologies for comparison
      % if ( ~isfield(stn,'monthly_nocs_srf') )
      %   stn = annocs(stn);
      % end;
      % if ( ~isfield(stn,'landy_sr') )
      %   stn = station_load_landy(stn);
      % end;

      if ( ~isfield(stn,bathyfld) )
        stn = get_ngdc_bathy_station(stn);
      end;
      if ( ~isfield(stn,slopefld) )
        stn = station_ngdc_offshore_slope(stn);
      end;
      if ( ~isfield(stn,bathorifld) )
        stn = station_optimal_isobath_orientation(stn);
      end;
      if ( ~isfield(stn,hfld) || ~isfield(stn,tufld) || ~isfield(stn,tvfld) )
        switch (TIDEPFX),
         case 'tmd_tide',	stn = station_tmd_tide(stn);
         case 'tpxo_tide',	stn = station_tmd_tide(stn);
         case 'ndbc',		warning('Tidal currents not available from NDBC!');
          stn.(tufld).date = stn.(hfld).date; stn.(tufld).date = repmat(0,size(stn.(tufld).date));
          stn.(tvfld).date = stn.(hfld).date; stn.(tvfld).date = repmat(0,size(stn.(tvfld).date));
          stn.(tspdfld).date = stn.(hfld).date; stn.(tspdfld).date = repmat(0,size(stn.(tspdfld).date));
          stn.(tdirfld).date = stn.(hfld).date; stn.(tdirfld).date = repmat(0,size(stn.(tdirfld).date));
         otherwise,		error('Unknown tide dataset "%s"',TIDEPFX);
        end;
      end;
      if ( ~isfield(stn,mhfld) )
        % Calculate the mean tide depth experienced by a watermass moving over
        % an M2 tidal ellipse centered on the coordinates of our station
        stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
      end;

      if ( ~isfield(stn,Ufld) )
        [Wix,Dix] = intersect_dates(stn.(Wfld).date,stn.(Dfld).date);
        stn.(Ufld).date = stn.(Wfld).date(Wix);
        stn.(Vfld).date = stn.(Wfld).date(Wix);
        [stn.(Ufld).data,stn.(Vfld).data] = spddir_to_uv(stn.(Wfld).data(Wix),stn.(Dfld).data(Dix));
      end;
      if ( ~isfield(stn,dsrfld) || ~isfield(stn,cfld) )
        switch (RAPFX),
         case 'erai',		stn = get_erai_station(stn);
         case 'ncep',		stn = get_ncep_station(stn,'narr');
         otherwise,		error('Unavailable gridded/reanalysis dataset "%s"',RAPFX);
        end;
        stn = station_heat_flux_term(stn,raq0fld,raqtfld,sfld,[],nanmean(stn.(mhfld).data));
      end;

      % Some datasets have dew temp, some have relhumid, some have spechumid!
      if ( ~isfield(stn,dfld) && ~isfield(stn,rhfld) && ~isfield(stn,qafld) )
        error('Found neither dew-point nor (rel/spec) humidity data (%s,%s,%s)!',dfld,rhfld,qafld);
      elseif ( isfield(stn,dfld) )
        if ( ~isfield(stn,rhfld) )
          disp([dfld '->' rhfld]);
          stn = station_dewp_to_relhumid(stn,afld,dfld,rhfld);
        end;
        if ( ~isfield(stn,qafld) )
          disp([rhfld '->' qafld]);
          stn = station_relhumid_to_spechumid(stn,afld,rhfld,qafld);
        end;
      elseif ( isfield(stn,rhfld) )
        disp([rhfld '->' dfld]);
        stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
        if ( ~isfield(stn,qafld) )
          disp([rhfld '->' qafld]);
          stn = station_relhumid_to_spechumid(stn,afld,rhfld,qafld);
        end;
      elseif ( isfield(stn,qafld) )
        disp([qafld '->' rhfld]);
        stn = station_spechumid_to_relhumid(stn,afld,qafld,rhfld);
        disp([rhfld '->' dfld]);
        stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
      end;

      if ( ~isfield(stn,whfld) )
        % If waves not from reanalysis, user must want model (WaveWatch III) or wind estimate
        switch (WAVEPFX),
         case 'ww3',		stn = get_ww3_station(stn);
         case 'ndbc',		stn = station_wind_to_wave(stn,Wfld,Dfld,wpfld,whfld,wdfld);
         otherwise,		error('Unknown wave source "%s"',WAVEPFX);
        end;
      end;
      if ( ~isfield(stn,ssufld) )
        stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wfld,Dfld,whfld,wpfld,wdfld);
      end;
      if ( ~isfield(stn,ufld) )
        switch (KMPFX),
          %function stn = get_fkeys_hycom(stn_or_stnm,mindt,maxdt,vars,flds,interpMethod,fkeyspath)
         case 'fkeys_hycom',	stn = get_fkeys_hycom(stn,[],[],[],[],'linear');
          %function stn = get_gom_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
         case 'gom_hycom',	stn = get_gom_hycom(stn);
         otherwise,		error('Unknown km-scale model "%s"',KMPFX);
        end;
        more off;
      end;
      if ( ~isfield(stn.(Tfld),'gradient_x') )
        % Calculate gradients and field Laplacians
        stn = calc_field_terms(stn,Tfld);
      end;
      % Spline-fit an hourly time series of mean currents to native data
      stn.(hufld) = interp_ts(stn.(ufld));
      stn.(hvfld) = interp_ts(stn.(vfld));
      stn.(netufld) = ts_op(stn.(tufld),stn.(hufld),'+');
      stn.(netvfld) = ts_op(stn.(tvfld),stn.(hvfld),'+');

      if ( ~isfield(stn,qeufld) )
        stn = calc_quasi_eulerian(stn,STOKESPFX,KMPFX,QEPFX);
      end;


      if ( ~isfield(stn,climq0fld) )
        switch (CLIMPFX),
         case 'daily_oaflux',	stn = station_load_oaflux(stn);
         otherwise,		error('Unknown flux climatology "%s"',CLIMPFX);
        end;
        stn = station_heat_flux_term(stn,climq0fld,climqtfld,sfld,[],nanmean(stn.(mhfld).data));
      end;

      sai_opts = [];
      sai_opts.kd = [0.15,0.45];
      if ( ~isfield(stn,asrfld) )
        % If absorbed short-wave not already present, user wants absorption calculation
        stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,sai_opts);
      end;
      if ( ~isfield(stn,lrfld) )
        % If reanalysis long-wave flux was not specified, user wants a bulk estimate
        %station_bulk_longwave(stn,afld,qfld,pfld,dsrfld,sfld,cfld,dlrf,ulrf,lrf)
        %stn = station_bulk_longwave(stn,afld,qafld,pfld,dsrfld,sfld,cfld,dlrfld,ulrfld,lrfld);
        stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);
      end;


      if ( ~isfield(stn,pblzfld) )
        % Default: 600m
        pblzfld = 600;
      end;

      % Fluxes WITH warm-layer adjustment
      % Assume max warm-layer depth somewhere near bottom boundary layer top
      max_wl = nanmax(stn.(mhfld).data) - 0.5;
      stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                              pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
                              Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,true,max_wl);
      % Algorithm sometimes returns complex numbers!
      stn.(qlhfld).data = real(stn.(qlhfld).data);
      stn.(qshfld).data = real(stn.(qshfld).data);
      stn.(qrhfld).data = real(stn.(qrhfld).data);
      stn.(q0fld).data = real(stn.(q0fld).data);
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],mhfld);

      stn.(qradfld) = ts_op(stn.(asrfld),stn.(lrfld),'+');
      stn.(qturfld) = ts_op(stn.(qlhfld),stn.(qshfld),'+');

      % Net flux without absorption calculation or benthic flux - for comparison
      stn.(sqradfld) = ts_op(stn.(srfld),stn.(lrfld),'+');
      stn.(sq0fld) = ts_op(stn.(sqradfld),stn.(qturfld),'+');
      stn = station_heat_flux_term(stn,sq0fld,sqtfld,sfld,[],mhfld);

      %%%
      %% Save results to MAT file for future reference
      disp(['Saving to ' matfname]);
      station = stn;
      save(matfname,'station');
      station = []; clear station;

    end; %if ( exist(matfname,'file') ) else


    % Ignore the along-shore component of model heat advection
    stn = station_cross_shore_advection(stn,bathorifld,...
                                        qeufld,qevfld,Tfld,kmtfld,...
                                        ['raw_' udTfld],udTfld,...
                                        qtfld,qtAdvfld);
    stn = station_calc_kdel2t(stn,K_theta,Tfld,...
                              ['raw_' kd2Tfld],kd2Tfld,...
                              qtAdvfld,dTfld);
    stn = station_heat_flux_term_inverse(stn,dTffld,...
                                         dTfld,sfld,[],mhfld);


    %%%
    %% Benthic Heat Exchanges

    % % IGNORE benthic heat flux for this site
    % stn.(bq0fld) = stn.(q0fld);

    % APPLY benthic heat flux model
    sbe_opts = [];
    stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld,sbe_opts);

    stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],mhfld);

    stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'+');
    stn = station_heat_flux_term(stn,bq0fld,bq0tfld,sfld,[],mhfld);

    stn.(bdTfld) = ts_op(stn.(dTfld),stn.(qbotfld),'+');
    stn = station_heat_flux_term_inverse(stn,bdTffld,bdTfld,sfld,[],mhfld);


    %%%
    %% Horizontal Convection

    bet = stn.(slopefld);
    % stn = station_horizontal_convection(stn,sfld,[],mhfld,tspdfld,bq0fld,bet,hc_opts);

    [tix,hix,tspdix,aix,Wix,q0ix,bq0ix,dTix] = ...
      intersect_all_dates([],stn.(sfld).date,stn.(mhfld).date,stn.(tspdfld).date,stn.(afld).date,stn.(Wfld).date,stn.(q0fld).date,stn.(bq0fld).date,stn.(bdTffld).date);
    t = stn.(sfld).data(tix);
    s = repmat(36,size(stn.(sfld).data(tix)));
    h = stn.(mhfld).data(hix);
    tspd = stn.(tspdfld).data(tspdix);
    at = stn.(afld).data(aix);
    W = stn.(Wfld).data(Wix);
    q0 = stn.(q0fld).data(q0ix);
    bq0 = stn.(bq0fld).data(bq0ix);
    dT = stn.(bdTffld).data(dTix);

    qstr=bdTffld;
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

      stnfld = ['hc_' fld];
      stn.(stnfld).date = dts;
      stn.(stnfld).data = dat;
    end;

    %%%% ??? DEBUG
    % Truncate start of all time series for display purposes
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
    fac = stn.hc_termFactor.data(startix:end);
    dTdt.date = stn.hc_dTdt.date(startix:end);
    dTdt.data = stn.hc_dTdt.data(startix:end);

    bt = stn.(btfld);
    [ig,startix] = min(abs(bt.date-dts(1)));
    [ig,endix] = min(abs(bt.date-dts(end)));
    bt.date = bt.date(startix:endix);
    bt.data = bt.data(startix:endix);
    %%%% ??? DEBUG

    dsr_lpfld = [asrfld '_24_hour_sum'];
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
    bigfh=figure; maxigraph; hold on;
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
    % chkann(stn,'daily',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX);


  end; %if ( doHeatBudget )

  %DEBUG:
  toc,
  set_more;

return;
