function stn = optim_q0(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%function stn = optim_q0(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%
% Optimize all parameters of reef ocean heat budget by graphing results for the caller
%
% Last Saved Time-stamp: <Fri 2012-03-30 12:30:53  Lew.Gramer>

  doPlot = true;
  doPrint = true;

  set_more off
  timenow,
  tic,

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

% %%%% ??? HACK HACK HACK
% disp('%%%% ??? HACK HACK HACK - forcing ERAI winds');
% Wfld='erai_wind_speed';
% Dfld='erai_wind_dir';
% %%%% ??? HACK HACK HACK

% %%%% ??? HACK HACK HACK
% disp('%%%% ??? HACK HACK HACK - forcing MLRF2 in situ PAR -> insolation');
% dsrfld='bic_surf_dsrf';
% %%%% ??? HACK HACK HACK

  %%%% ??? DEBUG
  if ( ~strcmp(sfld,[ISPFX '_sea_t']) )
    disp(['Using sea temperature ' sfld]);
  end;
  if ( ~strcmp(afld,[ISPFX '_air_t']) )
    disp(['Using air temperature ' afld]);
  end;
  %%%% ??? DEBUG


  datapath = get_thesis_path('../data');
  figspath = get_thesis_path('../figs');


  stn = get_station_from_station_name(stn_or_stnm);
  stnm = lower(stn.station_name);
  STNM = upper(stn.station_name);

  opts = station_heat_budget_options(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld);

  stn.commentstr = '';

  if ( ~isfield(stn,bathyfld) )
    stn = get_ngdc_bathy_station(stn);
  end;
  if ( ~isfield(stn,slopefld) )
    stn = station_ngdc_offshore_slope(stn);
  end;
  if ( ~isfield(stn,bathorifld) )
    stn = station_optimal_isobath_orientation(stn);
  end;

  warning('off','Ecoforecasts:mergedNonTS');
  stn = load_all_ndbc_data(stn);

  %%%
  %% Sea temperature (in situ NDBC, site-specific, remote sensing, or reanalysis)

  if ( ~isfield(stn,sfld) )
    if ( regexp(sfld,'misst_sst') )
      stn = get_misst_station(stn);
      stn.hourly_misst_sst = interp_ts(stn.misst_sst);
    elseif ( regexp(sfld,'avhrr_weekly_sst') )
      stn = get_avhrr_weekly_field(stn,true);
    elseif ( regexp(sfld,'erai_') )
      if ( ~strcmpi(ISPFX,'erai') && ~strcmpi(RAPFX,'erai') )
        stn = get_erai_station(stn);
      end;
    elseif ( regexp(sfld,'microcat') )
      stn = get_looe1_microcat(stn);
    elseif ( regexp(sfld,'adcp') )
      stn = get_looe1_adcp(stn);
    elseif ( regexp(sfld,'ndbc') )
      if ( ~strcmpi(ISPFX,'ndbc') )
        stn = load_all_ndbc_data(stn);
      end;
    elseif ( regexp(sfld,'^(sea|ct_|ctd_)') )
      if ( ~strcmpi(ISPFX,'icon') )
        stn = load_station_data(stn);
      end;
    else
      error(['Unknown sea temperature field ',sfld]);
    end;
  end;


  if ( ~isfield(stn,afld) || ~isfield(stn,sfld) || ~isfield(stn,Wfld) )
    switch (ISPFX),
     case 'ndbc',	%stn = load_all_ndbc_data(stn);
     case 'erai',	stn = get_erai_station(stn);
     case 'icon',	stn = load_station_data(stn);
     otherwise,		error('Unknown in situ dataset "%s"',ISPFX);
    end;
  end;
  warning('on','Ecoforecasts:mergedNonTS');

  stn = verify_variable(stn, afld);
  stn = verify_variable(stn, pfld);
  stn = verify_variable(stn, Wfld);
  stn = verify_variable(stn, sfld);

  if ( ~isfield(stn,qsfld) )
    % Saturated ("sea-surface") specific humidity [kg/kg]
    stn = station_relhumid_to_spechumid(stn,sfld,100,qsfld);
    % Adjustment for salinity with thanks to Stommel
    stn.(qsfld).data = 0.98 .* stn.(qsfld).data;
  end;


  stn = station_tmd_tide(stn);
  stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
  % Assume max warm-layer depth somewhere near bottom boundary layer top
  max_wl = nanmax(stn.(mhfld).data) - 0.5;


  if ( ~isfield(stn,dsfld) )
    stn.(dsfld).date = stn.(sfld).date(1:end-1);
    stn.(dsfld).data = diff(stn.(sfld).data);
    stn = filter_gaps(stn,sfld,dsfld,(1.5/24));
  end;
  if ( ~isfield(stn,dsffld) )
    stn = station_heat_flux_term_inverse(stn,dsffld,dsfld,sfld,[],mhfld);
  end;


  if ( ~isfield(stn,dsrfld) || ~isfield(stn,cfld) )
    switch (RAPFX),
     case 'erai',		stn = get_erai_station(stn);
     case 'ncep',		stn = get_ncep_station(stn,'narr');
     %%% Not Yet Implemented
     % case 'era40',		stn = get_era40_station(stn);
     % case 'cfsr',		stn = get_ncep_station(stn,'cfsr');
     otherwise,		error('Unavailable gridded/reanalysis dataset "%s"',RAPFX);
    end;
    stn = station_heat_flux_term(stn,raq0fld,raqtfld,sfld,[],nanmean(stn.(mhfld).data));
  end;

  if ( ~isfield(stn,dsrfld) )
    if ( regexp(dsrfld,'bic_surf') )
      x = load_station_data('mlrf1');
      stn.bic_surf_par = x.bic_surf_par;
      x=[]; clear x
      stn = station_par_to_insol(stn,'bic_surf_par','bic_surf_dsrf','bic_surf_usrf','bic_surf_srf');
    else
      error(['Unknown insolation field ',dsrfld]);
    end;
  end;

  if ( ~isfield(stn,whfld) )
    % If waves not from reanalysis, user must want model (WaveWatch III) or wind estimate
    switch (WAVEPFX),
     case 'ww3',	stn = get_ww3_station(stn);
     case 'ndbc',	stn = station_wind_to_wave(stn,Wfld,Dfld,wpfld,whfld,wdfld);
     case 'erai',	stn = get_erai_station(stn); %IF NOT ALSO OUR REANALYSIS DATASET
     otherwise,		error('Unknown wave source "%s"',WAVEPFX);
    end;
  end;

  % If we loaded this bogus "net" field for any reason above, recalculate it properly
  if ( isfield(stn,'erai_net_heat_flux') )
    stn.erai_turbulent_heat_flux = ts_op(stn.erai_latent_heat_flux,stn.erai_sensible_heat_flux,'+');
    stn.erai_radiative_heat_flux = ts_op(stn.erai_srf,stn.erai_lrf,'+');
    stn.erai_actual_net_heat_flux = ts_op(stn.erai_turbulent_heat_flux,stn.erai_radiative_heat_flux,'+');
  end;

  %% Low-pass filter winds for quasi-Eulerian currents
  stn = verify_variable(stn,Ulpfld);
  stn = verify_variable(stn,Vlpfld);
  stn.(Wlpfld).date = stn.(Ulpfld).date;
  stn.(Wlpfld).data = uv_to_spd(stn.(Ulpfld).data,stn.(Vlpfld).data);
  stn.(Dlpfld).date = stn.(Ulpfld).date;
  stn.(Dlpfld).data = uv_to_dir(stn.(Ulpfld).data,stn.(Vlpfld).data);

  stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wlpfld,Dlpfld,whfld,wpfld,wdfld);
  % Cross- and long-shore components of currents
  stn = station_reorient_vectors(stn,bathorifld,ssufld,ssvfld);

  stn.(netufld) = ts_op(stn.(tufld),stn.(ssufld),'+');
  stn.(netvfld) = ts_op(stn.(tvfld),stn.(ssvfld),'+');
  % Cross- and long-shore components of currents
  stn = station_reorient_vectors(stn,bathorifld,netufld,netvfld);


  %% Kilometer-scale Ocean Data

  if ( ~isfield(stn,ufld) )
    switch (KMPFX),
     case 'fkeys_hycom',
      stn = get_fkeys_hycom(stn);
      opts.km_scale_advection = true;
      opts.calculate_advection = get_opt(opts,'calculate_advection',true);
      opts.calculate_diffusion = true;
     case 'gom_hycom',
      stn = get_gom_hycom(stn);
      opts.km_scale_advection = true;
      opts.calculate_advection = get_opt(opts,'calculate_advection',true);
      opts.calculate_diffusion = true;
     case 'avhrr_weekly',
      if ( ~isfield(stn,Tfld) )
        disp('Loading AVHRR_WEEKLY SST instead of hydrodynamic model data...');
        stn = get_avhrr_weekly_field(stn,true);
      end;
      opts.km_scale_advection = false;
      opts.calculate_advection = get_opt(opts,'calculate_advection',true);
      opts.calculate_diffusion = true;
      disp('ONLY STOKES');
     case 'none',
      disp('Loading NO hydrodynamic model data...');
      opts.km_scale_advection = false;
      opts.calculate_advection = false;
      opts.calculate_diffusion = false;
      disp('ONLY STOKES');
     otherwise,
      error('Unknown km-scale data source "%s"',KMPFX);
    end;
    more off;
  end;


  opts.grid_interp_method = get_opt(opts,'grid_interp_method','linear');


  if ( ~isfield(stn,hkmtxfld) )
    warning('No km-scale sea temperature gradients');
    stn.(udTfld) = stn.(sfld);
    stn.(udTfld).data(:) = 0;
    stn.(udTffld) = stn.(udTfld);
    stn.(udTffld).data(:) = 0;

  else

    if ( ~isfield(stn,hkmtxsfld) )
      % Cross- and long-shore components of sea temperature gradient
      stn = station_reorient_vectors(stn,bathorifld,hkmtxfld,hkmtyfld);
    end;

    % stn.(udTxsfld) = ts_op(stn.(hkmtxsfld),stn.(netxsfld),'.*');
    stn.(udTxsfld) = ts_op(stn.(hkmtxsfld),stn.(ssxsfld),'.*');
    % Convert to units of [K/hr]
    stn.(udTxsfld).data = -(3600*stn.(udTxsfld).data);
    stn = station_heat_flux_term_inverse(stn,udTfxsfld,udTxsfld,sfld,[],mhfld);

    % stn.(udTlsfld) = ts_op(stn.(hkmtlsfld),stn.(netlsfld),'.*');
    stn.(udTlsfld) = ts_op(stn.(hkmtlsfld),stn.(sslsfld),'.*');
    % Convert to units of [K/hr]
    stn.(udTlsfld).data = -(3600*stn.(udTlsfld).data);
    stn = station_heat_flux_term_inverse(stn,udTflsfld,udTlsfld,sfld,[],mhfld);

    % Include both vector components of model heat advection
    stn.(udTfld) = ts_op(stn.(udTxsfld),stn.(udTlsfld),'+');
    stn = station_heat_flux_term_inverse(stn,udTffld,udTfld,sfld,[],mhfld);

  end;


  %% Surface fluxes

  % This result is NOT used - except to filter dates...
  stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,opts);

  stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);

  %%%% ??? DEBUG
  disp('Forcing reanalysis downward longwave');
  stn.(dlrfld) = stn.erai_dlrf;
  %stn.(ulrfld) = stn.erai_ulrf;
  stn.(lrfld) = ts_op(stn.erai_dlrf,stn.erai_ulrf,'-');

  stn.commentstr = [stn.commentstr,' (ERAI DLRF) '];
  %%%% ??? DEBUG


  doWarms = [false, true];
  kds = {...
      0.1,[0.05,0.15,45],[0.05,0.15,137],[0.05,0.15,228],[0.05,0.15,320],...
      0.2,[0.10,0.30,45],[0.10,0.30,137],[0.10,0.30,228],[0.10,0.30,320],...
      0.3,[0.20,0.40,45],[0.20,0.40,137],[0.20,0.40,228],[0.20,0.40,320],...
        };
  %cbds = {0,3.8e-5,3.8e-4,8.0e-4,16.0e-4,24.0e-4};
  cbds = {8.0e-4};
  advfacs = {...
      0,...
            };
  kths = {...
      0 ,...
      5, [0 ,10, 45],...
      10, [0 ,20, 45],...
  % LONF1 - related to tidal currents?
  %[0 ,5 ,183],...
         };


  switch ( STNM ),

   case 'FWYF1',
    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.035,0.220, 76], };
    % kths = { 0 };

    % % ERAI met, ERAI air_t, new Ppen - warm layer adjustment matches
    % % climatology closely, but grossly overestimates interannual budget
    % doWarms = [true];
    % kds = { ...
    %     [0.035,0.220, 91], ...
    %       };
    % kths = { ...
    %     [0,20, 45] ...
    %        };

    % % ERAI met, NDBC air_t, new Ppen - warm layer adjustment matches
    % % climatology closely, but grossly overestimates interannual budget
    % % doWarms = [true];
    % doWarms = [false];
    % kds = { ...
    %     [0.025,0.220, 91], ...
    %       };
    % kths = { ...
    %     [0,10, 45] ...
    %        };

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves
    doWarms = [true];
    kds = { ...
        [0.050,0.350,137] ...
        [0.050,0.350,183] ...
        [0.050,0.350,228] ...
        0.2 ...
        [0.100,0.350,137] ...
        [0.100,0.350,183] ...
        [0.100,0.350,228] ...
        0.225 ...
          };
    kths = { ...
        0 ...
        [0,10, 45] ...
        [0,20, 45] ...
           };

   case 'MLRF1',
    % % ERAI met, air_t, old Ppen
    % kds = { [0.045,0.375,45], };
    % kths = { 0 };

    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.035,0.190, 55], };
    % kths = { 0 };

    % % ERAI met, NDBC air_t, new Ppen - Kd optimized to close interannual budget
    % doWarms = [true];
    % kds = { ...
    %     [0.050,0.250,137], ...
    %       };
    % kths = { ...
    %     0, ...
    %     % [-10, -5,137] ...
    %     % [-12, -4,137] ...
    %     % [-15,  0,137] ...
    %        };

    % % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves
    % doWarms = [true];
    % kds = { ...
    %     [0.050,0.350,137] ...
    %     [0.050,0.350,183] ...
    %     [0.050,0.350,228] ...
    %     0.2 ...
    %     [0.100,0.350,137] ...
    %     [0.100,0.350,183] ...
    %     [0.100,0.350,228] ...
    %     0.225 ...
    %       };
    % kths = { ...
    %     0 ...
    %     [0,10, 45] ...
    %     [0,20, 45] ...
    %        };

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
    doWarms = [true];
    kds = { ...
        [0.050,0.350,137] ...
        [0.100,0.300, 45] ...
        [0.100,0.300,137] ...
        [0.100,0.300,183] ...
          };
    advfacs = { ...
        0,...
        1,...
              };
    kths = { ...
        0 ...
        [0,10, 45] ...
        [0,20, 45] ...
           };

   case 'LONF1',
    % ERAI met, NDBC air_t, old Ppen
    kds = { [0.20,0.37,228], };
    kths = { [0 ,5 ,183], };

   case 'SMKF1',
    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.045,0.375,45], };
    % kths = { 0 };

    % % ERAI met, NDBC air_t, higher Ppen - amplitude problem but good annual...
    % kds = { [0.045,0.375, 91], }

    % % ERAI met, NDBC air_t, new and improved higher Ppen
    % kds = { [0.060,0.400, 45], };
    % kths = { 0 };

    % ERAI met, NDBC air_t, new Ppen - Kd optimized to close interannual budget
    % Q0 w/Warm Layer also matches OAFlux annual amplitude! But 2005 TOO HOT
    doWarms = [true];
    % Satellite chlor_a from USF suggests a SEMI-annual Kd cycle
    % kds = { [0.060,0.250, 45, 183], };
    % kds = { [0.100,0.250, 45, 183], };
    kds = { ...
        % [0.060,0.250, 45, 183],...
        % [0.100,0.250, 45, 183],...
        % [0.060,0.250, 45],...
        [0.100,0.250, 45],...
        [0.100,0.300, 45],...
          };
    kths = { 0 };

   case 'LOOE1',

   case 'HAWK1',

   case 'SANF1',
    % ERAI met, NDBC air_t, old Ppen
    kds = { [0.20,0.37,228], };
    kths = { 0 };

   case 'DRYF1',
    opts.kd = get_opt(opts,'kd',[0.045,0.375,45]);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',8.0e-4);
    opts.K_theta = get_opt(opts,'K_theta',0);

    % ERAI met, NDBC air_t, old Ppen
    kds = { [0.20,0.40,320], };
    kths = { 0 };

   otherwise,
    error('Station %s options not implemented yet!',STNM);
  end;


  % default_cbd = 3.8e-4; %From literature
  default_cbd = 8.0e-4; %Best for all coral reef sites tried
  [ig,default_cbdix] = min(abs([cbds{:}] - default_cbd));

  for ix=1:numel(kds)
    stn.optim_q0.kdstrs{ix} = num2str(kds{ix},'%g,');
  end;
  for ix=1:numel(cbds)
    stn.optim_q0.cbdstrs{ix} = num2str(cbds{ix},'C_b_d=%g');
  end;
  for ix=1:numel(advfacs)
    stn.optim_q0.advfacstrs{ix} = num2str(advfacs{ix},'%g,');
  end;
  for ix=1:numel(kths)
    stn.optim_q0.kthstrs{ix} = num2str(kths{ix},'%g,');
  end;

  for doWarmix=1:numel(doWarms)
    doWarm = logical(doWarms(doWarmix));

    % tic,

    % Fluxes WITH or WITHOUT warm-layer adjustment
    % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
    %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
    %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);
    stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                            pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
                            Dfld,ssufld,ssvfld,wpfld,whfld,pblzfld,doWarm,max_wl);
                            % Dfld,ssufld,ssvfld,wpfld,whfld,600,doWarm,max_wl);
                            % % Dfld,netufld,netvfld,wpfld,whfld,pblzfld,doWarm,max_wl);

    % Algorithm sometimes returns complex numbers!
    stn.(qlhfld).data = real(stn.(qlhfld).data);
    stn.(qshfld).data = real(stn.(qshfld).data);
    stn.(qrhfld).data = real(stn.(qrhfld).data);
    stn.(qturfld) = ts_op(stn.(qlhfld),stn.(qshfld),'+');
    if ( isfield(stn,qrhfld) && is_valid_ts(stn.(qrhfld)) )
      stn.(qturfld) = ts_op(stn.(qturfld),stn.(qrhfld),'+');
    end;
    badix = find(~isfinite(stn.(qturfld).data));
    stn.(qturfld).date(badix) = [];
    stn.(qturfld).data(badix) = [];
    stn = station_heat_flux_term(stn,qturfld,qturtfld,sfld,[],mhfld);

    % Cross- and long-shore components of the wind stress
    stn = station_reorient_vectors(stn,bathorifld,tauxfld,tauyfld);

    % toc,

    % Just in case something above reset it!
    more off;

    disp(['Looping ' num2str(numel(kds)) ' attenuation options']);
    for kdix=1:numel(kds)
      kd = kds{kdix};
      opts.kd = kd;
      opts.kd_debug = false;

      stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,opts);

      stn.(qradfld) = ts_op(stn.(asrfld),stn.(lrfld),'+');
      badix = find(~isfinite(stn.(qradfld).data));
      stn.(qradfld).date(badix) = [];
      stn.(qradfld).data(badix) = [];
      stn = station_heat_flux_term(stn,qradfld,qradtfld,sfld,[],mhfld);
      stn.(qcoolfld) = ts_op(stn.(lrfld),stn.(qturfld),'+');
      stn = station_heat_flux_term(stn,qcoolfld,qcooltfld,sfld,[],mhfld);

      stn.(q0fld) = ts_op(stn.(qradfld),stn.(qturfld),'+');
      % badix = find(~isfinite(stn.(q0fld).data));
      % stn.(q0fld).date(badix) = [];
      % stn.(q0fld).data(badix) = [];
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],mhfld);

      % Net flux without absorption calculation or benthic flux - for comparison
      stn.(sqradfld) = ts_op(stn.(srfld),stn.(lrfld),'+');
      stn.(sq0fld) = ts_op(stn.(sqradfld),stn.(qturfld),'+');
      stn = station_heat_flux_term(stn,sq0fld,sqtfld,sfld,[],mhfld);

      for cbdix=1:numel(cbds)
        cbd = cbds{cbdix};
        opts.benthic_debug = false;
        opts.b_convective_coefficient = cbd;

        %% Benthic Heat Exchanges
        stn = station_benthic_exchange(stn,sfld,tufld,tvfld,qbfld,btfld,qbofld,opts);
        % badix = find(~isfinite(stn.(qbofld).data) | abs(stn.(qbofld).data)>2e3);
        badix = find(~isfinite(stn.(qbofld).data));
        stn.(qbofld).date(badix) = [];
        stn.(qbofld).data(badix) = [];
        stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],mhfld);

        stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'-');
        stn = station_heat_flux_term(stn,bq0fld,bq0tfld,sfld,[],mhfld);
        if ( isfield(stn,bq0lpfld) ); stn = rmfield(stn,bq0lpfld); end;
        stn = verify_variable(stn,bq0lpfld);

        disp(['(Re-)Looping ' num2str(numel(advfacs)) ' advection options']);
        for advix=1:numel(advfacs)
          advfac = advfacs{advix};
          opts.advection_factor = advfac;

          %% Km-scale Heat Advection
          if ( opts.calculate_advection )
            adv.date = stn.(udTfld).date;
            adv.data = advfac.*stn.(udTfld).data;
            stn.(qtAdvfld) = ts_op(stn.(bq0tfld),adv,'+');
            adv=[]; clear adv
          else
            stn.(qtAdvfld) = stn.(bq0tfld);
          end; %if ( opts.calculate_advection )
          stn = station_heat_flux_term_inverse(stn,qtAdvffld,qtAdvfld,sfld,[],mhfld);



          disp(['(Re-)(Re-)Looping ' num2str(numel(kths)) ' diffusion options']);
          for kthix=1:numel(kths)
            kth = kths{kthix};
            opts.K_theta = kth;


            %% Km-scale Heat Diffusion
            if ( opts.calculate_diffusion )
              stn = station_calc_kdel2t(stn,opts.K_theta,Tfld,...
                                        rawkd2Tfld,kd2Tfld,...
                                        qtAdvfld,bdTfld,opts.grid_interp_method,false);
              stn = station_heat_flux_term_inverse(stn,kd2Tffld,kd2Tfld,sfld,[],mhfld);
            else
              stn.(bdTfld) = stn.(qtAdvfld);
            end; %if ( isempty(opts.laplacian_climatology) ) else
            stn = station_heat_flux_term_inverse(stn,bdTffld,bdTfld,sfld,[],mhfld);


            %% Low-pass filtered total heat flux
            if ( isfield(stn,bdTflpfld) ); stn = rmfield(stn,bdTflpfld); end;
            stn = verify_variable(stn,bdTflpfld);


            %% Apply Horizontal Convection to total fluxes
            bet = stn.(slopefld);
            [tix,hix,tspdix,aix,Wix,sq0ix,bq0ix,dTix] = ...
                intersect_all_dates([],stn.(sfld).date,stn.(mhfld).date,stn.(tspdfld).date,stn.(afld).date,stn.(Wfld).date,stn.(sq0fld).date,stn.(bq0fld).date,stn.(bdTflpfld).date);
            dts = stn.(sfld).date(tix);
            t = stn.(sfld).data(tix);
            s = repmat(36,size(stn.(sfld).data(tix)));
            h = stn.(mhfld).data(hix);
            tspd = stn.(tspdfld).data(tspdix);
            at = stn.(afld).data(aix);
            W = stn.(Wfld).data(Wix);
            sq0 = stn.(sq0fld).data(sq0ix);
            bq0 = stn.(bq0fld).data(bq0ix);
            dT = stn.(bdTflpfld).data(dTix);

            % opts.hc_debug = get_opt(opts,'hc_debug',false);
            opts.hc_debug = false;
            % opts.hc_R = get_opt(opts,'hc_R',(1.00-0.08));
            opts.hc_scaling = get_opt(opts,'hc_scaling','US');
            opts.hc_max_onset_secs = get_opt(opts,'hc_max_onset_secs',12*3600);
            res = horizontal_convection(t,s,h,dT,bet,opts,dts,dT,W);
            stn.(hcdTdt).date = dts;
            stn.(hcdTdt).data = res.dTdt;
            stn = station_heat_flux_term_inverse(stn,hcdTdtf,hcdTdt,sfld,[],mhfld);

            % Horizontal convection diagnostics
            stn.(hcu).date = dts;
            stn.(hcu).data = res.u;
            stn.(hcdTdx).date = dts;
            stn.(hcdTdx).data = res.dTdx;
            stn = station_heat_flux_term_inverse(stn,hcdTdxf,hcdTdx,sfld,[],mhfld);
            stn.(hcdTdthc).date = dts;
            stn.(hcdTdthc).data = res.dTdthc;
            stn = station_heat_flux_term_inverse(stn,hcdTdthcf,hcdTdthc,sfld,[],mhfld);


%%%%??? DEBUG - TRY ELIMINATING HORIZONTAL CONVECTION!
% stn.(hcdTdt) = stn.(bq0tfld);
% stn.(hcdTdt) = stn.(bdTflpfld);
%%%%??? DEBUG

            % annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,[],[1:3]);

            if ( doPlot )

              %% Calculate and plot time series and errors from this option
              clear t q
              [tix,qix] = intersect_dates(stn.(sfld).date,stn.(hcdTdt).date);
              t.date = stn.(sfld).date(tix);	t.data = stn.(sfld).data(tix);
              q.date = stn.(hcdTdt).date(qix);	q.data = stn.(hcdTdt).data(qix);

%%%%??? DEBUG
begyr = 1996;
if ( strcmpi(KMPFX,'none') ); begyr = 1970; end;
endyr = 2011;
begdt = datenum(begyr,1,2);
enddt = datenum(endyr,1,1);
t.data(begdt>t.date|t.date>enddt)=[]; t.date(begdt>t.date|t.date>enddt)=[];
q.data(begdt>q.date|q.date>enddt)=[]; q.date(begdt>q.date|q.date>enddt)=[];
begyr = min(get_year(t.date));
endyr = max(get_year(t.date));
%%%%??? DEBUG


              % Evaluate total error for this option
              t0 = t.data(1);
              sq.date = q.date;
              sq.data = t0 + nancumsum(q.data) - q.data(1);
              stn.optim_q0.error(kdix,doWarmix,advix,kthix) = sqrt(sum((t.data-sq.data).^2));

              stn.optim_q0.doWarm(kdix,doWarmix,advix,kthix) = doWarm;
              stn.optim_q0.kd{kdix,doWarmix,advix,kthix} = kd;
              stn.optim_q0.cbd(kdix,doWarmix,advix,kthix) = cbd;
              stn.optim_q0.kth{kdix,doWarmix,advix,kthix} = kth;
              stn.optim_q0.q(kdix,doWarmix,advix,kthix) = q;
              stn.optim_q0.sq(kdix,doWarmix,advix,kthix) = sq;

              % Evaluate climatological error for this option
              if ( ~isfield(stn.optim_q0,'climt') )
                [cum,tid] = grp_ts(t.data,t.date,'daily');
                stn.optim_q0.climt.data = cum;
                stn.optim_q0.climt.date = tid;
              end;
              % [cum,tid] = grp_ts(q.data,q.date,'daily',@nansum);
              [cum,tid] = grp_ts(q.data,q.date,'daily',@nanmean);
              cum = 24*cum;
              stn.optim_q0.climq(kdix,doWarmix,advix,kthix).data = cum;
              stn.optim_q0.climq(kdix,doWarmix,advix,kthix).date = tid;
              t0 = stn.optim_q0.climt.data(1);
              sq.date = tid;
              sq.data = t0 + nancumsum(cum) - cum(1);
              stn.optim_q0.climsq(kdix,doWarmix,advix,kthix).date = sq.date;
              stn.optim_q0.climsq(kdix,doWarmix,advix,kthix).data = sq.data;
              stn.optim_q0.climerror(kdix,doWarmix,advix,kthix) = sqrt(sum((stn.optim_q0.climt.data-sq.data).^2));

              % Evaluate seasonal amplitude error for this option
              stn.optim_q0.minseassq = +Inf;
              stn.optim_q0.maxseassq = -Inf;
              yrs = unique(get_year(t.date));
              for yrix = 1:numel(yrs)
                yr = yrs(yrix);
                ix = find(get_year(t.date)==yr);
                if ( isempty(ix) )
                  warning('NO DATA for year %g??',yr);
                else
                  t0 = t.data(ix(1));
                  sq.date = q.date(ix);
                  sq.data = t0 + nancumsum(q.data(ix)) - q.data(ix(1));
                  stn.optim_q0.seasyear(kdix,doWarmix,advix,kthix,yrix) = yr;
                  stn.optim_q0.seassq(kdix,doWarmix,advix,kthix,yrix) = sq;
                  stn.optim_q0.minseassq = nanmin(stn.optim_q0.minseassq,nanmin(sq.data(:)));
                  stn.optim_q0.maxseassq = nanmax(stn.optim_q0.maxseassq,nanmax(sq.data(:)));
                  stn.optim_q0.seaserror(kdix,doWarmix,advix,kthix,yrix) = sqrt(sum((t.data(ix)-sq.data).^2));
                end; %if isempty(ix) else
              end; %for yrix

            end; %if doPlot

          end; %for kthix=1:numel(kths)

        end; %for advix=1:numel(advfacs)

      end; %for cbdix=1:numel(cbds)

    end; %for kdix=1:numel(kds)

  end; %for doWarmix=1:numel(doWarms)

  if ( doPlot )

    %{
    for doWarmix=1:numel(doWarms)
      doWarm = logical(doWarms(doWarmix));
      if ( doWarm );	doWarmStr = 'Warm Layer';
      else;		doWarmStr = 'No Warm Layer';
      end;

      for advix=1:numel(advfacs)
        advfac = advfacs{advix};
        advStr = [' Adv=',num2str(advfac,'%g,')];

        fmg; plot(squeeze(stn.optim_q0.error(:,doWarmix,advix,:))); titlename([STNM ' Total Error: ' doWarmStr advStr]);
        legend(stn.optim_q0.kthstrs);
        set(gca,'XTick',[1:numel(kds)],'XTickLabel',stn.optim_q0.kdstrs);
        ylim([nanmin(stn.optim_q0.error(:)),nanmax(stn.optim_q0.error(:))]);

        fmg; plot(squeeze(nanmedian(stn.optim_q0.seaserror(:,doWarmix,advix,:),4))); titlename([STNM ' Mdn Seas Err: ' doWarmStr advStr]);
        legend(stn.optim_q0.kthstrs);
        set(gca,'XTick',[1:numel(kds)],'XTickLabel',stn.optim_q0.kdstrs);
        ylim([nanmin(stn.optim_q0.seaserror(:)),nanmax(stn.optim_q0.seaserror(:))]);
      end;
    end;
    %}


    for advix=1:numel(advfacs)
      advfac = advfacs{advix};
      for kthix=1:numel(kths)
        kth = kths{kthix};

        for doWarmix=1:numel(doWarms)
          doWarm = logical(doWarms(doWarmix));
          if ( doWarm ); doWarmStr = 'Warm Layer';
          else;	         doWarmStr = 'No Warm Layer';
          end;

          fmg; plot_ts(stn.optim_q0.climt,squeeze(stn.optim_q0.climsq(:,doWarmix,advix,kthix)));
          datetick3('x',3);
          titlename([STNM ' Daily Clim: ' doWarmStr ' ' stn.optim_q0.cbdstrs{default_cbdix}]);
          legend({'T_s',stn.optim_q0.kdstrs{:}});
          xlim(stn.optim_q0.climt.date([1 end]));
          ylim([minmin([stn.optim_q0.climsq.data]),maxmax([stn.optim_q0.climsq.data])]);
          % ylim([21,35]);
          % ylim([15,50]);
          % ylim([15,35]);
          ylim([16,34]);
          ylabel('^oC');
          % appendtitlename([' (' strrep(KMPFX,'_','\_') ' K_\theta=' num2str(kth,'%g,') ')']);
          appendtitlename([' (' strrep(QEPFX,'_','\_') ' Adv=' num2str(advfac,'%g,') ' K_\theta=' num2str(kth,'%g,') ')']);
          appendtitlename([' (' num2str(begyr) '-' num2str(endyr) ')']);
        end; %for doWarmix=1:numel(doWarms)
        1;

        %{
        for doWarmix=1:numel(doWarms)
          doWarm = logical(doWarms(doWarmix));
          if ( doWarm ); doWarmStr = 'Warm Layer';
          else;	         doWarmStr = 'No Warm Layer';
          end;

          fmg; plot_ts(stn.(sfld),squeeze(stn.optim_q0.sq(:,doWarmix,advix,kthix))); titlename([STNM ' Time Series: ' doWarmStr]);
          legend({'T_s',stn.optim_q0.kdstrs{:}});
          ylim([minmin([stn.optim_q0.sq.data]),maxmax([stn.optim_q0.sq.data])]);
        end;

        yrs = unique(get_year(stn.optim_q0.sq(1,1,1,kthix).date));
        nrows = floor(sqrt(numel(yrs)));
        ncols = ceil(numel(yrs)/nrows);

        for doWarmix=1:numel(doWarms)
          doWarm = logical(doWarms(doWarmix));
          fmg;
          if (~doWarm)	suptitle([STNM ' Annual TS: No Warm Layer']);
          else		suptitle([STNM ' Annual TS: Warm Layer']);
          end;
          for yrix=1:numel(yrs)
            yr = yrs(yrix);
            if ( is_valid_ts(stn.optim_q0.seassq(1,doWarmix,advix,kthix,yrix)) )
              ix = find(get_year(stn.(sfld).date)==yr);
              t.date = stn.(sfld).date(ix);
              t.data = stn.(sfld).data(ix);

              subplot_tight(nrows,ncols,yrix);
              plot_ts(t,stn.optim_q0.seassq(:,doWarmix,advix,kthix,yrix));
              datetick('x',17,'keeplimits');
              % legend({'T_s',stn.optim_q0.kdstrs{:}});
              xlabel(num2str(yr));
              xlim([min(stn.optim_q0.seassq(kdix,doWarmix,advix,kthix,yrix).date),...
                    max(stn.optim_q0.seassq(kdix,doWarmix,advix,kthix,yrix).date)]);
              ylim([stn.optim_q0.minseassq,stn.optim_q0.maxseassq]);
              % ylim([minmin([stn.optim_q0.seassq.data]),maxmax([stn.optim_q0.seassq.data])]);
            end; %if is_valid_ts
          end; %for yrix
        end; %for doWarmix=1:numel(doWarms)
        1;
        %}

      end; %for kthix
    end; %for advix

    if ( doWarm );	doWarmStr = ' (WARM)';
    else;		doWarmStr = ' (NO WARM)';
    end;


    % Compare published climatologies to our (last) estimate
    stn = station_load_oaflux(stn);


    %%%%??? DEBUG
    % % Plot daily climatology comparisons
    % compare_flux_climatologies(stn,'daily',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,doWarmStr);
    % figfname= fullfile(figspath,[lower(stn.station_name) '_' hcdTdt '_chkann.']);
    % % print('-dtiff',[figfname 'tiff']);
    % % print('-dpng',[figfname 'png']);
    %%%%??? DEBUG

    %%%%??? DEBUG
    % % Show year-by-year comparison for whichever configuration we calculated last!
    % annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,doWarmStr,[],begyr);
    %%%%??? DEBUG

    %%%%??? DEBUG
    % ms1_clim(stn,'daily',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,[],[],[],[],1996);
    % figfname= fullfile(figspath,[lower(stn.station_name) '_' hcdTdt '_weekly.']);
    % % % print('-dtiff',[figfname 'tiff']);
    % % % print('-dpng',[figfname 'png']);
    %%%%??? DEBUG


    %%%%??? DEBUG
    % dys = datenum(lastyr,12,31) - datenum(firstyr,1,1) + 1;
    % plot_fluxes(stn,begyr,1,dys,{sfld,btfld},[],...
    %             {'ndbc_hfbulk_heat_flux_term',qtfld,dTfld,bdTfld,[bdTfld '_netqf']},[],...
    %             {'NDBC sea temperature','Modeled substrate temperature',...
    %              'Bulk Q_0',...
    %            '(Q_0 == \gammaQ_S_W + Q_L_W + Q_L_H + Q_S_H)/\rhoC_ph',...
    %            'K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + Q_0/\rhoC_ph',...
    %            'K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + (Q_0+Q_b)/\rhoC_ph',...
    %            'HC(K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + (Q_0+Q_b)/\rhoC_ph)',...
    %           });
    % appendtitlename(strrep(sprintf(' (%s,%s)', upper(RAPFX), upper(KMPFX)),'_','\_'));
    %%%%??? DEBUG


    % Compare monthly climatologies vs. flux implied by actual dTs vs. our (last) estimate
    stn = compare_monthly_flux_climatologies(stn,sfld,sq0fld);

    % Compare simple cumulative sums of daily climatology vs. our (last) estimate
    fmg;
    [dix,cix,six,bix,hix] = ...
        intersect_all_dates([],stn.(dsffld).date,stn.(climq0fld).date,stn.(sq0fld).date,stn.(bq0fld).date,stn.(hcdTdtf).date);
    plot(stn.(dsffld).date(dix(1):end),nancumsum(stn.(dsffld).data(dix(1):end)),'k');
    plot(stn.(climq0fld).date(cix(1):end),nancumsum(stn.(climq0fld).data(cix(1):end).*24),'c');
    plot(stn.(sq0fld).date(six(1):end),nancumsum(stn.(sq0fld).data(six(1):end)),'r-.');
    plot(stn.(bq0fld).date(bix(1):end),nancumsum(stn.(bq0fld).data(bix(1):end)),'m:');
    plot(stn.(hcdTdtf).date(hix(1):end),nancumsum(stn.(hcdTdtf).data(hix(1):end)),'b--');
    legend('Actual \rhoC_ph\partial_tT_s','OAFlux/ISCCP Q_0',...
           'G&M Q_0','G&M Q_0(\gamma)+Q_b','G&M \rhoC_ph\partial_tT_s',...
           'Location','Best');
    datetick3('x',20,'keeplimits');
    ylim([-3e6,+3e6]);    ylabel('W/m^2');
    titlename([STNM ' simple cumulative sums: Hourly fluxes K_d=' num2str(opts.kd) ' K\theta=' num2str(opts.K_theta) stn.commentstr doWarmStr ]);


    % Create one-day averages (which we will sub-sample) to compare with OAFlux climatology
    s1dfld = [sfld '_1_d_avg'];
    stn = verify_variable(stn,s1dfld);
    a1dfld = [afld '_1_d_avg'];
    stn = verify_variable(stn,a1dfld);
    sq01dfld = [sq0fld '_1_d_avg'];
    stn = verify_variable(stn,sq01dfld);
    sr1dfld = [srfld '_1_d_avg'];
    stn = verify_variable(stn,sr1dfld);
    asr1dfld = [asrfld '_1_d_avg'];
    stn = verify_variable(stn,asr1dfld);
    lr1dfld = [lrfld '_1_d_avg'];
    stn = verify_variable(stn,lr1dfld);
    qlh1dfld = [qlhfld '_1_d_avg'];
    stn = verify_variable(stn,qlh1dfld);
    qsh1dfld = [qshfld '_1_d_avg'];
    stn = verify_variable(stn,qsh1dfld);
    qcool1dfld = [qcoolfld '_1_d_avg'];
    stn = verify_variable(stn,qcool1dfld);
    udTffld1d = [udTffld '_1_d_avg'];
    stn = verify_variable(stn,udTffld1d);
    kd2Tffld1d = [kd2Tffld '_1_d_avg'];
    stn = verify_variable(stn,kd2Tffld1d);
    hcdTdthcf1d = [hcdTdthcf '_1_d_avg'];
    stn = verify_variable(stn,hcdTdthcf1d);


    [cix,six] = intersect_dates(stn.(climq0fld).date,stn.(sq01dfld).date);

if(0)
    fmg;
    boxplot_ts(stn.(climq0fld),'month','index',cix,...
               'title',[STNM,' OAFlux/ISCCP Q_0']);
    ylim([-1000,1000]); ylabel('W/m^2');

    fmg;
    boxplot_ts(stn.(sq01dfld),'month','index',six,...
               'title',[STNM,' Gramer&Mariano Sea-surface Q_0',doWarmStr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',sq01dfld,'.tiff']));
    end;


    [cix,six] = intersect_dates(stn.(climsrfld).date,stn.(sr1dfld).date);

    fmg;
    boxplot_ts(stn.(climsrfld),'month','index',cix,...
               'title',[STNM,' OAFlux/ISCCP Q_S_W']);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',climsrfld,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(sr1dfld),'month','index',six,...
               'title',[STNM,' Gramer&Mariano Q_S_W']);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',sr1dfld,'.tiff']));
    end;


    [cix,six] = intersect_dates(stn.(climqlhfld).date,stn.(qlh1dfld).date);

    fmg;
    boxplot_ts(stn.(climqlhfld),'month','index',cix,...
               'title',[STNM,' OAFlux Q_L_H']);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',climqlhfld,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(qlh1dfld),'month','index',six,...
               'title',[STNM,' Gramer&Mariano Q_L_H',doWarmStr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',qlh1dfld,'.tiff']));
    end;


    fmg;
    sh=boxplot_ts(stn.(sr1dfld),'month','allcolors','b');
    lh=boxplot_ts(stn.(qcool1dfld),'month','allcolors','r');
    ylim([-1000,1000]); ylabel('W/m^2');
    legend([sh(1),lh(1)], 'Q_S_W','Q_L_W+Q_L_H+Q_S_H', 'Location','South');
    titlename([STNM,' Gramer&Mariano Air-Sea fluxes (1d avg)',doWarmStr]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',sr1dfld,'-vs-',qcool1dfld,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(udTffld1d),'month','allcolors','k',...
               'title',[STNM,' Gramer&Mariano Km-scale Heat advection (1d avg)',doWarmStr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',udTffld1d,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(kd2Tffld1d),'month','allcolors','k',...
               'title',[STNM,' Gramer&Mariano Km-scale Heat diffusion (1d avg)',doWarmStr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',kd2Tffld1d,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(hcdTdthcf1d),'month','allcolors','k',...
               'title',[STNM,' Gramer&Mariano Horizontal Convection (1d avg)',doWarmStr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',hcdTdthcf1d,'.tiff']));
    end;
end;

    fmg;
    lh=boxplot_ts(stn.(sr1dfld),'month','allcol','b');
    rh=boxplot_ts(stn.(climsrfld),'month','allcol','r');
    legend([lh(1),rh(1)],'ERA-Interim','OAFlux/ISCCP');
    titlename([STNM,' Net insolation Q_S_W']);
    ylabel('Wm^-^2'); ylim([-1000,1000]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',sr1dfld,'-vs-',climsrfld,'.tiff']));
    end;

    fmg;
    lh=boxplot_ts(stn.(lr1dfld),'month','allcol','b');
    rh=boxplot_ts(stn.(climlrfld),'month','allcol','r');
    legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux/ISCCP');
    titlename([STNM,' Net longwave flux Q_L_W']);
    ylabel('Wm^-^2'); ylim([-1000,1000]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',lr1dfld,'-vs-',climlrfld,'.tiff']));
    end;

    fmg;
    lh=boxplot_ts(stn.(qlh1dfld),'month','allcol','b');
    rh=boxplot_ts(stn.(climqlhfld),'month','allcol','r');
    legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux');
    titlename([STNM,' Latent heat flux Q_L_H',doWarmStr]);
    ylabel('Wm^-^2'); ylim([-1000,1000]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',qlh1dfld,'-vs-',climqlhfld,'.tiff']));
    end;

    fmg;
    lh=boxplot_ts(stn.(qsh1dfld),'month','allcol','b');
    rh=boxplot_ts(stn.(climqshfld),'month','allcol','r');
    legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux');
    titlename([STNM,' Sensible heat flux Q_S_H',doWarmStr]);
    ylabel('Wm^-^2'); ylim([-1000,1000]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',qsh1dfld,'-vs-',climqshfld,'.tiff']));
    end;

    nd = ts_op(stn.(a1dfld),stn.(s1dfld),'-',@(x)(intersect_dates(x.date,stn.(climafld).date)));
    od = ts_op(stn.(climafld),stn.(climsfld),'-',@(x)(intersect_dates(x.date,stn.(a1dfld).date)));
    fmg; grpplot_ts(nd,@get_week,@nanmedian,0,'b'); grpplot_ts(od,@get_week,@nanmedian,0,'r');
    legend('In situ','OAFlux');
    titlename([STNM,' Air-Sea Temperature Difference']);
    ylabel('K'); ylim([-5,2]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-median-air-sea-temperature-',ISPFX,'-vs-',CLIMPFX,'.tiff']));
    end;

  end; %if ( doPlot )

  matfname = fullfile(datapath,[stnm '_' taufld '.mat']);
  if ( ~exist(matfname,'file') )
    date = stn.(taufld).date;
    data = stn.(taufld).data;
    save(matfname,'date','data');
    clear date data;
  end;

  toc,
  timenow,
  set_more;

return;
