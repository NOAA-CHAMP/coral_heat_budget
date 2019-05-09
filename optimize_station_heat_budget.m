function stn = optimize_station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,opts)
%function stn = optimize_station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,opts)
%
% "Optimize" all params of reef ocean heat budget by graphing results.  Most
% args are optional and specify fieldname prefixes: e.g.,. RAPFX=='erai' uses
% data fields from ECMWF Reanalysis-Interim (ERAI) for radiative fluxes and
% meteorological inputs; ISPFX=='ndbc' uses in situ sea/air temp. and winds.
% See the documentation of function STATION_HEAT_BUDGET for more details.
%
% SAMPLE usages for the optional SUBSTITUTE_FIELD_NAMES cell array arg:
%  {'sfld','hourly_misst_sst'} evaluate heat budget based on MISST instead of in situ
%  {'reanalysis_shortwave',true} use upward shortwave from reanalysis instead of bulk
%  {'override_light_kds',{{0.2,[.05,.25,91]}}} run heat budget with various Kd values
%  {'override_advection_factors',{{0,[0,.5,0]}}} ditto various advection weightings
%  {'comparisonsfld','ndbc_sea_t'} compare result with actual in situ sea temperature
%
% Specifying STN.(COMPARISONSFLD) plots a sea temperature that may be
% different from STN.(SFLD) for comparison with the heat budget results.
%
% OPTS is a struct, containing name-value pairs for any named option
% recognized by the code that follows, or any functions called below.
%
% Last Saved Time-stamp: <Sun 2014-09-28 18:43:27 Eastern Daylight Time gramer>


  % Save struct STN to a (*large*) MAT file before finishing?
  doSave = false;

  % Plot (multiple and) various results for the caller?
  doPlot = true;
  doBoxplots = doPlot;
  doClimplots = doPlot;
  doDiagplots = doPlot;
  doAnnSubs = doPlot;
  doScatterplots = doPlot;

  % Print resulting plots to individual TIFF files?
  doPrint = doPlot;


  doBoxplots = false;
  doClimplots = false;
  doDiagplots = false;
  % Comment out to enable year-by-year comparison plots
  doAnnSubs = false;
  doScatterplots = false;

  doPrint = false;


  set_more off
  timenow,
  tic,


  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  tmp.stn_or_stnm = stn_or_stnm;
  station_heat_budget_field_names;
  % HACK to preserve buoy numbers from FIX_VARNAMELENS (v.)
  stn_or_stnm = tmp.stn_or_stnm; clear tmp;

  % If COMPARISONFLD not set, compare sea temperature with its own heat budget 
  if ( ~exist('comparisonsfld','var') || isempty(comparisonsfld) )
    comparisonsfld = sfld;
  end;


  datapath = get_thesis_path('../data');
  figspath = get_thesis_path('../figs');


  stn = get_station_from_station_name(stn_or_stnm);
  stn.commentstr = '';
  stnm = lower(stn.station_name);
  STNM = upper(stn.station_name);

  % Save information about what inputs we used
  stn.RAPFX=RAPFX;
  stn.KMPFX=KMPFX;
  stn.ISPFX=ISPFX;
  stn.TIDEPFX=TIDEPFX;
  stn.WAVEPFX=WAVEPFX;
  stn.substitute_field_names = substitute_field_names;
  stn.QEPFX=QEPFX;
  stn.DTPFX=DTPFX;

  stn.sfld  = sfld;
  stn.qfld  = sqtfld;
  stn.bqfld = bq0tfld;
  stn.dtfld = bdTfld;
  stn.hcfld = hcdTdt;


  oldstr = '';
  if ( use_old_options )
    stn.commentstr = [stn.commentstr,' (OLD OPTIONS)'];
    oldstr = '-OLD';
  end;
  if ( use_old_erai_only_options )
    stn.commentstr = [stn.commentstr,' (OLD ERAI OPTIONS)'];
    oldstr = '-OLDERA';
  end;


  % Default heat budget options

  if ( ~isfield(stn,'opts') )
    stn.opts = [];
  end;
  if ( exist('opts','var') )
    for optscfld=fieldnames(opts)'
      stn.opts.(optscfld{:}) = opts.(optscfld{:});
    end;
  end;
  stn.opts = station_heat_budget_options(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names);


  % Stations where no high-resolution bathymetry is available
  switch ( stnm ),
   case '42003',
    disp(['** Using zero mean tide height anomaly (no bathymetry available) **']);
    bathyfld = 0;
  end;


  % Number of points to use in finite-difference templates for gradients
  npts = 5;
  switch ( STNM ),
   case {'SANF1','DRYF1'},
    % Topography around some stations will not allow larger findiff templates
    npts = 3;
  end;


  % Filter out dates flagged as "bad" in "<stnm>-bad-dates.csv" file?
  stn.opts.keep_bad_dates = get_opt(stn.opts,'keep_bad_dates',false);
  if ( use_old_erai_only_options )
    stn.opts.keep_bad_dates = true;
  end;

  stn.opts.default_salinity = get_opt(stn.opts,'default_salinity',35.5);
  sal = stn.opts.default_salinity;
  %% SENSITIVITY ANALYSIS 
  % sal = stn.opts.default_salinity - 0.5;
  % sal = stn.opts.default_salinity + 0.5;

  %%%% ??? DEBUG Try low-frequency ADCP currents for advection?
  stn.opts.force_adcp_currents = get_opt(stn.opts,'force_adcp_currents',false);
  if ( stn.opts.force_adcp_currents )
    stn.commentstr = [stn.commentstr,' (ADCP)'];
    %%%%DEBUG???
    % stn.opts.add_alongshore_advection = false;
  end;

  stn.opts.advection_factor_ts = get_opt(stn.opts,'advection_factor_ts',empty_ts);

  stn.opts.grid_interp_method = get_opt(stn.opts,'grid_interp_method','linear');
  disp(['** Using ',stn.opts.grid_interp_method,' grid interpolation **']);


  stn.opts.force_reanalysis_winds = get_opt(stn.opts,'force_reanalysis_winds',false);
  if ( stn.opts.force_reanalysis_winds )
    disp('%% Forcing use of ERAI winds **');
    stn.commentstr = [stn.commentstr,' (EWnd)'];
    Wfld='erai_wind_speed';
    Dfld='erai_wind_dir';
  end;

  stn.opts.force_bic_insolation = get_opt(stn.opts,'force_bic_insolation',false);
  if ( stn.opts.force_bic_insolation )
    disp('%% Forcing use of MLRF2 in situ PAR -> insolation **');
    stn.commentstr = [stn.commentstr,' (BIC)'];
    dsrfld='bic_surf_dsrf';
  end;

  if ( ~strcmp(sfld,[ISPFX '_sea_t']) )
    disp(['Using sea temperature ',sfld]);
  end;
  if ( ~strcmp(afld,[ISPFX '_air_t']) )
    disp(['Using air temperature ',afld]);
  end;
  if ( ~strcmp(sfld,comparisonsfld) )
    disp(['COMPARISON sea temperature for plots ',comparisonsfld]);
  end;

  if ( adjust_waves )
    disp(['** Will apply empirical adjustment to ',upper(WAVEPFX),' wave data **']);
    stn.commentstr = [stn.commentstr,' (Wv adj) '];
  end;
  if ( adjust_reanalysis )
    disp(['** Will apply empirical adjustment to ',upper(RAPFX),' reanalysis **']);
    stn.commentstr = [stn.commentstr,' (RA adj) '];
  end;

  if ( ignore_benthos )
    disp('** Ignoring absorption and benthic exchange **');
    stn.commentstr = [stn.commentstr,' (NO Q_b) '];
  end;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Load basic station data

  if ( ~isnumeric(bathyfld) && ~isfield(stn,bathyfld) )
    stn = get_ngdc_bathy_station(stn);
  end;
  if ( ~isfield(stn,slopefld) )
    stn = station_ngdc_offshore_slope(stn);
  end;
  if ( ~isfield(stn,bathorifld) )
    stn = station_optimal_isobath_orientation(stn);
  end;

  warning('OFF','Ecoforecasts:mergedNonTS');
  stn = load_all_ndbc_data(stn);
  warning('ON','Ecoforecasts:mergedNonTS');



  %%%
  %% Sea temperature (in situ NDBC, site-specific, remote sensing, or reanalysis)

  if ( ~isfield(stn,sfld) )
    if ( ~isempty(regexp(sfld,'misst')) )
      stn = get_misst_station(stn);
      stn.hourly_misst = interp_ts(stn.misst_sst);
    elseif ( ~isempty(regexp(sfld,'avhrr_weekly_sst')) )
      stn = get_avhrr_weekly_field(stn,true,stn.opts.grid_interp_method,npts,stn.opts.keep_bad_dates);
    elseif ( ~isempty(regexp(sfld,'avhrr_sst')) )
      stn = get_avhrr_field(stn,true,[],stn.opts.grid_interp_method);
      % % Re-interpolate gradient and Laplacian time series fields
      % stn.(hkmtxfld) = interp_ts(stn.(kmtxfld),inf);
      % stn.(hkmtyfld) = interp_ts(stn.(kmtyfld),inf);
      % stn.(hkmtlfld) = interp_ts(stn.(kmtlfld),inf);
      % Time series fields may contain different numbers of (non-NaN) points
      [stn.(hkmtxfld),stn.(hkmtyfld)] = intersect_tses(stn.(hkmtxfld),stn.(hkmtyfld));
    elseif ( ~isempty(regexp(sfld,'erai_')) )
      if ( ~strcmpi(ISPFX,'erai') && ~strcmpi(RAPFX,'erai') )
        stn = get_erai_station(stn);
        if ( adjust_reanalysis )
          stn = adjust_erai_station(stn);
        end;
      end;
    elseif ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      stn = get_looe1_microcat(stn);
      if ( ~isempty(regexp(sfld,'mc_')) && ~isfield(stn,sfld) )
        warning('Interpolating SFLD "%s" from STN.microcat_seatemp',sfld);
        stn.(sfld) = interp_ts(stn.microcat_seatemp);
      end;
    elseif ( ~isempty(regexp(sfld,'adcp')) || ~isempty(regexp(sfld,'ad_')) )
      stn = get_looe1_adcp(stn);
      if ( ~isempty(regexp(sfld,'ad_')) && ~isfield(stn,sfld) )
        warning('Copying SFLD "%s" from STN.adcp_seatemp',sfld);
        stn.(sfld) = stn.adcp_seatemp;
      end;
    elseif ( ~isempty(regexp(sfld,'ndbc')) )
      if ( ~strcmpi(ISPFX,'ndbc') )
        stn = load_all_ndbc_data(stn);
      end;
    %elseif ( ~isempty(regexp(sfld,'\<(sea|ct_|ctd_)')) )
    elseif ( ~isempty(regexp(sfld,'^(ct_|ctd_)')) )
      if ( ~strcmpi(ISPFX,'icon') )
        stn = load_station_data(stn);
      end;

    elseif ( ~isempty(regexp(sfld,'seatemp')) )
      switch ( stnm ),
       case {...
           'tavrk',...
           'consh',...
           'condp',...
           'bnpin',...
           'bnpmi',...
           'bnppa',...
           'bnpon',...
           'bnpnn',...
            },
        stns = get_lirman_thermistors;
        stn.seatemp = stns.(stnm).seatemp;
        stns=[]; clear stns;
       otherwise,
        warning('Not sure how to load STN.seatemp for site "%s"',stnm);
      end;
    end;

    if ( ~isempty(regexp(sfld,'merged_sea_t')) )
      disp('** Merging NDBC_SEA_T w/coincident CT_SHALLOW_SEATEMP **');
      % Want to "filter bad dates" out of order here, before merging...
      stn = station_filter_bad_dates(stn);
      x=load_station_data(stn.station_name);
      if ( ~isfield(stn,'ndbc_sea_t') || ~isfield(x,'ct_shallow_seatemp') )
        error('Cannot create MERGED_SEA_T: missing NDBC_SEA_T or CT_SHALLOW_SEATEMP?');
      end;
      dts=x.ct_shallow_seatemp.date;
      dat.merged_sea_t.date=datenum(get_year(dts),get_month(dts),get_dom(dts),get_hour(dts),0,0);
      dat.merged_sea_t.data=x.ct_shallow_seatemp.data;
      stn.merged_sea_t = stn.ndbc_sea_t;
      warning('OFF','Ecoforecasts:mergedNonTS');
      stn = merge_station_data(stn,dat);
      warning('ON','Ecoforecasts:mergedNonTS');
      x=[]; dts=[]; dat=[]; clear x dts dat
    end;

  end;


  if ( ~strcmp(sfld,comparisonsfld) && ~isfield(stn,comparisonsfld) )
    if ( ~isempty(regexp(comparisonsfld,'misst')) )
      stn = get_misst_station(stn);
      stn.hourly_misst = interp_ts(stn.misst_sst);
    elseif ( ~isempty(regexp(comparisonsfld,'avhrr_weekly_sst')) )
      stn = get_avhrr_weekly_field(stn,true,stn.opts.grid_interp_method,npts,stn.opts.keep_bad_dates);
    elseif ( ~isempty(regexp(comparisonsfld,'avhrr_sst')) )
      stn = get_avhrr_field(stn,true,[],stn.opts.grid_interp_method);
      % % Re-interpolate gradient and Laplacian time series fields
      % stn.(hkmtxfld) = interp_ts(stn.(kmtxfld),inf);
      % stn.(hkmtyfld) = interp_ts(stn.(kmtyfld),inf);
      % stn.(hkmtlfld) = interp_ts(stn.(kmtlfld),inf);
      % Time series fields may contain different numbers of (non-NaN) points
      [stn.(hkmtxfld),stn.(hkmtyfld)] = intersect_tses(stn.(hkmtxfld),stn.(hkmtyfld));
    elseif ( ~isempty(regexp(comparisonsfld,'microcat')) || ~isempty(regexp(comparisonsfld,'mc_')) )
      stn = get_looe1_microcat(stn);
      if ( ~isempty(regexp(comparisonsfld,'mc_')) && ~isfield(stn,comparisonsfld) )
        warning('Interpolating COMPARISONSFLD "%s" from STN.microcat_seatemp',comparisonsfld);
        stn.(comparisonsfld) = interp_ts(stn.microcat_seatemp);
      end;
    elseif ( ~isempty(regexp(comparisonsfld,'adcp')) || ~isempty(regexp(comparisonsfld,'ad_')) )
      stn = get_looe1_adcp(stn);
      if ( ~isempty(regexp(comparisonsfld,'ad_')) && ~isfield(stn,comparisonsfld) )
        warning('Copying COMPARISONSFLD "%s" from STN.adcp_seatemp',comparisonsfld);
        stn.(comparisonsfld) = stn.adcp_seatemp;
      end;
    elseif ( ~isempty(regexp(comparisonsfld,'ndbc')) )
      if ( ~strcmpi(ISPFX,'ndbc') )
        stn = load_all_ndbc_data(stn);
      end;
    %elseif ( ~isempty(regexp(comparisonsfld,'\<(sea|ct_|ctd_)')) )
    elseif ( ~isempty(regexp(comparisonsfld,'^(ct_|ctd_)')) )
      if ( ~strcmpi(ISPFX,'icon') )
        stn = load_station_data(stn);
      end;

    elseif ( ~isempty(regexp(comparisonsfld,'seatemp')) )
      switch ( stnm ),
       case {...
           'tavrk',...
           'consh',...
           'condp',...
           'bnpin',...
           'bnpmi',...
           'bnppa',...
           'bnpon',...
           'bnpnn',...
            },
        stns = get_lirman_thermistors;
        stn.seatemp = stns.(stnm).seatemp;
        stns=[]; clear stns;
       otherwise,
        warning('Not sure how to load STN.*seatemp* for site "%s"',stnm);
      end;
    end;
  end;

  if ( stn.opts.force_adcp_currents )
    if ( ~isfield(stn,'adcp_u') || ~isfield(stn,'adcp_v') )
      disp(['Processing ADCP currents ** from LOOE1 ** for advection']);
      x = get_looe1_adcp;
      stn.adcp_u = x.adcp_u;
      stn.adcp_v = x.adcp_v;
      x=[]; clear x;
    else
      disp(['Processing ADCP currents for advection']);
    end;
  end;


  %%%
  %% Meteorological variables (in situ, reanalysis, or a mix)

  warning('OFF','Ecoforecasts:mergedNonTS');
  if ( ~isfield(stn,afld) || ~isfield(stn,sfld) || ~isfield(stn,Wfld) )
    switch (ISPFX),
     case 'ndbc',	%stn = load_all_ndbc_data(stn);
     case 'erai',	stn = get_erai_station(stn);
      if ( adjust_reanalysis )
        stn = adjust_erai_station(stn);
      end;
     case 'icon',	stn = load_station_data(stn);
     case 'ncep',		stn = get_ncep_station(stn,'narr');
      disp('** Using default APBL height 600m with NCEP NARR **');
      stn.(pblzfld) = stn.(dsrfld);
      stn.(pblzfld).data(:) = 600;
     otherwise,		error('Unknown in situ dataset "%s"',ISPFX);
    end;
  end;
  warning('ON','Ecoforecasts:mergedNonTS');

  %% Apply additional processing (moving average, low-pass or QC filter, etc.)
  stn = verify_variable(stn, afld);
  stn = verify_variable(stn, pfld);
  stn = verify_variable(stn, Wfld);
  stn = verify_variable(stn, sfld);

  if ( ~isfield(stn,sfld) )
    error(['Unknown sea temperature field ',sfld]);
  end;
  if ( ~isfield(stn,afld) )
    error(['Unknown air temperature field ',afld]);
  end;


  %% Filter out periods when the heat budget does NOT work well...
  % These are generally dates when variability in in situ data is suspect

  %bad_years=[];
  stn.opts.bad_years = get_opt(stn.opts,'bad_years',[]);
  bad_years=stn.opts.bad_years;

  if ( isempty(bad_years) )
    % NOTE: Per-station defaults
    switch ( lower(stn.station_name) ),
     case 'looe1',
      if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      end;
     case 'mlrf1',
      if ( ~isempty(regexp(sfld,'sea_t')) )
      end;
     case 'smkf1',
      if ( ~isempty(regexp(sfld,'sea_t')) )
        bad_years=[2004,2008];
      end;
    end;
  end;

  if ( ~isempty(find(ismember(get_year(stn.(sfld).date),bad_years))) )
    warning('** Removing %s year(s) %s',sfld,num2str(bad_years(:)'));
    stn.commentstr = [stn.commentstr,' (sans ',num2str(bad_years(:)'),') '];
    for badyr = bad_years(:)';
      stn.(sfld).data(get_year(stn.(sfld).date)==badyr)=[];
      stn.(sfld).date(get_year(stn.(sfld).date)==badyr)=[];
    end;
  end;



  %% Additional QA/QC: Remove any user-defined "bad dates" from time series
  if ( ~stn.opts.keep_bad_dates )
    stn = station_filter_bad_dates(stn);
  end;

  if ( ~isfield(stn,qsfld) )
    % Saturated ("sea-surface") specific humidity [kg/kg]
    stn = station_relhumid_to_spechumid(stn,sfld,100,qsfld);
    % Adjustment for salinity with thanks to Stommel
    stn.(qsfld).data = 0.98 .* stn.(qsfld).data;
  end;

  %% Low-pass filter winds for quasi-Eulerian currents
  stn = verify_variable(stn,Ulpfld);
  stn = verify_variable(stn,Vlpfld);
  stn.(Wlpfld).date = stn.(Ulpfld).date;
  stn.(Wlpfld).data = uv_to_spd(stn.(Ulpfld).data,stn.(Vlpfld).data);
  stn.(Dlpfld).date = stn.(Ulpfld).date;
  stn.(Dlpfld).data = uv_to_dir(stn.(Ulpfld).data,stn.(Vlpfld).data);


  stn = station_tmd_tide(stn);

  h_override = get_opt(stn.opts,'override_tide_height',[]);
  if ( ~isempty(h_override) )
    stn.(hfld).data(:) = h_override;
  end;

  %stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld,tufld,tvfld);
  stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
% %%%%%%%%%%???DEBUG
% warning('EXPERIMENT WITH REDUCED DEPTH');
%   stn = station_mean_tide_height(stn,mhfld,stn.depth,hfld);
% %%%%%%%%%%???DEBUG


  mh_override = get_opt(stn.opts,'override_mean_tide_height',[]);
  if ( ~isempty(mh_override) )
    stn.(mhfld).data(:) = mh_override;
  end;

  %% SENSITIVITY ANALYSIS 
  % stn.(mhfld).data(:) = stn.depth;

  % Assume max warm-layer depth somewhere near bottom boundary layer top
  max_wl = nanmax(stn.(mhfld).data) - 0.5;


  %%%
  %% Radiative fluxes (in situ, reanalysis, or both)

  if ( ~isfield(stn,dsrfld) || ~isfield(stn,dlrfld) || ~isfield(stn,cfld) )
    switch (RAPFX),
     case 'erai',		stn = get_erai_station(stn);
      if ( adjust_reanalysis )
        stn = adjust_erai_station(stn);
      end;
     case 'ncep',		stn = get_ncep_station(stn,'narr');
      disp('** Using default APBL height 600m with NCEP NARR **');
      stn.(pblzfld) = stn.(dsrfld);
      stn.(pblzfld).data(:) = 600;
     %%% Not Yet Implemented
     % case 'era40',		stn = get_era40_station(stn);
     % case 'cfsr',		stn = get_ncep_station(stn,'cfsr');
     otherwise,		error('Unavailable gridded/reanalysis dataset "%s"',RAPFX);
    end;
    stn = station_heat_flux_term(stn,raq0fld,raqtfld,sfld,sal,nanmean(stn.(mhfld).data));
  end;


  %%%
  %% If caller wishes to use a different source than RAPFX for specific humidity
  warning('OFF','Ecoforecasts:mergedNonTS');
  if ( ~isfield(stn,qafld) )
    switch (qafld),
     case 'ndbc_spechumid',
      if ( ~isfield(stn,'ndbc_air_t') )
        stn = load_all_ndbc_data(stn);
        if ( ~isfield(stn,'ndbc_dew_t') )
          error('No dewpoint temperature to calculate ndbc_spechumid!');
        end;
        stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
        stn = station_relhumid_to_spechumid(stn,'ndbc_air_t','ndbc_relhumid',qafld);
      end;
     case 'ncep_spechumid',
      x = get_ncep_station(stnm,'narr');
      stn.(qafld) = x.(qafld);
      x=[]; clear x;
     otherwise,
      error('Unknown specific humidity field "%s"',qafld);
    end;
    disp(['** Using Specific Humidity field ',qafld,' **']);
  end;
  warning('ON','Ecoforecasts:mergedNonTS');


  %%%
  %% Actual sea temperature change, and implied flux
  if ( ~isfield(stn,dsfld) )
    stn.(dsfld).date = stn.(sfld).date(2:end);
    stn.(dsfld).data = diff(stn.(sfld).data);
    stn = filter_gaps(stn,sfld,dsfld,(1.5/24));
  end;
  if ( ~isfield(stn,dsffld) )
    stn = station_heat_flux_term_inverse(stn,dsffld,dsfld,sfld,sal,mhfld);
  end;


  %%%
  %% Change in daily mean sea temperature, and implied flux
  if ( ~isfield(stn,dsfld) )
    stn.(dsfld).date = stn.(sfld).date(2:end);
    stn.(dsfld).data = diff(stn.(sfld).data);
    stn = filter_gaps(stn,sfld,dsfld,(1.5/24));
  end;
  if ( ~isfield(stn,dsffld) )
    stn = station_heat_flux_term_inverse(stn,dsffld,dsfld,sfld,sal,mhfld);
  end;



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Optimization parameters
  %%  The ranges defined just below are for use with a new station, when we
  %%  are "searching blind" for the right mix of parameter values.

  % "Fudge factor" for Latent Heat Flux
  qlh_adj = 1.0;

  % Apply TOGA-COARE 3.0a Warm-layer adjustment for sea temp. sensor depth?
  doWarms = [false, true];
  % Attentuation coefficient/climatology for "penetrative" (PAR+NUV) insolation
  kds = {...
      0.1,[0.05,0.15,45],[0.05,0.15,137],[0.05,0.15,228],[0.05,0.15,320],...
      0.2,[0.10,0.30,45],[0.10,0.30,137],[0.10,0.30,228],[0.10,0.30,320],...
      0.3,[0.20,0.40,45],[0.20,0.40,137],[0.20,0.40,228],[0.20,0.40,320],...
        };

  % Bulk benthic heat coefficient for flow over sea-floor
  % % %cbds = {0,3.8e-5,2.9e-4,3.8e-4,8.0e-4,16.0e-4,24.0e-4};
  % % Cbd~CD^2, based on estimate CD~0.017 in Davis and Monismith (2011)
  % cbds = {2.9e-4};
  % Davis and Monismith (2011) CD times Cbh, where Cbh~1e-3
  cbds = {0.017*1e-2};

  % Heat advection fudge factor (0=no advection)
  advfacs = {...
      0,...
      1,...
            };
  % Heat diffusion fudge factor (0=no diffusion)
  kths = {...
      0 ,...
      5, [0 ,10, 45],...
      10, [0 ,20, 45],...
         };


  switch ( STNM ),

   case 'LKWF1',
    qlh_adj = 0.90;
    doWarms = [false];
    kds = { ...
        [0.050,0.400, 67] ...		% BEST SO FAR:  0.4,3.5,1.0
          };
    advfacs = { ...
        [0.00,1.0, 45] ...
              };
    kths = { ...
        [0,20, 45] ...
           };


   case 'FWYF1',
    if ( use_old_erai_only_options )
      stn.opts.Ppen = 0.473;
      cbds = {8e-4};
      doWarms = [false];
      kds = { ...
          [0.035,0.220, 76] ...
            };
      advfacs = { ...
          [0.0,1.0,137] ...
                };
      kths = { ...
          0 ...
             };
      % Essentially TURN OFF outlier removal
      qlh_err_cutoff = 100000;
      bq0_err_cutoff = 100000;
      hc_err_cutoff = 100.00;
    elseif ( use_old_options )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          %[0.025,0.300,110] ...     % Best for no waves 1987-2011 and(?) WW3 1999-2011
          [0.025,0.300,137] ...     % Best for ERAI waves 1993-2011
            };
      advfacs = { ...
          [0.00,1.00, 91] ...
                };
      kths = { ...
          [0,10, 91] ...
             };
    else
      doWarms = [true];
      kds = { ...
          [0.025,0.225,102] ...     % Best with hc_scaling=SU,hc_warming_factor=0.66
            };
      advfacs = { ...
          [0.00,0.50,102] ...
                };
      kths = { ...
          [0,20, 45] ...
             };
    end;


   case 'MLRF1',
    if ( use_old_erai_only_options )
      stn.opts.Ppen = 0.473;
      cbds = {8e-4};
      doWarms = [false];
      kds = { ...
          [0.035,0.190, 55] ...
          % [0.050,0.350, 91] ... %Best for new options
          % [0.045,0.375, 45] ... %ERAI met, air_t, old Ppen
            };
      advfacs = { ...
          1 ...
                };
      kths = { ...
          0 ...
             };
      % Essentially TURN OFF outlier removal
      qlh_err_cutoff = 100000;
      bq0_err_cutoff = 100000;
      hc_err_cutoff = 100.00;
    elseif ( use_old_options )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          % [0.075,0.275,110] ...     % Best for WW3 1999-2011?
          [0.050,0.350, 91] ... %Best so far
            };
%       kds = { ...
%           % [0.075,0.275,110] ...     % Best for WW3 1999-2011?
%           [0.050,0.350, 91] ... %Best so far
%           % [0.050,0.200,288] ... % Best match w/1m post-cleaning measurements
%           % [0.350,0.700,288] ... % Best(?) for mean depth 3.55m
%           % [0.350,0.700, 91] ... % Best(??) for mean depth 3.55m
%             };
% %%%%%%%%%%???DEBUG
% warning('EXPERIMENT WITH REDUCED DEPTH');
% %%%%%%%%%%???DEBUG
      advfacs = { ...
          [0.0,1.0, 45] ... %Best so far
                };
      kths = { ...
          [0,20, 45] ... %Best so far
             };
    else
      doWarms = [true];
      kds = { ...
          % %[.038,.250, 67] ...     % Best with hc_scaling=SU,hc_warming_factor=0.66
          % [.035,.250, 70] ...     % Best with TIDAL ADVECTION, hc_scaling=SU,hc_warming_factor=0.66
          [.035,.250, 69] ...     % Best with TIDAL ADVECTION, hc_scaling=SU,hc_warming_factor=0.66, consistent w/SMKF1
            };
      advfacs = { ...
          [0.00,1.0, 45] ...
                };
      kths = { ...
          [0,20, 45] ...
             };
    end;

   case 'LONF1',
    if ( use_old_erai_only_options )
      stn.opts.Ppen = 0.473;
      cbds = {8e-4};
      doWarms = [false];
      kds = { ...
          [0.200,0.370,228] ...
            };
      advfacs = { ...
          1 ...
                };
      kths = { ...
          [0,5,183] ...
             };
      % Essentially TURN OFF outlier removal
      qlh_err_cutoff = 100000;
      bq0_err_cutoff = 100000;
      hc_err_cutoff = 100.00;
    elseif ( use_old_options )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      %stn.opts.hc_scaling = 'SU';
      %stn.commentstr = [stn.commentstr,' (',stn.opts.hc_scaling,') '];
      doWarms = [false];
      kds = { ...
          % [0.700,1.200, 20] ... % Ended closest to balance
          [0.675,1.250, 10] ... % Best for ERAI waves adv=[0,1,91],K=[0,2,91]
            };
      advfacs = { ...
          [0.0,1.0, 91] ... % Best so far
                };
      kths = { ...
          [0, 2, 91] ... % We want SOME diffusion if we can have it
             };
    else
      doWarms = [false];
      kds = { ...
          %[.500,1.250, 45] ...     % Best with hc_scaling=SU,hc_warming_factor=0.66
          [.475,1.275, 45] ...     % Best with hc_scaling=SU,hc_warming_factor=0.66
            };
      advfacs = { ...
          [0.00,1.0, 45] ...
                };
      kths = { ...
          [0,20, 45] ...
             };
    end;


   case 'SMKF1',
    if ( use_old_erai_only_options )
      stn.opts.Ppen = 0.473;
      cbds = {8e-4};
      doWarms = [false];
      kds = { ...
          [0.045,0.375, 45] ...
            };
      advfacs = { ...
          1 ...
                };
      kths = { ...
          0 ...
             };
      % Essentially TURN OFF outlier removal
      qlh_err_cutoff = 100000;
      bq0_err_cutoff = 100000;
      hc_err_cutoff = 100.00;
    elseif ( use_old_options )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      %stn.opts.hc_scaling = 'SU';
      %stn.commentstr = [stn.commentstr,' (',stn.opts.hc_scaling,') '];
      doWarms = [true];
      kds = { ...
          % [0.100,0.400,110] ...    % "Best" for *all* HC scalings [US][SU]
          [0.100,0.600, 70] ...    % Best for NO GRADIENT case - Best so far?
            };
      advfacs = { ...
          [0.0,0.25, 91] ...
                };
      kths = { ...
          [0, 2, 91] ...
             };
    else
      doWarms = [true];
      kds = { ...
          [0.066,0.450, 69] ...
            };
      advfacs = { ...
          [0.0,1.00, 45] ...
                };
      kths = { ...
          [0,20, 45] ...
             };
    end;

   case 'SANF1',
    if ( use_old_options )
      % Number of points to use in finite-difference templates for gradients
      %%%%npts = 3; %%%%MUST BE SET ABOVE!
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      % So-so RMSE 5.5 (best was 2.1), but best *climatological* RMSE (0.8)
      kds = { ...
          ... %From *_options.m: [0.066,0.450, 69]
          [0.010,0.300, 81] ... %Best so far
            };
      advfacs = { ...
          0 ...
                };
      kths = { ...
          [0, 2, 91] ... %Best so far
             };
    else
      doWarms = [true];
      kds = { ...
          [0.015,0.150, 67] ...
            };
      advfacs = { ...
          0 ...
                };
      kths = { ...
          [0,20, 45] ...
             };
    end;

   case 'DRYF1',
    if ( use_old_options )
      % Number of points to use in finite-difference templates for gradients
      %%%%npts = 3; %%%%MUST BE SET ABOVE!
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      %%%%???DEBUG
      % stn.opts.hc_scaling = 'SU';
      % stn.commentstr = [stn.commentstr,' (',stn.opts.hc_scaling,') '];
      %%%%???DEBUG
      % qlh_adj = 1.1;
      % doWarms = [false];
      %%%%???DEBUG
      doWarms = [true];
      kds = { ...
        % [0.300,0.700,342] ...
        % [0.300,0.700,  0] ... %Best so far?
        % [0.250,0.750,  0] ...
        % [0.200,0.800,  0] ...
        [0.250,0.750,342] ... %Best so far?
        % [0.200,0.800,342] ...
            };
      advfacs = { ...
          [0.00,1.00, 45] ... %Best so far
                };
    kths = { ...
        % [0, 2, 45] ...
        [0,20, 45] ...
           };

    else
      qlh_adj = 0.90;
      doWarms = [false];
      kds = { ...
          [0.150,0.500,354] ...
            };
      advfacs = { ...
          [0.00,1.00, 45] ...
                };
      kths = { ...
          [0,20, 45] ...
             };
    end;


   case 'LOOE1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      if ( use_old_options )
        % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
        doWarms = [true];
        kds = { ...
            [0.050,0.400,137] ... %BEST FOR 0, 0, all years
              };
        advfacs = { ...
            [0.00,0.50, 45] ...
                  };
        kths = { ...
            0 ...                 %BEST
               };
      else
        doWarms = [true];
        kds = { ...
            [0.050,0.150, 45] ...
              };
        advfacs = { ...
            [0.00,0.25, 45] ...
                  };
        kths = { ...
            [0,20, 45] ...
               };
      end;
    elseif ( ~isempty(regexp(sfld,'adcp')) || ~isempty(regexp(sfld,'ad_')) )
      if ( use_old_options )
        % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
        doWarms = [true];
        kds = { ...
            [0.010,0.400,113] ... %BEST SO FAR
              };
        advfacs = { ...
            [0.00,0.50, 45] ...
                  };
        kths = { ...
            [0,20, 45] ...
               };
      else
        doWarms = [true];
        kds = { ...
            [0.025,0.200, 80] ...
              };
        advfacs = { ...
            [0.00,0.25, 45] ...
                  };
        kths = { ...
            [0,20, 45] ...
               };
      end;
    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;

   case 'HAWK1',
    error('Hawk Channel SFP station HAWK1 *not yet implemented*');

   case {...
       'NCORA',
       'NCORB',
       'NCOR1',
       'NCOR2',
       'NCOR3',
       'NCOR4',
        },
    error('NCORE inshore sites *not yet implemented*');

   case 'NCORC',
    error('NCORE site "C" *not yet implemented*');

   case {'KLGF1','NCORK'},
    error('NCORE Key Largo site *not yet implemented*');

   case {'MRTF1','NCORM'},
    error('NCORE Marathon site *not yet implemented*');

   case {...
       'BNPIN',...
       'BNPMI',...
       'BNPON',...
       'BNPNN',...
       'BNPPA',...
       'TAVRK',...
       'CONSH',...
       'CONDP',...
        },
    % stn.opts.hc_scaling = 'SS';
    % stn.opts.hc_scaling = 'SU';
    % stn.opts.hc_scaling = 'UU';
    % stn.commentstr = [stn.commentstr,' (',strrep(stn.opts.hc_scaling,'_','\_'),') '];
    doWarms = [false];
    kds = { ...
        % [0.700,1.200, 20] ... % Ended closest to balance
        [0.675,1.250, 10] ... % Best for ERAI waves adv=[0,1,91],K=[0,2,91]
          };
    advfacs = { ...
        [0.0,1.0, 91] ... % Best so far
              };
    kths = { ...
        %[0, 2, 91] ... % We want SOME diffusion if we can have it
        [0,20, 91] ... % Never tested
           };

   case '42003',
    doWarms = [true];
    kds = { ...
        0.2 ...
          };
    advfacs = { ...
        1 ...
              };
    kths = { ...
        20 ...
           };

   otherwise,
    error('Station %s options not implemented yet!',STNM);
  end;


  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% SENSITIVITY ANALYSIS 

  % doWarms(end+1) = ~doWarms(1);
  % cbds{end+1} = cbds{1}-1e-4;
  % cbds{end+1} = cbds{1}+1e-4;
  % kds{end+1} = kds{1}-0.025;
  % kds{end+1} = kds{1}+0.025;
  % advfacs{end+1} = 0;
  % advfacs{end+1} = 1;
  % advfacs{end+1} = [0,1,45];
  % kths{end+1} = 0;
  % kths{end+1} = 20;
  % kths{end+1} = { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) };
  %

  cbds    = get_opt(stn.opts,'override_benthic_cbd',cbds);
  doWarms = get_opt(stn.opts,'override_do_warm_adj',doWarms);
  kds     = get_opt(stn.opts,'override_light_kds',kds);
  advfacs = get_opt(stn.opts,'override_advection_factors',advfacs);
  kths    = get_opt(stn.opts,'override_k_thetas',kths);

  stn.opts.Ppen = get_opt(stn.opts,'override_Ppen',[]);
  qlh_adj = get_opt(stn.opts,'override_qlh_adj',qlh_adj);

  hcs = get_opt(stn.opts,'override_hc_scaling',stn.opts.hc_scaling);
  if ( ~strcmpi(stn.opts.hc_scaling,hcs) )
    stn.opts.hc_scaling = hcs;
    stn.commentstr = [stn.commentstr,' (really ',stn.opts.hc_scaling,') '];
  end;

  stn.opts.keep_err_outliers = get_opt(stn.opts,'keep_err_outliers',false);
  if ( stn.opts.keep_err_outliers )
    % Essentially TURN OFF outlier removal
    qlh_err_cutoff = 100000;
    bq0_err_cutoff = 100000;
    hc_err_cutoff = 100.00;
  end;

  %% SENSITIVITY ANALYSIS 
  %%%%%%%%%%%%%%%%%%%%%%%%%


  % % default_cbd = 3.8e-4; %From literature
  % default_cbd = 8.0e-4; %Highest convergent value tried
  % To coincide with estimate CD~0.017 in Davis and Monismith (2011)
  default_cbd = 2.9e-4;
  [ig,default_cbdix] = min(abs([cbds{:}] - default_cbd));

  for ix=1:numel(kds)
    if ( iscell(kds{ix}) )
      stn.optim.kdstrs{ix} = char(kds{ix}{2});
    elseif ( isnumeric(kds{ix}) )
      stn.optim.kdstrs{ix} = num2str(kds{ix},'%g,');
    elseif ( is_ts(kds{ix}) )
      stn.optim.kdstrs{ix} = 'time series';
    else
      stn.optim.kdstrs{ix} = 'unknown option type?';
    end;
  end;
  for ix=1:numel(cbds)
    if ( iscell(cbds{ix}) )
      stn.optim.cbdstrs{ix} = char(cbds{ix}{2});
    elseif ( isnumeric(cbds{ix}) )
      stn.optim.cbdstrs{ix} = num2str(cbds{ix},'C_b_d=%g');
    elseif ( is_ts(cbds{ix}) )
      stn.optim.cbdstrs{ix} = 'time series';
    else
      stn.optim.cbdstrs{ix} = 'unknown option type?';
    end;
  end;
  for ix=1:numel(advfacs)
    if ( iscell(advfacs{ix}) )
      stn.optim.advfacstrs{ix} = char(advfacs{ix}{2});
    elseif ( isnumeric(advfacs{ix}) )
      stn.optim.advfacstrs{ix} = num2str(advfacs{ix},'%g,');
    elseif ( is_ts(advfacs{ix}) )
      stn.optim.advfacstrs{ix} = 'time series';
    else
      stn.optim.advfacstrs{ix} = 'unknown option type?';
    end;
  end;
  for ix=1:numel(kths)
    if ( iscell(kths{ix}) )
      stn.optim.kthstrs{ix} = char(kths{ix}{2});
    elseif ( isnumeric(kths{ix}) )
      stn.optim.kthstrs{ix} = num2str(kths{ix},'%g,');
    elseif ( is_ts(kths{ix}) )
      stn.optim.kthstrs{ix} = 'time series';
    else
      stn.optim.kthstrs{ix} = 'unknown option type?';
    end;
  end;



  %%%
  %% Radiative fluxes ADDITIONAL PROCESSING (should follow "options" setting above)

  if ( ~isfield(stn,dsrfld) )
    if ( ~isempty(regexp(dsrfld,'bic_surf')) )
      disp('** Using BIC-derived downward shortwave radiation **');
      x = load_station_data('mlrf1');
      stn.bic_surf_par = x.bic_surf_par;
      x=[]; clear x
      stn = station_par_to_insol(stn,'bic_surf_par','bic_surf_dsrf','bic_surf_usrf','bic_surf_srf');
      stn.commentstr = [stn.commentstr,' (BIC DSRF) '];
    else
      error(['Unknown insolation field ',dsrfld]);
    end;
  end;
  %DEBUG:  disp(reanalysis_shortwave);
  %DEBUG:  disp(usrfld);
  if ( ~isfield(stn,usrfld) )
    if ( reanalysis_shortwave )
      error('Reanalysis upward shortwave field "%s" not found',usrfld);
    end;
    stn.opts.albedo = get_opt(stn.opts,'albedo',[]);
    if ( ~isempty(stn.opts.albedo) )
      disp('** Using options-specified albedo **');
      stn.(albfld) = stn.(dsrfld);
      stn.(albfld).data(:) = build_clim_opt(stn.opts.albedo,'albedo',stn.(dsrfld).date);
      stn.commentstr = [stn.commentstr,' (option Alb) '];
    else
      disp('** Using bulk upward shortwave radiation **');
      stn = station_bulk_albedo(stn,albfld,Wfld,cfld);
      stn.commentstr = [stn.commentstr,' (bulk USRF) '];
    end; %if isempty albedo else

    stn.opts.albedo_increment = get_opt(stn.opts,'albedo_increment',0.0);
    albinc = stn.opts.albedo_increment;
    %% SENSITIVITY ANALYSIS 
    % albinc = stn.opts.albedo_increment - 0.01;
    % albinc = stn.opts.albedo_increment + 0.01;
    if ( albinc > 0 )
      disp(['** Albedo + ',num2str(albinc),'%! **']);
      stn.(albfld).data = stn.(albfld).data + (albinc/100);
      stn.commentstr = [stn.commentstr,' (+',num2str(albinc),'%) '];
    end;

    stn.(usrfld) = ts_op(stn.(dsrfld),stn.(albfld),'.*');
    stn.(usrfld).data = -stn.(usrfld).data;
    % Add down- and upward shortwave fluxes together
    stn.(srfld) = ts_op(stn.(dsrfld),stn.(usrfld),'+');
  end;

  stn = station_heat_flux_term(stn,srfld,srtfld,sfld,sal,mhfld);


  if ( ~isfield(stn,dlrfld) )
      error(['Unknown downward longwave field ',dlrfld]);
  end;
  if ( ~isfield(stn,ulrfld) )
    if ( reanalysis_longwave )
      error('Reanalysis upward longwave field "%s" not found',ulrfld);
    end;
    disp('** Using bulk upward longwave radiation **');
    %%%% Calculated after STATION_HEAT_FLUX call(s) below instead...
    % %stn = station_bulk_ulr(stn,sfld,ulrfld,dlrfld);
    % stn = station_bulk_ulr(stn,sfld,ulrfld);
    stn.commentstr = [stn.commentstr,' (bulk ULRF) '];
  end;

  %%%% Calculated after STATION_HEAT_FLUX call(s) below...
  % % Add down- and upward longwave fluxes together
  % stn.(lrfld) = ts_op(stn.(dlrfld),stn.(ulrfld),'-');
  % stn = station_heat_flux_term(stn,lrfld,[lrfld,'_term'],sfld,sal,mhfld);



  if ( ~isfield(stn,whfld) )
    % If waves not from reanalysis, user must want model (WaveWatch III) or wind estimate
    switch (WAVEPFX),
     case 'erai',	stn = get_erai_station(stn); %IF NOT ALSO OUR REANALYSIS DATASET
      if ( adjust_reanalysis )
        stn = adjust_erai_station(stn);
      end;
     case 'ww3',	stn = get_ww3_station(stn);
     case 'ndbc',       stn = station_wind_to_wave(stn,Wfld,Dfld,wpfld,whfld,wdfld);
     otherwise,		error('Unknown wave source "%s"',WAVEPFX);
    end;
  end;

  if ( adjust_waves )
    switch (WAVEPFX),
     case 'erai',	stn = adjust_erai_station_waves(stn);
     case 'ww3',	stn = adjust_ww3_station_waves(stn);
     case 'ndbc',	
      stn.(wpfld).data = ((stn.(wpfld).data.*1.00) + 0.0);
      %stn.(whfld).data = ((stn.(whfld).data.*0.25) + 0.3);
      %stn.(whfld).data = ((stn.(whfld).data.*0.314) + 0.34);
      stn.(whfld).data = ((stn.(whfld).data.*0.15) + 0.25);
      stn.(wdfld).data = ((stn.(wdfld).data.*1.00) + 0.0);
     otherwise,		error('Unknown wave source "%s"',WAVEPFX);
    end;
  end;


  % If we loaded this bogus "net" field for any reason above, recalculate it properly
  if ( isfield(stn,'erai_net_heat_flux') )
    stn.erai_turbulent_heat_flux = ts_op(stn.erai_latent_heat_flux,stn.erai_sensible_heat_flux,'+');
    stn.erai_radiative_heat_flux = ts_op(stn.erai_srf,stn.erai_lrf,'+');
    stn.erai_actual_net_heat_flux = ts_op(stn.erai_turbulent_heat_flux,stn.erai_radiative_heat_flux,'+');
  end;
  % Or do the same calculation for NCEP NARR if it was loaded
  if ( isfield(stn,'ncep_net_heat_flux') )
    stn.ncep_turbulent_heat_flux = ts_op(stn.ncep_latent_heat_flux,stn.ncep_sensible_heat_flux,'+');
    stn.ncep_radiative_heat_flux = ts_op(stn.ncep_srf,stn.ncep_lrf,'+');
    stn.ncep_actual_net_heat_flux = ts_op(stn.ncep_turbulent_heat_flux,stn.ncep_radiative_heat_flux,'+');
  end;

  stn.opts.ocean_current_wind_multiplier = get_opt(stn.opts,'ocean_current_wind_multiplier',[]);
  if ( ~isempty(stn.opts.ocean_current_wind_multiplier) )
    disp(['%% Ocean current ',num2str(stn.opts.ocean_current_wind_multiplier),' x ',Wlpfld,' **']);
    stn.commentstr = [stn.commentstr,' (WFrc)'];
    stn.(sssfld) = ts_op(ts_fun(stn.(Wlpfld),@kts2mps),stn.opts.ocean_current_wind_multiplier,'*');
    stn.(ssdfld) = stn.(Dlpfld);
    stn = station_spddir_to_uv(stn,sssfld,ssdfld,ssufld,ssvfld,true);
  else
    % NOTE: STOKES_DRIFT (v. STATION_STOKES_DRIFT) converts wind speed WLPFLD from [Kts] to [m/s]
    stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wlpfld,Dlpfld,whfld,wpfld,wdfld);
  end;

  %% SENSITIVITY ANALYSIS 
  % stn.(sssfld).date=stn.(Wlpfld).date;  stn.(ssdfld).date=stn.(Dlpfld).date;
  % stn.(ssufld).date=stn.(Ulpfld).date;  stn.(ssvfld).date=stn.(Vlpfld).date;
  % stn.(ssufld).data=0.01.*stn.(Ulpfld).data; stn.(ssvfld).data=0.01.*stn.(Vlpfld).data;
  % stn.(ssufld).data=0.02.*stn.(Ulpfld).data; stn.(ssvfld).data=0.02.*stn.(Vlpfld).data;
  % stn.(ssufld).data=0.03.*stn.(Ulpfld).data; stn.(ssvfld).data=0.03.*stn.(Vlpfld).data;
  % stn.(ssufld).data=kts2mps(stn.(ssufld).data); stn.(ssvfld).data=kts2mps(stn.(ssvfld).data);
  % stn.(ssdfld).data=stn.(Dlpfld).data; stn.(sssfld).date=uv_to_spd(stn.(ssufld).data,stn.(ssvfld).data);

  % Cross- and long-shore components of currents
  stn = station_reorient_vectors(stn,bathorifld,ssufld,ssvfld);


  %% Kilometer-scale Ocean Data

  if ( ~isfield(stn,ufld) || ~isfield(stn,vfld) || ~isfield(stn,Tfld) )
    switch (KMPFX),
     case 'fkeys_hycom',
      stn = get_fkeys_hycom(stn,[],[],[],[],stn.opts.grid_interp_method);
      % Cross- and long-shore components of currents
      stn = station_reorient_vectors(stn,bathorifld,ufld,vfld);
      stn.opts.km_scale_advection = true;
      stn.opts.calculate_advection = get_opt(stn.opts,'calculate_advection',true);
      stn.opts.calculate_diffusion = true;
     case 'gom_hycom',
      stn = get_gom_hycom(stn,[],[],[],[],stn.opts.grid_interp_method);
      % Cross- and long-shore components of currents
      stn = station_reorient_vectors(stn,bathorifld,ufld,vfld);
      stn.opts.km_scale_advection = true;
      stn.opts.calculate_advection = get_opt(stn.opts,'calculate_advection',true);
      stn.opts.calculate_diffusion = true;
     case {'avhrr_weekly','avhrr'},
      if ( ~isfield(stn,Tfld) )
        if ( strcmpi(KMPFX,'avhrr') )
          disp('Loading AVHRR SST instead of hydrodynamic model data...');
          stn = get_avhrr_field(stn,true,[],stn.opts.grid_interp_method);
          % % Re-interpolate gradient and Laplacian time series fields
          % stn.(hkmtxfld) = interp_ts(stn.(kmtxfld),inf);
          % stn.(hkmtyfld) = interp_ts(stn.(kmtyfld),inf);
          % stn.(hkmtlfld) = interp_ts(stn.(kmtlfld),inf);
          % Time series fields may contain different numbers of (non-NaN) points
          [stn.(hkmtxfld),stn.(hkmtyfld)] = intersect_tses(stn.(hkmtxfld),stn.(hkmtyfld));
        else
          disp('Loading AVHRR_WEEKLY SST instead of hydrodynamic model data...');
          stn = get_avhrr_weekly_field(stn,true,stn.opts.grid_interp_method,npts,stn.opts.keep_bad_dates);
        end;
      end;
      % If we have in situ currents, TRY them for km-scale advection!
      if ( stn.opts.force_adcp_currents )
        if ( strcmpi([KMPFX '_u'],ufld) && strcmpi([KMPFX '_v'],vfld) )
          ufld = 'adcp_u_40_h_lp';
          vfld = 'adcp_v_40_h_lp';
        end;
        stn = verify_variable(stn,{ufld,vfld});
        hufld = [ufld];
        hvfld = [vfld];
        stn.opts.km_scale_advection = true;
        stn.opts.calculate_advection = get_opt(stn.opts,'calculate_advection',true);
        stn.opts.calculate_diffusion = true;
        disp('** (INSITU+STOKES) ADVECTION **');
      else
        stn.opts.km_scale_advection = false;
        stn.opts.calculate_advection = get_opt(stn.opts,'calculate_advection',true);
        stn.opts.calculate_diffusion = true;
        disp('** ONLY STOKES ADVECTION **');
      end;
     case 'none',
      disp('Loading NO kilometer-scale data...');
      stn.opts.km_scale_advection = false;
      stn.opts.calculate_advection = false;
      stn.opts.calculate_diffusion = false;
      disp('** NO STOKES ADVECTION **');
     otherwise,
      error('Unknown km-scale data source "%s"',KMPFX);
    end;
    more off;
  end;


  if ( isfield(stn,ufld) )
    if ( ~isfield(stn,hufld) )
      % Spline-fit an hourly time series of mean currents to native data
      stn.(hufld) = interp_ts(stn.(ufld));
      stn.(hvfld) = interp_ts(stn.(vfld));
    end;

    stn.(qeufld) = ts_op(stn.(hufld),stn.(ssufld),'+');
    stn.(qevfld) = ts_op(stn.(hvfld),stn.(ssvfld),'+');
  else
    stn.(qeufld) = stn.(ssufld);
    stn.(qevfld) = stn.(ssvfld);
  end;
  % Cross- and long-shore components of currents
  stn = station_reorient_vectors(stn,bathorifld,qeufld,qevfld);

  stn.(netufld) = ts_op(stn.(tufld),stn.(qeufld),'+');
  stn.(netvfld) = ts_op(stn.(tvfld),stn.(qevfld),'+');
  % Cross- and long-shore components of currents
  stn = station_reorient_vectors(stn,bathorifld,netufld,netvfld);


  if ( ~isfield(stn,Tfld) )
    if ( ~strcmpi(KMPFX,'none') )
      warning('No sea-surface temperature field %s',Tfld);
    end;
  elseif ( ~isfield(stn.(Tfld),'gradient_x') || ~isfield(stn,kmtxfld) )
    % Calculate gradients and field Laplacians - and interpolate to site
    stn = calc_field_terms(stn,Tfld,kmtfld,stn.opts.grid_interp_method,stn.lat,stn.lon,npts);
  end;

  % Precalculate full-scale advection (quality factor Fq is applied later)
  if ( ~isfield(stn,hkmtxfld) )
    warning('** No km-scale sea temperature gradients **');
    stn.(udTfld) = stn.(sfld);
    stn.(udTfld).data(:) = 0;
    stn.(udTffld) = stn.(udTfld);

  else

    if ( ~isfield(stn,hkmtxsfld) )
      stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',Inf);
      maxGrad = stn.opts.max_sst_gradient;
      badix = find(abs(stn.(hkmtxfld).data)>maxGrad | abs(stn.(hkmtyfld).data)>maxGrad);
      if ( ~isempty(badix) )
        disp(['** Zeroing unphysical gradients (|dSST/dx,dy|>',num2str(maxGrad*1e3),'K/km) **']);
        %tic,
        zeroix = [];
        %for ix=badix(:)'; zeroix = [zeroix ix-(24*7)+1:ix+(24*7)-1]; end;
        for ix=-(24*7)+1:(24*7)-1; zeroix = [zeroix badix(:)+ix]; end;
        zeroix = unique(zeroix);
        stn.(hkmtxfld).data(zeroix) = 0;
        stn.(hkmtyfld).data(zeroix) = 0;
        disp(['Zeroed ',num2str(numel(zeroix)),' points']);
        %toc,
      end;

      % Cross- and long-shore components of sea temperature gradient
      stn = station_reorient_vectors(stn,bathorifld,hkmtxfld,hkmtyfld);
    end;

    % Include tidal currents in advection estimate?
    stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',false);
    if( stn.opts.tidal_advection )
      stn.(udTxsfld) = ts_op(stn.(hkmtxsfld),stn.(netxsfld),'.*');
      stn.(udTlsfld) = ts_op(stn.(hkmtlsfld),stn.(netlsfld),'.*');
      stn.commentstr = [stn.commentstr,' (TidAdv) '];
    else
      % Exclude tidal currents from advection estimate!
      stn.(udTxsfld) = ts_op(stn.(hkmtxsfld),stn.(qexsfld),'.*');
      stn.(udTlsfld) = ts_op(stn.(hkmtlsfld),stn.(qelsfld),'.*');
    end;
    % Convert to units of [K/hr]
    stn.(udTxsfld).data = -(3600*stn.(udTxsfld).data);
    stn.(udTlsfld).data = -(3600*stn.(udTlsfld).data);
    stn = station_heat_flux_term_inverse(stn,udTfxsfld,udTxsfld,sfld,sal,mhfld);
    stn = station_heat_flux_term_inverse(stn,udTflsfld,udTlsfld,sfld,sal,mhfld);

    stn.opts.add_alongshore_advection = get_opt(stn.opts,'add_alongshore_advection',true);
    if( stn.opts.add_alongshore_advection )
      % Include both vector components of model heat advection
      stn.(udTfld) = ts_op(stn.(udTxsfld),stn.(udTlsfld),'+');
    else
      disp('** Only allowing cross-shore advection **');
      % Ignore along-shore component of heat advection - even 900m gridpoint
      % size at the reef crest may comprise a complex patchwork of deeper
      % (>30m) and shallower (<5m) water we can assume may form a barrier.
      stn.(udTfld) = stn.(udTxsfld);
      stn.commentstr = [stn.commentstr,' (\partial_x_su) '];
    end;

    stn.opts.max_heat_advection = get_opt(stn.opts,'max_heat_advection',Inf);
    maxAdv = stn.opts.max_heat_advection;
    if ( maxAdv < Inf)
      disp(['** Limiting QC''d advection to ',num2str(maxAdv),'K/hr **']);
      %%%%???DEBUG
      stn = qa_ts(stn,udTfld);
      stn.(udTfld).data(stn.(udTfld).data>maxAdv) = maxAdv;
      stn.(udTfld).data(stn.(udTfld).data<-maxAdv) = -maxAdv;
    end;

    stn = station_heat_flux_term_inverse(stn,udTffld,udTfld,sfld,sal,mhfld);

  end;


  %% Surface fluxes


  % This result is NOT used - except to filter dates...
  stn.opts.kd_debug = false;
  stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,stn.opts,asrdiagfld);


  for doWarmix=1:numel(doWarms)
    % Fluxes WITH or WITHOUT warm-layer adjustment
    doWarm = logical(doWarms(doWarmix));

    %tic,

    %% Upward and net longwave flux are now calculated further down using DTCOOL
    % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
    %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
    %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl,sal);
    warning('OFF','Ecoforecasts:Heat:NoRadiative');
    stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                            pfld,sfld,asrfld,[],TURPFX,dsrfld,dlrfld,rfld,...
                            Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl,sal);
                            % Dfld,ssufld,ssvfld,wpfld,whfld,pblzfld,doWarm,max_wl,sal);
                            % % Dfld,ssufld,ssvfld,wpfld,whfld,600,doWarm,max_wl,sal);
                            % % % Dfld,netufld,netvfld,wpfld,whfld,pblzfld,doWarm,max_wl,sal);
    warning('ON','Ecoforecasts:Heat:NoRadiative');

    % Algorithm sometimes returns complex numbers!
    stn.(qlhfld).data = real(stn.(qlhfld).data);
    if ( qlh_adj ~= 1.0 )
      disp(['** Latent heat flux "adjustment factor" ',num2str(qlh_adj)]);
      stn.commentstr = [stn.commentstr,' (Q_L_H\times',num2str(qlh_adj),') '];
      stn.(qlhfld).data = stn.(qlhfld).data.*qlh_adj;
    end;
    stn = station_heat_flux_term(stn,qlhfld,qlhtfld,sfld,sal,mhfld);
    stn.(qshfld).data = real(stn.(qshfld).data);
    stn = station_heat_flux_term(stn,qshfld,qshtfld,sfld,sal,mhfld);
    stn.(qturfld) = ts_op(stn.(qlhfld),stn.(qshfld),'+');
    if ( isfield(stn,qrhfld) && is_valid_ts(stn.(qrhfld)) )
      stn.(qrhfld).data = real(stn.(qrhfld).data);
      stn = station_heat_flux_term(stn,qrhfld,qrhtfld,sfld,sal,mhfld);
      stn.(qturfld) = ts_op(stn.(qturfld),stn.(qrhfld),'+');
    end;
    badix = find(~isfinite(stn.(qturfld).data));
    stn.(qturfld).date(badix) = [];
    stn.(qturfld).data(badix) = [];
    stn = station_heat_flux_term(stn,qturfld,qturtfld,sfld,sal,mhfld);

    % Cross- and long-shore components of the wind stress
    stn = station_reorient_vectors(stn,bathorifld,tauxfld,tauyfld);

    stn.opts.longwave_cool_skin = get_opt(stn.opts,'longwave_cool_skin',true);
    if ( ~stn.opts.longwave_cool_skin )
      disp(['** Longwave upward flux from bulk temperature ',sfld]);
      stn.(scoolfld) = stn.(sfld);
    elseif ( ~isfield(stn,scoolfld) )
      if ( ~isfield(stn,diagfld) || ~isfield(stn.(diagfld),'dtcool') )
        warning('Found no cool-skin temperature adjustment in STN.%s',diagfld);
        stn.(scoolfld) = stn.(sfld);
      else
        [six,diagix] = intersect_dates(stn.(sfld).date,stn.(diagfld).date);
        stn.(scoolfld).date = stn.(diagfld).date(diagix);
        stn.(scoolfld).data = stn.(sfld).data(six) - stn.(diagfld).dtcool(diagix);
      end;
    end;

    if ( ~reanalysis_longwave )
      % % Kraus & Businger (1994) formula
      % stn = station_bulk_ulr(stn,scoolfld,ulrfld,dlrfld);

      % Simple gray-body calculation of ULR
      %stn = station_bulk_ulr(stn,scoolfld,ulrfld);
      stn.opts.epsilon_water = get_opt(stn.opts,'epsilon_water',[]);
      stn = station_bulk_ulr(stn,scoolfld,ulrfld,[],[],stn.opts.epsilon_water);

      %% SENSITIVITY ANALYSIS 
      % stn = station_bulk_ulr(stn,scoolfld,ulrfld,0.95); stn.opts.b_epsw=0.95;
      % stn = station_bulk_ulr(stn,scoolfld,ulrfld,0.98); stn.opts.b_epsw=0.98;
    end;
    % Add down- and upward longwave fluxes together
    stn.(lrfld) = ts_op(stn.(dlrfld),stn.(ulrfld),'-');
    stn = station_heat_flux_term(stn,lrfld,lrtfld,sfld,sal,mhfld);

    %toc,

    % Just in case something above reset it!
    more off;

    disp(['Looping ' num2str(numel(kds)) ' attenuation options']);
    for kdix=1:numel(kds)
      kd = kds{kdix};
      stn.opts.kd = kd;
      stn.opts.kd_debug = true;

      %% SENSITIVITY ANALYSIS 
      % stn.opts.bottom_reflectance = 0.00;
      % stn.opts.bottom_reflectance = 0.17;
      % stn.opts.bottom_reflectance = 0.24;
      % stn.opts.bottom_reflectance = 1.00;

      stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,stn.opts,asrdiagfld);
      stn = station_heat_flux_term(stn,asrfld,asrtfld,sfld,sal,mhfld);


      stn.(qradfld) = ts_op(stn.(asrfld),stn.(lrfld),'+');
      badix = find(~isfinite(stn.(qradfld).data));
      stn.(qradfld).date(badix) = [];
      stn.(qradfld).data(badix) = [];
      stn = station_heat_flux_term(stn,qradfld,qradtfld,sfld,sal,mhfld);
      stn.(qcoolfld) = ts_op(stn.(lrfld),stn.(qturfld),'+');
      stn = station_heat_flux_term(stn,qcoolfld,qcooltfld,sfld,sal,mhfld);

      stn.(q0fld) = ts_op(stn.(qradfld),stn.(qturfld),'+');
      % badix = find(~isfinite(stn.(q0fld).data));
      % stn.(q0fld).date(badix) = [];
      % stn.(q0fld).data(badix) = [];
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,sal,mhfld);

      % Net flux without absorption calculation or benthic flux - for comparison
      stn.(sqradfld) = ts_op(stn.(srfld),stn.(lrfld),'+');
      stn.(sq0fld) = ts_op(stn.(sqradfld),stn.(qturfld),'+');
      stn = station_heat_flux_term(stn,sq0fld,sqtfld,sfld,sal,mhfld);

      disp(['(Re-)Looping ' num2str(numel(cbds)) ' benthic flux options']);
      for cbdix=1:numel(cbds)
        cbd = cbds{cbdix};
        stn.opts.benthic_debug = false;
        stn.opts.b_convective_coefficient = cbd;

        %% Benthic Heat Exchanges
        stn = station_benthic_exchange(stn,sfld,tufld,tvfld,qbfld,btfld,qbofld,stn.opts,qbdiagfld);
        % %%%%??? DEBUG: Use quasi-Eulerian instead of tidal currents
        % stn = station_benthic_exchange(stn,sfld,qeufld,qevfld,qbfld,btfld,qbofld,stn.opts,qbdiagfld);
        % badix = find(~isfinite(stn.(qbofld).data) | abs(stn.(qbofld).data)>2e3);
        badix = find(~isfinite(stn.(qbofld).data));
        stn.(qbofld).date(badix) = [];
        stn.(qbofld).data(badix) = [];
        stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,sal,mhfld);

        stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'+');
        % Try eliminating benthic exchange sub-model        
        if ( ignore_benthos )
          stn.(bq0fld) = stn.(sq0fld);
        end;
        stn = station_heat_flux_term(stn,bq0fld,bq0tfld,sfld,sal,mhfld);
        if ( ~strcmp(bq0lpfld,bq0fld) )
          if ( isfield(stn,bq0lpfld) ); stn = rmfield(stn,bq0lpfld); end;
          stn = verify_variable(stn,bq0lpfld);
        end;

        disp(['(Re-)(Re-)Looping ' num2str(numel(advfacs)) ' advection options']);
        for advix=1:numel(advfacs)
          advfac = advfacs{advix};
          stn.opts.advection_factor = advfac;

          %% Km-scale Heat Advection
          if ( stn.opts.calculate_advection )
            a = build_clim_opt(advfac,'advection_factor',stn.(udTfld).date);
            stn.opts.advection_factor_ts = stn.(udTfld);
            stn.opts.advection_factor_ts.data(:) = a;
            adv.date = stn.(udTfld).date;
            adv.data = a.*stn.(udTfld).data;
            stn.(fqudTfld) = adv;
            stn.(qtAdvfld) = ts_op(stn.(bq0tfld),adv,'+');
            adv=[]; clear adv
          else
            stn.(udTfld) = stn.(bq0tfld);
            stn.(udTfld).data(:) = 0;
            stn.(fqudTfld) = stn.(udTfld);
            stn.(qtAdvfld) = stn.(bq0tfld);
          end; %if ( stn.opts.calculate_advection )

          stn = station_heat_flux_term_inverse(stn,fqudTffld,fqudTfld,sfld,sal,mhfld);
          stn = station_heat_flux_term_inverse(stn,qtAdvffld,qtAdvfld,sfld,sal,mhfld);



          disp(['(Re-)(Re-)(Re-)Looping ' num2str(numel(kths)) ' diffusion options']);
          for kthix=1:numel(kths)
            kth = kths{kthix};
            stn.opts.K_theta = kth;

            %% Km-scale Heat Diffusion
            if ( stn.opts.calculate_diffusion )
              if ( ~isfield(stn,hkmtlfld) )
                disp(['** Recalcing from raw fields: ',rawkd2Tfld,' **']);
                stn = station_calc_kdel2t(stn,stn.opts.K_theta,Tfld,...
                                          rawkd2Tfld,kd2Tfld,...
                                          qtAdvfld,bdTfld,stn.opts.grid_interp_method,false);
              else
                K_theta = build_clim_opt(stn.opts.K_theta,'K_theta',stn.(hkmtlfld).date,false);
                stn.(kd2Tfld).date = stn.(hkmtlfld).date;
                stn.(kd2Tfld).data = K_theta .* stn.(hkmtlfld).data;
                stn.(bdTfld) = ts_op(stn.(qtAdvfld),stn.(kd2Tfld),'+');
              end;
            else
              stn.(kd2Tfld) = stn.(qtAdvfld);
              stn.(kd2Tfld).data(:) = 0;
              stn.(bdTfld) = stn.(qtAdvfld);
            end; %if ( isempty(stn.opts.laplacian_climatology) ) else
            stn = station_heat_flux_term_inverse(stn,kd2Tffld,kd2Tfld,sfld,sal,mhfld);
            stn = station_heat_flux_term_inverse(stn,bdTffld,bdTfld,sfld,sal,mhfld);


            %% Low-pass filtered total heat flux
            if ( ~strcmpi(bdTflpfld,bdTffld) )
              if ( isfield(stn,bdTflpfld) ); stn = rmfield(stn,bdTflpfld); end;
              stn = verify_variable(stn,bdTflpfld);
            end;


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
            disp(['** Horizontal convection fluxes from: ',bdTflpfld,' **']);

            % stn.opts.hc_debug = false;
            stn.opts.hc_debug = true;

            % stn.opts.hc_debug = get_opt(stn.opts,'hc_debug',false);
            % stn.opts.hc_R = get_opt(stn.opts,'hc_R',(1.00-0.08));
            %% SENSITIVITY ANALYSIS 
            % stn.opts.hc_R = (1.00-0.00);
            % stn.opts.hc_R = (1.00-0.11);
            % stn.opts.hc_R = (1.00-0.20);
            stn.opts.hc_scaling = get_opt(stn.opts,'hc_scaling','US');
            %% SENSITIVITY ANALYSIS 
            % stn.opts.hc_scaling = 'UU';
            % stn.opts.hc_scaling = 'SS';
            % stn.opts.hc_max_onset_secs = get_opt(stn.opts,'hc_max_onset_secs',12*3600);
            stn.opts.hc_max_onset_secs = get_opt(stn.opts,'hc_max_onset_secs',24*3600);
            % stn.opts.hc_max_onset_secs = get_opt(stn.opts,'hc_max_onset_secs',36*3600);
            hcres = horizontal_convection(t,s,h,dT,bet,stn.opts,dts,dT,W);
            if ( isempty(strfind(stn.commentstr,' wf=')) )
              stn.commentstr = [stn.commentstr,' ',strrep(stn.opts.hc_scaling,'_','\_'),...
                                ' wf=',num2str(hcres.hc_warming_factor)];
            end;
            
            stn.(hcdTdt).date = dts;
            stn.(hcdTdt).data = hcres.dTdt;
            stn = station_heat_flux_term_inverse(stn,hcdTdtf,hcdTdt,sfld,sal,mhfld);

            % Horizontal convection diagnostics
            stn.(hcu).date = dts;
            stn.(hcu).data = hcres.u;
            % Assuming fit in Fig. 10, panel (c) of Monismith et al (2006)
            stn.([hcu,'_SS']).date = dts;
            stn.([hcu,'_SS']).data = hcres.u_SS;
            % Assuming fit in Fig. 10, panel (a) of Monismith et al (2006)
            stn.([hcu,'_US']).date = dts;
            stn.([hcu,'_US']).data = hcres.u_US;
            % Assuming fit in Fig. 10, panel (f) of Monismith et al (2006)
            stn.([hcu,'_SU']).date = dts;
            stn.([hcu,'_SU']).data = hcres.u_SU;
            % Assuming fit in Fig. 10, panel (d) of Monismith et al (2006)
            stn.([hcu,'_UU']).date = dts;
            stn.([hcu,'_UU']).data = hcres.u_UU;

            stn.(hcdTdx).date = dts;
            stn.(hcdTdx).data = hcres.dTdx;
            stn = station_heat_flux_term_inverse(stn,hcdTdxf,hcdTdx,sfld,sal,mhfld);
            stn.(hcdTdthc).date = dts;
            stn.(hcdTdthc).data = hcres.dTdthc;
            stn = station_heat_flux_term_inverse(stn,hcdTdthcf,hcdTdthc,sfld,sal,mhfld);

            %toc,

            %% Calculate estimation errors for each term in the heat budget
            stn = station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names);
            %DEBUG:
            disp('PRE-QC'); dump_station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names);

            %%%%??? DEBUG - TRY ELIMINATING HORIZONTAL CONVECTION!
            % stn.(hcdTdt) = stn.(bdTfld);
            %%%%??? DEBUG - TRY ELIMINATING HORIZONTAL CONVECTION AND ADVECTION/DIFFUSION!
            % stn.(hcdTdt) = stn.(bq0tfld);
            %%%%??? DEBUG - TRY ELIMINATING EVERYTHING BUT Q0!
            % stn.(hcdTdt) = stn.(qtfld);


            %% 
            % Final QC: remove points with anomalous error est. FROM THIS ESTIMATE

            baddts = [];
            if ( ~exist('qlh_err_cutoff','var') || isempty(qlh_err_cutoff) )
              qlh_err_cutoff = 1000;
            end;
            stn = verify_variable(stn,qlhfld1derr);
            baddts = [baddts;...
                      stn.(qlhfld1derr).date(abs(stn.(qlhfld1derr).data)>qlh_err_cutoff)];
            if ( ~exist('bq0_err_cutoff','var') || isempty(bq0_err_cutoff) )
              bq0_err_cutoff = 1500;
            end;
            stn = verify_variable(stn,bq0fld1derr);
            baddts = [baddts;...
                      stn.(bq0fld1derr).date(abs(stn.(bq0fld1derr).data)>bq0_err_cutoff)];
            if ( ~exist('hc_err_cutoff','var') || isempty(hc_err_cutoff) )
              hc_err_cutoff = 13.00;
            end;
            stn = verify_variable(stn,hcdTdt1derr);
            baddts = [baddts;...
                      stn.(hcdTdt1derr).date(abs(stn.(hcdTdt1derr).data)>hc_err_cutoff)];

            baddays = unique(floor(baddts));
            disp(['** Days with anomalous errors removed: ',num2str(numel(baddays)),' **']);
            %DEBUG:            disp('** NOT REMOVED **'); if (0);
            if ( ~isempty(baddays) )

              %DEBUG:
              find_date_ranges(baddays,1);

              %DEBUG:
              disp(['** Budget days remaining: ',num2str(numel(stn.(hcdTdt).data)/24),' **']);
              %DEBUG:
              prctile(stn.(hcdTdt1derr).data,[99,99.9,99.99,99.999,99.9999,]),

              filtgood = @(x)(find(~ismember(floor(x.date),baddays)));

              stn.(srfld) = subset_ts(stn.(srfld),filtgood);
              stn.(asrfld) = subset_ts(stn.(asrfld),filtgood);
              stn.(lrfld) = subset_ts(stn.(lrfld),filtgood);
              stn.(qlhfld) = subset_ts(stn.(qlhfld),filtgood);
              stn.(qshfld) = subset_ts(stn.(qshfld),filtgood);

              stn.(q0fld) = subset_ts(stn.(q0fld),filtgood);
              stn.(qtfld) = subset_ts(stn.(qtfld),filtgood);
              stn.(sq0fld) = subset_ts(stn.(sq0fld),filtgood);
              stn.(sqtfld) = subset_ts(stn.(sqtfld),filtgood);
              stn.(bq0fld) = subset_ts(stn.(bq0fld),filtgood);
              stn.(bq0tfld) = subset_ts(stn.(bq0tfld),filtgood);

              stn.(qtAdvfld) = subset_ts(stn.(qtAdvfld),filtgood);
              stn.(qtAdvffld) = subset_ts(stn.(qtAdvffld),filtgood);
              stn.(fqudTfld) = subset_ts(stn.(fqudTfld),filtgood);
              stn.(fqudTffld) = subset_ts(stn.(fqudTffld),filtgood);
              stn.(kd2Tfld) = subset_ts(stn.(kd2Tfld),filtgood);
              stn.(kd2Tffld) = subset_ts(stn.(kd2Tffld),filtgood);

              stn.(bdTfld) = subset_ts(stn.(bdTfld),filtgood);
              stn.(bdTffld) = subset_ts(stn.(bdTffld),filtgood);
              stn.(hcdTdt) = subset_ts(stn.(hcdTdt),filtgood);
              stn.(hcdTdtf) = subset_ts(stn.(hcdTdtf),filtgood);


              stn.([srfld,'_err']) = subset_ts(stn.([srfld,'_err']),filtgood);
              stn.([asrfld,'_err']) = subset_ts(stn.([asrfld,'_err']),filtgood);
              stn.([lrfld,'_err']) = subset_ts(stn.([lrfld,'_err']),filtgood);
              stn.([qlhfld,'_err']) = subset_ts(stn.([qlhfld,'_err']),filtgood);
              stn.([qshfld,'_err']) = subset_ts(stn.([qshfld,'_err']),filtgood);

              stn.([q0fld,'_err']) = subset_ts(stn.([q0fld,'_err']),filtgood);
              stn.([sq0fld,'_err']) = subset_ts(stn.([sq0fld,'_err']),filtgood);
              stn.([bq0fld,'_err']) = subset_ts(stn.([bq0fld,'_err']),filtgood);

              stn.([fqudTfld,'_err']) = subset_ts(stn.([fqudTfld,'_err']),filtgood);
              stn.([fqudTffld,'_err']) = subset_ts(stn.([fqudTffld,'_err']),filtgood);
              stn.([kd2Tfld,'_err']) = subset_ts(stn.([kd2Tfld,'_err']),filtgood);
              stn.([kd2Tffld,'_err']) = subset_ts(stn.([kd2Tffld,'_err']),filtgood);
              stn.([qtAdvfld,'_err']) = subset_ts(stn.([qtAdvfld,'_err']),filtgood);
              stn.([qtAdvffld,'_err']) = subset_ts(stn.([qtAdvffld,'_err']),filtgood);

              stn.([hcdTdt,'_err']) = subset_ts(stn.([hcdTdt,'_err']),filtgood);

            end; %if ( ~isempty(baddays) )


            %DEBUG:
            disp(numel(stn.(hcdTdt).data));

            % Recalculate global errors and covariances after outlier removal
            stn = station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names);


            if ( doPlot )


              %% Calculate and plot time series and errors from this option
              clear t sqt dt q qe
              [t,sqt,bqt,dt,q,qe] = ...
                  intersect_tses(stn.(sfld),stn.(sqtfld),stn.(bq0tfld),stn.(bdTfld),stn.(hcdTdt),stn.([hcdTdt,'_err']));


              %% Limit accumulation period so different methods intercompare
              if ( ~exist('begyr','var' ) )
                begyr = get_year(t.date(1));
                %%%%??? DEBUG
                % % USF AVHRR 1km SST data (for Florida) before 1996 appears highly suspect
                % if ( ~isempty(regexp(KMPFX,'avhrr')) ); begyr = 1996; end;
                % if ( ~isempty(regexp(sfld,'avhrr')) ); begyr = 1996; end;
                % % Make periods for WAVEPFX='ww3','erai','ndbc' all match
                % begyr = 2000;
                % % Make periods for KMPFX='none','avhrr_weekly' both match
                % begyr = 1996;
                % % Consistent with LONF1, DRYF1 (latest) T_s deployments
                % begyr = 1993;
              end;

              if ( ~exist('endyr','var' ) )
                % Some datasets end a few hours before the New Year, so add 1 day
                endyr = get_year(t.date(end)+1);
                %%%%??? DEBUG
                % endyr = 2011;
                % % Consistent with SANF1, DRYF1 T_s failures
                % endyr = 2005;
                % % Consistent with SMKF1 T_s failure
                % endyr = 2007;
              end;

              begdt = datenum(begyr,1,2);
              enddt = datenum(endyr,1,1);
              t.data(begdt>t.date|t.date>enddt)=[];
              t.date(begdt>t.date|t.date>enddt)=[];
              sqt.data(begdt>sqt.date|sqt.date>enddt)=[];
              sqt.date(begdt>sqt.date|sqt.date>enddt)=[];
              bqt.data(begdt>bqt.date|bqt.date>enddt)=[];
              bqt.date(begdt>bqt.date|bqt.date>enddt)=[];
              dt.data(begdt>dt.date|dt.date>enddt)=[];
              dt.date(begdt>dt.date|dt.date>enddt)=[];
              q.data(begdt>q.date|q.date>enddt)=[];
              q.date(begdt>q.date|q.date>enddt)=[];
              qe.data(begdt>qe.date|qe.date>enddt)=[];
              qe.date(begdt>qe.date|qe.date>enddt)=[];
              begyr = min(get_year(t.date));
              endyr = max(get_year(t.date));


              % Evaluate total error for this option
              t0 = t.data(1);
              sq.date = q.date;
              sq.data = t0 + nancumsum(q.data) - q.data(1);
              stn.optim.error(kdix,doWarmix,advix,kthix) = sqrt(sum((t.data-sq.data).^2));

              ssqt.date = sqt.date;
              ssqt.data = t0 + nancumsum(sqt.data) - sqt.data(1);
              sbqt.date = bqt.date;
              sbqt.data = t0 + nancumsum(bqt.data) - bqt.data(1);
              sdt.date = dt.date;
              sdt.data = t0 + nancumsum(dt.data) - dt.data(1);

              stn.optim.doWarm(kdix,doWarmix,advix,kthix) = doWarm;
              stn.optim.kd{kdix,doWarmix,advix,kthix} = kd;
              stn.optim.cbd(kdix,doWarmix,advix,kthix) = cbd;
              stn.optim.kth{kdix,doWarmix,advix,kthix} = kth;
              stn.optim.q(kdix,doWarmix,advix,kthix) = q;
              stn.optim.sq(kdix,doWarmix,advix,kthix) = sq;

              % Evaluate daily climatological error for this option
              if ( ~isfield(stn.optim,'climt') )
                [cum,tid] = grp_ts(t.data,t.date,'daily',@nanmean,23);
                stn.optim.climt.date = tid;
                stn.optim.climt.data = cum;
              end;
              t0 = stn.optim.climt.data(1);

              % [cum,tid] = grp_ts(q.data,q.date,'daily',@nansum,23);
              [cum,tid] = grp_ts(q.data,q.date,'daily',@nanmean,23);
              cum = 24*cum;
              stn.optim.climq(kdix,doWarmix,advix,kthix).date = tid;
              stn.optim.climq(kdix,doWarmix,advix,kthix).data = cum;
              sq.date = tid;
              sq.data = t0 + nancumsum(cum) - cum(1);
              stn.optim.climsq(kdix,doWarmix,advix,kthix).date = sq.date;
              stn.optim.climsq(kdix,doWarmix,advix,kthix).data = sq.data;
              stn.optim.climerror(kdix,doWarmix,advix,kthix) = ...
                  sqrt(sum((stn.optim.climt.data-sq.data).^2)/numel(sq.data));
              %[cum,tid] = grp_ts(qe.data,qe.date,'daily',@(x)(nanmean(abs(x))),23);
              [cum,tid] = grp_ts(qe.data,qe.date,'daily',@(x)(nanmax(abs(x))),23);
              cum = 24*cum;
              % [cum,tid] = grp_ts(qe.data,qe.date,'daily',@(x)(nansum(abs(x))),23);
              stn.optim.climsq_err(kdix,doWarmix,advix,kthix).date = tid;
              stn.optim.climsq_err(kdix,doWarmix,advix,kthix).data = cum;

              [cum,tid] = grp_ts(sqt.data,sqt.date,'daily',@nanmean,23);
              cum = 24*cum;
              stn.optim.climsqt(kdix,doWarmix,advix,kthix).date = tid;
              stn.optim.climsqt(kdix,doWarmix,advix,kthix).data = cum;
              ssqt.date = tid;
              ssqt.data = t0 + nancumsum(cum) - cum(1);
              stn.optim.climssqt(kdix,doWarmix,advix,kthix).date = ssqt.date;
              stn.optim.climssqt(kdix,doWarmix,advix,kthix).data = ssqt.data;

              [cum,tid] = grp_ts(bqt.data,bqt.date,'daily',@nanmean,23);
              cum = 24*cum;
              stn.optim.climbqt(kdix,doWarmix,advix,kthix).date = tid;
              stn.optim.climbqt(kdix,doWarmix,advix,kthix).data = cum;
              sbqt.date = tid;
              sbqt.data = t0 + nancumsum(cum) - cum(1);
              stn.optim.climsbqt(kdix,doWarmix,advix,kthix).date = sbqt.date;
              stn.optim.climsbqt(kdix,doWarmix,advix,kthix).data = sbqt.data;

              [cum,tid] = grp_ts(dt.data,dt.date,'daily',@nanmean,23);
              cum = 24*cum;
              stn.optim.climdt(kdix,doWarmix,advix,kthix).date = tid;
              stn.optim.climdt(kdix,doWarmix,advix,kthix).data = cum;
              sdt.date = tid;
              sdt.data = t0 + nancumsum(cum) - cum(1);
              stn.optim.climsdt(kdix,doWarmix,advix,kthix).date = sdt.date;
              stn.optim.climsdt(kdix,doWarmix,advix,kthix).data = sdt.data;


              % Evaluate error vs. once-daily mean sea temperature for this option
              [cum,tid] = grp_ts(t.data,t.date,@floor,@nanmean,24);
              stn.optim.dayt.date = tid;
              stn.optim.dayt.data = cum; % We refilter above, so recalc every time
              stn.optim.daydt.date = tid(2:end);
              stn.optim.daydt.data = diff(cum); % We refilter above, so recalc every time
              stn.optim.daydt.data(diff(tid) > 1.1)=[]; stn.optim.daydt.date(diff(tid) > 1.1)=[];
              [cum,tid] = grp_ts(q.data,q.date,@floor,@nansum,24);
              stn.optim.dayq(kdix,doWarmix,advix,kthix).date = tid;
              stn.optim.dayq(kdix,doWarmix,advix,kthix).data = cum;
              t0 = stn.optim.dayt.data(1);
              daysq.date = tid;
              daysq.data = t0 + nancumsum(cum) - cum(1);
              stn.optim.daysq(kdix,doWarmix,advix,kthix).date = daysq.date;
              stn.optim.daysq(kdix,doWarmix,advix,kthix).data = daysq.data;
              % Daily mean temperature vs. accumulated fluxes - WHOLE RECORD
              [tix,qix] = intersect_dates(stn.optim.dayt.date,daysq.date);
              stn.optim.dayerror(kdix,doWarmix,advix,kthix) = ...
                  sqrt(sum((stn.optim.dayt.data(tix)-daysq.data(qix)).^2)/numel(daysq.data(qix)));
              % Change in daily mean temperature vs. daily sum of fluxes
              [tix,qix] = intersect_dates(stn.optim.daydt.date,stn.optim.dayq.date);
              stn.optim.dayrmse(kdix,doWarmix,advix,kthix) = ...
                  sqrt(sum((stn.optim.daydt.data(tix)-stn.optim.dayq.data(qix)).^2)/numel(qix));
              % [cum,tid] = grp_ts(qe.data,qe.date,@floor,@nansum,24);
              [cum,tid] = grp_ts(qe.data,qe.date,@floor,@(x)(nansum(abs(x))),24);
              stn.optim.daysq_err(kdix,doWarmix,advix,kthix).date = tid;
              stn.optim.daysq_err(kdix,doWarmix,advix,kthix).data = cum;

              % Evaluate seasonal amplitude error for this option
              stn.optim.minseassq = +Inf;
              stn.optim.maxseassq = -Inf;
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
                  stn.optim.seasyear(kdix,doWarmix,advix,kthix,yrix) = yr;
                  stn.optim.seassq(kdix,doWarmix,advix,kthix,yrix) = sq;
                  stn.optim.minseassq = nanmin(stn.optim.minseassq,nanmin(sq.data(:)));
                  stn.optim.maxseassq = nanmax(stn.optim.maxseassq,nanmax(sq.data(:)));
                  stn.optim.seaserror(kdix,doWarmix,advix,kthix,yrix) = sqrt(sum((t.data(ix)-sq.data).^2));
                end; %if isempty(ix) else
              end; %for yrix

            end; %if doPlot

          end; %for kthix=1:numel(kths)

        end; %for advix=1:numel(advfacs)

      end; %for cbdix=1:numel(cbds)

    end; %for kdix=1:numel(kds)

  end; %for doWarmix=1:numel(doWarms)


  %%%
  %% Calculate Sub-Grid Scale Heat Diffusion (as a residual - of LAST estimate!)
  if ( stn.opts.calculate_diffusion )
    if ( ~isfield(stn,Tfld) || ~isfield(stn,hkmtlfld) )
      warning('No sea-surface temperature field for SGS diffusion ("%s")',Tfld);
    else
      stn = station_calc_sgs_diffusion(stn,sfld,hcdTdt,Tfld,hkmtlfld,...
                                       sgskd2Tfld,sgskfld,sgsdTdt);
    end;
  end;


  % Report estimation errors for each heat budget term after outlier removal
  %DEBUG:
  disp('POST-QC'); dump_station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names);


  % Create fields for comparison of daily averages, sums, and changes
  stn.(dly_sfld).date=[]; stn.(dly_sfld).data=[];
  [stn.(dly_sfld).data,stn.(dly_sfld).date] = grp_ts(stn.(sfld).data,stn.(sfld).date,@floor,@nanmean,23);
  stn.(dly_dsfld).date = stn.(dly_sfld).date(2:end);
  stn.(dly_dsfld).data = diff(stn.(dly_sfld).data);
  stn = filter_gaps(stn,dly_sfld,dly_dsfld,1.5);
  stn = station_heat_flux_term_inverse(stn,dly_dsffld,dly_dsfld,sfld,sal,mhfld);
  stn.(dly_srfld) = oshb_daily_sum(stn.(srfld));
  stn.(dly_srtfld) = oshb_daily_sum(stn.(srtfld));
  stn.(dly_srfld) = oshb_daily_sum(stn.(srfld));
  stn.(dly_asrtfld) = oshb_daily_sum(stn.(asrtfld));
  stn.(dly_lrfld) = oshb_daily_sum(stn.(lrfld));
  stn.(dly_lrtfld) = oshb_daily_sum(stn.(lrtfld));
  stn.(dly_qlhfld) = oshb_daily_sum(stn.(qlhfld));
  stn.(dly_qlhtfld) = oshb_daily_sum(stn.(qlhtfld));
  stn.(dly_qshfld) = oshb_daily_sum(stn.(qshfld));
  stn.(dly_qshtfld) = oshb_daily_sum(stn.(qshtfld));
  stn.(dly_q0fld) = oshb_daily_sum(stn.(q0fld));
  stn.(dly_qtfld) = oshb_daily_sum(stn.(qtfld));
  stn.(dly_sq0fld) = oshb_daily_sum(stn.(sq0fld));
  stn.(dly_sqtfld) = oshb_daily_sum(stn.(sqtfld));
  stn.(dly_bq0fld) = oshb_daily_sum(stn.(bq0fld));
  stn.(dly_bq0tfld) = oshb_daily_sum(stn.(bq0tfld));
  stn.(dly_bdTfld) = oshb_daily_sum(stn.(bdTfld));
  stn.(dly_bdTffld) = oshb_daily_sum(stn.(bdTffld));
  stn.(dly_hcdTdthc) = oshb_daily_sum(stn.(hcdTdthc));
  stn.(dly_hcdTdthcf) = oshb_daily_sum(stn.(hcdTdthcf));
  stn.(dly_hcdTdt) = oshb_daily_sum(stn.(hcdTdt));
  stn.(dly_hcdTdtf) = oshb_daily_sum(stn.(hcdTdtf));

  % Create running one-day averages, e.g., to compare with climatologies
  stn = verify_variable(stn,s1dfld);
  stn = verify_variable(stn,a1dfld);
  stn = verify_variable(stn,sq01dfld);
  stn = verify_variable(stn,sr1dfld);
  stn = verify_variable(stn,asr1dfld);
  stn = verify_variable(stn,lr1dfld);
  stn = verify_variable(stn,qlh1dfld);
  stn = verify_variable(stn,qsh1dfld);
  stn = verify_variable(stn,qcool1dfld);
  stn = verify_variable(stn,fqudTffld1d);
  stn = verify_variable(stn,qbofld1d);
  stn = verify_variable(stn,kd2Tffld1d);
  stn = verify_variable(stn,hcdTdtf1d);
  stn = verify_variable(stn,hcdTdthcf1d);



  if ( doPlot )

    %{
    for doWarmix=1:numel(doWarms)
      doWarm = logical(doWarms(doWarmix));
      if ( doWarm );	doWarmStr = 'Warm Layer ';
      else;		doWarmStr = 'No Warm Layer ';
      end;

      for advix=1:numel(advfacs)
        advfac = advfacs{advix};
        advStr = [' Adv=',stn.optim.advfacstrs{advix},' '];

        fmg; plot(squeeze(stn.optim.error(:,doWarmix,advix,:))); titlename([STNM ' Total Error: ' doWarmStr advStr stn.commentstr]);
        legend(stn.optim.kthstrs);
        set(gca,'XTick',[1:numel(kds)],'XTickLabel',stn.optim.kdstrs);
        ylim([nanmin(stn.optim.error(:)),nanmax(stn.optim.error(:))]);

        fmg; plot(squeeze(nanmedian(stn.optim.seaserror(:,doWarmix,advix,:),4))); titlename([STNM ' Mdn Seas Err: ' doWarmStr advStr stn.commentstr]);
        legend(stn.optim.kthstrs);
        set(gca,'XTick',[1:numel(kds)],'XTickLabel',stn.optim.kdstrs);
        ylim([nanmin(stn.optim.seaserror(:)),nanmax(stn.optim.seaserror(:))]);
      end;
    end;
    %}


    % for cbdix=1:numel(cbds)
    %   cbd = cbds{cbdix};
    for advix=1:numel(advfacs)
      advfac = advfacs{advix};
      for kthix=1:numel(kths)
        kth = kths{kthix};

        for doWarmix=1:numel(doWarms)
          %doWarm = logical(doWarms(doWarmix));
          doWarm = stn.optim.doWarm(1,doWarmix,1,1);
          if ( doWarm ); doWarmStr = 'Warm Layer';
          else;	         doWarmStr = 'No Warm Layer';
          end;


          plot_daily_clim_heat_budget(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,doWarmix,advix,kthix);
          if ( doWarm )
            appendtitlename(' (WARM ');
          else
            appendtitlename(' (NO WARM ');
          end;
          appendtitlename([strrep(hcdTdt,'_','\_'),' Adv=',stn.optim.advfacstrs{advix},' K_\theta=',stn.optim.kthstrs{end},') ']);

          %{
          fmg;
          climsq = squeeze(stn.optim.climsq(:,doWarmix,advix,kthix));
          climsq_err = stn.optim.climsq_err(end,doWarmix,advix,kthix);
          climsq_minus_err = ts_op(climsq(end),climsq_err,'-');
          climsq_plus_err = ts_op(climsq(end),climsq_err,'+');
          % climsq_minus_err = ts_op(climsq(end),ts_op(climsq_err,24,'*'),'-');
          % climsq_plus_err = ts_op(climsq(end),ts_op(climsq_err,24,'*'),'+');
          % % climsq_minus_err.date = climsq_err.date;
          % % % climsq_minus_err.data = climsq(end).data - cumsum(climsq_err.data);
          % % climsq_minus_err.data = climsq(end).data - (climsq_err.data.*24);
          % % climsq_plus_err.date = climsq_err.date;
          % % % climsq_plus_err.data = climsq(end).data + cumsum(climsq_err.data);
          % % % plot_ts(stn.optim.climt,climsq,climsq_minus_err,'k^',climsq_plus_err,'kv');
          % % climsq_plus_err.data = climsq(end).data + (climsq_err.data.*24);
          lhs=plot_ts(stn.optim.climt,'k','LineWidth',3,climsq,'Color',[.5,.5,.5],'LineWidth',1.5,climsq_minus_err,'k:','LineWidth',1.5,climsq_plus_err,'k:','LineWidth',1.5);
          lhs(end) = [];

          % climssqt = squeeze(stn.optim.climssqt(end,doWarmix,advix,kthix));
          % climsbqt = squeeze(stn.optim.climsbqt(end,doWarmix,advix,kthix));
          % climsdt = squeeze(stn.optim.climsdt(end,doWarmix,advix,kthix));
          % lhs(end+1:end+3) = plot_ts(climssqt,'r--',climsbqt,'co',climsdt,'m-.');

          datetick3('x',3);
          titlename([STNM ' Daily Clim: ' doWarmStr ' ' stn.optim.cbdstrs{default_cbdix}]);
          % % legend({'T_s',stn.optim.kdstrs{:},[stn.optim.kdstrs{end},' \pm error']});
          kdstrs = strcat( stn.optim.kdstrs',{' Errs '},...
                           cellstr(num2str(stn.optim.dayrmse(:,doWarmix,advix,kthix),'%.1f')),{', '},...
                           cellstr(num2str(stn.optim.dayerror(:,doWarmix,advix,kthix),'%.1f')),{', '},...
                           cellstr(num2str(stn.optim.climerror(:,doWarmix,advix,kthix),'%.1f')) );
          %%%%DEBUG:
          for kdstrix=1:numel(kdstrs); disp(kdstrs{kdstrix}); end;
          % legend(lhs,{'T_s',kdstrs{:},[stn.optim.kdstrs{end},' \pm error']}, 'Location','South');
          legend(lhs,...
                 {'T_s',kdstrs{:},[stn.optim.kdstrs{end},' \pm error'],...
                  % 'Q_0/\rhoC_ph','(Q_0(\gamma)+Q_b)/\rhoC_ph','Non-HC',...
                 }, 'Location','South');
          xlim(stn.optim.climt.date([1 end]));
          ylim([minmin([stn.optim.climsq.data]),maxmax([stn.optim.climsq.data])]);
          ylim([16,34]);
          % ylim([16,50]);
          ylabel('^oC');
          % appendtitlename([' (' strrep(QEPFX,'_','\_') ' Adv=' stn.optim.advfacstrs{advix} ' K_\theta=' stn.optim.kthstrs{kthix} ')' stn.commentstr]);
          appendtitlename([' (' strrep(hcdTdt,'_','\_') ' Adv=' stn.optim.advfacstrs{advix} ' K_\theta=' stn.optim.kthstrs{kthix} ')' stn.commentstr]);
          appendtitlename([' (' num2str(begyr) '-' num2str(endyr) ')']);
          %}

          % Print out the daily climatology plot for the LAST set of options
          if ( advix==numel(advfacs) && kthix==numel(kths) && doWarmix==numel(doWarms) )
            if ( doPrint )
              print('-dtiff',fullfile(figspath,[stnm,oldstr,'-daily-clim-',hcdTdt,'.tiff']));
            end;
          end;

        end; %for doWarmix=1:numel(doWarms)

        %{
        for doWarmix=1:numel(doWarms)
          doWarm = logical(doWarms(doWarmix));
          if ( doWarm ); doWarmStr = 'Warm Layer';
          else;	         doWarmStr = 'No Warm Layer';
          end;

          fmg; plot_ts(stn.(sfld),squeeze(stn.optim.sq(:,doWarmix,advix,kthix))); titlename([STNM ' Time Series: ' doWarmStr stn.commentstr]);
          legend({'T_s',stn.optim.kdstrs{:}});
          ylim([minmin([stn.optim.sq.data]),maxmax([stn.optim.sq.data])]);
        end;

        yrs = unique(get_year(stn.optim.sq(1,1,1,kthix).date));
        nrows = floor(sqrt(numel(yrs)));
        ncols = ceil(numel(yrs)/nrows);

        for doWarmix=1:numel(doWarms)
          doWarm = logical(doWarms(doWarmix));
          fmg;
          if (~doWarm)	suptitle([STNM ' Annual TS: No Warm Layer' stn.commentstr]);
          else		suptitle([STNM ' Annual TS: Warm Layer' stn.commentstr]);
          end;
          for yrix=1:numel(yrs)
            yr = yrs(yrix);
            if ( is_valid_ts(stn.optim.seassq(1,doWarmix,advix,kthix,yrix)) )
              ix = find(get_year(stn.(sfld).date)==yr);
              t.date = stn.(sfld).date(ix);
              t.data = stn.(sfld).data(ix);

              subplot_tight(nrows,ncols,yrix);
              plot_ts(t,stn.optim.seassq(:,doWarmix,advix,kthix,yrix));
              datetick('x',17,'keeplimits');
              % legend({'T_s',stn.optim.kdstrs{:}});
              xlabel(num2str(yr));
              xlim([min(stn.optim.seassq(kdix,doWarmix,advix,kthix,yrix).date),...
                    max(stn.optim.seassq(kdix,doWarmix,advix,kthix,yrix).date)]);
              ylim([stn.optim.minseassq,stn.optim.maxseassq]);
              % ylim([minmin([stn.optim.seassq.data]),maxmax([stn.optim.seassq.data])]);
            end; %if is_valid_ts
          end; %for yrix
        end; %for doWarmix=1:numel(doWarms)
        1;
        %}

      end; %for kthix
    end; %for advix
    % end; %for cbdix


    % Append comments from *last* option tried, to STN.COMMENTSTR
    if ( doWarm )
      stn.commentstr = [stn.commentstr,' (WARM '];
    else
      stn.commentstr = [stn.commentstr,' (NO WARM '];
    end;
    %stn.commentstr = [stn.commentstr,strrep(QEPFX,'_','\_'),' Adv=',stn.optim.advfacstrs{advix},' K_\theta=',stn.optim.kthstrs{end},') '];
    stn.commentstr = [stn.commentstr,strrep(hcdTdt,'_','\_'),' Adv=',stn.optim.advfacstrs{advix},' K_\theta=',stn.optim.kthstrs{end},') '];

    if ( doScatterplots )
      % scatter_fit_ts(stn.(dly_bq0tfld),stn.optim.daydt,[],[],'(Q_0+Q_b)/\rhoC_ph','\Delta_1_dT_s',[],[],true);
      % axis([-3.5,3.5,-3.5,3.5]);
      % scatter_fit_ts(stn.(dly_bdTfld),stn.optim.daydt,[],[],'u^.\nabla_hT_s+K_H\nabla_h^2T_s+(Q_0+Q_b)/\rhoC_ph','\Delta_1_dT_s',[],[],true);
      % axis([-3.5,3.5,-3.5,3.5]);
      % scatter_fit_ts(stn.optim.dayq,stn.optim.daydt,[],[],'\partial_tT_s','\Delta_1_dT_s',[],[],true);
      % axis([-3.5,3.5,-3.5,3.5]);
      scatter_fit_ts_seasons(stn.(dly_bq0tfld),stn.optim.daydt,[],[],'\Sigma_1_d(Q_0+Q_b)/\rhoC_ph','\Delta_1_dT_s',[],[],true);
      axis([-3.5,3.5,-3.5,3.5]);
      scatter_fit_ts_seasons(stn.optim.dayq,stn.optim.daydt,[],[],'\Sigma_1_d\partial_tT_s','\Delta_1_dT_s',[],[],true);
      axis([-3.5,3.5,-3.5,3.5]);
    end;

    if ( doDiagplots || doAnnSubs )
      % Show year-by-year comparison SUBPLOTS for whichever configuration we ran last
      annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,[],[],begyr,[]);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-annsubs-',hcdTdt,'.tiff']));
      end;

      % % Show year-by-year comparison TIME SERIES for whichever configuration we ran last
      % plot_budget_years(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,[],[],begyr,[]);
      % ylim([-50,150]);
      % if ( doPrint )
      %   print('-dtiff',fullfile(figspath,[stnm,oldstr,'-annsubs_ts-',hcdTdt,'.tiff']));
      % end;
    end;


    % Compare published climatologies to our (last) estimate
    if ( (doClimplots || doBoxplots) && ~isfield(stn,climq0fld) )
      switch ( CLIMPFX )
       case 'daily_oaflux',
        stn = station_load_oaflux(stn);
       otherwise,
        error('Do not know how to load climatology "%s"',CLIMPFX);
      end;
    end;


    if ( doClimplots )


      %% Climatological Comparison figures


      % %%%% ??? DEBUG
      % ms1_clim(stn,'daily',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,[],[],[],[],begyr);
      % figfname= fullfile(figspath,[stnm,oldstr,'_',hcdTdt,'_weekly.']);
      % if ( doPrint )
      %   print('-dtiff',[figfname 'tiff']);
      %   print('-dpng',[figfname 'png']);
      % end;

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
      % %appendtitlename(strrep(sprintf(' (%s,%s) %s', upper(RAPFX), upper(KMPFX), stn.commentstr),'_','\_'));
      % appendtitlename(stn.commentstr);
      % %%%% ??? DEBUG

      % Compare simple cumulative sums of daily climatology vs. our (last) estimate
      compare_heat_budget_cumsums(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-simple-sum-',hcdTdtf,'.tiff']));
      end;


      % Compare monthly climatologies vs. flux implied by actual dTs vs. our (last) estimate
      stn = compare_monthly_flux_climatologies(stn,sfld,sq0fld,doPrint,[],oldstr);


      % Plot daily climatology comparisons for whichever configuration we ran last
      compare_flux_climatologies(stn,'daily',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-chkann-',hcdTdt,'.tiff']));
      end;
      compare_flux_climatologies(stn,'yearly',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-chkann-interann-',hcdTdt,'.tiff']));
      end;

    end; %if ( doClimplots )


    if ( doBoxplots )

      % Compare one-day averages (sub-sampled) with OAFlux climatology

      [cix,six] = intersect_dates(stn.(climq0fld).date,stn.(sq01dfld).date);

      fmg;
      boxplot_ts(stn.(climq0fld),'month','mean',true,'index',cix,...
                 'title',[STNM,' OAFlux/ISCCP Q_0',stn.commentstr]);
      ylim([-1000,1000]); ylabel('W/m^2');

      fmg;
      boxplot_ts(stn.(sq01dfld),'month','mean',true,'index',six,...
                 'title',[STNM,' Gramer&Mariano Sea-surface Q_0',stn.commentstr]);
      ylim([-1000,600]); ylabel('W/m^2');
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',sq01dfld,'-closeup.tiff']));
      end;
      ylim([-1000,1000]);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',sq01dfld,'.tiff']));
      end;


      [cix,six] = intersect_dates(stn.(climsrfld).date,stn.(sr1dfld).date);

      fmg;
      boxplot_ts(stn.(climsrfld),'month','mean',true,'index',cix,...
                 'title',[STNM,' OAFlux/ISCCP Q_S_W',stn.commentstr]);
      ylim([-1000,1000]); ylabel('W/m^2');
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',climsrfld,'.tiff']));
      end;

      fmg;
      boxplot_ts(stn.(sr1dfld),'month','mean',true,'index',six,...
                 'title',[STNM,' Gramer&Mariano Q_S_W',stn.commentstr]);
      ylim([-1000,1000]); ylabel('W/m^2');
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',sr1dfld,'.tiff']));
      end;


      [cix,six] = intersect_dates(stn.(climqlhfld).date,stn.(qlh1dfld).date);

      fmg;
      boxplot_ts(stn.(climqlhfld),'month','mean',true,'index',cix,...
                 'title',[STNM,' OAFlux Q_L_H',stn.commentstr]);
      ylim([-1000,1000]); ylabel('W/m^2');

      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',climqlhfld,'.tiff']));
      end;

      fmg;
      boxplot_ts(stn.(qlh1dfld),'month','mean',true,'index',six,...
                 'title',[STNM,' Gramer&Mariano Q_L_H',stn.commentstr]);
      ylim([-1000,1000]); ylabel('W/m^2');
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',qlh1dfld,'.tiff']));
      end;


      fmg;
      sh=boxplot_ts(stn.(sr1dfld),'month','mean',true,'allcolors','b');
      lh=boxplot_ts(stn.(qcool1dfld),'month','mean',true,'allcolors','r');
      ylim([-1000,1000]); ylabel('W/m^2');
      legend([sh(1),lh(1)], 'Q_S_W','Q_L_W+Q_L_H+Q_S_H', 'Location','South');
      titlename([STNM,' Gramer&Mariano Air-Sea fluxes (1d avg)',stn.commentstr]);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',sr1dfld,'-vs-',qcool1dfld,'.tiff']));
      end;

      fmg;
      boxplot_ts(stn.(fqudTffld1d),'month','mean',true,'allcolors','k',...
                 'title',[STNM,' Gramer&Mariano Km-scale Heat advection (1d avg)',stn.commentstr]);
      ylim([-1000,1000]); ylabel('W/m^2');
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',fqudTffld1d,'.tiff']));
      end;

      fmg;
      boxplot_ts(stn.(qbofld1d),'month','mean',true,'allcolors','k',...
                 'title',[STNM,' Gramer&Mariano Benthic Heat flux (1d avg)',stn.commentstr]);
      ylim([-1000,1000]); ylabel('W/m^2');
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',qbofld1d,'.tiff']));
      end;

      fmg;
      boxplot_ts(stn.(kd2Tffld1d),'month','mean',true,'allcolors','k',...
                 'title',[STNM,' Gramer&Mariano Km-scale Heat diffusion (1d avg)',stn.commentstr]);
      ylim([-1000,1000]); ylabel('W/m^2');
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',kd2Tffld1d,'.tiff']));
      end;

      fmg;
      boxplot_ts(stn.(hcdTdthcf1d),'month','mean',true,'allcolors','k',...
                 'title',[STNM,' Gramer&Mariano Horizontal Convection (1d avg)',stn.commentstr]);
      ylim([-1000,1000]); ylabel('W/m^2');
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',hcdTdthcf1d,'.tiff']));
      end;
      %end;

      fmg;
      lh=boxplot_ts(stn.(sr1dfld),'month','mean',true,'allcol','b');
      rh=boxplot_ts(stn.(climsrfld),'month','mean',true,'allcol','r');
      legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux/ISCCP');
      titlename([STNM,' Net insolation Q_S_W',stn.commentstr]);
      ylabel('Wm^-^2'); ylim([-1000,1000]);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',sr1dfld,'-vs-',climsrfld,'.tiff']));
      end;

      fmg;
      lh=boxplot_ts(stn.(lr1dfld),'month','mean',true,'allcol','b');
      rh=boxplot_ts(stn.(climlrfld),'month','mean',true,'allcol','r');
      legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux/ISCCP');
      titlename([STNM,' Net longwave flux Q_L_W',stn.commentstr]);
      ylabel('Wm^-^2'); ylim([-1000,1000]);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',lr1dfld,'-vs-',climlrfld,'.tiff']));
      end;

      fmg;
      lh=boxplot_ts(stn.(qlh1dfld),'month','mean',true,'allcol','b');
      rh=boxplot_ts(stn.(climqlhfld),'month','mean',true,'allcol','r');
      legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux');
      titlename([STNM,' Latent heat flux Q_L_H',stn.commentstr]);
      ylabel('Wm^-^2'); ylim([-1000,1000]);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',qlh1dfld,'-vs-',climqlhfld,'.tiff']));
      end;

      fmg;
      lh=boxplot_ts(stn.(qsh1dfld),'month','mean',true,'allcol','b');
      rh=boxplot_ts(stn.(climqshfld),'month','mean',true,'allcol','r');
      legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux');
      titlename([STNM,' Sensible heat flux Q_S_H',stn.commentstr]);
      ylabel('Wm^-^2'); ylim([-1000,1000]);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-boxplot-',qsh1dfld,'-vs-',climqshfld,'.tiff']));
      end;

      nd = ts_op(stn.(a1dfld),stn.(s1dfld),'-',@(x)(intersect_dates(x.date,stn.(climafld).date)));
      od = ts_op(stn.(climafld),stn.(climsfld),'-',@(x)(intersect_dates(x.date,stn.(a1dfld).date)));
      fmg; grpplot_ts(nd,@get_week,@nanmedian,0,'b'); grpplot_ts(od,@get_week,@nanmedian,0,'r');
      legend('In situ','OAFlux');
      titlename([STNM,' Air-Sea Temperature Difference']);
      ylabel('K'); ylim([-5,2]);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-median-air-sea-temperature-',ISPFX,'-vs-',CLIMPFX,'.tiff']));
      end;

    end; %if ( doBoxplots )


    if ( doDiagplots )

      % Compare magnitudes of advective heating vs. air-sea-benthic heating
      [bqt,udT] = intersect_tses(stn.(bq0tfld),stn.(udTfld));
      fmg; plot3(get_yearday(udT.date),udT.data,get_year(udT.date),'.',get_yearday(bqt.date),bqt.data,get_year(bqt.date),'r.'); datetick3; view(2); ylim([-.5,.5]); legend('u_s_f_c^.\nablaT_k_m','(Q_0+Q_b)/\rhoC_ph'); titlename(strrep([STNM,' ',udTfld,' VS. ',bq0tfld],'_','\_'));


      % Plot histograms of estimated representation error for terms of Q0
      fmg;
      pos=get(gcf,'Position'); pos(4)=pos(4)*0.40; set(gcf,'Position',pos);
      spt(1,4,1); hist(stn.([qlhfld,'_err']).data(stn.([qlhfld,'_err']).data<500),1000); xlim([0,80]); grid on;
      annotation('textbox',[0.20,0.85,.01,.01],'String','(a)','LineStyle','none','FontSize',14,'FontWeight','bold');
      spt(1,4,2); hist(stn.([qshfld,'_err']).data(stn.([qshfld,'_err']).data<500),1000); xlim([0,80]); grid on;
      annotation('textbox',[0.41,0.85,.01,.01],'String','(b)','LineStyle','none','FontSize',14,'FontWeight','bold');
      spt(1,4,3); hist(stn.([srfld,'_err']).data(0<stn.([srfld,'_err']).data&stn.([srfld,'_err']).data<500),1000); xlim([0,80]); grid on;
      annotation('textbox',[0.64,0.85,.01,.01],'String','(c)','LineStyle','none','FontSize',14,'FontWeight','bold');
      spt(1,4,4); hist(stn.([lrfld,'_err']).data(stn.([lrfld,'_err']).data<500),1000); xlim([0,80]); grid on;
      annotation('textbox',[0.86,0.85,.01,.01],'String','(d)','LineStyle','none','FontSize',14,'FontWeight','bold');
      suptitlename([STNM,' ',strrep(sq0fld,'_','\_'),' error histograms']);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,oldstr,'-error-histograms-',sq0fld,'.tiff']));
      end;


      % Plot balances between air-sea, aborption, and benthic exchange models
      %DEBUG:
      [sq,bq,bo,bl,dc]=intersect_tses(stn.(sq0fld),stn.(bq0fld),stn.(qbofld),stn.(asrdiagfld).lost_sr,stn.(qbdiagfld).qbcdd); fmg; plot(sq.date,cumsum(sq.data),bq.date,cumsum(bq.data),bo.date,cumsum(bo.data),bo.date,cumsum(bq.data+bl.data),bo.date,cumsum(bq.data+bl.data+dc.data)); datetick3; legend('Q_0','Q_0(\gamma)+Q_b','Q_b','Q_0(\gamma)+Q_b+(1-\gamma)Q_S_W','Q_0(\gamma)+Q_b+(1-\gamma)Q_S_W+Q_b_C_D^D'); titlename([STNM,' Flux Balances: ',strrep(bq0fld,'_','\_')]);

      if ( include_rain_flux )
        % Plot relative contributions to total error variance
        %q2 = ts_fun(stn.([q0fld,'_err']),@(x)(x.^2));
        bq2 = ts_fun(stn.([bq0fld,'_err']),@(x)(x.^2));
        aqsw2 = ts_fun(stn.([asrfld,'_err']),@(x)(x.^2));
        qlw2 = ts_fun(stn.([lrfld,'_err']),@(x)(x.^2));
        qlh2 = ts_fun(stn.([qlhfld,'_err']),@(x)(x.^2));
        qsh2 = ts_fun(stn.([qshfld,'_err']),@(x)(x.^2));
        qbo2 = ts_fun(stn.([qbofld,'_err']),@(x)(x.^2));
        % qcv2 = ts_op(bq2, ts_op(ts_op(aqsw2,qlw2,'+'), ts_op(qlh2,qsh2,'+'), '+'), '-');
        fmg; grpplot_ts(ts_op(aqsw2,bq2,'/'),[],@nanmean,[],'m'); grpplot_ts(ts_op(qlw2,bq2,'/'),[],@nanmean,[],'r'); grpplot_ts(ts_op(qlh2,bq2,'/'),[],@nanmean,[],'g'); grpplot_ts(ts_op(qsh2,bq2,'/'),[],@nanmean,[],'b'); grpplot_ts(ts_op(qbo2,bq2,'/'),[],@nanmean,[],'k'); ylim([-0.5,1.1]); legend('\sigma\gammaQ_S_W/\Sigma','\sigmaQ_L_W/\Sigma','\sigmaQ_L_H/\Sigma','\sigmaQ_S_H/\Sigma','\sigmaQ_b/\Sigma', 'Location','Best'); titlename([STNM,' error variance contribution (',strrep(q0fld,'_','\_'),')']);
        if ( doPrint )
          print('-dtiff',fullfile(figspath,[stnm,oldstr,'-relative-error-variances-',q0fld,'.tiff']));
        end;

        % % Plot individual components of 1-d flux error
        % fmg; grpplot_ts(stn.([asrfld,'_err']),[],[],[],'m'); grpplot_ts(stn.([lrfld,'_err']),[],[],[],'r'); grpplot_ts(stn.([qlhfld,'_err']),[],[],[],'g'); grpplot_ts(stn.([qshfld,'_err']),[],[],[],'b'); grpplot_ts(stn.([qbofld,'_err']),[],[],[],'c'); grpplot_ts(stn.([bq0fld,'_err']),[],[],[],'k'); legend('\sigma\gammaQ_S_W','\sigmaQ_L_W','\sigmaQ_L_H','\sigmaQ_S_H','\sigmaQ_b','\sigma(Q_0(\gamma)+Q_b)', 'Location','Best'); titlename([STNM,' individual error components (',strrep(q0fld,'_','\_'),')']);
        % if ( doPrint )
        %   print('-dtiff',fullfile(figspath,[stnm,oldstr,'-individual-errors-',bq0fld,'.tiff']));
        % end;

        % Plot individual components of total error
        % sother = ts_op( ts_op(stn.([fqudTffld,'_err']),stn.([kd2Tffld,'_err']),'+'), ...
        %                 stn.([hcdTdthcf,'_err']), '+');
        sother = stn.([hcdTdthcf,'_err']);
        fmg; grpplot_ts(stn.([asrfld,'_err']),[],[],[],'m'); grpplot_ts(stn.([lrfld,'_err']),[],[],[],'r'); grpplot_ts(stn.([qlhfld,'_err']),[],[],[],'g'); grpplot_ts(stn.([qshfld,'_err']),[],[],[],'b'); grpplot_ts(stn.([qbofld,'_err']),[],[],[],'y','Color',[.8,.8,.5]); grpplot_ts(stn.([bq0fld,'_err']),[],[],[],'k'); grpplot_ts(sother,[],[],[],'c'); grpplot_ts(stn.([hcdTdtf,'_err']),[],[],[],'k--','LineWidth',2); legend('\sigma\gammaQ_S_W','\sigmaQ_L_W','\sigmaQ_L_H','\sigmaQ_S_H','\sigmaQ_b','\sigma(Q_0(\gamma)+Q_b)','\sigma(HC)','\partial_tT', 'Location','Best'); titlename([STNM,' individual error components (',strrep(hcdTdt,'_','\_'),')']);
        ylim([0,150]); xlabel('Year-Month'); ylabel('Error [W/m^-^2]');
        if ( doPrint )
          print('-dtiff',fullfile(figspath,[stnm,oldstr,'-individual-errors-',hcdTdt,'.tiff']));
        end;
      end; %if ( include_rain_flux )

      % Diagnostics: Annual amplitude of horizontal convection
      fmg; boxplot_ts(stn.(hcdTdtf)); titlename([STNM,' \partial_tT^.\rhoC_ph']); ylim([-500,500]);
      fmg; boxplot_ts(stn.(hcu)); titlename([STNM,' u_H_C']); ylim([0,1]);
      fmg; boxplot_ts(stn.(hcdTdx)); titlename([STNM,' \nablaT_H_C']); ylim([-3e-4,3e-4]);

    end; %if ( doDiagplots )

  end; %if ( doPlot )


  if ( doClimplots )
    % Regress individual budget terms vs. observed variability (for Table 4)
    stn.budget_regression = ...
        regress_temp_budget_terms(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,doClimplots,true);
  end;


  if (doSave)
    matfname = fullfile(datapath,[stnm '_' taufld '.mat']);
    if ( ~exist(matfname,'file') )
      date = stn.(taufld).date;
      data = stn.(taufld).data;
      save(matfname,'date','data');
      clear date data;
    end;

    matfname = fullfile(datapath,[stnm '_' hcdTdt '.mat']);
    if ( ~exist(matfname,'file') )
      station=[]; clear station;
      station.station_name = stn.station_name;
      station.lon = stn.lon; station.lat = stn.lat; station.depth = stn.depth;
      flds = grepstruct(stn,['(',TURPFX,'|',QEPFX,')']);
      for cfld=flds'
        station.(cfld{:}) = stn.(cfld{:});
      end;
      disp(['Saving result to ',matfname]);
      save(matfname,'station');
      station=[]; clear station;
    end;
  end; %if (doSave)

  %% Print a brief report (as for Table 4 of Gramer & Mariano "... Heat Budget") 
  dump_robust_fit(stn, dly_dsfld, {...
      dly_srtfld,dly_asrtfld,dly_lrtfld,dly_qlhtfld,dly_qshtfld,dly_sqtfld,dly_bq0tfld,dly_bdTfld,dly_hcdTdt});

  if ( nargout < 1 )
    stn=[]; clear stn
  end;

  toc,
  timenow,
  set_more;

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newts = oshb_daily_sum(ts)
  newts.date = [];
  newts.data = [];
  [newts.data,newts.date] = grp_ts(ts.data,ts.date,@floor,@nansum,24);
return;
