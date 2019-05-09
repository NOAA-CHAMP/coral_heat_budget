1;
%%%% SCRIPT to be called from an M-function or the command line, to set:
%% Variable-name prefixes ("PFX") for various input and output datasets; AND,
%% All station struct fieldnames used to produce the ocean heat budget. Only
%% effect of this script is to set a large number of (mostly CHAR) variables.
%%
%% NOTE WELL: Calls FIX_VARNAMELENGTHS (v.) to fix variable name string
%%             lengths to meet MATLAB limitations. MAJOR SIDE EFFECTS...
%%
%% VARIABLES WHICH MAY BE PRESET TO CONTROL THIS SCRIPT:
%%   ISPFX = 'ndbc';
%%   SPREFIX = '';
%%   TIDEPFX = 'tmd_tide';
%%   RAPFX = 'erai';	%Or 'ncep'
%%   TURPFX = [SPREFIX ISPFX '_' RAPFX '_30a'];
%%   TUR30PFX = [SPREFIX ISPFX '_' RAPFX '_30'];
%%   CLIMPFX = 'daily_oaflux';
%%   WAVEPFX = 'erai';
%%   KMPFX = 'gom_hycom';	%OR 'fkeys_hycom' OR 'avhrr_weekly' OR 'avhrr' OR 'none'
%%   STOKESPFX = [WAVEPFX '_' ISPFX '_stokes'];
%%   QEPFX = [WAVEPFX '_gom'];
%%   Q0_LOWPASS = '_48_h_lp';
%%   QE_LOWPASS = '_72_h_lp';
%%   HCPFX = [TURPFX '_' QEPFX '_hc'];
%%   CLIMHCPFX = [CLIMPFX '_' HCPFX];
%%
%% Last Saved Time-stamp: <Sun 2013-07-21 15:33:00 Eastern Daylight Time gramer>


  %%%
  %% Do any initial fieldname substitutions - will also redo at end
  if ( ~exist('substitute_field_names','var') )
    substitute_field_names = {};
  end;
  for ix=1:2:numel(substitute_field_names)
    assignin('caller',substitute_field_names{ix},substitute_field_names{ix+1});
  end;

  %%%
  %% Variable-name prefixes ("PFX") for various input and output datasets

  % ISPFX - In situ (station) data
  if ( ~exist('ISPFX','var') || isempty(ISPFX) )
    ISPFX = 'ndbc';
  end;

  if ( ~exist('use_old_options','var') || isempty(use_old_options) )
    use_old_options = false;
  end;
  if ( ~exist('use_old_erai_only_options','var') || isempty(use_old_erai_only_options) )
    use_old_erai_only_options = false;
  end;
  if ( use_old_erai_only_options )
    warning('Reproducing old ISPFX==ERAI results');
    begyr = 1996;
    endyr = 2010;
    adjust_waves = false;
    adjust_reanalysis = false;
    reanalysis_shortwave = true;
    reanalysis_longwave = true;
  end;

  % Special handling - if sea temperature SFLD is NOT from source ISFPFX
  if ( ~exist('SPREFIX','var') || isempty(SPREFIX) )
    SPREFIX = '';
    if ( exist('sfld','var') && ~isempty(sfld) && ~strcmpi(sfld,[ISPFX '_sea_t'])  && ~strcmpi(sfld,[ISPFX '_seatemp']) )
      SPREFIX = [strrep(sfld,'_','') '_'];
    end;
  end;

  % TIDEPFX - Tidal heights and currents
  if ( ~exist('TIDEPFX','var') || isempty(TIDEPFX) )
    TIDEPFX = 'tmd_tide';
  end;

  % RAPFX - Atmospheric reanalysis
  if ( ~exist('RAPFX','var') || isempty(RAPFX) )
    RAPFX = 'erai';	%Or 'ncep'
  end;

  % WAVEPFX - Surface waves
  if ( ~exist('adjust_waves','var') )
    adjust_waves = true;
  end;
  if ( ~exist('WAVEPFX','var') || isempty(WAVEPFX) )
    % WAVEPFX = 'erai';
    % WAVEPFX = 'ndbc';
    WAVEPFX = 'ww3';
  end;

  % TURPFX - Turbulent surface fluxes
  if ( ~exist('include_rain_flux','var') )
    include_rain_flux = true;
  end;
  if ( include_rain_flux )
    % With warm-layer adjustment
    if ( ~exist('TURPFX','var') || isempty(TURPFX) )
      TURPFX = [SPREFIX ISPFX '_' RAPFX '_' WAVEPFX '_30a'];
    end;
    % WITHOUT warm-layer adjustment (just rain and cold-skin)
    if ( ~exist('TUR30PFX','var') || isempty(TUR30PFX) )
      TUR30PFX = [SPREFIX ISPFX '_' RAPFX '_' WAVEPFX '_30'];
    end;
  else
    % WITHOUT warm-layer adjustment or rain (just cold-skin)
    if ( ~exist('TURPFX','var') || isempty(TURPFX) )
      TURPFX = [SPREFIX ISPFX '_' RAPFX '_' WAVEPFX '_26'];
    end;
    % WITHOUT warm-layer adjustment, rain, or cold-skin
    if ( ~exist('TUR30PFX','var') || isempty(TUR30PFX) )
      TUR30PFX = [SPREFIX ISPFX '_' RAPFX '_' WAVEPFX '_20'];
    end;
  end;

  % Climatological net surface fluxes - for comparison 
  if ( ~exist('CLIMPFX','var') || isempty(CLIMPFX) )
    CLIMPFX = 'daily_oaflux';
  end;

  % KMPFX - Kilometer-scale ocean currents and sea temperature
  if ( ~exist('KMPFX','var') || isempty(KMPFX) )
    KMPFX = 'gom_hycom';	%OR 'fkeys_hycom' OR 'avhrr_weekly' OR 'avhrr' OR 'none'
  end;

  % STOKESPFX - Surface and near-surface currents (mean, Stokes drift)
  if ( ~exist('STOKESPFX','var') || isempty(STOKESPFX) )
    STOKESPFX = [WAVEPFX '_' ISPFX '_stokes'];
  end;

  % QEPFX - Quasi-Eulerian currents and heat advection
  if ( ~exist('QEPFX','var') || isempty(QEPFX) )
    switch (KMPFX),
     case 'fkeys_hycom',	QEPFX = [WAVEPFX '_fkeys'];
     case 'gom_hycom',		QEPFX = [WAVEPFX '_gom'];
     case 'avhrr_weekly',	QEPFX = [WAVEPFX '_avhrr'];
     case 'avhrr',		QEPFX = [WAVEPFX '_avhrrd'];
     case 'none',		QEPFX = [WAVEPFX '_none'];
     otherwise,			error('Unknown km-scale model "%s"',KMPFX);
    end;
  end;

  % DTPFX - Combined Heat Budget prefix - Rad/Tur + Adv/Diff
  if ( ~exist('DTPFX','var') || isempty(DTPFX) )
    switch (KMPFX),
     case 'fkeys_hycom',	DTPFX = [TURPFX '_fkeys'];
     case 'gom_hycom',		DTPFX = [TURPFX '_gom'];
     case 'avhrr_weekly',	DTPFX = [TURPFX '_avhrr'];
     case 'avhrr',		DTPFX = [TURPFX '_avhrrd'];
     case 'none',		DTPFX = [TURPFX '_none'];
     otherwise,			error('Unknown km-scale model "%s"',KMPFX);
    end;
  end;

  % Special options
  if ( ~exist('Q0_LOWPASS','var') )
    if ( use_old_options )
      Q0_LOWPASS = '_48_h_lp';
    else
      Q0_LOWPASS = '_24_h_lp';
    end;
  end;
  if ( ~exist('QE_LOWPASS','var') )
    QE_LOWPASS = '_72_h_lp';
  end;


  %%%
  %% All station struct fieldnames used to produce heat budget 

  bathyfld = 'ngdc_92m_bathy';
  slopefld = 'ngdc_offshore_slope';
  bathorifld = 'isobath_orientation';

  % Tide data (or model)
  hfld = [TIDEPFX '_i_depth'];
  mhfld = ['mean_' hfld];
  tufld = [TIDEPFX '_u'];
  tvfld = [TIDEPFX '_v'];
  tspdfld = [TIDEPFX '_speed'];
  tdirfld = [TIDEPFX '_dir'];

  % Sea temperature (normally, in situ from ISPFX source)
  if ( ~exist('sfld','var') || isempty(sfld) )
    sfld = [ISPFX '_sea_t'];
  end;
  dsfld = [sfld '_diff'];
  dsffld = [dsfld '_flux'];
  % "Cool-skin" temperature (from TOGA-COARE 3.0a or HFBULKTC)
  scoolfld = [sfld,'_cool'];

  % Meteorology (in situ)
  if ( ~exist('afld','var') || isempty(afld) )
    afld = [ISPFX '_air_t'];
  end;

  %%%% Air pressure and humidity
  %%%% In situ records often missing or incomplete. Air-sea flux algorithms
  %%%% are quite insensitive to normal ranges in these, so use reanalysis.
  % pfld = [ISPFX '_barom'];
  % dfld = [ISPFX '_dew_t'];
  % rhfld = [ISPFX '_relhumid'];
  % qafld = [ISPFX '_spechumid'];
  pfld = [RAPFX '_barom'];
  switch (RAPFX),
   case 'ncep',		dfld = [RAPFX '_dewp'];
   otherwise,		dfld = [RAPFX '_dew_t'];
  end;
  rhfld = [RAPFX '_relhumid'];
  qafld = [RAPFX '_spechumid'];

  qsfld = [SPREFIX ISPFX '_sea_spechumid'];

  switch (ISPFX),
   case 'ndbc',		WINDINFIX = '_wind1';
   otherwise,		WINDINFIX = '_wind';
  end;
  Wfld = [ISPFX WINDINFIX '_speed'];
  Dfld = [ISPFX WINDINFIX '_dir'];
  Ufld = [ISPFX WINDINFIX '_u'];
  Vfld = [ISPFX WINDINFIX '_v'];

  %%%% ??? DEBUG: Low-pass filter winds for quasi-Eulerian currents
  Ulpfld = [Ufld QE_LOWPASS];
  Vlpfld = [Vfld QE_LOWPASS];
  Wlpfld = [Wfld QE_LOWPASS];
  Dlpfld = [Dfld QE_LOWPASS];
  %%%% ??? DEBUG: Low-pass filter winds for quasi-Eulerian currents


  % Meteorology (gridded/reanalysis)

  if ( include_rain_flux )
    rfld = [RAPFX '_precip'];
  else
    rfld=[];
  end;

  cfld = [RAPFX '_cloud_cover'];
  pblzfld = [RAPFX '_pblz'];


  % Radiative surface fluxes

  % BULK or REANALYSIS for (*upward*) shortwave radiation?
  if ( ~exist('reanalysis_shortwave','var') )
    reanalysis_shortwave = false;
  end;
  if ( reanalysis_shortwave )
    dsrfld = [RAPFX '_dsrf'];
    usrfld = [RAPFX '_usrf'];
    albfld = [RAPFX '_albedo'];
    srfld  = [RAPFX '_srf'];
    asrfld = ['absorbed_' RAPFX '_srf'];
    gamfld = ['absorbed_' RAPFX '_gamma'];
    asrdiagfld = ['absorbed_' RAPFX '_diag'];
  else
    dsrfld = [RAPFX '_dsrf'];
    usrfld = [RAPFX '_' ISPFX '_usrf'];
    albfld = [RAPFX '_' ISPFX '_albedo'];
    srfld  = [RAPFX '_' ISPFX '_srf'];
    asrfld = ['absorbed_' srfld];
    gamfld = ['absorbed_' RAPFX '_' ISPFX '_gamma'];
    asrdiagfld = ['absorbed_' RAPFX '_' ISPFX '_diag'];
  end;

  % BULK or REANALYSIS for (*upward*) longwave radiation?
  if ( ~exist('reanalysis_longwave','var') )
    reanalysis_longwave = false;
  end;
  if ( reanalysis_longwave )
    dlrfld = [RAPFX '_dlrf'];
    ulrfld = [RAPFX '_ulrf'];
    lrfld = [RAPFX '_lrf'];
  else
    dlrfld = [RAPFX '_dlrf'];
    ulrfld = [SPREFIX RAPFX '_' ISPFX '_ulrf'];
    lrfld = [SPREFIX RAPFX '_' ISPFX '_lrf'];
  end;

  srtfld = [srfld '_term'];
  asrtfld = [asrfld '_term'];
  lrtfld = [lrfld '_term'];

  % Water-benthos fluxes
  if ( ~exist('ignore_benthos','var') )
    ignore_benthos = false;
  end;
  qbfld = ['b_' RAPFX '_srf'];
  btfld = ['b_' RAPFX '_t'];
  qbofld = ['b_' RAPFX '_qbo'];
  qbdiagfld = ['b_' RAPFX '_diag'];

  % Should we apply empirical corrections to reanalysis data?
  if ( ~exist('adjust_reanalysis','var') )
    adjust_reanalysis = true;
  end;
  if ( adjust_reanalysis )
    if ( include_rain_flux )
      rfld = [rfld '_adj'];
    end;
    dsrfld = [dsrfld '_adj'];
    dlrfld = [dlrfld '_adj'];
  end;


  % Turbulent and net surface fluxes
  qlhfld = [TURPFX '_latent_flux'];
  qshfld = [TURPFX '_sensible_flux'];
  qrhfld = [TURPFX '_rain_flux'];
  % Diagnostics and intermediate products of TOGA/COARE algorithm
  diagfld = [TURPFX '_cordiags'];

  qlhtfld = [qlhfld '_term'];
  qshtfld = [qshfld '_term'];
  qrhtfld = [qrhfld '_term'];

  adj_qlhfld = ['adj_' qlhfld];

  taufld = [TURPFX '_wind_stress'];
  tauxfld = [taufld '_u'];
  tauyfld = [taufld '_v'];


  % Absorbed short-wave + (BULK or REANALYSIS for long-wave?)
  qradfld = [SPREFIX RAPFX '_' ISPFX '_arf'];
  %qradfld = [RAPFX '_arf'];
  qradtfld = [qradfld '_term'];
  % Latent + Sensible + Rain fluxes
  qturfld = [TURPFX '_turbulent_flux'];
  qturtfld = [qturfld '_term'];
  % Long-wave (see above) + Latent + Sensible + Rain fluxes fluxes
  qcoolfld = [TURPFX '_cooling_flux'];
  qcooltfld = [qcoolfld '_term'];

  q0fld = [TURPFX '_net_flux'];
  q0lpfld = [q0fld Q0_LOWPASS];
  qtfld = [q0fld '_term'];

  % Fluxes without warm-layer adjustment
  qlh30fld = [TUR30PFX '_latent_flux'];
  qsh30fld = [TUR30PFX '_sensible_flux'];
  qrh30fld = [TUR30PFX '_rain_flux'];
  qtur30fld = [TUR30PFX '_turbulent_flux'];
  q030fld = [TUR30PFX '_net_flux'];
  qt30fld = [q030fld '_term'];
  q030lpfld = [q030fld Q0_LOWPASS];

  % Fluxes without absorption correction
  % (BULK or REANALYSIS for long-wave component of net radiation?)
  sqradfld = [SPREFIX RAPFX '_' ISPFX '_rf'];
  %sqradfld = [RAPFX '_rf'];
  sq0fld = ['simple_' q0fld];
  sqtfld = ['simple_' qtfld];
  % ... AND without warm-layer adjustment
  sq030fld = ['simple_' q030fld];
  sqt30fld = ['simple_' qt30fld];


  % Climatological inputs and net surface fluxes - for comparison 
  climsfld = [CLIMPFX '_seatemp'];
  climafld = [CLIMPFX '_air_t'];
  climsrfld = [CLIMPFX '_srf'];
  climasrfld = ['absorbed_' CLIMPFX '_srf'];
  climlrfld = [CLIMPFX '_lrf'];
  climevapfld = [CLIMPFX '_evap'];
  climqlhfld = [CLIMPFX '_latent_heat_flux'];
  climqshfld = [CLIMPFX '_sensible_heat_flux'];
  climq0fld = [CLIMPFX '_net_heat_flux'];
  climqtfld = [climq0fld '_term'];

  raevapfld = [RAPFX '_evap'];
  raqlhfld = [RAPFX '_latent_flux'];
  raqshfld = [RAPFX '_sensible_flux'];
  % This field name is normally proscribed by the Reanalysis code...
  %raq0fld = [RAPFX '_net_flux'];
  raq0fld = [RAPFX '_net_heat_flux'];
  raqtfld = [raq0fld '_term'];


  % Ocean processes
  whfld = [WAVEPFX '_sigwavehgt'];
  wpfld = [WAVEPFX '_peakwaveper'];
  wdfld = [WAVEPFX '_peakwavedir'];
  if ( adjust_waves )
    whfld = [WAVEPFX '_sigwavehgt_adj'];
    wpfld = [WAVEPFX '_peakwaveper_adj'];
    wdfld = [WAVEPFX '_peakwavedir_adj'];
  end;

  ufld = [KMPFX '_u'];
  vfld = [KMPFX '_v'];
  switch (KMPFX),
   case {'avhrr_weekly','avhrr'},
    Tfld = [KMPFX '_sst_field'];
    kmtfld = [KMPFX '_sst'];
   otherwise,
    Tfld = [KMPFX '_seatemp_field'];
    kmtfld = [KMPFX '_seatemp'];
  end;
  kmtxfld = [kmtfld '_x'];
  kmtyfld = [kmtfld '_y'];
  kmtlfld = [kmtfld '_l'];
  kmtxsfld = [kmtfld '_xshore'];
  kmtlsfld = [kmtfld '_lshore'];

  % Hourly fit to native data
  hufld = ['hourly_' ufld];
  hvfld = ['hourly_' vfld];
  hkmtfld = ['hourly_' kmtfld];
  hkmtxfld = ['hourly_' kmtxfld];
  hkmtyfld = ['hourly_' kmtyfld];
  hkmtlfld = ['hourly_' kmtlfld];
  hkmtxsfld = ['hourly_' kmtxsfld];
  hkmtlsfld = ['hourly_' kmtlsfld];


  % Near-bottom currents (mean + tide): used both for benthic heat exchange
  % and (friction velocity-dependent) horizontal convection calculations
  netufld = [TIDEPFX '_' KMPFX '_u'];
  netvfld = [TIDEPFX '_' KMPFX '_v'];

  netxsfld = [TIDEPFX '_' KMPFX '_xshore'];
  netlsfld = [TIDEPFX '_' KMPFX '_lshore'];


  % Surface and near-surface currents (mean, Stokes drift)
  ssufld = [STOKESPFX '_u'];
  ssvfld = [STOKESPFX '_v'];
  sssfld = [STOKESPFX '_speed'];
  ssdfld = [STOKESPFX '_dir'];

  ssxsfld = [STOKESPFX '_xshore'];
  sslsfld = [STOKESPFX '_lshore'];


  % Quasi-Eulerian currents and heat advection
  % Careful not to double-add advection by km-scale currents!
  qeufld = [QEPFX '_u'];
  qevfld = [QEPFX '_v'];
  qesfld = [QEPFX '_speed'];
  qedfld = [QEPFX '_dir'];

  qexsfld = [QEPFX '_xshore'];
  qelsfld = [QEPFX '_lshore'];

  udTfld = [QEPFX '_advected_heat'];
  udTffld = [udTfld '_flux'];
  rawudTfld = ['raw_' udTfld];
  rawudTffld = [rawudTfld '_flux'];

  % Result of applying "quality" factor
  fqudTfld = ['fq_',udTfld];
  fqudTffld = ['fq_',udTffld];

  udTxsfld = [QEPFX '_advected_heat_xshore'];
  udTfxsfld = [udTxsfld '_flux'];
  rawudTxsfld = ['raw_' udTxsfld];
  udTlsfld = [QEPFX '_advected_heat_lshore'];
  udTflsfld = [udTlsfld '_flux'];
  rawudTlsfld = ['raw_' udTlsfld];

  qtAdvfld = [DTPFX '_qtadv'];
  qtAdvffld = [qtAdvfld '_flux'];

  % Km-scale heat diffusion
  switch (KMPFX),
   case 'fkeys_hycom',		model_K_theta = 2.5;
   case 'gom_hycom',		model_K_theta = 20;
   case 'avhrr_weekly',		model_K_theta = 20;
   case 'avhrr',		model_K_theta = 20;
   case 'none',			model_K_theta = 20;
   otherwise,			error('Unknown km-scale model "%s"',KMPFX);
  end;
  kd2Tfld = [KMPFX '_diffused_heat'];
  kd2Tffld = [kd2Tfld '_flux'];
  rawkd2Tfld = ['raw_' kd2Tfld];

  % Total budget
  dTfld = [DTPFX '_dt'];
  dTffld = [dTfld '_flux'];
  dTflpfld = [dTffld Q0_LOWPASS];

  % Total vertical budget with benthic flux
  qbotfld = [qbofld '_term'];
  bq0fld = ['b_' q0fld];
  bq0tfld = [bq0fld '_term'];
  bq0lpfld = [bq0fld Q0_LOWPASS];

  bdTfld = ['b_' dTfld];
  bdTffld = [bdTfld '_flux'];
  bdTflpfld = [bdTffld Q0_LOWPASS];

  % Total budget with benthic flux and horizontal convection
  if ( ~exist('HCPFX','var') || isempty(HCPFX) )
    % HCPFX = 'hc';
    HCPFX = [DTPFX '_hc'];
  end;
  hcfactor = [HCPFX '_termFactor'];
  hcdTdt = [HCPFX '_dTdt'];
  hcdTdtf = [hcdTdt '_flux'];
  hcdTdx = [HCPFX '_dTdx'];
  hcdTdxf = [hcdTdx '_flux'];
  hcdTdthc = [hcdTdt 'hc'];
  hcdTdthcf = [hcdTdthc '_flux'];
  hcu = [HCPFX '_u'];


  % Sub-grid-scale heat diffusion (calculated as residual)
  sgskd2Tfld = [HCPFX '_sgs_diffused_heat'];
  % sgskfld = [HCPFX '_sgs_thermal_diffusivity'];
  sgskfld = [HCPFX '_sgs_K_theta'];
  sgsdTdt = [HCPFX '_sgs_final_budget'];


  % Climatological total heat budget terms - for comparison 
  climqtAdvfld = [CLIMPFX '_' QEPFX '_qtadv'];
  climdTfld = [CLIMPFX '_' QEPFX '_dt'];
  climdTffld = [climdTfld '_flux'];
  climbq0fld = ['b_' climq0fld];
  climbq0tfld = [climbq0fld '_term'];
  climbdTfld = ['b_' climdTfld];
  climbdTffld = [climbdTfld '_flux'];
  if ( ~exist('CLIMHCPFX','var') || isempty(CLIMHCPFX) )
    CLIMHCPFX = [CLIMPFX '_' HCPFX];
  end;
  climhcfactor = [CLIMHCPFX '_termFactor'];
  climhcdTdt = [CLIMHCPFX '_dTdt'];


  % For removal of budget outliers based on estimation error
  qlhfld1derr = [qlhfld,'_err_1_d_sum'];
  bq0fld1derr = [bq0fld,'_err_1_d_sum'];
  hcdTdt1derr = [hcdTdt,'_err_1_d_sum'];


  % Fields for comparison of daily averages, sums, and changes
  dly_sfld = [sfld '_dly'];
  dly_dsfld = [dly_sfld '_diff'];
  dly_dsffld = [dly_dsfld '_flux'];
  dly_srfld = [srfld '_dly'];
  dly_srtfld = [srtfld '_dly'];
  dly_asrfld = [asrfld '_dly'];
  dly_asrtfld = [asrtfld '_dly'];
  dly_lrfld = [lrfld '_dly'];
  dly_lrtfld = [lrtfld '_dly'];
  dly_qlhfld = [qlhfld '_dly'];
  dly_qlhtfld = [qlhtfld '_dly'];
  dly_qshfld = [qshfld '_dly'];
  dly_qshtfld = [qshtfld '_dly'];
  dly_q0fld = [q0fld '_dly'];
  dly_qtfld = [qtfld '_dly'];
  dly_sq0fld = [sq0fld '_dly'];
  dly_sqtfld = [sqtfld '_dly'];
  dly_bq0fld = [bq0fld '_dly'];
  dly_bq0tfld = [bq0tfld '_dly'];
  dly_bdTfld = [bdTfld '_dly'];
  dly_bdTffld = [bdTffld '_dly'];
  dly_hcdTdthc = [hcdTdthc '_dly'];
  dly_hcdTdthcf = [hcdTdthcf '_dly'];
  dly_hcdTdt = [hcdTdt '_dly'];
  dly_hcdTdtf = [hcdTdtf '_dly'];


  % Fields for climatological comparison of results
  s1dfld = [sfld '_1_d_avg'];
  a1dfld = [afld '_1_d_avg'];
  sq01dfld = [sq0fld '_1_d_avg'];
  sr1dfld = [srfld '_1_d_avg'];
  asr1dfld = [asrfld '_1_d_avg'];
  lr1dfld = [lrfld '_1_d_avg'];
  qlh1dfld = [qlhfld '_1_d_avg'];
  qsh1dfld = [qshfld '_1_d_avg'];
  qcool1dfld = [qcoolfld '_1_d_avg'];
  fqudTffld1d = [fqudTffld '_1_d_avg'];
  qbofld1d = [qbofld '_1_d_avg'];
  kd2Tffld1d = [kd2Tffld '_1_d_avg'];
  hcdTdtf1d = [hcdTdtf '_1_d_avg'];
  hcdTdthcf1d = [hcdTdthcf '_1_d_avg'];


  %%%
  %% Redo fieldname substitutions - in case they were overwritten above
  for ix=1:2:numel(substitute_field_names)
    assignin('caller',substitute_field_names{ix},substitute_field_names{ix+1});
  end;


  %%%
  %% In case any of the field name strings defined by the script are too long
  fix_varnamelengths;
