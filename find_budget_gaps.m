function gapflds = find_budget_gaps(stn,myDt)

  hfld = 'tmd_tide_i_depth';

  % Meteorology (in situ)
  afld = 'ndbc_air_t';
  sfld = 'ndbc_sea_t';
  pfld = 'ndbc_barom';
  dfld = 'ndbc_dew_t';
  Wfld = 'ndbc_wind1_speed';
  Dfld = 'ndbc_wind1_dir';
  Ufld = 'ndbc_wind1_u';
  Vfld = 'ndbc_wind1_v';

  if ( ~exist('RAPFX','var') || isempty(RAPFX) )
    RAPFX = 'erai';
    % RAPFX = 'ncep';
  end;

  % Meteorology (gridded/reanalysis)
  %rhfld = 'ndbc_relhumid';
  %qafld = 'ndbc_spechumid';
  %qsfld = 'ndbc_sea_spechumid';

  rhfld = [RAPFX '_relhumid'];
  qafld = [RAPFX '_spechumid'];
  qsfld = [RAPFX '_sea_spechumid'];

  rfld = [RAPFX '_precip'];
  cfld = [RAPFX '_cloud_cover'];
  alt_dfld = [RAPFX '_dew_t'];
  alt_rhfld = [RAPFX '_relhumid'];
  alt_qafld = [RAPFX '_spechumid'];
  alt_qsfld = [RAPFX '_sea_spechumid'];

  % Radiative surface fluxes
  dsrfld = [RAPFX '_dsrf'];
  usrfld = [RAPFX '_usrf'];
  srfld = [RAPFX '_srf'];
  asrfld = ['absorbed_' RAPFX '_srf'];

  dlrfld = [RAPFX '_ndbc_dlrf'];
  ulrfld = [RAPFX '_ndbc_ulrf'];
  lrfld = [RAPFX '_ndbc_lrf'];


  % Turbulent and net surface fluxes
  TURPFX = ['ndbc_' RAPFX '_30a'];

  qlhfld = [TURPFX '_latent_heat_flux'];
  qshfld = [TURPFX '_sensible_heat_flux'];
  qrhfld = [TURPFX '_rain_heat_flux'];

  q0fld = [TURPFX '_net_heat_flux'];
  qtfld = [TURPFX '_net_heat_flux_term'];

  % Ocean processes
  WAVEPFX = 'ww3';
  whfld = [WAVEPFX '_sigwavehgt'];
  wpfld = [WAVEPFX '_peakwaveper'];
  wdfld = [WAVEPFX '_peakwavedir'];

  if ( ~exist('KMPFX','var') || isempty(KMPFX) )
    KMPFX = 'fkeys_hycom';
  end;
  ufld = [KMPFX '_u'];
  vfld = [KMPFX '_v'];
  Tfld = [KMPFX '_seatemp_field'];

  STOKESPFX = [WAVEPFX '_ndbc_stokes'];
  ssufld = [STOKESPFX '_u'];
  ssvfld = [STOKESPFX '_v'];
  sssfld = [STOKESPFX '_speed'];
  ssdfld = [STOKESPFX '_dir'];

  % Quasi-Eulerian currents and heat advection
  % Careful not to double-add advection by km-scale currents!
  QEPFX = 'ww3_fkeys_qe';
  qeufld = [QEPFX '_u'];
  qevfld = [QEPFX '_v'];
  qesfld = [QEPFX '_speed'];
  qedfld = [QEPFX '_dir'];

  udTfld = [QEPFX '_advected_heat'];

  qtAdvfld = [TURPFX '_' QEPFX '_qtadv'];

  % Km-scale heat diffusion
  K_theta = 2.5;
  kd2Tfld = [KMPFX '_diffused_heat'];


  % Total budget
  dTfld = [TURPFX '_' QEPFX '_dt'];


  if ( ~exist('myDt','var') || isempty(myDt) )
    myDt = datenum(2007,6,1,0,0,0);
  end;
  mydelDt = 1;

  % flds = { qeufld,Tfld,udTfld,qtfld,qtAdvfld,kd2Tfld,qtAdvfld,dTfld,sfld,hfld };
  flds = {Wfld,afld,rhfld,pfld,sfld,asrfld,lrfld,dsrfld,dlrfld,rfld,Dfld,qeufld,qevfld,wpfld,whfld};

  gapflds = {};
  for fldix = 1:length(flds)
    fld = flds{fldix};
    ix = find(abs(myDt-stn.(fld).date)<mydelDt,1);
    if ( isempty(ix) )
      disp(fld);
      gapflds{end+1} = fld;
    end;
  end;


return;
