function testab(stn)


  %%%
  %% Variable-name prefixes ("PFX") for various input and output datasets

  % ISPFX - In situ (station) data
  ISPFX = 'ndbc';

  % TIDEPFX - Tidal heights and currents
  TIDEPFX = 'tmd_tide';

  % RAPFX - Atmospheric reanalysis
  if ( ~exist('RAPFX','var') || isempty(RAPFX) )
    RAPFX = 'erai';	%Or 'ncep'
  end;

  % TURPFX - Turbulent surface fluxes
  % With warm-layer adjustment
  TURPFX = [ISPFX '_' RAPFX '_30a'];
  % WITHOUT warm-layer adjustment (just cold-skin and rain)
  TUR30PFX = [ISPFX '_' RAPFX '_30'];

  % WAVEPFX - Surface waves
  WAVEPFX = 'ww3';

  % KMPFX - Kilometer-scale ocean currents and sea temperature
  if ( ~exist('KMPFX','var') || isempty(KMPFX) )
    KMPFX = 'fkeys_hycom';	%Or 'gom_hycom'
  end;

  % STOKESPFX - Surface and near-surface currents (mean, Stokes drift)
  STOKESPFX = [WAVEPFX '_' ISPFX '_stokes'];

  % QEPFX - Quasi-Eulerian currents and heat advection
  switch (KMPFX),
   case 'fkeys_hycom',	QEPFX = [WAVEPFX '_fkeys_qe'];
   case 'gom_hycom',	QEPFX = [WAVEPFX '_gom_qe'];
   otherwise,		error('Unknown km-scale model "%s"',KMPFX);
  end;



  %%%
  %% All station struct fieldnames used to produce heat budget 

  % Tide data (or model)
  hfld = [TIDEPFX '_i_depth'];
  tufld = [TIDEPFX '_u'];
  tvfld = [TIDEPFX '_v'];

  % Meteorology (in situ)
  afld = [ISPFX '_air_t'];
  sfld = [ISPFX '_sea_t'];

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

  qsfld = [ISPFX '_sea_spechumid'];

  Wfld = [ISPFX '_wind1_speed'];
  Dfld = [ISPFX '_wind1_dir'];
  Ufld = [ISPFX '_wind1_u'];
  Vfld = [ISPFX '_wind1_v'];


  % Meteorology (gridded/reanalysis)
  rfld = [RAPFX '_precip'];
  cfld = [RAPFX '_cloud_cover'];
  pblzfld = [RAPFX '_pblz'];


  % Radiative surface fluxes
  dsrfld = [RAPFX '_dsrf'];
  usrfld = [RAPFX '_usrf'];
  srfld = [RAPFX '_srf'];
  asrfld = ['absorbed_' RAPFX '_srf'];
  gamfld = ['absorbed_' RAPFX '_gamma'];

  % BULK or REANALYSIS for long-wave???
  dlrfld = [RAPFX '_' ISPFX '_dlrf'];
  ulrfld = [RAPFX '_' ISPFX '_ulrf'];
  lrfld = [RAPFX '_' ISPFX '_lrf'];
  % dlrfld = [RAPFX '_dlrf'];
  % ulrfld = [RAPFX '_ulrf'];
  % lrfld = [RAPFX '_lrf'];

  % Water-benthos fluxes
  qbfld = ['benthic_' RAPFX '_srf'];
  btfld = ['benthic_' RAPFX '_t'];
  qbofld = ['benthic_' RAPFX '_qbo'];


  % Turbulent and net surface fluxes
  qlhfld = [TURPFX '_latent_heat_flux'];
  qshfld = [TURPFX '_sensible_heat_flux'];
  qrhfld = [TURPFX '_rain_heat_flux'];

  q0fld = [TURPFX '_net_heat_flux'];
  qtfld = [q0fld '_term'];

  % Fluxes without warm-layer adjustment
  q030fld = [TUR30PFX '_net_heat_flux'];
  qt30fld = [q030fld '_term'];


  % Ocean processes
  whfld = [WAVEPFX '_sigwavehgt'];
  wpfld = [WAVEPFX '_peakwaveper'];
  wdfld = [WAVEPFX '_peakwavedir'];

  ufld = [KMPFX '_u'];
  vfld = [KMPFX '_v'];
  Tfld = [KMPFX '_seatemp_field'];

  % Hourly fit to native data
  hufld = ['hourly_' KMPFX '_u'];
  hvfld = ['hourly_' KMPFX '_v'];


  % Near-bottom currents (mean + tide): used both for benthic heat exchange
  % and (friction velocity-dependent) horizontal convection calculations
  netufld = [TIDEPFX '_' KMPFX '_u'];
  netvfld = [TIDEPFX '_' KMPFX '_v'];


  % Surface and near-surface currents (mean, Stokes drift)
  ssufld = [STOKESPFX '_u'];
  ssvfld = [STOKESPFX '_v'];
  sssfld = [STOKESPFX '_speed'];
  ssdfld = [STOKESPFX '_dir'];


  % Quasi-Eulerian currents and heat advection
  % Careful not to double-add advection by km-scale currents!
  qeufld = [QEPFX '_u'];
  qevfld = [QEPFX '_v'];
  qesfld = [QEPFX '_speed'];
  qedfld = [QEPFX '_dir'];

  udTfld = [QEPFX '_advected_heat'];

  qtAdvfld = [TURPFX '_' QEPFX '_qtadv'];

  % Km-scale heat diffusion
  switch (KMPFX),
   case 'fkeys_hycom',	K_theta = 2.5;
   case 'gom_hycom',	K_theta = 20;
   otherwise,		error('Unknown km-scale model "%s"',KMPFX);
  end;
  kd2Tfld = [KMPFX '_diffused_heat'];


  % Total budget
  dTfld = [TURPFX '_' QEPFX '_dt'];
  dTffld = [dTfld '_heat_flux'];

  % Total budget with benthic flux
  qbotfld = [qbofld '_term'];
  bq0fld = ['benthic_' q0fld];
  bq0tfld = [bq0fld '_term'];
  bdTfld = ['benthic_' dTfld];
  bdTffld = [bdTfld '_heat_flux'];

  % Total budget with benthic flux and horizontal convection
  HCPFX = 'hc_';


  % figure; maxigraph; plot_ts(stn.erai_dsrf,'go',stn.erai_srf,'ro',stn.absorbed_erai_srf,'b.');
  % axis([datenum(2005,4,1),datenum(2005,4,12),-200,1000]); datetick3('x',2,'keeplimits');
  % hold on;
  % % xlim([datenum(2005,4,1),datenum(2005,4,12)]); datetick3('x',2,'keeplimits');
  % % plot(stn.absorbed_erai_gamma.date,stn.absorbed_erai_gamma.data*100,'k');
  % titlename('Before');

  stn = station_absorbed_insolation(stn,asrfld,srfld,hfld,[],[],gamfld,qbfld,'aidiag');

  figure; maxigraph; plot_ts(stn.erai_srf,'ro',stn.absorbed_erai_srf,'b.');
  % axis([datenum(2005,4,1),datenum(2005,4,12),-200,1000]); datetick3('x',2,'keeplimits');
  hold on;
  % xlim([datenum(2005,4,1),datenum(2005,4,12)]); datetick3('x',2,'keeplimits');
  % plot(stn.absorbed_erai_gamma.date,stn.absorbed_erai_gamma.data*100,'m');
  % plot(stn.aidiag.sun_alt.date(stn.aidiag.sun_alt.data>-20),stn.aidiag.sun_alt.data(stn.aidiag.sun_alt.data>-20)*10,'yo');
  % plot(stn.aidiag.sun_alt_correction.date,stn.aidiag.sun_alt_correction.data*100,'ko');
  % plot(stn.aidiag.tau.date,stn.aidiag.tau.data*100,'ko','Color',[.5,.5,.5]);
  plot(stn.ndbc_sea_t.date,(stn.ndbc_sea_t.data-25).*600,'k-','color',[.5,.5,.5]);
  plot(stn.bic_surf_par.date+(0.5/24),stn.bic_surf_par.data.*0.473,'k-');
  titlename('FINAL');

return;
