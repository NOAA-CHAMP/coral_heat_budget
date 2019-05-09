function stn = recalc_gamma(stn)
%ERROR('Replaced by STATION_HEAT_BUDGET');
error('Replaced by STATION_HEAT_BUDGET');
%function stn = recalc_gamma(stn)
%
% Recalculate heat budget terms, including insolation absorption factor
% ("gamma") for the site in station struct STN.
%
% Last Saved Time-stamp: <Tue 2011-04-19 10:39:15  Lew.Gramer>

  set_more off

  % Tide data (or model)
  TIDEPFX = 'tmd_tide';
  hfld = [TIDEPFX '_i_depth'];
  tufld = [TIDEPFX '_u'];
  tvfld = [TIDEPFX '_v'];

  % Meteorology (in situ)
  afld = 'ndbc_air_t';
  sfld = 'ndbc_sea_t';
  pfld = 'erai_barom';
  dfld = 'erai_dew_t';
  rhfld = 'erai_relhumid';
  qafld = 'erai_spechumid';
  qsfld = 'erai_sea_spechumid';
  Wfld = 'ndbc_wind1_speed';
  Dfld = 'ndbc_wind1_dir';
  Ufld = 'ndbc_wind1_u';
  Vfld = 'ndbc_wind1_v';


  RAPFX = 'erai';

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
  dlrfld = [RAPFX '_ndbc_dlrf'];
  ulrfld = [RAPFX '_ndbc_ulrf'];
  lrfld = [RAPFX '_ndbc_lrf'];
  % dlrfld = [RAPFX '_dlrf'];
  % ulrfld = [RAPFX '_ulrf'];
  % lrfld = [RAPFX '_lrf'];

  % Water-benthos fluxes
  qbfld = ['benthic_' RAPFX '_srf'];
  btfld = ['benthic_' RAPFX '_t'];
  qbofld = ['benthic_' RAPFX '_qbo'];

  % Turbulent and net surface fluxes
  TURPFX = ['ndbc_' RAPFX '_30a'];

  qlhfld = [TURPFX '_latent_heat_flux'];
  qshfld = [TURPFX '_sensible_heat_flux'];
  qrhfld = [TURPFX '_rain_heat_flux'];

  q0fld = [TURPFX '_net_heat_flux'];
  qtfld = [q0fld '_term'];

  % Ocean processes
  WAVEPFX = 'ww3';
  whfld = [WAVEPFX '_sigwavehgt'];
  wpfld = [WAVEPFX '_peakwaveper'];
  wdfld = [WAVEPFX '_peakwavedir'];

  KMPFX = 'fkeys_hycom';
  %KMPFX = 'gom_hycom';
  ufld = [KMPFX '_u'];
  vfld = [KMPFX '_v'];
  Tfld = [KMPFX '_seatemp_field'];

  netufld = [TIDEPFX '_' KMPFX '_u'];
  netvfld = [TIDEPFX '_' KMPFX '_v'];

  STOKESPFX = [WAVEPFX '_ndbc_stokes'];
  ssufld = [STOKESPFX '_u'];
  ssvfld = [STOKESPFX '_v'];
  sssfld = [STOKESPFX '_speed'];
  ssdfld = [STOKESPFX '_dir'];

  % Quasi-Eulerian currents and heat advection
  % Careful not to double-add advection by km-scale currents!
  switch (KMPFX),
   case 'fkeys_hycom',	QEPFX = [WAVEPFX '_fkeys_qe'];
   case 'gom_hycom',	QEPFX = [WAVEPFX '_gom_qe'];
   otherwise,		error('Unknown km-scale model "%s"',KMPFX);
  end;
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

  % Total budget with benthic flux
  qbotfld = [qbofld '_term'];
  bq0fld = ['benthic_' q0fld];
  bdTfld = ['benthic_' dTfld];


  stn = station_absorbed_insolation(stn,asrfld,srfld,hfld,[],[],gamfld,qbfld);

  if ( ~isfield(stn,lrfld) )
    stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);
  end;

  % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
  %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
  %                         Dfld,qeufld,qevfld,wpfld,whfld,true);
  [swix,lwix,lhix,shix,rhix] = ...
      intersect_all_dates([],stn.(asrfld).date,stn.(lrfld).date,...
                          stn.(qlhfld).date,stn.(qshfld).date,stn.(qrhfld).date);

  stn.(q0fld).date = stn.(qlhfld).date(lhix);
  stn.(q0fld).data = stn.(asrfld).data(swix) + stn.(lrfld).data(lwix) ...
      + stn.(qlhfld).data(lhix) + stn.(qshfld).data(shix) + stn.(qrhfld).data(rhix);

  %station_heat_flux_term(stn,nffld,htfld,tfld,sfld,dfld)
  stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],hfld);


  %%%
  %% Eulerian (km + Stokes) Heat Advection, Km-scale Heat Diffusion
  stn = station_calc_udotdelt(stn,qeufld,qevfld,Tfld,kmtfld,...
                              ['raw_' udTfld],udTfld,...
                              qtfld,qtAdvfld);
  stn = station_calc_kdel2t(stn,K_theta,Tfld,...
                            ['raw_' kd2Tfld],kd2Tfld,...
                            qtAdvfld,dTfld);


  stn = station_heat_flux_term_inverse(stn,[dTfld '_heat_flux'],...
                                       dTfld,sfld,[],hfld);


  %%%
  %% Benthic heat exchanges

  % Spline-fit an hourly time series of mean currents to native data
  x.meanu = interp_ts(stn.(ufld));
  x.meanv = interp_ts(stn.(vfld));
  stn.(netufld) = ts_op(stn.(tufld),x.meanu,'+');
  stn.(netvfld) = ts_op(stn.(tvfld),x.meanv,'+');

  dTdt.date = stn.(sfld).date(1:end-1);
  dTdt.data = diff(stn.(sfld).data);
  gapix = find(diff(stn.(sfld).date) > (2/24));
  dTdt.date(gapix) = [];
  dTdt.data(gapix) = [];

  global hbl
  global hbls
  % hbls = logspace(-2,0,5);
  hbls=0.01;
  % hbls=0.03;
  % hbls=1.00;
  global hb
  global hbs
  % hbs = 2:8;
  hbs=3;
  % hbs=10;
  global Cbd
  global Cbds
  Cbds = [];
  % Cbds = logspace(-5,-3,20);
  Cbds=3.8e-4;
  % Cbds=3e-4;
  % Cbds=0;
  global resid
  resid = [];
  resid = repmat(nan,[length(hbls) length(hbs) length(Cbds)]);
  global maxerr
  maxerr = [];
  maxerr = repmat(nan,[length(hbls) length(hbs) length(Cbds)]);
  for hblix = 1:length(hbls)
   hbl = hbls(hblix);
   disp(hbl);
   for hix = 1:length(hbs)
    hb = hbs(hix);
    for ix = 1:length(Cbds)
      Cbd = Cbds(ix);
      stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld);
      stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],hfld);

      stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'+');
      stn = station_heat_flux_term(stn,bq0fld,[bq0fld '_term'],sfld,[],hfld);

      stn.(bdTfld) = ts_op(stn.(dTfld),stn.(qbotfld),'+');
      stn = station_heat_flux_term_inverse(stn,[bdTfld '_heat_flux'],bdTfld,sfld,[],hfld);

      resid(hblix,hix,ix) = error_ts(stn.(bdTfld),dTdt,'rmse');
      maxerr(hblix,hix,ix) = nanmax(cumsum(stn.(bdTfld).data));
    end;
   end;
  end;

  % stn = compare_heat_budgets(stn,sfld,RAPFX,KMPFX); appendtitlename(' (with benthic)');


  % Horizontal convection!
  R = (1.00-0.08);
  global cfacs
  cfacs = 0:0.1:4;
cfacs = 0.2;
  % cfac = 1.2;
  global wfacs
  wfacs = 0:0.1:4;
wfacs = 2.2;
  % wfac = 1.2;
  lagoff = 0;
  resid = [];
  resid = repmat(nan,[length(cfacs) length(wfacs)]);
  maxerr = [];
  maxerr = repmat(nan,[length(cfacs) length(wfacs)]);
  for cix = 1:length(cfacs)
   cfac = cfacs(cix);
   for wix = 1:length(wfacs)
    wfac = wfacs(wix);
    disp([cfac,wfac]);
    stn = station_horizontal_convection(stn,sfld,[],hfld,...
                                        [bdTfld '_heat_flux_24_hour_average'],...
                                        [bdTfld '_qvf'],[bdTfld '_qf'],R,cfac,wfac,...
                                        bdTfld,[bdTfld '_netqf'],lagoff);
    % stn = station_heat_flux_term_inverse(stn,[bdTfld '_netqf_heat_flux'],...
    %                                      [bdTfld '_netqf'],...
    %                                      sfld,[],hfld);
    goodix = find(~isnan(stn.([bdTfld '_netqf']).data));
    resid(cix,wix) = error_ts(stn.([bdTfld '_netqf']),dTdt,'rmse',goodix);
    maxerr(cix,wix) = max(abs(cumsum( stn.([bdTfld '_netqf']).data(goodix) )));
   end;
  end;

  % stn = compare_heat_budgets(stn,sfld,RAPFX,KMPFX);
  % appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
  % if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;

  % ms1_ann2;

  set_more;

return;
