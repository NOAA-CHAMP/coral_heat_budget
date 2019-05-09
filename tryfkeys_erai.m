function stn = tryfkeys_erai(stn,R,cfac,wfac,lagoff,doPlot,doAverages)
%function stn = tryfkeys_erai(stn,R,cfac,wfac,lagoff,doPlot,doAverages)
%
% Try various combinations of heat budget terms using FKEYS HYCOM 1km data
%
% Last Saved Time-stamp: <Tue 2011-04-19 10:38:20  Lew.Gramer>

  set_more off

  datapath = get_thesis_path('../data');

  if ( ischar(stn) )
    stnm = lower(stn);
    matfname = fullfile(datapath,[stnm '_ms.mat']);
    disp(['Loading ' matfname]);
    x = load(matfname,'station');
    stn = x.station;
    x = []; clear x;
  end;

  if ( ~exist('R','var') || isempty(R) )
    R = (1.00-0.08);
  end;
  if ( ~exist('cfac','var') || isempty(cfac) )
    cfac = 0.3;
  end;
  if ( ~exist('wfac','var') || isempty(wfac) )
    wfac = 0.3;
  end;
  if ( ~exist('lagoff','var') || isempty(lagoff) )
    lagoff = 0;
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = true;
  end;
  if ( ~exist('doAverages','var') || isempty(doAverages) )
    doAverages = false;
  end;

  if ( ~exist('K_theta','var') || isempty(K_theta) )
    % FKEYS equivalent Fickian heat dispersion
    K_theta = 2.5;
  end;

  % Load heat budget monthly/daily climatologies for comparison
  if ( doAverages )
    if ( ~isfield(stn,'monthly_nocs_srf') )
      stn = annocs(stn);
    end;
    if ( ~isfield(stn,'landy_sr') )
      stn = station_load_landy(stn);
    end;
    if ( ~isfield(stn,'daily_oaflux_net_heat_flux') )
      stn = station_load_oaflux(stn);
    end;
  end;


  if ( ~isfield(stn,'fkeys_hycom_u') )
    stn = get_fkeys_hycom(stn);
    set_more off
  end;

  % How do we estimate water depth?
  hfld = 'tmd_tide_i_depth';


  %%%
  % First, recalculate HC(Q_0) result each time also
  %%%
  PFX = 'ndbc_erai_30a_';

  origsrfld = 'erai_srf';
  srfld = 'erai_srf_absorbed';
  lrfld = 'erai_lrf';

  if ( ~isfield(stn,origsrfld) )
    stn = get_erai_station(stn);
  end;

  stn = station_absorbed_insolation(stn,srfld,origsrfld,hfld);

  lffld = [PFX 'latent_heat_flux'];
  sffld = [PFX 'sensible_heat_flux'];
  rffld = [PFX 'rain_heat_flux'];

  hffld = [PFX 'absorbed_heat_flux'];
  ntfld = [PFX 'absorbed_heat_flux_term'];
  [swix,lwix,lfix,sfix] = intersect_all_dates([], ...
                                              stn.(srfld).date, ...
                                              stn.(lrfld).date, ...
                                              stn.(lffld).date, ...
                                              stn.(sffld).date ...
                                              );
  stn.(hffld).date = stn.(srfld).date(swix);
  stn.(hffld).data = stn.(srfld).data(swix) + ...
      stn.(lrfld).data(lwix) + ...
      stn.(lffld).data(lfix) + ...
      stn.(sffld).data(sfix);
  % Include precipitation flux if we got it
  if ( isfield(stn,rffld) )
    [rfix,hfix] = intersect_dates(stn.(rffld).date,stn.(hffld).date);
    stn.(hffld).date = stn.(hffld).date(hfix);
    stn.(hffld).data = stn.(hffld).data(hfix) + stn.(rffld).data(rfix);
  end;
  % Get rid of any annoying imaginary parts?!
  stn.(hffld).data = real(stn.(hffld).data);

  % Calculate heat budget term (Q0/Cp*rho*h)
  stn = station_heat_flux_term(stn,hffld,ntfld,'ndbc_sea_t',[],hfld);

  % Do simple QA/QC
  badix = find(abs(stn.(hffld).data) > 2000);
  if ( ~isempty(badix) )
    warning('Deleting %d bad points from %s',length(badix),hffld);
    stn.(hffld).date(badix) = [];
    stn.(hffld).data(badix) = [];
  end;

  badix = find(abs(stn.(ntfld).data) > 10);
  if ( ~isempty(badix) )
    warning('Deleting %d bad points from %s',length(badix),ntfld);
    stn.(ntfld).date(badix) = [];
    stn.(ntfld).data(badix) = [];
  end;

  if ( ~strcmp(srfld,origsrfld) )
    orighffld = [PFX 'net_heat_flux'];
    origntfld = [PFX 'heat_flux_term'];

    [swix,lwix,lfix,sfix] = intersect_all_dates([], ...
                                                stn.(origsrfld).date, ...
                                                stn.(lrfld).date, ...
                                                stn.(lffld).date, ...
                                                stn.(sffld).date ...
                                                );
    stn.(orighffld).date = stn.(origsrfld).date(swix);
    stn.(orighffld).data = stn.(origsrfld).data(swix) + ...
        stn.(lrfld).data(lwix) + ...
        stn.(lffld).data(lfix) + ...
        stn.(sffld).data(sfix);
    % Include precipitation flux if we got it
    if ( isfield(stn,rffld) )
      [rfix,hfix] = intersect_dates(stn.(rffld).date,stn.(orighffld).date);
      stn.(orighffld).date = stn.(orighffld).date(hfix);
      stn.(orighffld).data = stn.(orighffld).data(hfix) + stn.(rffld).data(rfix);
    end;

    % Get rid of any annoying imaginary parts?!
    stn.(orighffld).data = real(stn.(orighffld).data);
    % Calculate heat budget term (Q0/Cp*rho*h)
    stn = station_heat_flux_term(stn,orighffld,origntfld,'ndbc_sea_t',[],hfld);
    % Do simple QA/QC
    badix = find(abs(stn.(orighffld).data) > 2000);
    if ( ~isempty(badix) )
      stn.(orighffld).date(badix) = [];
      stn.(orighffld).data(badix) = [];
    end;
    badix = find(abs(stn.(origntfld).data) > 10);
    if ( ~isempty(badix) )
      stn.(origntfld).date(badix) = [];
      stn.(origntfld).data(badix) = [];
    end;
  end;


  hafld = [hffld '_24_hour_average'];
  hcfld = ['netqf'];
  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],hfld,...
                                      hafld,'qvf','qf',R,cfac,wfac,...
                                      ntfld,hcfld,lagoff);


  %%%
  % Then, recalculate various advection terms
  %%%
  if ( ~isfield(stn,'fkeys_hycom_dt') )
    stn = station_calc_udotdelt(stn,'fkeys_hycom_u','fkeys_hycom_v',...
                                'fkeys_hycom_seatemp_field',...
                                'fkeys_hycom_seatemp',...
                                'daily_fkeys_hycom_advected_heat',...
                                'fkeys_hycom_advected_heat',...
                                ntfld,'fkeys_hycom_dt');
  end;

  stn = station_calc_udotdelt(stn,'fkeys_hycom_u','fkeys_hycom_v',...
                              'fkeys_hycom_seatemp_field',...
                              'fkeys_hycom_seatemp',...
                              'daily_fkeys_hycom_advected_heat',...
                              'fkeys_hycom_advected_heat',...
                              hcfld,'fkeys_hycom_dt_netqf');

  stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_dt_heat_flux','fkeys_hycom_dt',...
                                       'ndbc_sea_t',[],hfld);

  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],hfld,...
                                      'fkeys_hycom_dt_heat_flux_24_hour_average',...
                                      'fkeys_hycom_qvf','fkeys_hycom_qf',R,cfac,wfac,...
                                      'fkeys_hycom_dt','fkeys_hycom_netqf',lagoff);

  if ( ~isfield(stn,'ndbc_sea_t_implied_heat_flux') )
    stn.ndbc_sea_t_diff.date = stn.ndbc_sea_t.date(1:end-1);
    stn.ndbc_sea_t_diff.data = diff(stn.ndbc_sea_t.data);
    badix = find(diff(stn.ndbc_sea_t_diff.date) >= (1.1/24.0));
    stn.ndbc_sea_t_diff.date(badix) = [];
    stn.ndbc_sea_t_diff.data(badix) = [];
    stn = station_heat_flux_term_inverse(stn,'ndbc_sea_t_implied_heat_flux','ndbc_sea_t_diff',...
                                         'ndbc_sea_t',[],hfld);
  end;

  stn = station_heat_flux_term_inverse(stn,'netqf_heat_flux','netqf',...
                                       'ndbc_sea_t',[],hfld);
  stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_netqf_heat_flux','fkeys_hycom_netqf',...
                                       'ndbc_sea_t',[],hfld);



  if ( ~isfield(stn,'ww3_sigwavehgt') )
    stn = get_ww3_station(stn);
  end;
  if ( ~isfield(stn,'ww3_stokes_speed') )
    stn = station_stokes_drift(stn,'ww3_stokes_speed','ww3_stokes_dir','ww3_stokes_u','ww3_stokes_v','ndbc_wind1_speed','ndbc_wind1_dir','ww3_sigwavehgt','ww3_peakwaveper','ww3_peakwavedir');
  end;
  if ( ~isfield(stn,'fkeys_hycom_quasi_eulerian_speed') )
    stn = calc_quasi_eulerian(stn,'ww3_stokes','fkeys_hycom','fkeys_hycom_quasi_eulerian');
  end;
  if ( ~isfield(stn,'fkeys_hycom_qedt') )
    stn = station_calc_udotdelt(stn,'fkeys_hycom_quasi_eulerian_u','fkeys_hycom_quasi_eulerian_v',...
                                'fkeys_hycom_seatemp_field',...
                                'fkeys_hycom_seatemp',...
                                'daily_fkeys_hycom_quasi_eulerian_advected_heat',...
                                'fkeys_hycom_quasi_eulerian_advected_heat',...
                                ntfld,'fkeys_hycom_qedt');
    stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_qedt_heat_flux','fkeys_hycom_qedt',...
                                         'ndbc_sea_t',[],hfld);
  end;

  stn = station_calc_udotdelt(stn,'fkeys_hycom_quasi_eulerian_u','fkeys_hycom_quasi_eulerian_v',...
                              'fkeys_hycom_seatemp_field',...
                              'fkeys_hycom_seatemp',...
                              'daily_fkeys_hycom_quasi_eulerian_advected_heat',...
                              'fkeys_hycom_quasi_eulerian_advected_heat',...
                              'netqf','fkeys_hycom_qedt_netqf');

  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],hfld,...
                                      'fkeys_hycom_qedt_heat_flux_24_hour_average',...
                                      'fkeys_hycom_qeqvf','fkeys_hycom_qeqf',R,cfac,wfac,...
                                      'fkeys_hycom_qedt','fkeys_hycom_qenetqf',lagoff);


  % Now play with adding Laplacian heat diffusion term KH*del2(T)

  stn = station_calc_kdel2t(stn,K_theta,'fkeys_hycom_seatemp_field',...
                            'native_fkeys_hycom_diffused_heat',...
                            'fkeys_hycom_diffused_heat',...
                            'fkeys_hycom_qedt','fkeys_hycom_qelt');
  stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_qelt_heat_flux',...
                                       'fkeys_hycom_qelt',...
                                       'ndbc_sea_t',[],hfld);

  stn = station_calc_kdel2t(stn,K_theta,'fkeys_hycom_seatemp_field',...
                            'native_fkeys_hycom_diffused_heat',...
                            'fkeys_hycom_diffused_heat',...
                            'fkeys_hycom_qedt_netqf','fkeys_hycom_qelt_netqf');
  stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_qelt_netqf_heat_flux',...
                                       'fkeys_hycom_qelt_netqf',...
                                       'ndbc_sea_t',[],hfld);

  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],hfld,...
                                      'fkeys_hycom_qelt_heat_flux_24_hour_average',...
                                      'fkeys_hycom_qelqvf','fkeys_hycom_qelqf',R,cfac,wfac,...
                                      'fkeys_hycom_qelt','fkeys_hycom_qelnetqf',lagoff);
  stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_qelnetqf_heat_flux',...
                                       'fkeys_hycom_qelnetqf',...
                                       'ndbc_sea_t',[],hfld);



  disp('Adjusting K_H for budget residual...');
  stn = station_calc_sgs_diffusion(stn,'ndbc_sea_t', ...
                                   'fkeys_hycom_qelnetqf', ...
                                   'fkeys_hycom_seatemp_field', ...
                                   'fkeys_hycom_seatemp_laplacian', ...
                                   'fkeys_hycom_sgs_diffused_heat', ...
                                   'fkeys_hycom_sgs_thermal_diffusivity', ...
                                   'fkeys_hycom_absorbed_nn30a_final_budget');


  % Calculate daily and monthly averages for comparisons
  if ( doAverages )

    disp('Calculating daily means');
    dlyhffld = ['daily_' hffld]; dlycihffld = ['daily_ci_' hffld];
    dlydts = floor(stn.(hffld).date);
    stn.(dlyhffld).date = unique(dlydts);
    stn.(dlycihffld).date = unique(dlydts);
    [stn.(dlyhffld).data,stn.(dlycihffld).data] = ...
        grpstats(real(stn.(hffld).data),dlydts,{'mean','meanci'});

    dlydts = floor(stn.fkeys_hycom_qelt_heat_flux.date);
    stn.daily_fkeys_hycom_qelt_heat_flux.date = unique(dlydts);
    stn.daily_ci_fkeys_hycom_qelt_heat_flux.date = unique(dlydts);
    [stn.daily_fkeys_hycom_qelt_heat_flux.data,...
     stn.daily_ci_fkeys_hycom_qelt_heat_flux.data] = ...
        grpstats(real(stn.fkeys_hycom_qelt_heat_flux.data),dlydts,{'mean','meanci'});

    dlydts = floor(stn.fkeys_hycom_qelnetqf.date);
    stn.daily_fkeys_hycom_qelnetqf.date = unique(dlydts);
    stn.daily_ci_fkeys_hycom_qelnetqf.date = unique(dlydts);
    [stn.daily_fkeys_hycom_qelnetqf.data,...
     stn.daily_ci_fkeys_hycom_qelnetqf.data] = ...
        grpstats(real(stn.fkeys_hycom_qelnetqf.data),dlydts,{'mean','meanci'});

    dlydts = floor(stn.fkeys_hycom_qelnetqf_heat_flux.date);
    stn.daily_fkeys_hycom_qelnetqf_heat_flux.date = unique(dlydts);
    stn.daily_ci_fkeys_hycom_qelnetqf_heat_flux.date = unique(dlydts);
    [stn.daily_fkeys_hycom_qelnetqf_heat_flux.data,...
     stn.daily_ci_fkeys_hycom_qelnetqf_heat_flux.data] = ...
        grpstats(real(stn.fkeys_hycom_qelnetqf_heat_flux.data),dlydts,{'mean','meanci'});


    disp('Calculating monthly means');
    mlyhffld = ['monthly_' hffld]; mlycihffld = ['monthly_ci_' hffld];
    [yr,mo,dy] = datevec(stn.(hffld).date); mlydts = datenum(yr,mo,1);
    stn.(mlyhffld).date = unique(mlydts);
    stn.(mlycihffld).date = unique(mlydts);
    [stn.(mlyhffld).data,stn.(mlycihffld).data] = ...
        grpstats(real(stn.(hffld).data),mlydts,{'mean','meanci'});

    [yr,mo,dy] = datevec(stn.fkeys_hycom_qelt_heat_flux.date); mlydts = datenum(yr,mo,1);
    stn.monthly_fkeys_hycom_qelt_heat_flux.date = unique(mlydts);
    stn.monthly_ci_fkeys_hycom_qelt_heat_flux.date = unique(mlydts);
    [stn.monthly_fkeys_hycom_qelt_heat_flux.data,stn.monthly_ci_fkeys_hycom_qelt_heat_flux.data] = ...
        grpstats(real(stn.fkeys_hycom_qelt_heat_flux.data),mlydts,{'mean','meanci'});

    [yr,mo,dy] = datevec(stn.fkeys_hycom_qelnetqf.date); mlydts = datenum(yr,mo,1);
    stn.monthly_fkeys_hycom_qelnetqf.date = unique(mlydts);
    stn.monthly_ci_fkeys_hycom_qelnetqf.date = unique(mlydts);
    [stn.monthly_fkeys_hycom_qelnetqf.data,stn.monthly_ci_fkeys_hycom_qelnetqf.data] = ...
        grpstats(real(stn.fkeys_hycom_qelnetqf.data),mlydts,{'mean','meanci'});

    [yr,mo,dy] = datevec(stn.fkeys_hycom_qelnetqf_heat_flux.date); mlydts = datenum(yr,mo,1);
    stn.monthly_fkeys_hycom_qelnetqf_heat_flux.date = unique(mlydts);
    stn.monthly_ci_fkeys_hycom_qelnetqf_heat_flux.date = unique(mlydts);
    [stn.monthly_fkeys_hycom_qelnetqf_heat_flux.data,stn.monthly_ci_fkeys_hycom_qelnetqf_heat_flux.data] = ...
        grpstats(real(stn.fkeys_hycom_qelnetqf_heat_flux.data),mlydts,{'mean','meanci'});

  end; %if doAverages


  disp('Calculating GoM HYCOM results also...');
  stn = trygom(stn,R,cfac,wfac,lagoff,false,hffld,ntfld);
  set_more off


  % Do various interesting plots
  if ( doPlot )

    [firstyr,ig,ig] = datevec(stn.fkeys_hycom_dt.date(1));
    [lastyr,ig,ig] = datevec(stn.fkeys_hycom_dt.date(end));
    dys = datenum(lastyr,12,31) - datenum(firstyr,1,1) + 1;

    if (0)
      absflds = { ...
          'ndbc_erai_30a_wind_stress_30_day_maximum', ...
          'fkeys_hycom_speed_7_day_maximum', ...
                };
      accflds = { ...
          'fkeys_hycom_netqf', ...
          'fkeys_hycom_dt_netqf', ...
          'netqf', ...
                };
      fh = plot_fluxes(stn,firstyr,1,dys,[],absflds,accflds);
      appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
      if ( lagoff ~= 0 )
        appendtitlename(sprintf(' lag:%d', lagoff));
      end;

      plot_fluxes(stn,firstyr,1,dys,{'ndbc_sea_t'},[],{'ndbc_erai_30a_heat_flux_term','gom_hycom_dt','fkeys_hycom_dt','fkeys_hycom_qedt'},[],{'NDBC sea temperature','ERA-Interim/TOGA-COARE Q_0','GoM 4km HYCOM + Q_0','FKEYS 1km HYCOM + Q_0','FKEYS 1km HYCOM + Stokes + Q_0'});
      appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
      if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;
      % print('-dtiff','../figs/mlrf1-without-siphon.tiff');

      plot_fluxes(stn,firstyr,1,dys,{'ndbc_sea_t'},[],{'netqf','gom_hycom_dt_netqf','fkeys_hycom_dt_netqf','fkeys_hycom_netqf'},[],{'NDBC sea temperature','HC(Q_0) (thermal siphon)','GoM 4km HYCOM + HC(Q_0)','FKEYS 1km HYCOM + HC(Q_0)','HC(FKEYS 1km HYCOM + Q_0)'});
      appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
      if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;
      % print('-dtiff','../figs/mlrf1-with-siphon-2008.tiff');

      plot_fluxes(stn,firstyr,1,dys,{'ndbc_sea_t'},[],{'netqf','fkeys_hycom_dt','fkeys_hycom_netqf','fkeys_hycom_qedt','fkeys_hycom_qenetqf'},[],{'NDBC sea temperature','HC(Q_0) (thermal siphon)','FKEYS 1km HYCOM + Q_0','HC(FKEYS 1km HYCOM + Q_0)','FKEYS + Stokes + Q_0','HC(FKEYS + Stokes + Q_0)'});
      appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
      if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;
      % print('-dtiff','../figs/mlrf1-with-siphon-all-years.tiff');

      plot_fluxes(stn,firstyr,1,dys,{'ndbc_sea_t'},[],{'netqf','fkeys_hycom_netqf','fkeys_hycom_qenetqf'},[],{'NDBC sea temperature','HC(Q_0) (thermal siphon)','HC(FKEYS 1km HYCOM + Q_0)','HC(FKEYS 1km HYCOM + Stokes + Q_0)'});
      appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
      if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;
      % print('-dtiff','../figs/mlrf1-just-siphon-all-years.tiff');
    end; %if 0


    plot_fluxes(stn,firstyr,1,dys, ...
                {'ndbc_sea_t'},[],...
                {...
                    'netqf',...
                    ...
                    'gom_hycom_qenetqf',...
                    'fkeys_hycom_qenetqf',...
                    ...
                    'gom_hycom_qelnetqf',...
                    'fkeys_hycom_qelnetqf',...
                    ...
                    'fkeys_hycom_absorbed_nn30a_final_budget',...
                },[],...
                {...
                    'NDBC sea temperature',...
                    ...
                    'HC(Q_0) (thermal siphon)',...
                    ...
                    'HC(GoM 4km HYCOM u^.\nablaT + Stokes + Q_0)',...
                    'HC(FKEYS 1km HYCOM u^.\nablaT + Stokes + Q_0)',...
                    ...
                    'HC(GoM 4km HYCOM u^.\nablaT + \nabla^2T + Stokes + Q_0)',...
                    'HC(FKEYS 1km HYCOM u^.\nablaT + \nabla^2T + Stokes + Q_0)',...
                    ...
                    'HC(FKEYS 1km HYCOM u^.\nablaT + \nabla^2T + Stokes + Q_0)+\nabla^2T_S_G_S',...
                });
    appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
    if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;


    plot_fluxes(stn,firstyr,1,dys, ...
                {'ndbc_sea_t'},[],...
                {...
                    'ndbc_erai_30a_heat_flux_term',...
                    'ndbc_erai_30a_absorbed_heat_flux_term',...
                    ...
                    'gom_hycom_qedt_netqf',...
                    'fkeys_hycom_qedt_netqf',...
                    ...
                    'gom_hycom_qelt_netqf',...
                    'fkeys_hycom_qelt_netqf',...
                },[],...
                {...
                    'NDBC sea temperature',...
                    ...
                    'ERA-I/TOGA-COARE 3.0a Q_0',...
                    'ERA-I/TOGA-COARE 3.0a Q_0 w/\gamma',...
                    ...
                    'GoM 4km HYCOM u^.\nablaT + Stokes + HC(Q_0)',...
                    'FKEYS 1km HYCOM u^.\nablaT + Stokes + HC(Q_0)',...
                    ...
                    'GoM 4km HYCOM u^.\nablaT + \nabla^2T + Stokes + HC(Q_0)',...
                    'FKEYS 1km HYCOM u^.\nablaT + \nabla^2T + Stokes + HC(Q_0)',...
                });
    appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
    if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;


    if ( doAverages )
      % Plot successively more refined daily estimates vs. daily climatology
      figure;
      maxigraph;
      hold on;
      plot(stn.(dlyhffld).date,stn.(dlyhffld).data,'r.-');
      plot(stn.daily_fkeys_hycom_qelt_heat_flux.date,stn.daily_fkeys_hycom_qelt_heat_flux.data,'bo-');
      plot(stn.daily_fkeys_hycom_qelnetqf_heat_flux.date,stn.daily_fkeys_hycom_qelnetqf_heat_flux.data,'g^-');
      plot(stn.daily_oaflux_net_heat_flux.date,stn.daily_oaflux_net_heat_flux.data,'k-');
      xlim([stn.daily_fkeys_hycom_qelt_heat_flux.date(1),stn.daily_oaflux_net_heat_flux.date(end)]);
      datetick3;
      legend('TOGA/COARE 3.0a Q_0','FKEYS 1km HYCOM+Stokes+Q_0','HC(FKEYS 1km HYCOM+Stokes+Q_0)','OAFlux Q_0', 'Location','Best');
      titlename('Daily heat fluxes [W/m^2]');

      % Plot most refined estimates with daily 95% error bars vs. climatology +- turbulent errors
      lh=[];
      figure;
      maxigraph;
      hold on;
      lh(end+1)=plot(stn.(dlyhffld).date,stn.(dlyhffld).data,'r.-');
      plot(stn.(dlycihffld).date,stn.(dlycihffld).data,'r+');
      lh(end+1)=plot(stn.daily_fkeys_hycom_qelnetqf_heat_flux.date,stn.daily_fkeys_hycom_qelnetqf_heat_flux.data,'g^-');
      plot(stn.daily_ci_fkeys_hycom_qelnetqf_heat_flux.date,stn.daily_ci_fkeys_hycom_qelnetqf_heat_flux.data,'g+');
      lh(end+1)=plot(stn.daily_oaflux_net_heat_flux.date,stn.daily_oaflux_net_heat_flux.data,'k-');
      [ix1,ix2]=intersect_dates(stn.daily_oaflux_net_heat_flux.date,stn.daily_oaflux_latent_flux_err.date);
      plot(stn.daily_oaflux_net_heat_flux.date(ix1),[stn.daily_oaflux_net_heat_flux.data(ix1)-stn.daily_oaflux_sensible_flux_err.data(ix2)-stn.daily_oaflux_latent_flux_err.data(ix2),stn.daily_oaflux_net_heat_flux.data(ix1)+stn.daily_oaflux_sensible_flux_err.data(ix2)+stn.daily_oaflux_latent_flux_err.data(ix2)],'k+');
      xlim([stn.daily_fkeys_hycom_qelnetqf_heat_flux.date(1),stn.daily_oaflux_net_heat_flux.date(end)]);
      datetick3;
      legend(lh,'Q_O \pm95%','HC(FKEYS 1km HYCOM+Stokes+Q_0) \pm95%','OAFlux Q_0 \pmerr', 'Location','Best');
      titlename('Daily heat fluxes [W/m^2]');


      % Plot successively more refined monthly estimates vs. monthly climatologies
      lh=[];
      figure;
      maxigraph;
      hold on;
      lh(end+1)=plot(stn.(mlyhffld).date,stn.(mlyhffld).data,'r.-');
      lh(end+1)=plot(stn.monthly_fkeys_hycom_qelt_heat_flux.date,stn.monthly_fkeys_hycom_qelt_heat_flux.data,'bo-');
      lh(end+1)=plot(stn.monthly_fkeys_hycom_qelnetqf_heat_flux.date,stn.monthly_fkeys_hycom_qelnetqf_heat_flux.data,'g^-');
      plot(stn.monthly_ci_fkeys_hycom_qelnetqf_heat_flux.date,stn.monthly_ci_fkeys_hycom_qelnetqf_heat_flux.data,'g+');
      lh(end+1)=plot(stn.landy_net_heat_flux.date,stn.landy_net_heat_flux.data,'k:', 'Color',[.5,.5,.5]);
      lh(end+1)=plot(stn.monthly_nocs_net_heat_flux.date,stn.monthly_nocs_net_heat_flux.data,'k-');
      xlim(stn.monthly_fkeys_hycom_qelt_heat_flux.date([1 end]));
      datetick3;
      legend(lh,'TOGA/COARE 3.0a Q_0','FKEYS 1km HYCOM+Stokes+Q_0','HC(FKEYS 1km HYCOM+Stokes+Q_0)\pm95%','Large&Yeager Q_0','NOC Southampton Q_0', 'Location','Best');
      titlename('Monthly heat fluxes [W/m^2]');
    end; %if doAverages


  end; %if doPlot


  set_more;

return;
