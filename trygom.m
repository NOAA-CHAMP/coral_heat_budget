function stn = trygom(stn,R,cfac,wfac,lagoff,doPlot,nffld,htfld)
%function stn = trygom(stn,R,cfac,wfac,lagoff,doPlot,nffld,htfld)
%
% Try various combinations of heat budget terms using Gulf of Mexico HYCOM 4km data
%
% Last Saved Time-stamp: <Tue 2011-04-19 10:37:20  Lew.Gramer>

  datapath = get_thesis_path('../data');

  if ( ischar(stn) )
    stnm = stn;
    matfname = fullfile(datapath,sprintf('%s_ms.mat',stnm));
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
  if ( ~exist('nffld','var') || isempty(nffld) )
    nffld = 'ndbc_ncep_30a_net_heat_flux';
  end;
  if ( ~exist('htfld','var') || isempty(htfld) )
    htfld = 'ndbc_ncep_30a_heat_flux_term';
  end;


  % Load heat budget monthly climatologies for comparison
  if ( ~isfield(stn,'monthly_nocs_srf') )
    stn = annocs(stn);
  end;
  if ( ~isfield(stn,'landy_sr') )
    stn = station_load_landy(stn);
  end;


  if ( ~isfield(stn,'gom_hycom_u') )
    stn = get_gom_hycom(stn);
    set_more off
  end;

  badix = find(abs(stn.(nffld).data) > 2000);
  if ( ~isempty(badix) )
    warning('Deleting %d bad points from %s',length(badix),nffld);
    stn.(nffld).date(badix) = [];
    stn.(nffld).data(badix) = [];
  end;

  badix = find(abs(stn.(htfld).data) > 10);
  if ( ~isempty(badix) )
    warning('Deleting %d bad points from %s',length(badix),htfld);
    stn.(htfld).date(badix) = [];
    stn.(htfld).data(badix) = [];
  end;

  if ( ~isfield(stn,'gom_hycom_dt') )
    stn = station_calc_udotdelt(stn,'gom_hycom_u','gom_hycom_v',...
                                'gom_hycom_seatemp_field',...
                                'gom_hycom_seatemp',...
                                'daily_gom_hycom_advected_heat',...
                                'gom_hycom_advected_heat',...
                                htfld,'gom_hycom_dt');
  end;


  stn = station_calc_udotdelt(stn,'gom_hycom_u','gom_hycom_v',...
                              'gom_hycom_seatemp_field',...
                              'gom_hycom_seatemp',...
                              'daily_gom_hycom_advected_heat',...
                              'gom_hycom_advected_heat',...
                              'netqf','gom_hycom_dt_netqf');


  if ( ~isfield(stn,'gom_hycom_dt_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'gom_hycom_dt_heat_flux','gom_hycom_dt',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;

  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],'tmd_tide_i_depth',...
                                      'gom_hycom_dt_heat_flux_24_hour_average',...
                                      'gom_hycom_qvf','gom_hycom_qf',R,cfac,wfac,...
                                      'gom_hycom_dt','gom_hycom_netqf',lagoff);

  if ( ~isfield(stn,'ndbc_sea_t_implied_heat_flux') )
    stn.ndbc_sea_t_diff.date = stn.ndbc_sea_t.date(1:end-1);
    stn.ndbc_sea_t_diff.data = diff(stn.ndbc_sea_t.data);
    badix = find(diff(stn.ndbc_sea_t_diff.date) >= (1.1/24.0));
    stn.ndbc_sea_t_diff.date(badix) = [];
    stn.ndbc_sea_t_diff.data(badix) = [];
    stn = station_heat_flux_term_inverse(stn,'ndbc_sea_t_implied_heat_flux','ndbc_sea_t_diff',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;

  stn = station_heat_flux_term_inverse(stn,'netqf_heat_flux','netqf',...
                                       'ndbc_sea_t',[],'tmd_tide_i_depth');

  stn = station_heat_flux_term_inverse(stn,'gom_hycom_netqf_heat_flux','gom_hycom_netqf',...
                                       'ndbc_sea_t',[],'tmd_tide_i_depth');



  if ( ~isfield(stn,'ww3_sigwavehgt') )
    stn = get_ww3_station(stn);
  end;
  if ( ~isfield(stn,'ww3_stokes_speed') )
    stn = station_stokes_drift(stn,'ww3_stokes_speed','ww3_stokes_dir','ww3_stokes_u','ww3_stokes_v','ndbc_wind1_speed','ndbc_wind1_dir','ww3_sigwavehgt','ww3_peakwaveper','ww3_peakwavedir');
  end;
  if ( ~isfield(stn,'gom_hycom_quasi_eulerian_speed') )
    stn = calc_quasi_eulerian(stn,'ww3_stokes','gom_hycom','gom_hycom_quasi_eulerian');
  end;
  if ( ~isfield(stn,'gom_hycom_qedt') )
    stn = station_calc_udotdelt(stn,'gom_hycom_quasi_eulerian_u','gom_hycom_quasi_eulerian_v',...
                                'gom_hycom_seatemp_field',...
                                'gom_hycom_seatemp',...
                                'daily_gom_hycom_quasi_eulerian_advected_heat',...
                                'gom_hycom_quasi_eulerian_advected_heat',...
                                htfld,'gom_hycom_qedt');
  end;

  stn = station_calc_udotdelt(stn,'gom_hycom_quasi_eulerian_u','gom_hycom_quasi_eulerian_v',...
                              'gom_hycom_seatemp_field',...
                              'gom_hycom_seatemp',...
                              'daily_gom_hycom_quasi_eulerian_advected_heat',...
                              'gom_hycom_quasi_eulerian_advected_heat',...
                              'netqf','gom_hycom_qedt_netqf');

  if ( ~isfield(stn,'gom_hycom_qedt_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'gom_hycom_qedt_heat_flux','gom_hycom_qedt',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;

  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],'tmd_tide_i_depth',...
                                      'gom_hycom_qedt_heat_flux_24_hour_average',...
                                      'gom_hycom_qeqvf','gom_hycom_qeqf',R,cfac,wfac,...
                                      'gom_hycom_qedt','gom_hycom_qenetqf',lagoff);


  % Now play with adding heat diffusion term
  if ( ~isfield(stn,'gom_hycom_qelt_heat_flux') )
    stn = station_calc_kdel2t(stn,20,'gom_hycom_seatemp_field',...
                              'daily_gom_hycom_diffused_heat','gom_hycom_diffused_heat',...
                              'gom_hycom_qedt','gom_hycom_qelt');
    stn = station_heat_flux_term_inverse(stn,'gom_hycom_qelt_heat_flux','gom_hycom_qelt',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;

  stn = station_calc_kdel2t(stn,20,'gom_hycom_seatemp_field',...
                            'daily_gom_hycom_diffused_heat','gom_hycom_diffused_heat',...
                            'gom_hycom_qedt_netqf','gom_hycom_qelt_netqf');
  stn = station_heat_flux_term_inverse(stn,'gom_hycom_qelt_netqf_heat_flux',...
                                       'gom_hycom_qelt_netqf',...
                                       'ndbc_sea_t',[],'tmd_tide_i_depth');

  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],'tmd_tide_i_depth',...
                                      'gom_hycom_qelt_heat_flux_24_hour_average',...
                                      'gom_hycom_qelqvf','gom_hycom_qelqf',R,cfac,wfac,...
                                      'gom_hycom_qelt','gom_hycom_qelnetqf',lagoff);


  % Do various interesting plots
  if ( doPlot )

    [firstyr,ig,ig] = datevec(stn.gom_hycom_dt.date(1));
    dys = floor(now) - datenum(firstyr,1,1) + 1;

    plot_fluxes(stn,2008,1,366,{'ndbc_sea_t'},[],{htfld,'gom_hycom_dt','flkeys_hycom_dt','flkeys_hycom_qedt'},[],{'NDBC sea temperature','NCEP NARR/TOGA-COARE Q_0','GoM 4km HYCOM + Q_0','FlKeys 1km HYCOM + Q_0','FlKeys 1km HYCOM + Stokes + Q_0'});
    appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
    if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;
    % print('-dtiff','../figs/mlrf1-without-siphon.tiff');

    plot_fluxes(stn,2008,1,366,{'ndbc_sea_t'},[],{'netqf','gom_hycom_dt_netqf','flkeys_hycom_dt_netqf','flkeys_hycom_netqf'},[],{'NDBC sea temperature','HC(Q_0) (thermal siphon)','GoM 4km HYCOM + HC(Q_0)','FlKeys 1km HYCOM + HC(Q_0)','HC(FlKeys 1km HYCOM + Q_0)'});
    appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
    if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;
    % print('-dtiff','../figs/mlrf1-with-siphon-2008.tiff');

    plot_fluxes(stn,firstyr,1,dys,{'ndbc_sea_t'},[],{'netqf','gom_hycom_dt','gom_hycom_netqf','gom_hycom_qedt','gom_hycom_qenetqf'},[],{'NDBC sea temperature','HC(Q_0) (thermal siphon)','GoM 4km HYCOM + Q_0','HC(GoM 4km HYCOM + Q_0)','GoM + Stokes + Q_0','HC(GoM + Stokes + Q_0)'});
    appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
    if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;
    % print('-dtiff','../figs/mlrf1-with-siphon-all-years.tiff');

    plot_fluxes(stn,firstyr,1,dys,{'ndbc_sea_t'},[],{'netqf','gom_hycom_netqf','gom_hycom_qenetqf'},[],{'NDBC sea temperature','HC(Q_0) (thermal siphon)','HC(GoM 4km HYCOM + Q_0)','HC(GoM 4km HYCOM + Stokes + Q_0)'});
    appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
    if (lagoff); appendtitlename(sprintf(' lag:%d', lagoff)); end;
    % print('-dtiff','../figs/mlrf1-just-siphon-all-years.tiff');

  end;

return;
