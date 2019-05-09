function stn = tryflkeys(stn,R,cfac,wfac,lagoff)
%function stn = tryflkeys(stn,R,cfac,wfac,lagoff)

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

  % Load heat budget monthly climatologies for comparison
  if ( ~isfield(stn,'monthly_nocs_srf') )
    stn = annocs(stn);
  end;
  if ( ~isfield(stn,'landy_sr') )
    stn = station_load_landy(stn);
  end;


  if ( ~isfield(stn,'flkeys_hycom_u') )
    stn = get_flkeys_hycom(stn);
  end;

  nffld = 'ndbc_ncep_30a_net_heat_flux';
  badix = find(abs(stn.(nffld).data) > 2000);
  if ( ~isempty(badix) )
    warning('Deleting %d bad points from %s',length(badix),nffld);
    stn.(nffld).date(badix) = [];
    stn.(nffld).data(badix) = [];
  end;

  htfld = 'ndbc_ncep_30a_heat_flux_term';
  badix = find(abs(stn.(htfld).data) > 10);
  if ( ~isempty(badix) )
    warning('Deleting %d bad points from %s',length(badix),htfld);
    stn.(htfld).date(badix) = [];
    stn.(htfld).data(badix) = [];
  end;

  if ( ~isfield(stn,'flkeys_hycom_dt') )
    stn = station_calc_udotdelt(stn,'flkeys_hycom_u','flkeys_hycom_v',...
                                'flkeys_hycom_seatemp_field',...
                                'flkeys_hycom_seatemp',...
                                'daily_flkeys_hycom_advected_heat',...
                                'flkeys_hycom_advected_heat',...
                                htfld,'flkeys_hycom_dt');
  end;

  if ( ~isfield(stn,'flkeys_hycom_dt_netqf') )
    stn = station_calc_udotdelt(stn,'flkeys_hycom_u','flkeys_hycom_v',...
                                'flkeys_hycom_seatemp_field',...
                                'flkeys_hycom_seatemp',...
                                'daily_flkeys_hycom_advected_heat',...
                                'flkeys_hycom_advected_heat',...
                                'netqf','flkeys_hycom_dt_netqf');
  end;

  if ( ~isfield(stn,'flkeys_hycom_dt_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'flkeys_hycom_dt_heat_flux','flkeys_hycom_dt',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;

  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],'tmd_tide_i_depth',...
                                      'flkeys_hycom_dt_heat_flux_24_hour_average',...
                                      'flkeys_hycom_qvf','flkeys_hycom_qf',R,cfac,wfac,...
                                      'flkeys_hycom_dt','flkeys_hycom_netqf',lagoff);

  if ( ~isfield(stn,'ndbc_sea_t_implied_heat_flux') )
    stn.ndbc_sea_t_diff.date = stn.ndbc_sea_t.date(1:end-1);
    stn.ndbc_sea_t_diff.data = diff(stn.ndbc_sea_t.data);
    badix = find(diff(stn.ndbc_sea_t_diff.date) >= (1.1/24.0));
    stn.ndbc_sea_t_diff.date(badix) = [];
    stn.ndbc_sea_t_diff.data(badix) = [];
    stn = station_heat_flux_term_inverse(stn,'ndbc_sea_t_implied_heat_flux','ndbc_sea_t_diff',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;

  if ( ~isfield(stn,'netqf_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'netqf_heat_flux','netqf',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;
  if ( ~isfield(stn,'flkeys_hycom_netqf_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'flkeys_hycom_netqf_heat_flux','flkeys_hycom_netqf',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;



  if ( ~isfield(stn,'ww3_sigwavehgt') )
    stn = get_ww3_station(stn);
  end;
  if ( ~isfield(stn,'ww3_stokes_speed') )
    stn = station_stokes_drift(stn,'ww3_stokes_speed','ww3_stokes_dir','ww3_stokes_u','ww3_stokes_v','ndbc_wind1_speed','ndbc_wind1_dir','ww3_sigwavehgt','ww3_peakwaveper','ww3_peakwavedir');
  end;
  if ( ~isfield(stn,'flkeys_hycom_quasi_eulerian_speed') )
    stn = calc_quasi_eulerian(stn,'ww3_stokes','flkeys_hycom','flkeys_hycom_quasi_eulerian');
  end;
  if ( ~isfield(stn,'flkeys_hycom_qedt') )
    stn = station_calc_udotdelt(stn,'flkeys_hycom_quasi_eulerian_u','flkeys_hycom_quasi_eulerian_v',...
                                'flkeys_hycom_seatemp_field',...
                                'flkeys_hycom_seatemp',...
                                'daily_flkeys_hycom_quasi_eulerian_advected_heat',...
                                'flkeys_hycom_quasi_eulerian_advected_heat',...
                                htfld,'flkeys_hycom_qedt');
  end;
  if ( ~isfield(stn,'flkeys_hycom_qedt_netqf') )
    stn = station_calc_udotdelt(stn,'flkeys_hycom_quasi_eulerian_u','flkeys_hycom_quasi_eulerian_v',...
                                'flkeys_hycom_seatemp_field',...
                                'flkeys_hycom_seatemp',...
                                'daily_flkeys_hycom_quasi_eulerian_advected_heat',...
                                'flkeys_hycom_quasi_eulerian_advected_heat',...
                                'netqf','flkeys_hycom_qedt_netqf');
  end;

  if ( ~isfield(stn,'flkeys_hycom_qedt_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'flkeys_hycom_qedt_heat_flux','flkeys_hycom_qedt',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;

  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],'tmd_tide_i_depth',...
                                      'flkeys_hycom_qedt_heat_flux_24_hour_average',...
                                      'flkeys_hycom_qeqvf','flkeys_hycom_qeqf',R,cfac,wfac,...
                                      'flkeys_hycom_qedt','flkeys_hycom_qenetqf',lagoff);




  [firstyr,ig,ig] = datevec(stn.flkeys_hycom_dt.date(1));
  dys = floor(now) - datenum(firstyr,1,1) + 1;
%%%% DEBUG
firstyr = 2008;
dys = 365;
%%%% DEBUG
  absflds = { ...
      'ndbc_ncep_30a_wind_stress_30_day_maximum', ...
      'flkeys_hycom_speed_7_day_maximum', ...
            };
  accflds = { ...
      'flkeys_hycom_netqf', ...
      'flkeys_hycom_dt_netqf', ...
      'netqf', ...
            };
  fh = plot_fluxes(stn,firstyr,1,dys,[],absflds,accflds);
  appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
  if ( lagoff ~= 0 )
    appendtitlename(sprintf(' lag:%d', lagoff));
  end;

%   plot_fluxes(stn,2004,109,6.5*365);
%   plot_fluxes(stn,2006,45,4.7*365);
%   figure(fh);

return;
