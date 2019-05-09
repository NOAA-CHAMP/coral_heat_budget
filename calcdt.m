function station = calcdt(stnm)
%function station = calcdt(stnm)
%
% Calculate dT/dt as the sum of net surface heat flux (Q0) and advective heat
% flux terms, per eqn. (1) of Gramer, MS Thesis Proposal, 2009.
%
% Last Saved Time-stamp: <Wed 2011-07-06 06:25:44  lew.gramer>

  datapath = get_thesis_path('../data');

  matfname = fullfile(datapath, sprintf('%s_dt.mat',stnm));

  if ( exist(matfname,'file') )
    disp(['Loading from ' matfname]);
    load(matfname,'station');

  else

    % Fields to KEEP from old "_thesis.mat" file
    flds = { ...
        'station_name','lon','lat',...
        ...
        'licor_surf_par','bic_surf_par','bic_shallow_par',...
        'tide','wave_height','WaveHgt','wave_i_depth','IDepth','wave_tide_height','TidalHgt',...
        ...
        'ndbc_wind1_dir','ndbc_wind1_speed','ndbc_wind1_u','ndbc_wind1_v',...
        'ndbc_air_t','ndbc_barom','ndbc_sea_t','ndbc_tide',...
        'ndbc_sigwavehgt','ndbc_avgwaveper','ndbc_avgwavedir',...
        'ndbc_relhumid','ndbc_spechumid',...
        ...
        'ncep_wind_dir','ncep_wind_speed','ncep_air_t','ncep_barom','ncep_sea_t',...
        'ncep_cloud_cover','ncep_relhumid','ncep_spechumid',...
        'ncep_par','ncep_parW',...
        'ncep_dsrf','ncep_usrf','ncep_srf',...
        'ncep_dlrf','ncep_ulrf','ncep_lrf',...
        'ncep_latent_heat_flux','ncep_sensible_heat_flux',...
        'ncep_net_heat_flux','ncep_heat_flux_term',...
        ...
        'sat_net_heat_flux','sat_heat_flux_term','sat_wind_stress',...
        ...
        'nocs_net_heat_flux','nocs_heat_flux_term',...
        ...
        'ndbc_hfbulk_net_heat_flux','ndbc_hfbulk_heat_flux_term',...
        'ndbc_bulk_net_heat_flux','ndbc_bulk_heat_flux_term','ndbc_bulk_wind_stress',...
        'ndbc_bulk_rh_net_heat_flux','ndbc_bulk_rh_heat_flux_term','ndbc_bulk_rh_wind_stress',...
        'ndbc_bulk_26_net_heat_flux','ndbc_bulk_26_heat_flux_term','ndbc_bulk_26_wind_stress',...
        'ndbc_bulk_30_net_heat_flux','ndbc_bulk_30_heat_flux_term','ndbc_bulk_30_wind_stress',...
        'ndbc_bulk_30a_net_heat_flux','ndbc_bulk_30a_heat_flux_term','ndbc_bulk_30a_wind_stress',...
        ...
        'ndbc_ncep_rh_net_heat_flux','ndbc_ncep_rh_heat_flux_term','ndbc_ncep_rh_wind_stress',...
        'ndbc_ncep_26_net_heat_flux','ndbc_ncep_26_heat_flux_term','ndbc_ncep_26_wind_stress',...
        'ndbc_ncep_30_net_heat_flux','ndbc_ncep_30_heat_flux_term','ndbc_ncep_30_wind_stress',...
        'ndbc_ncep_30a_net_heat_flux','ndbc_ncep_30a_heat_flux_term','ndbc_ncep_30a_wind_stress',...
        ...
        'ndbc_hfbulk_latent_heat_flux','ndbc_hfbulk_latent_flux_term',...
        'ndbc_hfbulk_sensible_heat_flux','ndbc_hfbulk_sensible_flux_term',...
        ...
        'ndbc_ncep_30a_latent_heat_flux','ndbc_ncep_30a_latent_flux_term',...
        'ndbc_ncep_30a_sensible_heat_flux','ndbc_ncep_30a_sensible_flux_term',...
        'ndbc_ncep_30a_shortwave_flux_term','ndbc_ncep_30a_longwave_flux_term',...
        ...
        'bic_ncep_net_heat_flux','bic_ncep_heat_flux_term','bic_ncep_wind_stress',...
        'bic_ncep_rh_net_heat_flux','bic_ncep_rh_heat_flux_term','bic_ncep_rh_wind_stress',...
        'bic_ncep_26_net_heat_flux','bic_ncep_26_heat_flux_term','bic_ncep_26_wind_stress',...
        'bic_ncep_30_net_heat_flux','bic_ncep_30_heat_flux_term','bic_ncep_30_wind_stress',...
        'bic_ncep_30a_net_heat_flux','bic_ncep_30a_heat_flux_term','bic_ncep_30a_wind_stress',...
           };
%         'ncep_bulk_net_heat_flux','ncep_bulk_heat_flux_term',...
%         'ncep_bulk_net_heat_flux','ncep_bulk_heat_flux_term','ncep_bulk_wind_stress',...
%         ...

    oldmatfname = fullfile(datapath, sprintf('%s_thesis.mat',stnm));
    if ( exist(oldmatfname,'file') )
      disp(['Reloading from ' oldmatfname]);
      olds = load(oldmatfname);
    else
      warning('RECALCULATING Q0 FROM SCRATCH: NO "%s" FILE!',oldmatfname);
      olds.(stnm) = calcq0(stnm);
      save(oldmatfname,'-struct','olds',stnm);
    end;

    for fldix = 1:length(flds)
      fld = flds{fldix};
      if ( isfield(olds.(stnm),fld) )
        station.(fld) = olds.(stnm).(fld);
      end;
    end;
    olds = [];
    clear olds;

    station = get_ww3_station(station);

    disp(['Saving to ' matfname]);
    save(matfname,'station');

  end;


  % Calculate heat advection U * grad(T) :

  % Extract and estimate point ocean currents, and 1km satellite SST fields


  % HACK workaround for now-fixed bug in UV_TO_DIR
  station.ww3_peakwavedir.data(isnan(station.ww3_peakwavedir.data)) = 0;


  station = verify_variable(station,'ndbc_wind1_u');
  station = verify_variable(station,'ndbc_wind1_v');

  % PFX = 'ndbc_ncep_30a_';
  PFX = '';

  % Fields to REMOVE here so we can recalculate them below
  newflds = { ...
      [PFX 'stokes_drift_speed'],[PFX 'stokes_drift_dir'],...
      [PFX 'stokes_drift_u'],[PFX 'stokes_drift_v'],...
      'ndbc_wind1_wvper','ndbc_wind1_wvhgt','ndbc_wind1_wvdir',...
      [PFX 'wind_stokes_drift_speed'],[PFX 'wind_stokes_drift_dir'],...
      [PFX 'wind_stokes_drift_u'],[PFX 'wind_stokes_drift_v'],...
      'global_hycom_u','global_hycom_v',...
      [PFX 'stokes_drift_u_1_day_lowpass'],...
      [PFX 'stokes_drift_v_1_day_lowpass'],...
      [PFX 'quasi_eulerian_u'],[PFX 'quasi_eulerian_v'],...
      [PFX 'advected_heat'],...
      'ndbc_hfbulk_heat_flux_term',...
      'ndbc_bulk_rh_dt','ndbc_bulk_26_dt','ndbc_bulk_30_dt','ndbc_bulk_30a_dt',...
      'ndbc_ncep_rh_dt','ndbc_ncep_26_dt','ndbc_ncep_30_dt','ndbc_ncep_30a_dt',...
      'bic_ncep_dt','bic_ncep_rh_dt','bic_ncep_26_dt','bic_ncep_30_dt','bic_ncep_30a_dt',...
            };
  for fldix = 1:length(newflds)
    fld = newflds{fldix};
    if ( isfield(station,fld) )
      station = rmfield(station,fld);
    end;
  end;


  station = ...
      station_stokes_drift(station,...
                           [PFX 'stokes_drift_speed'],[PFX 'stokes_drift_dir'],...
                           [PFX 'stokes_drift_u'],[PFX 'stokes_drift_v'],...
                           'ndbc_wind1_speed','ndbc_wind1_dir',...
                           'ww3_sigwavehgt','ww3_peakwaveper','ww3_peakwavedir');

  station = ...
      station_wind_to_wave(station,...
                           'ndbc_wind1_speed','ndbc_wind1_dir',...
                           'ndbc_wind1_wvper','ndbc_wind1_wvhgt','ndbc_wind1_wvdir');

  station = ...
      station_stokes_drift(station,...
                           [PFX 'wind_stokes_drift_speed'],[PFX 'wind_stokes_drift_dir'],...
                           [PFX 'wind_stokes_drift_u'],[PFX 'wind_stokes_drift_v'],...
                           'ndbc_wind1_speed','ndbc_wind1_dir',...
                           'ndbc_wind1_wvhgt','ndbc_wind1_wvper','ndbc_wind1_wvdir');

  rawstokesu = [PFX 'stokes_drift_u'];
  rawstokesv = [PFX 'stokes_drift_v'];
  % rawstokesu = [PFX 'wind_stokes_drift_u'];
  % rawstokesv = [PFX 'wind_stokes_drift_v'];

  % stokesu = [rawstokesu '_1_day_lowpass'];
  % stokesv = [rawstokesv '_1_day_lowpass'];
  % station = verify_variable(station, stokesu);
  % station = verify_variable(station, stokesv);
  stokesu = [rawstokesu];
  stokesv = [rawstokesv];


  hycomfname = fullfile(datapath,[stnm '_hycom.mat']);
  if ( exist(hycomfname,'file') )
    disp(['Loading from presaved ' hycomfname]);
    load(hycomfname);

    dts = stn.global_hycom_u.date(1):(1/24):stn.global_hycom_u.date(end);
    station.global_hycom_u.date = dts';
    station.global_hycom_u.data = interp1(stn.global_hycom_u.date,...
                                          stn.global_hycom_u.data,...
                                          dts,'pchip',nan)';
    station.global_hycom_v.date = dts';
    station.global_hycom_v.data = interp1(stn.global_hycom_v.date,...
                                          stn.global_hycom_v.data,...
                                          dts,'pchip',nan)';
    stn = []; clear stn;
  else
    warning('Substituting all-zero surface currents for missing "%s"!', hycomfname);
    station.global_hycom_u.date = station.(rawstokesu).date;
    station.global_hycom_u.data = repmat(0,size(station.global_hycom_u.date));
    station.global_hycom_v = station.global_hycom_u;
  end;


  [hix,wix] = intersect_dates(station.global_hycom_u.date,...
                              station.(stokesu).date);

  euleru = [PFX 'quasi_eulerian_u'];
  eulerv = [PFX 'quasi_eulerian_v'];
  eulers = [PFX 'quasi_eulerian_speed'];
  eulerd = [PFX 'quasi_eulerian_dir'];
  station.(euleru).date = station.global_hycom_u.date(hix);
  station.(euleru).data = station.global_hycom_u.data(hix) + ...
      station.(stokesu).data(wix);
  station.(eulerv).date = station.global_hycom_v.date(hix);
  station.(eulerv).data = station.global_hycom_v.data(hix) + ...
      station.(stokesv).data(wix);
  station.(eulers).date = station.(euleru).date;
  station.(eulers).data = uv_to_spd(station.(euleru).data,station.(eulerv).data);
  station.(eulerd).date = station.(euleru).date;
  station.(eulerd).data = uv_to_dir_curr(station.(euleru).data,station.(eulerv).data);

  station = get_avhrr_weekly_field(station);

  UdotdelT = [PFX 'advected_heat'];
  station = station_advect_field(station,UdotdelT,...
                                 euleru,eulerv,...
                                 'avhrr_weekly_sst');


  % Extract or generate values for varying depth 'h'
  station = station_tides_to_idepths(station);

  idepthfld = 'ndbc_tide_i_depth';
  if ( ~isfield(station,idepthfld) )
    idepthfld = 'tide_i_depth';
    if ( ~isfield(station,idepthfld) )
      idepthfld = 'wave_i_depth';
      if ( ~isfield(station,idepthfld) )
        if ( isfield(station,'IDepth') )
          station.(idepthfld) = station.IDepth;
          station = rmfield(station,'IDepth');
        else
          station = station_tmd_tide(station,station.ndbc_sea_t.date);
          idepthfld = 'tmd_tide_i_depth';
          if ( ~isfield(station,idepthfld) )
            idepthfld = 'tpxo_tide_i_depth';
            if ( ~isfield(station,idepthfld) )
              idepthfld = [];
            end;
          end;
        end;
      end;
    end;
  end;



  for q0basename = { ...
      'ndbc_hfbulk',...
      'ndbc_bulk_rh','ndbc_bulk_26','ndbc_bulk_30','ndbc_bulk_30a',...
      'ndbc_ncep_26','ndbc_ncep_30','ndbc_ncep_30a',...
      'bic_ncep','bic_ncep_rh','bic_ncep_26','bic_ncep_30','bic_ncep_30a',...
                   }

    nffld = [q0basename{:} '_net_heat_flux'];
    htfld = [q0basename{:} '_heat_flux_term'];
    dtfld = [q0basename{:} '_dt'];

    if ( isfield(station,nffld) )
      if (isfield(station,htfld)); station = rmfield(station,htfld);      end;
      station = station_heat_flux_term(station,nffld,htfld,'ndbc_sea_t',[],idepthfld);
      for termfld = {'latent','sensible'}
        lffld = [q0basename{:} '_' termfld{:} '_heat_flux'];
        ltfld = [q0basename{:} '_' termfld{:} '_flux_term'];
        if ( isfield(station,lffld) )
          if (isfield(station,ltfld)); station = rmfield(station,ltfld);    end;
          station = station_heat_flux_term(station,lffld,ltfld,'ndbc_sea_t',[],idepthfld);
        end;
      end;
    end;

    if ( isfield(station,htfld) )
      %DEBUG:      disp(htfld);
      [ix1,ix2] = intersect_dates(station.(htfld).date, ...
                                  station.(UdotdelT).date);
      station.(dtfld).date = station.(htfld).date(ix1);
      station.(dtfld).data = station.(UdotdelT).data(ix2) + ...
          station.(htfld).data(ix1);
    else
      %DEBUG:
      disp(['No field ' htfld]);
    end;
  end;

  %DEBUG:
  if (exist('idepthfld') && ~isempty(idepthfld)); disp(['Used depth time series ' idepthfld]); end;

  % station = verify_variable(station,'ndbc_sea_t_1_day_decimate');
  % station = verify_variable(station,'ndbc_ncep_26_heat_flux_term_1_day_decimate');
  % station = verify_variable(station,'ndbc_ncep_26_dt_1_day_decimate');
  % station = verify_variable(station,'ndbc_ncep_30_heat_flux_term_1_day_decimate');
  % station = verify_variable(station,'ndbc_ncep_30_dt_1_day_decimate');
  % station = verify_variable(station,'ndbc_ncep_30a_heat_flux_term_1_day_decimate');
  % station = verify_variable(station,'ndbc_ncep_30a_dt_1_day_decimate');

  % station = verify_variable(station,'ndbc_sea_t_1_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_26_heat_flux_term_1_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_26_dt_1_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_30_heat_flux_term_1_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_30_dt_1_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_30a_heat_flux_term_1_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_30a_dt_1_day_lowpass');

  % station = verify_variable(station,'ndbc_sea_t_7_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_26_heat_flux_term_7_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_26_dt_7_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_30_heat_flux_term_7_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_30_dt_7_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_30a_heat_flux_term_7_day_lowpass');
  % station = verify_variable(station,'ndbc_ncep_30a_dt_7_day_lowpass');

% %%%% ??? DEBUG
%   disp(['NOT (Re-)saving to ' matfname]);
  disp(['(Re-)saving to ' matfname]);
  save(matfname,'station');

return;
