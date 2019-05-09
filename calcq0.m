function stn = calcq0(stnm_or_stn)
%function stn = calcq0(stnm_or_stn)
%
% Calculate or extract surface turbulent (QLH, QSH) and radiative (QSW, QLW)
% heat fluxes on struct STN or station named STNM.
%
% CALLS: LOAD_STATION_DATA, LOAD_ALL_NDBC_DATA, GET_NCEP_STATION,
% STATION_DEWP_TO_RELHUMID, STATION_RELHUMID_TO_SPECHUMID, STATION_PAR_TO_INSOL,
% STATION_BULK_LONGWAVE, STATION_NDBC_HFBULK, STATION_HEAT_FLUX.
%
% Last Saved Time-stamp: <Fri 2011-01-21 11:32:42  lew.gramer>

  timenow,

  datapath = get_thesis_path('../data');

  % Temporarily turn off block-headed MATLAB MORE behavior
  set_more off;

  if ( isstruct(stnm_or_stn) )
    stn = stnm_or_stn;
    stnm = stn.station_name;
  elseif ( ischar(stnm_or_stn) )
    stnm = stnm_or_stn;
    try
      stn = load_station_data(stnm);
    catch
      warning('NO raw station data - will use NDBC only.');
      stn = [];
      stn.station_name = stnm;
    end;
    if ( strcmpi(stnm,'mlrf1') )
      warning('Off','LoadStationData:WrongName');
      stn2 = load_station_data('mlrf2');
      warning('On','LoadStationData:WrongName');
      flds = fieldnames(stn2);
      flds = flds(strmatch('bic_',flds));
      for ix = 1:length(flds)
        stn.(flds{ix}) = stn2.(flds{ix});
      end;
      stn2 = []; clear stn2;
    end;
  else
    error('STNM_OR_STN must either be station struct or station name string!');
  end;

  if ( ~isfield(stn,'ndbc_air_t') )
    stn = load_all_ndbc_data(stn);
  end;

  if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
  end;


  if ( ~isfield(stn,'ncep_lrf') )
    stn = get_ncep_station(stn, 'narr');
  end;

  if ( strcmpi(stnm,'mlrf1') || strcmpi(stnm,'sanf1') || strcmpi(stnm,'dryf1') || strcmpi(stnm,'plsf1') )
    disp('Getting RH from SMKF1');
    smkf1 = load_all_ndbc_data([], 'smkf1');
    stn.ndbc_dew_t = smkf1.ndbc_dew_t;
    smkf1 = []; clear smkf1;
  end;
  if ( strcmpi(stnm,'fwyf1') )
    disp('Getting RH from LKWF1');
    lkwf1 = load_all_ndbc_data([], 'lkwf1');
    stn.ndbc_dew_t = lkwf1.ndbc_dew_t;
    lkwf1 = []; clear lkwf1;
  end;


  [stn,Q0_factor] = ...
      station_heat_flux(stn,'ncep_wind_speed','ncep_air_t','ncep_relhumid',...
                        'ncep_barom','ncep_sea_t','ncep_srf','ncep_lrf',...
                        'ncep_bulk');


  stn.ncep_heat_flux_term.date = stn.ncep_net_heat_flux.date;
  stn.ncep_heat_flux_term.data = (stn.ncep_net_heat_flux.data ./ Q0_factor) .* (60*60);

  stn.ncep_srf_term.date = stn.ncep_srf.date;
  stn.ncep_srf_term.data = (stn.ncep_srf.data ./ Q0_factor) .* (60*60);
  stn.ncep_lrf_term.date = stn.ncep_lrf.date;
  stn.ncep_lrf_term.data = (stn.ncep_lrf.data ./ Q0_factor) .* (60*60);
  stn.ncep_latent_flux_term.date = stn.ncep_latent_heat_flux.date;
  stn.ncep_latent_flux_term.data = (stn.ncep_latent_heat_flux.data ./ Q0_factor) .* (60*60);
  stn.ncep_sensible_flux_term.date = stn.ncep_sensible_heat_flux.date;
  stn.ncep_sensible_flux_term.data = (stn.ncep_sensible_heat_flux.data ./ Q0_factor) .* (60*60);


  if ( isfield(stn,'ndbc_relhumid') && (length(find(isfinite(stn.ndbc_relhumid.data))) > 0) )
    rhfld = 'ndbc_relhumid';
    shfld = 'ndbc_spechumid';
  elseif ( isfield(stn,'ndbc_dew_t') && (length(find(isfinite(stn.ndbc_dew_t.data))) > 0) )
    rhfld = 'ndbc_relhumid';
    stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t',rhfld);
    shfld = 'ndbc_spechumid';
  else
    rhfld = 'ncep_relhumid';
    shfld = 'ncep_spechumid';
    warning('No in situ Relative Humidity! Using NCEP...');
  end;
  if ( ~isfield(stn,shfld) )
    stn = station_relhumid_to_spechumid(stn,'ndbc_air_t',rhfld,shfld);
  end;


  licor_valid_idx = [];
  if ( isfield(stn,'licor_surf_par') && ~isempty(stn.licor_surf_par.date) )
    licor_valid_dates_fname = fullfile(datapath,[stnm '_valid_licor_surf_par.dat']);
    if ( ~exist(licor_valid_dates_fname,'file') )
      warning('Unable to use LICOR_SURF_PAR: No file "%s" found!', licor_valid_dates_fname);
    else
      licor_valid_dates = [];
      lvds = load(licor_valid_dates_fname);
      for ix = 1:size(lvds,1)
        dts = datenum(lvds(ix,1),lvds(ix,2),lvds(ix,3)) ...
              : (1/24) ...
              : datenum(lvds(ix,4),lvds(ix,5),lvds(ix,6),23,59,59);
        licor_valid_dates(end+1:end+length(dts),1) = dts;
      end;
      [licor_valid_idx,ig] = intersect_dates(stn.licor_surf_par.date,licor_valid_dates);

      stn = station_par_to_insol(stn,'licor_surf_par','licor_surf_dsrf',...
                                 'licor_surf_usrf','licor_surf_srf',licor_valid_idx);
    end;
  end;

  bic_valid_idx = [];
  if ( isfield(stn,'bic_surf_par') && ~isempty(stn.bic_surf_par.date) )
    bic_valid_dates_fname = fullfile(datapath,[stnm '_valid_bic_surf_par.dat']);
    if ( ~exist(bic_valid_dates_fname,'file') )
      warning('Using *all* of BIC_SURF_PAR: No file "%s" found!', bic_valid_dates_fname);
      %%%% ??? HACK
      bic_valid_idx = 1:length(stn.bic_surf_par.date);
    else
      bic_valid_dates = [];
      lvds = load(bic_valid_dates_fname);
      for ix = 1:size(lvds,1)
        dts = datenum(lvds(ix,1),lvds(ix,2),lvds(ix,3)) ...
              : (1/24) ...
              : datenum(lvds(ix,4),lvds(ix,5),lvds(ix,6),23,59,59);
        bic_valid_dates(end+1:end+length(dts),1) = dts;
      end;
      [bic_valid_idx,ig] = intersect_dates(stn.bic_surf_par.date,bic_valid_dates);
    end;

    stn = station_par_to_insol(stn,'bic_surf_par','bic_surf_dsrf',...
                               'bic_surf_usrf','bic_surf_srf',bic_valid_idx);
  end;

  if ( isempty(licor_valid_idx) && isempty(bic_valid_idx) )
    insfld = 'ncep_dsrf';
    uswfld = 'ncep_usrf';
    swfld = 'ncep_srf';
    warning('No in situ PAR! Using NCEP insolation...');
  elseif ( isempty(bic_valid_idx) )
    insfld = 'licor_surf_dsrf';
    uswfld = 'licor_surf_usrf';
    swfld = 'licor_surf_srf';
  else
    insfld = 'bic_surf_dsrf';
    uswfld = 'bic_surf_usrf';
    swfld = 'bic_surf_srf';
  end;


  stn = station_bulk_longwave(stn,'ndbc_air_t',shfld, ...
                              'ndbc_barom','ncep_cloud_cover','ndbc_sea_t', ...
                              'ncep_cloud_cover','ndbc_bulk_rh_dlrf', ...
                              'ndbc_bulk_rh_ulrf','ndbc_bulk_rh_lrf');

  stn = station_ndbc_hfbulk(stn,'ncep_srf','ndbc_bulk_rh_lrf');



  if ( ~strcmp(shfld,'ncep_spechumid') || ~strcmp(insfld,'ncep_dsrf') )

    if ( ~strcmp(shfld,'ncep_spechumid') )
      [stn,Q0_factor] = ...
          station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                            'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
                            'ndbc_bulk_rh');
      [stn,Q0_factor] = ...
          station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                            'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lrf',...
                            'ndbc_ncep_rh');

      [stn,Q0_factor] = ...
          station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                            'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
                            'ndbc_bulk_26','ncep_dsrf','ndbc_bulk_rh_dlrf');
      [stn,Q0_factor] = ...
          station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                            'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lrf',...
                            'ndbc_ncep_26','ncep_dsrf','ncep_dlrf');

      [stn,Q0_factor] = ...
          station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                            'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
                            'ndbc_bulk_30','ncep_dsrf','ndbc_bulk_rh_dlrf','ncep_precip');
      [stn,Q0_factor] = ...
          station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                            'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lrf',...
                            'ndbc_ncep_30','ncep_dsrf','ncep_dlrf','ncep_precip');

      [stn,Q0_factor] = ...
          station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                            'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
                            'ndbc_bulk_30a','ncep_dsrf','ndbc_bulk_rh_dlrf','ncep_precip',...
                            'ndbc_wind1_dir','default','default','default','default','default',true);
      [stn,Q0_factor] = ...
          station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                            'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lrf',...
                            'ndbc_ncep_30a','ncep_dsrf','ncep_dlrf','ncep_precip',...
                            'ndbc_wind1_dir','default','default','default','default','default',true);
    end;


    stn = station_bulk_longwave(stn,'ndbc_air_t',shfld,...
                                'ndbc_barom',insfld,'ndbc_sea_t', ...
                                'ndbc_cloud_cover','ndbc_bulk_dlrf',...
                                'ndbc_bulk_ulrf','ndbc_bulk_lrf');

    [stn,Q0_factor] = ...
        station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                          'ndbc_barom','ndbc_sea_t',swfld,'ndbc_bulk_lrf',...
                          'ndbc_bulk');


%     [stn,Q0_factor] = ...
%         station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ncep_relhumid',...
%                           'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
%                           'ncep_bulk_rh');

%     [stn,Q0_factor] = ...
%         station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ncep_relhumid',...
%                           'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
%                           'ncep_bulk_30','ncep_dsrf','ndbc_bulk_rh_dlrf','ncep_precip');
  end;

  % Last but not least, try it with satellite insolation data!
  if ( ~isfield(stn,'sat_par') )
    stn = get_satpar_insol(stn);
  end;
  stn = station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                          'ndbc_barom','ndbc_sea_t','sat_net_insol',...
                          'sat_net_longwave','sat');


  s = annocs(stnm);
  flds = fieldnames(s);
  for ix = 1:length(flds)
    fld = flds{ix};
    stn.(fld) = s.(fld);
  end;
  s = []; clear s;


  % Reset MORE state to previous value
  set_more;

  timenow,

return;
