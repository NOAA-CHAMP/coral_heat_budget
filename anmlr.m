function stn = anmlr(stn)

error('This was a run-once testing script!');

    datapath = get_thesis_path('../data');

    if ( ~strcmpi(stn.station_name,'mlrf1') )
      error('This should only be run for station MLRF1!');
    end;

    warning('Off','LoadStationData:WrongName');
    stn2 = load_station_data('mlrf2');
    warning('On','LoadStationData:WrongName');
    flds = fieldnames(stn2);
    flds = flds(strmatch('bic_',flds));
    for ix = 1:length(flds)
      stn.(flds{ix}) = stn2.(flds{ix});
    end;
    stn2 = []; clear stn2;

    warning('Using *all* of BIC_SURF_PAR!');
    bic_valid_idx = 1:length(stn.bic_surf_par.date);

    stn = station_par_to_insol(stn,'bic_surf_par','bic_surf_dsrf',...
                               'bic_surf_usrf','bic_surf_srf',bic_valid_idx);

    [stn,Q0_factor] = ...
        station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ncep_relhumid',...
                          'ndbc_barom','ndbc_sea_t','bic_surf_srf','ncep_lrf',...
                          'bic_ncep');

    rhfld = 'ndbc_relhumid';
    [stn,Q0_factor] = ...
        station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                          'ndbc_barom','ndbc_sea_t','bic_surf_srf','ncep_lrf',...
                          'bic_ncep_rh');
    [stn,Q0_factor] = ...
        station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                          'ndbc_barom','ndbc_sea_t','bic_surf_srf','ncep_lrf',...
                          'bic_ncep_26','bic_surf_dsrf','ncep_dlrf');
%     [stn,Q0_factor] = ...
%         station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
%                           'ndbc_barom','ndbc_sea_t','bic_surf_srf','ncep_lrf',...
%                           'bic_ncep_30','bic_surf_dsrf','ncep_dlrf','ncep_precip');
%     [stn,Q0_factor] = ...
%         station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
%                           'ndbc_barom','ndbc_sea_t','bic_surf_srf','ncep_lrf',...
%                           'bic_ncep_30a','bic_surf_dsrf','ncep_dlrf','ncep_precip',...
%                           'ndbc_wind1_dir','default','default','default','default','default',true);

    matfname = fullfile(datapath, 'mlrf1_dt.mat');
    disp(['(Re-)saving to ' matfname]);
    station = stn;
    save(matfname,'station');
    station = []; clear station;

return;
