function stn = foo(stn_or_stnm)

  %%%
  %% Call SCRIPT to set:
  %% Variable-name prefixes ("PFX") for various input and output datasets; AND,
  %% All station struct fieldnames used to produce heat budget 
  station_heat_budget_field_names;

  stn = get_station_from_station_name(stn_or_stnm);

    if ( ~isfield(stn,afld) )
      warning('off','Ecoforecasts:mergedNonTS');
      switch (ISPFX),
       case 'ndbc',		stn = load_all_ndbc_data(stn);
       otherwise,		error('Unknown in situ dataset "%s"',ISPFX);
      end;
      warning('on','Ecoforecasts:mergedNonTS');
    end;
    if ( ~isfield(stn,Ufld) )
      [Wix,Dix] = intersect_dates(stn.(Wfld).date,stn.(Dfld).date);
      stn.(Ufld).date = stn.(Wfld).date(Wix);
      stn.(Vfld).date = stn.(Wfld).date(Wix);
      [stn.(Ufld).data,stn.(Vfld).data] = spddir_to_uv(stn.(Wfld).data(Wix),stn.(Dfld).data(Dix));
    end;
    if ( ~isfield(stn,whfld) )
      % If waves not from reanalysis, user must want model (WaveWatch III) or wind estimate
      switch (WAVEPFX),
       case 'ww3',		stn = get_ww3_station(stn);
       case 'ndbc',		stn = station_wind_to_wave(stn,Wfld,Dfld,wpfld,whfld,wdfld);
       case 'erai',		stn = get_erai_station(stn); %IF NOT ALSO OUR REANALYSIS DATASET
       otherwise,		error('Unknown wave source "%s"',WAVEPFX);
      end;
    end;

      %%%% ??? DEBUG: Low-pass filter winds for quasi-Eulerian currents
      Ulpfld = [Ufld '_30_hour_lowpass'];
      Vlpfld = [Vfld '_30_hour_lowpass'];
      Wlpfld = [Wfld '_30_hour_lowpass'];
      Dlpfld = [Dfld '_30_hour_lowpass'];
      stn = verify_variable(stn,Ulpfld);
      stn = verify_variable(stn,Vlpfld);
      stn.(Wlpfld).date = stn.(Ulpfld).date;
      stn.(Wlpfld).data = uv_to_spd(stn.(Ulpfld).data,stn.(Vlpfld).data);
      stn.(Dlpfld).date = stn.(Ulpfld).date;
      stn.(Dlpfld).data = uv_to_dir(stn.(Ulpfld).data,stn.(Vlpfld).data);
      stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wlpfld,Dlpfld,whfld,wpfld,wdfld);

return;


  [six,aix,qaix,qsix,Wix,hix,qlhix,qshix,swix,lwix] = ...
      intersect_all_dates([], stn.ndbc_sea_t.date,stn.ndbc_air_t.date,...
                          stn.ndbc_spechumid.date,stn.ndbc_sea_spechumid.date,...
                          stn.ndbc_wind1_speed.date,...
                          stn.tmd_tide_i_depth.date,...
                          stn.ndbc_ncep_30a_latent_heat_flux.date,...
                          stn.ndbc_ncep_30a_sensible_heat_flux.date,...
                          stn.ncep_dsrf.date,stn.ncep_dlrf.date);

  dts = stn.ndbc_sea_t.date(six);

return;

1;

dlon = 0.02;
dlat = 0.02;

lon = -80:dlon:-79;
lat = 24:dlat:25;

lonp=[-79.20 -79.30];
latp=[+24.13 +24.17];
% lonp=[-79.2 -79.2];
% latp=[24.5 24.1];
% lonp=[-79.2 -79.3];
% latp=[24.5 24.5];

figure;
hold on;
xlim([min(lon) max(lon)]);
ylim([min(lat) max(lat)]);

plot(lonp,  latp, 'bo');


Y = lonp + (i .* latp);
XI = linspace(1, length(Y), 1000);
y = interp1(Y, XI, 'pchip');
lonl = real(y);
latl = imag(y);

plot(lonl,  latl, 'r.');


lons = round(lonl ./ (dlon/2)) .* (dlon/2);
lats = round(latl ./ (dlat/2)) .* (dlat/2);

uniqpts = find(diff(lons) | diff(lats));
lons = lons([1 uniqpts]);
lats = lats([1 uniqpts]);

plot(lons,  lats, 'gs');
