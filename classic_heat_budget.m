function stn = classic_heat_budget(stn_or_stnm,RAPFX,KMPFX)
%ERROR('Replaced by STATION_HEAT_BUDGET');
error('Replaced by STATION_HEAT_BUDGET');
%function stn = classic_heat_budget(stn_or_stnm,RAPFX,KMPFX)
%
% Calculate all components of a reef heat budget EXCEPT horizontal convection
% and substrate (benthic) flux.
%
% Last Saved Time-stamp: <Tue 2011-04-19 10:40:03  Lew.Gramer>

  set_more off
  tic,

  datapath = get_thesis_path('../data');
  figspath = get_thesis_path('../figs');

  %%%
  %% All station struct fieldnames used to produce heat budget 

  % Tide data (or model)
  TIDEPFX = 'tmd_tide';
  hfld = [TIDEPFX '_i_depth'];
  tufld = [TIDEPFX '_u'];
  tvfld = [TIDEPFX '_u'];

  % Meteorology (in situ)
  afld = 'ndbc_air_t';
  sfld = 'ndbc_sea_t';
  %%%% DEBUG???
  % pfld = 'ndbc_barom';
  pfld = 'erai_barom';
  dfld = 'ndbc_dew_t';
  rhfld = 'ndbc_relhumid';
  qafld = 'ndbc_spechumid';
  qsfld = 'ndbc_sea_spechumid';
  Wfld = 'ndbc_wind1_speed';
  Dfld = 'ndbc_wind1_dir';
  Ufld = 'ndbc_wind1_u';
  Vfld = 'ndbc_wind1_v';


  if ( ~exist('RAPFX','var') || isempty(RAPFX) )
    RAPFX = 'erai';
    % RAPFX = 'ncep';
  end;

  % Meteorology (gridded/reanalysis)
  rfld = [RAPFX '_precip'];
  cfld = [RAPFX '_cloud_cover'];
  pblzfld = [RAPFX '_pblz'];
  alt_dfld = [RAPFX '_dew_t'];
  alt_rhfld = [RAPFX '_relhumid'];
  alt_qafld = [RAPFX '_spechumid'];
  alt_qsfld = [RAPFX '_sea_spechumid'];

  % DEBUG: ??? For now, ALWAYS use reanalysis humidities
  rhfld = alt_rhfld;
  qafld = alt_qafld;
  qsfld = alt_qsfld;


  % Radiative surface fluxes
  dsrfld = [RAPFX '_dsrf'];
  usrfld = [RAPFX '_usrf'];
  srfld = [RAPFX '_srf'];
  asrfld = ['absorbed_' RAPFX '_srf'];
  gamfld = ['absorbed_' RAPFX '_gamma'];

  dlrfld = [RAPFX '_ndbc_dlrf'];
  ulrfld = [RAPFX '_ndbc_ulrf'];
  lrfld = [RAPFX '_ndbc_lrf'];

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
  switch (KMPFX),
   case 'fkeys_hycom',	K_theta = 2.5;
   case 'gom_hycom',	K_theta = 20;
   otherwise,		error('Unknown km-scale model "%s"',KMPFX);
  end;
  kd2Tfld = [KMPFX '_diffused_heat'];


  % Total budget
  dTfld = [TURPFX '_' QEPFX '_dt'];


  if ( ischar(stn_or_stnm) )
    stnm = stn_or_stnm;
    stn.station_name = stnm;
  else
    stn = stn_or_stnm;
  end;
  if ( ~isfield(stn,'station_name') )
    error('Station struct STN had no .station_name field!');
  end;

  matfname = fullfile(datapath,[lower(stn.station_name) '_' dTfld '.mat']);

  if ( exist(matfname,'file') )

    disp(['Loading from ' matfname]);
    load(matfname,'station');
    flds = fieldnames(station);
    for fldix=1:length(flds)
      fld = flds{fldix};
      stn.(fld) = station.(fld);
    end;
    station = []; clear station;

  else


    %%%
    %% Station coordinates and bathymetry

    if ( ~isfield(stn,'lon') )
      [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
    end;


    %%%
    %% Meteorology (in situ and gridded/reanalysis)

    if ( ~isfield(stn,afld) )
      stn = load_all_ndbc_data(stn);
    end;
    if ( ~isfield(stn,hfld) )
      stn = station_tmd_tide(stn);
    end;
    if ( ~isfield(stn,Ufld) )
      [Wix,Dix] = intersect_dates(stn.(Wfld).date,stn.(Dfld).date);
      stn.(Ufld).date = stn.(Wfld).date(Wix);
      stn.(Vfld).date = stn.(Wfld).date(Wix);
      [stn.(Ufld).data,stn.(Vfld).data] = spddir_to_uv(stn.(Wfld).data(Wix),stn.(Dfld).data(Dix));
    end;
    if ( ~isfield(stn,dsrfld) || ~isfield(stn,cfld) )
      switch (RAPFX),
       case 'erai',		stn = get_erai_station(stn);
       case 'ncep',		stn = get_ncep_station(stn,'narr');
        %%% Not Yet Implemented
        % case 'era40',		stn = get_era40_station(stn);
        % case 'cfsr',		stn = get_ncep_station(stn,'cfsr');
       otherwise,			error('Unknown gridded/reanalysis dataset "%s"',RAPFX);
      end;
    end;
    if ( ~isfield(stn,dfld) || ~isfield(stn.(dfld),'date') || numel(stn.(dfld).date)<(numel(stn.(afld).date)/3) )
      if ( isfield(stn,rhfld) )
        disp([rhfld '->' dfld]);
        stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
      elseif ( isfield(stn,qafld) )
        disp([qafld '->' dfld]);
        stn = station_spechumid_to_relhumid(stn,afld,qafld,rhfld);
        stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
      elseif ( isfield(stn,alt_rhfld) )
        dfld = alt_dfld;
        rhfld = alt_rhfld;
        qafld = alt_qafld;
        qsfld = alt_qsfld;
        % Some reanalysis have dew_t, some have relhumid!
        if ( ~isfield(stn,dfld) )
          disp([rhfld '->' dfld]);
          stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
        end;
      else
        error('Found no humidity/dew-point data (%s,%s,%s,%s)!',dfld,alt_dfld,rhfld,alt_rhfld);
      end;
    end;
    if ( ~isfield(stn,rhfld) )
      stn = station_dewp_to_relhumid(stn,afld,dfld,rhfld);
    end;
    if ( ~isfield(stn,qafld) )
      stn = station_relhumid_to_spechumid(stn,afld,rhfld,qafld);
    end;
    if ( ~isfield(stn,qsfld) )
      stn.(qsfld).date = stn.(sfld).date;
      % Stommel's seawater salinity adjustment to saturated specific humidity
      stn.(qsfld).data = 0.98 .* relhumid_to_spechumid(stn.(sfld).data,100);
    end;


    %%%
    %% Waves and Stokes Drift

    if ( ~isfield(stn,whfld) )
      % If waves not from reanalysis, user must want model (WaveWatch III) or wind estimate
      switch (WAVEPFX),
       case 'ww3',		stn = get_ww3_station(stn);
       case 'ndbc',		stn = station_wind_to_wave(stn,Wfld,Dfld,wpfld,whfld,wdfld);
       otherwise,			error('Unknown wave source "%s"',WAVEPFX);
      end;
    end;
    if ( ~isfield(stn,ssufld) )
      %station_stokes_drift(stn,ssfld,sdfld,sufld,svfld,wsfld,wdfld,hsfld,tpfld,tdfld)
      stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wfld,Dfld,whfld,wpfld,wdfld);
    end;  


    %%%
    %% Kilometer-scale Ocean Data

    if ( ~isfield(stn,ufld) )
      switch (KMPFX),
       case 'fkeys_hycom',	stn = get_fkeys_hycom(stn);
       case 'gom_hycom',	stn = get_gom_hycom(stn);
       otherwise,		error('Unknown km-scale model "%s"',KMPFX);
      end;
      more off
    end;
    if ( ~isfield(stn,Tfld) )
      stn = calc_field_terms(stn,Tfld);
    end;


    %%%
    %% Quasi-Eulerian (km + surface) Currents

    if ( ~isfield(stn,qeufld) )
      stn = calc_quasi_eulerian(stn,STOKESPFX,KMPFX,QEPFX);
    end;


    %%%
    %% Radiative Fluxes

    if ( ~isfield(stn,asrfld) )
      % If absorbed short-wave not already present, user wants abs. calculation
      stn = station_absorbed_insolation(stn,asrfld,srfld,hfld,[],[],gamfld,qbfld);
    end;
    if ( ~isfield(stn,lrfld) )
      % If reanalysis long-wave flux was not specified, user wants a bulk estimate
      %station_bulk_longwave(stn,afld,qfld,pfld,dsrfld,sfld,cfld,dlrf,ulrf,lrf)
      %stn = station_bulk_longwave(stn,afld,qafld,pfld,dsrfld,sfld,cfld,dlrfld,ulrfld,lrfld);
      stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);
    end;


    %%%
    %% Turbulent and Net Fluxes

    if ( ~isfield(stn,pblzfld) )
      % Default: 600m
      pblzfld = 600;
    end;

    if ( ~isfield(stn,q0fld) )
      %station_heat_flux(stn,wfld,afld,qfld,pfld,tfld,sfld,lfld,PFX,dsfld,dlfld,prfld,wdfld,oufld,ovfld,wpfld,whfld,pblzfld,doWarm,doPlot)
      stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                              pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
                              Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,true);
      %station_heat_flux_term(stn,nffld,htfld,tfld,sfld,dfld)
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],hfld);
    end;


    disp(['Saving to ' matfname]);
    station = stn;
    save(matfname,'station');
    station = []; clear station;


    %%%
    %% Eulerian (km + Stokes) Heat Advection, Km-scale Heat Diffusion
    if ( ~isfield(stn,dTfld) )
      stn = station_calc_udotdelt(stn,qeufld,qevfld,Tfld,kmtfld,...
                                  ['raw_' udTfld],udTfld,...
                                  qtfld,qtAdvfld);
      stn = station_calc_kdel2t(stn,K_theta,Tfld,...
                                ['raw_' kd2Tfld],kd2Tfld,...
                                qtAdvfld,dTfld);
      stn = station_heat_flux_term_inverse(stn,[dTfld '_heat_flux'],...
                                           dTfld,sfld,[],hfld);
    end;

    disp(['RE-saving to ' matfname]);
    station = stn;
    save(matfname,'station');
    station = []; clear station;


  end; %if ( exist(matfname,'file') ) else

  toc,
  set_more;

  stn = compare_heat_budgets(stn,sfld,RAPFX,KMPFX);
  % ylim([0 400]);
  figfname= fullfile(figspath,[lower(stn.station_name) '_' dTfld '.tiff']);
  print('-dtiff',figfname);

return;
