function stn = station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,forceRecalc)
%function stn = station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,forceRecalc)
%
% Calculate all components of a reef heat budget: air-sea flux, advection
% and diffusion, horizontal convection, and substrate (benthic) flux. The
% various optional *PFX string args specify datasets from which Reanalysis
% (RA), kilometer-scale advection/diffusion (KM), in situ data (IS), tide
% height/currents (TIDE), and surface waves (WAVE) are extracted. Optional
% cell array arg SUBSTITUTE_FIELD_NAMES specifies alternate field names (see
% below). By default, save result of basic heat budget in a .MAT file, and
% reloads it on succeeding calls. If optional FORCERECALC==True, recalculate
% budget, i.e., don't reload results from any previously saved .MAT file.
%
% SAMPLE usages for the optional SUBSTITUTE_FIELD_NAMES cell array arg:
%  {'sfld','hourly_misst_sst'} evaluate heat budget based on MISST instead of in situ
%  {'comparisonsfld','ndbc_sea_t'} compare result with actual in situ sea temperature
%  {'reanalysis_shortwave',true} use upward shortwave from reanalysis instead of bulk
%
% Last Saved Time-stamp: <Tue 2012-07-31 16:40:18  lew.gramer>

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

  %%%% ??? DEBUG
  if ( ~strcmp(sfld,[ISPFX '_sea_t']) )
    disp(['Using sea temperature ' sfld]);
  end;
  if ( ~strcmp(afld,[ISPFX '_air_t']) )
    disp(['Using air temperature ' afld]);
  end;
  %%%% ??? DEBUG


  if ( ischar(stn_or_stnm) )
    stnm = stn_or_stnm;
    stn.station_name = stnm;
  else
    stn = stn_or_stnm;
  end;
  if ( ~isfield(stn,'station_name') )
    error('Station struct STN had no .station_name field!');
  end;

  if ( ~exist('forceRecalc','var') || isempty(forceRecalc) )
    forceRecalc = false;
  end;

  set_more off
  tic,

  datapath = get_thesis_path('../data');
  figspath = get_thesis_path('../figs');

  stn.commentstr = '';

  %%%
  %% Heat Budget Options structure
  if ( isfield(stn,'station_heat_budget_options') )
    opts = stn.station_heat_budget_options;
  else
    opts = station_heat_budget_options(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names);
  end;


  %%%
  %% Check for presaved results

  matfname = fullfile(datapath,[lower(stn.station_name) '_' TURPFX '_' QEPFX '_heat_budget.mat']);

  if ( ~forceRecalc && exist(matfname,'file') )

    disp(['Loading from ' matfname]);
    load(matfname,'station');
    flds = fieldnames(station);
    for fldix=1:length(flds)
      fld = flds{fldix};
      if ( ~isfield(stn,fld) )
        stn.(fld) = station.(fld);
      else
        disp(['Field already in STN: ' fld]);
      end;
    end;
    stn.commentstr = station.commentstr;
    station = []; clear station;


  else

    disp(['Calculating... Will save to ' matfname]);

    %%%
    %% Station coordinates and bathymetry

    if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
      [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
    end;
    if ( ~isfield(stn,bathyfld) )
      stn = get_ngdc_bathy_station(stn);
    end;
    if ( ~isfield(stn,slopefld) )
      stn = station_ngdc_offshore_slope(stn);
    end;
    if ( ~isfield(stn,bathorifld) )
      stn = station_optimal_isobath_orientation(stn);
    end;


    %%%
    %% Sea temperature (in situ NDBC, site-specific, remote sensing, or reanalysis)

    if ( ~isfield(stn,sfld) )
      if ( regexp(sfld,'misst_sst') )
        stn = get_misst_station(stn);
        stn.hourly_misst_sst = interp_ts(stn.misst_sst);
      elseif ( regexp(sfld,'avhrr_weekly_sst') )
        stn = get_avhrr_weekly_field(stn,true);
        stn.hourly_avhrr_weekly_sst = interp_ts(stn.avhrr_weekly_sst);
      elseif ( regexp(sfld,'erai_') )
        if ( ~strcmpi(ISPFX,'erai') && ~strcmpi(RAPFX,'erai') )
          stn = get_erai_station(stn);
        end;
      elseif ( regexp(sfld,'microcat') )
        stn = get_looe1_microcat(stn);
      elseif ( regexp(sfld,'adcp') )
        stn = get_looe1_adcp(stn);
      elseif ( regexp(sfld,'ndbc') )
        if ( ~strcmpi(ISPFX,'ndbc') )
          stn = load_all_ndbc_data(stn);
        end;
      elseif ( regexp(sfld,'^(sea|ct_|ctd_)') )
        if ( ~strcmpi(ISPFX,'icon') )
          stn = load_station_data(stn);
        end;
      else
        error(['Unknown sea temperature field ',sfld]);
      end;
    end;


    %%%
    %% Meteorology (in situ and gridded/reanalysis)

    if ( ~isfield(stn,afld) )
      warning('off','Ecoforecasts:mergedNonTS');
      switch (ISPFX),
       case 'ndbc',		stn = load_all_ndbc_data(stn);
       case 'erai',		stn = get_erai_station(stn);
       case 'icon',		stn = load_station_data(stn);
       otherwise,		error('Unknown in situ dataset "%s"',ISPFX);
      end;
      warning('on','Ecoforecasts:mergedNonTS');

      % E.g., if SFLD is 'ndbc_sea_t_12_hour_lowpass'
      stn = verify_variable(stn, afld);
      stn = verify_variable(stn, pfld);
      stn = verify_variable(stn, Wfld);
      stn = verify_variable(stn, sfld);
    end;

    if ( ~isfield(stn,hfld) || ~isfield(stn,tufld) || ~isfield(stn,tvfld) )
      switch (TIDEPFX),
       case 'tmd_tide',		stn = station_tmd_tide(stn);
       case 'tpxo_tide',	stn = station_tmd_tide(stn);
       case {'ndbc','icon'},	warning('In situ tidal currents not available!');
        stn.(tufld).date = stn.(hfld).date; stn.(tufld).data = repmat(0,size(stn.(tufld).date));
        stn.(tvfld).date = stn.(hfld).date; stn.(tvfld).data = repmat(0,size(stn.(tvfld).date));
        stn.(tspdfld).date = stn.(hfld).date; stn.(tspdfld).data = repmat(0,size(stn.(tspdfld).date));
        stn.(tdirfld).date = stn.(hfld).date; stn.(tdirfld).data = repmat(0,size(stn.(tdirfld).date));
       otherwise,		error('Unknown tide dataset "%s"',TIDEPFX);
      end;
    end;
    if ( ~isfield(stn,mhfld) )
      % % Calculate the mean tide depth experienced by a watermass moving over
      % % an M2 tidal ellipse centered on the coordinates of our station
      % stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld,tufld,tvfld);

      % Just assume (small) cross-shore tidal ellipse to estimate mean depth
      stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
    end;

    if ( ~isfield(stn,dsfld) )
      stn.(dsfld).date = stn.(sfld).date(1:end-1);
      stn.(dsfld).data = diff(stn.(sfld).data);
      stn = filter_gaps(stn,sfld,dsfld,(1.5/24));
    end;
    if ( ~isfield(stn,dsffld) )
      stn = station_heat_flux_term_inverse(stn,dsffld,dsfld,sfld,[],mhfld);
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
       otherwise,		error('Unavailable gridded/reanalysis dataset "%s"',RAPFX);
      end;
      stn = station_heat_flux_term(stn,raq0fld,raqtfld,sfld,[],nanmean(stn.(mhfld).data));
    end;

    % Some datasets have dew temp, some have relhumid, some have spechumid!
    if ( ~isfield(stn,dfld) && ~isfield(stn,rhfld) && ~isfield(stn,qafld) )
      error('Found neither dew-point nor (rel/spec) humidity data (%s,%s,%s)!',dfld,rhfld,qafld);
    elseif ( isfield(stn,dfld) )
      if ( ~isfield(stn,rhfld) )
        disp([dfld '->' rhfld]);
        stn = station_dewp_to_relhumid(stn,afld,dfld,rhfld);
      end;
      if ( ~isfield(stn,qafld) )
        disp([rhfld '->' qafld]);
        stn = station_relhumid_to_spechumid(stn,afld,rhfld,qafld);
      end;
    elseif ( isfield(stn,rhfld) )
      disp([rhfld '->' dfld]);
      stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
      if ( ~isfield(stn,qafld) )
        disp([rhfld '->' qafld]);
        stn = station_relhumid_to_spechumid(stn,afld,rhfld,qafld);
      end;
    elseif ( isfield(stn,qafld) )
      disp([qafld '->' rhfld]);
      stn = station_spechumid_to_relhumid(stn,afld,qafld,rhfld);
      disp([rhfld '->' dfld]);
      stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
    end;

    if ( ~isfield(stn,qsfld) )
      % Saturated ("sea-surface") specific humidity [kg/kg]
      stn = station_relhumid_to_spechumid(stn,sfld,100,qsfld);
      % Adjustment for salinity with thanks to Stommel
      stn.(qsfld).data = 0.98 .* stn.(qsfld).data;
    end;


    %%%
    %% Waves and Stokes Drift

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
    stn = verify_variable(stn,Ulpfld);
    stn = verify_variable(stn,Vlpfld);
    stn.(Wlpfld).date = stn.(Ulpfld).date;
    stn.(Wlpfld).data = uv_to_spd(stn.(Ulpfld).data,stn.(Vlpfld).data);
    stn.(Dlpfld).date = stn.(Ulpfld).date;
    stn.(Dlpfld).data = uv_to_dir(stn.(Ulpfld).data,stn.(Vlpfld).data);
    %%%% ??? DEBUG: Low-pass filter winds for quasi-Eulerian currents

    %%%% ??? DEBUG: Low-pass filter winds for quasi-Eulerian currents
    if ( ~isfield(stn,ssufld) )
      % stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wfld,Dfld,whfld,wpfld,wdfld);
      stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wlpfld,Dlpfld,whfld,wpfld,wdfld);
    end;
    stn.commentstr = [stn.commentstr ' (QE~72hlp) '];
    %%%% ??? DEBUG: Low-pass filter winds for quasi-Eulerian currents


    %%%
    %% Kilometer-scale Ocean Data

    if ( ~isfield(stn,ufld) )
      switch (KMPFX),
       case 'fkeys_hycom',
        %function stn = get_fkeys_hycom(stn_or_stnm,mindt,maxdt,vars,flds,interpMethod,fkeyspath)
        stn = get_fkeys_hycom(stn,[],[],[],[],'linear');
       case 'gom_hycom',
        %function stn = get_gom_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
        stn = get_gom_hycom(stn);
       case 'avhrr_weekly',
        if ( ~isfield(stn,Tfld) )
          disp('Loading AVHRR_WEEKLY SST instead of hydrodynamic model data...');
          stn = get_avhrr_weekly_field(stn,true);
        end;
        opts.km_scale_advection = false;
        disp('ONLY STOKES');
       case 'none',
        disp('Loading NO hydrodynamic model data...');
        opts.km_scale_advection = false;
        disp('ONLY STOKES');
       otherwise,
        error('Unknown km-scale data source "%s"',KMPFX);
      end;
      more off;
    end;



    %%%
    %% Near-bottom (km + tide) Currents

    opts.km_scale_advection = get_opt(opts,'km_scale_advection',true);
    if ( ~opts.km_scale_advection )
      %%%% ??? DEBUG: Assume no km-scale advection at all
      stn.(ufld) = stn.(tufld);    stn.(ufld).data(:) = 0;
      stn.(vfld) = stn.(tvfld);    stn.(vfld).data(:) = 0;
      stn.(hufld) = stn.(tufld);   stn.(hufld).data(:) = 0;
      stn.(hvfld) = stn.(tvfld);   stn.(hvfld).data(:) = 0;
      stn.commentstr = [stn.commentstr ' (NO KM-ADVEC) '];

    else
      % Spline-fit an hourly time series of mean currents to native data
      stn.(hufld) = interp_ts(stn.(ufld));
      stn.(hvfld) = interp_ts(stn.(vfld));

    end;

    stn.(netufld) = ts_op(stn.(tufld),stn.(hufld),'+');
    stn.(netvfld) = ts_op(stn.(tvfld),stn.(hvfld),'+');

    %%%% ??? DEBUG: Low-pass filter km-scale currents for advection
    % hulpfld = [ufld '_90_day_lowpass'];
    % hvlpfld = [vfld '_90_day_lowpass'];
    % stn = verify_variable(stn,hulpfld);
    % stn = verify_variable(stn,hvlpfld);
    % stn.(netufld) = ts_op(stn.(tufld),stn.(hulpfld),'+');
    % stn.(netvfld) = ts_op(stn.(tvfld),stn.(hvlpfld),'+');
    % stn.commentstr = [stn.commentstr ' (90dLP KM-ADVEC) '];
    %%%% ??? DEBUG: Low-pass filter km-scale currents for advection


    %%%
    %% Quasi-Eulerian (km + surface) Currents

    if ( ~isfield(stn,qeufld) )
      % % stn = calc_quasi_eulerian(stn,STOKESPFX,KMPFX,QEPFX);

      stn = station_calc_quasi_eulerian(stn,ssufld,ssvfld,hufld,hvfld,QEPFX);

      % stn = station_calc_quasi_eulerian(stn,ssufld,ssvfld,netufld,netvfld,QEPFX);
      % stn.commentstr = [stn.commentstr ' (SS+=tide) '];

      % %%%% ??? DEBUG: Low-pass filter km-scale currents for advection
      % stn = station_calc_quasi_eulerian(stn,ssufld,ssvfld,hulpfld,hvlpfld,QEPFX);
    end;


    %%%
    %% Climatological surface fluxes for comparison 
    if ( ~isfield(stn,climq0fld) )
      switch (CLIMPFX),
       case 'daily_oaflux',	stn = station_load_oaflux(stn);
       otherwise,		error('Unknown flux climatology "%s"',CLIMPFX);
      end;
      stn = station_heat_flux_term(stn,climq0fld,climqtfld,sfld,[],nanmean(stn.(mhfld).data));
    end;


    %%%
    %% Radiative Fluxes

    % Follows seasonal pattern in Kd from MLRF1 PAR sensors, and Hu pers. comm.
    opts.kd = get_opt(opts,'kd',[0.1,0.25,274]);
    stn.commentstr = [stn.commentstr ' (Kd=' num2str(opts.kd,'%g,') ') '];

    if ( ~isfield(stn,asrfld) )
      % If absorbed short-wave not already present, user wants absorption calculation
      stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,opts);
    end;
    if ( ~isfield(stn,lrfld) )
      % If reanalysis long-wave flux not specified (or not available), do bulk estimate
      %station_bulk_longwave(stn,afld,qfld,pfld,dsrfld,sfld,cfld,dlrf,ulrf,lrf)
      %stn = station_bulk_longwave(stn,afld,qafld,pfld,dsrfld,sfld,cfld,dlrfld,ulrfld,lrfld);
      stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);
    end;


    %%%
    %% Turbulent and Net Fluxes

    if ( ~isfield(stn,pblzfld) )
      % Default: 600m
      pblzfld = 600;
      disp(['Using constant (Atmospheric) Planetary Boundary Layer height ' num2str(pblzfld)]);
    end;

    if ( ~isfield(stn,q0fld) )
      %station_heat_flux(stn,wfld,afld,qfld,pfld,tfld,sfld,lfld,PFX,dsfld,dlfld,prfld,
      %                  wdfld,oufld,ovfld,wpfld,whfld,pblzfld,doWarm,doPlot,max_pwp)

      % Fluxes WITH warm-layer adjustment
      % Assume max warm-layer depth somewhere near bottom boundary layer top
      max_wl = nanmax(stn.(mhfld).data) - 0.5;
      stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                              pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
                              Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,true,max_wl);
      % Algorithm sometimes returns complex numbers!
      stn.(qlhfld).data = real(stn.(qlhfld).data);
      stn.(qshfld).data = real(stn.(qshfld).data);
      stn.(qrhfld).data = real(stn.(qrhfld).data);
      stn.(q0fld).data = real(stn.(q0fld).data);

      % Cross- and long-shore components of the wind stress
      stn = station_reorient_vectors(stn,bathorifld,tauxfld,tauyfld);

      %%%% ??? DEBUG
      % Algorithm does not handle hurricanes very well!
      badix = find(-1500 > stn.(qlhfld).data | stn.(qlhfld).data > 100);
      if ( ~isempty(badix) )
        warning('Removing %d extreme QLH estimates! (Summary below)',length(badix));
        nansummary(stn.(qlhfld).data(badix)),
        stn.(q0fld) = ts_op(stn.(q0fld),stn.(qlhfld),'-');
        stn.(q0fld) = ts_op(stn.(q0fld),stn.(qshfld),'-');
        stn.(q0fld) = ts_op(stn.(q0fld),stn.(qrhfld),'-');

        stn.(qlhfld).date(badix) = [];      stn.(qlhfld).data(badix) = [];
        stn.(qshfld).date(badix) = [];      stn.(qshfld).data(badix) = [];
        stn.(qrhfld).date(badix) = [];      stn.(qrhfld).data(badix) = [];

        stn.(q0fld) = ts_op(stn.(q0fld),stn.(qlhfld),'+');
        stn.(q0fld) = ts_op(stn.(q0fld),stn.(qshfld),'+');
        stn.(q0fld) = ts_op(stn.(q0fld),stn.(qrhfld),'+');

        stn.(taufld).date(badix) = [];      stn.(taufld).data(badix) = [];
        stn.(tauxfld).date(badix) = [];     stn.(tauxfld).data(badix) = [];
        stn.(tauyfld).date(badix) = [];     stn.(tauyfld).data(badix) = [];
      end;
      %%%% ??? DEBUG


      stn.(qradfld) = ts_op(stn.(asrfld),stn.(lrfld),'+');
      stn.(qturfld) = ts_op(stn.(qlhfld),stn.(qshfld),'+');
      if ( isfield(stn,qrhfld) && is_valid_ts(stn.(qrhfld)) )
        stn.(qturfld) = ts_op(stn.(qturfld),stn.(qrhfld),'+');
      end;
      stn.(qcoolfld) = ts_op(stn.(lrfld),stn.(qturfld),'+');

      stn = station_heat_flux_term(stn,qradfld,qradtfld,sfld,[],mhfld);
      stn = station_heat_flux_term(stn,qturfld,qturtfld,sfld,[],mhfld);
      stn = station_heat_flux_term(stn,qcoolfld,qcooltfld,sfld,[],mhfld);
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],mhfld);

      % Fluxes without warm-layer adjustment - for comparison
      stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                              pfld,sfld,asrfld,lrfld,TUR30PFX,dsrfld,dlrfld,rfld,...
                              Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,false);
      % Algorithm sometimes returns complex numbers!
      stn.(qlh30fld).data = real(stn.(qlh30fld).data);
      stn.(qsh30fld).data = real(stn.(qsh30fld).data);
      stn.(qrh30fld).data = real(stn.(qrh30fld).data);
      stn.(q030fld).data = real(stn.(q030fld).data);
      stn = station_heat_flux_term(stn,q030fld,qt30fld,sfld,[],mhfld);
      stn.(qtur30fld) = ts_op(stn.(qlh30fld),stn.(qsh30fld),'+');

      % Net flux without absorption calculation or benthic flux - for comparison
      stn.(sqradfld) = ts_op(stn.(srfld),stn.(lrfld),'+');
      stn.(sq0fld) = ts_op(stn.(sqradfld),stn.(qturfld),'+');
      stn.(sq030fld) = ts_op(stn.(sqradfld),stn.(qtur30fld),'+');
      stn = station_heat_flux_term(stn,sq0fld,sqtfld,sfld,[],mhfld);
      stn = station_heat_flux_term(stn,sq030fld,sqt30fld,sfld,[],mhfld);
    end;


    %%%
    %% Save results to MAT file for future reference
    disp(['Saving to ' matfname]);
    station = stn;
    save(matfname,'station');
    station = []; clear station;
  end; %if ( exist(matfname,'file') ) else

  % EXPERIMENTING: Regardless of saved MAT file, ALWAYS redo these calculations
  if (1)

    % Check result if we ignore warm-layer adjustment
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',true);
    if ( ~opts.do_warm_layer )
      q0fld = q030fld;
      qtfld = qt30fld;
      sq0fld = sq030fld;
      sqtfld = sqt30fld;
      stn.commentstr = [stn.commentstr ' (NO WARM) '];
    end;

    % Adjust latent heat flux by a fixed factor, as in Monismith et al (2006)
    qlhadj = 1.00;
    if ( qlhadj ~= 1.00 )
      %%%% ??? DEBUG
      stn.commentstr = [stn.commentstr ' Q_L_H*' num2str(qlhadj) ' '];
      stn.(adj_qlhfld).date = stn.(qlhfld).date;
      stn.(adj_qlhfld).data = stn.(qlhfld).data .* qlhadj;
      [swix,lwix,lhix,shix,rhix] = ...
          intersect_all_dates([],stn.(asrfld).date,stn.(lrfld).date,...
                              stn.(adj_qlhfld).date,stn.(qshfld).date,stn.(qrhfld).date);
      stn.(q0fld).date = stn.(qlhfld).date(lhix);
      stn.(q0fld).data = stn.(asrfld).data(swix) + stn.(lrfld).data(lwix) ...
          + stn.(adj_qlhfld).data(lhix) + stn.(qshfld).data(shix) + stn.(qrhfld).data(rhix);
      % Algorithm sometimes returns complex numbers!
      stn.(q0fld).data = real(stn.(q0fld).data);
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],mhfld);
    end;


    %%%
    %% Eulerian (km-scale + Stokes) Heat Advection

    switch (KMPFX),
     case 'avhrr_weekly',
      opts.gradient_climatology = [];
      % opts.laplacian_climatology = [];
      disp('ONLY STOKES, AVHRR_WEEKLY SST GRADIENTS, and FIELD LAPLACIAN');
     case 'none',
      disp('ONLY STOKES');
    end;

    opts.add_alongshore_advection = get_opt(opts,'add_alongshore_advection',true);
    opts.gradient_climatology = get_opt(opts,'gradient_climatology',[]);

    opts.grid_interp_method = get_opt(opts,'grid_interp_method','linear');
    disp(['Field gradient interpolation method ' opts.grid_interp_method]);

    if ( ~isfield(stn,Tfld) )
      warning('No sea-surface temperature field %s',Tfld);
    elseif ( ~isfield(stn.(Tfld),'gradient_x') || ~isfield(stn,kmtxfld) )
      % Calculate gradients and field Laplacians - and interpolate to site
      stn = calc_field_terms(stn,Tfld,kmtfld,opts.grid_interp_method);
    end;

    if ( isempty(opts.gradient_climatology) )
      if ( opts.add_alongshore_advection )
        % Include both vector components of model heat advection
        stn = station_calc_udotdelt(stn,qeufld,qevfld,Tfld,kmtfld,...
                                    rawudTfld,udTfld,...
                                    qtfld,qtAdvfld,opts.grid_interp_method);
      else
        % Ignore the along-shore component of model heat advection
        stn.commentstr = [stn.commentstr ' ^.XSHORE (' num2str(stn.(bathorifld)) 'o) '];
        stn = station_cross_shore_advection(stn,bathorifld,...
                                            qeufld,qevfld,Tfld,kmtfld,...
                                            rawudTfld,udTfld,...
                                            qtfld,qtAdvfld,opts.grid_interp_method);
      end; %if ( opts.add_alongshore_advection ) else

    else
      stn.commentstr = [stn.commentstr ' \partial_x_sT (' num2str(opts.gradient_climatology,'%g ') ') '];

      %DEBUG:
      grclim = [hkmtxsfld '_clim'];

      stn.(grclim) = stn.(qtfld);
      stn.(grclim).data(:) = build_clim_opt(opts.gradient_climatology,grclim,stn.(grclim).date);

      if ( ~isfield(stn,qexsfld) )
        stn = station_reorient_vectors(stn,bathorifld,qeufld,qevfld,qexsfld,qelsfld);
      end;

      stn.(udTfld) = ts_op(stn.(qexsfld),stn.(grclim),'*');
      stn.(qtAdvfld) = ts_op(stn.(udTfld),stn.(qtfld),'+');

    end; %if ( isempty(opts.gradient_climatology) ) else


    %%%
    %% Km-scale Heat Diffusion

    opts.K_theta = get_opt(opts,'K_theta',model_K_theta);
    if ( isnumeric(opts.K_theta) )
      stn.commentstr = [stn.commentstr ' K_\theta (' num2str(opts.K_theta,'%g ') ') '];
    elseif ( all(isfield(opts.K_theta,{'func','arg'})) )
      stn.commentstr = [stn.commentstr ' K_\theta~' strrep(opts.K_theta.arg,'_','\_') ' '];
      opts.K_theta = opts.K_theta.func(stn.(opts.K_theta.arg));
    end;

    opts.laplacian_climatology = get_opt(opts,'laplacian_climatology',[]);

    if ( isempty(opts.laplacian_climatology) )

      %DEBUG:      stn.commentstr,
      stn = station_calc_kdel2t(stn,opts.K_theta,Tfld,...
                                rawkd2Tfld,kd2Tfld,...
                                qtAdvfld,dTfld,opts.grid_interp_method);
      stn = station_heat_flux_term_inverse(stn,kd2Tffld,kd2Tfld,sfld,[],mhfld);

    else
      stn.commentstr = [stn.commentstr ' \nabla^2T (' num2str(opts.laplacian_climatology,'%g ') ') '];

      stn.(hkmtlfld) = stn.(qtfld);
      stn.(hkmtlfld).data(:) = build_clim_opt(opts.laplacian_climatology,hkmtlfld,stn.(hkmtlfld).date);

      k = build_clim_opt(opts.K_theta,'opts.K_theta',stn.(hkmtlfld).date);

      stn.(kd2Tfld).date = stn.(hkmtlfld).date;
      stn.(kd2Tfld).data = 3600 .* k .* stn.(hkmtlfld).data;
      stn.(dTfld) = ts_op(stn.(kd2Tfld),stn.(qtAdvfld),'+');

    end; %if ( isempty(opts.laplacian_climatology) ) else

    stn = station_heat_flux_term_inverse(stn,dTffld,...
                                         dTfld,sfld,[],mhfld);



    %%%
    %% Benthic Heat Exchanges

    stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld,opts);
    stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],mhfld);

    stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'-');
    stn = station_heat_flux_term(stn,bq0fld,bq0tfld,sfld,[],mhfld);

    stn.(bdTfld) = ts_op(stn.(dTfld),stn.(qbotfld),'-');
    stn = station_heat_flux_term_inverse(stn,bdTffld,bdTfld,sfld,[],mhfld);


    %%%
    %% Quality Control
    %% Budget should skip gaps in sea temperature time series
    stn = filter_gaps(stn,sfld,sq0fld,(5/24),1);
    stn = filter_gaps(stn,sfld,q0fld,(5/24),1);
    stn = filter_gaps(stn,sfld,dTfld,(5/24),1);
    stn = filter_gaps(stn,sfld,dTffld,(5/24),1);
    stn = filter_gaps(stn,sfld,bq0fld,(5/24),1);
    stn = filter_gaps(stn,sfld,bdTfld,(5/24),1);
    stn = filter_gaps(stn,sfld,bdTffld,(5/24),1);


    %%%
    %% Horizontal Convection

    bet = stn.(slopefld);
    % stn = station_horizontal_convection(stn,sfld,[],mhfld,tspdfld,bq0fld,bet,opts,dts);

    % qstr=q0fld;
    % qstr=sq0fld;
    % qstr=dTffld;
    qstr=bdTffld;

    stn.q0 = stn.(qstr);

    % dTlpfld = ['q0'];				hc_basis='';
    % dTlpfld = ['q0_12_hour_lowpass'];		hc_basis=':12hlp';
    dTlpfld = ['q0_24_hour_lowpass'];		hc_basis=':24hlp';
    % dTlpfld = ['q0_24_hour_average'];		hc_basis=':24hma';
    % dTlpfld = ['q0_30_hour_lowpass'];		hc_basis=':30hlp';
    % dTlpfld = ['q0_30_hour_average'];		hc_basis=':30hma';
    % dTlpfld = ['q0_40_hour_lowpass'];		hc_basis=':40hlp';
    % dTlpfld = ['q0_72_hour_lowpass'];		hc_basis=':72hlp';

    stn = verify_variable(stn,dTlpfld);
    qstr=[qstr hc_basis];

    [tix,hix,tspdix,aix,Wix,sq0ix,bq0ix,dTix] = ...
      intersect_all_dates([],stn.(sfld).date,stn.(mhfld).date,stn.(tspdfld).date,stn.(afld).date,stn.(Wfld).date,stn.(sq0fld).date,stn.(bq0fld).date,stn.(dTlpfld).date);
    dts = stn.(sfld).date(tix);
    t = stn.(sfld).data(tix);
    s = repmat(36,size(stn.(sfld).data(tix)));
    h = stn.(mhfld).data(hix);
    tspd = stn.(tspdfld).data(tspdix);
    at = stn.(afld).data(aix);
    W = stn.(Wfld).data(Wix);
    sq0 = stn.(sq0fld).data(sq0ix);
    bq0 = stn.(bq0fld).data(bq0ix);
    dT = stn.(dTlpfld).data(dTix);

    % q0lp = stn.(q0lpfld).data(q0lpix);
    % Wlp = stn.(Wlpfld).data(Wlpix);


    q=dT;

    opts.hc_debug = get_opt(opts,'hc_debug',true);
    opts.R = get_opt(opts,'hc_R',(1.00-0.08));

    % opts.hc_scaling = get_opt(opts,'hc_scaling','Seasonal');
    % opts.hc_scaling = get_opt(opts,'hc_scaling','SS');
    opts.hc_scaling = get_opt(opts,'hc_scaling','US');
    % opts.hc_scaling = get_opt(opts,'hc_scaling','SU');
    % opts.hc_scaling = get_opt(opts,'hc_scaling','UU');
    maxhrs = 12;
    opts.hc_max_onset_secs = get_opt(opts,'hc_max_onset_secs',maxhrs*3600);
    %%%% ??? DEBUG
    stn.commentstr = [stn.commentstr ' (HC' hc_basis ' ' opts.hc_scaling ' maxT: ' num2str(maxhrs) 'h) '];
    % res = horizontal_convection(t,s,h,q,bet,opts,dts,q0lp,Wlp);
    res = horizontal_convection(t,s,h,q,bet,opts,dts,q,W);

    flds = fieldnames(res);
    for fldix = 1:length(flds)
      fld = flds{fldix};
      dat = res.(fld);
      res.(fld) = [];
      res.(fld).date = dts;
      res.(fld).data = dat;

      stnfld = [HCPFX '_' fld];
      stn.(stnfld).date = dts;
      stn.(stnfld).data = dat;
    end;
    stn = station_heat_flux_term_inverse(stn,hcdTdtf,hcdTdt,sfld,[],mhfld);
    stn = station_heat_flux_term_inverse(stn,hcdTdthcf,hcdTdthc,sfld,[],mhfld);


    %%%
    %% Calculate Sub-Grid Scale Heat Diffusion (as a residual)
    if ( ~isfield(stn,Tfld) )
      warning('No sea-surface temperature field for SGS diffusion ("%s")',Tfld);
    else
      stn = station_calc_sgs_diffusion(stn,sfld,hcdTdt,Tfld,hkmtlfld,...
                                       sgskd2Tfld,sgskfld,sgsdTdt);

      % sgsklp = [sgskfld '_30_day_lowpass'];
      % stn = verify_variable(stn,sgsklp);
      % fmg;
      % plot_ts(stn.(sgskfld),stn.(sgsklp));
      % legend('K_\theta_s_g_s','K_\theta_s_g_s 30dLP');
    end;


    %%%
    %% Plot several diagnostic figures

    %%%% ??? DEBUG
    % Truncate start of all time series for display purposes
    % % Start date of FKEYS HYCOM
    % startdt = datenum(2004,1,1);
    % % Start date of GOM HYCOM
    % startdt = datenum(2003,1,1);
    % % Start date of WW3
    % startdt = datenum(1999,7,1);
    % Start date of in situ data
    startdt = stn.(bdTfld).date(1);
    finaldt = stn.(bdTfld).date(end);

    startix = find(dts>=startdt,1);

    %%%% DEBUG???
    ix = find(~isfinite(sq0(startix:end)),1,'last');
    if ( ~isempty(ix) )
      startix = startix + ix;
    end;

    dts = dts(startix:end);
    t = t(startix:end);
    h = h(startix:end);
    tspd = tspd(startix:end);
    at = at(startix:end);
    W = W(startix:end);
    sq0 = sq0(startix:end);
    bq0 = bq0(startix:end);
    dT = dT(startix:end);
    fac = stn.(hcfactor).data(startix:end);
    dTdt.date = stn.(hcdTdt).date(startix:end);
    dTdt.data = stn.(hcdTdt).data(startix:end);

    bt = stn.(btfld);
    [ig,startix] = min(abs(bt.date-dts(1)));
    [ig,endix] = min(abs(bt.date-dts(end)));
    bt.date = bt.date(startix:endix);
    bt.data = bt.data(startix:endix);
    %%%% ??? DEBUG

    dsr_lpfld = [asrfld '_24_hour_sum'];
    stn = verify_variable(stn,dsr_lpfld);
    [ig,dsr_lpix] = intersect_dates(dts,stn.(dsr_lpfld).date);

    % lhf_lpfld = [qlhfld '_24_hour_sum'];
    % stn = verify_variable(stn,lhf_lpfld);
    % [ig,lhf_lpix] = intersect_dates(dts,stn.(lhf_lpfld).date);

    udT_lpfld = [udTfld];
    stn = verify_variable(stn,udT_lpfld);
    [ig,udT_lpix] = intersect_dates(dts,stn.(udT_lpfld).date);


%%%% ??? DEBUG
if (0)
    %%%% ??? DEBUG
    T0 = t(1);
    bigfh=figure, maxigraph; hold on;
    plot(dts,t,dts,at,'k:',bt.date,bt.data,'m--',dts,T0+cumsum(sq0.*fac),dts,T0+cumsum(bq0.*fac),dts,T0+cumsum(dT.*fac),dTdt.date,T0+cumsum(dTdt.data));
    % plot(stn.(dsr_lpfld).date(dsr_lpix),(stn.(dsr_lpfld).data(dsr_lpix)./1000),stn.(lhf_lpfld).date(lhf_lpix),-(stn.(lhf_lpfld).data(lhf_lpix)./1000),'y:','Color',[.8,.8,0]);
    % plot(stn.(dsr_lpfld).date(dsr_lpix),(stn.(dsr_lpfld).data(dsr_lpix)./1000),stn.(udT_lpfld).date(udT_lpix),T0+cumsum(stn.(udT_lpfld).data(udT_lpix)),'y:','Color',[.8,.8,0]);

%%%% DEBUG???
    % plot(stn.(dsr_lpfld).date(dsr_lpix),T0+(stn.(dsr_lpfld).data(dsr_lpix).*fac),stn.(udT_lpfld).date(udT_lpix),T0+cumsum(stn.(udT_lpfld).data(udT_lpix)),'y:','Color',[.8,.8,0]);
%%%% DEBUG???

    % plot(dts,h,dts,tspd,'k:',dts,W,'o','Color',[0,.8,.2]);
    % plot(dts,h,dts,tspd,'k:','Color',[0,.8,.2]);
    datetick3('x',2,'keeplimits');
    % legend('T_s','T_a','T_b','T_0+Q_0/\rhoC_ph','T_0+(Q_0(\gamma)+Q_b)/\rhoC_ph','T_0+\partial_tT_k_m',['T_0+\partial_tT ' opts.hc_scaling],'T_0+\Sigma_1_d Q_S_W\times1h/\rho^.C_p^.h','T_0+\Sigmau^.\nablaT_k_m','h_t_i_d_e','SPD_t_i_d_e','W', 'Location','SouthWest'); %'Best');
    % legend('T_s','T_a','T_b','T_0+Q_0/\rhoC_ph','T_0+(Q_0(\gamma)+Q_b)/\rhoC_ph','T_0+\partial_tT_k_m',['T_0+\partial_tT ' opts.hc_scaling],'T_0+\Sigma_1_d Q_S_W\times1h/\rho^.C_p^.h','T_0+\Sigmau^.\nablaT_k_m','h_t_i_d_e','SPD_t_i_d_e', 'Location','SouthWest'); %'Best');
    legend('T_s','T_a','T_b','T_0+Q_0/\rhoC_ph','T_0+(Q_0(\gamma)+Q_b)/\rhoC_ph','T_0+\partial_tT_k_m',['T_0+\partial_tT ' opts.hc_scaling],'T_0+\Sigma_1_d Q_S_W\times1h/\rho^.C_p^.h','T_0+\Sigmau_q_e^.\nablaT_k_m', 'Location','SouthWest'); %'Best');
    titlename([ stn.commentstr stn.station_name ' ' strrep(qstr,'_','\_') ]);
    disp(stn.commentstr);
    grid on;
    % axis('tight');
    % axis('tight'); ylim([-50 150]);
    % axis([startdt,startdt+365,0,60]);
    % axis([startdt,startdt+93,0,26]);
    % axis([startdt,startdt+5,0,26]);
    % axis([startdt,datenum(2008,12,31),-100,100]);
    axis([startdt,finaldt,0,60]);
    %%%% ??? DEBUG
    datetick3('x',2,'keeplimits');
    figfname= fullfile(figspath,[lower(stn.station_name) '_' hcdTdt '.']);
    % print('-dtiff',[figfname 'tiff']);
    % print('-dpng',[figfname 'png']);
end;

    % % Plot daily climatology comparisons
    % compare_flux_climatologies(stn,'daily',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,stn.commentstr);
    % figfname= fullfile(figspath,[lower(stn.station_name) '_' hcdTdt '_chkann.']);
    % % print('-dtiff',[figfname 'tiff']);
    % % print('-dpng',[figfname 'png']);


    % [stn,hls,hxs]=multiplot_station(stn,{'absorbed_erai_srf','ndbc_erai_30a_latent_flux','ndbc_erai_30a_net_flux','benthic_erai_qbo','benthic_ndbc_erai_30a_net_flux','ndbc_air_t','ndbc_sea_t'},[],[],{'SW','LH','Q_0','BH','Q_0+Q_b','T_a','T_s'},[datenum(2005,4,1),datenum(2005,5,1)],{[-1000,1000],[-1000,1000],[-1000,1000],[-1000,1000],[-1000,1000],[16,28],[16,28]});
    % axes(hxs(1));
    % hold on; plot(stn.erai_srf.date,stn.erai_srf.data,'ro','MarkerSize',1.5);


%%%% ??? DEBUG
% warning('All graphs commented out!');
if (1)
    % Simple year-by-year comparison (SUBPLOTs) of truth and heat budget
    % annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,stn.commentstr,6:9,2002,2010,'ndbc_sea_t');
    annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,stn.commentstr,[],1993);
    subplots_set('ylim',[5,45]);
    figfname= fullfile(figspath,[lower(stn.station_name) '_' hcdTdt '_annsubs.']);
    % print('-dtiff',[figfname 'tiff']);
    % print('-dpng',[figfname 'png']);
end;

%%%% ??? DEBUG
if (1)
    ms1_clim(stn,'daily',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,stn.commentstr,[],[],[],1996);
    figfname= fullfile(figspath,[lower(stn.station_name) '_' hcdTdt '_weekly.']);
    % % print('-dtiff',[figfname 'tiff']);
    % % print('-dpng',[figfname 'png']);
end;

%%%% ??? DEBUG
if (0)
    begyr=1998; endyr=1998; mos=1:6;
    annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,stn.commentstr,mos,begyr,endyr);
    T0 = stn.(sfld).data(find(stn.(sfld).date>=datenum(begyr,mos(1),1),1));
    advPlotFld = udTfld;
    % advPlotFld = qtAdvfld;
    begix = find(stn.(advPlotFld).date>=datenum(begyr,mos(1),1),1);
    plot(stn.(advPlotFld).date(begix:end),...
         T0+cumsum(stn.(advPlotFld).data(begix:end)),'m-.');
    ylim([6,30]);
end;

if (0)
    % fmg; plot_ts(stn.(q0fld),stn.(bdTffld),stn.(qturfld),stn.(qradfld),stn.([hcdTdt,'hc_flux'])); legend('Q_0','Q_0+Q_b','Tur','Rad','HC'); xlim(datenum(1998,[1,6],1)); datetick3('x',26,'keeplimits');
    fmg; plot_ts(stn.(qtfld),stn.(bdTfld),stn.(qturtfld),stn.(qradtfld),stn.([hcdTdt,'hc']),stn.(sssst)); legend('Q_0','Q_0+Q_b','Tur','Rad','HC','u^.\nablaT'); xlim(datenum(1998,[1,6],1)); datetick3('x',26,'keeplimits');
    titlename(strrep([qstr],'_','\_'));
end;

    %%%
    %% Save results to MAT file for future reference
    disp(['NOT RE-Saving to ' matfname]);
    % station = stn;
    % save(matfname,'station');
    % station = []; clear station;

  end; %IF(1) %if ( exist(matfname,'file') ) else


  toc,
  set_more;


  % stn = compare_heat_budgets(stn,sfld,RAPFX,KMPFX);
  % ylim([5 35]);
  % figfname= fullfile(figspath,[lower(stn.station_name) '_' dTfld '.tiff']);
  % print('-dtiff',figfname);

  if (nargout < 1)
    stn = [];    clear stn;
  end;

return;
