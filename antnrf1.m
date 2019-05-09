function stn = antnrf1(stn,doBudget)
%function stn = antnrf1(stn,doBudget)
%
% Analyze all data for the "TNRF1" joint GLERL/FIO/AOML experimental reef
% monitoring station, and calculate heat fluxes and comparisons from it. If
% optional arg DOBUDGET is True (DEFAULT), calculate heat budget for TNRF1.
%
% Last Saved Time-stamp: <Wed 2011-07-06 06:25:20  lew.gramer>


  set_more off;

  datapath = get_thesis_path('../data');
  tnrf1path = fullfile(datapath,'tnrf1');

  if ( ~exist('doBudget','var') || isempty(doBudget) )
    doBudget = true;
  end;

  matfname = fullfile(datapath,'tnrf1.mat');

  if ( exist(matfname,'file') )
    disp(['Loading presave MAT file ' matfname]);
    x=load(matfname,'stn');
    flds = fieldnames(x.stn);
    for ix=1:length(flds)
      fld = flds{ix};
      stn.(fld) = x.stn.(fld);
    end;
    x=[]; clear x;


  else

    if ( ~isfield(stn,'station_name') )
      stn.station_name = 'TNRF1';
    end;
    if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
      [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
    end;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Proxy longer (QC'd) met. from LONF1,SMKF1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Using regional sites as proxies for some meteorological data');
    lonf1 = load_all_ndbc_data([],'lonf1');
    flds = grepstruct(lonf1,'ndbc_');
    for fld=flds; stn.(fld{:}) = lonf1.(fld{:}); end;
    lonf1 = []; clear lonf1;

    smkf1 = load_all_ndbc_data([],'smkf1');
    % One-hour lag seems to give best coherence and fit!
    stn.ndbc_dew_t.date = smkf1.ndbc_dew_t.date(1:end-1);
    stn.ndbc_dew_t.data = smkf1.ndbc_dew_t.data(2:end);
    smkf1 = []; clear smkf1;

    stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');

    bbfname = fullfile(datapath,'tnrf1.bb');
    if ( exist(bbfname,'file') )
      disp(sprintf('Loading "fact" data file "%s"', bbfname));
      stn = load_bb_data(bbfname, stn);
    end;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load in situ WXT and ADCP (and NXIC?) data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clear any previouus ADCP data in STN
    stn.adcp_u.date = []; stn.adcp_u.data = [];
    stn.adcp_v.date = []; stn.adcp_v.data = [];
    stn.adcp_good.date = []; stn.adcp_good.data = [];


    % Get all known dates of in situ data from TNRF1

    for dt = [ datenum(2009,7,20):datenum(2009,9,09) ...
               datenum(2010,5,27):floor(now) ]
      %DEBUG:    disp(datestr(dt));


      % WXT520 data

      datfname = fullfile(tnrf1path,'metwxtaver',['metwxtaver' datestr(dt,'yyyymmdd') '.dat']);
      warning('off','ImportMetWxtAver:NoFile');
      result = importmetwxtaver(datfname);
      warning('on','ImportMetWxtAver:NoFile');
      if ( ~isempty(result) )
        stn = merge_station_data(stn,result);
      end;
      result = []; clear result;



      % ADCP data

      datfname = fullfile(tnrf1path,'adcpcurr',['adcpcurr' datestr(dt,'yyyymmdd') '.dat']);

      % Stupid IMPORTDATA does not handle zero-size files well
      fexists = false;
      if ( exist(datfname,'file') )
        fid = fopen(datfname,'r');
        if ( fid > 0 )
          fread(fid,1);
          fexists = ~feof(fid);
          fclose(fid);
        end;
      end;

      if ( fexists )
        dat = importadcp(datfname);
      else
        % warning('No ADCPCURR file "%s" found.',datfname);
        continue;
      end;

      result.adcp_i_depth.date = dat.currprof.date;
      result.adcp_i_depth.data = dat.currprof.i_depth;

      result.adcp_seatemp.date = dat.currprof.date;
      result.adcp_seatemp.data = dat.currprof.seatemp;

      result.adcp_barotropic_u.date = dat.currprof.date;
      result.adcp_barotropic_u.data = dat.currprof.u;

      result.adcp_barotropic_v.date = dat.currprof.date;
      result.adcp_barotropic_v.data = dat.currprof.v;

      stn = merge_station_data(stn,result);
      result = []; clear result;

      % Have to do profiles by hand
      nprofs = length(dat.currprof.date);
      nlvls = size(dat.currprof.uprof,2);
      stn.adcp_u.date(end+1:end+nprofs,1) = dat.currprof.date;
      stn.adcp_u.data(end+1:end+nprofs,1:nlvls) = dat.currprof.uprof;
      stn.adcp_v.date(end+1:end+nprofs,1) = dat.currprof.date;
      stn.adcp_v.data(end+1:end+nprofs,1:nlvls) = dat.currprof.vprof;
      stn.adcp_good.date(end+1:end+nprofs,1) = dat.currprof.date;
      stn.adcp_good.data(end+1:end+nprofs,1:nlvls) = dat.currprof.profgood;

      dat = []; clear dat;
    end;

    % Do some bin averaging on ADCP data
    %DEBUG:  maxigraph(figure); plot(stn.adcp_v.date,[nanmean(stn.adcp_v.data(:,[19:21]),2) , nanmean(stn.adcp_v.data(:,[36:38]),2) , stn.adcp_i_depth.data-5.05]); datetick3; legend('V19-21','V36-38','H');

    highbins = 36:38;
    goods = nanmean(stn.adcp_good.data(:,[highbins]),2);
    goodix = find(goods > 50);
    stn.adcp_high_u.date = stn.adcp_u.date(goodix);
    stn.adcp_high_u.data = nanmean(stn.adcp_u.data(goodix,[highbins]),2);
    stn.adcp_high_v.date = stn.adcp_v.date(goodix);
    stn.adcp_high_v.data = nanmean(stn.adcp_v.data(goodix,[highbins]),2);
    stn = verify_variable(stn,'adcp_high_u_40_hour_lowpass');
    stn = verify_variable(stn,'adcp_high_v_40_hour_lowpass');

    midbins = 19:21;
    goods = nanmean(stn.adcp_good.data(:,[midbins]),2);
    goodix = find(goods > 50);
    stn.adcp_mid_u.date = stn.adcp_u.date(goodix);
    stn.adcp_mid_u.data = nanmean(stn.adcp_u.data(goodix,[midbins]),2);
    stn.adcp_mid_v.date = stn.adcp_v.date(goodix);
    stn.adcp_mid_v.data = nanmean(stn.adcp_v.data(goodix,[midbins]),2);
    stn = verify_variable(stn,'adcp_mid_u_40_hour_lowpass');
    stn = verify_variable(stn,'adcp_mid_v_40_hour_lowpass');

    lowbins = 5:7;
    goods = nanmean(stn.adcp_good.data(:,[lowbins]),2);
    goodix = find(goods > 50);
    stn.adcp_low_u.date = stn.adcp_u.date(goodix);
    stn.adcp_low_u.data = nanmean(stn.adcp_u.data(goodix,[lowbins]),2);
    stn.adcp_low_v.date = stn.adcp_v.date(goodix);
    stn.adcp_low_v.data = nanmean(stn.adcp_v.data(goodix,[lowbins]),2);
    stn = verify_variable(stn,'adcp_low_u_40_hour_lowpass');
    stn = verify_variable(stn,'adcp_low_v_40_hour_lowpass');


    [ix1,ix2] = intersect_dates(stn.adcp_low_u_40_hour_lowpass.date,stn.adcp_high_u_40_hour_lowpass.date);
    stn.adcp_shear_u_40_hour_lowpass.date = stn.adcp_low_u_40_hour_lowpass.date(ix1);
    stn.adcp_shear_u_40_hour_lowpass.data = stn.adcp_low_u_40_hour_lowpass.data(ix1) ...
        - stn.adcp_high_u_40_hour_lowpass.data(ix2);
    stn.adcp_shear_v_40_hour_lowpass.date = stn.adcp_low_v_40_hour_lowpass.date(ix1);
    stn.adcp_shear_v_40_hour_lowpass.data = stn.adcp_low_v_40_hour_lowpass.data(ix1) ...
        - stn.adcp_high_v_40_hour_lowpass.data(ix2);


    disp(['Saving to MAT file ' matfname]);
    save(matfname,'stn');

  end; %if ( exist(matfname,'file') ) else


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Load model, reanalysis, and derived data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ( ~isfield(stn,'ncep_lrf') )
    stn = get_ncep_station(stn, 'narr');
  end;

  if ( ~isfield(stn,'ww3_peakwavedir') )
    stn = get_ww3_station(stn);
  end;

  if ( ~isfield(stn,'stokes_drift_speed') )
    stn = ...
        station_stokes_drift(stn,...
                             'stokes_drift_speed','stokes_drift_dir',...
                             'stokes_drift_u','stokes_drift_v',...
                             'wxt_wspeed','wxt_wdir',...
                             'ww3_sigwavehgt','ww3_peakwaveper','ww3_peakwavedir');
  end;
  if ( ~isfield(stn,'global_hycom_u') )
    stn = load_global_hycom(stn);
  end;

  if ( ~isfield(stn,'quasi_eulerian_u') )
    stn = calc_quasi_eulerian(stn);
  end;


  UdotdelT = 'advected_heat';
  if ( ~isfield(stn,UdotdelT) )
    stn = get_avhrr_weekly_field(stn);

    stn = station_advect_field(stn,UdotdelT,...
                               'quasi_eulerian_u','quasi_eulerian_v',...
                               'avhrr_weekly_sst');
  end;



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate air-sea heat fluxes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % stn = station_relhumid_to_spechumid(stn,'ncep_air_t','ncep_relhumid','ncep_spechumid');
  % stn = station_relhumid_to_spechumid(stn,'wxt_air_t','wxt_relhumid','wxt_spechumid');

  if ( doBudget )

    more off;

    stn = ...
        station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid',...
                          'ndbc_barom','adcp_seatemp','ncep_srf','ncep_lrf',...
                          'ndbc_ncep_30a','ncep_dsrf','ncep_dlrf','ncep_precip',...
                          'ndbc_wind1_dir','default','default','default','default','default',true);
                          % 'ndbc_wind1_dir','default','default','ww3_peakwaveper','ww3_sigwavehgt',true);
    stn = station_heat_flux_term(stn,'ndbc_ncep_30a_net_heat_flux', ...
                                 'ndbc_ncep_30a_heat_flux_term', ...
                                 'adcp_seatemp',[],'adcp_i_depth');


    stn = ...
        station_heat_flux(stn,'wxt_wspeed','wxt_air_t','wxt_relhumid',...
                          'wxt_barom','adcp_seatemp','ncep_srf','ncep_lrf',...
                          'wxt_ncep_30a','ncep_dsrf','ncep_dlrf','wxt_precip',...
                          'wxt_wdir','default','default','default','default','default',true);
                          % 'wxt_wdir','default','default','ww3_peakwaveper','ww3_sigwavehgt',true);
    stn = station_heat_flux_term(stn,'wxt_ncep_30a_net_heat_flux', ...
                                 'wxt_ncep_30a_heat_flux_term', ...
                                 'adcp_seatemp',[],'adcp_i_depth');


    stn = ...
        station_heat_flux(stn,'ncep_wind_speed','ncep_air_t','ncep_relhumid',...
                          'ncep_barom','adcp_seatemp','ncep_srf','ncep_lrf',...
                          'adcp_ncep_30a','ncep_dsrf','ncep_dlrf','ncep_precip',...
                          'ncep_wind_dir','default','default','default','default','default',true);
                          % 'ncep_wind_dir','default','default','ww3_peakwaveper','ww3_sigwavehgt',true);
    stn = station_heat_flux_term(stn,'adcp_ncep_30a_net_heat_flux', ...
                                 'adcp_ncep_30a_heat_flux_term', ...
                                 'adcp_seatemp',[],'adcp_i_depth');

  end; %if doBudget

  set_more;

return;
