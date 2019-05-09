function stn = get_erai_station_2(stn_or_stnms,forceExtract,yrs)
%function stn = get_erai_station_2(stn_or_stnms,forceExtract,yrs)
%
% Extract subsets of reanalysis and forecast data from European Centre for
% Medium-range Weather Forecasting (ECMWF) Interim Reanalysis "ERA-Interim",
% for a station or list of station names STN_OR_STNMS, for 1987-present. If
% struct STN is given, it is returned with new 3-hourly fields .raw_erai_*
% and hourly interpolated fields .erai_* added for all variables; if a name
% string or cellstr of names is given instead, return the last station *that
% required extraction* in STN. If all stations were previously extracted, the
% last-named is loaded from its MAT file (unless optional arg FORCEEXTRACT).
%
% For full documentation on ECMWF 'ERA Interim' forecast and reanalysis
% products and GRiB file oddities, see: http://www.ecmwf.int/research/era
%
% TO DOWNLOAD NEW DATA to the GRiB file source directory ERAPATH (see code),
% try running the script ERAPATH/get_ERA_Interim.pl *on a Linux-like host*!
% Sorry, but that is the only option offered by our friends across the Pond.
% For documentation on the dataset which that script downloads, please see:
%   http://data.ecmwf.int/data/
%
% N.B., SIDE EFFECT: Data for any stations which are extracted are saved to
% MAT files with names like [station_name '_erai.mat'] in the DATAPATH. NOTE:
% If a value is passed in for arg YRS, that side effect may case station .MAT
% files to be INCOMPLETE; if .MAT files already exist, arg YRS is IGNORED!
%
% CALLING EXAMPLES:
%
%  >> % Extract list of stations with coords known by GET_STATION_COORDS (v.)
%  >> stn = get_erai_station({'srvi2','lppr1','lciy2'});  % Return LCIY2 data
%
%  >> % Extract and return a brand new station from ERAI at a given location
%  >> stn.station_name='flower'; stn.lon=-93.62; stn.lat=27.93;
%  >> stn = get_erai_station(stn);
%
%  >> % Station was previously extracted from GRiB files: reload from MAT file
%  >> stn = get_erai_station('mlrf1');
%
%  >> % Station is (re-)extracted from GRiB files, but only for years 2000-2016
%  >> stn = get_erai_station('mlrf1',true,2000:2016);
%
% Last Saved Time-stamp: <Fri 2016-09-16 15:26:13 Eastern Daylight Time gramer>

  set_more off

  datapath = get_thesis_path('../data');
  erapath = fullfile(datapath,'ECMWF');

  if ( iscellstr(stn_or_stnms) )
    stnms = stn_or_stnms;
  elseif ( ischar(stn_or_stnms) )
    stnms = {stn_or_stnms};
  elseif ( isstruct(stn_or_stnms) )
    stn = stn_or_stnms;
    if ( ~isfield(stn,'station_name') )
      error('First arg was a struct without .station_name field!');
    end;
    stnms = {stn.station_name};
  else
    error('First arg must be a station struct, name string, or cellstr of names!');
  end;

  if ( ~exist('forceExtract','var') || isempty(forceExtract) )
    forceExtract = false;
  end;

  % Conversion factor from W/m^2 to micromole quanta/m^2.s
  % See Morel and Smith, 1974
  UMOL_PER_WATT = 4.1513469579;
  % NOTE: Dye (JGR-A, 2004) recommends 4.56 mumol/J instead!

  % Do not re-extract stations we have already saved
  badix = [];
  if ( ~forceExtract )
    for stix=1:length(stnms)
      stnm = stnms{stix};
      matfname = fullfile(datapath,[stnm '_erai.mat']);
      if ( exist(matfname,'file') )
        disp(['Found ',matfname]);
        badix(end+1) = stix;
      end;
    end;
  end;

  % Just return the LAST of the named stations we're extracting.

  % If we are not extracting ANY, just return last-named station.
  if ( length(badix) == length(stnms) )
    load(matfname,'station');
    flds = fieldnames(station);
    for ix = 1:length(flds)
      fld = flds{ix};
      stn.(fld) = station.(fld);
    end;
    station = []; clear station;

    %%%% DEBUG??? Post-correct a spline-fitting error
    stn.erai_cloud_cover.data(stn.erai_cloud_cover.data>1) = 1;
    %%%% DEBUG???

    %%%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%%
    set_more;
    return;
  end;

  tic,

  stnms(badix) = [];

  result = [];
  stn.station_name = stnms{end};
  if ( isfield(stn,'lon') && isfield(stn,'lat') )
    result(length(stnms)).lon = stn.lon;
    result(length(stnms)).lat = stn.lat;
  end;



  dds = { ...
      { 'an', ...
        { ...
            'Albedo', ... %                        true      Albedo @ surface
            'Charnock', ... %                      true      Charnock @ surface
            'Mean_sea_level_pressure', ... % Pa    true      Mean_sea_level_pressure @ surface
            'N10_metre_U_wind_component', ... % m s-1 true   N10_metre_U_wind_component @ surface
            'N10_metre_V_wind_component', ... % m s-1 true   N10_metre_V_wind_component @ surface
            'N2_metre_dewpoint_temperature', ... % K true    N2_metre_dewpoint_temperature @ surface
            'N2_metre_temperature', ... % K        true      N2_metre_temperature @ surface
            'Sea_surface_temperature', ... % K     true      Sea_surface_temperature @ surface
            'Skin_temperature', ... % K            true      Skin_temperature @ surface
            'Surface_pressure', ... % Pa           true      Surface_pressure @ surface
            'Total_cloud_cover', ... %             true      Total_cloud_cover @ surface
        },...
        { ...
            'none', ... %                        true      Albedo @ surface
            'none', ... %                      true      Charnock @ surface
            'mult 0.01', ... % Pa    true      Mean_sea_level_pressure @ surface
            'none', ... % m s-1 true   N10_metre_U_wind_component @ surface
            'none', ... % m s-1 true   N10_metre_V_wind_component @ surface
            'add -273.14', ... % K true    N2_metre_dewpoint_temperature @ surface
            'add -273.14', ... % K        true      N2_metre_temperature @ surface
            'add -273.14', ... % K     true      Sea_surface_temperature @ surface
            'add -273.14', ... % K            true      Skin_temperature @ surface
            'mult 0.01', ... % Pa           true      Surface_pressure @ surface
            'none', ... %             true      Total_cloud_cover @ surface
        }, ...
        { ...
            'albedo', ... %                        true      Albedo @ surface
            'charnock', ... %                      true      Charnock @ surface
            'barom', ... % Pa    true      Mean_sea_level_pressure @ surface
            'wind_u', ... % m s-1 true   N10_metre_U_wind_component @ surface
            'wind_v', ... % m s-1 true   N10_metre_V_wind_component @ surface
            'dew_t', ... % K true    N2_metre_dewpoint_temperature @ surface
            'air_t', ... % K        true      N2_metre_temperature @ surface
            'sea_t', ... % K     true      Sea_surface_temperature @ surface
            'sst', ... % K            true      Skin_temperature @ surface
            'barom_surf', ... % Pa           true      Surface_pressure @ surface
            'cloud_cover', ... %             true      Total_cloud_cover @ surface
        }, ...
      }, ...
      { 'fc', ...
        { ...
            'Boundary_layer_height', ... % m       true      Boundary_layer_height @ surface
            'Convective_precipitation', ... % m    true      Convective_precipitation @ surface
            'Downward_UV_radiation_at_the_surface', ... % w m-2 s true Downward_UV_radiation_at_the_surface @ surface
            'East-West_surface_stress', ... % N m-2 s true   East-West_surface_stress @ surface
            'Evaporation', ... %    m of water     true      Evaporation @ surface
            'Instantaneous_surface_heat_flux', ... % W m-2 true Instantaneous_surface_heat_flux @ surface
            'Instantaneous_X_surface_stress', ... % N m-2 true Instantaneous_X_surface_stress @ surface
            'Instantaneous_Y_surface_stress', ... % N m-2 true Instantaneous_Y_surface_stress @ surface
            'N10_metre_wind_gust', ... % m s-1     true      N10_metre_wind_gust @ surface
            'North-South_surface_stress', ... % N m-2 s true North-South_surface_stress @ surface
            'Photosynthetically_active_radiation_at_the_surface', ... % w m-2 s true Photosynthetically_active_radiation_at_the_surface @ surface
            'Runoff', ... %         m              true      Runoff @ surface
            'Sunshine_duration', ... % s           true      Sunshine_duration @ surface
            'Surface_latent_heat_flux', ... % W m-2 s true   Surface_latent_heat_flux @ surface
            'Surface_net_solar_radiation_clear_sky', ... % W m-2 s true Surface_net_solar_radiation_clear_sky @ surface
            'Surface_net_thermal_radiation_clear_sky', ... % W m-2 s true Surface_net_thermal_radiation_clear_sky @ surface
            'Surface_sensible_heat_flux', ... % W m-2 s true Surface_sensible_heat_flux @ surface
            'Surface_solar_radiation_downwards', ... % W m-2 s true Surface_solar_radiation_downwards @ surface
            'Surface_solar_radiation', ... % W m-2 s true    Surface_solar_radiation @ surface
            'Surface_thermal_radiation_downwards', ... % W m-2 s true Surface_thermal_radiation_downwards @ surface
            'Surface_thermal_radiation', ... % W m-2 s true  Surface_thermal_radiation @ surface
            'Total_precipitation', ... % m         true      Total_precipitation @ surface
        },...
        { ...
            'none', ... % m       true      Boundary_layer_height @ surface
            'deaccum_sum', ... % m    true      Convective_precipitation @ surface
            'deaccum', ... % w m-2 s true Downward_UV_radiation_at_the_surface @ surface
            'deaccum', ... % N m-2 s true   East-West_surface_stress @ surface
            'deaccum', ... %    m of water     true      Evaporation @ surface
            'mult -1', ... % W m-2 true Instantaneous_surface_heat_flux @ surface
            'none', ... % N m-2 true Instantaneous_X_surface_stress @ surface
            'none', ... % N m-2 true Instantaneous_Y_surface_stress @ surface
            'none', ... % m s-1     true      N10_metre_wind_gust @ surface
            'deaccum', ... % N m-2 s true North-South_surface_stress @ surface
            'deaccum', ... % w m-2 s true Photosynthetically_active_radiation_at_the_surface @ surface
            'deaccum', ... %         m              true      Runoff @ surface
            'deaccum', ... % s           true      Sunshine_duration @ surface
            'deaccum', ... % W m-2 s true   Surface_latent_heat_flux @ surface
            'deaccum', ... % W m-2 s true Surface_net_solar_radiation_clear_sky @ surface
            'deaccum', ... % W m-2 s true Surface_net_thermal_radiation_clear_sky @ surface
            'deaccum', ... % W m-2 s true Surface_sensible_heat_flux @ surface
            'deaccum', ... % W m-2 s true Surface_solar_radiation_downwards @ surface
            'deaccum', ... % W m-2 s true    Surface_solar_radiation @ surface
            'deaccum', ... % W m-2 s true Surface_thermal_radiation_downwards @ surface
            'deaccum', ... % W m-2 s true  Surface_thermal_radiation @ surface
            'deaccum_sum', ... % m         true      Total_precipitation @ surface
        }, ...
        { ...
            'pblz', ... % m       true      Boundary_layer_height @ surface
            'conv_precip', ... % m    true      Convective_precipitation @ surface
            'uvW', ... % w m-2 s true Downward_UV_radiation_at_the_surface @ surface
            'wind_stress_u', ... % N m-2 s true   East-West_surface_stress @ surface
            'evap', ... %    m of water     true      Evaporation @ surface
            'net_heat_flux', ... % W m-2 true Instantaneous_surface_heat_flux @ surface
            'instant_wind_stress_u', ... % N m-2 true Instantaneous_X_surface_stress @ surface
            'instant_wind_stress_v', ... % N m-2 true Instantaneous_Y_surface_stress @ surface
            'wind_gust', ... % m s-1     true      N10_metre_wind_gust @ surface
            'wind_stress_v', ... % N m-2 s true North-South_surface_stress @ surface
            'parW', ... % w m-2 s true Photosynthetically_active_radiation_at_the_surface @ surface
            'runoff', ... %         m              true      Runoff @ surface
            'sunshine_duration', ... % s           true      Sunshine_duration @ surface
            'latent_heat_flux', ... % W m-2 s true   Surface_latent_heat_flux @ surface
            'clear_sky_srf', ... % W m-2 s true Surface_net_solar_radiation_clear_sky @ surface
            'clear_sky_lrf', ... % W m-2 s true Surface_net_thermal_radiation_clear_sky @ surface
            'sensible_heat_flux', ... % W m-2 s true Surface_sensible_heat_flux @ surface
            'dsrf', ... % W m-2 s true Surface_solar_radiation_downwards @ surface
            'srf', ... % W m-2 s true    Surface_solar_radiation @ surface
            'dlrf', ... % W m-2 s true Surface_thermal_radiation_downwards @ surface
            'lrf', ... % W m-2 s true  Surface_thermal_radiation @ surface
            'precip', ... % m         true      Total_precipitation @ surface
        }, ...
      }, ...
      { 'wv', ...
        { ...
            'Mean_wave_direction', ... %  degrees  true      Mean_wave_direction @ msl
            'Mean_wave_period', ... % s            true      Mean_wave_period @ msl
            'Significant_wave_height', ... % m     true      Significant_wave_height @ msl
        },...
        { ...
            'none', ... %   true      Mean_wave_direction @ msl
            'none', ... % s            true      Mean_wave_period @ msl
            'none', ... % m     true      Significant_wave_height @ msl
        }, ...
        { ...
            'peakwavedir', ... %   true      Mean_wave_direction @ msl
            'peakwaveper', ... % s            true      Mean_wave_period @ msl
            'sigwavehgt', ... % m     true      Significant_wave_height @ msl
        }, ...
      }, ...
        };

  lats = [];
  lons = [];

  for stix=1:length(stnms)
    result(stix).station_name = stnms{stix};
  end;

  for ddix=1:length(dds)
    dd = dds(ddix);
    dd = dd{:};
    flds = dd{4};
    for fldix = 1:length(flds)
      fld = ['raw_erai_' flds{fldix}];
      for stix=1:length(stnms)
        result(stix).(fld) = struct('date',[],'data',[]);
      end;
    end;
  end;

  [thisyear,thismonth,ig] = datevec(now);
  if ( ~exist('yrs','var') || isempty(yrs) )
    yrs = 1987:thisyear;
    %DEBUG:    yrs = 2007;
    %DEBUG:    yrs = 2013;
    %DEBUG:    yrs = 2012:2013;
    %DEBUG:    yrs = 2000:2013;
  else
    warning('All MAT files will contain years: %s',num2str(yrs)); 
  end;

  for yr = yrs(:)'

    %DEBUG:
    disp(yr);

    switch (yr),
     %DEBUG:     case 2009,		mos = 4:6;
     %DEBUG:     case 2013,		mos = 2:4;
     case thisyear,	mos = 1:(thismonth-1);
     otherwise,		mos = 1:12;
    end;

    for mo = mos(:)'

      for ddix=1:length(dds)
        dd = dds(ddix);
        dd = dd{:};
        ds = dd{1};
        vars = dd{2};
        cvts = dd{3};
        flds = dd{4};

        fname = fullfile(erapath,sprintf('ERA_Interim_%04d%02d_%s.grib',yr,mo,ds));
        if ( ~exist(fname,'file') )
          warning('GRiB file not found! "%s"',fname);
          continue;
        end;

        % nc = mDataset(fname);
        %% Many screens full of Java warnings we can't turn off...
        [spewage,nc] = evalc('mDataset(fname)');
        if ( isempty(nc) )
          warning('Invalid GRiB file "%s"',fname);
          clear nc;
          continue;
        end;
        if ( isempty(lats) )
          %DEBUG:
          disp('Extracting grid coordinates');
          lons = cast( nc{'lon'}(:,:,:), 'double' );
          lats = cast( nc{'lat'}(:,:,:), 'double' );
          lons(lons > 180) = lons(lons > 180) - 360;
          dlat = min(diff(unique(lats(:))));
          dlon = min(diff(unique(lons(:))));
          badix = [];
          for stix=1:length(stnms)
            stnm = stnms{stix};
            result(stix).station_name = stnm;
            if ( ~isfield(result(stix),'lon') || ~isfield(result(stix),'lat') ...
                 || isempty(result(stix).lon) || isempty(result(stix).lat) )
              try
                [result(stix).lon,result(stix).lat,result(stix).depth] = get_station_coords(stnm);
              catch
                warning('Station "%s": found no coordinates! Skipping...',stnm);
                badix(end+1) = stix;
                continue;
              end;
            end;
            [laterr,result(stix).erai_latix] = min(abs(result(stix).lat-lats));
            [lonerr,result(stix).erai_lonix] = min(abs(result(stix).lon-lons));
            if ( laterr > (1.1*dlat) || lonerr > (1.1*dlon) )
              warning('Station "%s" outside GRiB file region [%g,%g,%g,%g]! Skipping...',...
                      stnm,min(lats),max(lats),min(lons),max(lons));
              badix(end+1) = stix;
            end; %if laterr
          end; %for stix
          stnms(badix) = [];
          result(badix) = [];
          if ( isempty(result) )
            close(nc); clear nc;
            error('None of the named stations fall within GRiB file region!');
          end;
        end; %if isempty(lats)

        hrs = cast( nc{'time'}(:,:,:), 'double');
        dts = datenum(yr,mo,1) + (hrs./24);
        ndts = numel(dts);
        for varix = 1:length(vars)
          var = vars{varix};
          fld = ['raw_erai_' flds{varix}];
          cvt = cvts{varix};
          dat = cast( nc{var}(:,:,:), 'double' );
          dat = get_erai_station_apply_cvt(cvt,dts,dat);

          for stix=1:length(stnms)
            result(stix).(fld).date(end+1:end+ndts,1) = dts(:);

            %% For simple nearest-neighbor assignment - DOES NOT HANDLE NaNs!
            % result(stix).(fld).data(end+1:end+ndts,1) ...
            %     = dat(:,result(stix).erai_latix,result(stix).erai_lonix);

            if ( is_erai_ocean_var(flds{varix}) )
              % Do what we have to do to get a value
              result(stix).(fld).data(end+1:end+ndts,1) = ...
                  interp_field(lats,lons,dat,result(stix).lat,result(stix).lon,@nanmean,'warn',true);
            else
              % Use bilinear (2-D) spatial interpolation
              %% Stupid INTERP2/INTERP3 is just too hard to use
              result(stix).(fld).data(end+1:end+ndts,1) = ...
                  interp_field(lats,lons,dat,result(stix).lat,result(stix).lon,'linear','warn');
            end; %if is_ocean_var else
          end; %for stix

          dat = []; clear dat;
        end; %for varix

        close(nc); clear nc;
        hrs = []; clear hrs;
        dts = []; clear dts;

      end; %for ddix

    end; %for mo

  end; %for yr

  % Special handling for raw surface wave data
  for stix=1:length(stnms)
    % The following code is identical to that in GET_WW3_STATION (v.)...
    % Calculate "u" and "v" peak wavenumber vector components by weighting
    % PEAKWAVEDIR with putative (deep-water) wave celerity: thus. the lower
    % frequency waves dominate when averaged in with shorter wavelengths.
    % (As for any vector, interpolating PEAKWAVEDIR means special handling!)
    g = sw_g(result(stix).lat,0);
    c = g .* (result(stix).raw_erai_peakwaveper.data) ./ (2*pi);
    % % *OR* just do a simple vector averaging with all weights == 1??
    % c = 1;
    result(stix).raw_erai_peakwave_speed.date = result(stix).raw_erai_peakwavedir.date;
    result(stix).raw_erai_peakwave_speed.data = c;

    result(stix).raw_erai_peakwave_u.date = result(stix).raw_erai_peakwavedir.date;
    result(stix).raw_erai_peakwave_v.date = result(stix).raw_erai_peakwavedir.date;
    [result(stix).raw_erai_peakwave_u.data,...
     result(stix).raw_erai_peakwave_v.data] = ...
        spddir_to_uv(result(stix).raw_erai_peakwave_speed.data,result(stix).raw_erai_peakwavedir.data);
  end;

  % Calculate some additional flux fields from those given
  for stix=1:length(stnms)
    result(stix).raw_erai_usrf.date = result(stix).raw_erai_dsrf.date;
    result(stix).raw_erai_usrf.data = ...
        result(stix).raw_erai_dsrf.data - result(stix).raw_erai_srf.data;
    result(stix).raw_erai_ulrf.date = result(stix).raw_erai_dlrf.date;
    result(stix).raw_erai_ulrf.data = ...
        result(stix).raw_erai_dlrf.data - result(stix).raw_erai_lrf.data;
  end;

  % Now spline-interpolate raw 3-hourly data into hourly fields
  rawflds = grepstruct(result(1),'raw_erai_');
  for fldix = 1:length(rawflds)
    rawfld = rawflds{fldix};
    fld = rawfld(5:end);

    rawdts = result(1).(rawfld).date;
    dts = [rawdts(1):(1/24):rawdts(end)]';
    ndts = length(dts);
    for stix=1:length(stnms)
      rawdat = result(stix).(rawfld).data;
      result(stix).(fld).date(1:ndts,1) = dts(:);
      if ( any(isfinite(rawdat)) )
        result(stix).(fld).data(1:ndts,1) = interp1(rawdts,rawdat,dts(:),'spline');
      else
        warning('Field %s.%s is all NaN!',result(stix).station_name,rawfld);
        result(stix).(fld).data(1:ndts,1) = NaN;
      end;
      result(stix) = filter_gaps(result(stix),rawfld,fld,[],(6/24),[],nan);
    end;
  end;

  % Very simplistic QA:
  % Certain spline-interpolated fields should never be negative
  for stix=1:length(stnms)
    result(stix).erai_cloud_cover.data(result(stix).erai_cloud_cover.data<0) = 0;
    result(stix).erai_cloud_cover.data(result(stix).erai_cloud_cover.data>1) = 1;
    result(stix).erai_conv_precip.data(result(stix).erai_conv_precip.data<0) = 0;
    result(stix).erai_uvW.data(result(stix).erai_uvW.data<0) = 0;
    % (Evaporation should never be positive!)
    result(stix).erai_evap.data(result(stix).erai_evap.data>0) = 0;
    result(stix).erai_parW.data(result(stix).erai_parW.data<0) = 0;
    result(stix).erai_runoff.data(result(stix).erai_runoff.data<0) = 0;
    result(stix).erai_sunshine_duration.data(result(stix).erai_sunshine_duration.data<0) = 0;
    result(stix).erai_clear_sky_srf.data(result(stix).erai_clear_sky_srf.data<0) = 0;
    result(stix).erai_dsrf.data(result(stix).erai_dsrf.data<0) = 0;
    result(stix).erai_usrf.data(result(stix).erai_usrf.data<0) = 0;
    result(stix).erai_srf.data(result(stix).erai_srf.data<0) = 0;
    result(stix).erai_dlrf.data(result(stix).erai_dlrf.data<0) = 0;
    result(stix).erai_ulrf.data(result(stix).erai_ulrf.data<0) = 0;
    result(stix).erai_precip.data(result(stix).erai_precip.data<0) = 0;
  end;

  % Additional calculations for some (raw and) interpolated fields
  for stix=1:length(stnms)
    % Dew P -> Rel Humid -> Spec Humid
    result(stix).erai_relhumid.date = result(stix).erai_dew_t.date;
    result(stix).erai_relhumid.data = ...
        dewp_to_relhumid(result(stix).erai_air_t.data,result(stix).erai_dew_t.data);
    result(stix).erai_spechumid.date = result(stix).erai_relhumid.date;
    result(stix).erai_spechumid.data = ...
        relhumid_to_spechumid(result(stix).erai_air_t.data,result(stix).erai_relhumid.data);

    % Wind U,V [m/s] -> Wind Speed [kts] -> Wind Dir
    result(stix).erai_wind_speed.date = result(stix).erai_wind_u.date;
    result(stix).erai_wind_speed.data = ...
        mps2kts(uv_to_spd(result(stix).erai_wind_u.data,result(stix).erai_wind_v.data));
    result(stix).erai_wind_dir.date = result(stix).erai_wind_u.date;
    result(stix).erai_wind_dir.data = ...
        uv_to_dir(result(stix).erai_wind_u.data,result(stix).erai_wind_v.data);

    % We want rain rate in mm/hr, not m/hr!
    result(stix).raw_erai_conv_precip.data = result(stix).raw_erai_conv_precip.data*1e3;
    result(stix).raw_erai_precip.data = result(stix).raw_erai_precip.data*1e3;
    result(stix).erai_conv_precip.data = result(stix).erai_conv_precip.data*1e3;
    result(stix).erai_precip.data = result(stix).erai_precip.data*1e3;

    % Photosynthetically Active Radiation [W/m^2] -> micro-mol quanta/m^2/s
    result(stix).erai_par.date = result(stix).erai_parW.date;
    result(stix).erai_par.data = result(stix).erai_parW.data .* UMOL_PER_WATT;

    % Vector-interpolated peak wave direction
    % This overwrites the interpolated garbage calculated above
    result(stix).erai_peakwavedir.date = result(stix).erai_peakwave_u.date;
    result(stix).erai_peakwavedir.data = ...
        uv_to_dir(result(stix).erai_peakwave_u.data,result(stix).erai_peakwave_v.data);
  end;  

  % Finally, save all the station data we extracted as MAT files
  for stix=1:length(stnms)
    stnm = stnms{stix};
    matfname = fullfile(datapath,[stnm '_erai.mat']);
    station = result(stix);
    disp(['Saving to ' matfname]);
    save(matfname,'station');
    station = []; clear station;
  end;

  % And return data for the last station we extracted in struct STN
  flds = fieldnames(result(end));
  for ix = 1:length(flds)
    fld = flds{ix};
    stn.(fld) = result(end).(fld);
  end;

  result = []; clear result;

  %%%% DEBUG??? Post-correct a spline-fitting error
  stn.erai_cloud_cover.data(stn.erai_cloud_cover.data>1) = 1;
  %%%% DEBUG???

  toc,
  set_more;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   PRIVATE FUNCTIONS   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dat = get_erai_station_apply_cvt(cvt,dts,dat)
  [op,valstr] = strtok(cvt);
  val = str2num(valstr);
  switch ( op ),
   case 'add',
    dat = dat + val;
   case 'mult',
    dat = dat .* val;
   case 'deaccum',
    % This is an oddity of some ECMWF forecast products: the value in the
    % GRiB file is actually an accumulated quantity, reset every 12 hours
    resetix = find(get_hour(dts) == 3 | get_hour(dts) == 15);
    resets = dat(resetix,:,:);
    dat(2:end,:,:) = diff(dat)./(3*3600);
    dat(resetix,:,:) = resets./(3*3600);
   case 'deaccum_sum',
    % See comment for 'deaccum' above: This option is an ODDITY OF AN ODDITY
    % for precipitation data, in that the value in the GRiB file is actually
    % an accumulated TOTAL quantity (i.e, *not* per sec), reset every 12 hrs
    resetix = find(get_hour(dts) == 3 | get_hour(dts) == 15);
    resets = dat(resetix,:,:);
    dat(2:end,:,:) = diff(dat)./(3);
    dat(resetix,:,:) = resets./(3);
   case 'none',
   otherwise,
    error('Unrecognized conversion string "%s"',cvt);
  end;
return;


%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%

function yorn = is_erai_ocean_var(fld)
  yorn = ismember(fld,{'sea_t','peakwavedir','peakwaveper','sigwavehgt'});
return;
