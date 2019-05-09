function [stn,clim] = station_load_oaflux(stn_or_stnm)
%function [stn,clim] = station_load_oaflux(stn_or_stnm)
%
% Load flux terms (net, short-/longwave radiative, sensible/latent, wind
% speed in m/s, evaporation, specific humidity, air and sea temp.) from daily
% OAflux v.3 heat budget of Yu and Weller (2007), for station named by string
% STN.station_name or STNM. Return station struct STN with new '.oaflux_*'
% time series fields added, and CLIM struct with just '.oaflux_*' fields.
%
% Last Saved Time-stamp: <Fri 2011-11-04 17:11:40  Lew.Gramer>

  datapath = get_thesis_path('../data');
  oafluxpath = fullfile(datapath,'oaflux');

  if ( ischar(stn_or_stnm) )
    stn.station_name = stn_or_stnm;
  elseif ( isfield(stn_or_stnm,'station_name') )
    stn = stn_or_stnm;
  else
    error('First arg must be station name or struct with field .station_name!');
  end;

  if ( ~isfield(stn,'lat') )
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
  end;

  %DEBUG:  tic,

  clim = [];

  matfname = fullfile(datapath,sprintf('%s_oaflux.mat',lower(stn.station_name)));

  if ( exist(matfname,'file') )
    disp(sprintf('Loading station OAFlux climatology from "%s"',matfname));
    load(matfname,'clim');

  else
    disp('Subsetting station OAFlux climatology from raw netCDF...');

    %  NC file rootname  NC variable names    Factor  STN field names
    fnames = { ...
    %daily mean net surface heat flux, positive downward: W/m^2
        'qnet',		{'qnet'},		1,   {'net_heat_flux'}; ...
    %daily mean net surface fullsky shortwave radiation flux, positive downward: W/m^2
        'sw_isccp',	{'nswrs'},		1,   {'srf'}; ...
    %daily mean net surface fullsky longwave radiation flux, positive upward: W/m^2 -> pos. down
        'lw_isccp',	{'nlwrs'},		-1,  {'lrf'}; ...
    %daily mean net surface fullsky longwave radiation flux, positive upward: W/m^2 -> pos. down
        'sh_oaflux',	{'shtfl','err'},	-1,  {'sensible_heat_flux','sensible_flux_err'}; ...
    %daily mean surface latent heat flux, positive upward: W/m^2 -> pos. down
        'lh_oaflux',	{'lhtfl','err'},	-1,  {'latent_heat_flux','latent_flux_err'}; ...
    %daily mean evaporation rate: cm/yr -> mm/hr
        'evapr_oaflux',	{'evapr','err'},	.00114155, {'evap','evap_err'}; ...
    %daily mean neutral wind speed at 10m: m/s -> kts
        'ws_oaflux',	{'wnd10','err'},	1.9438, {'wind_speed','wind_speed_err'}; ...
    %daily mean air temperature at 2m: degree C
        'ta_oaflux',	{'tmp2m','err'},	1,   {'air_t','air_t_err'}; ...
    %daily mean sea surface temperature: degree C
        'ts_oaflux',	{'tmpsf','err'},	1,   {'seatemp','seatemp_err'}; ...
    %daily mean specific humidity at 2m: g/Kg -> Kg/Kg
        'qa_oaflux',	{'hum2m','err'},	1e3, {'spechumid','spechumid_err'}; ...
             };

    yrs = [1987:get_year(now)];
    %DEBUG:    yrs = [2000];
    %DEBUG:    yrs = [2010];

    for ix = 1:length(fnames)
      fname = fnames{ix,1};
      flds  = fnames{ix,4};
      for fldix = 1:length(flds)
        fld = ['oaflux_' flds{fldix}];
        clim.(fld) = struct('date',[],'data',[]);
      end;
    end;

    lats = [];
    lons = [];
    nc = [];
    % NETCDF JAVA datasets left open forever make MATLAB unhappy: use TRY-CATCH 
    try
      for yr = yrs(:)'

        %DEBUG:
        disp(yr);

        dts = datenum(yr,1,1):1:datenum(yr,12,31);

        for ix = 1:length(fnames)
          fname = fnames{ix,1};
          fac   = fnames{ix,3};

          vars  = fnames{ix,2};
          flds  = fnames{ix,4};

          ncfname = fullfile(oafluxpath,sprintf('%s_%04d.nc',fname,yr));
          %DEBUG:          disp(ncfname);
          if ( ~exist(ncfname,'file') )
            warning('Skipping missing "%s"',ncfname);
            continue;
          end;
          nc = mDataset(ncfname);

          for fldix = 1:length(flds)
            var = vars{fldix};
            fld = ['oaflux_' flds{fldix}];
            %DEBUG:  nj_info(nc),
            if ( isempty(lats) )
              lats = nc{'lat'}(:,:);
              lons = nc{'lon'}(:,:);
              % Idiotic negative West longitudes
              lons(lons>180) = lons(lons>180) - 360;
              [ig,lonix]=min(abs(stn.lon-lons));
              [ig,latix]=min(abs(stn.lat-lats));
              clim.oaflux_lonix = lonix;
              clim.oaflux_latix = latix;
            end; %if isempty(lats)
            dat = cast(nc{var}(:,clim.oaflux_latix,clim.oaflux_lonix),'double');

            % First time through each file, figure out how many days it contains
            if ( fldix==1 )
              ndays = length(dat);
              dts = dts(1:ndays);
            end;
            clim.(fld).date(end+1:end+ndays,1) = dts(:);
            clim.(fld).data(end+1:end+ndays,1) = fac .* dat(:);
          end; %for flds

          close(nc); nc = [];

        end; %for fnames
      end;  %for yrs

    catch
      if ( ~isempty(nc) )
        close(nc); nc = [];
        warning('Failure on file "%s"',ncfname);
        rethrow(lasterror);
      end;
    end; %try

    disp(sprintf('Saving station OAFlux data to "%s"',matfname));
    save(matfname,'clim');

  end; %if exist(matfname) else

  flds = fieldnames(clim);
  for fldix = 1:length(flds)
    fld = flds{fldix};
    stnfld = ['daily_' fld];
    stn.(stnfld) = clim.(fld);
  end;

  %DEBUG:  toc,

return;
