function stn = get_fkeys_hycom(stn_or_stnm,mindt,maxdt,vars,flds,interpMethod,fkeyspath)
%function stn = get_fkeys_hycom(stn_or_stnm,mindt,maxdt,vars,flds,interpMethod,fkeyspath)
%
% Add fields for ocean surface U and V current components, sea temperature T,
% salinity S, and time series of 17x17 grid-point field for U, V, T, and S,
% from non-assimilative RSMAS FKEYS HYCOM (1/100 degree) hydrodynamic model
% of Kourafalou and Kang, to struct STN.  If station name string STNM (five
% characters, e.g., 'mlrf1') is given instead of a struct, or if struct has
% no valid .lon and .lat fields (e.g., -80.38, 25.01), GET_STATION_COORDS is
% called with the station name, to try to retrieve the station's location.
%
% Optional VARS is a string or cellstr specifying which variables to extract.
% If VARS is specified, optional FLDS may also be specified for the corresp.
% list of struct variable names to be added to STN, e.g., 'seatemp' for 't'.
% DEFAULT VARS and corresponding FLDS cell arrays are as follows:
%  vars = { 'u', 'v', 'temperature', 'salinity', 'u' ,      'v' ,      'temperature' ,  'salinity',       };
%  flds = { 'u', 'v', 'seatemp',     'salinity', 'u_field', 'v_field', 'seatemp_field', 'salinity_field', };
%
% INTERP_FIELD is called with INTERPMETHOD (DEFAULT 'linear') to do subset.
%
% NOTE: This differs from the similar function GET_FLKEYS_HYCOM, in that we
% assume all FKEYS HYCOM data are accessible as local netCDF files in a dir
% FKEYSPATH off the THESIS path (DEFAULT: [DATAPATH '/hycom/FKEYS']).
%
% CALLS: MDATASET (netCDF-Java); INTERP_FIELD, GET_STATION_COORDS (Ecoforecasts).
%
% Last Saved Time-stamp: <Wed 2011-09-07 14:38:23  lew.gramer>

  set_more off

  datapath = get_thesis_path('../data');

  if ( ~exist('vars','var') || isempty(vars) )
    vars = { 'u', 'v', 'temperature', 'salinity', 'u' ,      'v' ,      'temperature' ,  'salinity',       };
    flds = { 'u', 'v', 'seatemp',     'salinity', 'u_field', 'v_field', 'seatemp_field', 'salinity_field', };
  end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = lower(vars);
  end;

  interpMethodDefaulted = false;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
    % If user previously subset this site, use interpMethod specified in that
    % initial run as our default for this site as well (below)!
    interpMethodDefaulted = true;
  end;
  if ( interpMethod(1) == '*' )
    interpMethod = interpMethod(2:end);
  end;

  if ( ~exist('fkeyspath','var') || isempty(fkeyspath) )
    fkeyspath = fullfile(datapath,'hycom','FKEYS');
  end;

  % What are the first and last days with available model data?
  %%%% ??? NOTE: These are currently hard-wired: will need modifying later!
  zerodt = datenum(2004,1,1);
  lastdt = datenum(2009,1,1) - (0.5/24);
  % Model outputs are once-per-six hours time series
  alldts = zerodt:(6/24):lastdt;



  if ( ischar(stn_or_stnm) )
    stn.station_name = stn_or_stnm;
  elseif ( isstruct(stn_or_stnm) )
    stn = stn_or_stnm;
    if ( ~isfield(stn,'station_name') )
      error('Station STRUCT with no station_name field!');
    end;
  else
    error('First arg must either be station STRUCT or station name string!');
  end;

  if ( ~isfield(stn,'lon') )
    if ( isfield(stn,'station_name') )
      [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
    else
      error('Station specified without coordinates or name!');
    end;
  end;


  matfname = fullfile(datapath,[lower(stn.station_name) '_fkeys_hycom.mat']);

  % If we already did this before, just load MAT file and subset
  if ( exist(matfname,'file') )
    disp(['Reloading from MAT file ' matfname]);
    load(matfname,'station');

  else

    disp('Loading original data from local netCDF files...');

    % Grid-point radii for 'field' 3-D time-series (e.g., seatemp_field)
    % Choose one axis slightly larger to guard against transposition errors
    yrad = 8;
    xrad = 9;

    if ( ~exist('mindt','var') || isempty(mindt) )
      mindt = zerodt;
    end;
    if ( ~exist('maxdt','var') || isempty(maxdt) )
      maxdt = lastdt;
    end;

    ixes = find(mindt <= alldts & alldts <= maxdt);
    alldts = alldts(ixes);

    % First time ever calling this function, we want to get ALL the data
    [begyr,ig,ig] = datevec(zerodt);
    [endyr,ig,ig] = datevec(lastdt);


    if ( isfield(stn,'station_name') )
      station.station_name = stn.station_name;
    end;
    station.lon = stn.lon; station.lat = stn.lat;
    station.fkeys_hycom_interp_method = lower(interpMethod);
    for vix = 1:length(vars)
      fld = ['fkeys_hycom_' flds{vix}];
      % Preallocate for fast execution
      station.(fld) = [];
      if ( ~isempty(strfind(fld,'_field')) )
        station.(fld).date = repmat(nan,[length(alldts),1]);
        station.(fld).field = repmat(nan,[length(alldts),(2*yrad)+1,(2*xrad)+1]);
      else
        station.(fld).date = repmat(nan,[length(alldts),1]);
        station.(fld).data = repmat(nan,[length(alldts),1]);
      end;
    end;

    lats = [];
    lons = [];

    for yr=begyr:endyr
      %DEBUG:
      disp(yr); tic,

      if ( yr>=2008 )
        expt = '902';
      else
        expt = '303';
      end;

      jds = 1:365;
      % Will we never need to deal with years 1900 or 2100??
      if ( mod(yr,4) == 0 )
        jds = 1:366;
      end;

      for jd=jds(:)'
        for hr=0:6:18

          dt = datenum(yr,1,1,hr,0,0) + jd - 1;
          [ig,dtix] = min(abs(alldts - dt));

          % First make sure we have ALL variables for this date/time
          fpatt = fullfile(fkeyspath,sprintf('%s_archv.%04d_%03d_%02d_3z*.nc',expt,yr,jd,hr));
          flist = dir(fpatt);
          if ( isempty(flist) )
            warning('Missing all variables for date/time "%s"!',fpatt);
            continue;
          elseif ( length(flist) < length(unique(vars)) )
            warning('Missing some variables for date/time "%s"!',fpatt);
            continue;
          end;

          for vix = 1:length(vars)
            var = vars{vix};
            varl = var(1);
            fld = ['fkeys_hycom_' flds{vix}];
            %DEBUG:            disp(fld);

            station.(fld).date(dtix,1) = dt;

            %303_archv.2004_001_06_3zs.nc
            fbasename = sprintf('%s_archv.%04d_%03d_%02d_3z%s.nc',expt,yr,jd,hr,varl);
            fname = fullfile(fkeyspath,fbasename);
            if ( ~exist(fname,'file') )
              warning('Missing file "%s"!',fname);
              continue;
            end;
            %DEBUG:            disp(fname);
            nc = mDataset(fname);
            if ( isempty(nc) )
              warning('Unable to open "%s"!',fname);
              continue;
            end;

            % NJTBX-2.0 Toolbox does not handle unclosed Datasets very well
            try
              if ( isempty(lats) || isempty(lons) )
                lats = cast(nc{'Latitude'}(:),'double');
                lons = cast(nc{'Longitude'}(:),'double');

                % [yix,xix] = query_fkeys_hycom_indices(station.lon,station.lat);
                [yerr,yix] = min(abs(lats - station.lat));
                [xerr,xix] = min(abs(lons - station.lon));
                if ( yerr > max(diff(unique(lats)))+0.01 || xerr > max(diff(unique(lons)))+0.01 )
                  close(nc); clear nc
                  error('Station "%s" at %g,%g is outside FKEYS HYCOM domain!',...
                        station.station_name,station.lat,station.lon);
                end;

                lat = lats(yix-1:yix+1);
                lon = lons(xix-2:xix+2);
              end;

              if ( ~isempty(strfind(fld,'_field')) )
                if ( ~isfield(station.(fld),'lat') || ~isfield(station.(fld),'lon') )
                  station.(fld).lat = lats(yix-yrad:yix+yrad);
                  station.(fld).lon = lons(xix-xrad:xix+xrad);
                end;

                station.(fld).field(dtix,:,:) = ...
                    squeeze(cast(nc{var}(1,1,yix-yrad:yix+yrad,xix-xrad:xix+xrad),'double'));

              else
                datsz = getShape(nc{var});
                % 2-D data elements
                if ( length(datsz) == 3 )
                  dat = squeeze(cast(nc{var}(1,yix-1:yix+1,xix-2:xix+2),'double'));
                % 3-D data elements
                else
                  dat = squeeze(cast(nc{var}(1,1,yix-1:yix+1,xix-2:xix+2),'double'));
                end;
                % % More X-Y vs. row-column confusion: stupid MATLAB...
                % station.(fld).data(dtix,1) = ...
                %     interp2(lon,lat,dat,station.lon,station.lat,interpMethod);
                station.(fld).data(dtix,1) = ...
                    interp_field(lat,lon,dat,station.lat,station.lon,interpMethod,'warn');
              end;

            catch
              close(nc); clear nc;
              rethrow(lasterror);
            end;

            close(nc); clear nc;

          end; %for vix = 1:length(vars)
        end; %for hr=0:6:18
      end; %for jd=1:366

      %DEBUG:
      toc,
    end; %for yr=begyr:endyr

    % Remove dates for which we were missing any data files (u, v, T, or S)
    for vix = 1:length(vars)
      fld = ['fkeys_hycom_' flds{vix}];
      badix = find(~isfinite(station.(fld).date));
      station.(fld).date(badix) = [];
      if ( ~isempty(strfind(fld,'_field')) )
        station.(fld).field(badix,:,:) = [];
      else
        station.(fld).data(badix) = [];
      end;
    end;

    % Calculate speed and direction from model U and V currents
    if ( isfield(station,'fkeys_hycom_u') && isfield(station,'fkeys_hycom_v') )
      station.fkeys_hycom_speed.date = station.fkeys_hycom_u.date;
      station.fkeys_hycom_speed.data = uv_to_spd(station.fkeys_hycom_u.data,station.fkeys_hycom_v.data);
      station.fkeys_hycom_dir.date = station.fkeys_hycom_u.date;
      station.fkeys_hycom_dir.data = uv_to_dir_curr(station.fkeys_hycom_u.data,station.fkeys_hycom_v.data);
    end;

    disp(['Saving to MAT file ' matfname]);
    save(matfname,'station');

  end;

  % Limit our result fields to only those dates we requested
  for vix = 1:length(vars)
    fld = ['fkeys_hycom_' flds{vix}];

    % Gross quality control
    if ( ~isfield(station,fld) )
      warning('No field "%s" found after load!',fld);
    elseif ( ~isempty(strfind(fld,'_field')) )
      % Make sure lon and lat are row vectors as for other models
      station.(fld).lon = station.(fld).lon(:);
      station.(fld).lat = station.(fld).lat(:);

      % 40 just happens to be an upper bound for all of u,v,T,S
      ixes = find(-10 >= station.(fld).field | station.(fld).field >= 40);
      station.(fld).field(ixes) = nan;

      ixes = find(alldts(1) <= station.(fld).date & station.(fld).date <= alldts(end));
      station.(fld).date = station.(fld).date(ixes);
      station.(fld).field = station.(fld).field(ixes,:,:);
    else
      % 40 just happens to be an upper bound for all of u,v,T,S
      ixes = find(-10 <= station.(fld).data & station.(fld).data <= 40);
      station.(fld).date = station.(fld).date(ixes);
      station.(fld).data = station.(fld).data(ixes);

      ixes = find(alldts(1) <= station.(fld).date & station.(fld).date <= alldts(end));
      station.(fld).date = station.(fld).date(ixes);
      station.(fld).data = station.(fld).data(ixes);
    end;
  end;

  % If we are requesting a different interpolation method than that used to
  % save our MAT file, then re-call INTERP_FIELD on our four "_field" fields
  if ( ~isfield(station,'fkeys_hycom_interp_method') || ...
       ~strcmpi(station.fkeys_hycom_interp_method,interpMethod) )

    % If someone previously subset this site, and this caller did not specify
    % interpMethod, use interpMethod specified in that initial run as default:
    % This complexity is useful for those sites where 'linear' would always
    % get us into trouble, e.g., those very near land-masked pixels.
    if ( interpMethodDefaulted && isfield(station,'fkeys_hycom_interp_method') )
      disp([mfilename ': Using original site INTERPMETHOD ''' station.fkeys_hycom_interp_method '''']);
    else
      warning('Reinterpolating time series fields using %s',upper(interpMethod));
      for vix = 1:length(vars)
        fld = ['fkeys_hycom_' flds{vix}];
        if ( isempty(strfind(fld,'_field')) )
          fnm = [fld '_field'];
          if ( isfield(station,fld) && isfield(station,fnm) )
            f = station.(fnm);
            station.(fld).data = ...
                interp_field(f.lat,f.lon,f.field,station.lat,station.lon,interpMethod,'warn');
          end;
        end;
      end;
    end;

  end;

  allflds = grepstruct(station,'fkeys_hycom');
  for fldix = 1:length(allflds)
    fld = allflds{fldix};
    stn.(fld) = station.(fld);
  end;
  station = []; clear station;

  set_more;

return;
