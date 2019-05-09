function stn = get_flkeys_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
%function stn = get_flkeys_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
%
% Add fields for ocean surface U and V current components, sea temperature,
% salinity, mixed-layer depth, and 9x9 grid-point sea temperature field, from
% the assimilative NRL Gulf of Mexico HYCOM (1/25 degree) hydrodynamic model
% of P. Hogan et al., to struct STN.  If a station name string STNM (five
% characters, e.g., 'mlrf1') is given instead of a struct, or if struct has
% no valid .lon and .lat fields (e.g., -80.38, 25.01), GET_STATION_COORDS is
% called with the station name, to try to retrieve the station's location.
%
% Optional VARS is a string or cellstr specifying which variables to extract.
% No online catalog seems to be available, but see the URL templates below
% for a hard-won list of all the variables available from this dataset.
% If VARS is specified, optional FLDS may also be specified for the corresp.
% list of struct variable names to be added to STN, e.g., 'seatemp'.
% DEFAULT VARS and corresponding FLDS cell arrays are as follows:
%   vars = { 'u', 'v', 'mld', 'qtot',          'temperature', 'salinity', 'u' ,      'v' ,      'temperature'    };
%   flds = { 'u', 'v', 'mld', 'net_heat_flux', 'seatemp',     'salinity', 'u_field', 'v_field', 'seatemp_field' };
%
% DEFAULT for BASEURL:
%  'http://tds.hycom.org/thredds/dodsC/flkeys'
% If you have another ocean model with a THREDDS interface, specify it here.
%
% TO DO - NEW BASE URL:
%  'http://tds.hycom.org/thredds/dodsC/flkeys'
%
% SAMPLE URL:
%  http://tds.hycom.org/thredds/dodsC/flkeys.ascii?Latitude[244],Longitude[298]
%  (ORIG) http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1/2003
%
% CALLS: MDATASET (netCDF-Java), GET_STATION_COORDS, QUERY_FLKEYS_HYCOM_INDICES,
%        QUERY_FLKEYS_HYCOM_STATION_FIELD_COORDS.
%
% Last Saved Time-stamp: <Fri 2013-04-26 13:48:38 Eastern Daylight Time gramer>

  set_more off;

  datapath = get_thesis_path('../data');

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

    if ( ~exist('vars','var') || isempty(vars) )
      vars = { 'u', 'v', 'temperature', 'salinity', 'u' ,      'v' ,      'temperature',  };
      flds = { 'u', 'v', 'seatemp',     'salinity', 'u_field', 'v_field', 'seatemp_field' };
    end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = lower(vars);
  end;


  % What are the first and last days with available model data?
  %%%% ??? NOTE: These are currently hard-wired: will need modifying later!
  zerodt = datenum(2008,1,1);
  lastdt = datenum(2009,1,1) - (0.5/24);
  % Model outputs are once-per-six hours time series
  alldts = zerodt:(6/24):lastdt;


  matfname = fullfile(datapath,[lower(stn.station_name) '_flkeys_hycom.mat']);

  % If we already did this before, just load MAT file and subset
  if ( exist(matfname,'file') )
    disp(['Reloading from MAT file ' matfname]);
    x = load(matfname,'stn');
    allflds = grepstruct(x.stn,'flkeys_hycom');
    for fldix = 1:length(allflds)
      fld = allflds{fldix};
      stn.(fld) = x.stn.(fld);
    end;
    x = []; clear x;

  else

    disp('Loading original data from online netCDF...');

    % Grid-point radii for time-series fields (e.g., seatemp_field)
    yrad = 4;
    xrad = 4;

    if ( ~exist('mindt','var') || isempty(mindt) )
      mindt = zerodt;
    end;
    if ( ~exist('maxdt','var') || isempty(maxdt) )
      maxdt = lastdt;
    end;

    if ( ~exist('baseurl','var') || isempty(baseurl) )
      % %% "Legacy" THREDDS interface
      % %% baseurl = 'http://tds.hycom.org/opendap/nph-dods/datasets/hycom/GOMl0.04/expt_20.1';

      % % New THREDDS interface
      % baseurl = 'http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1';

      % New THREDDS interface
      baseurl = 'http://tds.hycom.org/thredds/dodsC/flkeys';
    end;

    ixes = find(mindt <= alldts & alldts <= maxdt);
    alldts = alldts(ixes);

    % First time ever calling this function, we want to get ALL the data
    [begyr,ig,ig] = datevec(zerodt);
    [endyr,ig,ig] = datevec(lastdt);


    [yix,xix] = query_flkeys_hycom_indices(stn.lon,stn.lat);

    for yr=begyr:endyr

      dts = [];

      % url = sprintf('%s/%04d',baseurl,yr);
      url = baseurl;
      %DEBUG:
      disp(url);
      nc = mDataset(url);
      if ( isempty(nc) )
        warning('Invalid URL "%s"!',url);
        continue;
      end;

      for vix = 1:length(vars)
        var = vars{vix};
        fld = ['flkeys_hycom_' flds{vix}];
        %DEBUG:
        disp(fld);

        if ( isempty(dts) )
          dts = cast(getTimes(nc{var}),'double');
          ndts = numel(dts);
        end;
        if ( yr == begyr )
          stn.(fld) = struct('date',[],'data',[]);
        end;

        stn.(fld).date(end+1:end+ndts,1) = dts';
        % if ( strcmp(flds{vix},'seatemp_field') )
        if ( ~isempty(strfind(flds{vix},'_field')) )
          stn.(fld).data(end+1:end+ndts,1:(2*yrad)+1,1:(2*xrad)+1) = ...
              squeeze(cast(nc{var}(1:end,1,yix-yrad:yix+yrad,xix-xrad:xix+xrad),'double'));
        else
          datsz = getShape(nc{var});
          % 2-D data elements
          if ( length(datsz) == 3 )
            stn.(fld).data(end+1:end+ndts,1) = ...
                squeeze(cast(nc{var}(1:end,yix,xix),'double'));
            % 3-D data elements
          else
            stn.(fld).data(end+1:end+ndts,1) = ...
                squeeze(cast(nc{var}(1:end,1,yix,xix),'double'));
          end;
        end;
      end;

      close(nc); clear nc;

    end;

    % Calculate speed and direction from model U and V currents
    if ( isfield(stn,'flkeys_hycom_u') && isfield(stn,'flkeys_hycom_v') )
      stn.flkeys_hycom_speed.date = stn.flkeys_hycom_u.date;
      stn.flkeys_hycom_speed.data = uv_to_spd(stn.flkeys_hycom_u.data,stn.flkeys_hycom_v.data);
      stn.flkeys_hycom_dir.date = stn.flkeys_hycom_u.date;
      stn.flkeys_hycom_dir.data = uv_to_dir_curr(stn.flkeys_hycom_u.data,stn.flkeys_hycom_v.data);
    end;

    % Modify STN.flkeys_hycom_*_field structs: change field .data to
    % .field, and add Nx1 fields .lon and .lat, before saving to MAT.
    stn = query_flkeys_hycom_station_field_coords(stn,'flkeys_hycom_u_field');
    stn = query_flkeys_hycom_station_field_coords(stn,'flkeys_hycom_v_field');
    stn = query_flkeys_hycom_station_field_coords(stn,'flkeys_hycom_seatemp_field');

    disp(['Saving to MAT file ' matfname]);
    save(matfname,'stn');

  end;

  % Limit our result fields to only those dates we requested
  for vix = 1:length(vars)
    fld = ['flkeys_hycom_' flds{vix}];

    % Gross quality control
    if ( ~isfield(stn,fld) )
      warning('No field "%s" found after load!',fld);
    % elseif ( strcmp(flds{vix},'seatemp_field') )
    elseif ( ~isempty(strfind(flds{vix},'_field')) )
      ixes = find(-4 >= stn.(fld).field | stn.(fld).field >= 40);
      stn.(fld).field(ixes) = nan;

      ixes = find(alldts(1) <= stn.(fld).date & stn.(fld).date <= alldts(end));
      stn.(fld).date = stn.(fld).date(ixes);
      stn.(fld).field = stn.(fld).field(ixes,:,:);
    else
      ixes = find(-3000 <= stn.(fld).data & stn.(fld).data <= 3000);
      stn.(fld).date = stn.(fld).date(ixes);
      stn.(fld).data = stn.(fld).data(ixes);

      ixes = find(alldts(1) <= stn.(fld).date & stn.(fld).date <= alldts(end));
      stn.(fld).date = stn.(fld).date(ixes);
      stn.(fld).data = stn.(fld).data(ixes);
    end;
  end;

  set_more;

return;
