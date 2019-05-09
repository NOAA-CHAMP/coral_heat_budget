function stn = get_ncep_station(stn, dataset)
%function stn = get_ncep_station(stn, dataset)
%
% Load NCEP reanalysis MAT files (saved e.g., by INCREMENTALLY_SAVE_NCEP_NARR
% or similar, after calls to e.g., QUERY_NCEP_NARR_SUBSET) from reanalysis
% DATASET, and append data from them to individual station data struct STN.
%
% Last Saved Time-stamp: <Thu 2010-08-05 12:52:25  Lew.Gramer>

  if ( isstruct(stn) )
    stnm = lower(stn.station_name);
  elseif ( ischar(stn) )
    stnm = lower(stn);
    stn = [];
  else
    error('First arg must be a station STRUCT or name STRING!');
  end;

  datapath = get_thesis_path('../data');

  ncepmatfname = fullfile(datapath,sprintf('%s_ncep_%s.mat',stnm,dataset));
  if ( exist(ncepmatfname,'file') )
    dat = load(ncepmatfname);
    disp(['Loading pre-merged ' ncepmatfname]);
    stn = get_ncep_station_merge_fields(stn, dat.station, dataset);
    dat = []; clear dat;

    if (isfield(stn,'ncep_air_t'))
      if ( nanmean(stn.ncep_air_t.data) > 273 )
        % Convert units [K]->[oC], [Pa]->[hPa]
        stn.ncep_air_t.data = stn.ncep_air_t.data - 273.15;
        if (isfield(stn,'ncep_sea_t')); stn.ncep_sea_t.data = stn.ncep_sea_t.data - 273.15; end;
        if (isfield(stn,'ncep_barom')); stn.ncep_barom.data = stn.ncep_barom.data ./ 100; end;

        disp(['REsaving converted units to ' ncepmatfname '...']);
        station = stn;
        save(ncepmatfname,'station');
        station = []; clear station;
      end;
    end;

  else

    ncep = [];
    anySoFar = false;

    disp(['Loading year-month NCEP ' dataset ' .MAT files...']);
    for yr = 1987:2011
      % Accident of history - sometimes I saved whole years, sometimes individual months
      fname = fullfile(datapath,'ncep',sprintf('ncep_%s_%04d.mat',dataset,yr));
      if ( exist(fname,'file') )
        anySoFar = true;
        disp(fname);
        load(fname);
        ncep = get_ncep_station_mash_all_fields(ncep, dat.(stnm));
        dat = []; clear dat;
      else
        for mo = 1:12
          dt = datenum(yr,mo,1);
          % NCEP Reanalyses do not allow us to look into the future
          if (dt > now)
            break;
          end;
          fname = fullfile(datapath,'ncep',sprintf('ncep_%s_%04d_%02d.mat',dataset,yr,mo));
          if ( exist(fname,'file') )
            anySoFar = true;
            disp(fname);
            load(fname);
            ncep = get_ncep_station_mash_all_fields(ncep, dat.(stnm));
            dat = []; clear dat;
          else
            if ( anySoFar )
              warning('No file found! "%s"', fname);
            else
              % Ignore error for dates prior to the FIRST data of this dataset
            end;
          end; %if exist(fname) else
        end; %for mo
      end; %if exist(fname) else
    end; %for yr

    % Interpolate all fields to hourly values
    stn = get_ncep_station_merge_fields(stn, ncep, dataset);

    % Convert units [K]->[oC], [Pa]->[hPa]
    if (isfield(stn,'ncep_air_t')); stn.ncep_air_t.data = stn.ncep_air_t.data - 273.15; end;
    if (isfield(stn,'ncep_sea_t')); stn.ncep_sea_t.data = stn.ncep_sea_t.data - 273.15; end;
    if (isfield(stn,'ncep_barom')); stn.ncep_barom.data = stn.ncep_barom.data ./ 100; end;

    % Calculate some derived fields

    stn.ncep_par.date = stn.ncep_dsrf.date;
    stn.ncep_parW.date = stn.ncep_dsrf.date;
    [stn.ncep_par.data, stn.ncep_parW.data] = ...
        insol_to_par(stn.ncep_dsrf.data);

    [ix1,ix2] = intersect_dates(stn.ncep_dsrf.date,stn.ncep_usrf.date);
    stn.ncep_srf.date = stn.ncep_dsrf.date(ix1);
    stn.ncep_srf.data = stn.ncep_dsrf.data(ix1) - stn.ncep_usrf.data(ix2);
    stn.ncep_lrf.date = stn.ncep_dlrf.date;
    stn.ncep_lrf.data = stn.ncep_dlrf.data - stn.ncep_ulrf.data;

    [ix1,ix2] = intersect_dates(stn.ncep_srf.date,stn.ncep_lrf.date);
    stn.ncep_net_heat_flux.date = stn.ncep_lrf.date(ix2);
    stn.ncep_net_heat_flux.data = stn.ncep_srf.data + stn.ncep_lrf.data + ...
        stn.ncep_latent_heat_flux.data + stn.ncep_sensible_heat_flux.data;

    stn = calc_ncep_wind(stn);

    disp(['Saving to ' ncepmatfname '...']);
    station = stn;
    save(ncepmatfname,'station');
    station = []; clear station;

  end; %if exist(ncepmatfname) else

return;


%%%%%%%%%%
%%%%%%%%%% PRIVATE FUNCTIONS
%%%%%%%%%%


% Just combine all fields as loaded from each year-month MAT file, into one
% master struct NCEP for this station - we will merge with STN later
function ncep = get_ncep_station_mash_all_fields(ncep,dat)
  if ( isempty(ncep) )
    ncep = dat;
  else
    flds = fieldnames(dat);
    for fldix = 1:length(flds)
      fld = flds{fldix};
      if ( ~isfield(ncep, fld) )
        ncep.(fld).date = [];
        ncep.(fld).data = [];
      end;
      if ( isfield(dat.(fld),'data') )
        ndat = numel(dat.(fld).data);
        ncep.(fld).date(end+1:end+ndat,1) = dat.(fld).date;
        ncep.(fld).data(end+1:end+ndat,1) = dat.(fld).data;
      end;
    end;
  end;
return;


% Merge data from master NCEP struct with station struct STN
function stn = get_ncep_station_merge_fields(stn,ncep,dataset)

  flds = fieldnames(ncep);
  flds = flds(strmatch('ncep_',flds));
  % Another accident of history
  if ( any(~strmatch(flds, ['ncep_' dataset '_'])) )
    newflds = strrep(flds,'ncep_',['ncep_' dataset '_' ]);
  else
    newflds = flds;
  end;
  for fldix = 1:length(flds)
    fld = flds{fldix};
    newfld = newflds{fldix};
    if ( ~isfield(stn, newfld) )
      stn.(newfld).date = [];
      stn.(newfld).data = [];
    end;
    dts = ncep.(fld).date(1):(1/24):ncep.(fld).date(end);
    ndat = numel(dts);
    stn.(newfld).date(end+1:end+ndat,1) = dts;
    goodix = find(isfinite(ncep.(fld).data));
    stn.(newfld).data(end+1:end+ndat,1) = ...
        interp1(ncep.(fld).date(goodix),ncep.(fld).data(goodix),dts);
  end;

return;
