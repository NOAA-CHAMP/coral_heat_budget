function stn = get_ww3_station(stn)
%function stn = get_ww3_station(stn)
%
% Load NOAA WaveWatch III MAT files (saved e.g., by EXTRACT_WW3_MULTI and
% EXTRACT_WW3 or similar calling, e.g., netCDF-Java Toolbox and LOADDAP),
% and append data from them to fields in single station data struct STN. If
% STN is a string, assume it is a station 5-character code.
%
% Last Saved Time-stamp: <Thu 2012-03-22 11:00:17  Lew.Gramer>

  set_more off;

  if ( isstruct(stn) )
    stnm = lower(stn.station_name);
  elseif ( ischar(stn) )
    stnm = lower(stn);
    stn = [];
  else
    error('First arg must be a station STRUCT or name STRING!');
  end;

  datapath = get_thesis_path('../data');
  %ww3path = '\\cygnus\gramer\backup\RSMAS\data';
  ww3path = fullfile(datapath,'ww3');

  ww3matfname = fullfile(datapath,sprintf('%s_ww3.mat',stnm));
  if ( exist(ww3matfname,'file') )
    disp(['Loading pre-merged ' ww3matfname]);
    dat = load(ww3matfname);
    stn = get_ww3_station_merge_fields(stn, dat.station);
    dat = []; clear dat;

  else

    %DEBUG:
    tic,

    % List of datasets to merge together, in order from WORST (most likely to
    % be overwritten by overlapping data from later datasets) to BEST
    datasets = {'wna', 'at_4m'};

    % Identical for both datasets
    vars = {'tp','dp','hs'};
    flds = { 'peakwaveper','peakwavedir','sigwavehgt', };

    for dix = 1:length(datasets)
      dataset = datasets{dix};

      ww3 = [];
      anySoFar = false;

      disp(['Loading year-month WW3 ' dataset ' .MAT files...']);
      yrs = 1987:get_year(now);
      for yr = yrs(:)'
        mos = 1:12;
        if ( yr == get_year(now) )
          mos = 1:(get_month(now));
        end;
        for mo = mos(:)'
          dt = datenum(yr,mo,1);
          if (dt > now)
            break;
          end;
          fname = fullfile(ww3path,sprintf('ww3.%s.%04d%02d.mat',dataset,yr,mo));
          if ( exist(fname,'file') )
            anySoFar = true;
            disp(['Merging ' fname]);
            load(fname);
            if ( ~isfield(dat,stnm) )
              warning('Ecoforecasts:WW3:NoStation',...
                      'Station %s not found in "%s"! Rerunning EXTRACTION...',...
                      stnm,fname);
              dat.(stnm).lon = stn.lon;
              dat.(stnm).lat = stn.lat;
              [dat.(stnm).lonix,dat.(stnm).latix] = ...
                  gridnbhd_km(dat.lon, dat.lat,dat.(stnm).lon, dat.(stnm).lat, 0);
              for vix = 1:length(vars)
                var = vars{vix};
                fld = ['ww3_' dataset '_' flds{vix}];
                dat.(stnm).(fld).date = dat.date(:);
                dat.(stnm).(fld).data = squeeze(dat.(var)(:,dat.(stnm).latix,dat.(stnm).lonix));
                %DEBUG:            {stnm,fld,nanmean(dat.(stnm).(fld).data)}
              end; %for vix
              disp(['REsaving ' fname]);
              save(fname,'dat');  
            end; %if ( ~isfield(dat,stnm) )
            ww3 = get_ww3_station_mash_all_fields(ww3, dat.(stnm),dataset);
            dat = []; clear dat;
          elseif ( anySoFar && dix > 1 )
            disp(['Missing year-month ' fname]);
          end; %if exist(fname) elseif
        end; %for mo
      end; %for yr

      if ( ~anySoFar )

        warning('Ecoforecasts:WW3:NoFiles','No files found for dataset "%s"!', dataset);

      else

        % Also calculate peak wave vector components "U" and "V"
        [Tix,Dix] = intersect_dates(ww3.ww3_peakwaveper.date,ww3.ww3_peakwavedir.date);
        T = ww3.ww3_peakwaveper.data(Tix);
        D = ww3.ww3_peakwavedir.data(Dix);
        % Calculate "u" and "v" peak wavenumber vector components by weighting
        % PEAKWAVEDIR with putative (deep-water) wave celerity: thus. the lower
        % frequency waves dominate when averaged in with shorter wavelengths.
        % (As for any vector, interpolating PEAKWAVEDIR means special handling!)
        if ( isfield(stn,'lat') ); lat=stn.lat; else lat=25; end;
        g = sw_g(lat,0);
        c = g .* T ./ (2*pi);
        % % *OR* just do a simple vector averaging with all weights == 1??
        % c = 1;
        [u,v] = spddir_to_uv(c,D);
        ww3.ww3_peakwave_u.date = ww3.ww3_peakwavedir.date(Dix);
        ww3.ww3_peakwave_u.data = u;
        ww3.ww3_peakwave_v.date = ww3.ww3_peakwavedir.date(Dix);
        ww3.ww3_peakwave_v.data = v;

        % Interpolate all fields to hourly values
        stn = get_ww3_station_merge_fields(stn, ww3, dataset);

      end; %if ~anySoFar else

    end; %for dataset


    disp(['Saving to ' ww3matfname '...']);
    station = stn;
    save(ww3matfname,'station');
    station = []; clear station;

    %DEBUG:
    toc,

  end; %if exist(ww3matfname) else

  set_more;

return;


%%%%%%%%%%
%%%%%%%%%% PRIVATE FUNCTIONS
%%%%%%%%%%


%%% GET_WW3_STATION_MASH_ALL_FIELDS
% Just combine all fields as loaded from each year-month MAT file, into one
% master struct WW3 for this dataset and station - will merge with STN later
function ww3 = get_ww3_station_mash_all_fields(ww3,dat,dataset)

  flds = fieldnames(dat);
  if ( ~exist('dataset','var') || isempty(dataset) )
    newflds = flds;
  else
    newflds = strrep(flds,['ww3_' dataset '_' ],'ww3_');
  end;

  % Peak Period below 2s is NOT REAL - just means a model was spinning up!
  wvper = ['ww3_' dataset '_peakwaveper'];
  if ( isfield(dat,wvper) )
    badix = find(dat.(wvper).data<2);
    if ( ~isempty(badix) )
      nwvper = length(dat.(wvper).data);
      for fldix = 1:length(flds)
        fld = flds{fldix};
        if ( is_ts(dat.(fld)) )
          if ( length(dat.(fld).data) ~= nwvper )
            error('Fields differ in length?? %s vs. %s',fld,wvper);
          end;
          dat.(fld).date(badix) = [];
          dat.(fld).data(badix) = [];
        end;
      end;
    end;
  end;

  for fldix = 1:length(flds)
    fld = flds{fldix};
    newfld = newflds{fldix};
    if ( ~isfield(ww3, newfld) )
      ww3.(newfld).date = [];
      ww3.(newfld).data = [];
    end;
    if ( isfield(dat.(fld),'data') )
      ndat = numel(dat.(fld).data);
      ww3.(newfld).date(end+1:end+ndat,1) = dat.(fld).date;
      ww3.(newfld).data(end+1:end+ndat,1) = dat.(fld).data;
    end;
  end;

return;


%%% GET_WW3_STATION_MERGE_FIELDS
function stn = get_ww3_station_merge_fields(stn,ww3,dataset)
%function stn = get_ww3_station_merge_fields(stn,ww3,dataset)
% Merge data from master WW3 struct with station struct STN:
% Copy all fields at native time resolution from extracted struct WW3 to new
% struct STN, with "raw_" prefix added to them. Remove any DATASET infixes,
% merging overlapping dates by choosing more recent dataset FIRST. Finally,
% interpolate to hourly time series with new field names.

  flds = fieldnames(ww3);
  flds = flds(strmatch('ww3_',flds));
  if ( ~exist('dataset','var') || isempty(dataset) )
    newflds = flds;
  else
    newflds = strrep(flds,['ww3_' dataset '_' ],'ww3_');
  end;

  newdat = [];
  for fldix = 1:length(flds)
    fld = flds{fldix};
    newfld = newflds{fldix};
    rawfld = ['raw_' newfld];

    if ( isempty(find(isfinite(ww3.(fld).data))) )
      warning('Ecoforecasts:WW3:NoValidMerge','No valid values in field "%s"!',fld);
      newdat.(rawfld) = struct('date',[],'data',[]);
      newdat.(newfld) = struct('date',[],'data',[]);
    else
      newdat.(rawfld) = ww3.(fld);
      if ( strcmp(newfld,'ww3_sigwavehgt') )
        % NOTE: Negative values are *never* valid for raw wave model fields
        % Also note: Not quite sure where these negative values crept in...
        newdat.(rawfld).data(newdat.(rawfld).data<eps) = 0;
      end;
      newdat.(newfld) = interp_ts(newdat.(rawfld));
    end;
  end;

  % Some basic quality control (SHOULD have been done during initial extraction!)
  if ( isfield(newdat,'ww3_peakwaveper') )
    % Peak Period below 2s is NOT real - just means a model was spinning up!
    baddts = newdat.ww3_peakwaveper.date(newdat.ww3_peakwaveper.data<2);
    if ( ~isempty(baddts) )
      [badix,ig] = intersect_dates(newdat.ww3_sigwavehgt.date,baddts);
      % newdat.ww3_sigwavehgt.data(badix) = nan;
      newdat.ww3_sigwavehgt.date(badix) = [];
      newdat.ww3_sigwavehgt.data(badix) = [];
      [badix,ig] = intersect_dates(newdat.ww3_peakwavedir.date,baddts);
      % newdat.ww3_peakwavedir.data(badix) = nan;
      newdat.ww3_peakwavedir.date(badix) = [];
      newdat.ww3_peakwavedir.data(badix) = [];
      [badix,ig] = intersect_dates(newdat.ww3_peakwaveper.date,baddts);
      % newdat.ww3_peakwaveper.data(badix) = nan;
      newdat.ww3_peakwaveper.date(badix) = [];
      newdat.ww3_peakwaveper.data(badix) = [];

      if ( isfield(newdat,'ww3_peakwave_u') && isfield(newdat,'ww3_peakwave_v') )
        [badix,ig] = intersect_dates(newdat.ww3_peakwave_u.date,baddts);
        % newdat.ww3_peakwave_u.data(badix) = nan;
        newdat.ww3_peakwave_u.date(badix) = [];
        newdat.ww3_peakwave_u.data(badix) = [];
        [badix,ig] = intersect_dates(newdat.ww3_peakwave_v.date,baddts);
        % newdat.ww3_peakwave_v.data(badix) = nan;
        newdat.ww3_peakwave_v.date(badix) = [];
        newdat.ww3_peakwave_v.data(badix) = [];
      end;
    end;
  end;

  % As for any vector, interpolating PEAKWAVEDIR means special handling!
  if ( isfield(newdat,'ww3_peakwave_u') && isfield(newdat,'ww3_peakwave_v') )
    if ( numel(find(isfinite(newdat.ww3_peakwave_u.data))) < 1 || ...
         numel(find(isfinite(newdat.ww3_peakwave_v.data))) < 1 )

      warning('Ecoforecasts:WW3:NoUV','No valid peakwave_u/v found for "%s"!', dataset);
      newdat.ww3_peakwavedir = struct('date',[],'data',[]);

    else
      [uix,vix] = intersect_dates(newdat.ww3_peakwave_u.date,newdat.ww3_peakwave_v.date);
      u = newdat.ww3_peakwave_u.data(uix);
      v = newdat.ww3_peakwave_v.data(vix);
      D = uv_to_dir(u,v);
      newdat.ww3_peakwavedir.date = newdat.ww3_peakwave_u.date(uix);
      newdat.ww3_peakwavedir.data = D;
    end;
  end;

  if ( ~isempty(newdat) )
    warning('off','Ecoforecasts:mergedNonTS');
    stn = merge_station_data(stn,newdat);
    warning('on','Ecoforecasts:mergedNonTS');
  end;
  newdat = [];
  clear newdat;

return;
