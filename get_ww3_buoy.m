function stn = get_ww3_buoy(stn)
%function stn = get_ww3_buoy(stn)
%
% Load NOAA WaveWatch III MAT files (saved e.g., by EXTRACT_WW3_MULTI and
% EXTRACT_WW3 or similar calling, e.g., netCDF-Java Toolbox and LOADDAP),
% and extract data for a new buoy location (i.e., a monitoring site not
% previously subsetted by the EXTRACT* functions), and append those data to
% fields in single station data struct STN. If STN is a string, assume it is
% a buoy 5-number code, and *prepend a 'b' to it* for use as a field name.
%
% Last Saved Time-stamp: <Thu 2011-01-13 13:12:38  lew.gramer>

  set_more off;

  if ( isstruct(stn) )
    stnm = lower(stn.station_name);
  elseif ( ischar(stn) )
    stnm = lower(stn);
    stn = [];
  else
    error('First arg must be a station STRUCT or name STRING!');
  end;

  buoynm = stnm;
  if ( ~isempty(str2num(stnm(1))) )
    stnm = [ 'b' stnm ];
  end;

  datapath = get_thesis_path('../data');
  arcdatapath = fullfile(datapath,'ww3');
  %arcdatapath = '\\cygnus\gramer\backup\RSMAS\data';

  % tp - Peak Wave Period
  % dp - Peak Wave Direction
  % hs - Significant Wave Height
  vars = { 'tp',         'dp',         'hs', };
  flds = { 'peakwaveper','peakwavedir','sigwavehgt', };

  % List of datasets to merge together, in order from WORST (most likely to
  % be overwritten by overlapping data from later datasets) to BEST
  datasets = {'wna', 'at_4m'};

  ww3matfname = fullfile(datapath,sprintf('%s_ww3.mat',buoynm));
  if ( exist(ww3matfname,'file') )
    dat = load(ww3matfname);
    disp(['Loading pre-merged ' ww3matfname]);
    stn = get_ww3_station_merge_fields(stn, dat.station);
    dat = []; clear dat;

  else

    %DEBUG:
    tic,

    for dix = 1:length(datasets)
      dataset = datasets{dix};

      ww3 = [];
      anySoFar = false;

      disp(['Loading year-month WW3 ' dataset ' .MAT files...']);
      for yr = 1987:2011
        for mo = 1:12
          dt = datenum(yr,mo,1);
          if (dt > now)
            break;
          end;
          fname = fullfile(arcdatapath,sprintf('ww3.%s.%04d%02d.mat',dataset,yr,mo));
          if ( exist(fname,'file') )
            anySoFar = true;
            disp(['Merging from ' fname]);
            load(fname);
            if ( ~isfield(dat,stnm) )
              [dat.(stnm).lon,dat.(stnm).lat,ig] = get_station_coords(buoynm);
              [dat.(stnm).lonix,dat.(stnm).latix] = ...
                  gridnbhd_km(dat.lon, dat.lat,dat.(stnm).lon, dat.(stnm).lat, 0);
              for vix = 1:length(vars)
                var = vars{vix};
                fld = ['ww3_' flds{vix}];
                dat.(stnm).(fld).date = dat.date(:);
                dat.(stnm).(fld).data = squeeze(dat.(var)(:,dat.(stnm).latix,dat.(stnm).lonix));
              end; %for vix
            end;
            ww3 = get_ww3_station_mash_all_fields(ww3, dat.(stnm),dataset);
            dat = []; clear dat;
          end; %if exist(fname) else
        end; %for mo
      end; %for yr

      if ( ~anySoFar )

        warning('No files found for dataset "%s"!', dataset);

      else

        % Also calculate peak wave vector components "U" and "V"
        [Tix,Dix] = intersect_dates(ww3.ww3_peakwaveper.date,ww3.ww3_peakwavedir.date);
        T = ww3.ww3_peakwaveper.data(Tix);
        D = ww3.ww3_peakwavedir.data(Dix);
        % Calculate "u" and "v" peak wavenumber vector components by weighting
        % PEAKWAVEDIR with putative (deep-water) wave celerity: thus. the lower
        % frequency waves dominate when averaged in with shorter wavelengths.
        % (As for any vector, interpolating PEAKWAVEDIR means special handling!)
        g = grv(25);
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


% Just combine all fields as loaded from each year-month MAT file, into one
% master struct WW3 for this dataset and station - will merge with STN later
function ww3 = get_ww3_station_mash_all_fields(ww3,dat,dataset)

  flds = fieldnames(dat);
  if ( ~exist('dataset','var') || isempty(dataset) )
    newflds = flds;
  else
    newflds = strrep(flds,['ww3_' dataset '_' ],'ww3_');
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


% Merge data from master WW3 struct with station struct STN
function stn = get_ww3_station_merge_fields(stn,ww3,dataset)

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
    newdat.(newfld).date = ...
        [ww3.(fld).date(1):(1/24):ww3.(fld).date(end)]';
    goodix = find(isfinite(ww3.(fld).data));
    newdat.(newfld).data = ...
        interp1(ww3.(fld).date(goodix),ww3.(fld).data(goodix),newdat.(newfld).date);
  end;

  % As for any vector, interpolating PEAKWAVEDIR means special handling!
  if ( isfield(newdat,'ww3_peakwave_u') && isfield(newdat,'ww3_peakwave_v') )
    [uix,vix] = intersect_dates(newdat.ww3_peakwave_u.date,newdat.ww3_peakwave_v.date);
    u = newdat.ww3_peakwave_u.data(uix);
    v = newdat.ww3_peakwave_v.data(vix);
    D = uv_to_dir(u,v);
    newdat.ww3_peakwavedir.date = newdat.ww3_peakwave_u.date(uix);
    newdat.ww3_peakwavedir.data = D;
  end;

  if ( ~isempty(newdat) )
    stn = merge_station_data(stn,newdat);
  end;
  newdat = [];
  clear newdat;

return;
