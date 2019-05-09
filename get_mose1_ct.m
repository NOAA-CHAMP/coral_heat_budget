function stn = get_mose1_ct(stn)
%function stn = get_mose1_ct(stn)
%
% Load multiple raw CT ASCII data files from Moser Channel (SFP) deployment,
% and concatenate into a group of time series fields in STN for sea
% temperature, conductiviy, salinity, and density (assuming p=4.9dbar).
%
% Last Saved Time-stamp: <Sat 2010-11-20 21:34:23 Eastern Standard Time gramer>

  %moserpath = 'C:\Documents and Settings\gramer\My Documents\RSMAS\Coastal\thesis\data\SFP';
  datapath = get_thesis_path('../data');
  moserpath = fullfile(datapath,'SFP');

  stn.ct_fnames = { ...
      'MOSER_CHANNEL_CT_2004.asc', ...
      'MOSER_CHANNEL_CT_2005.asc', ...
      'MOSER_CHANNEL_CT_2006.asc', ...
      'MOSER_CHANNEL_CT_2007.asc', ...
      'MOSER_CHANNEL_CT_2008.asc', ...
           };

  matfname = fullfile(datapath,'mose1_ct.mat');
  if ( exist(matfname,'file') )
    result = load(matfname,'stn');

    if ( isfield(result.stn,'ct_fnames') ...
         && all(size(result.stn.ct_fnames)==size(stn.ct_fnames)) ...
         && all(strcmp(result.stn.ct_fnames,stn.ct_fnames)) )

      disp(['Reloading previously parsed data from MAT file: ' matfname]);
      flds = fieldnames(result.stn);
      for fldix = 1:length(flds)
        fld = flds{fldix};
        stn.(fld) = result.stn.(fld);
      end;
      result = []; clear result;
      %%%%%%%%%%%%%%%% EARLY RETURN
      return;
      %%%%%%%%%%%%%%%% EARLY RETURN

    end;

  end;


  disp(['Parsing raw ASCII CT data...']);

  stn.ct_file_begdates = [];
  stn.ct_file_enddates = [];
  stn.ct_seatemp = struct('date',[],'data',[]);
  stn.ct_cond = struct('date',[],'data',[]);
  stn.ct_salin = struct('date',[],'data',[]);
  stn.ct_dens = struct('date',[],'data',[]);

  for fnameix=1:length(stn.ct_fnames)
    fname = fullfile(moserpath,stn.ct_fnames{fnameix});
    %DEBUG:    disp(fname);

    rawdat = load(fname);

    dy = rawdat(:,4);
    mo = rawdat(:,5);
    yr = rawdat(:,6);
    hr = rawdat(:,7);
    mn = rawdat(:,8);

    c = rawdat(:,11);
    t = rawdat(:,12);
    s = rawdat(:,13);

    % Convert conductivity from Siemens/m to mS/cm for consistency
    c = c * 10;

    yr(yr < 70) = yr(yr < 70) + 2000;
    yr(yr < 100) = yr(yr < 100) + 1900;
    dts = datenum(yr,mo,dy,hr,mn,0);

    % Basic QA/QC
    badix = find(datenum(1980,1,1)>dts | dts>datenum(2020,1,1) ...
                 | 2>t | t>38 | 40>c | c>70 | 20>s | s>40);
    if ( ~isempty(badix) )
      %DEBUG:
      disp(['Removing ' num2str(length(badix)) ' bad points: ' fname]);
      dts(badix) = [];
      c(badix) = [];
      t(badix) = [];
      s(badix) = [];
    end;

    if ( isempty(dts) )
      warning('No valid data found?? In "%s"',fname);
      continue;
    end;

    stn.ct_file_begdates(end+1) = dts(1);
    stn.ct_file_enddates(end+1) = dts(end);

    ndts = length(dts);


    % Calculate salinity in PSU and density using (*obsolete*) SeaWater toolkit
    % Instrument at ~16 feet ~ 4.9m ~ 4.9dbar
    p = 4.9;

    d = sw_dens(s,t,p);

    % Append data from this file to STN struct fields
    stn.ct_seatemp.date(end+1:end+ndts,1) = dts(:);
    stn.ct_seatemp.data(end+1:end+ndts,1) = t(:);
    stn.ct_cond.date(end+1:end+ndts,1) = dts(:);
    stn.ct_cond.data(end+1:end+ndts,1) = c(:);
    stn.ct_salin.date(end+1:end+ndts,1) = dts(:);
    stn.ct_salin.data(end+1:end+ndts,1) = s(:);
    stn.ct_dens.date(end+1:end+ndts,1) = dts(:);
    stn.ct_dens.data(end+1:end+ndts,1) = d(:);
  end;


  % For files with overlapping dates, pick old file over new one
  badix = [];
  dts = stn.ct_seatemp.date;
  overlapix = find(diff(dts) <= 0);
  for ix = 1:length(overlapix)
    lastgood = dts(overlapix);
    endix = find((dts-(0.4/24) > lastgood),1);
    dupix = (overlapix+1):(endix-1);
    badix = unique([badix(:) ; dupix(:)]);
  end;

  if ( ~isempty(badix) )
    %DEBUG:
    disp(['Removing ' num2str(length(badix)) ' overlapping points']);

    stn.ct_seatemp.date(badix) = [];
    stn.ct_seatemp.data(badix) = [];
    stn.ct_cond.date(badix) = [];
    stn.ct_cond.data(badix) = [];
    stn.ct_salin.date(badix) = [];
    stn.ct_salin.data(badix) = [];
    stn.ct_dens.date(badix) = [];
    stn.ct_dens.data(badix) = [];
  end;

  disp(['Saving parsed CT data to MAT file: ' matfname]);
  save(matfname,'stn');

return;
