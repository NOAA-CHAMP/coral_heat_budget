function stn = get_looe1_microcat(stn)
%function stn = get_looe1_microcat(stn)
%
% Load multiple raw MicroCat ASCII data files from Looe Key Reef (SFP)
% deployment, and concatenate into a group of time series fields in STN for
% sea temperature, conductiviy, salinity, and density (assuming p=4.9dbar).
%
% See file THESIS/data/SFP/SERIAL_NUMBERS.txt for dates and serial #s.
%
% Last Saved Time-stamp: <Fri 2012-03-23 15:44:40  Lew.Gramer>

  %looepath = 'C:\Documents and Settings\gramer\My Documents\RSMAS\Coastal\thesis\data\SFP';
  datapath = get_thesis_path('../data');
  looepath = fullfile(datapath,'SFP');

  if ( ~exist('stn','var') )
    stn = [];
  end;

  matfname = fullfile(datapath,'looe1_microcat.mat');
  if ( exist(matfname,'file') )
    disp(['Reloading pre-saved file ' matfname]);
    load(matfname,'result');

  else

    disp(['Parsing raw ASCII MicroCat data...']);

    result.microcat_fnames = { ...
        '1116JUL06.ASC', ...
        '2938LooeFEB07(No_Real_Time).asc', ...
        'looeNov07.asc', ...
        '0147JUN08.ASC', ...
        '072JAN2509.asc', ...
        'MC 143 Looe (Recov May 27 2009).asc', ...
        'MC 0152 Looe Recovered May 04 2010.asc', ...
             };

    result.microcat_file_begdates = [];
    result.microcat_file_enddates = [];
    result.microcat_seatemp = struct('date',[],'data',[]);
    result.microcat_cond = struct('date',[],'data',[]);
    result.microcat_salin = struct('date',[],'data',[]);
    result.microcat_dens = struct('date',[],'data',[]);

    for fnameix=1:length(result.microcat_fnames)
      fname = fullfile(looepath,result.microcat_fnames{fnameix});
      %DEBUG:    disp(fname);

      fid = fopen(fname,'r');
      if ( fid < 0 )
        warning('Skipping unopenable file "%s"',fname);
        continue;
      end;
      while ( isempty(ferror(fid)) && ~feof(fid) )
        ln = fgetl(fid);
        lasthdr = strmatch('start sample number =',ln);
        if ( ~isempty(lasthdr) ); break; end;
      end;
      if ( isempty(lasthdr) )
        warning('End of header not found in "%s"',fname);
        fclose(fid);
        continue;
      end;
      C = textscan(fid,'%n%n%s%s', 'Delimiter',',');
      fclose(fid);

      t = C{1};
      c = C{2};
      % Convert conductivity from Siemens/m to mS/cm for consistency
      c = c * 10;

      dt = strcat(C{3}, {' '}, C{4});
      % Filter out misformatted dates that DATENUM would choke on
      % Example of a good date: '03-13-2006 23:00:01'
      matchix = ...
          regexp(dt,...
                 '^[01][0-9][-][0-3][0-9][-][12][90][0-9][0-9] [0-2][0-9]:[0-9][0-9]:[0-9][0-9]$');
      goodix = find(~cellfun(@isempty,matchix));
      t = t(goodix);
      c = c(goodix);
      dt = dt(goodix);
      dts = datenum(dt);

      % Basic QA/QC
      badix = find(datenum(1980,1,1)>dts | dts>datenum(2020,1,1) ...
                   | 2>t | t>38 | 40>c | c>70);
      if ( ~isempty(badix) )
        %DEBUG:
        disp(['Removing ' num2str(length(badix)) ' bad points: ' fname]);
        dts(badix) = [];
        t(badix) = [];
        c(badix) = [];
      end;

      if ( isempty(dts) )
        warning('No valid data found?? In "%s"',fname);
        continue;
      end;

      result.microcat_file_begdates(end+1) = dts(1);
      result.microcat_file_enddates(end+1) = dts(end);

      % If we get a WEIRD FILE with HYPERSAMPLING - group into half-hourly means
      if ( mean(diff(dts)) < (0.4/24) )
        halfhours = round(dts/(0.5/24))*(0.5/24);
        t = grpstats(t,halfhours);
        c = grpstats(c,halfhours);
        dts = unique(halfhours);
      end;

      ndts = length(dts);


      % Calculate salinity in PSU and density using (*obsolete*) SeaWater toolkit
      % Instrument at ~16 feet ~ 4.9m ~ 4.9dbar
      p = 4.9;

      cndr = c ./ sw_c3515();
      s = sw_salt(cndr,t,p);
      d = sw_dens(s,t,p);

      % Append data from this file to RESULT struct fields
      result.microcat_seatemp.date(end+1:end+ndts,1) = dts(:);
      result.microcat_seatemp.data(end+1:end+ndts,1) = t(:);
      result.microcat_cond.date(end+1:end+ndts,1) = dts(:);
      result.microcat_cond.data(end+1:end+ndts,1) = c(:);
      result.microcat_salin.date(end+1:end+ndts,1) = dts(:);
      result.microcat_salin.data(end+1:end+ndts,1) = s(:);
      result.microcat_dens.date(end+1:end+ndts,1) = dts(:);
      result.microcat_dens.data(end+1:end+ndts,1) = d(:);
    end;


    % For files with overlapping dates, pick old file over new one
    badix = [];
    dts = result.microcat_seatemp.date;
    overlapix = find(diff(dts) <= 0);
    for ix = 1:length(overlapix)
      lastgood = dts(overlapix);
      endix = find((dts-(0.4/24) > lastgood),1);
      dupix = (overlapix+1):(endix-1);
      badix = unique([badix(:) ; dupix(:)]);
    end;

    %DEBUG:
    disp(['Removing ' num2str(length(badix)) ' overlapping points']);

    result.microcat_seatemp.date(badix) = [];
    result.microcat_seatemp.data(badix) = [];
    result.microcat_cond.date(badix) = [];
    result.microcat_cond.data(badix) = [];
    result.microcat_salin.date(badix) = [];
    result.microcat_salin.data(badix) = [];
    result.microcat_dens.date(badix) = [];
    result.microcat_dens.data(badix) = [];

    disp(['Saving result to file ' matfname]);
    save(matfname,'result');

  end; %if exist(matfname) else


  flds = fieldnames(result);
  for ix = 1:length(flds)

    fld = flds{ix};
    stn.(fld) = result.(fld);

    % Select an hourly subset for intermittent time series: preferrable to,
    % e.g., hourly spline fit for some analyses, as it uses only REAL DATA
    if ( is_valid_ts(stn.(fld)) )
      hfld = ['hourly_' fld];
      dts = [stn.(fld).date(1):(1/24):stn.(fld).date(end)]';
      actualdts = interp1(stn.(fld).date,stn.(fld).date,dts,'nearest');
      stn.(hfld).date = dts;
      stn.(hfld).data = interp1(stn.(fld).date,stn.(fld).data,dts,'nearest');
      dupix = find(diff(actualdts)<eps) + 1;
      stn.(hfld).date(dupix) = [];
      stn.(hfld).data(dupix) = [];
    end;

  end;

  result = []; clear result;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METADATA FROM EACH OF THE ABOVE FILES - FOR REFERENCE PURPOSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % '1116JUL06.ASC', ...
% start time =  03-13-2006  23:00:01
% sample interval = 1800 seconds
% start sample number = 1
%  22.3072, 0.00092, 03-13-2006, 23:00:01
%  ...
%  22.3206,-0.00033, 08-01-2006, 14:30:00

% % '2938LooeFEB07(No_Real_Time).asc', ...
% start time =  08-03-2006  16:00:01
% sample interval = 16 seconds
% start sample number = 1
%  26.4566, 0.00235, 08-03-2006, 16:00:01
%  ...
%  30.0561, 5.93225, 09-09-2006, 03:38:28

% % 'looeNov07.asc', ...
% start time =  02-16-2007  00:00:01
% sample interval = 1800 seconds
% start sample number = 1
%  19.6257, 0.00002, 02-16-2007, 00:00:01
%  ...
%  21.7218, 0.03043, 10-21-2007, 14:00:01

% % '0147JUN08.ASC', ...
% start time =  10-13-2007  00:00:01
% sample interval = 1800 seconds
% start sample number = 1
%  27.7854, 0.00251, 10-13-2007, 00:00:01
%  ...
%  21.4713, 0.00453, 06-12-2008, 01:00:01

% % '072JAN2509.asc', ...
% start time =  06-10-2008  16:00:01
% sample interval = 1800 seconds
% start sample number = 1
%  20.5166,-0.00002, 06-10-2008, 16:00:01
%  ...
%  17.7067, 0.00063, 01-26-2009, 00:00:01

% % 'MC 143 Looe (Recov May 27 2009).asc', ...
% start time =  01-25-2009  15:00:01
% sample interval = 1800 seconds
% start sample number = 2
%  22.8152,-0.00001, 01-25-2009, 14:49:05
%  ...
%  29.9290, 0.01506, 05-27-2009, 20:00:01

% % 'MC 0152 Looe Recovered May 04 2010.asc', ...
% start time =  05-27-2009  00:00:01
% sample interval = 1800 seconds
% start sample number = 1
%  25.0952, 0.00028, 05-27-2009, 00:00:01
%  ...
%  23.0020, 0.00467, 05-04-2010, 23:30:00
