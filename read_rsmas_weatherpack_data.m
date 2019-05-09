function stn = read_rsmas_weatherpack_data(stn_or_stnm)
%function stn = read_rsmas_weatherpack_data(stn_or_stnm)
%
% Load RSMAS Weatherpack data (see thesis/data/RSMAS/README.txt)
%
% Last Saved Time-stamp: <Tue 2011-12-27 12:02:16  lew.gramer>

  set_more off;
  tic,

  datapath = get_thesis_path('../data');
  rsmaspath = fullfile(datapath,'RSMAS');

  % Weatherpak Header
  % DATE          YY/MM/DD UTC
  % TIME          HH:MM:SS
  % WSPD#1        Wind speed #1 (m/s)
  % WDIR#1        Wind Direction #1
  % SD#1          
  % GUST#1        Wind gust #1 (m/s)
  % WSPD#2        Wind speed #2 (m/s)
  % WDIR#2        Wind direction #2
  % SD#2          
  % GUST#2        Wind gust #2 (m/s)
  % AVETEMP       Temperature (°C)
  % INSTTEMP      Temperature (°C)
  % AVERH         Humidity (%)
  % AVEBARO       Barometric Pressure
  % PIR           Incoming LW (W/m2)
  % PIRCASE               
  % PIRDOME               
  % PSP           Incoming SW (W/m2)
  % RAINRATE      Rain rate (mm/hr)
  % RAINACC       Accumulated Rain (mm)
  % TB-RR         
  % TB-RA         
  % VBATT         
  flds = {'','',...
          'wind1_speed','wind1_dir','','wind1_gust',...
          'wind2_speed','wind2_dir','','wind2_gust',...
          'air_t','','relhumid','barom',...
          'dlrf','','','dsrf','precip','precip_acc',...
          '','','',...
         };

  if ( nargin < 1 || isempty(stn_or_stnm) )
    stn = [];
  else
    stn = get_station_from_station_name(stn_or_stnm);
  end;

  matfname = fullfile(datapath,'rsmas_weatherpack.mat');
  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname,'result');
    %DEBUG:    disp('Loaded');

  else

    disp('Extracting from raw Weatherpack files...');
    result = [];

    firstyear = 2003;    firstjday = 237;

    % % Assume we have downloaded up through yesterday
    % lastyear = get_year(now-1);    lastjday = get_jday(now-1);

    % ??? Assume we have downloaded through 2011 Dec 07
    lastyear  = 2011;    lastjday  = 341;

    %DEBUG:    firstyear = 2011;    firstjday = 300;
    %DEBUG:    lastyear  = 2003;    lastjday  = 365;

    for yr = firstyear:lastyear

      %DEBUG:
      disp(yr);

      if ( mod(yr,4) == 0 )
        begjd=1; endjd=366;
      else
        begjd=1; endjd=365;
      end;
      if ( yr == firstyear )        begjd=firstjday;      end;
      if ( yr == lastyear )         endjd=lastjday;       end;

      for jd = begjd:endjd
        fname = fullfile(rsmaspath,sprintf('zena.dat.%d.%03d',yr,jd));
        if ( ~exist(fname,'file') )
          warning('Missing %s',fname);
          continue;
        end;

        fid=fopen(fname,'r');
        if ( fid < 1 )
          warning('Unable to open %s',fname);
          continue;
        end;
        C=textscan(fid,['%[^ ] %[^ ] %f %f %f %f %f %f %f %f %f %f %f %f %f ' ...
                        '%f %f %f %f %f %f %f %f\n']);
        fclose(fid);

        if ( numel(C) ~= numel(flds) )
          C=[]; clear C
          warning('Invalid format %s',fname);
          continue;
        end;

        dts = datenum(strcat(C{1},C{2}),'yy/mm/ddHH:MM:SS');
        for fldix=1:numel(flds)
          fld = flds{fldix};
          if ( ~isempty(fld) )
            fld = ['rsmas_',fld];
            rawfld = ['raw_',fld];
            dat = C{fldix};
            res.(rawfld).date = dts(:);
            res.(rawfld).data = dat(:);
          end;
        end;
        dts=[]; clear dts
        C=[]; clear C

        result = merge_station_data(result,res);
        res=[]; clear res

      end; %for jd = begjd:endjd

    end; %for yr = firstyear:lastyear

    % Calculate hourly averages of all fields
    for fldix=1:numel(flds)
      fld = flds{fldix};
      if ( ~isempty(fld) )
        fld = ['rsmas_',fld];
        rawfld = ['raw_',fld];
        %DEBUG:      disp(fld);
        result.(fld) = interp_ts(result.(rawfld));
      end;
    end;

    disp(['Saving ',matfname]);
    save(matfname,'result');

  end; %if ( exist(matfname,'file') ) else


  % Transfer results to our return struct
  stn = merge_station_data(stn,result);
  result = []; clear result
  %DEBUG:  disp('Merged');

  if ( ~isfield(stn,'station_name') )
    stn.station_name = 'rsmas';
  end;

  toc,
  set_more;

return;
