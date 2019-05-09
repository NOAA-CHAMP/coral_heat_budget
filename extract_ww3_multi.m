function extract_ww3_multi(dataset, region, yrs)
%function extract_ww3_multi(dataset, region, yrs)
%
% Extract NOAA WaveWatch III multi-grid wave model output from GRiB2 files
%
% For details, see NOAA NCEP Web pages, e.g., http://nomad5.ncep.noaa.gov
%
% DEFAULTS: DATASET = 'multi_1.at_4m' (Atlantic coast 4-minute lat/lon resn.)
%           REGION = 'fknms'          (Florida Keys NMS, [24 27 -84 -79])
%           YRS = 2005:NOW            (All available WW3 Multi years)
%
% Last Saved Time-stamp: <Thu 2012-03-08 12:28:13  lew.gramer>

  set_more off;

  datapath = get_thesis_path('../data');
  ww3path = fullfile(datapath, 'ww3');

  if ( ~exist('dataset','var') || isempty(dataset) )
    % Multi-mesh model: Atlantic 4-minute (lat/lon) resolution grid
    dataset = 'multi_1.at_4m';
  end;
  datafld = strrep(dataset,'multi_1.','');

  if ( ~exist('region','var') || isempty(region) )
    region = 'fknms';
  end;

  [thisyear,thismonth,ig] = datevec(now);
  if ( ~exist('yrs','var') || isempty(yrs) )
    yrs = 2005:thisyear;
    %DEBUG:    yrs = 2007;
    %DEBUG:    yrs = 2011:2012;
  end;

  switch ( region ),
   case 'fknms',
    minlat = 24;
    maxlat = 27;
    minlon = -84;
    maxlon = -79;
   otherwise,
    error('Unknown region argument "%s"!', region);
  end;


  stncfg = {};
  cfgfname = fullfile(datapath,['ww3-' datafld '.cfg']);
  if ( exist(cfgfname,'file') )
    fid = fopen(cfgfname,'r');
    stncfg = textscan(fid,'%s%d%d', 'Delimiter',',', 'CommentStyle','#'); 
    fclose(fid);
  end;

  %ftp://polar.ncep.noaa.gov/pub/history/waves/multi_1.at_4m.hs.200711.grb2

  % NOTE: If we add any variables here later, be sure to change the IF
  % statement 'if ( ~isfield(dat,stnm) )' below to read 'if (1)'!
  % tp - Peak Wave Period
  % dp - Peak Wave Direction
  % hs - Significant Wave Height
  vars = { 'tp',         'dp',         'hs', };
  varnms = { ...
      'Primary_wave_mean_period', ...
      'Primary_wave_direction', ...
      'Significant_height_of_combined_wind_waves_and_swell', ...
           };
  flds = { 'peakwaveper','peakwavedir','sigwavehgt', };
  %DEBUG:  vars = { 'hs', };  flds = { 'sigwavehgt', };

  f = [];
  doFTP = false;
  %DEBUG:
  doFTP = true;
  if ( doFTP )
    fhost = 'polar.ncep.noaa.gov';
    fuser = 'anonymous';
    fpawd = 'lew.gramer@noaa.gov';
    fdir = '/pub/history/waves';

    f = ftp(fhost,fuser,fpawd);
    binary(f);
    cd(f, fdir);
  end;


  for yr = yrs(:)'

    switch ( yr ),
     %DEBUG:     case 2007,		mos = 11;
     %DEBUG:     case 2010,		mos = 10:12;
     %DEBUG:     case 2011,		mos = 11:12;
     case thisyear,	mos = 1:thismonth;
     otherwise,		mos = 1:12;
    end;

    for mo = mos(:)'

      matfname = fullfile(ww3path, ...
                          sprintf('ww3.%s.%04d%02d.mat', datafld, yr, mo));
      if ( exist(matfname,'file') )
        disp(['Loading pre-saved MAT file ' matfname]);
        load(matfname,'dat');

      else
        %DEBUG:        tic,

        dat = [];
        clear dat;
        dat = [];

        for vix = 1:length(vars)

          var = vars{vix};
          varnm = varnms{vix};

          fbasename = sprintf('%s.%s.%04d%02d.grb2', dataset, var, yr, mo);
          fname = fullfile(ww3path,fbasename);
                           
          if ( ~isempty(f) )
            if ( ~exist(fname,'file') )
              %DEBUG:
              disp(['Attempting to FTP ' fbasename]);
              try,
                mget(f, fbasename, ww3path);
              catch
                warning('Unable to FTP "%s"',fbasename);
              end;
            end;
          end;

          if ( ~exist(fname,'file') )
            warning('"%s" not available to process',fname);
            continue;
          end;

          %DEBUG:
          disp(fname);

          nc = mDataset(fname);
          if ( isempty(nc) )
            warning('"%s" NJTBX open failed',fname);
            % if ( ~isempty(f) )
            %     delete(fname);
            % end;
            continue;
          end;

          if ( ~exist('lats','var') )
            lats = double( nc{'lat'}(1:end) );
            lons = double( nc{'lon'}(1:end) );
            lons(lons > 180) = lons(lons > 180) - 360;

            % Why do these NCEP meteo-wonks insist on adding a record for
            % the first hour of next month into this month's data file??
            hrs = double( nc{'time'}(1:end-1) );

            % Latitude indices are reversed (naturally!)
            [minlonix,minlatix] = gridnbhd_km(lons,lats,minlon,maxlat,0);
            [maxlonix,maxlatix] = gridnbhd_km(lons,lats,maxlon,minlat,0);
          end;

          if ( ~isfield(dat,'n') )
            dat.n = length(hrs);
            dat.minlonix = minlonix;
            dat.minlatix = minlatix;
            dat.date = datenum(yr,mo,1) + (hrs/24.0);
            dat.lon = lons(minlonix:maxlonix);
            dat.lat = lats(minlatix:maxlatix);
            dat.nlon = length(dat.lon);
            dat.nlat = length(dat.lat);
          end;

          for ix = 1:length(hrs)
            rawdat = double( squeeze(nc{varnm}(ix,1:end,1:end)) );
            dat.(var)(ix,1:dat.nlat,1:dat.nlon) = rawdat(minlatix:maxlatix,minlonix:maxlonix);
            rawdat = []; clear rawdat;
          end;
          dat.(var)(dat.(var) > 360) = nan;

          close(nc);
          clear nc;

          % if ( ~isempty(f) )
          %     delete(fname);
          % end;

        end; %for vix

        if ( ~isfield(dat,vars{1}) )
          warning('Skipping year-month %d-%d!', yr, mo);
          %DEBUG:          toc,
          continue;
        end;

        disp(['Pre-saving to ' matfname]);
        save(matfname,'dat');  
        %DEBUG:        toc,
      end; %if exist(matfname) else

      %DEBUG:      id=62; figure; contourf(squeeze(dat.hs(id,:,:))); maxigraph; caxis([0 3]); colorbar; hold on; set(gca,'yd','rev'); title(['sigwavehgt ' datestr(dat.date(id))]);
      %DEBUG:      id=81; figure; contourf(squeeze(dat.hs(id,:,:))); maxigraph; caxis([0 3]); colorbar; hold on; set(gca,'yd','rev'); title(['sigwavehgt ' datestr(dat.date(id))]);

      %%%%DEBUG:      toc,

      % Subset our world list of stations to those inside our BBOX
      if ( ~exist('stns','var') || isempty(stns) )
        stns = get_all_station_metadata;
        XV = [ min(dat.lon) max(dat.lon) max(dat.lon) min(dat.lon) ];
        YV = [ max(dat.lat) max(dat.lat) min(dat.lat) min(dat.lat) ];
        goodix = find( inside(stns.lons, stns.lats, XV, YV) );
        clear XV YV;
        stns.codes  = stns.codes(goodix);
        stns.lons   = stns.lons(goodix);
        stns.lats   = stns.lats(goodix);
        stns.depths = stns.depths(goodix);
      end;

      %HACK ALERT - should only need to be done ONCE!
      %DEBUG:      minlonix,      minlatix,
      if ( ~isfield(dat,'minlonix') )
        dat.minlonix = 226;
        dat.minlatix = 301;
      end;

      for stix = 1:length(stns.lons)

        stnm = lower(stns.codes{stix});
        % NOTE: If we ever add any variables to the list above, change the
        % following line temporarily to simply read "if (1)"!
        if ( ~isfield(dat,stnm) )
          %DEBUG:
          disp(['Adding ' stnm ' to ' matfname]);
          dat.(stnm).lon = stns.lons(stix);
          dat.(stnm).lat = stns.lats(stix);
          ididx = [];
          if ( ~isempty(stncfg) )
            ididx = find(strcmpi(stnm,stncfg{1}));
          end;
          % [lonix,latix] = gridnbhd_km(dat.lon, dat.lat,...
          %                             dat.(stnm).lon, dat.(stnm).lat, 0);
          if ( ~isempty(ididx) )
            latix = double(stncfg{2}(ididx));
            lonix = double(stncfg{3}(ididx));
            lonix = lonix - dat.minlonix + 1;
            latix = latix - dat.minlatix + 1;
          else
            % [lonix,latix] = gridnbhd_km(lons, lats,...
            %                             dat.(stnm).lon, dat.(stnm).lat, 0);
            [lonix,latix] = gridnbhd_km(dat.lon, dat.lat,...
                                        dat.(stnm).lon, dat.(stnm).lat, 0);
            %DEBUG:            lonix, latix,
            %DEBUG:            fprintf('%s,%d,%d\n',stnm,latix,lonix);
          end;

          % dat.(stnm).lonix = lonix - minlonix + 1;
          % dat.(stnm).latix = latix - minlatix + 1;
          dat.(stnm).lonix = lonix;
          dat.(stnm).latix = latix;

          for vix = 1:length(vars)
            var = vars{vix};
            fld = ['ww3_' datafld '_' flds{vix}];
            dat.(stnm).(fld).date = dat.date(:);
            dat.(stnm).(fld).data = squeeze(dat.(var)(:,dat.(stnm).latix,dat.(stnm).lonix));
            %DEBUG:            {stnm,fld,nanmean(dat.(stnm).(fld).data)}
          end; %for vix

          %DEBUG:        plot(dat.(stnm).lonix,dat.(stnm).latix,'k*'); text(dat.(stnm).lonix,dat.(stnm).latix,stnm);
        end; %if ( ~isfield(dat,stnm) )

      end; %for stix

      %DEBUG:      figure; plot(dat.date,[dat.lkwf1.ww3_wna_sigwavehgt.data,dat.fwyf1.ww3_wna_sigwavehgt.data,dat.mlrf1.ww3_wna_sigwavehgt.data,dat.lonf1.ww3_wna_sigwavehgt.data,dat.smkf1.ww3_wna_sigwavehgt.data, ]); maxigraph; datetick3; legend('LKWF1','FWYF1','MLRF1','LONF1','SMKF1'); title('sigwavehgt');

      disp(['Saving (again) to ' matfname]);
      save(matfname,'dat');  

    end; %for mo

  end; %for yr

  if ( ~isempty(f) )
    close(f);
  end;

  set_more;

return;
