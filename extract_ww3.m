function extract_ww3(dataset, region, yrs)
%function extract_ww3(dataset, region, yrs)
%
% Extract NOAA WaveWatch III wave model output from GRiB files
%
% For details, see NOAA NCEP Web pages, e.g., http://nomad5.ncep.noaa.gov
%
% DEFAULTS: DATASET = 'wna'  (Western North Atlantic)
%           REGION = 'fknms' (Florida Keys NMS, [24 27 -84 -79])
%           YRS = 1999:2007  (All available WW3 years)
%
% Last Saved Time-stamp: <Wed 2010-12-22 15:48:09  lew.gramer>

  set_more off;

  datapath = get_thesis_path('../data');
  ww3path = fullfile(datapath, 'ww3');

  if ( ~exist('dataset','var') || isempty(dataset) )
    dataset = 'wna';
  end;
  if ( ~exist('region','var') || isempty(region) )
    region = 'fknms';
  end;
  if ( ~exist('yrs','var') || isempty(yrs) )
    yrs = 1999:2007;
  end;

  switch ( region ),
   case 'fknms',
    minlat = 24;
    maxlat = 27;
    minlon = 360-84;
    maxlon = 360-79;
   otherwise,
    error('Unknown region argument "%s"!', region);
  end;


  stncfg = {};
  cfgfname = fullfile(datapath,['ww3-' dataset '.cfg']);
  if ( exist(cfgfname,'file') )
    fid = fopen(cfgfname,'r');
    stncfg = textscan(fid,'%s%d%d', 'Delimiter',',', 'CommentStyle','#'); 
    fclose(fid);
  end;

  %ftp://polar.ncep.noaa.gov/pub/history/waves/wna.hs.200711.grb

  % NOTE: If we add any variables here later, be sure to change the IF
  % statement 'if ( ~isfield(dat,stnm) )' below to read 'if (1)'!
  % tp - Peak Wave Period
  % dp - Peak Wave Direction
  % hs - Significant Wave Height
  vars = { 'tp',         'dp',         'hs', };
  flds = { 'peakwaveper','peakwavedir','sigwavehgt', };
  %DEBUG:  vars = { 'hs', };  flds = { 'sigwavehgt', };

  f = [];
  doFTP = false;
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
     %DEBUG:     case 1999, mos = 7;
     %DEBUG:     case 2007, mos = 11;
     case 1999, mos = 7:12;
     case 2007, mos = 1:11;
     otherwise, mos = 1:12;
    end;

    for mo = mos(:)'

      matfname = fullfile(ww3path, ...
                          sprintf('ww3.%s.%04d%02d.mat', dataset, yr, mo));
      if ( exist(matfname,'file') )
        disp(['Loading pre-saved MAT file ' matfname]);
        load(matfname,'dat');

      else
        %DEBUG:
        tic,

        dat = [];
        clear dat;
        dat = [];

        for vix = 1:length(vars)

          var = vars{vix};

          fbasename = sprintf('%s.%s.%04d%02d.grb', dataset, var, yr, mo);
          fname = fullfile(ww3path,fbasename);
          %DEBUG:
          disp(fname);

          if ( ~isempty(f) )
            if ( ~exist(fname,'file') )
              %DEBUG:
              disp(['Attempting to FTP ' fbasename]);
              mget(f, fbasename, ww3path);
            end;
          end;

          if ( ~exist(fname,'file') )
            warning('"%s" not available to process',fname);
            continue;
          end;

          x = read_grib(fname,-1, 'ScreenDiag',0);
          %DEBUG:            size(x),

          if ( ~exist('lats','var') )
            nlon = x(1).gds.Ni;
            nlat = x(1).gds.Nj;

            dlon = x(1).gds.Di;
            %dlat = x(1).gds.Dj;
            % Using dlat=gds.Di because some of these stupid hacker-
            % produced "GRIB" files apparently HAVE NO gds.Dj field!
            dlat = x(1).gds.Di;

            lon1 = x(1).gds.Lo1; lon2 = x(1).gds.Lo2;
            lat1 = x(1).gds.La1; lat2 = x(1).gds.La2;

            lons = lon1 :  dlon : lon2;
            lons(lons > 180) = lons(lons > 180) - 360;
            % Latitude indices are reversed (naturally!)
            lats = lat1 : -dlat : lat2;

            minlon(minlon < 0) = minlon(minlon < 0) + 360;
            maxlon(maxlon < 0) = maxlon(maxlon < 0) + 360;
            minlonix = round( (minlon - lon1) / dlon ) + 1;
            maxlonix = round( (maxlon - lon1) / dlon ) + 1;

            minlatix = round( (lat1 - maxlat) / dlat ) + 1;
            maxlatix = round( (lat1 - minlat) / dlat ) + 1;
          end;

          if ( ~isfield(dat,'n') )
            dat.n = length(x);
            dat.minlonix = minlonix;
            dat.minlatix = minlatix;
            dat.lon = lons(minlonix:maxlonix);
            dat.lat = lats(minlatix:maxlatix);
            dat.nlon = length(dat.lon);
            dat.nlat = length(dat.lat);
            %% Nitwits - 199907 file has stimes of e.g., 11-Jul-2099 21:00:00!
            %for ix = 1:length(x)
            %  dat.date(ix) = datenum(x(ix).stime);
            %end;
            dat.date = [datenum(yr,mo,1) : (3/24) : (datenum(yr,(mo+1),1)-(1/24))]';
          end;

          % Why do these NCEP meteo-wonks insist on adding a record for the
          % first hour of next month into this month's data file??
          for ix = 1:(length(x)-1)
            rawdat = reshape(x(ix).fltarray, [nlon nlat])';
            dat.(var)(ix,1:dat.nlat,1:dat.nlon) = rawdat(minlatix:maxlatix,minlonix:maxlonix);
            rawdat = []; clear rawdat;
          end;
          dat.(var)(dat.(var) > 360) = nan;

          x = [];
          clear x;

          % if ( ~isempty(f) )
          %     delete(fname);
          % end;

        end; %for vix

        disp(['Pre-saving to ' matfname]);
        save(matfname,'dat');  
        %DEBUG:
        toc,
      end; %if exist(matfname) else

      %DEBUG:      id=62; figure; contourf(squeeze(dat.hs(id,:,:))); maxigraph; caxis([0 3]); colorbar; hold on; set(gca,'yd','rev'); title(['sigwavehgt ' datestr(dat.date(id))]);
      %DEBUG:      id=81; figure; contourf(squeeze(dat.hs(id,:,:))); maxigraph; caxis([0 3]); colorbar; hold on; set(gca,'yd','rev'); title(['sigwavehgt ' datestr(dat.date(id))]);

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
      % if ( ~isfield(dat,'minlonix') )
      %   dat.minlonix = 58;
      %   dat.minlatix = 94;
      % end;

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
            fld = ['ww3_' dataset '_' flds{vix}];
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
