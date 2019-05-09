function stns = get_pathfinder_sst(stns_or_stnms)
%function stns = get_pathfinder_sst(stns_or_stnms)
%
% Use OPeNDAP interface (via netCDF-Java MATLAB Toolbox) to subset Pathfinder
% v5.0 4km AVHRR night-time SST for station(s) identified by structs or name
% strings STNS_OR_STNMS. By DEFAULT download "Best SST" (quality 4 or higher),
% and only for full years 1995 to 2009 incl. If structs STNS are passed in,
% must contain either coordinates in fields .lon and .lat, or recognized 5 char
% station codes in .station_name. Array of structs STNS is returned with
% new once-daily SST time series fields STN(:).pathfinder_sst.
%
% Pathfinder v5.0 product documentation:
%  http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/doc/avhrr_pathfinder_sst.html
%
% Sample data access URLs:
%  http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/night/04km/2000/bsst/2000001.s04d1pfv50-bsst-16b.hdf
%  http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/night/04km/2009/bsst/2009001.s04d1pfrt-bsst.hdf
% As usual, remote sensing wonks change filenames seemingly at random,
% resulting in some of the ugly "contents.html" parsing contained below.
%
% Last Saved Time-stamp: <Thu 2011-03-31 07:08:29  Lew.Gramer>
%
%error('This function downloads a "base SST" which the dataset PIs do not recommend we use!');
error('This function downloads a "base SST" which the dataset PIs do not recommend we use!');

  set_more off;

  if ( isempty(stns_or_stnms) )
    error('First arg was empty!');
  elseif ( ischar(stns_or_stnms) )
    for ix=1:size(stns_or_stnms,1)
      stns(ix).station_name = stns_or_stnms(ix,:);
    end;
  elseif ( iscellstr(stns_or_stnms) )
    for ix=1:numel(stns_or_stnms)
      stns(ix).station_name = stns_or_stnms{ix};
    end;
  elseif ( isstruct(stn_or_stnm) )
    for ix=1:numel(stns_or_stnms)
      stns(ix) = stns_or_stnms(ix);
    end;
  else
    error('First arg must be a station struct(s) or station name string(s) (5-char)!');
  end;
  clear stns_or_stnms;

  if ( ~isfield(stns,'station_name') && isfield(stns,'name') )
    for ix=1:numel(stns)
      stns(ix).station_name = stns(ix).name;
    end;
  end;

  if ( ~isfield(stns,'lon') );    stns(1).lon = [];  end;
  if ( ~isfield(stns,'lat') );    stns(1).lat = [];  end;
  badix = [];
  for ix=1:numel(stns)
    if ( isempty(stns(ix).lon) || isempty(stns(ix).lat) )
      if ( isfield(stns(ix),'station_name') && ~isempty(stns(ix).station_name) )
        try
          [stns(ix).lon,stns(ix).lat,stns(ix).depth] = get_station_coords(stns(ix).station_name);
        catch
        end;
      end;
      if ( isempty(stns(ix).lon) || isempty(stns(ix).lat) )
        badix = [badix ix];
      end;
    end;
  end;
  stns(badix) = [];
  if ( isempty(stns) )
    error('No station struct with coordinates or .station_name field!');
  end;


  yrs = 1995:2009;
  %DEBUG:  yrs = 2003

  % Each year has 73 "pentads", regardless of leap years
  alldts = [];
  for yr=yrs(:)'
    alldts = [ alldts datenum(yr,1,0)+1:5:365 ];
  end;
  alldts = alldts';

  for ix=1:numel(stns)
    if ( isfield(stns(ix),'station_name') && ~isempty(stns(ix).station_name) )
      stations(ix).station_name = stns(ix).station_name;
    end;
    stations(ix).lon = stns(ix).lon;
    stations(ix).lat = stns(ix).lat;
    % Preallocate for speed
    stations(ix).pathfinder_sst.date = repmat(nan,[length(alldts) 1]);
    stations(ix).pathfinder_sst.data = repmat(nan,[length(alldts) 1]);
  end;



  baseurl = 'http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/night/04km';

  missedDays = 0;

  for yr = yrs(:)'

    %DEBUG:
    tic, disp(yr);

    % Regular Expression pattern to extract year-days from HDF filenames
    jdpatt = sprintf('%04d([0-3][0-9][0-9])[.]',yr);

    % Get data filenames for this year (also grabs extra day on leap-years!)
    % Sample URL: http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/night/04km/2005/bsst/contents.html
    url = sprintf('%s/%04d/bsst/contents.html',baseurl,yr);
    s = urlread(url);
    fnames = regexpi(s,'[^a-zA-Z0-9.]([12][90][a-zA-Z0-9.-]*[.]hdf)[^a-zA-Z0-9.]','tokens');

    for fnameix = 1:length(fnames)
      fname = fnames{fnameix}{:};
      %DEBUG:      disp(fname);

      jdstr = regexp(fname,jdpatt,'tokens');
      if ( isempty(jdstr) )
        warning('Found malformed HDF filename "%s"!', fname);
        continue;
      end;
      jd = str2double(jdstr{:});

      %DEBUG:      if (jd>6); break; end;

      dt = datenum(yr,1,1) + jd - 1;
      dtix = find(alldts == dt);
      if ( isempty(dtix) )
        warning('BAD DATE %g ("%s")?!',dt,datestr(dt));
        continue;
      end;


      url = sprintf('%s/%04d/bsst/%s',baseurl,yr,fname);

      %DEBUG:
      disp(url);

      nc = mDataset(url);
      if ( isempty(nc) )
        %DEBUG:
        disp(['BAD DATA: Year ' num2str(yr) ', Day ' num2str(jd)]);
        missedDays = missedDays + 1;
        continue;
      end;

      try
        if ( ~exist('lats','var') || isempty(lats) )
          lats = cast(nc{'lat'}(:),'double');
          lons = cast(nc{'lon'}(:),'double');
          for ix=1:numel(stations)
            [ig,stations(ix).pathfinder_latix]=min(abs(lats - stations(ix).lat));
            [ig,stations(ix).pathfinder_lonix]=min(abs(lons - stations(ix).lon));
          end;
        end;

        % Pixel-to-SST conversion handled by njTbx!
        for ix=1:numel(stations)
          sst = cast(nc{'bsst'}(stations(ix).pathfinder_latix,stations(ix).pathfinder_lonix),'double');
          sst(sst < 1.0) = nan;
          stations(ix).pathfinder_sst.date(dtix,1) = dt;
          stations(ix).pathfinder_sst.data(dtix,1) = sst;
        end;

      catch
        if ( exist('nc','var') && ~isempty(nc) )
          close(nc);
        end;
        clear nc;
        rethrow(lasterror);
      end;
      close(nc); clear nc;

      % Give JPL server a breather between dates!
      pause(1);
    end;

    %DEBUG:
    toc,

  end; %for yr

  if ( missedDays > 0 )
    warning('HDF files for %d days in 1995-2009 not found!',missedDays);
  end;

  set_more;

return;
