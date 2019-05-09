function [sst, LONS, LATS] = get_usf_sst_3d(bbox_or_stanm,yr,wk,dataset,region)
%function [sst,LONS,LATS] = get_usf_sst_3d(bbox_or_stanm,yr,wk,dataset,region)
%
% Download SST mean, anomaly or synoptic image and return SST field. If first
% arg BBOX_OR_STANM is a string, interpret as a SEAKEYS (or ICON) station
% code (e.g., 'fwyf1'), and choose bounding-box accordingly. Otherwise  this
% arg should be a 4-element vector containing [minlon maxlon minlat maxlat].
%
% NOTE: Based on the prior function "ANSST.m" from the same directory.
%
% Last Saved Time-stamp: <Sun 2013-06-09 18:44:06 Eastern Daylight Time gramer>

  % Store/retrieve data in this M-file's local directory
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = get_thesis_path('../data');
  end;
  avhrrpath = fullfile(datapath,'avhrr');

  if ( ~exist('bbox_or_stanm', 'var') || isempty(bbox_or_stanm) )
    % Default is to show the whole Florida Reef Tract
    bbox_or_stanm = [ -84 -79 23 26 ];
  end;
  if ( ~exist('dataset', 'var') || isempty(dataset) )
    dataset = 'mean';
  end;
  if ( ~exist('region', 'var') || isempty(region) )
    region = 'florida';
  end;

  switch ( dataset )
   case 'mean',
    minval = +5.0; maxval = +35.0;
   case 'anomaly',
    minval = -10.0; maxval = +10.0;
   case 'synoptic',
    minval = +0.0; maxval = +40.0;
    error('SYNOPTIC IMAGE ACCESS NOT YET READY!');
   otherwise,
    error('Unrecognized dataset "%s"!', dataset);
  end;

  switch ( region )
   case 'florida',
    minlon = -91; maxlon = -79; dlon = (1/110);
    minlat =  22; maxlat =  31; dlat = (1/110);
   otherwise,
    error('Unrecognized region "%s"!', region);
  end;

  boxradius = 40;

  if ( ~ischar(bbox_or_stanm) )
    bbox = bbox_or_stanm;
    miny = (bbox(1) - minlon) / dlon;
    maxy = (bbox(2) - minlon) / dlon;
    minx = (maxlat - bbox(3)) / dlat;
    maxx = (maxlat - bbox(4)) / dlat;
    stn_x = minx + round((maxx+minx)/2);
    stn_y = miny + round((maxy+miny)/2);

  else
    stanm = bbox_or_stanm;

    [stn_lat, stn_lon] = station_latlon(stanm);

    stn_x = round((maxlat - stn_lat) / dlat);
    stn_y = round((stn_lon - minlon) / dlon);

    minx = stn_x - boxradius; maxx = stn_x + boxradius;
    miny = stn_y - boxradius; maxy = stn_y + boxradius;
  end;


  bminlat = maxlat - (maxx*dlat);
  bmaxlat = maxlat - (minx*dlat);
  bminlon = (miny*dlon) + minlon;
  bmaxlon = (maxy*dlon) + minlon;

  % If user wants to get our lat/lon grid, build it
  LONS = [];
  LATS = [];
  if ( nargout > 1 )
    lons = bminlon:dlon:bmaxlon;
    if (abs(lons(end) - bmaxlon) > eps); lons(end+1) = bmaxlon; end;
    lats = bmaxlat:(-dlat):bminlat;
    if (abs(lats(end) - bminlat) > eps); lats(end+1) = bminlat; end;

    [LONS, LATS] = meshgrid(lons, lats);
  end;


  sstdims = [maxx-minx+1, maxy-miny+1];

  sst = repmat(nan, sstdims);


  if ( ischar(yr) && strcmpi(yr, 'bounds') )
    % User just needs size() of matrix that'd be returned by a "real" call
    % for this region (and may want to get LAT/LON grid at the same time...)
    return;
  end;


  % Sample URL (mean and anomaly datasets):
  % http://www.imars.usf.edu/MERGED_SST/products/anomaly/husf/images/florida/2000.12/weekly/all.2000344.2000350.sst.day_night.florida.mean.png

  % Sample URL (synoptic images):
  % http://www.imars.usf.edu/husf_avhrr/products/images/florida/2009.03/n15.20090303.2211.florida.true.png


  dy2 = (wk * 7);
  dy1 = dy2 - 6;
  dt = datenum(yr, 1, 1) + dy1 - 1;
  [ig, mo, ig] = datevec(dt);

  if ( dy2 >= 364 )
    % If final "week" is >= 7 days, always include up to last day of year
    dy2 = datenum(yr,12,31) - datenum(yr,1,1) + 1;
  end;


  fname = sprintf('all.%04d%03d.%04d%03d.sst.day_night.%s.%s.png', ...
                  yr, dy1, yr, dy2, region, dataset);

  % No peaking into the future!
  if ( dt > now + 7 )
    fprintf(2, '\n');
    warning('Cannot see the future! SKIPPED date %04d.%03d.%04d.%03d', ...
            yr, dy1, yr, dy2);
    return;
  end;


  url = sprintf('http://www.imars.usf.edu/MERGED_SST/products/anomaly/husf/images/%s/%04d.%02d/weekly/%s', ...
                region, yr, mo, fname);

  % Save images for future reference
  fpath = fullfile(avhrrpath, fname);
  if ( exist(fpath, 'file') )
    try
      sstbytes = imread(fpath);
    catch
      % Last ditch effort - delete corrupt file and load from URL anyway
      fprintf(2, '\n');
      warning('Corrupt cached file? Trying URL "%s"...', url);
      delete(fpath);
    end;
  end;
  if ( ~exist(fpath, 'file') )
    [fpath, fstatus] = urlwrite(url, fpath);
    if ( fstatus == 0 )
      fprintf(2, '\n');
      warning('SKIPPED! Failed to download "%s"...', url);
      return;
    end;
    sstbytes = imread(fpath);
  end;

  % Or, if we didn't want to save them
  %sstbytes = imread(url);

  % Calculate SST from 1-byte PNG (colormap) values
  sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;
  % Assume that BOTH 254 and 255 are cloud mask values(??)
  sst(sstbytes >= 254) = nan;
  % Just in case our assumptions about cloud mask are wrong!
  sst(minval > sst | sst > maxval) = nan;


  % Subset to user-specified bounding box - default is Straits of Florida
  sst = sst(minx:maxx, miny:maxy);

return;
