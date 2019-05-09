function [LONS, LATS, sst, chl, tsm] = anmodis(dtstr, plotflgs, boundingbox)
%function [LONS, LATS, sst, chl, tsm] = anmodis(dtstr, plotflgs, boundingbox)
%
% Map/plot Aqua MODIS data for the Straits of Florida, for the year, jday and
% hour/minute specified by string DTSTR (e.g., '2008320.2300'). Optional
% second arg PLOTFLGS is a logical vector: 1:plot SST field; 2:plot chl_a;
% 3:plot TSM. Scalar value assumed to apply to all three. DEFAULT: [0 0 0].
% Optional third arg BOUNDINGBOX ([lonW lonE latS latN]) is applied to all
% plots. DEFAULT: boundingbox = [ -83.00 -79.80 +24.00 +26.00 ]. BOUNDINGBOX
% may also be the 5-char code of an ICON station, e.g., 'SRVI2', Salt River.
% 
% Returned values: plaid arrays (see 'meshgrid') of LONS, LATS, sst (1km oC),
% chl_a (1km chlorophyll a, micro-g/m^2), tsm (total suspended matter, mg/l).
%
% Last Saved Time-stamp: <Mon 2010-10-18 11:21:09 Eastern Daylight Time gramer>


  LONS = [];
  LATS = [];
  sst = [];
  chl = [];

  if ( ~exist('dtstr', 'var') )
    error('anmodis:NoArg', ...
          'Must specify a yearjday[.hour[minute]] string (first arg)!');
  end;


  % Parse date string for year, julian day, hour and minute (if any)
  res = sscanf(dtstr, '%04d%03d.%02d%02d');
  jyear = res(1);
  jday = res(2);
  if (length(res) < 3); hr=[]; else; hr=res(3); end;
  if (length(res) < 4); mn=[]; else; mn=res(4); end;
  if (~isempty(mn)); mnstr=sprintf('%02d',mn); else; mnstr='[0-9][0-9]'; end;

  if ( ~exist('plotflgs', 'var') || isempty(plotflgs) )
    plotflgs = logical([0 0 0]);
  end;
  if ( isscalar(plotflgs) )
    plotflgs = logical([plotflgs plotflgs plotflgs]);
  end;
  plotsst = plotflgs(1);
  plotchl = plotflgs(2);
  plottsm = plotflgs(3);

  % Default: Include reasonable area surrounding WERA footprint
  if ( ~exist('boundingbox', 'var') || isempty(boundingbox) )
    % WERA domain
    %boundingbox = [ -82.50 -79.10 +24.40 +25.90 ];
    % Wider FRT
    boundingbox = [ -83.00 -79.80 +24.00 +26.00 ];
  end;

  if ( ischar(boundingbox) )
    [stnlon,stnlat] = get_station_coords(boundingbox);
    switch ( lower(boundingbox(1:4)) ),
     case 'cmrc', region = 'windward';
     case 'srvi', region = 'ias';
     case 'lppr', region = 'ias';
     case 'dbjm', region = 'windward';
     case 'lciy', region = 'ias';
     % otherwise, Default region='florida' is set down below
    end;
    dl = 0.4 + eps;
    boundingbox = [ (stnlon-dl) (stnlon+dl) (stnlat-dl) (stnlat+dl) ];
  end;

  if ( boundingbox(1) > boundingbox(2) || ...
       boundingbox(3) > boundingbox(4) || ...
       boundingbox(3) < -90 || boundingbox(4) > 90 )
    error('anmodis:BadBounds', ...
          'Bounding Box should be [minlon,maxlon,minlat,maxlat]!');
  end;

  if ( ~exist('region', 'var') || isempty(region) )
    region = 'florida';
  end;

  % Region dataset boundaries

  switch ( region )
   case 'flbay',
    % Index ranges per USF website: 0:791 0:1319
    vars = { 'EV_250_RefSB_Band00', 'EV_250_RefSB_Band01' };
    minlat = 24.4;  minlon = -83;
    maxlat = 26.2;  maxlon = -80;
    dlat = 0.00227272727272727;  dlon = 0.00227272727272727;
   case 'florida',
    % Index ranges per USF website: 0:989 0:1319
    minlat = 22;  minlon = -91;
    maxlat = 31;  maxlon = -79;
    dlat = 0.00909090909090909;  dlon = 0.00909090909090909;
   case 'ias',
    % Index ranges per USF website: 0:999 0:1499
    % Map Limit: "8, -98, 33, -58"
    minlat = 08;  minlon = -98;
    maxlat = 33;  maxlon = -58;
    dlat = 0.02500000000000000;  dlon = 0.02666666666666667;
   case 'windward',
    % Index ranges per USF website: 0:989 0:989
    % Map Limit: "16, -80, 25, -71"
    minlat = 16;  minlon = -80;
    maxlat = 25;  maxlon = -71;
    dlat = 0.00909090909090909;  dlon = 0.00909090909090909;
  end;


  % The subset (or intersection) of interest to us

  lon1 = boundingbox(1);
  lon2 = boundingbox(2);
  lat1 = boundingbox(3);
  lat2 = boundingbox(4);



  % NOTE: LOWER RIGHT corner of image would be:
  %  latix1 = 689;  latix2 = 989;  lonix1 = 1019;  lonix2 = 1319;
  lonix1 = floor((lon1 - minlon) / dlon);
  lonix2 = ceil((lon2 - minlon) / dlon);
  latix1 = floor((maxlat - lat2 + dlat) / dlat);
  latix2 = ceil((maxlat - lat1 + dlat) / dlat);


  % Shift user coords to be at the appropriate gridpoints
  lon1 = minlon + (lonix1 * dlon);
  lon2 = minlon + (lonix2 * dlon);
  lat1 = maxlat - (latix2 * dlat) - dlat;
  lat2 = maxlat - (latix1 * dlat) - dlat;

  lons = lon1:dlon:lon2;
  if (abs(lons(end) - lon2) > eps); lons(end+1) = lon2; end;
  lats = lat1:dlat:lat2;
  if (abs(lats(end) - lat2) > eps); lats(end+1) = lat2; end;


  % Find full filenames (incl. seconds of timestamp) for sst & chl_a/tsm data

  %http://www.imars.usf.edu/dods-bin/nph-dods/modis/level3/husf/florida/2009/003/1km/pass/final/MODIS.2009003.082513.florida.seadas_sst.hdf

  baseurl = sprintf(['http://www.imars.usf.edu/dods-bin/nph-dods/modis/level3/' ...
                     'husf/%s/%04d/%03d/1km/pass/final'], ...
                    region, jyear, jday);
  intermediate_baseurl = sprintf(['http://www.imars.usf.edu/dods-bin/nph-' ...
                      'dods/modis/level3/husf/%s/%04d/%03d/1km/pass/intermediate'], ...
                                 region, jyear, jday);

  baseurl = ...
      sprintf('http://cyclops.marine.usf.edu/modis/level3/husf/%s/%04d/%03d/1km/pass/final',...
              region,jyear,jday);
  intermediate_baseurl = ...
      sprintf('http://cyclops.marine.usf.edu/modis/level3/husf/%s/%04d/%03d/1km/pass/intermediate',...
              region,jyear,jday);

  fnames = urlread(baseurl);

  % Extract exactly two matching filenames from USF's directory listing
  % RECENT:
  %  >MODIS.2009003.160747.florida.seadas.hdf<
  %  >MODIS.2009003.160747.florida.seadas_sst.hdf<
  patt = sprintf('>MODIS[.]%04d%03d[.]%02d%s[0-9][0-9]*[.]%s[.](seadas|seadas_sst)[.]hdf<', ...
                 jyear, jday, hr, mnstr, region);
  begix = regexp(fnames, patt);
  historical_data = false;

  % VERY RECENT (no final version yet!)
  if ( isempty(begix) )
    warning('anmodis:TryingIntermediateData', ...
            'Found no matches for final "%s/%s" data! Trying INTERMEDIATE data...', ...
            baseurl, patt);
    intermediate_fnames = urlread(intermediate_baseurl);
    begix = regexp(intermediate_fnames, patt);
    if ( ~isempty(begix) )
      baseurl = intermediate_baseurl;
      fnames = intermediate_fnames;
    end;
  end;

  % HISTORICAL:
  %  >MODIS.2003135.154510.florida.sst.hdf<
  %  >MODIS.2003135.154510.florida.chl.hdf<
  if ( isempty(begix) )
    warning('anmodis:TryingHistoricalData', ...
            'Found no matches for recent "%s/%s" data! Trying historical data...', ...
            baseurl, patt);
    patt = sprintf('>MODIS[.]%04d%03d[.]%02d%s[0-9][0-9]*[.]%s[.](chl|sst)[.]hdf<', ...
                   jyear, jday, hr, mnstr, region);
    begix = regexp(fnames, patt);
    historical_data = true;
  end;

  if ( numel(begix) < 2 )
    error('anmodis:TooFewFiles', ...
          'Found less than two matches for "%s/%s"!', ...
          baseurl, patt);
  elseif ( numel(begix) > 2 )
    for ix = begix
      endix = strfind(fnames(ix:end), '<');
      disp(fnames(ix+1:ix+endix(1)-2));
    end;
    error('anmodis:TooManyFiles', ...
          'Found too many matches for "%s/%s"!\n\n%s', ...
          baseurl, patt, ...
          'Try making your time-string more specific...');
  end;

  % Figure out which filename is which
  for ix = begix
    endix = strfind(fnames(ix:end), '<');
    fname = fnames(ix+1:ix+endix(1)-2);
    if ( strfind(fname, 'sst') )
      sstfname = fname;
    else
      chlfname = fname;
    end;
  end;



  % Build Web data query string from filenames and our bounding box

  %http://www.imars.usf.edu/dods-bin/nph-dods/modis/level3/husf/florida/2009/003/1km/pass/final/MODIS.2009003.160747.florida.seadas_sst.hdf.ascii?seadas_sst[330:332][1167:1169],l2_flags[330:332][1167:1169],cloud_mask[330:332][1167:1169]
  %http://www.imars.usf.edu/dods-bin/nph-dods/modis/level3/husf/florida/2008/222/1km/pass/final/MODIS.2008222.153656.florida.seadas.hdf.ascii?chlor_a[330:332][1167:1169],tsm_clark[330:332][1167:1169]

  % DEBUG:  disp('Querying:'); tic;

  nrows = latix2 - latix1 + 1;
  ncols = lonix2 - lonix1 + 1;
  ixstr = sprintf('[%d:%d][%d:%d]', latix1, latix2, lonix1, lonix2);


  %
  % SST
  %

  if ( historical_data )
    querysst = sprintf('%s.ascii?sst%s', ...
                       sstfname, ixstr);
  else
    querysst = sprintf('%s.ascii?seadas_sst%s,l2_flags%s,cloud_mask%s', ...
                       sstfname, ixstr, ixstr, ixstr);
  end;
  % DEBUG:  disp('SST:'); toc;
  sst_s = query_dods(baseurl, querysst, nrows, ncols);
  if ( isempty(sst_s) )
    querysst = sprintf('%s.ascii?seadas_sst%s,cloud_mask%s', ...
                       sstfname, ixstr, ixstr);
    sst_s = query_dods(baseurl, querysst, nrows, ncols);
    if ( isempty(sst_s) )
      querysst = sprintf('%s.ascii?seadas_sst%s', ...
                         sstfname, ixstr);
      sst_s = query_dods(baseurl, querysst, nrows, ncols);
    end;
  end;
  % DEBUG:  disp('Done'); toc;


  %
  % CHL
  %

  if ( historical_data )
    querychl = sprintf('%s.ascii?chlor_a_2%s', ...
                       chlfname, ixstr);
  else
    querychl = sprintf('%s.ascii?chlor_a%s,l2_flags%s,tsm_clark%s', ...
                       chlfname, ixstr, ixstr, ixstr);
  end;
  % DEBUG:  disp('CHL/TSM:'); toc;
  chl_s = query_dods(baseurl, querychl, nrows, ncols);
  if ( isempty(chl_s) )
    querychl = sprintf('%s.ascii?chlor_a%s', ...
                       chlfname, ixstr);
    chl_s = query_dods(baseurl, querychl, nrows, ncols);
  end;
  % DEBUG:  disp('Done'); toc;


  %
  % Grid data
  %

  [LONS, LATS] = meshgrid(lons, lats);

  sst = repmat(nan, [length(lats) length(lons)]);
  chl = repmat(nan, [length(lats) length(lons)]);
  tsm = repmat(nan, [length(lats) length(lons)]);


  % Do any preprocessing required for each raw dataset

  % With the passage of time, USF DODS has served various datasets
  sstdat = [];
  if ( isempty(sst_s) )
    warning('No SST data found!');
  elseif ( isfield(sst_s, 'seadas_sst') )
    sstdat = sst_s.seadas_sst;
    sstinter = 0;
    sstslope = 0.005;
  else
    % Clear two "quality bits" at positions 1 and 2
    sstbit = cast(sst_s.sst, 'uint16');
    sstbit = bitset(sstbit, 1, 0);
    sstbit = bitset(sstbit, 2, 0);
    sstdat = cast(sstbit, 'double');
    sstinter = -300;
    sstslope = 0.01;
  end;

  chldat = [];
  if ( isempty(chl_s) )
    warning('No CHLOROPHYLL or TSM data found!');
  elseif ( isfield(chl_s, 'chlor_a') )
    chldat = chl_s.chlor_a;
    chlinter = 0;
    chlslope = 1.0;
  else
    % Clear two "quality bits" at positions 7 and 8
    chlbit = cast(chl_s.chlor_a_2, 'uint16');
    chlbit = bitset(chlbit, 7, 0);
    chlbit = bitset(chlbit, 8, 0);
    chldat = cast(chlbit, 'double');
    chlinter = -5;
    chlslope = 0.002;
  end;


  % Process flags and masks - replace "bad" values with NaN
  % (See: http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html)

  % long_name: "Level-2 Processing Flags\\000"
  % f01_name: "ATMFAIL\\000"	  % f02_name: "LAND\\000"
  % f03_name: "BADANC\\000"	  % f04_name: "HIGLINT\\000"	
  % f05_name: "HILT\\000"	  % f06_name: "HISATZEN\\000"	
  % f07_name: "COASTZ\\000"	  % f08_name: "NEGLW\\000"	
  % f09_name: "STRAYLIGHT\\000"	  % f10_name: "CLDICE\\000"	
  % f11_name: "COCCOLITH\\000"	  % f12_name: "TURBIDW\\000"
  % f13_name: "HISOLZEN\\000"	  % f14_name: "HITAU\\000"
  % f15_name: "LOWLW\\000"	  % f16_name: "CHLFAIL\\000"
  % f17_name: "NAVWARN\\000"	  % f18_name: "ABSAER\\000"
  % f19_name: "TRICHO\\000"	  % f20_name: "MAXAERITER\\000"
  % f21_name: "MODGLINT\\000"	  % f22_name: "CHLWARN\\000"
  % f23_name: "ATMWARN\\000"	  % f24_name: "DARKPIXEL\\000"
  % f25_name: "SEAICE\\000"	  % f26_name: "NAVFAIL\\000"
  % f27_name: "FILTER\\000"	  % f28_name: "SSTWARN\\000"
  % f29_name: "SSTFAIL\\000"	  % f30_name: "HIPOL\\000"
  % f31_name: "SPARE\\000"	  % f32_name: "OCEAN\\000"


  % Take whatever flags and masks we can get for this era
  flgs = repmat(1, size(sstdat));
  landmask = repmat(0, size(sstdat));

  if ( isfield(sst_s, 'cloud_mask') )
    flgs = flgs & (sst_s.cloud_mask == 0);
  end;
  if ( isfield(sst_s, 'l2_flags') )
    landmask = landmask | bitget(sst_s.l2_flags,1);
    flgs = flgs & (~landmask) & (~bitget(sst_s.l2_flags,2));
    % SSTFAIL:
    %        (~bitget(sst_s.l2_flags,29));
  end;
  % NOTE WELL: chl borrows cloud mask from SST data!
  if ( isfield(chl_s, 'l2_flags') )
    landmask = landmask | bitget(chl_s.l2_flags,1);
    flgs = flgs & (~landmask) & (~bitget(chl_s.l2_flags,2));
    % CHLFAIL:
    %        (~bitget(chl_s.l2_flags,16));
    % TSM-uncontaminated CHL values:
    %        (chl_s.tsm_clark < 10);
  end;


  % NOTE: 'flgs' is an OR of the sst and chl flags, if any

  if ( isempty(flgs) )
    error('Found no data, or else flags were all bad!');
  end;

  sst(flgs) = (sstinter + (sstdat(flgs) .* sstslope));
  sst(5 > sst | sst > 40) = nan;

  chl(flgs) = (chlinter + (chldat(flgs) .* chlslope));
  chl(0 > chl | chl > 10) = nan;

  if ( isfield(chl_s, 'tsm_clark') )
    tsm(flgs) = chl_s.tsm_clark(flgs);
    tsm(0 > tsm | tsm > 100) = nan;
  end;


  %
  % Plot/Maps
  %

  % Sea Surface Temperature
  if ( plotsst )
    sstclim = [max(20.0, min(sst(:))) min(32.0, max(sst(:)))];
    sstlim = mean(sst(:)) - std(sst(:));
    fh = anmodis_plot_prop(sstfname, LONS, LATS, sst, sstclim, sstlim, ...
                           'Sea Surface Temperature (^oC)', 'SST', 'default');
  end;

  % Chlorophyll /a/
  if ( plotchl )
    %[0.0 3.0], 1.5,
    %[0.0 0.6], 0.4,
    fh = anmodis_plot_prop(chlfname, LONS, LATS, chl, [0.0 3.0], 1.5, ...
                           'Chlorophyll \it{a} (mg^.m^-^3)', 'Chl', 'winter');
  end;

  % Total Suspended Matter
  if ( plottsm )
    fh = anmodis_plot_prop(chlfname, LONS, LATS, tsm, [0.0 3.0], 1.0, ...
                           'Total Suspended Matter (mg^.l^-^1)', 'TSM', 'summer');
  end;


return;


%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%

function [fh, peakixes] = anmodis_plot_prop(fn,LONS,LATS,prp,prplim,prppeaks,ttl,ylbl,cmap)

  if (~exist('cmap','var')||isempty(cmap)); cmap='default'; end;

  % Store/retrieve data in this M-file's local directory
  if ( ~exist('figspath', 'var') || isempty(figspath) )
    figspath = get_thesis_path('../figs');
  end;

  fh = [];
  peakixes = [];
  if ( all(isnan(prp)) )
    warning('No data for "%s"! Skipping plot...', ttl);
    return;
  end;


  fh = figure;
  set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% COMMENTED OUT FOR SPEED
%   ah(1) = subplot('position', [0.10 0.24 0.88 0.67]);
  ah(1) = subplot('position', [0.05 0.05 0.90 0.90]);
  hold on;

  % Assign colormap
  colormap(ah(1), cmap);

  % Surface-contour the property onto a map of Straits of Florida
  pcolor(LONS, LATS, prp);
  % CAREFUL: Non-interp maps look awful - but this IS a smoothing operation!
  shading('interp');
  set_pcolor_cursor;

  % Show (rough) h=-30 and h=-220 isobaths of the bottom topography
  map_sofla([LONS(1,1) LONS(end,end) LATS(1,1) LATS(end,end)], [-30 -220]);
  % Plot locations of all SEAKEYS and ICON stations
  plot_seakeys_icon;
  view(2);

  cbh = colorbar('EastOutside');

  set(ah(1), 'xlim', [LONS(1,1) LONS(end,end)]);
  set(ah(1), 'ylim', [LATS(1,1) LATS(end,end)]);
  if ( ~isempty(prplim) )
    set(ah(1), 'clim', prplim);
  else
    set(ah(1), 'clim', 'default');
  end;


  % Outline contiguous gridpoints containing values of interest
  if ( ~isempty(prppeaks) )
    peakixes = outline_peaks(fh, LONS, LATS, prp, prppeaks);
  end;

  ttlstr = sprintf('%s - %s', fn, ttl);
  th = title(strrep(ttlstr, '_', '\_'));
  set(fh, 'Name', ttlstr);
  figfname = fullfile(figspath, sprintf('%s.%s.png', fn, ylbl));
%   print('-dpng', figfname);

% COMMENTED OUT FOR SPEED
%   % Boxplot zonal statistics of property, below each longitude
%   ah(2) = subplot('position', [0.10 0.05 0.88 0.16]);
%   boxplot(prp);
%   ylabel(sprintf('%s', ylbl));
%   set(ah(2), 'XTickLabel', []);

return;
