function [LON,LAT,VAL] = query_modis_subset(fld,dts,bbox)
%function [LON,LAT,VAL] = query_modis_subset(fld,dts,bbox)
%
% Query USF MODIS dataset for a subset of data field FLD, for the given
% vector of datenums DTS, within the bounding box BBOX. BBOX is in the form
% [lonW lonE latS latN] (DEFAULT: [ -83.00 -79.80 +24.00 +26.00 ]). BBOX may
% also be the 5-char code of an ICON station, e.g., 'SRVI2', Salt River.
% 
% RETURNS: Plaid arrays (v. 'meshgrid') for LONS, LATS, and FLD: 'sst' [oC],
% 'chl' (chlorophyll a) [micro-g/m^2], 'tsm' (total susp. matter) [mg/l]. VAL
% is a matrix NxXxY, X the number of (1km or 4km) longitudinal gridpoints in
% BBOX, Y the number of latitudinal gridpoints, and N the number of synoptic
% fields with *sufficient good data in BBOX* during the given date range.
%
% Last Saved Time-stamp: <Mon 2010-03-29 22:17:11 Eastern Daylight Time gramer>

error('UNDER CONSTRUCTION');

  LONS = [];
  LATS = [];
  VAL = [];

  if ( ~exist('fld', 'var') || ~ischar(fld) )
    error('query_modis_subset:NoArg1', ...
          'First arg must be a variable name string: sst, chl, or tsm!');
  end;
  if ( ~exist('dts', 'var') || ~isnumeric(dts) )
    error('query_modis_subset:NoArg2', ...
          'Second arg must be valid date range, e.g., [datenum(y1,m1,d1):datenum(y2,m2,d2)]!');
  end;

  switch (fld),
   case 'sst',
    dataset = 'seadas_sst';
    historical_dataset = 'seadas_sst';
   case 'chl',
    dataset = 'seadas';
    historical_dataset = 'chl';
   case 'tsm',
    dataset = 'seadas';
    historical_dataset = 'tsm_clark';
   otherwise,
    error('Variable "%s" not yet supported: Should be sst|chl|tsm!',fld);
  end;


  % Default: Include reasonable area surrounding WERA footprint
  if ( ~exist('bbox', 'var') || isempty(bbox) )
    % WERA domain
    %bbox = [ -82.50 -79.10 +24.40 +25.90 ];
    % Wider FRT
    bbox = [ -83.00 -79.80 +24.00 +26.00 ];
  end;

  if ( ischar(bbox) )
    [stnlon,stnlat] = get_station_coords(bbox);
    switch ( lower(bbox(1:4)) ),
     case 'cmrc', region = 'windward';
     case 'srvi', region = 'ias';
     case 'lppr', region = 'ias';
     case 'dbjm', region = 'windward';
     case 'lciy', region = 'ias';
     % otherwise, Default region='florida' is set down below
    end;
    dl = 0.4 + eps;
    bbox = [ (stnlon-dl) (stnlon+dl) (stnlat-dl) (stnlat+dl) ];
  end;

  if ( bbox(1) > bbox(2) || ...
       bbox(3) > bbox(4) || ...
       bbox(3) < -90 || bbox(4) > 90 )
    error('query_modis_subset:BadBounds', ...
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

  lon1 = bbox(1);
  lon2 = bbox(2);
  lat1 = bbox(3);
  lat2 = bbox(4);



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


  jdays = floor(dts);
  ujdays = unique(jdays);

  for jdix = 1:length(ujdays)

    ujd = ujdays(jdix);
    [jyear,ig,ig] = datevec(ujd);
    jday = ujd - datenum(jyear,1,1) + 1;

    % Find full filenames (incl. seconds of timestamp) for sst & chl_a/tsm data
    %http://www.imars.usf.edu/dods-bin/nph-dods/modis/level3/husf/florida/2009/003/1km/pass/final/MODIS.2009003.082513.florida.seadas_sst.hdf
    baseurl = sprintf(['http://www.imars.usf.edu/dods-bin/nph-dods/modis/level3/' ...
                       'husf/%s/%04d/%03d/1km/pass/final'], ...
                      region, jyear, jday);
    intermediate_baseurl = sprintf(['http://www.imars.usf.edu/dods-bin/nph-' ...
                        'dods/modis/level3/husf/%s/%04d/%03d/1km/pass/intermediate'], ...
                                   region, jyear, jday);
    fnames = urlread(baseurl);

    % First be sure we have some FINAL data - if not, use INTERMEDIATE
    patt = sprintf('>MODIS[.]%04d%03d[.][0-9][0-9][0-9][0-9][0-9][0-9]*[.]%s[.][^<]*[.]hdf<', ...
                     jyear, jday, region);
    begix = regexp(fnames, patt);

    % VERY RECENT (no final version yet!)
    if ( isempty(begix) )
      intermediate_fnames = urlread(intermediate_baseurl);
      begix = regexp(intermediate_fnames, patt);
      if ( ~isempty(begix) )
        warning('query_modis_subset:UsingIntermediateData', ...
                'Found no matches for final "%s/%s" data! Using INTERMEDIATE data...', ...
                baseurl, patt);
        baseurl = intermediate_baseurl;
        fnames = intermediate_fnames;
      else
        warning('query_modis_subset:NoIntermediateData', ...
                'Found no matches for FINAL "%s/..." OR INTERMEDIATE "%s/%s" data!', ...
                baseurl, intermediate_baseurl, patt);
        continue;
      end;
    end;


    % Get all files/fields for hours matching those in DTS
    % NOTE: This implies we can get at most one field per hour

    dtix = find(jdays == ujd);
    for ix = 1:length(dtix)

      jd = jdays(dtix(ix));

      [ig,ig,ig,hr] = datevec(jd);
      mnstr='[0-9][0-9]';

      % Extract exactly two matching filenames from USF's directory listing
      % RECENT:
      %  >MODIS.2009003.160747.florida.seadas.hdf<
      %  >MODIS.2009003.160747.florida.seadas_sst.hdf<
      patt = sprintf('>MODIS[.]%04d%03d[.]%02d%s[0-9][0-9]*[.]%s[.](seadas|seadas_sst)[.]hdf<', ...
                     jyear, jday, hr, mnstr, region);
      begix = regexp(fnames, patt);
      historical_data = false;

      % HISTORICAL:
      %  >MODIS.2003135.154510.florida.sst.hdf<
      %  >MODIS.2003135.154510.florida.chl.hdf<
      if ( isempty(begix) )
        patt = sprintf('>MODIS[.]%04d%03d[.]%02d%s[0-9][0-9]*[.]%s[.](chl|sst)[.]hdf<', ...
                       jyear, jday, hr, mnstr, region);
        begix = regexp(fnames, patt);
        historical_data = true;
        if ( isempty(begix) )
          % warning('query_modis_subset:NoHistoricalData', ...
          %         'Found no matches for recent "%s/%s" data! Trying historical data...', ...
          %         baseurl, patt);
          continue;
        end;
      end;

      % This should never happen - logic bug!
      if ( numel(begix) < 2 )
        error('query_modis_subset:TooFewFiles', ...
              'Found less than two matches for "%s/%s"!', ...
              baseurl, patt);
      elseif ( numel(begix) > 2 )
        for ix = begix
          endix = strfind(fnames(ix:end), '<');
          disp(fnames(ix+1:ix+endix(1)-2));
        end;
        error('query_modis_subset:TooManyFiles', ...
              'Found too many matches for "%s/%s"!\n\n%s', ...
              baseurl, patt, ...
              'Verify that HTML catalog file is well-formed...');
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

      %DEBUG:  disp('Querying:'); tic;

      nrows = latix2 - latix1 + 1;
      ncols = lonix2 - lonix1 + 1;
      ixstr = sprintf('[%d:%d][%d:%d]', latix1, latix2, lonix1, lonix2);

      %
      % SST
      %
      if ( strcmp(fld,'sst') )
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
      elseif ( strcmp(fld,'chl') )
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
      % TSM
      %
      elseif ( strcmp(fld,'tsm') )
        if ( historical_data )
          warning('No historical data for TSM: "%s"',
        else
          querychl = sprintf('%s.ascii?chlor_a%s,l2_flags%s,tsm_clark%s', ...
                             chlfname, ixstr, ixstr, ixstr);

      %
      % Grid data
      %

      if ( isempty(LONS) )
        [LONS, LATS] = meshgrid(lons, lats);
      end;

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


    end; %for ix = 1:length(dtix)

  end; %for jdix = 1:length(ujdays)


return;


%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%

