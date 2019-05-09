function fld = get_ww3_multi_region(dataset,region,dts,bbox)
%function fld = get_ww3_multi_region(dataset,region,dts,bbox)
%
% Subset a region from NOAA WaveWatch III model output (GRB2 files FTP'd
% under local data directory). Default DTS is 2005-02-01 to 2015-03-31.
%
% Datasets: 'multi_1.at_4m', 'multi_1.wc_4m', 'multi_1.ep_10m', 'multi_1.glo_30m'
% Regions: 'amsam', 'cnmi', 'fknms', 'hi', 'nwhi', 'pr_vi'
%
% For dataset details, see NOAA NCEP, e.g., http://nomad5.ncep.noaa.gov
%
% Last Saved Time-stamp: <Wed 2019-03-06 12:06:51 Eastern Standard Time gramer>

  set_more off;

  fld = [];

  datapath = get_heat_budget_path('../data');
  ww3path = fullfile(datapath, 'ww3');

  if ( ~exist('dataset','var') || isempty(dataset) )
    % See assignment of DEFAULT_DATASET below
    dataset = [];
  end;

  if ( ~exist('region','var') || isempty(region) )
    region = 'pr_vi';
  end;

  switch ( region ),
   case 'amsam',
    minlat = -16.0;
    maxlat = -10.0;
    minlon = -172.0;
    maxlon = -168.0;
    % Multi-mesh model: Equatorial Pacific 10-minute resolution grid
    default_dataset = 'multi_1.ep_10m';
   case 'cnmi',
    minlat = 12.0;
    maxlat = 16.5;
    minlon = 143.0;
    maxlon = 147.0;
    % Multi-mesh model: Equatorial Pacific 10-minute resolution grid
    default_dataset = 'multi_1.ep_10m';
   case 'fknms',
    minlat = 24.0;
    maxlat = 27.0;
    minlon = -84.0;
    maxlon = -79.0;
    % Multi-mesh model: Atlantic 4-minute (lat/lon) resolution grid
    default_dataset = 'multi_1.at_4m';
   case 'hi',
    minlat = 17.5998;
    maxlat = 23.3332;
    minlon = -153.0;
    maxlon = -149.0;
    % Multi-mesh model: West Coast 4-minute resolution grid
    default_dataset = 'multi_1.wc_4m';
   case 'nwhi',
    minlat = 15.0;
    maxlat = 27.0;
    minlon = 194.0;
    maxlon = 211.0;
    % Multi-mesh model: Equatorial Pacific 10-minute resolution grid
    default_dataset = 'multi_1.ep_10m';
   case 'pr_vi',
    minlat = 16.5;
    maxlat = 19.5;
    minlon = -68.0;
    maxlon = -63.0;
    % Multi-mesh model: Atlantic 4-minute (lat/lon) resolution grid
    default_dataset = 'multi_1.at_4m';
   case 'berm',
    minlat = 30.0;
    maxlat = 34.0;
    minlon = -67.0;
    maxlon = -62.0;
    %% Multi-mesh model: Atlantic 4-minute (lat/lon) resolution grid
    %default_dataset = 'multi_1.at_4m';
    % OOPS! The AT_4min is U.S. Shelf only! Gotta use Global 30min
    default_dataset = 'multi_1.glo_30m';
   otherwise,
    error('Unknown region argument "%s"!',region);
  end;

  if ( isempty(dataset) )
    dataset = default_dataset;
  end;
  datafld = strrep(dataset,'multi_1.','');

  if ( ~exist('dts','var') || isempty(dts) )
    %%dts = datenum(2005,2,1):(3/24):datenum(2015,4,1)-(1/24);
    %dts = datenum(2005,2,1):(3/24):datenum(2015,8,1)-(1/24);
    %dts = datenum(2005,2,1):(3/24):datenum(2016,5,1)-(1/24);
    %dts = datenum(2005,2,1):(3/24):datenum(2017,4,1)-(1/24);
    switch ( region ),
     case {'amsam','cnmi'},
      dts = datenum(2005,2,1):(3/24):datenum(2017,6,1)-(1/24);
     otherwise,
      dts = datenum(2005,2,1):(3/24):datenum(2018,7,1)-(1/24);
    end;
  end;
  if ( ~exist('bbox','var') || isempty(bbox) )
    bbox = [minlon,maxlon,minlat,maxlat];
  end;

  savematfname = sprintf('ww3_%s_%s_%06.3f_%06.3f_%06.3f_%06.3f_%s_%s.mat',datafld,region,bbox(:),datestr(dts(1),'yyyymmddHHMM'),datestr(dts(end),'yyyymmddHHMM'));
  savematfname = fullfile(datapath,savematfname);
  %DEBUG:
  disp(savematfname);

  if ( exist(savematfname,'file') )
    disp(['Loading ',savematfname]);
    load(savematfname);

  else
    disp('Extracting data');
    dat=[]; clear dat
    yrmos = unique(get_yearmonth(dts));
    for yrmo = yrmos(:)'
      [mo,yr] = get_month(yrmo);
      matfname = sprintf('ww3.%s.%s.%04d%02d.mat',datafld,region,yr,mo);
      if ( ~exist(fullfile(ww3path,matfname),'file') )
        error('MISSING MONTH %d OF YEAR %d',mo,yr);
      else
        %DEBUG:
        disp(matfname);
        load(fullfile(ww3path,matfname));
        if ( isempty(fld) )
          % If our BBOX touches any gridsquares - encompass it
          lolonix = floor(interp1(dat.lon,1:numel(dat.lon),bbox(1),'linear','extrap'));
          hilonix = ceil(interp1(dat.lon,1:numel(dat.lon),bbox(2),'linear','extrap'));
          lonix = lolonix:hilonix; lonix(1>lonix|lonix>numel(dat.lon))=[];
          lolatix = floor(interp1(dat.lat,1:numel(dat.lat),bbox(4),'linear','extrap'));
          hilatix = ceil(interp1(dat.lat,1:numel(dat.lat),bbox(3),'linear','extrap'));
          latix = lolatix:hilatix; latix(1>latix|latix>numel(dat.lat))=[];
          if ( isempty(lonix) || isempty(latix) )
            error('Bounding box enclosed no data!');
          end;
          fld.lon = dat.lon(lonix);
          fld.lat = dat.lat(latix);
          % For some reason, the very first day of the first file has bad data!
          fld.tp.date = dat.date(2:end);
          fld.tp.field = dat.tp(2:end,latix,lonix);
          fld.dp.date = dat.date(2:end);
          fld.dp.field = dat.dp(2:end,latix,lonix);
          fld.hs.date = dat.date(2:end);
          fld.hs.field = dat.hs(2:end,latix,lonix);
        else
          ndts = numel(dat.date);
          fld.tp.date(end+1:end+ndts) = dat.date;
          fld.tp.field(end+1:end+ndts,:,:) = dat.tp(:,latix,lonix);
          fld.dp.date(end+1:end+ndts) = dat.date;
          fld.dp.field(end+1:end+ndts,:,:) = dat.dp(:,latix,lonix);
          fld.hs.date(end+1:end+ndts) = dat.date;
          fld.hs.field(end+1:end+ndts,:,:) = dat.hs(:,latix,lonix);
        end; %if isempty(fld)
        if 0;
          [u,v]=spddir_to_uv(dat.tp,dat.dp);
          fmg; contourf(dat.lon,dat.lat,squeeze(prctile(dat.hs,93))); cbh=colorbar; ylabel(cbh,'Wave height'); quiver(dat.lon,dat.lat,squeeze(prctile(u,93)),squeeze(prctile(v,93))); titlename(datestr(dat.date(2)));
        end;
        dat=[]; clear dat
      end; %if ~exist(matfname) else
    end; %for yrmo = yrmos(:)'

    disp(['Saving ',savematfname]);
    %save(savematfname,'fld');
    save(savematfname,'fld','-v7.3');

  end; %if ( exist(savematfname,'file') ) else

  if 0;
    dtix = 2;
    [u,v]=spddir_to_uv(fld.tp.field,fld.dp.field);
    fmg; contourf(fld.lon,fld.lat,squeeze(fld.hs.field(dtix,:,:))); set(gca,'CLim',[0,1.5]); cbh=colorbar; ylabel(cbh,'Wave height'); quiver(dat.lon,fld.lat,squeeze(u(dtix,:,:)),squeeze(v(dtix,:,:))); titlename(datestr(fld.hs.date(dtix)));
  end;

  set_more;

return;
