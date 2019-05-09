function [LONS,LATS,U,V,fh,zetar,divr,owp] = anwera(hr,plotflgs,boundingbox,doDownload)
%function [LONS,LATS,U,V,fh,zetar,divr,owp] = anwera(hr,plotflgs,boundingbox,doDownload)
%
% Map/plot WERA data for the Straits of Florida, for the year, jday and
% hour/minute specified by the 'hr' string (e.g., '20083202300').  Optional
% second arg is a logical vector, values as follows: 1: plot vector field;
% 2: curl; 3: divergence; 4: Okubo-Weiss parameter field. DEFAULT [0 0 0 1].
% Optional third arg: [lonW lonE latS latN] boundingbox for all plots. If
% optional fourth arg DODOWNLOAD is FALSE, only local files are checked.
%
% Returned values: plaid arrays (see 'meshgrid') of LONS, LATS, U and V
% current vector components, 'fh' vector of figure handles for all plots,
% 'zetar' plaid array of CURLs, 'divr' DIV array, 'owp' Okubo-Weiss array.
%
% Last Saved Time-stamp: <Thu 2013-04-18 12:40:01 Eastern Daylight Time gramer>

  % Export figures to a location relative to this M-file's local directory
  figspath = get_thesis_path('../figs');

  datapath = get_thesis_path('../data');
  werapath = fullfile(datapath,'wera');

  LONS = [];
  LATS = [];
  U = [];
  V = [];

  fh = nan;
  zetar = [];
  divr = [];
  owp = [];

  if ( ~exist('hr', 'var') )
    error('anwera:NoArg', ...
          'Must specify a year-jday-hour-minute string (first arg)!');
  end;

  if ( ~exist('plotflgs', 'var') || isempty(plotflgs) )
    plotflgs = logical([0 0 0 1]);
  end;
  if ( isscalar(plotflgs) )
    plotflgs = logical([plotflgs plotflgs plotflgs plotflgs]);
  end;
  plotUV = plotflgs(1);
  plotzetar = plotflgs(2);
  plotdivr = plotflgs(3);
  plotowp = plotflgs(4);

  % Default: Include reasonable area surrounding WERA footprint
  if ( ~exist('boundingbox', 'var') || isempty(boundingbox) )
    boundingbox = [ -80.50 -79.10 +24.80 +25.90 ];
  end;

  % Default: Download data files as needed from UM Web site
  if ( ~exist('doDownload', 'var') || isempty(doDownload) )
    doDownload = true;
  end;


  avernm = [hr '.aver.dat'];
  accnm = [hr '.acc.dat'];

  % averpath = sprintf('%s/%s', datapath, avernm);
  % accpath = sprintf('%s/%s', datapath, accnm);
  averpath = fullfile(werapath, avernm);
  accpath = fullfile(werapath, accnm);


  % % OLD BASE URL
  % baseurl = 'http://iwave.rsmas.miami.edu/wera/efs/data';
  % Download data files from RSMAS WERA website if necessary
  baseurl = 'http://iwave.rsmas.miami.edu/wera/efs/data_efsn';

  dy = hr(1:7);
  averurl = sprintf('%s/%s/%s', baseurl, dy, avernm);
  accurl = sprintf('%s/%s/%s', baseurl, dy, accnm);

  if ( ~exist(averpath, 'file') )
    if ( doDownload )
      [respath, result] = urlwrite(averurl, averpath);
      if ( ~exist(respath, 'file') || result ~= 1 )
        warning('anwera:NoURL', ...
                'No AVER file could be downloaded at URL:\n"%s"!', averurl);
        return;
      end;
    else
      warning('anwera:NoLocalAver', ...
              'No local AVER file, and DODOWNLOAD was FALSE:\n"%s"!', averurl);
      return;
    end;
  end;
  if ( ~exist(accpath, 'file') )
    if ( doDownload )
      [respath, result] = urlwrite(accurl, accpath);
      if ( ~exist(respath, 'file') || result ~= 1 )
        warning('anwera:NoURL', ...
                'No ACC file could be downloaded at URL:\n"%s"!', accurl);
        return;
      end;
    else
      warning('anwera:NoLocalAcc', ...
              'No local ACC file, and DODOWNLOAD was FALSE:\n"%s"!', averurl);
      return;
    end;
  end;


  aver = load(averpath);
  acc = load(accpath);

  if ( isempty(aver) || size(aver,1) ~= size(acc,1) )
    warning('anwera:NoData', ...
            'NO DATA\n%s empty or data sizes mismatch!\nCorrupted files?', ...
            averpath);
    return;
  end;

  lon = aver(:,1);
  lat = aver(:,2);
  drn = 90 - aver(:,3);
  spd = aver(:,4);
  n = aver(:,5);

  u = spd .* sind(drn);
  v = spd .* cosd(drn);

  % Use only gridpoints with adequate sample size and accuracy
  goodix = (n >= 4);
  goodix = (acc(goodix, 3) < 10.0);

  if ( numel(acc(goodix)) < 3 )
    warning('anwera:BadData', ...
            'BAD DATA\n%s contains too few points!\nIs file incomplete?', ...
            averpath);
    return;
  end;

  lon = lon(goodix);
  lat = lat(goodix);
  drn = drn(goodix);
  spd = spd(goodix);
  u = u(goodix);
  v = v(goodix);


  LAT = unique(lat);
  LON = unique(lon);
  dlon = min(abs(diff(LON)));
  dlat = min(abs(diff(LAT)));

  if ( numel(LAT) < 2 || numel(LON) < 2 )
    warning('anwera:BadCoords', ...
            'AVER file contained too few coordinates at URL:\n"%s"!', averurl);
    return;
  end;

  % If user gave us a scalar boundingbox - use data limits for this file!
  if ( numel(boundingbox) < 4 )
    disp('Using data limits');
    minlon = min(LON);
    maxlon = max(LON);
    minlat = min(LAT);
    maxlat = max(LAT);
  else
    minlon = boundingbox(1);
    maxlon = boundingbox(2);
    minlat = boundingbox(3);
    maxlat = boundingbox(4);
  end;

  lats = minlat:dlat:maxlat;
  if (abs(lats(end) - maxlat) > eps); lats(end+1) = maxlat; end;
  lons = minlon:dlon:maxlon;
  if (abs(lons(end) - maxlon) > eps); lons(end+1) = maxlon; end;

  [LONS, LATS] = meshgrid(lons, lats);
  SPD = repmat(nan, [length(lats) length(lons)]);
  DRN = repmat(nan, [length(lats) length(lons)]);
  U = repmat(nan, [length(lats) length(lons)]);
  V = repmat(nan, [length(lats) length(lons)]);
  for ix = 1:length(u)
    latix = find(abs(lats - lat(ix)) < dlat);
    lonix = find(abs(lons - lon(ix)) < dlon);
    SPD(latix, lonix) = spd(ix);
    DRN(latix, lonix) = drn(ix);
    U(latix, lonix) = u(ix);
    V(latix, lonix) = v(ix);
  end;


  % Plot actual current vectors, color coded for magnitude
  if ( plotUV )
    fh = figure;
    set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    hold on;
    contourf(lons, lats, SPD, 'LineStyle', 'none');
    quiver(lon, lat, u, v, 0.5);
    map_sofla([LONS(1,1) LONS(end,end) LATS(1,1) LATS(end,end)], [-30 -200 -350]);
    % Idiotic matlab - have to contour AGAIN to get colorbar right
    contour(lons, lats, SPD, 'LineStyle', 'none');
    cbh = colorbar('East');
    % set(get(cbh,'Title'),'string','cm/s')
    set(get(cbh,'YLabel'),'string','cm/s');
    % dcmh = datacursormode(fh);
    % set(dcmh, 'UpdateFcn', @latlonzU_select_cb);
    title(sprintf('%s - WERA 2km surface currents', hr));
    pfname = fullfile(figspath, sprintf('WERA.UV.%s', hr));
    print('-dtiff', '-r300', [pfname '.tiff']);
  end;


  % Curl, divergence, strain, stress, deformation, and Okubo-Weiss parameter
  [dx, dy] = degrees_to_meters(LONS(1,1), LONS(1,2), LATS(1,1), LATS(2,1));
  [zetar, divr, stnr, stsr, dfrm, owp] = okubo_weiss(U, V, dx, dy);

  % % Sanity check - is larger-scale circulation ~ non-divergent?
  % meandivr = mean(divr(~isnan(divr)));
  % fprintf(1, '%s mean divergence: %g\n', hr, meandivr);


  % Plot/Maps

  % Relative vorticity
  if ( plotzetar )
    fh = anwera_plot_prop(hr, LONS, LATS, SPD, U, V, zetar, ...
                          [-3e-3 +3e-3], 0, '\nabla\timesU', '\zeta_R');
    pfname = fullfile(figspath, sprintf('WERA.zetar.%s', hr));
    print('-dtiff', '-r300', [pfname '.tiff']);
  end;

  % Divergence
  if ( plotdivr )
    fh = anwera_plot_prop(hr, LONS, LATS, SPD, U, V, divr, ...
                          [-0.005 +0.010], 0, '\nabla^.U', 'DIV');
    pfname = fullfile(figspath, sprintf('WERA.divr.%s', hr));
    print('-dtiff', '-r300', [pfname '.tiff']);
  end;

  % Okubo-Weiss parameter (deformation-squared minus relvort-squared)
  if ( plotowp )
    [fh, peakixes] = ...
        anwera_plot_prop(hr, LONS, LATS, SPD, U, V, owp, ...
                         [-15e-9 +15e-9], -2e-9, ...
                         'Okubo-Weiss Parameter', 'O-WP');
    pfname = fullfile(figspath, sprintf('WERA.OWP.%s', hr));
    print('-dtiff', '-r300', [pfname '.tiff']);
  end;

return;


function [fh, peakixes] = anwera_plot_prop(hr, LONS, LATS, SPD, U, V, prp, prplim, prppeaks, ttl, ylbl)

  fh = figure;
  set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% COMMENTED OUT FOR SPEED
%   ah(1) = subplot('position', [0.10 0.24 0.88 0.67]);
  ah(1) = gca;
  hold on;

  % Surface-contour the property onto a map of Straits of Florida
  % surf(LONS, LATS, prp);
  contourf(LONS, LATS, prp);
  map_sofla([LONS(1,1) LONS(end,end) LATS(1,1) LATS(end,end)], [-30 -200 -350]);
  view(2);
  cbh = colorbar('East');
  if ( exist('ylbl','var') && ~isempty(ylbl) )
    % set(get(cbh,'Title'),'string',ylbl);
    set(get(cbh,'YLabel'),'string',ylbl);
  end;
  set(gca, 'xlim', [-80.50 -79.10]);
  set(gca, 'ylim', [+24.80 +25.90]);
  if ( ~isempty(prplim) )
    set(gca, 'clim', prplim);
  end;
  colormap(cool);

% COMMENTED OUT FOR SPEED
%   % Quiver-plot maximum current vector along each latitude
%   [Vmax, VmaxIx] = max(SPD, [], 2);
%   % Stupid f'in MATLAB; stupid f'in FORTRAN
%   for ix = 1:length(VmaxIx)
%     qlon(ix) = LONS(ix, VmaxIx(ix));
%     qlat(ix) = LATS(ix, VmaxIx(ix));
%     qu(ix) = U(ix, VmaxIx(ix));
%     qv(ix) = V(ix, VmaxIx(ix));
%   end;
%   % Add a scale vector
%   qlon(end+1) = -79.3;
%   qlat(end+1) = 24.9;
%   qu(end+1) = 100;
%   qv(end+1) = 0;
%   % quiver(qlon, qlat, qu, qv, 0.25);
%   % text(qlon(end), qlat(end)+0.05, '1 m/s');

  % Contour - outline contiguous gridpoints containing extreme values
  if ( ~isempty(prppeaks) )
    peakixes = outline_peaks(fh, LONS, LATS, prp, prppeaks);
  end;

  title(sprintf('%s - %s', hr, ttl));

% COMMENTED OUT FOR SPEED
%   % Boxplot zonal statistics of property, below each longitude
%   ah(2) = subplot('position', [0.10 0.05 0.88 0.16]);
%   boxplot(prp);
%   ylabel(sprintf('%s', ylbl));
%   set(ah(2), 'XTickLabel', []);

return;


function output_txt = latlonzU_select_cb(obj,event_obj)
  pos = get(event_obj,'Position');
  output_txt = {['LON: ',sprintf('%f', pos(1))],...
                ['LAT: ',sprintf('%f', pos(2))]};
  hdl = get(event_obj, 'Target');
  if ( ~isempty(hdl) )
  end;
return;
