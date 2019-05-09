function [fhs,zetas,divs,owps,maxowps1,maxowps2,dts,wspd,wdir] = animwera(daystrs, plotflgs)
%function [fhs,zetas,divs,owps,maxowps1,maxowps2,dts,wspd,wdir] = animwera(daystrs, plotflgs)
% 
% Plot curl (zetas), divergence (divs), Okubo-Weiss parameter (owps), and
% wind forcing at FWYF1 SEAKEYS Station (wspd, wdir), as well as subset
% fields (maxowps[12]) associated with 2-D WERA surface current fields. Arg
% 'daystrs' is a string array or cell array of strings specifying which
% 24-hour periods of the WERA data should be downloaded (if not present in
% the './data' dir), analyzed and plotted. OPTIONAL arg 'plotflgs' can be a
% logical vector specifying which plots to produce for each hour's worth of
% WERA vectors. (See 'help anwera' for the meaning of each vector element.)
% 
% Sample call (including subsequent call to 'reviewanim' - see help):
%  [f,z,d,o,mo1,mo2,dts,wspd,wdir] = animwera(num2str([2008183:2008185]'));
%  reviewanim(f);
% 
% Last Saved Time-stamp: <Fri 2009-04-17 17:00:45 Eastern Daylight Time gramer>
% 

  more_status = get(0, 'More');
  more('off');

  if ( ~exist('daystrs', 'var') || isempty(daystrs) )
    error('animwera:Usage', 'Missing day-strings argument!');
  end;
  if ( ~exist('plotflgs', 'var') || isempty(plotflgs) )
    % Use ANWERA's default plot flags by default
    plotflgs = [];
  end;

  hrs = 0:23;
  fhs = repmat(nan, [length(hrs) size(daystrs,1)]);
  zetas = [];
  divs = [];
  owps = [];
  maxowps1 = repmat(nan, [length(hrs) size(daystrs,1)]);
  maxowps2 = repmat(nan, [length(hrs) size(daystrs,1)]);
  maxzetas1 = repmat(nan, [length(hrs) size(daystrs,1)]);
  maxzetas2 = repmat(nan, [length(hrs) size(daystrs,1)]);
  dts = repmat(nan, [length(hrs) size(daystrs,1)]);

  for di = 1:size(daystrs, 1)
    daystr = daystrs(di, 1:size(daystrs, 2));
    disp(daystr);
    yr = str2num(daystr(1:4));
    day = str2num(daystr(5:7));
    for hi = 1:length(hrs)
      hr = hrs(hi);
      fname = sprintf('%s%02d00', daystr, hr);
      [LONS, LATS, U, V, fh, z, d, o] = anwera(fname, plotflgs);
      fhs(hi, di) = fh;
      if ( ~isempty(z) && ~isempty(d) && ~isempty(o) )

        % roi1Ix = find(-79.90 <= LONS & LONS <= -79.80 & ...
        %               +25.44 <= LATS & LATS <= +25.55);
        % Encompass the Fowey Rocks FWYF1 SEAKEYS station
        roi1Ix = find(-80.10 <= LONS & LONS <= -79.80 & ...
                      +25.55 <= LATS & LATS <= +25.65);
        maxzetas1(hi, di) = max(max(z(roi1Ix)));
        maxowps1(hi, di) = max(max(o(roi1Ix)));

        roi2Ix = find(-80.02 <= LONS & LONS <= -79.95 & ...
                      +25.45 <= LATS & LATS <= +25.52);
        maxzetas2(hi, di) = max(max(z(roi2Ix)));
        maxowps2(hi, di) = max(max(o(roi2Ix)));

        dts(hi, di) = datenum(yr, 1, 1, hr, 0, 0) + day - 1;

        % if ( isempty(zetas) )
        %   zetas = repmat(nan, [length(hrs) size(daystrs,1) size(z)]);
        %   divs = repmat(nan, [length(hrs) size(daystrs,1) size(z)]);
        %   owps = repmat(nan, [length(hrs) size(daystrs,1) size(z)]);
        % end;
        % zetas(hi, di, :, :) = z;
        % divs(hi, di, :, :) = d;
        % owps(hi, di, :, :) = o;
      end;
    end;
  end;

  % maxzetasd1 = maxzetas1 - mean(maxzetas1(~isnan(maxzetas1(:))));
  % maxzetasd2 = maxzetas2 - mean(maxzetas2(~isnan(maxzetas2(:))));

  fh = figure;
  set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
  % ah(1) = subplot('position', [0.10 0.24 0.88 0.67]);
  ah(1) = subplot('position', [0.10 0.34 0.88 0.67]);
  hold on;
  % plot(dts(:), [maxzetasd1(:) , maxzetasd2(:)]);
  % plot(dts(:), maxzetasd1(:));
  % title('\zeta time series');
  plot(dts(:), maxowps1(:));
  datetick;
  grid minor;
  title('Okubo-Weiss time series');

  % Display local wind forcing beneath time series
  [wdts,wspd] = load_g2_data('FWYF1-WIND1-SPEED.csv', min(dts(:)), max(dts(:)));
  [wdts,wdir] = load_g2_data('FWYF1-WIND1-DIR.csv', min(dts(:)), max(dts(:)));
  if ( numel(wspd) < 10 && numel(wdir) < 10 )
    warning('No WIND data found for this period!');

  else
    wu = wspd .* sind(wind2currdir(wdir));
    wv = wspd .* cosd(wind2currdir(wdir));
    if ( numel(wdts) ~= numel(dts) )
      wu = interp1(wdts, wu, dts, [], nan);
      wv = interp1(wdts, wv, dts, [], nan);
    end;
    [taux, tauy] = wstress(wu, wv, 43);		% FWYF1 anemometer at 43m
    taux = taux .* 0.1; tauy = tauy .* 0.1;	% dyn/cm^2 => N/m^2

    % ah(2) = subplot('position', [0.10 0.05 0.88 0.16]);
    ah(2) = subplot('position', [0.10 0.05 0.88 0.26]);
    taux(1:3:end) = nan; tauy(1:3:end) = nan;
    feather(ah(2), taux, tauy);
    xlim(ah(2), [1 length(wdts)]); ylim([-1 1]);
    set(ah(2), 'XTickLabel', []);
    xlabel('FWYF1 wind stress (\tau: N/m^2)');

  end;

  more(more_status);

return;
