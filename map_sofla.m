function [isobath, transects] = map_sofla(boundingbox, isodepth)
%function [isobath, transects] = map_sofla(boundingbox, isodepth)
%
% Draw a map of the South Florida and Bahamas coastlines, together with an
% isobath (at depth 'isodepth', DEFAULT -350m), suitable for plotting WERA
% current vectors or contour fields upon. Returns coordinates of inner
% (Florida near-shore) isobath, and coordinates of periodic 'transect' lines
% perpendicular to that isobath along its whole length.
% DEFAULT boundingbox = [-80.50 -79.10 +24.80 +25.90];
%
% Last Saved Time-stamp: <Tue 2010-09-14 13:02:57  lew.gramer>

  global sofla_coast;
  global sofla_topo;
  global sofla_topo_lats;
  global sofla_topo_lons;

  if ( ~exist('boundingbox', 'var') || isempty(boundingbox) )
    boundingbox = [-80.50 -79.10 +24.80 +25.90];
  end;
  if ( ~exist('isodepth', 'var') || isempty(isodepth) )
    isodepth = -350;
  end;
  if ( isscalar(isodepth) && ~isnan(isodepth) )
    isodepth = [ isodepth isodepth ];
  end;

  isobath = [];
  transects = [];

  hold on;

  % Load and draw medium-resolution coastline

  if ( ~exist('sofla_coast', 'var') || isempty(sofla_coast) )
    disp('Reloading coastline');
    % sofla_coast = load('sofla_coast_low.dat');
    sofla_coast = load('sofla_coast_medium.dat');
  end;

  %line(sofla_coast(:,1), sofla_coast(:,2), ...
  %     'Color', [.4 .3 .2]);
  fill(sofla_coast(:,1), sofla_coast(:,2), [.4 .3 .2]);


  % Load and draw requested isobath from Smith & Sandwell bottom topography

  if ( ~exist('sofla_topo', 'var') || isempty(sofla_topo) )
    % Smith & Sandwell annoyingly already prints something out
    % disp('Reloading topo...');

    [sofla_topo, lats, lons] = ...
        mygrid_sand([boundingbox(3) boundingbox(4) boundingbox(1) boundingbox(2)]);
    % Stupid longitude conventions
    lons(lons > 180) = lons(lons > 180) - 360;

    [sofla_topo_lons, sofla_topo_lats] = meshgrid(lons, lats);
  end;

  if ( numel(sofla_topo) < 4 )
    warning('map_sofla:NotEnoughTopo', ...
            'Too few topo points in bounding box: no isobaths rendered');
    return;
  end;

  if ( all(isnan(isodepth)) )
    disp('Isobars not plotted per user request...');
    return;
  end;


  [cs, ch] = contour(sofla_topo_lons, sofla_topo_lats, sofla_topo, ...
                     isodepth, 'LineColor', [.6 .5 .4]);
  if ( isempty(cs) )
    warning('map_sofla:NoContourLines', ...
            'No isobath contours found in bounding box!');
    return;
  end;

  clabel(cs, ch, 'Color', [.6 .5 .4]);

  set(gca, 'xlim', [boundingbox(1) boundingbox(2)]);
  set(gca, 'ylim', [boundingbox(3) boundingbox(4)]);

  if ( nargout > 0 )
    % Find lat/lon positions (nearest Florida) for desired isobath
    isobath = cs(1:2, 2:(cs(2,1) + 1));

    if ( nargout > 1 )
      % Subset gradient vectors at each interior point of isobath
      %[topx, topy] = gradient(sofla_topo);
      %lnix = ismember(sofla_topo_lons, isobath(1,:));
      %ltix = ismember(sofla_topo_lats, isobath(2,:));
      %topx = topx(lnix & ltix);
      %topy = topy(lnix & ltix);
      %transects(1,:,:) = isobath(:,2:end-1) - [topx ; topy];
      %transects(2,:,:) = isobath(:,2:end-1) + [topx ; topy];

      % Calculate secant and normal slopes at each interior point
      sct = isobath(:,3:end) - isobath(:,1:end-2);
      nrm = [sct(2,:) ; -sct(1,:)];
      % Calculate straight transects orthogonal to our isobath
      transects = repmat(nan, [2 2 length(isobath(1,1:end-2))]);
      transects(1,:,:) = isobath(:,2:end-1) - (2.*nrm);
      transects(2,:,:) = isobath(:,2:end-1) + (2.*nrm);

      %scts = [];
      %line(scts(:,1,:), scts(:,2,:));
    end;

  end;

return;
