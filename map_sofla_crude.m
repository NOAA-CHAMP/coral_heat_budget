function [isobath, transects, sct, nrm] = map_sofla(boundingbox, isodepth)

  global sofla_coast;
  global sofla_topo_low;
  global sofla_topo;

  if ( ~exist('boundingbox', 'var') || isempty(boundingbox) )
    boundingbox = [-80.50 -79.10 +24.80 +25.90];
  end;
  if ( ~exist('isodepth', 'var') || isempty(isodepth) )
    isodepth = -350;
  end;


  hold on;

  % Load and draw medium-resolution coastline

  if ( ~exist('sofla_coast', 'var') || isempty(sofla_coast) )
    disp('Reloading coast...');
    sofla_coast = load('sofla_coast_low.dat');
  end;

  line(sofla_coast(:,1), sofla_coast(:,2), ...
       'Color', [.4 .3 .2]);


  % Load and draw requested isobath from crude bottom topography

  llong = floor(boundingbox(1)+.5)-.5;
  rlong = ceil(boundingbox(2)+.5)-.5;
  blat = max(floor(boundingbox(3)+.5),-89)-.5;
  tlat = min(ceil(boundingbox(4)+.5),90)-.5;
  lgs = (llong:rlong);
  lts = (blat:tlat);
  lgsi = lgs(1):0.02:lgs(end);
  ltsi = lts(1):0.02:lts(end);
  [lgi, lti] = meshgrid(lgsi,ltsi);

  if ( ~exist('sofla_topo', 'var') || isempty(sofla_topo) )
    disp('Reloading topo...');
    load('topo', 'topo');

    if ( rlong<0 )
      topo = topo(lts+90.5, lgs+360.5);
    elseif ( llong < 0 && rlong >= 0 )
      topo = topo(lts+90.5, [(360.5+llong:end) (1:rlong+0.5)]);
    else
      topo = topo(lts+90.5, lgs+.5);
    end;

    sofla_topo_low = topo;
    clear topo;

    [lg,lt] = meshgrid(lgs,lts);
    sofla_topo = interp2(lg, lt, sofla_topo_low, lgi, lti, '*linear');

  end;

  [cs, hs] = contour(lgi, lti, sofla_topo, [isodepth isodepth], ...
                     'LineColor', [.6 .5 .4]);
  % clabel(cs, hs);

  isobath = repmat(nan, [2, length(ltsi)]);
  for ltix = 1:length(ltsi)
    zix = find(abs(sofla_topo(ltix,:) - isodepth) <= 10, 1);
    if ( ~isempty(zix) && lgsi(zix) < -79.5 )
      isobath(1, ltix) = lgsi(zix);
      isobath(2, ltix) = ltsi(ltix);
    end;
  end;


  % Calculate secant and normal slopes at each interior point
  sct = isobath(:,3:end) - isobath(:,1:end-2);
  nrm = -[sct(2,:) ; sct(1,:)];
  transects = repmat(nan, [2 2 length(isobath(1,1:end-2))]);
  transects(1,:,:) = isobath(:,2:end-1) - (2.*nrm);
  transects(2,:,:) = isobath(:,2:end-1) + (2.*nrm);

return;
