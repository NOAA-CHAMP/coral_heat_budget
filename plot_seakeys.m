function SEAKEYS = plot_seakeys(ah)
%function SEAKEYS = plot_seakeys(ah)
%
% Plot positions and labels of SEAKEYS stations on the existing map-plot in
% axes 'ah': station locations outside the current lat/lon limits of 'ah' are
% not plotted, and the xlim and ylim of 'ah' are not changed. If 'ah' is not
% specified, the current axes ('gca') is used. NOTE: Temporarily changes the
% 'hold' status to 'on' ('NextPlot' axes property to 'add') in 'ah'. After
% plotting, previous hold status should be restored.
%
% Last Saved Time-stamp: <Sat 2009-07-04 16:14:54 Eastern Daylight Time gramer>

  if ( ~exist('ah', 'var') || isempty(ah) )
    ah = gca;
  end;

  SEAKEYS = get_seakeys;

  axes(ah);
  hold_status = get(ah, 'NextPlot');
  hold(ah, 'on');
  xl = xlim(ah); yl = ylim(ah);

  % Complicated logic to find all stations "inside" our current plot axes
  tol = 0.1;
  dlon = abs(xl(2) - xl(1));
  dlat = abs(yl(2) - yl(1));
  minlon = xl(1) - (tol * dlon); maxlon = xl(2) + (tol * dlon);
  minlat = yl(1) - (tol * dlat); maxlat = yl(2) + (tol * dlat);
  inix = find(minlon <= SEAKEYS.lons & SEAKEYS.lons <= maxlon & ...
              minlat <= SEAKEYS.lats & SEAKEYS.lats <= maxlat);

  plot(SEAKEYS.lons(inix), SEAKEYS.lats(inix), 'pr');
  text(SEAKEYS.lons(inix), SEAKEYS.lats(inix), ...
       strcat(' \leftarrow ', SEAKEYS.codes(inix)), ...
       'Color', 'r', 'VerticalAlignment', 'middle');

  % Restore original plot axes size after we are done
  xlim(ah, xl); ylim(ah, yl);
  set(ah, 'NextPlot', hold_status);

return;
