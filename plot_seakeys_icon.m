function STATIONS = plot_seakeys_icon(ah,stnms)
%function STATIONS = plot_seakeys_icon(ah,stnms)
%
% Plot positions and labels of any SEAKEYS and ICON stations on the existing
% map-plot in axes 'ah': station locations outside the current lat/lon limits
% of 'ah' are not plotted, and the xlim and ylim of 'ah' are not changed. If
% 'ah' is not specified, the current axes ('gca') is used. NOTE: Temporarily
% changes the 'hold' status to 'on' ('NextPlot' axes property to 'add') in
% 'ah'. After plotting, previous hold status should be restored.
%
% If optional 'stnms' cell array of strings is given, only plot stations
% whose names (5-char codes) appear in that cellstr.
%
% Last Saved Time-stamp: <Tue 2009-07-14 09:50:32 Eastern Daylight Time Lew.Gramer>

  if ( ~exist('ah', 'var') || isempty(ah) )
    ah = gca;
  end;

  STATIONS = get_seakeys_icon;

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
  inix = find(minlon <= STATIONS.lons & STATIONS.lons <= maxlon & ...
              minlat <= STATIONS.lats & STATIONS.lats <= maxlat);
  if ( nargin > 1 && ~isempty(stnms) )
    inix = intersect(inix, find(ismember(STATIONS.codes, stnms)));
  end;

  plot(STATIONS.lons(inix), STATIONS.lats(inix), 'pr');
  text(STATIONS.lons(inix), STATIONS.lats(inix), ...
       strcat(' \leftarrow ', STATIONS.codes(inix)), ...
       'Color', 'r', 'VerticalAlignment', 'middle');

  % Restore original plot axes size after we are done
  xlim(ah, xl); ylim(ah, yl);
  set(ah, 'NextPlot', hold_status);

return;
