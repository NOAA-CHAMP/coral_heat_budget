function peakixes = outline_peaks(fh, LONS, LATS, prp, cutoffs)
%function peakixes = outline_peaks(fh, LONS, LATS, prp, cutoffs)
%
% Use MATLAB 'contour' to outline peak values in property plot 'fh'.
%
% Last Saved Time-stamp: <Tue 2009-02-03 18:44:15 Eastern Standard Time gramer>

  peakixes = {};

  if ( numel(cutoffs) == 1 )
    %%%% FORGOT WHY THIS abs() WAS IN HERE??
    % prp = abs(prp);
    cutoffs = [cutoffs cutoffs];
  end;

  figure(fh);
  [cs, ch] = contour(LONS, LATS, prp, cutoffs, 'r');

  % Parse out the matrix of contour vertices (see HELP CONTOURC)
  ix = 1;
  curcon = 1;
  while ( curcon <= size(cs, 2) )
    nverts = cs(2, curcon);

    % Points and lines are not polygons... Stupid MATLAB.
    if ( nverts >= 3 )
      verts = cs(:, curcon+1:curcon+nverts);

      % Sometimes CONTOURC annoyingly repeats vertices multiple times
      uniqix = find((diff(verts(1,:)) ~= 0) | (diff(verts(2,:)) ~= 0));
      % End vertex is always dupliate
      verts = verts(:, uniqix);

      insideidx = inside(LONS, LATS, verts(1,:), verts(2,:));
      peakixes{ix} = find(insideidx > 0);

      ix = ix + 1;
    end;

    curcon = curcon + nverts + 1;
  end;

return;
