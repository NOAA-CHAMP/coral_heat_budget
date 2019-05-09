function showssts(ssts, urls)
%function showssts(ssts, urls)
% 
% Display a figure with subplots (cv.) for each SST field in 3-D matrix
% (or alternatively cell array of matrices) 'ssts'. Optionally use USF URL
% for each SST field, to build an appropriate title for its SST plot.
% 
% Last Saved Time-stamp: <Thu 2009-07-02 14:35:33 Eastern Daylight Time Lew.Gramer>

  if ( ~exist('urls','var') )
    urls = {};
  end;

  nplts = min(size(ssts,1), 36);
  nrows = min(ceil(sqrt(nplts)), 6);
  ncols = ceil(nplts/nrows);

  minx =    661; maxx =     741;
  miny =   1048; maxy =    1128;
  minval = +0.0; maxval = +40.0;

  figure;
  set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

  minmask = +Inf;
  for ix = 1:nplts
    if ( iscell(ssts) )
      sst = ssts{ix};
    else
      sst = squeeze(ssts(ix,:,:));
    end;
    minmask = min(minmask, length(find(isnan(sst(:)))));
    subplot(nrows, ncols, ix);

    pcolor(flipud(sst));
    shading('flat');
    % What did we plot here - relative or absolute SST??
    % set(gca, 'clim', [-3 +3]);
    % set(gca, 'clim', [22 33]);
    set_pcolor_cursor;

    if ( numel(urls) >= ix && length(urls{ix}) > 0 )
      [ig,sat,yr,mo,dy,hr,mn] = parseusfurl(urls{ix});
      title(sprintf('%s %04d-%02d-%02d %02d:%02d', sat, yr, mo, dy, hr, mn));
    end;
  end;

  disp(minmask);

return;
