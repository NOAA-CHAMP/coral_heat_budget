function showurls(urls)

  nplts = min(length(urls), 36);
  nrows = min(ceil(sqrt(nplts)), 6);
  ncols = ceil(nplts/nrows);

  minx =    661; maxx =     741;
  miny =   1048; maxy =    1128;
  minval = +0.0; maxval = +40.0;

  figure;
  set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

  for ix = 1:nplts
    url = urls{ix};
    subplot(nrows, ncols, ix);

    sstbytes = imread(url);
%     sstbytes = sstbytes(minx:maxx, miny:maxy);
    sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;
    sst(sstbytes >= 254) = nan;
    sst(minval > sst | sst > maxval) = nan;

    pcolor(flipud(sst));
    shading('flat');
    h = datacursormode(gcf);
    set(h, 'UpdateFcn', @xyzc_select_cb);
  end;

return;
