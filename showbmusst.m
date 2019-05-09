function fh = showbmusst(stanm, sm, n, dts)
%function fh = showbmusst(stanm, sm, n, dts)
%
% Produce a subplot display showing the synoptic SST images within list of
% datenums 'dts', that correspond to Best Matching Unit(s) 'n' from a Self
% Organizing Map (SOM) analysis. Arg 'n' is an integer scalar or vector.
%
% Last Saved Time-stamp: <Thu 2010-09-16 14:11:02  lew.gramer>

  % Store/retrieve data in this M-file's local directory
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = get_thesis_path('../data');
  end;
  avhrrpath = fullfile(datapath,'avhrr');

  fh = [];

  [dts, sortix] = sort(dts);

  bmus = sm.bmus(sortix);
  synix = find(ismember(bmus, n));

  nimgs = length(synix);

  maxrowcol = 5;

  nplts = min(nimgs, (maxrowcol^2));
  nrows = min(floor(sqrt(nimgs)), maxrowcol);
  ncols = ceil(nplts/nrows);

  minval = +0.0; maxval = +40.0;

  for imgix = 1:(maxrowcol^2):nimgs

   fh(end+1) = figure;
   set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
   for ix = 1:nplts
    if ( (ix + imgix) > nimgs )
      break;
    end;
    myids = synix(ix + (imgix-1));

    subplot(nrows,ncols,ix);

    [yr,mo,dy,hr,mn,sc] = datevec(dts(myids));
    ppatt = fullfile( avhrrpath, ...
                      sprintf('%s-*.%04d%02d%02d.%02d%02d.florida.true.png', ...
                              stanm, yr, mo, dy, hr, round(mn+(sc/60))) );
    dir_st = dir(ppatt);
    if ( length(dir_st) < 1 )
      warning('Skipped! Pattern "%s" not matched.', ppatt);
      continue;
    end;
    fname = dir_st(1).name;
    pname = fullfile(avhrrpath, fname);

    sstbytes = imread(pname);

    % Calculate SST from 1-byte PNG (colormap) values
    sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;
    % Assume that BOTH 254 and 255 are cloud mask values(??)
    sst(sstbytes >= 254) = nan;
    % Just in case our assumptions about cloud mask are wrong!
    sst(minval > sst | sst > maxval) = nan;

    % Correct for diurnal cycle by removing synoptic median
    sst = sst - nanmean(sst(:));

    pcolor(flipud(sst));
    shading('flat');
    set(gca,'clim',[-3 +3]);
    set_pcolor_cursor;
    title(sprintf( '%d %04d%02d%02d %02d%02d', ...
                   bmus(myids), yr, mo, dy, hr, mn ));
   end;

  end;

return;
