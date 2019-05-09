function showdailysst(dailyssts,N)

  if ( ~exist('N','var') )
    N = [];
  end;

  nplts = min(length(dailyssts), 42);
  nrows = min(ceil(sqrt(nplts)), 5);
  ncols = ceil(nplts/nrows);

  figure;

  for ix=1:nplts

    subplot(nrows,ncols,ix);

    if ( ~isempty(dailyssts{ix}) )
      pcolor(flipud(dailyssts{ix}));
      shading('flat');
      set(gca,'clim',[20 32]);
      set_pcolor_cursor;
      if ( numel(N) >= ix )
        title(sprintf('N = %d', N{ix}));
      end;
    end;

  end;

return;