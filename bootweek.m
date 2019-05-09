function [m,s] = bootweek(stn,fld,nboot)
%function [m,s] = bootweek(stn,fld,nboot)
%
% Bootstrap estimate the weekly mean and stdev of STN.(FLD).
% Returns 3x52 mean and stdev., each with 95% confidence interval
%
% Last Saved Time-stamp: <Mon 2010-10-04 12:43:21 Eastern Daylight Time gramer>

  set_more off;

  if (~exist('nboot','var') || isempty(nboot) )
    nboot = 300;
  end;

  for wk = 1:52
    %DEBUG:    disp(wk);
    %DEBUG:    tic,
    ix = find(get_week(stn.(fld).date) == wk);
    dat = real(stn.(fld).data(ix));
    m(2,wk) = median(bootstrp(nboot, @nanmean, dat));
    m([1 3],wk) = bootci(nboot, @nanmean, dat);
    if ( nargout > 1 )
      s(2,wk) = median(bootstrp(nboot, @nanstd, dat));
      s([1 3],wk) = s(2,wk);
    end;
    %DEBUG:    toc,
  end;

  set_more;

return;
