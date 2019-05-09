function [m,s] = bootmon(stn,fld,nboot)
%function [m,s] = bootmon(stn,fld,nboot)
%
% Bootstrap estimate the monthly mean and stdev of STN.(FLD).
% Returns 3x12 mean and stdev., each with 95% confidence interval
%
% Last Saved Time-stamp: <Mon 2010-10-04 14:30:59 Eastern Daylight Time gramer>

  set_more off;

  if (~exist('nboot','var') || isempty(nboot) )
    nboot = 1000;
  end;

  for mo = 1:12
    %DEBUG:    disp(mo);
    %DEBUG:    tic,
    ix = find(get_month(stn.(fld).date) == mo);
    dat = real(stn.(fld).data(ix));
    m(2,mo) = median(bootstrp(nboot, @nanmean, dat));
    m([1 3],mo) = bootci(nboot, @nanmean, dat);
    if ( nargout > 1 )
      s(2,mo) = median(bootstrp(nboot, @nanstd, dat));
      % TOO SLOW! Comment out for now...
      % s([1 3],mo) = bootci(nboot, @nanstd, dat);
      s([1 3],mo) = s(2,mo);
    end;
    %DEBUG:    toc,
  end;

  set_more;

return;
