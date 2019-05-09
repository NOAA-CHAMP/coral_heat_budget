function [m,s] = bootyearday(stn,fld,nboot)
%function [m,s] = bootyearday(stn,fld,nboot)
%
% Bootstrap estimate the hourly mean and stdev of STN.(FLD).
% Returns 3x8760 mean and stdev., each with 95% confidence interval
%
% Last Saved Time-stamp: <Fri 2013-03-22 23:24:43 Eastern Daylight Time gramer>

  set_more off;

  if (~exist('nboot','var') || isempty(nboot) )
    nboot = 30;
  end;

  %DEBUG:
  tic,
  for hrix = 1:8760
    %DEBUG:    disp(hr);
    %DEBUG:    tic,
    ix = find(round(get_yearday(stn.(fld).date)*24)+1 == hrix);
    dat = real(stn.(fld).data(ix));
    m(2,hrix) = median(bootstrp(nboot, @nanmean, dat));
  end;
  %DEBUG:
  toc,
  %DEBUG:
  tic,
  for hrix = 1:8760
    ix = find(round(get_yearday(stn.(fld).date)*24)+1 == hrix);
    dat = real(stn.(fld).data(ix));
    m([1 3],hrix) = bootci(nboot, @nanmean, dat);
    if ( nargout > 1 )
      s(2,hrix) = median(bootstrp(nboot, @nanstd, dat));
      s([1 3],hrix) = s(2,hrix);
    end;
    %DEBUG:    toc,
  end;
  %DEBUG:
  toc,

  set_more;

return;
