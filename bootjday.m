function [m,s] = bootjday(stn,fld,nboot)
%function [m,s] = bootjday(stn,fld,nboot)
%
% Bootstrap estimate the daily mean and stdev of STN.(FLD).
% Returns 3x365 mean and stdev., each with 95% confidence interval
%
% Last Saved Time-stamp: <Fri 2013-03-22 18:43:41 Eastern Daylight Time gramer>

  set_more off;

  if (~exist('nboot','var') || isempty(nboot) )
    nboot = 300;
  end;

  %DEBUG:
  tic,
  for jd = 1:365
    %DEBUG:    disp(jd);
    %DEBUG:    tic,
    ix = find(get_jday_no_leap(stn.(fld).date) == jd);
    dat = real(stn.(fld).data(ix));
    m(2,jd) = median(bootstrp(nboot, @nanmean, dat));
  end;
  %DEBUG:
  toc,
  %DEBUG:
  tic,
  for jd = 1:365
    ix = find(get_jday_no_leap(stn.(fld).date) == jd);
    dat = real(stn.(fld).data(ix));
    m([1 3],jd) = bootci(nboot, @nanmean, dat);
    if ( nargout > 1 )
      s(2,jd) = median(bootstrp(nboot, @nanstd, dat));
      s([1 3],jd) = s(2,jd);
    end;
    %DEBUG:    toc,
  end;
  %DEBUG:
  toc,

  set_more;

return;
