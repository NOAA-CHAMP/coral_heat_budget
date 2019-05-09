function app_ii_disp(str,ts,fn)
%function app_ii_disp(str,ts,fn)
% Display median and IQR for time series TS subset by FN
  if ( ~exist('fn','var') || isempty(fn) )
    %DEBUG:    disp('NO FILTER');
    disp({str,nanmedian(ts.data),iqr(ts.data),});
  else
    %DEBUG:    disp(char(fn));
    disp({str,nanmedian(ts.data(fn(ts))),iqr(ts.data(fn(ts))),});
  end;
return;
