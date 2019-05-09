function app_ii_fit_disp(ts1,ts2,fn1,fn2,ttl1,ttl2)
%function app_ii_fit_disp(ts1,ts2,fn1,fn2,ttl1,ttl2)
% Display bias, slope, RMSE for robust fit between TS1 and TS2 subset by FN1 and FN2

  [B,S]=scatter_fit_ts(ts1,ts2,fn1,fn2,ttl1,ttl2,'none');
  disp(sprintf('%45s \t %0.2g, %0.2g, %0.2g',[ttl1,' vs. ',ttl2],B(1),B(2),S.s));

return;
