function [dts,x] = mocum(cs)
%function mocum(cs)
%
% Accumulate monthly mean values CS.momean over multiple years: return
% datenum vector DTS and cumulative vector X; plot results as PLOT(DTS,X).
% The input CS may be the result of calling, e.g., ANCLIM. (v.)

  y = [cs.momean(:,:)]'; y = y(:)';
  d = 1:length(y);
  z = interp1(d(~isnan(y)),y(~isnan(y)),d);
  z(isnan(z)) = 0;
  x = cumsum(z);
  x(isnan(y)) = nan;

  dts = datenum(cs.yrs(1),1,1):30:datenum(cs.yrs(end),12,31);
  dts = dts(1:length(x));

  figure;
  plot(dts,x);
  maxigraph;
  datetick3;  

return;
