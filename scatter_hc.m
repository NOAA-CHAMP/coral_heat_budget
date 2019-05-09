function [t,f] = scatter_hc(stn,tfld,ffld,tfix,ffix)
%function [t,f] = scatter_hc(stn,tfld,ffld,tfix,ffix)
%
% Linearly regress the rate of change in STN.(TFLD) (determined using an
% N-point cenetered finite difference), against the values in STN.(FFLD).
% Useful for comparing sea temperature change to heat flux, for example!
%
% Last Saved Time-stamp: <Thu 2010-07-22 14:37:30  Lew.Gramer>

  stn = verify_variable(stn,tfld);
  stn = verify_variable(stn,ffld);


  f.date = stn.(ffld).date;
  f.data = real(stn.(ffld).data);

  % Centered N-point finite difference
  t = findiff_ts(stn.(tfld),5);

  if ( ~exist('tfix','var') || isempty(tfix) )
    tfix = 1:length(t.data);
  end;
  if ( ~exist('ffix','var') || isempty(ffix) )
    ffix = 1:length(f.data);
  end;

  tlbl = strrep(tfld,'_','\_');
  xlbl = ['d(' tlbl ')/dt [K/hr]'];
  ylbl = [strrep(ffld,'_','\_') ' [K/hr]'];
  [B,Stats,fh,fix1,fix2] = scatter_fit_ts(t,f,tfix,ffix,xlbl,ylbl,[],true);

  [ix1,ix2]=intersect_dates(t.date,f.date);
  x = t.data(ix1);
  y = f.data(ix2);
  Bs = bootstrp(100,@robustfit,x,y);
  B = mean(Bs);
  Bmin = min(Bs);
  Bmax = max(Bs);

  figure;
  % plot(x,y,'.', 'MarkerSize',3, 'Color',[.4 .4 .4]);
  plot(x,y,'.', 'MarkerSize',6, 'Color',[.8 .8 .8]);
  plot(x,(B(1)+(B(2).*x)),'-', 'Color',[.0 .0 .0]);
  plot(x,(Bmin(1)+(Bmin(2).*x)),'-', 'Color',[.0 .0 .0]);
  plot(x,(Bmax(1)+(Bmax(2).*x)),'-', 'Color',[.0 .0 .0]);
  plot(x,x,'--', 'Color',[.5 .5 .5]);
  maxigraph;

  multiplot_datetick({stn.(tfld).date,t.date,f.date,t.date(fix1)}, ...
                     {stn.(tfld).data,t.data,f.data,Stats.resid}, ...
                     'Time series and residuals',[],{tlbl,xlbl,ylbl,'resid'});

return;
