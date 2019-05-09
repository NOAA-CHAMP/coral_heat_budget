function [t,f] = scatter_hc(stn,tfld,ffld,tfix,ffix)
%function [t,f] = scatter_hc(stn,tfld,ffld,tfix,ffix)

  stn = verify_variable(stn,tfld);
  stn = verify_variable(stn,ffld);


  f.date = stn.(ffld).date;
  f.data = real(stn.(ffld).data);

  t.date = stn.(tfld).date(1:end-1);
  t.data = diff(stn.(tfld).data);

  badix = find(abs(t.data) > 0.3);
  t.date(badix) = [];
  t.data(badix) = [];

  if ( ~exist('tfix','var') || isempty(tfix) )
    tfix = 1:length(t.data);
  end;
  if ( ~exist('ffix','var') || isempty(ffix) )
    ffix = 1:length(f.data);
  end;

  scatter_fit_ts(t,f,tfix,ffix,'dT/dt','Heat flux','resid',true);

return;
