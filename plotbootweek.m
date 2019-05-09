function plotbootweek(stats,absfld,accfld)
%function plotbootweek(stats,absfld,accfld)
%
% Plot bootstrap estimates of weekly mean for, e.g., sea temperature and
% cumulative heat flux. STATS is a struct built by, e.g., BOOTWEEKALL (v.).
%
% Last Saved Time-stamp: <Tue 2010-07-13 07:59:10  Lew.Gramer>

  if ( ~exist('absfld','var') || isempty(absfld) )
    absfld = 'ndbc_sea_t';
  end;
  if ( ~exist('accfld','var') || isempty(accfld) )
    accfld = 'netqf';
  end;

  % Need to multiply each weekly mean by hours/week, to accumulate
  hrswks = 24*7;

  figure;
  maxigraph;
  plot(1:52, ...
       [stats.(absfld).m(2,:) ; ...
        stats.(absfld).m(2,52) + (cumsum(stats.(accfld).m(2,:)).*hrswks)]);
  legend('T','dT/dt');

return;
