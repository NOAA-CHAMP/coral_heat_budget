function plotbootmon(stats,absfld,accfld)
%function plotbootmon(stats,absfld,accfld)
%
% Plot bootstrap estimates of monthly mean for, e.g., sea temperature and
% cumulative heat flux. STATS is a struct built by, e.g., BOOTMONALL (v.).
%
% Last Saved Time-stamp: <Sun 2010-07-11 14:38:13  Lew.Gramer>

  if ( ~exist('absfld','var') || isempty(absfld) )
    absfld = 'ndbc_sea_t';
  end;
  if ( ~exist('accfld','var') || isempty(accfld) )
    accfld = 'netqf';
  end;

  % Need to multiply each monthly mean by hours/month, to accumulate
  dymos = [31 28 31 30 31 30 31 31 30 31 30 31];
  hrmos = 24*repmat(dymos, [1 (length(stats.(accfld).m(2,:))/12)]);

  figure;
  maxigraph;
  plot(1:12, ...
       [stats.(absfld).m(2,:) ; ...
        stats.(absfld).m(2,12) + (cumsum(stats.(accfld).m(2,:)).*hrmos)]);
  legend('T','dT/dt');

return;
