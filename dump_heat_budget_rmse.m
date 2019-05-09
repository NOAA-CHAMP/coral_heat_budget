function s = dump_heat_budget_rmse(s)
%function s = dump_heat_budget_rmse(s)
%
% See CALC_HEAT_BUDGET_RMSE.

  %fmg; hist(s.dtbqse.data,5000); xlim([-.025,.025]);

  disp([s.N,s.Nc]);
  disp([median(diff(s.raw_ts.date)),max(diff(s.raw_ts.date))]);
  disp([median(diff(s.raw_td.date)),max(diff(s.raw_td.date))]);
  disp([s.qrmse, s.qrmsec]);
  disp([s.bqrmse,s.bqrmsec]);
  disp([s.dtrmse,s.dtrmsec]);
  disp([s.hcrmse,s.hcrmsec,s.optim.climerror]);

return;
