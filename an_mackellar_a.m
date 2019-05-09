1;

% MLRF1:
% Mid-day sea-surface albedo % (A)
% Summary: 3.47258      5.17058      6.02169      6.65006      19.9141      6.15611      1.48948
% Mid-day net insolation loss from water column, QSW*(1-gamma), Wm^-2
% Summary: 0.05738663      8.105462      28.06719       85.9336       262.865      54.02119      58.91194
% Mid-day total reflected insolation QSW*(A+(1-gamma)), Wm^-2
% Summary: 1.832059      36.87412      58.48329      114.1491      299.0979      82.82573       61.5217
% Mid-day total reflectivity % (A+1-gamma, compare MacKellar et al 2012)
% Summary: 5.76042       8.0772      11.8285      23.9222      42.6132       16.213       9.7201

% LONF1:
% Mid-day sea-surface albedo % (A)
% Summary: 3.46411      5.14318      5.99055      6.59074      16.7824      6.12505      1.49435
% Mid-day net insolation loss from water column, QSW*(1-gamma), Wm^-2
% Summary: 0.05242881      6.961466       24.4326      57.06756      165.9721      36.33613      35.35472
% Mid-day total reflected insolation QSW*(A+(1-gamma)), Wm^-2
% Summary: 1.222632      35.11221      52.81502      87.91299       204.583      65.02883      38.95395
% Mid-day total reflectivity % (A+1-gamma, compare MacKellar et al 2012)
% Summary: 5.72558      8.33793      10.9678      16.4286      26.7538      12.5406      4.89957

% DRYF1:
% Mid-day sea-surface albedo % (A)
% Summary: 3.46111      4.99952       5.9931      6.67666      22.7439       6.1591      1.73011
% Mid-day net insolation loss from water column, QSW*(1-gamma), Wm^-2
% Summary: 0.561308        32.919      78.15355      147.8576      292.8016      94.93276      72.31569
% Mid-day total reflected insolation QSW*(A+(1-gamma)), Wm^-2
% Summary: 2.004511      60.23679      105.8117      179.0734      328.4014      122.9419      76.05132
% Mid-day total reflectivity % (A+1-gamma, compare MacKellar et al 2012)
% Summary: 8.78027        14.99      21.7253      30.9961      42.9104      23.0177       8.3025


disp([upper(stn.station_name),':']);
disp('Mid-day sea-surface albedo % (A)');
a=ts_op(stn.erai_ndbc_srf,subset_ts(stn.erai_dsrf_adj,@ts_florida_midday),'/');
nansummary((1-a.data)*100),

disp('Mid-day net insolation loss from water column, QSW*(1-gamma), Wm^-2');
ag = ts_op(stn.erai_ndbc_srf,subset_ts(stn.absorbed_erai_ndbc_srf,@ts_florida_midday),'-');
nansummary(ag.data),

disp('Mid-day total reflected insolation QSW*(A+(1-gamma)), Wm^-2');
ag = ts_op(stn.erai_dsrf_adj,subset_ts(stn.absorbed_erai_ndbc_srf,@ts_florida_midday),'-');
nansummary(ag.data),

disp('Mid-day total reflectivity % (A+1-gamma, compare MacKellar et al 2012)');
ta = ts_op(ag,stn.erai_dsrf_adj,'/');
nansummary(ta.data*100), 
