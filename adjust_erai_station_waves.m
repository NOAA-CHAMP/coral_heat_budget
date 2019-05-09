function stn = adjust_erai_station_waves(stn)
%function stn = adjust_erai_station_waves(stn)
%
% Apply empirical adjustments to ERA-Interim reanalysis SURFACE WAVE data to better match
% in situ data from Haus et al. studies in March-June 2005.
%
% Last Saved Time-stamp: <Fri 2012-03-09 20:52:45  Lew.Gramer>


  % Create new fields for adjusted values
  stn.erai_peakwaveper_adj = stn.erai_peakwaveper;
  stn.erai_sigwavehgt_adj = stn.erai_sigwavehgt;
  stn.erai_peakwavedir_adj = stn.erai_peakwavedir;


  %%%%
  %% Significant wave height

  % ERAI significant wave height error AFTER empirical corrections using
  % Haus study sites: Site 3 slope 1.04, bias 0.01, r^2=0.46, RMSE~0.21;
  % Site 4: 0.87, -0.08, 0.35, 0.23; Site 7: 0.97, -0.00, 0.46, 0.21.

  % Regression against all seasons - piecewise by wave PERIOD
  [chopix,ig] = intersect_dates(stn.erai_sigwavehgt.date,stn.erai_peakwaveper.date(stn.erai_peakwaveper.data<=6.5));
  a = -0.173; 
  b = +0.978;  % RMSE ~ 0.24[m], r^2 ~ 0.43
  stn.erai_sigwavehgt_adj.data(chopix) = ( (stn.erai_sigwavehgt.data(chopix) .* b) + a );

  [rollix,ig] = intersect_dates(stn.erai_sigwavehgt.date,stn.erai_peakwaveper.date(stn.erai_peakwaveper.data>6.5));
  % % Linear ROBUST
  % a = +0.386; 
  % b = +0.220;  % RMSE ~ 0.17[m], r^2 ~ 0.40
  % stn.erai_sigwavehgt_adj.data(rollix) = ( (stn.erai_sigwavehgt.data(rollix) .* b) + a );

  % Power2 CFIT
  a = -0.129; 
  b = -2.813;
  c = +0.833;  % RMSE ~ 0.14[m], r^2 ~ 0.60
  stn.erai_sigwavehgt_adj.data(rollix) = ( ( a .* (stn.erai_sigwavehgt.data(rollix) .^ b) ) + c );

  % Sanity check
  stn.erai_sigwavehgt_adj.data(stn.erai_sigwavehgt_adj.data<0) = 0;

return;
