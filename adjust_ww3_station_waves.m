function stn = adjust_ww3_station_waves(stn)
%function stn = adjust_ww3_station_waves(stn)
%
% Apply empirical adjustments to NOAA WaveWatch III model SURFACE WAVE data
% to better match in situ data from Haus et al. studies in March-June 2005.
%
% Last Saved Time-stamp: <Fri 2012-03-09 20:53:32  Lew.Gramer>


  % Create new fields for adjusted values
  stn.ww3_peakwaveper_adj = stn.ww3_peakwaveper;
  stn.ww3_sigwavehgt_adj = stn.ww3_sigwavehgt;
  stn.ww3_peakwavedir_adj = stn.ww3_peakwavedir;


  % % Initial adjustment based on simplistic linear regressions
  % stn.ww3_peakwaveper_adj.data = ((stn.ww3_peakwaveper.data.*1.00) + 0.00);
  % stn.ww3_sigwavehgt_adj.data = ((stn.ww3_sigwavehgt.data.*0.63) + 0.25);  %r^2~0.54, RMSE~0.2
  % stn.ww3_peakwavedir_adj.data = ((stn.ww3_peakwavedir.data.*1.00) + 0.00);

  % LINEAR ADJUSTED vs. IN SITU Hs FIT:
  %  Wv3: +0.000+1.00x r2~0.54 RMSE~0.197
  %  Wv4: -0.007+0.74x r2~0.29 RMSE~0.248
  %  Wv7: -0.073+1.03x r2~0.52 RMSE~0.190



  %%%%
  %% Significant wave height

  % LINEAR ADJUSTED vs. IN SITU FIT:
  %  Wv3: +0.006+1.01x r2~0.73 RMSE~0.147
  %  Wv4: -0.069+0.84x r2~0.44 RMSE~0.216
  %  Wv7: -0.052+1.03x r2~0.74 RMSE~0.146

  % WW3 significant wave height error AFTER empirical corrections using
  % Haus study sites: Site 3 slope 1.01, bias 0.01, r^2=0.72, RMSE~0.15;
  % Site 4: 0.84, -0.07, 0.44, 0.22; Site 7: 1.03, -0.05, 0.74, 0.15.

  % Piecewise breakpoint in Peak Wave Period [s]
  cutoff_per = 6.5;

  % Regression against all seasons - piecewise by wave PERIOD
  [chopix,ig] = intersect_dates(stn.ww3_sigwavehgt.date,stn.ww3_peakwaveper.date(stn.ww3_peakwaveper.data<=cutoff_per));
  a = +0.115; 
  b = +0.892;  % RMSE ~ 0.15[m], r^2 ~ 0.77
  stn.ww3_sigwavehgt_adj.data(chopix) = ( (stn.ww3_sigwavehgt.data(chopix) .* b) + a );

  [rollix,ig] = intersect_dates(stn.ww3_sigwavehgt.date,stn.ww3_peakwaveper.date(stn.ww3_peakwaveper.data>cutoff_per));
  % % Linear ROBUST
  % a = +0.330; 
  % b = +0.336;  % RMSE ~ 0.15[m], r^2 ~ 0.54
  % stn.ww3_sigwavehgt_adj.data(rollix) = ( (stn.ww3_sigwavehgt.data(rollix) .* b) + a );

  % Power2 CFIT
  a = +0.432; 
  b = +0.723;
  c = +0.243;  % RMSE ~ 0.15[m], r^2 ~ 0.54
  stn.ww3_sigwavehgt_adj.data(rollix) = ( ( a .* (stn.ww3_sigwavehgt.data(rollix) .^ b) ) + c );


  % Sanity check
  stn.ww3_sigwavehgt_adj.data(stn.ww3_sigwavehgt_adj.data<0) = 0;

return;
