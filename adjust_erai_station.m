function stn = adjust_erai_station(stn)
%function stn = adjust_erai_station(stn)
%
% Apply empirical adjustments to ERA-Interim reanalysis data to better match
% SEAKEYS in situ data throughout the year.
%
% Last Saved Time-stamp: <Sat 2012-06-02 17:53:23  Lew.Gramer>


  % Create new fields for adjusted values
  stn.raw_erai_dsrf_adj = stn.raw_erai_dsrf;
  stn.raw_erai_dlrf_adj = stn.raw_erai_dlrf;
  stn.raw_erai_precip_adj = stn.raw_erai_precip;
  stn.erai_dsrf_adj = stn.erai_dsrf;
  stn.erai_dlrf_adj = stn.erai_dlrf;
  stn.erai_precip_adj = stn.erai_precip;

  %%%%
  %% Downward shortwave radiative flux (insolation)
  % % Regression against all seasons
  % a = -10.00; 
  % b = .7877;
  % stn.raw_erai_dsrf_adj.data = ( (stn.raw_erai_dsrf.data .* b) + a );
  % stn.erai_dsrf_adj.data = ( (stn.erai_dsrf.data .* b) + a );

  % % Regression of ONE-DAY AVERAGES
  % as = [+10.22,-60.06,-24.46, -6.48];
  % bs = [0.9533,1.1551,1.0997,1.1056];
  % for seasix = 1:4
  %   a = as(seasix);
  %   b = bs(seasix);
  %   ix = find(get_season(stn.raw_erai_dsrf.date)==seasix);
  %   stn.raw_erai_dsrf_adj.data(ix,1) = ( (stn.raw_erai_dsrf.data(ix) .* b) + a );
  %   ix = find(get_season(stn.erai_dsrf.date)==seasix);
  %   stn.erai_dsrf_adj.data(ix,1) = ( (stn.erai_dsrf.data(ix) .* b) + a );
  % end;

  % % Regression of hourly mid-day (TS_FLORIDA_MIDDAY) values
  % as = [181.53,253.92,205.95,165.41];
  % bs = [0.7990,0.7825,0.8486,0.7965];
  % for seasix = 1:4
  %   a = as(seasix);
  %   b = bs(seasix);
  %   ix = find(get_season(stn.raw_erai_dsrf.date)==seasix);
  %   stn.raw_erai_dsrf_adj.data(ix,1) = ( (stn.raw_erai_dsrf.data(ix) .* b) + a );
  %   ix = find(get_season(stn.erai_dsrf.date)==seasix);
  %   stn.erai_dsrf_adj.data(ix,1) = ( (stn.erai_dsrf.data(ix) .* b) + a );
  % end;

  % % Regression of hourly daylight (GET_DAYLIGHT) values
  % as = [85.127,81.437,78.319,69.435];
  % bs = [0.7390,0.7601,0.8036,0.7665];
  % for seasix = 1:4
  %   a = as(seasix);
  %   b = bs(seasix);
  %   ix = find(get_season(stn.raw_erai_dsrf.date)==seasix);
  %   stn.raw_erai_dsrf_adj.data(ix,1) = ( (stn.raw_erai_dsrf.data(ix) .* b) + a );
  %   ix = find(get_season(stn.erai_dsrf.date)==seasix);
  %   stn.erai_dsrf_adj.data(ix,1) = ( (stn.erai_dsrf.data(ix) .* b) + a );
  % end;

  % % Regression of DAILY AVERAGE against all seasons (FWYF1)
  % a = -21.56; 
  % b = 0.9555;
  % Regression of DAILY AVERAGE against all seasons (VKAF1)
  a = 0.45; 
  b = 0.9091;

  ix = find(stn.raw_erai_dsrf_adj.data > 100);
  stn.raw_erai_dsrf_adj.data(ix) = ( (stn.raw_erai_dsrf.data(ix) .* b) + a );

  ix = find(stn.erai_dsrf_adj.data > 100);
  stn.erai_dsrf_adj.data(ix) = ( (stn.erai_dsrf.data(ix) .* b) + a );

  % Sanity check
  stn.raw_erai_dsrf_adj.data(stn.raw_erai_dsrf_adj.data<0) = 0;
  stn.erai_dsrf_adj.data(stn.erai_dsrf_adj.data<0) = 0;


  %%%%
  %% Downward longwave radiative flux
  as = [ -8.48,+30.61,104.55,-34.83];
  bs = [1.0185,0.9160,0.7503,1.0835];
  for seasix = 1:4
    a = as(seasix);
    b = bs(seasix);
    ix = find(get_season(stn.raw_erai_dlrf.date)==seasix);
    stn.raw_erai_dlrf_adj.data(ix,1) = ( (stn.raw_erai_dlrf.data(ix) .* b) + a );
    ix = find(get_season(stn.erai_dlrf.date)==seasix);
    stn.erai_dlrf_adj.data(ix,1) = ( (stn.erai_dlrf.data(ix) .* b) + a );
  end;

  % Sanity check
  stn.raw_erai_dlrf_adj.data(stn.raw_erai_dlrf_adj.data<0) = 0;
  stn.erai_dlrf_adj.data(stn.erai_dlrf_adj.data<0) = 0;


  %%%%
  %% Precipitation
  stn.raw_erai_precip_adj.data(stn.raw_erai_precip_adj.data<0) = 0;
  stn.erai_precip_adj.data(stn.erai_precip_adj.data<0) = 0;


return;
