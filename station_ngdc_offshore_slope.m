function stn = station_ngdc_offshore_slope(stn_or_stnm)
%function stn = station_ngdc_offshore_slope(stn_or_stnm)
%
% Estimate cross-shore slope (rise / run) of each SEAKEYS station in the
% direction of greatest descent (normally, offshore). Result was the average
% of using two sources: NGDC 3-arcsecond resolution Coastal Relief Model, and
% Google Earth bathymetry (combination of Scripps, NOAA, US Navy, NGA, and
% GEBCO, the General Bathymetric Chart of the Oceans).
%
% Last Saved Time-stamp: <Wed 2013-04-10 17:07:18 Eastern Daylight Time gramer>

  stn = get_station_from_station_name(stn_or_stnm);

  switch ( lower(stn.station_name) ),
   case 'lkwf1', stn.ngdc_offshore_slope = 0.02;
   case 'fwyf1', stn.ngdc_offshore_slope = 0.04;
   case 'mlrf1', stn.ngdc_offshore_slope = 0.03;
   case 'smkf1', stn.ngdc_offshore_slope = 0.02;
   case 'looe1', stn.ngdc_offshore_slope = 0.04;
   case 'sanf1', stn.ngdc_offshore_slope = 0.02;

   % case 'lonf1', stn.ngdc_offshore_slope = 0.003;
   % case 'dryf1', stn.ngdc_offshore_slope = 0.003;
   % case 'plsf1', stn.ngdc_offshore_slope = 0.007;

   % Rechecked with 2km search radius in NGDC bathymetry
   case 'lonf1', stn.ngdc_offshore_slope = 0.002;
   case 'dryf1', stn.ngdc_offshore_slope = 0.004;
   case 'plsf1', stn.ngdc_offshore_slope = 0.007;

   % NCORE C and Key Largo, estimated from NGDC bathymetry 2013 Apr 10
   case 'ncorc', stn.ngdc_offshore_slope = 0.02;
   case 'klgf1', stn.ngdc_offshore_slope = 0.02;

   case 'bnpin', stn.ngdc_offshore_slope = 0.004;
   case 'bnpmi', stn.ngdc_offshore_slope = 0.007;
   case 'bnpon', stn.ngdc_offshore_slope = 0.020;
   case 'bnpnn', stn.ngdc_offshore_slope = 0.007;
   case 'bnppa', stn.ngdc_offshore_slope = 0.042;

   case '42003', stn.ngdc_offshore_slope = 0.0001;

   otherwise,    error('Unrecognized station name "%s"!',stn.station_name);
  end;

return;
