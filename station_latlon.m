function [stn_lat, stn_lon] = station_latlon(stanm)
%function [stn_lat, stn_lon] = station_latlon(stanm)
%
% error('This function has been replaced by GET_STATION_COORDS.m!');

error('This function has been replaced by GET_STATION_COORDS.m!');

    switch ( lower(stanm) )
     case 'fwyf1',
      stn_lat = 25.590; stn_lon = -80.100;
     case 'cryf1',
      stn_lat = 25.222; stn_lon = -80.212;
     case 'mlrf1',
      stn_lat = 25.010; stn_lon = -80.380;
     case 'lonf1',
      stn_lat = 24.840; stn_lon = -80.860;
     case 'tnrf1',
      stn_lat = 24.750; stn_lon = -80.78333;
     case 'smkf1',
      stn_lat = 24.630; stn_lon = -81.110;
     case 'amsf1',
      stn_lat = 24.525; stn_lon = -81.520;
     case 'sanf1',
      stn_lat = 24.460; stn_lon = -81.880;
     case 'plsf1',
      stn_lat = 24.690; stn_lon = -82.770;
     otherwise,
      error('Unrecognized station name "%s"!', stanm);
    end;

return;
