1;

% OPeNDAP/DODS Data URL:  http://nomads.ncep.noaa.gov:9090/dods/ofs/ofs20100503/hourly/rtofs_nowcast_atl 
% Description: Real-time Ocean Forecast System - Altantic Sector Native GRID Beta Version
% Longitude: -105.04000000000°E to 67.07750000000°E  (1300 points, avg. res. 0.132°)  
% Latitude: -25.20000000000°N to 76.43220000000°N  (1500 points, avg. res. 0.068°)  
% Time: 00Z02MAY2010 to 23Z02MAY2010  (24 points, avg. res. 0.042 days)  
% Variables: (total of 9)
%  mixhtsfc ** surface mixed layer depth [m] 
%  mntsfsfc ** surface montgomery stream function [m2/s] 
%  salinhl1p5 ** 1 hybrid level - 2 hybrid level 3-d salinity [-] 
%  sshgsfc ** surface sea surface height relative to geoid [m] 
%  ubaroocn ** entire ocean (considered as a single layer) barotropic u velocity [m/s] 
%  uogrdhl1p5 ** 1 hybrid level - 2 hybrid level u-component of current [m/s] 
%  vbaroocn ** entire ocean (considered as a single layer) barotropic v velocity [m/s] 
%  vogrdhl1p5 ** 1 hybrid level - 2 hybrid level v-component of current [m/s] 
%  wtmpchl1p5 ** 1 hybrid level - 2 hybrid level 3-d temperature [degc] 

minlon = -105.0400;
dlon = 0.1325;
maxlon = 67.0775;

minlat = -25.2000;
dlat = 0.0678;
maxlat = 76.4322;

dtime = (1/24);

x = loaddap('http://nomads.ncep.noaa.gov:9090/dods/ofs/ofs20100429/hourly/rtofs_nowcast_atl.dods?mixhtsfc[0:0][114:118][277:281]');
