function [dx,dy] = degrees_to_meters(lon1, lon2, lat1, lat2)

  dx = sw_dist([lat1 lat1], [lon1 lon2], 'km') * 1e3;
  dy = sw_dist([lat1 lat2], [lon1 lon1], 'km') * 1e3;

return;
