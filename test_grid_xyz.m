1;

m = csvread('../data/isobaths-mlrf1.csv');

dlon = 8.3300e-004;
dlat = 8.3300e-004;

lons = min(m(:,1)):dlon:max(m(:,1));
lats = min(m(:,2)):dlat:max(m(:,2));

[LON,LAT] = meshgrid(lons, lats);

Z = griddata(m(:,1),m(:,2),m(:,3), LON,LAT);
Z = cast(Z, 'single');

figure;
contour(LON,LAT,Z,[-2:-2:-10 -15:-5:-30 -100:-50:-200]);
