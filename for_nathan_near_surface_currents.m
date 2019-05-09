1;

nc = mDataset('../data/hycom/matthieu.lehenaff/archv.2014_001_12.nc');
u = cast(nc{'u'}(:,:,:,:),'double');
v = cast(nc{'v'}(:,:,:,:),'double');
w = cast(nc{'w_velocity'}(:,:,:,:),'double');

z = cast(nc{'Depth'}(:),'double');
lat = cast(nc{'Latitude'}(:),'double');
lon = cast(nc{'Longitude'}(:),'double');
close(nc); clear nc

[LON,Z,LAT] = meshgrid(lon(200:201),-z,lat(200:201));
fmg; quiver3(LON,LAT,Z,u(1:3,200:201,200:201),v(1:3,200:201,200:201),w(1:3,200:201,200:201),0.001); view(3);
