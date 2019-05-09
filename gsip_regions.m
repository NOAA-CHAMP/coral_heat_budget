1;

regions = {'asam','freef','ecarib','gbr'};
for ix = 1:length(regions)

  region = regions{ix};
  dx = 360/4096;

  switch ( lower(region) )
   % case 'world',  lonoff =   -7;  latoff =   -7;  lonlen = 4096;  latlen = 2048;
   case 'asam',   lonoff = 2018;  latoff =  687;  lonlen =  256;  latlen =  256;
   case 'freef',  lonoff = 2986;  latoff = 1109;  lonlen =  256;  latlen =  256;
   case 'ecarib', lonoff = 3190;  latoff = 1029;  lonlen =  256;  latlen =  256;
   case 'gbr',    lonoff = 1529;  latoff =  653;  lonlen =  256;  latlen =  256;
   otherwise,     error('Unrecognized region "%s"!', region);
  end;

  lonoff = lonoff + 8;
  latoff = latoff + 8;

  lon = (([lonoff:(lonoff+lonlen-1)] - 1) * dx) + (dx/2);
  lat = (([latoff:(latoff+latlen-1)] - 1) * dx) - 90 + (dx/2);

  lon(lon >= 180) = lon(lon >= 180) - 360;

  region,
  lon([1 end]),
  lat([1 end]),

end;
