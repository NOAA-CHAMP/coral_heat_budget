1;

% Cobble up some pretend data
dep = linspace(0,2000,40)';
% Is my guess a good one - your working in the IND again, right? :)
lat = linspace(-40,+15,382)';
lon = linspace(+40,+130,320)';

dens = repmat(linspace(25,33,40)', [1 382 320]);
sals = 36 + randn([40 382 320]);

% Find depth CLOSEST TO our desired isopycnal at each grid point
[nearest_depth,nearestix] = min(abs(dens-27),[],1);

% Find salinity at that nearest depth, at each grid point
nearest_salin = squeeze(sals(nearestix));

% Results are two 382x320 matrices - one of depth, one of salinity
