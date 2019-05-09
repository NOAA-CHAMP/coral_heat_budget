1;

boundingbox = [ -80.50 -79.10 +24.80 +25.90 ];

lon = boundingbox(1):0.02:boundingbox(2);
lat = boundingbox(3):0.02:boundingbox(4);

[LONS, LATS] = meshgrid(lon, lat);

[U, V] = meanflow(LONS, LATS, 0.0, 0.0);

% % Impose a Florida Current-like flow field over whole "Straits"
% [U, V] = eddy(LONS, LATS, U, V, -2.0, 3, -81.0, 26.0, 350);

% % Coastal counter-current inshore of the "Florida Current"
% [U, V] = eddy(LONS, LATS, U, V, 0.5, 'sin', -81.0, 26.0, 220);

% Cyclonic (in NHem) "frontal eddy"
[U, V] = eddy(LONS, LATS, U, V, -0.25, 0.5, -79.90, 25.25, 20);

% % Anti-cyclonic (in NHem) "frontal eddy"
% [U, V] = eddy(LONS, LATS, U, V, +0.25, 0.5, -80.00, 25.10, 20);

% Internal wave train(s) propagating offshore of the coast
for lambda = 18:1:22
  [U, V] = cwave(LONS, LATS, U, V, 0.10, -79.5, [], lambda);
end;

% Add stochastic noise to the whole vector field
[U, V] = randflow(LONS, LATS, U, V, 0.01);


% Calculate vector-field properties

[dx,dy] = degrees_to_meters(LONS(1,1), LONS(1,2), LATS(1,1), LATS(2,1));
[zetar, divr, strn, strs, dfrm, owp] = okubo_weiss(U, V, dx, dy);


% Plot results

figure;
quiver(LONS, LATS, U, V, 'color', 'black');
title('Currents');

fh = figure;
hold on;
% surf(LONS, LATS, zetar);
contourf(LONS, LATS, zetar);
colorbar;
% quiver(LONS, LATS, U, V, 'color', 'black');
xlim(boundingbox(1:2)); ylim(boundingbox(3:4));
btlim = 0e-5;
uplim = 10e-5;
% uplim = mean(zetar(:)) - 0.2 * std(zetar(:));
set(gca, 'clim', [btlim uplim]);
title('\zeta');
peakixes = outline_peaks(fh, LONS, LATS, zetar, [5e-5]);

% figure;
% hold on;
% % surf(LONS, LATS, strn);
% contourf(LONS, LATS, strn);
% colorbar;
% % quiver(LONS, LATS, U, V, 'color', 'black');
% xlim(boundingbox(1:2)); ylim(boundingbox(3:4));
% title('Normal Strain');

% figure;
% hold on;
% % surf(LONS, LATS, strs);
% contourf(LONS, LATS, strs);
% colorbar;
% % quiver(LONS, LATS, U, V, 'color', 'black');
% xlim(boundingbox(1:2)); ylim(boundingbox(3:4));
% title('Stress Strain');

% figure;
% hold on;
% % surf(LONS, LATS, dfrm);
% contourf(LONS, LATS, dfrm);
% colorbar;
% % quiver(LONS, LATS, U, V, 'color', 'black');
% xlim(boundingbox(1:2)); ylim(boundingbox(3:4));
% title('Deformation');

fh = figure;
hold on;
% surf(LONS, LATS, owp);
contourf(LONS, LATS, owp);
colorbar;
% quiver(LONS, LATS, U, V, 'color', 'black');
xlim(boundingbox(1:2)); ylim(boundingbox(3:4));
btlim = -10e-9;
uplim = -2e-9;
% uplim = mean(owp(:)) - 0.2 * std(owp(:));
set(gca, 'clim', [btlim uplim]);
title('OWP');
peakixes = outline_peaks(fh, LONS, LATS, owp, [uplim]);
