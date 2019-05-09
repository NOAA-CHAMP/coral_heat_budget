1;

map_sofla([-84 -79 24 27]);

disp('Hit Enter to begin mapping, then Enter again to stop...');
pause;

[X, Y] = ginput;

line(X, Y);

dist = 0;
for ix = 2:length(X)
    [d, ig] = sw_dist(Y([ix, ix-1]), X([ix, ix-1]), 'km');
    dist = dist + d;
end;

disp(dist);
