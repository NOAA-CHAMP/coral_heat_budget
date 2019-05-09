1;

bbox = [-85.6 -79.5 24.0 26.7];

map_sofla(bbox, [-30 -200]);

plot_seakeys_icon(gca, {'LKWF1','FWYF1','MLRF1','LONF1','SMKF1','SANF1','DRYF1','42003'});

x = cosd(bbox(3));
y = abs(bbox(3)-bbox(4)) / abs(bbox(1)-bbox(2));
set(gcf,'units','normalized','outerposition',[0 0 x y]);

title('Region of proposed research');

print('-dpng', '-loose', '../figs/region-map.png');
