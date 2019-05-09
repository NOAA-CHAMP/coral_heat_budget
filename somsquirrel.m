1;

somsst;

% h = som_show(sm);
h = som_show(sm, 'umat', 'all');
dat = sst; dat(~isfinite(dat)) = 0;
sd = som_data_struct(dat');
bmus = som_bmus(sm, sd);
som_trajectory(bmus);

close all;

h = som_show(sm);
