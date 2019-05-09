1;

stanm = 'smkf1';
year = 2005;
sea1 = 'wet';
sea2 = 'dry';

load(sprintf('%s-sst-sm-%04d-%s.mat', stanm, year, sea1));
load(sprintf('%s-sst-%04d-%s.mat', stanm, year, sea2), 'sst', 'sstdims');

mapdims = [4 4];

% If a pixel is NaN for every image in our dataset, it's LAND
landmask = all(isnan(sst), 1)';

[dlen ig] = size(sst);
sTrain = som_train_struct(sm,'dlen',dlen,'phase','finetune');
sTrain = som_set(sTrain,'data_name','sst','algorithm','batch');
sTrain.trainlen = sTrain.trainlen*4;
sm = som_batchtrain(sm,sst,sTrain,'tracking',1,'mask',double(~landmask));

sm.hits = sm.hits + som_hits(sm,sst);

sm.codebook(1:size(sm.codebook,1), landmask) = nan;

figure;
set(gcf, 'units','normalized', 'outerposition',[0 0 1 1]);
for ix = 1:(mapdims(1)*mapdims(2))
  subplot(mapdims(1), mapdims(2), ix);
  cdbk = reshape(sm.codebook(ix,:), sstdims);
  pcolor(flipud(cdbk));
  shading('flat');
  % set(gca,'clim',[22 32]);
  set(gca,'clim',[-2 +2]);
  set_pcolor_cursor;
  title(sprintf( 'N = %d', sm.hits(ix) ));
end;
ax = suplabel(sprintf('%s SOM "Modes" [%d %d], USF AVHRR SST %04d %s then %s', ...
                      upper(stanm), mapdims(1), mapdims(2), year, sea1, sea2), 't');
pfname = sprintf('../figs/%s-SOM-AVHRR-%04d-%s-then-%s-%dx%d.tiff', ...
                stanm, year, sea1, sea2, mapdims(1), mapdims(2));
print('-dtiff', '-r300', pfname);
