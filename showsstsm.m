1;

stanm = 'fwyf1';
year = 2005;
period = 'spring';

fprintf('%s AVHRR SST SOM - %04d %s\n', stanm, year, period);

pathroot = '.';
if ( ~exist('figspath', 'var') || isempty(figspath) )
  figspath = fullfile(pathroot, '../figs', '');
end;
clear pathroot;

switch (lower(period(1:3)))
 case {'win'}, indts = datenum(year,01,01):datenum(year,03,31);
 case {'spr'}, indts = datenum(year,04,01):datenum(year,06,30);
 case {'sum'}, indts = datenum(year,07,01):datenum(year,09,30);
 case {'aut'}, indts = datenum(year,10,01):datenum(year,12,31);
 case {'dry'}, indts = datenum(year,01,01):datenum(year,06,30);
 case {'wet'}, indts = datenum(year,07,01):datenum(year,12,31);
 case {'ann'}, indts = datenum(year,01,01):datenum(year,12,31);
 otherwise,    error('Unrecognized seasonal period "%s"!', period);
end;

fname = sprintf('%s-ssts-%04d-%s.mat', stanm, year, period);
% [dts,ssts] = anavhrr(indts, stanm, 'florida');
% save(fname, 'ssts', 'dts');
% load(fname, 'ssts', 'dts');

fname = sprintf('%s-sst-%04d-%s.mat', stanm, year, period);
% sstdims = [size(ssts,2), size(ssts,3)];
% sst = reshape(ssts, [size(ssts,1) (sstdims(1)*sstdims(2))]);
% sst = sst - repmat(nanmean(sst,2), [1 size(sst,2)]);
% save(fname, 'sst', 'sstdims');
% load(fname, 'sst', 'sstdims');
load(fname, 'sstdims');


clear ssts; % This sucker is huge - can't leave it lying around...


mapdims = [4 4];

fname = sprintf('%s-sst-sm-%04d-%s.mat', stanm, year, period);
% % If a pixel is NaN for every image in our dataset, it's LAND
% landmask = all(isnan(sst), 1)';
% sm = som_make(sst, 'msize',mapdims, 'randinit', 'shape','sheet', ...
%               'mask',double(~landmask), 'training','long', 'neigh','ep');
% [sm.bmus,sm.qerrs] = som_bmus(sm, sst);
% sm.hits = som_hits(sm,sst);
% save(fname, 'sm');
load(fname, 'sm');

clear sst; % Ditto...

% Not sure why somtoolbox doesn't do this automatically...
sm.codebook(1:size(sm.codebook,1), ~logical(sm.mask)) = nan;

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
  title(sprintf( 'N = %d (%.0f%%)', sm.hits(ix), (sm.hits(ix)*100/length(sm.bmus)) ));
end;
ax = suplabel(sprintf('%s SOM "Modes" [%d %d], USF AVHRR SST %04d %s', ...
                      upper(stanm), mapdims(1), mapdims(2), year, period), 't');
pfname = fullfile(figspath, ...
                  sprintf( '%s-SOM-AVHRR-%04d-%s-%dx%d', ...
                           stanm, year, period, mapdims(1), mapdims(2) ));
print('-dtiff', '-r300', [pfname '.tiff']);

fprintf('SOM used %d synoptic SST fields\n', length(sm.bmus));

clear cdbk ix ax pfname;
