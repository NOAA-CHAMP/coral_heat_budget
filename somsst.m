1;

format short; format compact;

more_status = get(0, 'More'); more('off');


datapath = '../data';

figspath = '../figs';


stanm = 'smkf1';

dataset = 'mean';


thisweek = floor((now - datenum(2009,1,1) + 1)/7);
yrs = repmat(1994:2008', [52 1]);
yrs = [repmat(1993, [1 14]) yrs(:)' repmat(2009, [1 thisweek])];
wks = [39:52 repmat(1:52, [1 15]) 1:thisweek];


% period = 'Summer';
period = 'Annual';

fprintf('Period is %s\n', period);
switch ( lower(period) )
 case 'annual',
  % noop
 case 'winter',
  seasonix = find(ismember(wks, [1:12 50:52]));
  yrs = yrs(seasonix);
  wks = wks(seasonix);
 case 'spring',
  seasonix = find(ismember(wks, [11:25]));
  yrs = yrs(seasonix);
  wks = wks(seasonix);
 case 'summer',
  seasonix = find(ismember(wks, [24:38]));
  yrs = yrs(seasonix);
  wks = wks(seasonix);
 case 'autumn',
  seasonix = find(ismember(wks, [37:51]));
  yrs = yrs(seasonix);
  wks = wks(seasonix);
 otherwise,
  error('Unrecognized seasonal period "%s"!', period);
end;


[emptysst, LONS, LATS] = ansst(stanm, 'bounds');
sstdims = size(emptysst);
clear emptysst;

if ( ~exist('sst', 'var') )

  sstfname = sprintf('%s-sst-%s-%04d-%04d-%s-%04d-%04d.mat', ...
                     stanm, dataset, yrs([1 end]), lower(period), ...
                     sstdims([1 2]));
  sstfname = fullfile(datapath, sstfname);

  if ( exist(sstfname, 'file') )
    load(sstfname, 'sst');
    fprintf('Reloading SST weekly composites from %s\n', sstfname);

  else
    fprintf('Loading weekly SSTs');

    oldyr = 0;
    for ix = 1:length(wks)
      yr = yrs(ix);
      wk = wks(ix);
      if ( yr ~= oldyr )
        fprintf('\n\t%04d:', yr);
        oldyr = yr;
      end;

      fprintf(' %02d', wk);
      x = ansst(stanm, yr, wk, dataset);
      % Turn matrix of pixels into a vector of variables
      sst(ix, 1:numel(x)) = x(:);
    end;
    fprintf('\n');
    save(sstfname, 'sst');

  end;

  % A week when every pixel was NaN is a missing or bad image - remove it
  rmrows = all(isnan(sst), 2);
  sst(rmrows, :) = [];


  % If a pixel is NaN for every week in our dataset, it's probably LAND
  landmask = all(isnan(sst), 1)';


  % 'goodsst' is 'sst' with all NaNs (masked values) set to 0
  goodsst = sst;
  goodsst(~isfinite(sst)) = 0;

  meandim = 'none';

%   % Remove the time mean from SSTs (a la Mariano et al, 2006)
%   meandim = 'space'; % i.e., the mean is a FUNCTION OF x,y
%   tsmean = repmat( mean(goodsst,1), [size(goodsst,1) 1] );
%   sst = sst - tsmean;
%   goodsst = goodsst - tsmean;

%   % Remove the spatial mean from SSTs (a la Mariano et al, 2006)
%   meandim = 'time'; % i.e., the mean is a FUNCTION OF time
%   tsmean = repmat( mean(goodsst,2), [1 size(goodsst,2)] );
%   sst = sst - tsmean;
%   goodsst = goodsst - tsmean;


%   % Remove spatial AND temporal mean from SSTs (a la Weisberg and He, 2006)
%   meandim = 'space(time)';
%   tsmean = repmat( mean(goodsst,1), [size(goodsst,1) 1] );
%   sst = sst - tsmean;
%   goodsst = goodsst - tsmean;
%   tsmean = repmat( mean(goodsst,2), [1 size(goodsst,2)] );
%   sst = sst - tsmean;
%   goodsst = goodsst - tsmean;


end;


% mapdims = [8 3];
% mapdims = [6 3];
% mapdims = [4 3];
mapdims = [3 3];
% mapdims = [3 2];

nmaps = mapdims(2) * mapdims(1);


if ( ~exist('sm', 'var') )

  somfname = sprintf('%s-som-sst-%s-%04d-%04d-%s-%04d-%04d-%d-%d.mat', ...
                     stanm, dataset, yrs([1 end]), lower(period), ...
                     sstdims([1 2]), mapdims([1 2]));
  somfname = fullfile(datapath, somfname);

  if ( exist(somfname, 'file') )

    load(somfname, 'sm', 'results', 'bmus', 'qerrs');
    fprintf('Reloaded SOM result from %s\n', somfname);


  else

    % Random initialization - slower runs, but LESS MEMORY (and NaNs are OK)
    sm = som_make(sst, 'msize',mapdims, 'randinit', 'shape','sheet', ...
                  'mask',double(~landmask), 'training','long', 'neigh','ep');

    % % Linear init. - faster SOM runs, but MORE MEMORY (and must coddle NaNs)
    % sm = som_make(goodsst, 'msize',mapdims, 'lininit', 'shape','sheet', ...
    %               'mask',double(~landmask), 'training','long', 'neigh','ep');

    results = sm.codebook;

    % Reset all landmasks to NaN in our 'results'
    results(1:size(results,1), landmask) = nan;

    % Calculate the best-matched Unit or "mode" for each week of data
    fprintf('Calculating Best-Match-Units\n');
    [bmus, qerrs] = som_bmus(sm, sst);
    % [bmus, qerrs] = som_bmus(sm, goodsst);

    save(somfname, 'sm', 'results', 'bmus', 'qerrs');
    fprintf('SOM completed and saved\n');

  end;

end;


fprintf('Plotting SOM "modes"\n');

% figure;
% % Do not try to show all 'variables' - there are thousands of them!
% som_show(sm, 'umat', 'all');
% som_trajectory(bmus);


minc = min(results(isfinite(results)));
maxc = max(results(isfinite(results)));

% Calculate number of matched weeks for each "mode"
wksmatched = som_hits(sm, sst);
pctmatched = wksmatched * (100 / length(bmus));

% Sort BMUs by number of weeks, from most-matched to least
[ign, mode_order] = sort(wksmatched, 'descend');


%
% Plot SOM results in order from most "hits" (matching samples) to least
%

figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
for orderix = 1:length(mode_order)
  % Plot "modes" from most frequently matched to least
  ix = mode_order(orderix);
  subplot(mapdims(2), mapdims(1), orderix);

  pcolor(LONS, LATS, reshape(results(ix, :), sstdims));
  shading('flat');
  % For pcolor and similar, let user see SST values when in datacursor mode
  set_pcolor_cursor;
%   set(gca, 'zlim', [minc maxc]);
%   set(gca, 'clim', [minc maxc]);

  % contourf(LONS, LATS, reshape(results(ix, :), sstdims));
  % set(gca, 'clim', [minc maxc]);

  title(sprintf('#%d (%d wks - %0.1f%%)', ...
                ix, wksmatched(ix), pctmatched(ix)));
end;
% suptitle or suplabel, from MATLAB online contrib
ax = suplabel(sprintf('%s SOM "Modes" vs. mean(%s), 1993-2009 %s (%d weeks)', ...
                      upper(stanm), meandim, period, size(sst,1)), 't');
% set(ax,'clim',[minc maxc]);
% colorbar('peer', ax);
print('-dpng', fullfile(figspath, ...
                        sprintf('%s-SOM-%s-%dx%d.png', stanm, period, ...
                                mapdims(1), mapdims(2))));

%
% Compute Principal Components / EOFs
%

fprintf('Doing complementary Principal Component Analysis\n');

% [v, d, z] = eof(goodsst');
covm = cov(goodsst');
[evecs, evals, explained] = pcacov(covm);
pcadat = goodsst.' * evecs;
% PCA scaling is arbitrary - this "fix" suggested by Arthur Mariano
pcadat = -pcadat';

% (For now, use same scale as SOM outputs) NOT
minc = min(pcadat(:));
maxc = max(pcadat(:));

pcadat(~isfinite(sst)) = nan;

figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
% Try to show the same number of PCA modes as we did SOM maps
totpct = 0;
ncomps = min(size(pcadat,1), (mapdims(2)*mapdims(1)));
for ix = 1:ncomps
  subplot(mapdims(2), mapdims(1), ix);

  pcolor(LONS, LATS, reshape(pcadat(ix, :), sstdims));
  shading('flat');
  % For pcolor and similar, let user see SST values when in datacursor mode
  h = datacursormode(gcf);
  set(h, 'UpdateFcn', @xyzc_select_cb);
%   set(gca, 'zlim', [minc maxc]);
%   set(gca, 'clim', [minc maxc]);

  % contourf(LONS, LATS, reshape(pcadat(ix, :), sstdims));
  % set(gca, 'clim', [minc maxc]);

  title(sprintf('#%d (%0.1f%% %s)', ix, explained(ix), '\sigma^2'));
  totpct = totpct + explained(ix);
end;
ax = suplabel(sprintf('%s PCA Modes vs. mean(%s), 1993-2009 %s (%3.1f%% of %s)', ...
                      upper(stanm), meandim, period, totpct, '\sigma^2'), 't');
% set(ax,'clim',[minc maxc]);
% colorbar('peer', ax);
print('-dpng', fullfile(figspath, ...
                        sprintf('%s-PCA-%s-%dx%d.png', stanm, period, ...
                                mapdims(1), mapdims(2))));

%%%% DEBUG: Do not do these 'clear' calls for now...
% clear goodsst;
% clear covm;
% clear pcadat;

more(more_status);
clear more_status;
