function [yrdts,yrdat,meandts,meandat] = anann(stn,fld)
%function [yrdts,yrdat,meandts,meandat] = anann(stn,fld)
%
% Analyze (with ANOVA1, MULTCOMPARE, and interannual ranges, qv.) the annual
% data distributions of multi-year time series STN.(FLD). Only analyze data
% from complete days, and from matching year days of all years; only analyze
% years with a minimum total sample count (6000 of 8760 possible hourly data)
% and sample density (no data gaps longer than 45 days).
%
% Last Saved Time-stamp: <Sat 2010-03-20 14:37:32 Eastern Daylight Time gramer>

  [yrs,mos,dys] = datevec(stn.(fld).date);
  newyrix = find(diff(yrs) > 0);

  % Use only complete years of data
  yrs = yrs(newyrix(1)+1:newyrix(end));
  mos = mos(newyrix(1)+1:newyrix(end));
  dys = dys(newyrix(1)+1:newyrix(end));

  uyrs = unique(yrs);

  for yrix = 1:length(uyrs)
    yr = uyrs(yrix);
    ix = find(yrs == yr);
    % We only want to compare years that aren't too screwed up
    if ( length(ix) > 6000 && max(diff(stn.(fld).date(ix))) < 45 )
      % For good years, use a strict hourly (interpolated) time series
      yrdts{yrix} = stn.(fld).date(ix(1)):(1/24):stn.(fld).date(ix(end));
      [yry,yrm,yrd,yrhr] = datevec(yrdts{yrix});
      firstix = find((yrhr == 0), 1, 'first');
      lastix = find((yrhr == 23), 1, 'last');

      yrdts{yrix} = yrdts{yrix}(firstix:lastix);
      [yry,yrm,yrd,yrhr] = datevec(yrdts{yrix});
      yrjds{yrix} = floor(yrdts{yrix}) - datenum(yry,1,1) + 1;
      yrdat{yrix} = interp1(stn.(fld).date(ix),stn.(fld).data(ix),yrdts{yrix});
      meandts(yrix) = yr;
      meandat(yrix) = nanmean(yrdat{yrix});
    else
      yrdts(yrix) = {[]};
      yrdat(yrix) = {[]};
      meandts(yrix) = yr;
      meandat(yrix) = nan;
    end;
  end;

  badix = find(isnan(meandat));
  yrdts(badix) = [];
  yrjds(badix) = [];
  yrdat(badix) = [];
  meandts(badix) = [];
  meandat(badix) = [];

  if ( numel(meandat) < 2 )
    error('After filtering, no good years were left to compare!');
  end;

  % Only compare means over MATCHING year-days of all good years
  jds = 1:366;
  for ix = 1:numel(yrdts)
    jds = intersect(jds,unique(yrjds{ix}));
  end;
  goodix = find(ismember(yrjds{1},jds));
  yrmat(1,1:length(goodix)) = yrdat{1}(goodix);
  for ix = 1:numel(yrdts)
    goodix = find(ismember(yrjds{ix},jds));
    yrmat(ix,:) = yrdat{ix}(goodix);
  end;

  yrmat = yrmat';

  % Do analysis of variance (and show BOXPLOT) for all good years
  [P,ANOVATAB,STATS] = anova1(yrmat, cellstr(num2str(meandts')), 'on');
  titlename(['Boxplot: Annual means']);
  maximize_graph;
  % print('-dpng', fullfile(figspath,['all-asam-boxplot-' ttltag '.png']));
  % ylim(sstlims);
  P,

  figure;
  [CMPS,MNS,fh,NMS] = multcompare(STATS);
  maximize_graph(fh);
  % xlim(sstlims);
  titlename([stn.station_name ' (Click station to test)']);
  % print('-dpng', fullfile(figspath,['all-asam-multcompare-' ttltag '.png']));
  [NMS num2cell(MNS)],

  figure;
  hold on;
  envelope = [ nanmin(yrmat(:)) nanmax(yrmat(:)) ];
  envelope(3) = envelope(2) - envelope(1);
  plot(jds, [ nanmin(yrmat) ; nanmean(yrmat); nanmax(yrmat) ; ...
              envelope(2) + envelope(3) + (nanmax(yrmat) - nanmin(yrmat)) ]);
  hold off;
  maximize_graph;
  datetick3;
  legend('Min','Mean','Max','Range');
  titlename(sprintf('%s inter-annual ranges: %d to %d ',...
                    stn.station_name,uyrs(1),uyrs(end)));
  % print('-dpng', fullfile(figspath,['all-asam-ranges-' ttltag '.png']));

return;
