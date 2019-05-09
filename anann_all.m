function anann_all(stn,fld)
%function anann_all(stn,fld)
%
% Analyze (with ANOVA1 and MULTCOMPARE, qv.) the annual data distributions of
% multi-year time series STN.(FLD). Only analyze data from years with at
% least 6000 of a possible 8760 hourly samples, but otherwise use all
% available data (unlike ANANN, qv., which uses only matching year-days).
%
% Last Saved Time-stamp: <Sat 2010-03-20 22:26:52 Eastern Daylight Time gramer>

  dat = stn.(fld).data;
  [yr,mo,dy] = datevec(stn.(fld).date);

  uyrs = unique(yr);

  badyrs = [];
  for yrix = 1:length(uyrs)
    ix = find(yr == uyrs(yrix));
    % Only do annual mean comparison on years that aren't too screwed up
    if ( length(ix) < 6000 && max(diff(stn.(fld).date(ix))) < 45 )
      badyrs = [badyrs(:)' uyrs(yrix)];
    end;
  end;

  dat = dat(~ismember(yr,badyrs));
  yr = yr(~ismember(yr,badyrs));

  yrstr = num2str(yr);

  [P,ANOVATAB,STATS] = anova1(dat, yrstr, 'on');
  titlename([stn.station_name ' boxplot: Annual means']);
  maximize_graph;
  % print('-dpng', fullfile(figspath,['all-boxplot-' ttltag '.png']));
  P,

  figure;
  [CMPS,MNS,fh,NMS] = multcompare(STATS);
  maximize_graph(fh);
  % xlim(sstlims);
  titlename([stn.station_name ' (Click station to test)']);
  % print('-dpng', fullfile(figspath,['all-multcompare-' ttltag '.png']));
  [NMS num2cell(MNS)],

return;
