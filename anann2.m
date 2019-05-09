function anann2(stn,fld)

  dat = stn.(fld).data;
  [yr,mo,dy] = datevec(stn.(fld).date);

  dat = dat(yr ~= 2004 & yr ~= 2008);
  yr = yr(yr ~= 2004 & yr ~= 2008);

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
