function dat = test_ncep(dataset,yr,mo,fld,stnm,stfld)

  datapath = get_thesis_path('../data');

  fname = fullfile(datapath, sprintf('ncep_%s_%04d_%02d.mat', dataset, yr, mo));

  load(fname);

  [dat.(stnm).lon,dat.(stnm).lat,dat.(stnm).depth] = get_station_coords(stnm);


  ts = 1:12;

  rowix = 1:size(dat.(fld),2);
  colix = 1:size(dat.(fld),3);

  x = dat.(fld)(ts,rowix,colix);

  cmin = nanmin(x(:));
  cmax = nanmax(x(:));

  plotrows = ceil(sqrt(numel(ts)));
  plotcols = floor(sqrt(numel(ts)));

  figure;
  for ix=ts
    subplot(plotrows,plotcols,ix);
    pcolor(dat.lon(rowix,colix),dat.lat(rowix,colix),squeeze(dat.(fld)(ix,rowix,colix)));
    caxis([cmin cmax]);
    title(num2str(ix));
    set_pcolor_cursor;
    hold on;
    plot(dat.(stnm).lon,dat.(stnm).lat,'sb');
    hold off;
  end;
  maxigraph;

  % figure;
  % colorbar;
  % caxis([cmin cmax]);

  figure;
  plot(dat.(stnm).(stfld).date(ts), dat.(stnm).(stfld).data(ts));
  maxigraph;
  ylim([cmin cmax]);
  titlename(strrep(sprintf('%s.%s', stnm,stfld),'_','\_'));
  datetick3;

return;
