1;

  datapath = get_ecoforecasts_path('data');

  sanf1 = load_bb_data(fullfile(datapath,'sanf1.bb'),[],'sanf1');
  sanf1 = merge_station_data(sanf1,sanf1);
  smkf1 = load_bb_data(fullfile(datapath,'smkf1.bb'),[],'smkf1');
  smkf1 = merge_station_data(smkf1,smkf1);
  mlrf1 = load_bb_data(fullfile(datapath,'mlrf1.bb'),[],'mlrf1');
  mlrf1 = merge_station_data(mlrf1,mlrf1);
  fwyf1 = load_bb_data(fullfile(datapath,'fwyf1.bb'),[],'fwyf1');
  fwyf1 = merge_station_data(fwyf1,fwyf1);

  per = 'season';
  fmg;
  h1=boxplot_ts(sanf1.amodis_chlor_a,per,'allcol','m','mean',1);
  h2=boxplot_ts(smkf1.amodis_chlor_a,per,'allcol','k','mean',1);
  h3=boxplot_ts(mlrf1.amodis_chlor_a,per,'allcol','b','mean',1);
  h4=boxplot_ts(fwyf1.amodis_chlor_a,per,'allcol','r','mean',1);
  legend([h1(1),h2(1),h3(1),h4(1)],'SANF1','SMKF1','MLRF1','FWYF1');
  ylim([0,3]);
  set(gca,'xticklabel',{'JFM','AMJ','JAS','OND'});

  per=@(dts)(ceil(get_month(dts)./2));
  fmg;
  h1=boxplot_ts(sanf1.amodis_chlor_a,per,'allcol','m','mean',1);
  h2=boxplot_ts(smkf1.amodis_chlor_a,per,'allcol','k','mean',1);
  h3=boxplot_ts(mlrf1.amodis_chlor_a,per,'allcol','b','mean',1);
  h4=boxplot_ts(fwyf1.amodis_chlor_a,per,'allcol','r','mean',1);
  legend([h1(1),h2(1),h3(1),h4(1)],'SANF1','SMKF1','MLRF1','FWYF1');
  ylim([0,3]);
  set(gca,'xticklabel',{'Jan-Feb','Mar-Apr','May-Jun','Jul-Aug','Sep-Oct','Nov-Dec'});
