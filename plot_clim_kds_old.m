1;

stnms={};
kds={};

stnms{end+1}='FWYF1';
% Original result as of Nov 2012
kds{end+1}=          [0.025,0.300,137];
% % Result as of Apr 2013
% kds{end+1}=          [0.025,0.225,102];

stnms{end+1}='MLRF1';
kds{end+1}=          [0.050,0.350, 91];
%kds{end+1}=          [.035,.250, 69];

stnms{end+1}='LONF1';
kds{end+1}=          [0.675,1.250, 10];
%kds{end+1}=          [.475,1.275, 45];

stnms{end+1}='SMKF1';
kds{end+1}=          [0.100,0.600, 70];
%kds{end+1}=          [0.066,0.450, 69];

stnms{end+1}='LOOE1';
kds{end+1}=          [0.010,0.400,113];
%kds{end+1}=          [0.025,0.200, 80];

stnms{end+1}='SANF1';
kds{end+1}=        [0.010,0.300, 81];
%kds{end+1}=        [0.015,0.150, 67];

stnms{end+1}='DRYF1';
kds{end+1}=        [0.250,0.750,342];
%kds{end+1}=        [0.150,0.500,354];


jds=[0:365];

clear kdclims
for ix=1:length(kds)
  kdclims(ix,1:length(jds)) = build_clim_opt(kds{ix},'Kd',jds);
end;
clear ix


doPlot=true;
doSave=true;
% doPlot=false;
% doSave=false;

if (doPlot)
  fmg;
  plot(jds,kdclims);
  legend(stnms, 'Location','North');
  titlename('Diffuse attenuation coefficient - climatology by site');
  xlabel('Year-Day');
  ylabel('K_d (m^-^1)');
  xlim([0,366]);
  ylim([0,1.4]);
  print('-dtiff',fullfile(get_thesis_path('../figs'),'Kd_by_site_OLD.tif'));
end;
if (doSave)
  save(fullfile(get_thesis_path('../data'),'Kd_by_site_OLD.mat'),'stnms','kds','jds','kdclims');
end;
