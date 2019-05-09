1;

diary off
more off;

figspath = get_thesis_path('../figs');
doPrint = false;

stnm='lonf1';

%%%%???DEBUG
diary(fullfile(get_thesis_path('../src'),[stnm,'-sensitivity.log']));


% @ROUND: Daily mean and climatology
cf=@floor;

RAPFX='erai';
KMPFX='avhrr_weekly';
ISPFX='ndbc';
TIDEPFX='tpxo_tide';
WAVEPFX='erai';
subs={};
%subs={'sfld','hourly_misst_sst', 'comparisonsfld','ndbc_sea_t'};


disp('*************************************************');
disp([upper(stnm),' SENSITIVITY ANALYSIS: ',datestr(now)]);
disp('*************************************************');

s.station_name = stnm;

%function stn = optimize_station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,opts)

disp('******************** CONTROL ********************');
fname=[stnm,'-sensitivity-','control'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

disp('******************** Albedo ********************');
fname=[stnm,'-sensitivity-','reanalysis_shortwave'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,{subs{:},'reanalysis_shortwave',true}),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

for alb = [0.04,0.11];
  fname=[stnm,'-sensitivity-','albedo-',num2str(alb)],
  fld=regexprep(fname,'[.-]','_');
  s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('albedo',alb)),cf);
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;
end;

disp('******************** Attenuation Kd ********************');
%% WAS:           [0.675,1.250, 10] ... % Best for ERAI waves adv=[0,1,91],K=[0,2,91]
%for kd = {[0.675,1.25,55],[0.675,1.35, 10],0.96};
% IS:          [.475,1.275, 45] ...     % Best with hc_scaling=SU,hc_warming_factor=0.66
for kd = {[.475,1.275, 91],[.475,1.375, 45],0.88};
  fname=[stnm,'-sensitivity-','kd',num2str(kd{:},'-%g')],
  fld=regexprep(fname,'[.-]','_');
  s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_light_kds',{kd})),cf);
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;
end;

disp('******************** Bottom reflectivity Ab ********************');
for Ab = [0.24,0.0,1.0];
  fname=[stnm,'-sensitivity-','Ab',num2str(Ab,'-%g')],
  fld=regexprep(fname,'[.-]','_');
  s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('bottom_reflectance',Ab)),cf);
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;
end;

disp('******************** Water emissivity epsilon_water ********************');
for epsw = [0.95,0.98];
  fname=[stnm,'-sensitivity-','eps',num2str(epsw,'-%g')],
  fld=regexprep(fname,'[.-]','_');
  s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('epsilon_water',epsw)),cf);
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;
end;

disp('******************** (No) Warm-layer adjustment ********************');
fname=[stnm,'-sensitivity-','do-warm','-true'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_do_warm_adj',{true})),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

disp('******************** Meteorology ********************');
fname=[stnm,'-sensitivity-','meteorology','-erai-erai-ndbc_sea_t'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,'erai',TIDEPFX,WAVEPFX,{subs{:},'sfld','ndbc_sea_t'}),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

fname=[stnm,'-sensitivity-','meteorology','-narr-ndbc-ndbc_sea_t'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,'ncep',KMPFX,ISPFX,TIDEPFX,WAVEPFX,{subs{:},'adjust_reanalysis',false}),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

fname=[stnm,'-sensitivity-','meteorology','-narr-narr-ndbc_sea_t'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,'ncep',KMPFX,'ncep',TIDEPFX,WAVEPFX,{subs{:},'sfld','ndbc_sea_t','adjust_reanalysis',false}),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

disp('******************** Wave model ********************');
fname=[stnm,'-sensitivity-','waves','-ww3'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,'ww3',subs),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

% disp('*************************************************');
% disp(['RESTARTING ',upper(stnm),' SENSITIVITY ANALYSIS: ',datestr(now)]);
% disp('*************************************************');

fname=[stnm,'-sensitivity-','waves','-ndbc'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,'ndbc',subs),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

disp('******************** Surface current ********************');
fname=[stnm,'-sensitivity-','currents','-1-pct-of-wind'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('ocean_current_wind_multiplier',0.01)),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

fname=[stnm,'-sensitivity-','currents','-2-pct-of-wind'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('ocean_current_wind_multiplier',0.02)),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

disp('******************** Advection factor Fq ********************');
for advfac = {0.0,1.0,[0.0,1.0, 91]};
  fname=[stnm,'-sensitivity-','advfac',num2str(advfac{:},'-%g')],
  fld=regexprep(fname,'[.-]','_');
  s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_advection_factors',{advfac})),cf);
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;
end;

% disp('*************************************************');
% disp(['FINISHING ',upper(stnm),' SENSITIVITY ANALYSIS: ',datestr(now)]);
% disp('*************************************************');

disp('******************** Ktheta ********************');
for kth = {20,0,'ndbc_wind1_speed'};
  if ( ~ischar(kth{:}) )
    fname=[stnm,'-sensitivity-','ktheta',num2str(kth{:},'-%g')],
  else
    fname=[stnm,'-sensitivity-','ktheta','-',kth{:}],
    stn = load_all_ndbc_data([],stnm);
    kth = { { stn.ndbc_wind1_speed,@(W)(min(20,((W./35).^2).*20)) } };
  end;
  fld=regexprep(fname,'[.-]','_');
  s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_k_thetas',{kth})),cf);
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;
  stn=[]; clear stn
end;

% disp('*************************************************');
% disp(['RESTARTING ',upper(stnm),' SENSITIVITY ANALYSIS: ',datestr(now)]);
% disp('*************************************************');

disp('******************** Advection+Diffusion ********************');
fname=[stnm,'-sensitivity-','adv-dif','-gom_hycom'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,'gom_hycom',ISPFX,TIDEPFX,WAVEPFX,subs,struct('grid_interp_method','nearest')),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

fname=[stnm,'-sensitivity-','adv-dif','-fkeys_hycom'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,'fkeys_hycom',ISPFX,TIDEPFX,WAVEPFX,subs,struct('grid_interp_method','nearest')),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

 % Use HYCOM sea temperature gradients, but empirical currents
fname=[stnm,'-sensitivity-','adv-dif','-gom_hycom-empirical'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,'gom_hycom',ISPFX,TIDEPFX,WAVEPFX,subs,struct('grid_interp_method','nearest','km_scale_advection',false)),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

disp('******************** Horizontal Convection ********************');
fname=[stnm,'-sensitivity-','hc_scaling','-SS'],
fld=regexprep(fname,'[.-]','_');
disp('HC=SS (steady momentum and thermal)');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('hc_scaling','SS')),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

fname=[stnm,'-sensitivity-','hc_scaling','-UU'],
fld=regexprep(fname,'[.-]','_');
disp('HC=UU (unsteady momentum and thermal)');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('hc_scaling','UU')),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

fname=[stnm,'-sensitivity-','hc_scaling','-US'],
fld=regexprep(fname,'[.-]','_');
disp('HC=US (unsteady momentum, steady thermal)');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('hc_scaling','US')),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

for R = [0.00,0.11,0.20];
  fname=[stnm,'-sensitivity-','hc_R',num2str(R,'-%g')],
  fld=regexprep(fname,'[.-]','_');
  s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('hc_R',(1.0-R))),cf);
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;
end;

disp('******************** Salinity ********************');
for sal = [35.2,36.2];
  fname=[stnm,'-sensitivity-','default_salinity',num2str(sal,'-%g')],
  fld=regexprep(fname,'[.-]','_');
  s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('default_salinity',sal)),cf);
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;
end;

disp('******************** Water Depth ********************');
fname=[stnm,'-sensitivity-','depth','-tmd_tide'],
fld=regexprep(fname,'[.-]','_');
disp('TMD_TIDE (GoM solution)');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,'tmd_tide',WAVEPFX,subs),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

stn=get_station_from_station_name(stnm);
fname=[stnm,'-sensitivity-','depth','-fixed-sensor'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_mean_tide_height',stn.depth)),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;

fname=[stnm,'-sensitivity-','depth','-mean-sensor'],
fld=regexprep(fname,'[.-]','_');
s.(fld)=calc_heat_budget_rmse(optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_tide_height',stn.depth)),cf);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,[fname,'.tif'])); end;
stn=[]; clear stn;

% %   cbds    = get_opt(stn.opts,'override_benthic_cbd',cbds);
% %   doWarms = get_opt(stn.opts,'override_do_warm_adj',doWarms);
% %   kds     = get_opt(stn.opts,'override_light_kds',kds);
% %   advfacs = get_opt(stn.opts,'override_advection_factors',advfacs);
% %   kths    = get_opt(stn.opts,'override_k_thetas',kths);

% %   stn.opts.Ppen = get_opt(stn.opts,'override_Ppen',[]);
% %   qlh_adj = get_opt(stn.opts,'override_qlh_adj',qlh_adj);

% %   hcs = get_opt(stn.opts,'override_hc_scaling',stn.opts.hc_scaling);


%% Create a summary report and plot of all total heat-budget RMSEs
anhcrmsec;

save(fullfile(get_thesis_path('../src'),[stnm,'-sensitivity-results.mat']),'s');


diary off;
more on;
