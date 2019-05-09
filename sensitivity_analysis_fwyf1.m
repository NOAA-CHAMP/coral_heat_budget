1;

doDiary=false;
doPrint=false;

diary off
more off;

figspath = get_thesis_path('../figs');

stnm='fwyf1';

RAPFX='erai';
KMPFX='avhrr_weekly';
ISPFX='ndbc';
TIDEPFX='tpxo_tide';
WAVEPFX='erai';
subs={};
%subs={'sfld','hourly_misst_sst', 'comparisonsfld','ndbc_sea_t'};


if (doDiary); diary([stnm,'-sensitivity.log']); end;


disp('*************************************************');
disp([upper(stnm),' SENSITIVITY ANALYSIS: ',datestr(now)]);
disp('*************************************************');

%function stn = optimize_station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,opts)

disp('******************** CONTROL ********************');
fname=[stnm,'-sensitivity-','control','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

disp('******************** Albedo ********************');
fname=[stnm,'-sensitivity-','reanalysis_shortwave','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,{subs{:},'reanalysis_shortwave',true});
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

for alb = [0.04,0.11];
  fname=[stnm,'-sensitivity-','albedo-',num2str(alb),'.tif'],
  optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('albedo',alb));
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
end;

disp('******************** Attenuation Kd ********************');
for kd = {[0.025,0.300, 91],[0.025,0.250,137],0.15};
  fname=[stnm,'-sensitivity-','kd',num2str(kd{:},'-%g'),'.tif'],
  optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_light_kds',{kd}));
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
end;

%%%%%???DEBUG
%keyboard;

disp('******************** Bottom reflectivity Ab ********************');
for Ab = [0.17,0.0,1.0];
  fname=[stnm,'-sensitivity-','Ab',num2str(Ab,'-%g'),'.tif'],
  optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('bottom_reflectance',Ab));
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
end;

% disp('*************************************************');
% disp(['RESTARTING ',upper(stnm),' SENSITIVITY ANALYSIS: ',datestr(now)]);
% disp('*************************************************');

disp('******************** Water emissivity epsilon_water ********************');
for epsw = [0.95,0.98];
  fname=[stnm,'-sensitivity-','eps',num2str(epsw,'-%g'),'.tif'],
  optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('epsilon_water',epsw));
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
end;

disp('******************** (No) Warm-layer adjustment ********************');
fname=[stnm,'-sensitivity-','do-warm','-false','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_do_warm_adj',{false}));
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

disp('******************** Meteorology ********************');
fname=[stnm,'-sensitivity-','meteorology','-erai-erai-ndbc_sea_t','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,'erai',TIDEPFX,WAVEPFX,{subs{:},'sfld','ndbc_sea_t'});
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

fname=[stnm,'-sensitivity-','meteorology','-narr-ndbc-ndbc_sea_t','.tif'],
optimize_station_heat_budget(stnm,'ncep',KMPFX,ISPFX,TIDEPFX,WAVEPFX,{subs{:},'adjust_reanalysis',false});
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

fname=[stnm,'-sensitivity-','meteorology','-narr-narr-ndbc_sea_t','.tif'],
optimize_station_heat_budget(stnm,'ncep',KMPFX,'ncep',TIDEPFX,WAVEPFX,{subs{:},'sfld','ndbc_sea_t','adjust_reanalysis',false});
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

disp('******************** Wave model ********************');
fname=[stnm,'-sensitivity-','waves','-ww3','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,'ww3',subs);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

% disp('*************************************************');
% disp(['RESTARTING ',upper(stnm),' SENSITIVITY ANALYSIS: ',datestr(now)]);
% disp('*************************************************');

fname=[stnm,'-sensitivity-','waves','-ndbc','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,'ndbc',subs);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

disp('******************** Surface current ********************');
fname=[stnm,'-sensitivity-','currents','-1-pct-of-wind','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('ocean_current_wind_multiplier',0.01));
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

fname=[stnm,'-sensitivity-','currents','-2-pct-of-wind','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('ocean_current_wind_multiplier',0.02));
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

disp('******************** Advection factor Fq ********************');
for advfac = {0,1,[0.0,0.5, 91]};
  fname=[stnm,'-sensitivity-','advfac',num2str(advfac{:},'-%g'),'.tif'],
  optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_advection_factors',{advfac}));
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
end;

% disp('*************************************************');
% disp(['FINISHING ',upper(stnm),' SENSITIVITY ANALYSIS: ',datestr(now)]);
% disp('*************************************************');

disp('******************** Ktheta ********************');
for kth = {10,0,'ndbc_wind1_speed'};
  if ( ~ischar(kth{:}) )
    fname=[stnm,'-sensitivity-','ktheta',num2str(kth{:},'-%g'),'.tif'],
  else
    fname=[stnm,'-sensitivity-','ktheta','-',kth{:},'.tif'],
    stn = load_all_ndbc_data([],stnm);
    kth = { { stn.ndbc_wind1_speed,@(W)(min(10,((W./35).^2).*10)) } };
  end;
  optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_k_thetas',{kth}));
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
  stn=[]; clear stn
end;

% disp('*************************************************');
% disp(['RESTARTING ',upper(stnm),' SENSITIVITY ANALYSIS: ',datestr(now)]);
% disp('*************************************************');

disp('******************** Advection+Diffusion ********************');
fname=[stnm,'-sensitivity-','adv-dif','-gom_hycom','.tif'],
optimize_station_heat_budget(stnm,RAPFX,'gom_hycom',ISPFX,TIDEPFX,WAVEPFX,subs);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

fname=[stnm,'-sensitivity-','adv-dif','-fkeys_hycom','.tif'],
optimize_station_heat_budget(stnm,RAPFX,'fkeys_hycom',ISPFX,TIDEPFX,WAVEPFX,subs);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

 % Use HYCOM sea temperature gradients, but empirical currents
fname=[stnm,'-sensitivity-','adv-dif','-gom_hycom-empirical','.tif'],
optimize_station_heat_budget(stnm,RAPFX,'gom_hycom',ISPFX,TIDEPFX,WAVEPFX,subs,struct('km_scale_advection',false));
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

disp('******************** Horizontal Convection ********************');
fname=[stnm,'-sensitivity-','hc_scaling','-SS','.tif'],
disp('HC=SS (steady thermal)');
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('hc_scaling','SS'));
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

fname=[stnm,'-sensitivity-','hc_scaling','-UU','.tif'],
disp('HC=UU (unsteady momentum)');
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('hc_scaling','UU'));
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

fname=[stnm,'-sensitivity-','hc_scaling','-SU','.tif'],
disp('HC=SU (steady thermal, unsteady momentum)');
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('hc_scaling','SU'));
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

for R = [0.00,0.11,0.20];
  fname=[stnm,'-sensitivity-','hc_R',num2str(R,'-%g'),'.tif'],
  optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('hc_R',(1.0-R)));
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
end;

disp('******************** Salinity ********************');
for sal = [35.0,36.0];
  fname=[stnm,'-sensitivity-','default_salinity',num2str(sal,'-%g'),'.tif'],
  optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('default_salinity',sal));
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
end;

disp('******************** Water Depth ********************');
fname=[stnm,'-sensitivity-','depth','-tmd_tide','.tif'],
disp('TMD_TIDE (GoM solution)');
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,'tmd_tide',WAVEPFX,subs);
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

stn=get_station_from_station_name(stnm);
fname=[stnm,'-sensitivity-','depth','-fixed-sensor','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_mean_tide_height',stn.depth));
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

fname=[stnm,'-sensitivity-','depth','-mean-sensor','.tif'],
optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_tide_height',stn.depth));
xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
stn=[]; clear stn;

% %   cbds    = get_opt(stn.opts,'override_benthic_cbd',cbds);
% %   doWarms = get_opt(stn.opts,'override_do_warm_adj',doWarms);
% %   kds     = get_opt(stn.opts,'override_light_kds',kds);
% %   advfacs = get_opt(stn.opts,'override_advection_factors',advfacs);
% %   kths    = get_opt(stn.opts,'override_k_thetas',kths);

% %   stn.opts.Ppen = get_opt(stn.opts,'override_Ppen',[]);
% %   qlh_adj = get_opt(stn.opts,'override_qlh_adj',qlh_adj);

% %   hcs = get_opt(stn.opts,'override_hc_scaling',stn.opts.hc_scaling);

diary off;
more on;
