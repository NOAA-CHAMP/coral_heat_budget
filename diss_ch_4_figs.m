1;


doPrint = true;

stnm = 'mlrf1';

if ( ~exist('stn','var') || ~isfield(stn,'qfld') )
  stn=[]; clear stn;
  load(fullfile(get_thesis_path('../data'),['sensitivity_analysis_',stnm,'.mat']));
  stn = s.([stnm,'_sensitivity_control']);
  s=[]; clear s;
end;
if ( ~isfield(stn,'ndbc_air_t') )
  y=load_all_ndbc_data([],stnm);
  stn.ndbc_air_t=y.ndbc_air_t;
  stn.ndbc_wind1_dir=y.ndbc_wind1_dir;
  stn.ndbc_wind1_speed=y.ndbc_wind1_speed;
  y=[]; clear y;
end;

clear s;

s.ts = 'ndbc_sea_t';
 s.ts_a = 'T_s';
s.ta = 'ndbc_air_t';
 s.ta_a = 'T_a';
s.ws = 'ndbc_wind1_speed';
 s.ws_a = 'U_1_0';
s.wu = 'ndbc_wind1_u';
 s.wu_a = 'U_1_0^x';
s.wv = 'ndbc_wind1_v';
 s.wv_a = 'U_1_0^y';

s.tsavg = [s.ts,'_3_d_avg'];
 s.tsavg_a = ['\mu_3_d',s.ts_a];
s.taavg = [s.ta,'_3_d_avg'];
 s.taavg_a = ['\mu_3_d',s.ta_a];
s.wsavg = [s.ws,'_3_d_avg'];
 s.wsavg_a = ['\mu_3_d',s.ws_a];

s.tsvar = [s.ts,'_1_d_dev_3_d_avg'];
 s.tsvar_a = ['\mu_3_d\sigma_1_d',s.ts_a];
s.tavar = [s.ta,'_1_d_dev_3_d_avg'];
 s.tavar_a = ['\mu_3_d\sigma_1_d',s.ta_a];
% s.wsvar = [s.ws,'_1_d_dev_3_d_avg'];
%  s.wsvar_a = ['\mu_3_d\sigma_1_d',s.ws_a];
% s.wsvar = [s.wu,'_3_d_dev_0_d_asof_add_',s.wv,'_3_d_dev'];
%  s.wsvar_a = ['\sigma_3_d',s.wu_a,'+\sigma_3_d',s.wv_a];
s.wsvar = [s.wu,'_7_d_dev_0_d_asof_sum_',s.wv,'_7_d_dev'];
 s.wsvar_a = ['\sigma_7_d',s.wu_a,'+\sigma_7_d',s.wv_a];

s.tsanom = [s.ts,'_3_d_avg_0_d_asof_diff_',s.ts,'_1_d_min'];
 s.tsanom_a = ['anom^m^i^n_3_d',s.ts_a];

s.q  = stn.qfld;
 s.q_a = ['Q_0/\rhoC_ph'];
s.bq = stn.bqfld;
 s.bq_a = ['[Q_0+Q_b]/\rhoC_ph'];
s.dt = stn.dtfld;
 s.dt_a = ['u\bullet\nablaT_s+K_h\nabla^2T_s+',s.bq];
s.hc = stn.hcfld;
 s.hc_a = ['\partial_tT_s'];


s.qsum  = [s.q,'_3_d_sum'];
 s.qsum_a = ['\Sigma_3_d',s.q_a];
s.bqsum = [s.bq,'_3_d_sum'];
 s.bqsum_a = ['\Sigma_3_d',s.bq_a];
s.dtsum = [s.dt,'_3_d_sum'];
 s.dtsum_a = ['\Sigma_3_d',s.dt_a];
s.hcsum = [s.hc,'_3_d_sum'];
 s.hcsum_a = ['\Sigma_3_d',s.hc_a];

stn = verify_variable(stn,{s.tsavg,s.taavg,s.wsavg,s.tsvar,s.tavar,s.wsvar,...
                    s.tsanom,s.qsum,s.bqsum,s.dtsum,s.hcsum});


%% NOTE: Taking absolute value of SUM here
s.aqsum  = [s.qsum,'_abs'];
 s.aqsum_a = ['|',s.qsum_a,'|'];
s.abqsum = [s.bqsum,'_abs'];
 s.abqsum_a = ['|',s.bqsum_a,'|'];
s.adtsum = [s.dtsum,'_abs'];
 s.adtsum_a = ['|',s.dtsum_a,'|'];
s.ahcsum = [s.hcsum,'_abs'];
 s.ahcsum_a = ['|',s.hcsum_a,'|'];

stn.(s.aqsum)  = ts_fun(stn.(s.qsum),@abs);
stn.(s.abqsum) = ts_fun(stn.(s.bqsum),@abs);
stn.(s.adtsum) = ts_fun(stn.(s.dtsum),@abs);
stn.(s.ahcsum) = ts_fun(stn.(s.hcsum),@abs);


%% NOTE: Taking sum of ABSOLUTE VALUES here
s.qa  = [s.q,'_abs'];
 s.qa_a = ['|',s.q,'|'];
s.bqa = [s.bq,'_abs'];
 s.bqa_a = ['|',s.bq,'|'];
s.dta = [s.dt,'_abs'];
 s.dta_a = ['|',s.dt,'|'];
s.hca = [s.hc,'_abs'];
 s.hca_a = ['|',s.hc,'|'];

stn.(s.qa)  = ts_fun(stn.(s.q),@abs);
stn.(s.bqa) = ts_fun(stn.(s.bq),@abs);
stn.(s.dta) = ts_fun(stn.(s.dt),@abs);
stn.(s.hca) = ts_fun(stn.(s.hc),@abs);

s.qasum  = [s.qa,'_3_d_sum'];
 s.qasum_a = ['\Sigma_3_d',s.qa_a];
s.bqasum = [s.bqa,'_3_d_sum'];
 s.bqasum_a = ['\Sigma_3_d',s.bqa_a];
s.dtasum = [s.dta,'_3_d_sum'];
 s.dtasum_a = ['\Sigma_3_d',s.dta_a];
s.hcasum = [s.hca,'_3_d_sum'];
 s.hcasum_a = ['\Sigma_3_d',s.hca_a];

stn = verify_variable(stn,{s.qasum,s.bqasum,s.dtasum,s.hcasum});



if (0)
  multiplot_station(stn,{s.ts,s.tsvar,s.ahcsum,s.aqsum,},...
                    [upper(stn.station_name),' Sea Temperature and Variability Metrics'],...
                    'Date',{s.ts_a,s.tsvar_a,s.ahcsum_a,s.aqsum_a,},...
                    [],{'default',[0,1],[0,1],[0,2]},10,{'k-','b-','g-','r-','c-'});
  if ( doPrint )
    print('-dtiff',fullfile(get_thesis_path('../DISS'),[lower(stn.station_name),'_',mfilename,'_',s.hcsum,'-TS-ALL.tif']));
  end;
end;

if (0)
  multiplot_station(stn,{s.ts,s.tsvar,s.ahcsum,s.aqsum,},...
                    [upper(stn.station_name),' Sea Temperature and Variability Metrics'],...
                    'Date',{s.ts_a,s.tsvar_a,s.ahcsum_a,s.aqsum_a,},...
                    [],{'default',[0,1],[0,1],[0,2]},10,{'k-','b-','g-','r-','c-'});
  xlim(datenum(2006,[5,9],[2,29]));
  datetick('x',2,'keeplimits');
  % % Arrow between flux and metric
  %tah=annotation('textarrow',[0.2904,0.2904],[0.4140,0.5665]);
  % Arrow between metric and sea temperature
  tah=annotation('textarrow',[0.2904,0.2904],[0.6650,0.7525]);
  set(tah,'String',{'Eddy Passage'},'LineWidth',2,'Color','r','FontName','TimesNewRoman','FontSize',12);
  if ( doPrint )
    print('-dtiff',fullfile(get_thesis_path('../DISS'),[lower(stn.station_name),'_',mfilename,'_',s.hcsum,'-TS-2006-eddy.tif']));
  end;
end;

if (0)
  multiplot_station(stn,{s.tsvar,s.ahcsum,s.aqsum,s.tavar,s.wsvar,},...
                    [upper(stn.station_name),' Sea Temperature and Variability Metrics'],...
                    'Date',{s.tsvar_a,s.ahcsum_a,s.aqsum_a,s.tavar_a,s.wsvar_a,},...
                    [],{[0,1],[0,1],[0,1],[0,1],[0,1]},10,{'k-','b-','g-','r-','c-','m-'});
  if ( doPrint )
    print('-dtiff',fullfile(get_thesis_path('../DISS'),[lower(stn.station_name),'_',mfilename,'_',s.hcsum,'.tif']));
  end;
end;

if (0)
  multiplot_station(stn,{s.tsvar,s.ahcsum,s.aqsum,s.tavar,s.wsvar,},...
                    [upper(stn.station_name),' Sea Temperature and Variability Metrics'],...
                    'Date',{s.tsvar_a,s.ahcsum_a,s.aqsum_a,s.tavar_a,s.wsvar_a,},...
                    [],{[0,1],[0,1],[0,1],[0,1],[0,1]},10,{'k-','b-','g-','r-','c-','m-'});
  xlim(datenum(2006,[5,9],[2,29]));
  datetick('x',2,'keeplimits');
  % tah=annotation('textarrow',[0.2904,0.2904],[0.6603,0.8128]);
  tah=annotation('textarrow',[0.2904,0.2904],[0.6086,0.7611]);
  set(tah,'String',{'Eddy Passage'},'LineWidth',2,'Color','r','FontName','TimesNewRoman','FontSize',12);
  if ( doPrint )
    print('-dtiff',fullfile(get_thesis_path('../DISS'),[lower(stn.station_name),'_',mfilename,'_',s.hcsum,'-2006-eddy.tif']));
  end;
end;


if (1)
  [ig,hl,ax,fh]=...
      multiplot_station(stn,{s.ts,s.tsvar,s.ahcsum,s.aqsum,s.tavar,s.wsvar,},...
                        [upper(stn.station_name),' Sea Temperature and Variability Metrics'],...
                        'Date',{s.ts_a,s.tsvar_a,s.ahcsum_a,s.aqsum_a,s.tavar_a,s.wsvar_a,},...
                        [],{[25,32],[0,1],[0,1],[0,2],[0,4],[0,20]},10,{'k-','b-','g-','r-','c-','m-'});
  ig=[]; clear ig
  xlim(datenum(2006,[5,9],[2,29]));
  datetick('x',2,'keeplimits');

  set(hl{4},'Color',[0.0,0.5,0.0]);
  set(hl{2},'Color',[0.5,0.5,0.5]);
  set([hl{:}],'LineWidth',1.5);
  set(ax,'FontSize',14);
  xlh=get(ax,'YLabel'); set([xlh{:}],'FontSize',16);

  % % Arrow between flux and metric
  %tah=annotation('textarrow',[0.2904,0.2904],[0.4140,0.5665]);
  % Arrow between metric and sea temperature
  tah=annotation('textarrow',[0.2904,0.2904],[0.6650,0.7525]);
  set(tah,'String',{'Eddy Passage'},'LineWidth',2,'Color','r','FontName','TimesNewRoman','FontSize',12);
  if ( doPrint )
    print('-dtiff',fullfile(get_thesis_path('../DISS'),[lower(stn.station_name),'_',mfilename,'_',s.hcsum,'-TS-2006-eddy-TALL.tif']));
  end;
end;

clear tah
