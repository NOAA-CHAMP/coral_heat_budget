function res = regress_temp_budget_terms(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,doPlot,saveFile)
%function res = regress_temp_budget_terms(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,doPlot,saveFile)
%
% Regress heat budget terms vs. variability in DAILY MEAN sea temperature.
% If DOPLOT (DEFAULT: false), plots ROBUSTFIT linear regressions; displays a
% report using DISP of simple R^2, RMSE, Scatter Index, and of R^2, RMSE, SI,
% bias, and slope from regression; if SAVEFILE (DEFAULT: false), also output
% to the CSV file FULLFILE(DATAPATH,[STNM,'-',MNM,'-',hcdTdt,'.csv']).
%
% See STATION_HEAT_BUDGET for the meaning of the other arguments.
%
% Last Saved Time-stamp: <Tue 2012-07-31 17:01:41  lew.gramer>


  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = false;
  end;
  if ( ~exist('saveFile','var') || isempty(saveFile) )
    saveFile = false;
  end;

  res = [];

  if ( strcmpi(ISPFX,'erai') )
    warning('Special ERAI options');
    begyr = 1996;
    adjust_waves = false;
    adjust_reanalysis = false;
    reanalysis_shortwave = true;
    reanalysis_longwave = true;
  end;

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;


  datapath = get_thesis_path('../data');

  mnm = lower(mfilename);
  MNM = upper(mfilename);

  stnm = lower(stn.station_name);
  STNM = upper(stn.station_name);

  cumfun = @floor; minN=24; cumstr='_1_d';
  %cumfun = @get_jday; minN=23; cumstr='_1_d_,_c_l_i_m';
  %cumfun = @get_week; minN=23*6; cumstr='_1_w_,_c_l_i_m';
  %cumfun = @get_yearweek; minN=24*7; cumstr='_1_w';
  %cumfun = @get_month; minN=23*25; cumstr='_1_m_,_c_l_i_m';
  %cumfun = @get_yearmonth; minN=24*28; cumstr='_m_o_n_t_h';

  disp([MNM,' Intersecting']);
  [raw_t,raw_qsw,raw_aqsw,raw_qlw,raw_qlh,raw_qsh,raw_sq0,raw_bq0,raw_bdt,raw_hcdt] = ...
      intersect_tses(stn.(sfld),...
                     stn.([srfld,'_term']),...
                     stn.([asrfld,'_term']),...
                     stn.([lrfld,'_term']),...
                     stn.([qlhfld,'_term']),...
                     stn.([qshfld,'_term']),...
                     stn.([sq0fld,'_term']),...
                     stn.([bq0fld,'_term']),...
                     stn.(bdTfld),...
                     stn.(hcdTdt));


  disp([MNM,' Accumulating']);
  [t.data,t.date] = grp_ts(raw_t.data,raw_t.date,cumfun,@nanmean,minN);
  dt.date = t.date(2:end); dt.data = diff(t.data);

  [qsw.data,qsw.date] = grp_ts(raw_qsw.data,raw_qsw.date,cumfun,@nansum,minN);
  [aqsw.data,aqsw.date] = grp_ts(raw_aqsw.data,raw_aqsw.date,cumfun,@nansum,minN);
  [qlw.data,qlw.date] = grp_ts(raw_qlw.data,raw_qlw.date,cumfun,@nansum,minN);
  [qlh.data,qlh.date] = grp_ts(raw_qlh.data,raw_qlh.date,cumfun,@nansum,minN);
  [qsh.data,qsh.date] = grp_ts(raw_qsh.data,raw_qsh.date,cumfun,@nansum,minN);

  [sq0.data,sq0.date] = grp_ts(raw_sq0.data,raw_sq0.date,cumfun,@nansum,minN);
  [bq0.data,bq0.date] = grp_ts(raw_bq0.data,raw_bq0.date,cumfun,@nansum,minN);
  [bdt.data,bdt.date] = grp_ts(raw_bdt.data,raw_bdt.date,cumfun,@nansum,minN);
  [hcdt.data,hcdt.date] = grp_ts(raw_hcdt.data,raw_hcdt.date,cumfun,@nansum,minN);


  disp([MNM,' Regressing']);

  [t,dt,qsw,aqsw,qlw,qlh,qsh,sq0,bq0,bdt,hcdt] = ...
      intersect_tses(t,dt,qsw,aqsw,qlw,qlh,qsh,sq0,bq0,bdt,hcdt);

  N = numel(dt.data);
  % nrm = nanmean(abs(dt.data));
  nrm = nanvar(dt.data);
  mn = nanmean(dt.data);


  if ( doPlot )
    fh = [];
  else
    fh = 'none';
  end;

  ix=0;

  ix=ix+1; res.nm{ix}='Q_S_W/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qsw,dt,[],[],[STNM,' \Sigma',cumstr,res.nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  res.rR2(ix)=Stats.regress_stats(1); res.rRMSE(ix)=Stats.s; res.rSI(ix)=Stats.s./nrm; res.rbias(ix)=B(1); res.rslope(ix)=B(2);
  [res.SSE(ix),res.R2(ix),res.RMSE(ix),res.SI(ix)] = r_t_b_t_calc(qsw,dt);

  ix=ix+1; res.nm{ix}='\gammaQ_S_W/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(aqsw,dt,[],[],[STNM,' \Sigma',cumstr,res.nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  res.rR2(ix)=Stats.regress_stats(1); res.rRMSE(ix)=Stats.s; res.rSI(ix)=Stats.s./nrm; res.rbias(ix)=B(1); res.rslope(ix)=B(2);
  [res.SSE(ix),res.R2(ix),res.RMSE(ix),res.SI(ix)] = r_t_b_t_calc(aqsw,dt);

  ix=ix+1; res.nm{ix}='Q_L_W/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qlw,dt,[],[],[STNM,' \Sigma',cumstr,res.nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  res.rR2(ix)=Stats.regress_stats(1); res.rRMSE(ix)=Stats.s; res.rSI(ix)=Stats.s./nrm; res.rbias(ix)=B(1); res.rslope(ix)=B(2);
  [res.SSE(ix),res.R2(ix),res.RMSE(ix),res.SI(ix)] = r_t_b_t_calc(qlw,dt);

  ix=ix+1; res.nm{ix}='Q_L_H/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qlh,dt,[],[],[STNM,' \Sigma',cumstr,res.nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  res.rR2(ix)=Stats.regress_stats(1); res.rRMSE(ix)=Stats.s; res.rSI(ix)=Stats.s./nrm; res.rbias(ix)=B(1); res.rslope(ix)=B(2);
  [res.SSE(ix),res.R2(ix),res.RMSE(ix),res.SI(ix)] = r_t_b_t_calc(qlh,dt);

  ix=ix+1; res.nm{ix}='Q_S_H/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qsh,dt,[],[],[STNM,' \Sigma',cumstr,res.nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  res.rR2(ix)=Stats.regress_stats(1); res.rRMSE(ix)=Stats.s; res.rSI(ix)=Stats.s./nrm; res.rbias(ix)=B(1); res.rslope(ix)=B(2);
  [res.SSE(ix),res.R2(ix),res.RMSE(ix),res.SI(ix)] = r_t_b_t_calc(qsh,dt);

  ix=ix+1; res.nm{ix}='Q_0/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(sq0,dt,[],[],[STNM,' \Sigma',cumstr,res.nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  res.rR2(ix)=Stats.regress_stats(1); res.rRMSE(ix)=Stats.s; res.rSI(ix)=Stats.s./nrm; res.rbias(ix)=B(1); res.rslope(ix)=B(2);
  [res.SSE(ix),res.R2(ix),res.RMSE(ix),res.SI(ix)] = r_t_b_t_calc(sq0,dt);

  ix=ix+1; res.nm{ix}='(Q_0(\gamma)+Q_b)/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(bq0,dt,[],[],[STNM,' \Sigma',cumstr,res.nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  res.rR2(ix)=Stats.regress_stats(1); res.rRMSE(ix)=Stats.s; res.rSI(ix)=Stats.s./nrm; res.rbias(ix)=B(1); res.rslope(ix)=B(2);
  [res.SSE(ix),res.R2(ix),res.RMSE(ix),res.SI(ix)] = r_t_b_t_calc(bq0,dt);

  ix=ix+1; res.nm{ix}='(Q_0(\gamma)+Q_b)/\rhoC_ph+u_s_f_c^.\nablaT_k_m+K_h\nabla^2T_k_m';
  [B,Stats]=scatter_fit_ts(bdt,dt,[],[],[STNM,' \Sigma',cumstr,res.nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  res.rR2(ix)=Stats.regress_stats(1); res.rRMSE(ix)=Stats.s; res.rSI(ix)=Stats.s./nrm; res.rbias(ix)=B(1); res.rslope(ix)=B(2);
  [res.SSE(ix),res.R2(ix),res.RMSE(ix),res.SI(ix)] = r_t_b_t_calc(bdt,dt);

  ix=ix+1; res.nm{ix}='\partial_tT_s';
  [B,Stats]=scatter_fit_ts(hcdt,dt,[],[],[STNM,' \Sigma',cumstr,res.nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  res.rR2(ix)=Stats.regress_stats(1); res.rRMSE(ix)=Stats.s; res.rSI(ix)=Stats.s./nrm; res.rbias(ix)=B(1); res.rslope(ix)=B(2);
  [res.SSE(ix),res.R2(ix),res.RMSE(ix),res.SI(ix)] = r_t_b_t_calc(hcdt,dt);

  disp(MNM);
  disp('Daily_Temperature_Change,Mean,Variance');
  disp({sfld,mn,nrm});
  % fprintf(1,'Daily_Temperature_Change,Mean,Variance\n');
  % fprintf(1,'%s,%g,%g\n',sfld,mn,nrm);
  disp('R^2  RMSE   SI    Regr.: R^2    RMSE    SI    A    B');
  % fprintf(1,'%s,R2,RMSE,SI,rR2,rRMSE,rSI,rA,rB\n',MNM);
  for ix=1:numel(res.R2)
    disp(res.nm{ix});
    disp([res.R2(ix),res.RMSE(ix),res.SI(ix),NaN,res.rR2(ix),res.rRMSE(ix),res.rSI(ix),res.rbias(ix),res.rslope(ix)]);
    % fprintf(1,'%s,%g,%g,%g,R:%g,%g,%g,%g,%g\n',res.nm{ix},res.R2(ix),...
    %         res.RMSE(ix),res.SI(ix),res.rR2(ix),res.rRMSE(ix),res.rSI(ix),res.rbias(ix),res.rslope(ix));
  end;

  if ( saveFile )
    fname = fullfile(datapath,[stnm,'-',mnm,'-',hcdTdt,'.csv']);
    if ( exist(fname,'file') )
      disp(['Overwriting regression results on ',fname]);
    else
      disp(['Writing regression results to ',fname]);
    end;
    fid = fopen(fname,'w+');
    if ( fid<2 )
      warning('UNABLE TO WRITE "%s"',fname);
    else
      fprintf(fid,'Daily_Temperature_Change,Mean,Variance\n');
      fprintf(fid,'%s,%g,%g\n',sfld,mn,nrm);
      fprintf(fid,'%s,R2,RMSE,SI,rR2,rRMSE,rSI,rA,rB\n',MNM);
      for ix=1:numel(res.R2)
        fprintf(fid,'%s,%g,%g,%g,%g,%g,%g,%g,%g\n',res.nm{ix},res.R2(ix),...
                res.RMSE(ix),res.SI(ix),res.rR2(ix),res.rRMSE(ix),res.rSI(ix),res.rbias(ix),res.rslope(ix));
      end;
      fclose(fid);
    end;
  end;    


return;

%%%%%%%%%%%%%%%%%%%%%
%% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%

function [sse,r2,rmse,si] = r_t_b_t_calc(trm,dt)
%function [sse,r2,rmse,si] = r_t_b_t_calc(trm,dt)
  N = numel(dt.data);
  nrm = nanvar(dt.data);
  mn = nanmean(dt.data);

  sse = sum((trm.data-dt.data).^2);
  r2 = 1-(sse./(N.*nrm));
  % % May not apply
  % ssr = sum((trm.data-mn).^2);
  % r2 = ssr./(N.*nrm);

  rmse = sqrt(sse./N);
  si = rmse./nrm;
return;
