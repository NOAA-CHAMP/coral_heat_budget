function stn = phd_ch3_parinsoldose(stn,RAPFX)
%function stn = phd_ch3_parinsoldose(stn,RAPFX)
%
% Calculate PAR and insolation above surface and near bottom, from various
% sources at MLRF1: in situ data; in situ with climatological attenuation;
% combinations of ECMWF Reanalysis-Interim (ERAI), in situ, and climatology.
%
% Last Saved Time-stamp: <Thu 2013-06-27 18:32:59 Eastern Daylight Time gramer>

  doPrint = false;
  %doPrint = true;

  disp(['Print figures? ',num2str(doPrint)]);

  % adjust_reanalysis = true;
  adjust_reanalysis = false;
  disp(['Adjust reanalysis? ',num2str(adjust_reanalysis)]);
  if ( adjust_reanalysis )
    adjstr = 'Adj ';
  else
    adjstr = '';
  end;

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

  %% IN SITU DATA

  if ( ~exist('stn','var') || ~isfield(stn,'station_name') || ~strcmpi(stn.station_name,'mlrf1') )
    stn = []; clear stn;
    stn = get_station_from_station_name('mlrf1');
    stn = load_all_ndbc_data(stn);
    stn = load_station_data(stn);

    % Parcel out only PAR measurements taken between 18 hours and 7 days
    % immediately following a known station (and BIC) cleaning date.
    stn = load_cleaning_dates(stn);
    sfkeepix = [];
    shkeepix = [];
    for ix=1:length(stn.cleaning_date)
      dtdif = stn.bic_surf_par.date - stn.cleaning_date(ix);
      % sfkeepix = union(sfkeepix,find( (18/24) < dtdif & dtdif < 7 ));
      % Try one full day to two weeks after a cleaning instead!
      sfkeepix = union(sfkeepix,find( 1 < dtdif & dtdif < 14 ));

      dtdif = stn.bic_shallow_par.date - stn.cleaning_date(ix);
      % shkeepix = union(shkeepix,find( (18/24) < dtdif & dtdif < 7 ));
      % Try one full day to two weeks after a cleaning instead!
      shkeepix = union(shkeepix,find( 1 < dtdif & dtdif < 14 ));
    end;
    stn.bic_surf_par_clean = subset_ts(stn.bic_surf_par,sfkeepix);
    stn.bic_shallow_par_clean = subset_ts(stn.bic_shallow_par,shkeepix);

    stn = station_tmd_tide(stn);

    % stn = get_ngdc_bathy_station(stn);
    % %stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld,tufld,tvfld);
    % stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
    % % Experiment - compare with underwater PAR measurements at ~1m depth
    % stn = station_mean_tide_height(stn,mhfld,1,hfld);
    % Assume mean depth just the site itself: we are more interested in the
    % depth experienced by the BENTHOS, i.e., the actual site depth at MLRF1
    stn = station_mean_tide_height(stn,mhfld,stn.depth,hfld);

    stn.opts = [];
% %% PRIOR MLRF1
%     stn.opts.kd = [0.050,0.350, 91];
%% ACTUAL MLRF1
    stn.opts.kd = [0.035,0.250, 70];
% %%%%??? DEBUG EXPERIMENT
%     stn.opts.kd = [0.20,0.70, 91];
%     stn.opts.kd = [0.20,0.50, 91];
%     stn.opts.kd = [0.30,0.50,288];
%     stn.opts.kd = [0.05,0.30,288];
% %% BEST match w/1m post-cleaning measurements
%     stn.opts.kd = [0.05,0.20,288];
%     stn.opts.kd = [0.050,0.350,288];
  end;

  stn.bic_surf_par_daily_dose = par_dose(stn.bic_surf_par);
  if (0)
    fmg; boxplot_ts(stn.bic_surf_par_daily_dose,'month');
    titlename([upper(stn.station_name),' daily BIC sfc. PAR [mol/m^2/day]']);
  end;
  stn = station_par_to_insol(stn,'bic_surf_par','bic_surf_insol');
  stn.bic_surf_insol_daily_dose = par_dose(stn.bic_surf_insol);
  if (1-1)
    fmg; boxplot_ts(stn.bic_surf_insol_daily_dose,'month');
    % ylim([0,60]);
    ylim([0,30]);
    titlename([upper(stn.station_name),' daily BIC above-water insolation [MJ/m^2]',phd_ch3_parinsoldose_title_str(stn.bic_surf_insol_daily_dose)]);
    if ( doPrint )
      print('-dtiff',fullfile(get_thesis_path('../figs'),'mlrf1-daily-bic_surf_insol.tif'));
    end;
  end;


  stn.bic_surf_par_clean_daily_dose = par_dose(stn.bic_surf_par_clean);
  if (0)
    fmg; boxplot_ts(stn.bic_surf_par_clean_daily_dose,'month');
    titlename([upper(stn.station_name),' daily BIC clean sfc. PAR [mol/m^2/day]']);
  end;
  stn = station_par_to_insol(stn,'bic_surf_par_clean','bic_surf_insol_clean');
  stn.bic_surf_insol_clean_daily_dose = par_dose(stn.bic_surf_insol_clean);
  if (1-1)
    fmg; boxplot_ts(stn.bic_surf_insol_clean_daily_dose,'month');
    ylim([0,60]);
    titlename([upper(stn.station_name),' daily BIC clean sfc. insol. [MJ/m^2]']);
  end;


  stn.bic_shallow_par_daily_dose = par_dose(stn.bic_shallow_par);
  if (0)
    fmg; boxplot_ts(stn.bic_shallow_par_daily_dose,'month');
    titlename([upper(stn.station_name),' daily BIC btm. PAR [mol/m^2/day]']);
  end;
  stn = station_par_to_insol(stn,'bic_shallow_par','bic_shallow_insol');
  stn.bic_shallow_insol_daily_dose = par_dose(stn.bic_shallow_insol);
  if (1-1)
    fmg; boxplot_ts(stn.bic_shallow_insol_daily_dose,'month');
    % ylim([0,10]);
    ylim([0,60]);
    titlename([upper(stn.station_name),' daily BIC btm. insol. [MJ/m^2]']);
  end;


  stn.bic_shallow_par_clean_daily_dose = par_dose(stn.bic_shallow_par_clean);
  if (0)
    fmg; boxplot_ts(stn.bic_shallow_par_clean_daily_dose,'month');
    titlename([upper(stn.station_name),' daily BIC clean btm. PAR [mol/m^2/day]']);
  end;
  stn = station_par_to_insol(stn,'bic_shallow_par_clean','bic_shallow_insol_clean');
  stn.bic_shallow_insol_clean_daily_dose = par_dose(stn.bic_shallow_insol_clean);
  if (1-1)
    fmg; boxplot_ts(stn.bic_shallow_insol_clean_daily_dose,'month');
    % ylim([0,10]);
    ylim([0,60]);
    titlename([upper(stn.station_name),' daily BIC clean btm. insol. [MJ/m^2]']);
  end;


  stn.bic_climopt_kd.date = stn.bic_surf_insol.date;
  stn.bic_climopt_kd.data = build_clim_opt([0.05,0.25,137],'Kd',stn.bic_surf_insol.date);

  stn.bic_climopt_insol.date = stn.bic_surf_insol.date;
  stn.bic_climopt_insol.data = stn.bic_surf_insol.data.*exp(2.7.*stn.bic_climopt_kd.data);

  stn.bic_climopt_par.date = stn.bic_surf_par.date;
  stn.bic_climopt_par.data = stn.bic_surf_par.data.*exp(2.7.*stn.bic_climopt_kd.data);

  stn.bic_climopt_par_daily_dose = par_dose(stn.bic_climopt_par);
  if (0)
    fmg; boxplot_ts(stn.bic_climopt_par_daily_dose,'month');
    titlename([upper(stn.station_name),' climatological btm. PAR [mol/m^2/day]']);
  end;


  %% REANALYSIS

  disp(['Comparisons with ',RAPFX]);

  switch (RAPFX),
   case 'erai',
    stn = get_erai_station(stn);
    stn = adjust_erai_station(stn);
   case 'ncep',
    stn = get_ncep_station(stn,'narr');
   otherwise,
    error('Unknown reanalysis prefix %s',RAPFX);
  end;


  %% Q_S_W^I

  if (1-1+0)
    scatter_fit_ts(stn.(dsrfld),stn.bic_surf_insol,[],[],strrep(dsrfld,'_','\_'),'BIC Q_S_W^I',[],[],true);
    scatter_fit_ts(stn.(dsrfld),stn.bic_surf_insol_clean,[],[],strrep(dsrfld,'_','\_'),'Cleaned BIC Q_S_W^I',[],[],true);
  end;
  if (1-1+0)
    scatter_fit_ts(stn.(dsrfld),stn.bic_surf_insol,@ts_boreal_warm,@ts_florida_midday,[adjstr,upper(RAPFX),' Q_S_W^I'],'BIC Q_S_W^I',[],[],true);
    scatter_fit_ts(anomalize_ts(stn.(dsrfld)),anomalize_ts(stn.bic_surf_insol),@ts_boreal_warm,@ts_florida_midday,[adjstr,upper(RAPFX),' Q_S_W^I anomalized'],'BIC Q_S_W^I anomalized',[],[],true);
    stn.([dsrfld,'_anom'])=ts_op(stn.(dsrfld),nanmean(stn.(dsrfld).data),'-');
    stn.bic_surf_insol_anom=ts_op(stn.bic_surf_insol,nanmean(stn.bic_surf_insol.data),'-');
    scatter_fit_ts(stn.([dsrfld,'_anom']),stn.bic_surf_insol_anom,@ts_boreal_warm,@ts_florida_midday,[adjstr,upper(RAPFX),' Q_S_W^I anom.'],'BIC Q_S_W^I anom.',[],[],true);
  end;

  stn.([dsrfld,'_daily_dose']) = par_dose(stn.(dsrfld));

  if (1)
    scatter_fit_ts(stn.([dsrfld,'_daily_dose']),stn.bic_surf_insol_daily_dose,[],[],[adjstr,upper(RAPFX),' Q_S_W^I dose'],'BIC Q_S_W^I dose',[],[],true);
    axis([0,27,0,27]);
    if ( doPrint )
      print('-dtiff',fullfile(get_thesis_path('../figs'),['mlrf1-',[dsrfld,'_daily_dose'],'-scatter-bic_surf_insol_daily_dose','.tif']));
    end;

    [ig,ig,ax] = scatter_fit_ts_seasons(stn.([dsrfld,'_daily_dose']),stn.bic_surf_insol_daily_dose,[],[],[adjstr,upper(RAPFX),' Q_S_W^I dose'],'BIC Q_S_W^I dose',[],[],true);
    subplots_set('xlim',[0,27],'ylim',[0,27]);
    legend(ax(1),'Location','Best'); legend(ax(2),'Location','Best'); legend(ax(3),'Location','Best'); legend(ax(4),'Location','Best'); 
    if ( doPrint )
      print('-dtiff',fullfile(get_thesis_path('../figs'),['mlrf1-',[dsrfld,'_daily_dose'],'-scatter-seasons-bic_surf_insol_daily_dose','.tif']));
    end;
  end;


  if (1)
    scatter_fit_ts(stn.([dsrfld,'_daily_dose']),stn.bic_surf_insol_clean_daily_dose,[],[],[adjstr,upper(RAPFX),' Q_S_W^I dose'],'Clean BIC Q_S_W^I dose',[],[],true);
    axis([0,27,0,27]);
    if ( doPrint )
      print('-dtiff',fullfile(get_thesis_path('../figs'),['mlrf1-',[dsrfld,'_daily_dose'],'-scatter-bic_surf_insol_clean_daily_dose','.tif']));
    end;

    [ig,ig,ax] = scatter_fit_ts_seasons(stn.([dsrfld,'_daily_dose']),stn.bic_surf_insol_clean_daily_dose,[],[],[adjstr,upper(RAPFX),' Q_S_W^I dose'],'Clean BIC Q_S_W^I dose',[],[],true);
    subplots_set('xlim',[0,27],'ylim',[0,27]);
    legend(ax(1),'Location','Best'); legend(ax(2),'Location','Best'); legend(ax(3),'Location','Best'); legend(ax(4),'Location','Best'); 
    if ( doPrint )
      print('-dtiff',fullfile(get_thesis_path('../figs'),['mlrf1-',[dsrfld,'_daily_dose'],'-scatter-seasons-bic_surf_insol_clean_daily_dose','.tif']));
    end;
  end;


  %% Q_S_W^O

  disp('** Using bulk upward shortwave radiation **');
  stn = station_bulk_albedo(stn,albfld,Wfld,cfld);
  stn.(usrfld) = ts_op(stn.(dsrfld),stn.(albfld),'.*');
  stn.(usrfld).data = -stn.(usrfld).data;
  % Add down- and upward shortwave fluxes together
  stn.(srfld) = ts_op(stn.(dsrfld),stn.(usrfld),'+');


  %% Q_S_W

  stn.([srfld,'_daily_dose']) = par_dose(stn.(srfld));
  if (1-1)
    fmg; boxplot_ts(stn.([srfld,'_daily_dose']),'month');
    % ylim([0,60]);
    ylim([0,30]);
    titlename([upper(stn.station_name),' daily RA insolation [MJ/m^2]',phd_ch3_parinsoldose_title_str(stn.([srfld,'_daily_dose']))]);
    if ( doPrint )
      print('-dtiff',fullfile(get_thesis_path('../figs'),['mlrf1-daily-',lower(RAPFX),'_insol.tif']));
    end;
  end;

  %% \gammaQ_S_W

  stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,stn.opts,asrdiagfld);
  stn.b_dsrf = stn.(asrdiagfld).bottom_dsr;
  stn.b_dsrf_daily_dose = par_dose(stn.b_dsrf);
  if (1-1)
    fmg; boxplot_ts(stn.b_dsrf_daily_dose,'month');
    % ylim([0,10]);
    ylim([0,60]);
    titlename([upper(stn.station_name),' daily ',adjstr,upper(RAPFX),' btm. insol. [MJ/m^2]']);
  end;

  stn.bic_surf_usrf = ts_op(stn.bic_surf_insol,stn.(albfld),'.*');
  stn.bic_surf_usrf.data = -stn.bic_surf_usrf.data;
  % Add down- and upward shortwave fluxes together
  stn.bic_surf_srf = ts_op(stn.bic_surf_insol,stn.bic_surf_usrf,'+');

  % stn.b_bic_surf_dsrf = ts_op(stn.bic_surf_srf,stn.(asrdiagfld).tau,'*');
  stn.b_bic_surf_dsrf = ts_op(stn.bic_surf_srf,stn.(asrdiagfld).tau_scatter,'*');
  stn.b_bic_surf_dsrf.data = stn.(asrdiagfld).Ppen .* stn.b_bic_surf_dsrf.data;
  stn.b_bic_surf_dsrf_daily_dose = par_dose(stn.b_bic_surf_dsrf);
  if (1-1)
    fmg; boxplot_ts(stn.b_bic_surf_dsrf_daily_dose,'month');
    % ylim([0,10]);
    ylim([0,60]);
    titlename([upper(stn.station_name),' daily BIC/',adjstr,upper(RAPFX),' btm. insol. [MJ/m^2]']);
  end;

  if (1)
    fh=fmg;

    mos = unique(get_month(stn.bic_shallow_insol_clean_daily_dose.date));
    boxplot_ts(stn.bic_shallow_insol_clean_daily_dose,'month','allcolors','k','positions',mos);
    th=text(5.5,22,'Clean 1m measured'); set(th,'FontWeight','bold');

    mos = unique(get_month(stn.b_bic_surf_dsrf_daily_dose.date));
    boxplot_ts(stn.b_bic_surf_dsrf_daily_dose,'month','allcolors',[.5,.5,.5],'positions',mos);
    th=text(5.5,2,'Bottom estimated'); set(th,'FontWeight','bold','Color',[.5,.5,.5]);

    % ylim([0,10]);
    % ylim([0,60]);
    ylim([0,30]);

    titlename([upper(stn.station_name),' daily underwater insolation [MJ/m^2]',phd_ch3_parinsoldose_title_str(stn.bic_shallow_insol_clean_daily_dose)]);
    if ( doPrint )
      print('-dtiff',fullfile(get_thesis_path('../figs'),'mlrf1-daily-underwater-insols.tif'));

      %% TO GET FIGURE RIGHT, I HAD TO DO THIS NONSENSE BY HAND! MATLAB sucks...
      %set(gca,'fonts',18)
      %set(gca,'fonts',18)
      %cs=get(gco,'child'); for cix=1:numel(cs); o=cs(cix); x=get(o); if (isfield(x,'FontSize')); disp({cix,get(o,'FontSize')}); set(o,'FontSize',18); end; end;
      %cs=get(gco,'child'); for cix=1:numel(cs); o=cs(cix); x=get(o); if (isfield(x,'FontSize')); disp({cix,get(o,'FontSize')}); set(o,'FontSize',18); end; end;
      %saveas(fh,fullfile(get_thesis_path('../figs'),'mlrf1-daily-underwater-insols-SAVE-AS.tif'),'tiff');
    end;
  end;

  if (1)
    fh=fmg;

    mos = unique(get_month(stn.bic_shallow_insol_daily_dose.date));
    boxplot_ts(stn.bic_shallow_insol_daily_dose,'month','allcolors','k','positions',mos);
    th=text(5.5,22,'All 1m measured'); set(th,'FontWeight','bold');

    mos = unique(get_month(stn.b_bic_surf_dsrf_daily_dose.date));
    boxplot_ts(stn.b_bic_surf_dsrf_daily_dose,'month','allcolors',[.5,.5,.5],'positions',mos);
    th=text(5.5,2,'Bottom estimated'); set(th,'FontWeight','bold','Color',[.5,.5,.5]);

    ylim([0,30]);
    titlename([upper(stn.station_name),' daily underwater insolation [MJ/m^2]',phd_ch3_parinsoldose_title_str(stn.bic_shallow_insol_daily_dose)]);
  end;

  if (1-1)
    scatter_fit_ts(stn.b_bic_surf_dsrf_daily_dose,stn.bic_shallow_insol_clean_daily_dose,[],[],'Sfc. BIC \Sigma_1_dQ_b_S_W^I','Clean u/w BIC',[],[],true);
  end;

return;


%%%%%%%%%%%%%%%%%%%%
%% PRIVATE FUNCTIONS

function str = phd_ch3_parinsoldose_title_str(ts)
  yr = unique(get_year(ts.date));
  str = [' ',num2str(yr(1)),'-',num2str(yr(end)),' (df=',num2str(numel(ts.data)),')'];
return;
