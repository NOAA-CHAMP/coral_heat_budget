function stn = phd_ch2_parinsoldose(stn)
%function stn = phd_ch2_parinsoldose(stn)
%
% Calculate PAR and insolation above surface and near bottom, from various
% sources at MLRF1: in situ data; in situ with climatological attenuation;
% combinations of ECMWF Reanalysis-Interim (ERAI), in situ, and climatology.
%
% Last Saved Time-stamp: <Fri 2012-12-07 11:16:24 Eastern Standard Time gramer>

  % doPrint = false;
  doPrint = true;


  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

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
%% ACTUAL MLRF1
    stn.opts.kd = [0.050,0.350, 91];
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
  if (1)
    fmg; boxplot_ts(stn.bic_surf_insol_daily_dose,'month');
    % ylim([0,60]);
    ylim([0,30]);
    titlename([upper(stn.station_name),' daily BIC above-water insolation [MJ/m^2]',phd_ch2_parinsoldose_title_str(stn.bic_surf_insol_daily_dose)]);
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


  disp('Comparisons with ERAI');

  stn = get_erai_station(stn);
  stn = adjust_erai_station(stn);

  if (0)
    scatter_fit_ts(stn.erai_dsrf,stn.bic_surf_insol)
    scatter_fit_ts(stn.erai_dsrf,stn.bic_surf_insol,@ts_boreal_warm,@ts_florida_midday,'ERAI Q_S_W^I','BIC Q_S_W^I',[],[],true);
    scatter_fit_ts(anomalize_ts(stn.erai_dsrf),anomalize_ts(stn.bic_surf_insol),@ts_boreal_warm,@ts_florida_midday,'ERAI Q_S_W^I','BIC Q_S_W^I',[],[],true);
    stn.erai_dsrf_anom=ts_op(stn.erai_dsrf,nanmean(stn.erai_dsrf.data),'-');
    stn.bic_surf_insol_anom=ts_op(stn.bic_surf_insol,nanmean(stn.bic_surf_insol.data),'-');
    scatter_fit_ts(stn.erai_dsrf_anom,stn.bic_surf_insol_anom,@ts_boreal_warm,@ts_florida_midday,'ERAI Q_S_W^I anom.','BIC Q_S_W^I anom.',[],[],true);
  end;

  disp('** Using bulk upward shortwave radiation **');
  stn = station_bulk_albedo(stn,albfld,Wfld,cfld);
  stn.(usrfld) = ts_op(stn.(dsrfld),stn.(albfld),'.*');
  stn.(usrfld).data = -stn.(usrfld).data;
  % Add down- and upward shortwave fluxes together
  stn.(srfld) = ts_op(stn.(dsrfld),stn.(usrfld),'+');

  stn.([srfld,'_daily_dose']) = par_dose(stn.(srfld));
  if (1)
    fmg; boxplot_ts(stn.([srfld,'_daily_dose']),'month');
    % ylim([0,60]);
    ylim([0,30]);
    titlename([upper(stn.station_name),' daily ERAI insolation [MJ/m^2]',phd_ch2_parinsoldose_title_str(stn.([srfld,'_daily_dose']))]);
    if ( doPrint )
      print('-dtiff',fullfile(get_thesis_path('../figs'),'mlrf1-daily-erai_insol.tif'));
    end;
  end;

  stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,stn.opts,asrdiagfld);
  stn.b_dsrf = stn.(asrdiagfld).bottom_dsr;
  stn.b_dsrf_daily_dose = par_dose(stn.b_dsrf);
  if (1-1)
    fmg; boxplot_ts(stn.b_dsrf_daily_dose,'month');
    % ylim([0,10]);
    ylim([0,60]);
    titlename([upper(stn.station_name),' daily ERAI btm. insol. [MJ/m^2]']);
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
    titlename([upper(stn.station_name),' daily BIC/ERAI btm. insol. [MJ/m^2]']);
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

    titlename([upper(stn.station_name),' daily underwater insolation [MJ/m^2]',phd_ch2_parinsoldose_title_str(stn.bic_shallow_insol_clean_daily_dose)]);
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


  if (1-1)
    scatter_fit_ts(stn.b_bic_surf_dsrf_daily_dose,stn.bic_shallow_insol_clean_daily_dose,[],[],'Sfc. BIC \Sigma_1_dQ_b_S_W^I','Clean u/w BIC',[],[],true);
  end;

return;


%%%%%%%%%%%%%%%%%%%%
%% PRIVATE FUNCTIONS

function str = phd_ch2_parinsoldose_title_str(ts)
  yr = unique(get_year(ts.date));
  str = [' ',num2str(yr(1)),'-',num2str(yr(end)),' (df=',num2str(numel(ts.data)),')'];
return;
