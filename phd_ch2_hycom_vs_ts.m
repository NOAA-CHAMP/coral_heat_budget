1;

tic,

figspath = get_thesis_path('../figs');
cstnms={'fwyf1','mlrf1','lonf1','smkf1','sanf1'};
for cstnm=cstnms;

  stnm=lower(cstnm{:});
  STNM=upper(cstnm{:});
  disp(STNM);
  stn=get_station_from_station_name(stnm);
  stn=load_all_ndbc_data(stn);
  stn=verify_variable(stn,{'ndbc_sea_t_1_d_avg','ndbc_sea_t_7_d_avg'});


  % Gulf of Mexico HYCOM
  interpMethod = [];
  switch (stnm),
    % case 'lkwf1',   interpMethod={@nanmean,2,2};
    % case 'dryf1',   interpMethod={@nanmean,2,2};
  end;
  stn=get_gom_hycom(stn,[],[],[],[],interpMethod);

  scatter_fit_ts_seasons(stn.gom_hycom_seatemp,stn.ndbc_sea_t,[],[],'GOM HYCOM',[STNM,' T_s'],[],[],true);
  subplots_set('xlim',[14,34],'ylim',[14,34]);
  print('-dtiff',fullfile(figspath,[stnm,'-gom_hycom_seatemp-scatter-ndbc_sea_t.tif']));

  scatter_fit_ts_seasons(stn.gom_hycom_seatemp,stn.ndbc_sea_t_1_d_avg,[],[],'GOM HYCOM',[STNM,' \mu_1_dT_s'],[],[],true);
  subplots_set('xlim',[14,34],'ylim',[14,34]);
  print('-dtiff',fullfile(figspath,[stnm,'-gom_hycom_seatemp-scatter-ndbc_sea_t_1_d_avg.tif']));


  % Florida Keys HYCOM
  interpMethod = [];
  switch (stnm),
    % These cause an error in GET_FKEYS_HYCOM - code needs updating slightly
   case 'lkwf1',   interpMethod={@nanmean,2,2};
   case 'dryf1',   interpMethod={@nanmean,2,2};
  end;
  stn=get_fkeys_hycom(stn,[],[],[],[],interpMethod);

  scatter_fit_ts_seasons(stn.fkeys_hycom_seatemp,stn.ndbc_sea_t,[],[],'FKEYS HYCOM',[STNM,' T_s'],[],[],true);
  subplots_set('xlim',[14,34],'ylim',[14,34]);
  print('-dtiff',fullfile(figspath,[stnm,'-fkeys_hycom_seatemp-scatter-ndbc_sea_t.tif']));

  scatter_fit_ts_seasons(stn.fkeys_hycom_seatemp,stn.ndbc_sea_t_1_d_avg,[],[],'FKEYS HYCOM',[STNM,' \mu_1_dT_s'],[],[],true);
  subplots_set('xlim',[14,34],'ylim',[14,34]);
  print('-dtiff',fullfile(figspath,[stnm,'-fkeys_hycom_seatemp-scatter-ndbc_sea_t_1_d_avg.tif']));


  stn=[]; clear stn;
end;

clear figspath cstnms cstnm stnm ans
toc,
