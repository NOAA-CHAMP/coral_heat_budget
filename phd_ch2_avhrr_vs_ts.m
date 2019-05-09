1;

figspath = get_thesis_path('../figs');
cstnms={'lkwf1','fwyf1','mlrf1','nfbf1','lonf1','smkf1','sanf1','dryf1','42003'};
for cstnm=cstnms;
  stnm=lower(cstnm{:});
  STNM=upper(cstnm{:});
  disp(STNM);
  tic,
    stn=get_station_from_station_name(stnm);
    stn=load_all_ndbc_data(stn);
    stn=verify_variable(stn,{'ndbc_sea_t_1_d_avg','ndbc_sea_t_7_d_avg'});

    interpMethod = []; npts=5;
	switch (stnm),
      case 'lkwf1',   interpMethod={@nanmean,2,2}; npts=3;
      case 'sanf1',   interpMethod='nearest';      npts=3;
      case 'dryf1',   interpMethod={@nanmean,2,2}; npts=3;
    end;
    stn=get_avhrr_weekly_field(stn,false,interpMethod,npts);

    scatter_fit_ts_seasons(stn.raw_avhrr_weekly_sst,stn.ndbc_sea_t,[],[],'Raw AVHRR Weekly',[STNM,' T_s'],[],[],true);
    subplots_set('xlim',[18,32],'ylim',[18,32]);
    print('-dtiff',fullfile(figspath,[stnm,'-raw_avhrr_weekly_sst-scatter-ndbc_sea_t.tif']));

    scatter_fit_ts_seasons(stn.raw_avhrr_weekly_sst,stn.ndbc_sea_t_1_d_avg,[],[],'Raw AVHRR Weekly',[STNM,' \mu_1_dT_s'],[],[],true);
    subplots_set('xlim',[18,32],'ylim',[18,32]);
    print('-dtiff',fullfile(figspath,[stnm,'-raw_avhrr_weekly_sst-scatter-ndbc_sea_t_1_d_avg.tif']));

    scatter_fit_ts_seasons(stn.raw_avhrr_weekly_sst,stn.ndbc_sea_t_7_d_avg,[],[],'Raw AVHRR Weekly',[STNM,' \mu_7_dT_s'],[],[],true);
    subplots_set('xlim',[18,32],'ylim',[18,32]);
    print('-dtiff',fullfile(figspath,[stnm,'-raw_avhrr_weekly_sst-scatter-ndbc_sea_t_7_d_avg.tif']));

    stn=[]; clear stn;
  toc,
end;
clear figspath cstnms cstnm stnm ans
