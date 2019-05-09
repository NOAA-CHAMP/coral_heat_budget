1;

if ( ~exist('s','var') )
  s=[];
end;

%for perfun=[@floor,@get_yearmonth]; for seasfilt=[@ts_isfinite,@ts_boreal_warm,@ts_boreal_cool]; for filterBadDates = [false,true];
for perfun=[@floor]; for seasfilt=[@ts_isfinite,@ts_boreal_warm,@ts_boreal_cool]; for filterBadDates = [true];

diary off;
more off;

switch (char(perfun))
 case 'floor';               cumfun = @nanmean;      minN = 24;      maxgaps=(25/24);
 case 'get_yearweek';        cumfun = @nanmean;      minN = 24*4.5;  maxgaps=8;
 case 'get_yearmonth';       cumfun = @nanmean;      minN = 24*20;   maxgaps=32;
 otherwise,                  error('Unsupported perfun');
end;

%seasfilt = [];
%seasfilt = @ts_boreal_warm;
%seasfilt = @ts_boreal_cool;

%filterBadDates = false;
%filterBadDates = true;

if ( filterBadDates ); filtstr = 'filtered'; else filtstr = 'unfiltered'; end;

fname = fullfile(get_thesis_path('../data'),[mfilename,'-looe1-',char(perfun),'-',char(seasfilt),'-',filtstr,'-results.log']);
diary(fname);

disp(mfilename);
timenow;
disp(char(perfun));
disp(char(seasfilt));
if ( filterBadDates )
  disp('WILL FILTER "BAD DATES"');
else
  disp('Dates unfiltered');
end;

sacs=0;
sars=0;
esacs=0;
esars=0;
swcs=0;
swrs=0;
eswcs=0;
eswrs=0;
sw2cs=0;
sw2rs=0;
esw2cs=0;
esw2rs=0;
cs=0;
rs=0;
windrs=0;
rhrs = 0;

snm='looe1';

for ctsnm={'mc_seatemp','ad_seatemp'};
  tsnm = ctsnm{:};

  disp('====================');
  disp(upper(tsnm));
  disp('====================');

  if ( ~isfield(s,snm) )
    s.(snm) = get_station_from_station_name(snm);
    s.(snm) = load_all_ndbc_data(s.(snm));
    s.(snm) = station_optimal_isobath_orientation(s.(snm));
    s.(snm) = verify_variable(s.(snm),'ndbc_wind1_u');
    s.(snm) = verify_variable(s.(snm),'ndbc_wind1_v');
    s.(snm) = station_reorient_vectors(s.(snm),s.(snm).isobath_orientation,...
                                       'ndbc_wind1_u','ndbc_wind1_v');

    x = get_looe1_microcat;
    s.(snm).mc_seatemp = interp_ts(x.microcat_seatemp);
    x=[]; clear x;

    x = get_looe1_adcp;
    s.(snm).ad_seatemp = x.adcp_seatemp;
    s.(snm).ad_x.date = x.adcp_x.date;
    s.(snm).ad_x.data = x.adcp_x.data;
    s.(snm).ad_l.date = x.adcp_l.date;
    s.(snm).ad_l.data = x.adcp_l.data;
    s.(snm).ad_speed.date = x.adcp_speed.date;
    s.(snm).ad_speed.data = x.adcp_speed.data;
    x=[]; clear x;

    s.(snm).ndbc_wind2.date = s.(snm).ndbc_wind1_xshore.date;
    s.(snm).ndbc_wind2.data = ((s.(snm).ndbc_wind1_xshore.data .^ 2) + (s.(snm).ndbc_wind1_lshore.data .^ 2));
    s.(snm).ndbc_wind.date = s.(snm).ndbc_wind2.date;
    s.(snm).ndbc_wind.data = sqrt(s.(snm).ndbc_wind2.data);
  end;

  [six,aix]=intersect_dates(s.(snm).(tsnm).date,s.(snm).ndbc_air_t.date);
  s.(snm).([tsnm,'_air_diff']).date=s.(snm).(tsnm).date(six);
  s.(snm).([tsnm,'_air_diff']).data=s.(snm).(tsnm).data(six)-s.(snm).ndbc_air_t.data(aix);

  if ( ~isfield(s.(snm),'erai_air_t') )
    x = get_erai_station(snm);
    disp('Adjusting ERAI radiation, waves');
    x = adjust_erai_station(x);
    x = adjust_erai_station_waves(x);
    s.(snm).erai_air_t = x.erai_air_t;
    s.(snm).erai_relhumid = x.erai_relhumid;
    s.(snm).erai_spechumid = x.erai_spechumid;
    s.(snm).erai_wind_speed = x.erai_wind_speed;
    s.(snm).erai_wind_dir = x.erai_wind_dir;
    s.(snm).erai_dsrf = x.erai_dsrf;
    s.(snm).erai_dlrf = x.erai_dlrf;
    s.(snm).erai_precip = x.erai_precip;
    s.(snm).erai_dsrf_adj = x.erai_dsrf_adj;
    s.(snm).erai_dlrf_adj = x.erai_dlrf_adj;
    s.(snm).erai_precip_adj = x.erai_precip_adj;
    s.(snm).erai_sigwavehgt = x.erai_sigwavehgt;
    s.(snm).erai_sigwavehgt_adj = x.erai_sigwavehgt_adj;
    x=[]; clear x;
    s.(snm) = verify_variable(s.(snm),'erai_wind_u');
    s.(snm) = verify_variable(s.(snm),'erai_wind_v');
    s.(snm) = station_reorient_vectors(s.(snm),s.(snm).isobath_orientation,...
                                       'erai_wind_u','erai_wind_v');

    s.(snm).erai_wind2.date = s.(snm).erai_wind_xshore.date;
    s.(snm).erai_wind2.data = ((s.(snm).erai_wind_xshore.data .^ 2) + (s.(snm).erai_wind_lshore.data .^ 2));
    s.(snm).erai_wind.date = s.(snm).erai_wind2.date;
    s.(snm).erai_wind.data = sqrt(s.(snm).erai_wind2.data);
  end;
  if ( ~isfield(s.(snm),'avhrr_weekly_sst') )
    switch (snm)
     case 'sanf1',
      s.(snm) = get_avhrr_weekly_field(s.(snm),true,'nearest',3,~filterBadDates);
     case 'dryf1',
      s.(snm) = get_avhrr_weekly_field(s.(snm),true,{'nanmean',2,2},3,~filterBadDates);
     otherwise,
      s.(snm) = get_avhrr_weekly_field(s.(snm),true,'linear',5,~filterBadDates);
    end;

    s.(snm) = rmfield(s.(snm), grepstruct(s.(snm),'_field'));
    s.(snm) = rmfield(s.(snm), grepstruct(s.(snm),'raw_'));

    s.(snm) = station_reorient_vectors(s.(snm),s.(snm).isobath_orientation,...
                                       'hourly_avhrr_weekly_sst_x','hourly_avhrr_weekly_sst_y');
  end;

  [s.(snm).ts_per.data,s.(snm).ts_per.date]=grp_ts(s.(snm).(tsnm).data,s.(snm).(tsnm).date,perfun,cumfun,minN);

  [s.(snm).sst_per.data,s.(snm).sst_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst.data,s.(snm).hourly_avhrr_weekly_sst.date,perfun,cumfun,minN);
  [s.(snm).sst_xs_per.data,s.(snm).sst_xs_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_xshore.data,s.(snm).hourly_avhrr_weekly_sst_xshore.date,perfun,cumfun,minN);
  [s.(snm).sst_ls_per.data,s.(snm).sst_ls_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_lshore.data,s.(snm).hourly_avhrr_weekly_sst_lshore.date,perfun,cumfun,minN);
  try,
    [s.(snm).sst_l_per.data,s.(snm).sst_l_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_l.data,s.(snm).hourly_avhrr_weekly_sst_l.date,perfun,cumfun,minN);
    [s.(snm).sst_dl_per.data,s.(snm).sst_dl_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_dl.data,s.(snm).hourly_avhrr_weekly_sst_dl.date,perfun,cumfun,minN);
  catch,
  end;

  [s.(snm).ta_per.data,s.(snm).ta_per.date]=grp_ts(s.(snm).ndbc_air_t.data,s.(snm).ndbc_air_t.date,perfun,cumfun,minN);
  [s.(snm).W_kts_per.data,s.(snm).W_kts_per.date]=grp_ts(s.(snm).ndbc_wind1_speed.data,s.(snm).ndbc_wind1_speed.date,perfun,cumfun,minN);
  [s.(snm).U_per.data,s.(snm).U_per.date]=grp_ts(s.(snm).ndbc_wind1_xshore.data,s.(snm).ndbc_wind1_xshore.date,perfun,cumfun,minN);
  [s.(snm).V_per.data,s.(snm).V_per.date]=grp_ts(s.(snm).ndbc_wind1_lshore.data,s.(snm).ndbc_wind1_lshore.date,perfun,cumfun,minN);
  [s.(snm).W_per.data,s.(snm).W_per.date]=grp_ts(s.(snm).ndbc_wind.data,s.(snm).ndbc_wind.date,perfun,cumfun,minN);
  [s.(snm).W2_per.data,s.(snm).W2_per.date]=grp_ts(s.(snm).ndbc_wind2.data,s.(snm).ndbc_wind2.date,perfun,cumfun,minN);


  [s.(snm).e_ta_per.data,s.(snm).e_ta_per.date]=grp_ts(s.(snm).erai_air_t.data,s.(snm).erai_air_t.date,perfun,cumfun,minN);
  [s.(snm).e_qa_per.data,s.(snm).e_qa_per.date]=grp_ts(s.(snm).erai_spechumid.data,s.(snm).erai_spechumid.date,perfun,cumfun,minN);
  [s.(snm).e_W_kts_per.data,s.(snm).e_W_kts_per.date]=grp_ts(s.(snm).erai_wind_speed.data,s.(snm).erai_wind_speed.date,perfun,cumfun,minN);
  [s.(snm).e_U_per.data,s.(snm).e_U_per.date]=grp_ts(s.(snm).erai_wind_xshore.data,s.(snm).erai_wind_xshore.date,perfun,cumfun,minN);
  [s.(snm).e_V_per.data,s.(snm).e_V_per.date]=grp_ts(s.(snm).erai_wind_lshore.data,s.(snm).erai_wind_lshore.date,perfun,cumfun,minN);
  [s.(snm).e_W_per.data,s.(snm).e_W_per.date]=grp_ts(s.(snm).erai_wind.data,s.(snm).erai_wind.date,perfun,cumfun,minN);
  [s.(snm).e_W2_per.data,s.(snm).e_W2_per.date]=grp_ts(s.(snm).erai_wind2.data,s.(snm).erai_wind2.date,perfun,cumfun,minN);

  % NOTE USE OF "@nansum" FOR INSOLATION!
  [s.(snm).e_qswi_per.data,s.(snm).e_qswi_per.date]=grp_ts(s.(snm).erai_dsrf_adj.data,s.(snm).erai_dsrf_adj.date,perfun,@nansum,minN);
  [s.(snm).e_qlwi_per.data,s.(snm).e_qlwi_per.date]=grp_ts(s.(snm).erai_dlrf_adj.data,s.(snm).erai_dlrf_adj.date,perfun,cumfun,minN);

  [s.(snm).e_hs_per.data,s.(snm).e_hs_per.date]=grp_ts(s.(snm).erai_sigwavehgt_adj.data,s.(snm).erai_sigwavehgt_adj.date,perfun,cumfun,minN);


  [s.(snm).ad_xs_per.data,s.(snm).ad_xs_per.date]=grp_ts(s.(snm).ad_x.data,s.(snm).ad_x.date,perfun,cumfun,minN);
  [s.(snm).ad_ls_per.data,s.(snm).ad_ls_per.date]=grp_ts(s.(snm).ad_l.data,s.(snm).ad_l.date,perfun,cumfun,minN);


  % Calculate one-day differences with gaps removed
  s.(snm).ts_per_dif.date=s.(snm).ts_per.date(2:end); s.(snm).ts_per_dif.data=diff(s.(snm).ts_per.data);
  s.(snm).ta_per_dif.date=s.(snm).ta_per.date(2:end); s.(snm).ta_per_dif.data=diff(s.(snm).ta_per.data);
  s.(snm).e_ta_per_dif.date=s.(snm).e_ta_per.date(2:end); s.(snm).e_ta_per_dif.data=diff(s.(snm).e_ta_per.data);
  s.(snm).e_qa_per_dif.date=s.(snm).e_qa_per.date(2:end); s.(snm).e_qa_per_dif.data=diff(s.(snm).e_qa_per.data);
  s.(snm).sst_per_dif.date=s.(snm).sst_per.date(2:end); s.(snm).sst_per_dif.data=diff(s.(snm).sst_per.data);
  % REMOVE GAPS...
  s.(snm) = filter_gaps(s.(snm),'ts_per','ts_per_dif',maxgaps);
  s.(snm) = filter_gaps(s.(snm),'ta_per','ta_per_dif',maxgaps);
  s.(snm) = filter_gaps(s.(snm),'e_ta_per','e_ta_per_dif',maxgaps);
  s.(snm) = filter_gaps(s.(snm),'e_qa_per','e_qa_per_dif',maxgaps);
  s.(snm) = filter_gaps(s.(snm),'sst_per','sst_per_dif',maxgaps);


  %% REMOVE ANOMALY???
  % s.(snm).U_per = anomalize_daily_mean_ts(s.(snm).U_per);
  % s.(snm).V_per = anomalize_daily_mean_ts(s.(snm).V_per);
  % s.(snm).W_per = anomalize_daily_mean_ts(s.(snm).W_per);
  % s.(snm).W2_per = anomalize_daily_mean_ts(s.(snm).W2_per);
  % s.(snm).e_hs_per = anomalize_daily_mean_ts(s.(snm).e_hs_per);
  % s.(snm).sst_xs_per = anomalize_daily_mean_ts(s.(snm).sst_xs_per);
  % s.(snm).sst_ls_per = anomalize_daily_mean_ts(s.(snm).sst_ls_per);
  % try,
  %   s.(snm).sst_l_per = anomalize_daily_mean_ts(s.(snm).sst_l_per);
  %   s.(snm).sst_dl_per = anomalize_daily_mean_ts(s.(snm).sst_dl_per);
  % catch,
  % end;
  % s.(snm).ts_per_dif = anomalize_daily_mean_ts(s.(snm).ts_per_dif);
  % s.(snm).ta_per_dif = anomalize_daily_mean_ts(s.(snm).ta_per_dif);
  % s.(snm).e_ta_per_dif = anomalize_daily_mean_ts(s.(snm).e_ta_per_dif);
  % s.(snm).e_qa_per_dif = anomalize_daily_mean_ts(s.(snm).e_qa_per_dif);
  % s.(snm).sst_per_dif = anomalize_daily_mean_ts(s.(snm).sst_per_dif);


  [c,r,p,ix]=cov_ts(s.(snm).ta_per_dif,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Dper Tair vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  sacs = sacs + c;
  sars = sars + r;

  [c,r,p,ix]=cov_ts(s.(snm).e_ta_per_dif,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' ERAI Dper Tair vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  esacs = esacs + c;
  esars = esars + r;


  if ( ismember(snm,{'lonf1','smkf1'}) )
    s.(snm) = station_dewp_to_relhumid(s.(snm),'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
    s.(snm) = station_relhumid_to_spechumid(s.(snm),'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
    [s.(snm).qa_per.data,s.(snm).qa_per.date]=grp_ts(s.(snm).ndbc_spechumid.data,s.(snm).ndbc_spechumid.date,perfun,cumfun,minN);
    s.(snm).qa_per_dif.date=s.(snm).qa_per.date(2:end); s.(snm).qa_per_dif.data=diff(s.(snm).qa_per.data);
    s.(snm) = filter_gaps(s.(snm),'qa_per','qa_per_dif',maxgaps);

    [c,r,p,ix]=cov_ts(s.(snm).qa_per_dif,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper qa vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  end;

  [c,r,p,ix]=cov_ts(s.(snm).e_qa_per_dif,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' ERAI Dper qa vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);


  [c,r,p,ix]=cov_ts(s.(snm).U_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Per U vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

  [c,r,p,ix]=cov_ts(s.(snm).V_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Per V vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

  [c,r,p,ix]=cov_ts(s.(snm).W_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Per Wind vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  swcs = swcs + c;
  swrs = swrs + r;

  [c,r,p,ix]=cov_ts(s.(snm).e_W_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' ERAI per Wind vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  eswcs = eswcs + c;
  eswrs = eswrs + r;

  [c,r,p,ix]=cov_ts(s.(snm).W2_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Per Wind2 vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  sw2cs = sw2cs + c;
  sw2rs = sw2rs + r;

  [c,r,p,ix]=cov_ts(s.(snm).e_W2_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' ERAI per Wind2 vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  esw2cs = esw2cs + c;
  esw2rs = esw2rs + r;


  if (0)
  end;

  [c,r,p,ix]=cov_ts(s.(snm).ndbc_air_t,s.(snm).(tsnm),seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Tair vs. Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  % sacs = sacs + c;
  % sars = sars + r;

  [c,r,p,ix]=cov_ts(s.(snm).ndbc_wind1_speed,s.(snm).(tsnm),seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Wind vs. Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  % swcs = swcs + c;
  % swrs = swrs + r;

  [c,r,p,ix]=cov_ts(s.(snm).ndbc_wind1_speed,s.(snm).([tsnm,'_air_diff']),seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Wind vs. Tsea-Tair: R=' num2str(r)]);
  cs = cs + c;
  rs = rs + r;

  ug = s.(snm).ndbc_air_t;
  ug.date(ug.data<1 | ~isfinite(ug.data)) = [];
  ug.data(ug.data<1 | ~isfinite(ug.data)) = [];
  ug.data = (6./ug.data);
  [c,r,p,ix]=cov_ts(s.(snm).ndbc_wind2,ug,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' U^2 + V^2 vs. Ug^2 ~ 6/Ta: C=' num2str(c)]);
  % cs = cs + c;


  if ( ismember(snm,{'lonf1','smkf1'}) )
    s.(snm).expcdew = s.(snm).ndbc_dew_t;
    s.(snm).expcdew.data = exp(0.07*s.(snm).expcdew.data);
    s.(snm).expmcta = s.(snm).ndbc_air_t;
    s.(snm).expmcta.data = exp(-0.07*s.(snm).expmcta.data);
    % [c,r,p,ix]=cov_ts(s.(snm).([tsnm,'_air_diff']),s.(snm).ndbc_wind1_speed,seasfilt);
    [c,r,p,ix]=cov_ts(s.(snm).expcdew,s.(snm).ndbc_wind1_speed,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' exp(C*Tdew) vs. exp(-C*Ta): R=' num2str(r) ' (',num2str(p),')']);
    windrs = windrs + r;

    expta = s.(snm).ndbc_air_t;
    expta.data = exp(0.06*expta.data);
    [c,r,p,ix]=cov_ts(s.(snm).ndbc_relhumid,expta,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' RH vs. exp(0.06*Ta): R=' num2str(r) ' (',num2str(p),')']);
    rhrs = rhrs + r;
  end;

  [c,r,p,ix]=cov_ts(s.(snm).e_qswi_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Adj ERAI QSWI vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

  [c,r,p,ix]=cov_ts(s.(snm).e_qlwi_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Adj ERAI QLWI vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

  [c,r,p,ix]=cov_ts(s.(snm).e_hs_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Adj ERAI Hs vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

  [c,r,p,ix]=cov_ts(s.(snm).sst_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' AVHRR SST vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

  [c,r,p,ix]=cov_ts(s.(snm).sst_per_dif,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

  [c,r,p,ix]=cov_ts(s.(snm).sst_xs_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST XS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

  [c,r,p,ix]=cov_ts(s.(snm).sst_ls_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST LS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

  try,
    [c,r,p,ix]=cov_ts(s.(snm).sst_l_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST L vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).sst_dl_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST DL vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  catch,
  end;

  [c,r,p,ix]=cov_ts(s.(snm).ad_xs_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Dper ADCP XS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
  [c,r,p,ix]=cov_ts(s.(snm).ad_ls_per,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' Dper ADCP LS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);


end;
sacs=sacs/6,
sars=sars/6,
esacs=esacs/6,
esars=esars/6,
swcs=swcs/6,
swrs=swrs/6,
eswcs=eswcs/6,
eswrs=eswrs/6,
sw2cs=sw2cs/6,
sw2rs=sw2rs/6,
esw2cs=esw2cs/6,
esw2rs=esw2rs/6,
cs=cs/6,
rs=rs/6,
windrs=windrs/2,
rhrs=rhrs/2,

clear c r six aix ug expcd expmct expta

timenow;

more on;

diary off;

end; end; end;
