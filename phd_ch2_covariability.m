1;

doDiary = false;
doCovs = false;
doMultiRegress = true;

if ( ~exist('s','var') )
  s=[];
end;

% stnms={'fwyf1','mlrf1','lonf1','smkf1','looe1','sanf1','dryf1'};
stnms={'mlrf1','lonf1'};

% %for perfun=[@floor,@get_yearmonth]; for seasfilt=[@ts_isfinite,@ts_boreal_warm,@ts_boreal_cool]; for filterBadDates = [false,true];
% for perfun=[@floor]; for cseasfilt={@ts_isfinite,@ts_boreal_warm,@ts_boreal_cool}; for filterBadDates = [true];
for perfun=[@floor]; for cseasfilt={@ts_isfinite}; for filterBadDates = [true];

seasfilt = cseasfilt{:};

diary off;
more off;

switch (char(perfun))
 case 'floor';               cumfun = @nanmean;      minN = 24;      maxgaps=(25/24);
 case 'get_yearweek';        cumfun = @nanmean;      minN = 24*4.5;  maxgaps=8;
 case 'get_yearmonth';       cumfun = @nanmean;      minN = 24*20;   maxgaps=32;
 otherwise,                  error('Unsupported perfun');
end;

%seasfilt = @ts_isfinite;
%seasfilt = @ts_boreal_warm;
%seasfilt = @ts_boreal_cool;

switch (char(seasfilt))
 case 'ts_isfinite';         per = 'year';
 case 'ts_boreal_warm';      per = 'warm';
 case 'ts_boreal_cool';      per = 'cool';
 otherwise,                  error('Unsupported seasfilt');
end;

%filterBadDates = false;
%filterBadDates = true;

if ( filterBadDates ); filtstr = 'filtered'; else filtstr = 'unfiltered'; end;

if ( doDiary )
  fname = fullfile(get_thesis_path('../data'),[mfilename,'-',char(perfun),'-',char(seasfilt),'-',filtstr,'-results.log']);
  diary(fname);
end;

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

for snmc=stnms
  snm=snmc{:};

  if ( ~isfield(s,snm) )
    s.(snm) = get_station_from_station_name(snm);
    s.(snm) = load_all_ndbc_data(s.(snm));
    s.(snm) = station_optimal_isobath_orientation(s.(snm));
    s.(snm) = verify_variable(s.(snm),'ndbc_wind1_u');
    s.(snm) = verify_variable(s.(snm),'ndbc_wind1_v');
    s.(snm) = station_reorient_vectors(s.(snm),s.(snm).isobath_orientation,...
                                       'ndbc_wind1_u','ndbc_wind1_v');

    [six,aix]=intersect_dates(s.(snm).ndbc_sea_t.date,s.(snm).ndbc_air_t.date);
    s.(snm).ndbc_sea_air_diff.date=s.(snm).ndbc_sea_t.date(six);
    s.(snm).ndbc_sea_air_diff.data=s.(snm).ndbc_sea_t.data(six)-s.(snm).ndbc_air_t.data(aix);

    s.(snm).ndbc_wind2.date = s.(snm).ndbc_wind1_xshore.date;
    s.(snm).ndbc_wind2.data = ((s.(snm).ndbc_wind1_xshore.data .^ 2) + (s.(snm).ndbc_wind1_lshore.data .^ 2));
    s.(snm).ndbc_wind.date = s.(snm).ndbc_wind2.date;
    s.(snm).ndbc_wind.data = sqrt(s.(snm).ndbc_wind2.data);
  end;
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

  [s.(snm).(['ts_',per]).data,s.(snm).(['ts_',per]).date]=grp_ts(s.(snm).ndbc_sea_t.data,s.(snm).ndbc_sea_t.date,perfun,cumfun,minN);

  [s.(snm).(['sst_',per]).data,s.(snm).(['sst_',per]).date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst.data,s.(snm).hourly_avhrr_weekly_sst.date,perfun,cumfun,minN);
  [s.(snm).(['sst_xs_',per]).data,s.(snm).(['sst_xs_',per]).date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_xshore.data,s.(snm).hourly_avhrr_weekly_sst_xshore.date,perfun,cumfun,minN);
  [s.(snm).(['sst_ls_',per]).data,s.(snm).(['sst_ls_',per]).date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_lshore.data,s.(snm).hourly_avhrr_weekly_sst_lshore.date,perfun,cumfun,minN);
  try,
    [s.(snm).(['sst_l_',per]).data,s.(snm).(['sst_l_',per]).date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_l.data,s.(snm).hourly_avhrr_weekly_sst_l.date,perfun,cumfun,minN);
    [s.(snm).(['sst_dl_',per]).data,s.(snm).(['sst_dl_',per]).date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_dl.data,s.(snm).hourly_avhrr_weekly_sst_dl.date,perfun,cumfun,minN);
  catch,
  end;

  [s.(snm).(['ta_',per]).data,s.(snm).(['ta_',per]).date]=grp_ts(s.(snm).ndbc_air_t.data,s.(snm).ndbc_air_t.date,perfun,cumfun,minN);
  [s.(snm).(['W_kts_',per]).data,s.(snm).(['W_kts_',per]).date]=grp_ts(s.(snm).ndbc_wind1_speed.data,s.(snm).ndbc_wind1_speed.date,perfun,cumfun,minN);
  [s.(snm).(['U_',per]).data,s.(snm).(['U_',per]).date]=grp_ts(s.(snm).ndbc_wind1_xshore.data,s.(snm).ndbc_wind1_xshore.date,perfun,cumfun,minN);
  [s.(snm).(['V_',per]).data,s.(snm).(['V_',per]).date]=grp_ts(s.(snm).ndbc_wind1_lshore.data,s.(snm).ndbc_wind1_lshore.date,perfun,cumfun,minN);
  [s.(snm).(['W_',per]).data,s.(snm).(['W_',per]).date]=grp_ts(s.(snm).ndbc_wind.data,s.(snm).ndbc_wind.date,perfun,cumfun,minN);
  [s.(snm).(['W2_',per]).data,s.(snm).(['W2_',per]).date]=grp_ts(s.(snm).ndbc_wind2.data,s.(snm).ndbc_wind2.date,perfun,cumfun,minN);


  [s.(snm).(['e_ta_',per]).data,s.(snm).(['e_ta_',per]).date]=grp_ts(s.(snm).erai_air_t.data,s.(snm).erai_air_t.date,perfun,cumfun,minN);
  [s.(snm).(['e_qa_',per]).data,s.(snm).(['e_qa_',per]).date]=grp_ts(s.(snm).erai_spechumid.data,s.(snm).erai_spechumid.date,perfun,cumfun,minN);
  [s.(snm).(['e_W_kts_',per]).data,s.(snm).(['e_W_kts_',per]).date]=grp_ts(s.(snm).erai_wind_speed.data,s.(snm).erai_wind_speed.date,perfun,cumfun,minN);
  [s.(snm).(['e_U_',per]).data,s.(snm).(['e_U_',per]).date]=grp_ts(s.(snm).erai_wind_xshore.data,s.(snm).erai_wind_xshore.date,perfun,cumfun,minN);
  [s.(snm).(['e_V_',per]).data,s.(snm).(['e_V_',per]).date]=grp_ts(s.(snm).erai_wind_lshore.data,s.(snm).erai_wind_lshore.date,perfun,cumfun,minN);
  [s.(snm).(['e_W_',per]).data,s.(snm).(['e_W_',per]).date]=grp_ts(s.(snm).erai_wind.data,s.(snm).erai_wind.date,perfun,cumfun,minN);
  [s.(snm).(['e_W2_',per]).data,s.(snm).(['e_W2_',per]).date]=grp_ts(s.(snm).erai_wind2.data,s.(snm).erai_wind2.date,perfun,cumfun,minN);

  % NOTE USE OF "@nansum" FOR INSOLATION!
  [s.(snm).(['e_qswi_',per]).data,s.(snm).(['e_qswi_',per]).date]=grp_ts(s.(snm).erai_dsrf_adj.data,s.(snm).erai_dsrf_adj.date,perfun,@nansum,minN);
  [s.(snm).(['e_qlwi_',per]).data,s.(snm).(['e_qlwi_',per]).date]=grp_ts(s.(snm).erai_dlrf_adj.data,s.(snm).erai_dlrf_adj.date,perfun,cumfun,minN);

  [s.(snm).(['e_hs_',per]).data,s.(snm).(['e_hs_',per]).date]=grp_ts(s.(snm).erai_sigwavehgt_adj.data,s.(snm).erai_sigwavehgt_adj.date,perfun,cumfun,minN);

  % Calculate one-day differences with gaps removed
  s.(snm).(['ts_',per,'_dif']).date=s.(snm).(['ts_',per]).date(2:end); s.(snm).(['ts_',per,'_dif']).data=diff(s.(snm).(['ts_',per]).data);
  s.(snm).(['ta_',per,'_dif']).date=s.(snm).(['ta_',per]).date(2:end); s.(snm).(['ta_',per,'_dif']).data=diff(s.(snm).(['ta_',per]).data);
  s.(snm).(['e_ta_',per,'_dif']).date=s.(snm).(['e_ta_',per]).date(2:end); s.(snm).(['e_ta_',per,'_dif']).data=diff(s.(snm).(['e_ta_',per]).data);
  s.(snm).(['e_qa_',per,'_dif']).date=s.(snm).(['e_qa_',per]).date(2:end); s.(snm).(['e_qa_',per,'_dif']).data=diff(s.(snm).(['e_qa_',per]).data);
  s.(snm).(['sst_',per,'_dif']).date=s.(snm).(['sst_',per]).date(2:end); s.(snm).(['sst_',per,'_dif']).data=diff(s.(snm).(['sst_',per]).data);
  % REMOVE GAPS...
  s.(snm) = filter_gaps(s.(snm),['ts_',per],['ts_',per,'_dif'],maxgaps);
  s.(snm) = filter_gaps(s.(snm),['ta_',per],['ta_',per,'_dif'],maxgaps);
  s.(snm) = filter_gaps(s.(snm),['e_ta_',per],['e_ta_',per,'_dif'],maxgaps);
  s.(snm) = filter_gaps(s.(snm),['e_qa_',per],['e_qa_',per,'_dif'],maxgaps);
  s.(snm) = filter_gaps(s.(snm),['sst_',per],['sst_',per,'_dif'],maxgaps);


  %% REMOVE ANOMALY???
  % s.(snm).(['U_',per]) = anomalize_daily_mean_ts(s.(snm).(['U_',per]));
  % s.(snm).(['V_',per]) = anomalize_daily_mean_ts(s.(snm).(['V_',per]));
  % s.(snm).(['W_',per]) = anomalize_daily_mean_ts(s.(snm).(['W_',per]));
  % s.(snm).(['W2_',per]) = anomalize_daily_mean_ts(s.(snm).(['W2_',per]));
  % s.(snm).(['e_hs_',per]) = anomalize_daily_mean_ts(s.(snm).(['e_hs_',per]));
  % s.(snm).(['sst_xs_',per]) = anomalize_daily_mean_ts(s.(snm).(['sst_xs_',per]));
  % s.(snm).(['sst_ls_',per]) = anomalize_daily_mean_ts(s.(snm).(['sst_ls_',per]));
  % try,
  %   s.(snm).(['sst_l_',per]) = anomalize_daily_mean_ts(s.(snm).(['sst_l_',per]));
  %   s.(snm).(['sst_dl_',per]) = anomalize_daily_mean_ts(s.(snm).(['sst_dl_',per]));
  % catch,
  % end;
  % s.(snm).(['ts_',per,'_dif']) = anomalize_daily_mean_ts(s.(snm).(['ts_',per,'_dif']));
  % s.(snm).(['ta_',per,'_dif']) = anomalize_daily_mean_ts(s.(snm).(['ta_',per,'_dif']));
  % s.(snm).(['e_ta_',per,'_dif']) = anomalize_daily_mean_ts(s.(snm).(['e_ta_',per,'_dif']));
  % s.(snm).(['e_qa_',per,'_dif']) = anomalize_daily_mean_ts(s.(snm).(['e_qa_',per,'_dif']));
  % s.(snm).(['sst_',per,'_dif']) = anomalize_daily_mean_ts(s.(snm).(['sst_',per,'_dif']));

  if ( ismember(snm,{'lonf1','smkf1'}) )
    s.(snm) = station_dewp_to_relhumid(s.(snm),'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
    s.(snm) = station_relhumid_to_spechumid(s.(snm),'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
    [s.(snm).(['qa_',per]).data,s.(snm).(['qa_',per]).date]=grp_ts(s.(snm).ndbc_spechumid.data,s.(snm).ndbc_spechumid.date,perfun,cumfun,minN);
    s.(snm).(['qa_',per,'_dif']).date=s.(snm).(['qa_',per]).date(2:end); s.(snm).(['qa_',per,'_dif']).data=diff(s.(snm).(['qa_',per]).data);
    s.(snm) = filter_gaps(s.(snm),['qa_',per],['qa_',per,'_dif'],maxgaps);
  end;


  if ( doMultiRegress )

    [ts,ta,qa,W2,U,V,sstx,sstl] = ...
        intersect_tses(s.(snm).(['ts_',per,'_dif']),...
                       s.(snm).(['ta_',per,'_dif']),...
                       s.(snm).(['e_qa_',per,'_dif']),...
                       s.(snm).(['W2_',per]),...
                       s.(snm).(['U_',per]),...
                       s.(snm).(['V_',per]),...
                       s.(snm).(['sst_xs_',per]),...
                       s.(snm).(['sst_ls_',per]));

    % R2~0.368
    X = [ts.data,ta.data,qa.data,W2.data,U.data,V.data,sstx.data,sstl.data];
    % % R2~0.348
    % X = [ts.data,ta.data,qa.data,W2.data];
    % % R2~0.330
    % X = [ts.data,ta.data,W2.data];
    % % R2~0.294
    % X = [ts.data,ta.data,qa.data];
    % % R2~0.302
    % X = [ts.data,ta.data,U.data,V.data];
    % % R2~0.360
    % X = [ts.data,ta.data,W2.data,U.data,V.data];
    % % R2~0.280
    % X = [ts.data,ta.data,sstx.data,sstl.data];
    % % R2~0.276
    % X = [ts.data,ta.data];

    % Normalize
    x=(X-repmat(nanmean(X),[size(X,1),1]))./repmat(nanstd(X),[size(X,1),1]);
    % Bartlett test
    [n,p,c] = barttest(x,0.05);
    % Robust multiple regression
    [B,STATS] = robustfit(x(:,2:end),x(:,1));

    % Simple multiple regression
    x(:,end+1) = 1;
    [BR,BRINT,RR,RRINT,STATSR] = regress(x(:,1),x(:,2:end),0.05);
    STATSR(1),

    % % Auto-correlation analysis
    % lags=0:10;
    % fld=['ts_',per,'_dif'];
    % lc=[];
    % for lix=1:numel(lags);
    %   lag=lags(lix);
    %   [ix,lagix]=intersect_dates(s.(snm).(fld).date,s.(snm).(fld).date+lag);
    %   lc(lix)=corr2(s.(snm).(fld).data(ix),s.(snm).(fld).data(lagix));
    %   disp([round(lag),lc(lix)]);
    % end;
    % fmg; plot(0:max(lags),lc(0<=lags));

  end; %if ( doMultiRegress )



  if ( doCovs )

    [c,r,p,ix]=cov_ts(s.(snm).(['ta_',per,'_dif']),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper Tair vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    sacs = sacs + c;
    sars = sars + r;

    [c,r,p,ix]=cov_ts(s.(snm).(['e_ta_',per,'_dif']),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' ERAI Dper Tair vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    esacs = esacs + c;
    esars = esars + r;


    if ( ismember(snm,{'lonf1','smkf1'}) )
      [c,r,p,ix]=cov_ts(s.(snm).(['qa_',per,'_dif']),s.(snm).(['ts_',per,'_dif']),seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' Dper qa vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    end;

    [c,r,p,ix]=cov_ts(s.(snm).(['e_qa_',per,'_dif']),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' ERAI Dper qa vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);


    [c,r,p,ix]=cov_ts(s.(snm).(['U_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Per U vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).(['V_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Per V vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).(['W_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Per Wind vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    swcs = swcs + c;
    swrs = swrs + r;

    [c,r,p,ix]=cov_ts(s.(snm).(['e_W_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' ERAI per Wind vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    eswcs = eswcs + c;
    eswrs = eswrs + r;

    [c,r,p,ix]=cov_ts(s.(snm).(['W2_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Per Wind2 vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    sw2cs = sw2cs + c;
    sw2rs = sw2rs + r;

    [c,r,p,ix]=cov_ts(s.(snm).(['e_W2_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' ERAI per Wind2 vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    esw2cs = esw2cs + c;
    esw2rs = esw2rs + r;


    if (0)
    end;

    [c,r,p,ix]=cov_ts(s.(snm).ndbc_air_t,s.(snm).ndbc_sea_t,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Tair vs. Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    % sacs = sacs + c;
    % sars = sars + r;

    [c,r,p,ix]=cov_ts(s.(snm).ndbc_wind1_speed,s.(snm).ndbc_sea_t,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Wind vs. Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    % swcs = swcs + c;
    % swrs = swrs + r;

    [c,r,p,ix]=cov_ts(s.(snm).ndbc_wind1_speed,s.(snm).ndbc_sea_air_diff,seasfilt);
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
      % [c,r,p,ix]=cov_ts(s.(snm).ndbc_sea_air_diff,s.(snm).ndbc_wind1_speed,seasfilt);
      [c,r,p,ix]=cov_ts(s.(snm).expcdew,s.(snm).ndbc_wind1_speed,seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' exp(C*Tdew) vs. exp(-C*Ta): R=' num2str(r) ' (',num2str(p),')']);
      windrs = windrs + r;

      expta = s.(snm).ndbc_air_t;
      expta.data = exp(0.06*expta.data);
      [c,r,p,ix]=cov_ts(s.(snm).ndbc_relhumid,expta,seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' RH vs. exp(0.06*Ta): R=' num2str(r) ' (',num2str(p),')']);
      rhrs = rhrs + r;
    end;

    [c,r,p,ix]=cov_ts(s.(snm).(['e_qswi_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Adj ERAI QSWI vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).(['e_qlwi_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Adj ERAI QLWI vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).(['e_hs_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Adj ERAI Hs vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).(['sst_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' AVHRR SST vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).(['sst_',per,'_dif']),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).(['sst_xs_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST XS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).(['sst_ls_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST LS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    try,
      [c,r,p,ix]=cov_ts(s.(snm).(['sst_l_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST L vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

      [c,r,p,ix]=cov_ts(s.(snm).(['sst_dl_',per]),s.(snm).(['ts_',per,'_dif']),seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST DL vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    catch,
    end;

  end; %if ( doCovs )

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
