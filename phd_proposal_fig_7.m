1;

%cstnms={'lkwf1'};
%cstnms={'lkwf1','mlrf1','lonf1','smkf1','42003'};
cstnms={'lkwf1','fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1','42003'};

for cstnm=cstnms
  stnm = cstnm{:};
  disp(stnm);

  stn=[]; clear stn;

if ( ~exist('stn','var') )
    % if ( ~exist('stnm','var') )
    %     %stnm='fwyf1';
    %     %stnm='mlrf1';
    %     stnm='lonf1';
    %     %stnm='dryf1';
    % end;
    stn = get_station_from_station_name(stnm);
    stn = load_all_ndbc_data(stn);
    % stn = qa_ts(stn,'ndbc_sea_t');
    % stn = qa_ts(stn,'ndbc_air_t');
    % stn = qa_ts(stn,'ndbc_wind1_speed');

    [vl,dr] = intersect_tses(stn.ndbc_wind1_speed,stn.ndbc_wind1_dir);

    stn.ndbc_wind1_speed_mps.date = vl.date;
    stn.ndbc_wind1_speed_mps.data = kts2mps(vl.data);
    stn.ndbc_wind1_u_mps.date = stn.ndbc_wind1_speed_mps.date;
    stn.ndbc_wind1_v_mps.date = stn.ndbc_wind1_speed_mps.date;
    [stn.ndbc_wind1_u_mps.data,stn.ndbc_wind1_v_mps.data] = ...
        spddir_to_uv(stn.ndbc_wind1_speed_mps.data,dr.data);
end;

%stn = verify_variable(stn,{'ndbc_wind1_u_mps_7_d_dev','ndbc_wind1_v_mps_7_d_dev'});
%    'ndbc_wind1_u_mps_3_d_avg_7_d_dev_0_d_asof_sum_ndbc_wind1_v_mps_3_d_avg_7_d_dev',...

flds = {...
    'ndbc_sea_t',...
    'ndbc_sea_t_1_d_dev_3_d_avg',...
    'ndbc_air_t',...
    'ndbc_air_t_1_d_dev_3_d_avg',...
    'ndbc_wind1_u_mps',...
    'ndbc_wind1_v_mps',...
    'ndbc_wind1_speed_mps_3_d_avg',...
    'ndbc_wind1_u_mps_7_d_dev_0_d_asof_sum_ndbc_wind1_v_mps_7_d_dev',...
    };
ylbls = {...
    'T_s [^oC]',...
    '\mu_3_d(\sigma_1_dT_s) [^oC]',...
    'T_a [^oC]',...
    '\mu_3_d(\sigma_1_dT_a) [^oC]',...
    'W_U [m/s]',...
    'W_V [m/s]',...
    '\mu_3_d(W) [m/s]',...
    '\sigma_7_d(W_U)+\sigma_7_d(W_V)',...
    };

stn = verify_variable(stn,flds);

%    datenum([1987,2012],[6,6],[1,1]),...
%    {[0,40],[0,2],[0,40],[0,4],[-50,+50],[-50,+50],[0,20],[0,20]});
multiplot_station(stn,flds,[],'',ylbls,...
    datenum([1987,2011],[6,12],[1,31]),...
    {[5,35],[0,2],[5,35],[0,4],[-30,+30],[-30,+30],[0,30],[0,30]});
datetick3('x',10,'keeplimits');

if (1)
  fignm = fullfile(get_thesis_path('../figs'),[lower(mfilename),'-',stnm,'.tif']);
  disp(['Printing ',fignm]);
  print('-dtiff',fignm);
end;

end;
