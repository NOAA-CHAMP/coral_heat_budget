function stn = get_avhrr_weekly_reoriented(stn_or_stnm,varargin)
%function stn = get_avhrr_weekly_reoriented(stn_or_stnm,varargin)

  stn = get_station_from_station_name(stn_or_stnm);
  stn = station_optimal_isobath_orientation(stn);
  stn = load_all_ndbc_data(stn);
  stn = get_avhrr_weekly_field(stn,true,varargin{:});
  stn = station_reorient_vectors(stn,'isobath_orientation','hourly_avhrr_weekly_sst_x','hourly_avhrr_weekly_sst_y');
  stn = station_reorient_field(stn,'isobath_orientation','avhrr_weekly_sst_field');

  stn = verify_variable(stn,'ndbc_wind1_u_24_h_lp');
  stn = verify_variable(stn,'ndbc_wind1_v_24_h_lp');

  stn.Ulp.date = stn.avhrr_weekly_sst_field.date;
  stn.Ulp.data = interp1(stn.ndbc_wind1_u_24_h_lp.date,stn.ndbc_wind1_u_24_h_lp.data,stn.Ulp.date);
  stn.Vlp.date = stn.avhrr_weekly_sst_field.date;
  stn.Vlp.data = interp1(stn.ndbc_wind1_v_24_h_lp.date,stn.ndbc_wind1_v_24_h_lp.data,stn.Vlp.date);

  if (0)
    fmg;
    subplot_tight(2,1,1);
    boxplot_ts(stn.hourly_avhrr_weekly_sst_xshore);
    xlabel('Cross-shore');
    ylim([-1,1]*1e-3);
    titlename([stn.station_name,'\partial_x_s_,_l_sT_A_V_H_R_R']);
    subplot_tight(2,1,2);
    boxplot_ts(stn.hourly_avhrr_weekly_sst_lshore);
    xlabel('Along-shore');
    ylim([-1,1]*1e-3);
    grid on;
  end;


  maxVlp = 5;

  for minUlp=5:-5:-20;
    maxUlp = minUlp + 5;
    % % for seas=1:4
    % for seas=4:4
    for mon=11:11
      % for cfld={'field','gradient_xshore','gradient_lshore','laplacian'}
      for cfld={'gradient_xshore'}
        fld = cfld{:};
        fmg;
        % ix = find(get_season(stn.avhrr_weekly_sst_field.date)==seas ...
        ix = find(get_year(stn.avhrr_weekly_sst_field.date)>1996 ...
                  & get_month(stn.avhrr_weekly_sst_field.date)==mon ...
                  & minUlp<=stn.Ulp.data & stn.Ulp.data<maxUlp & abs(stn.Vlp.data)<maxVlp);
        contourf(stn.avhrr_weekly_sst_field.lon,stn.avhrr_weekly_sst_field.lat,...
                 squeeze(nanmean(stn.avhrr_weekly_sst_field.(fld)(ix,:,:),1)));
        if ( strcmp(fld,'field') );	caxis([24,27]); %caxis([16,34]);
        else;				caxis([-3,3]*1e-4); end;
        colorbar;
        plot(stn.lon,stn.lat,'wp');
        titlename([stn.station_name,' month ',num2str(mon),' ',strrep(fld,'_','\_'),...
                   ' U_2_4_h_l_p\epsilon [',num2str([minUlp,maxUlp]),') |V_2_4_h_l_p|<',...
                   num2str(maxVlp),' N=',num2str(numel(ix))]);
      end;
    end;
  end;

return;
