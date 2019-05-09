1;

doPrint = false;
doRegression = false;

figspath = get_thesis_path('../figs');

fld = 'amodis_kd_par';
% fld = 'amodis_kd_par_anom';
% fld = 'amodis_kd_par_simple_anom';
fldstr = strrep(fld,'_','\_');

for cst={'fwyf1','mlrf1','lonf1','smkf1','sanf1'};
  st=cst{:};
  stns.(st) = read_usf_kd_par(st);

  if ( doRegression )
    stns.(st) = load_all_ndbc_data(stns.(st));
    stns.(st) = verify_variable(stns.(st),{'ndbc_wind1_speed_30_d_lp','ndbc_wind1_u_30_d_lp','ndbc_wind1_v_30_d_lp'});
    stns.(st) = station_optimal_isobath_orientation(stns.(st));

    % stns.(st) = station_reorient_vectors(stns.(st),'isobath_orientation',...
    %                                      'ndbc_wind1_u_30_d_lp','ndbc_wind1_v_30_d_lp',...
    %                                      'ndbc_wind1_xshore_30_d_lp','ndbc_wind1_lshore_30_d_lp');
    stns.(st) = station_reorient_vectors(stns.(st),'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');
    stns.(st) = verify_variable(stns.(st),{'ndbc_wind1_xshore_30_d_lp','ndbc_wind1_lshore_30_d_lp'});
  end;
end;


if ( 1 )
  cs = {'k','b','r','m','y'};
  sts={'fwyf1','mlrf1','lonf1','smkf1','sanf1'};

  fmg;
  for stix=1:numel(sts)
    st=sts{stix};
    c = cs{stix};
    plot_ts(stns.(st).amodis_kd_par_max,[c,'--']);
    lhs=plot_ts(stns.(st).(fld),[c,'-']);
    plot_ts(stns.(st).amodis_kd_par_min,[c,':']);
    lh(stix) = lhs(1);
    clear ig lhs
  end;
  legend(lh, upper(sts), 'Location','Best');
  if ( ~isempty(regexp(fld,'anom')) )
    ylim([-0.4,+0.4]);
  else
    ylim([0.0,1.0]);
  end;
  titlename([fldstr,' time series']);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,['USF_',fld,'_time_series.tiff']));
  end;

  fmg;
  for stix=1:numel(sts)
    st=sts{stix};
    c = cs{stix};
    grpplot_ts(stns.(st).amodis_kd_par_max,@get_month,@nanmedian,[],[c,'--']);
    [ig,lhs]=grpplot_ts(stns.(st).(fld),@get_month,@nanmedian,[],[c,'-']);
    grpplot_ts(stns.(st).amodis_kd_par_min,@get_month,@nanmedian,[],[c,':']);
    lh(stix) = lhs(1);
    clear ig lhs
  end;
  legend(lh, upper(sts), 'Location','Best');
  if ( ~isempty(regexp(fld,'anom')) )
    ylim([-0.4,+0.4]);
  else
    ylim([0.0,1.0]);
  end;
  titlename([fldstr,' climatologies']);
  if ( doPrint )
    print('-dtiff',fullfile(figspath,['USF_',fld,'_climatologies.tiff']));
  end;
end;


if ( doRegression )

  subfn = [];
  % subfn = @ts_jfm;
  % subfn = @ts_amj;
  % subfn = @ts_jas;
  % subfn = @ts_ond;

  for cst={'fwyf1','mlrf1','lonf1','smkf1','sanf1'};
    st=cst{:};
    if ( ~isempty(regexp(fld,'anom')) )
      axlm = [-20,+10,-0.3,+0.3];
    else
      axlm = [-20,+10,0.1,0.4];
      if ( strcmpi(st,'lonf1') ); axlm = [-20,+10,0.5,0.8]; end;
    end;

    for csubfn = {@ts_jfm,@ts_amj,@ts_jas,@ts_ond};
      subfn = csubfn{:};

      ix = [];
      if ( isa(subfn,'function_handle') )
        ix = subfn(stns.(st).(fld));
      end;

      fh = fmg;
      spt(2,2,1);
      scatter_fit_ts(stns.(st).ndbc_wind1_u_30_d_lp,stns.(st).(fld),...
                     [],ix,'U^x_3_0_d_l_p',fldstr,fh);
      axis(axlm); legend('Location','SouthEast');
      spt(2,2,2);
      scatter_fit_ts(stns.(st).ndbc_wind1_v_30_d_lp,stns.(st).(fld),...
                     [],ix,'U^y_3_0_d_l_p',fldstr,fh);
      axis(axlm); legend('Location','SouthEast');
      spt(2,2,3);
      scatter_fit_ts(stns.(st).ndbc_wind1_xshore_30_d_lp,stns.(st).(fld),...
                     [],ix,'U^x^s_3_0_d_l_p',fldstr,fh);
      axis(axlm); legend('Location','SouthEast');
      spt(2,2,4);
      scatter_fit_ts(stns.(st).ndbc_wind1_lshore_30_d_lp,stns.(st).(fld),...
                     [],ix,'U^l^s_3_0_d_l_p',fldstr,fh);
      axis(axlm); legend('Location','SouthEast');
      suptitlename([upper(st),' K_d vs. wind components ',upper(char(subfn))]);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[lower(st),'_',fld,'_scatter_wind_',lower(char(subfn)),'.tiff']));
      end;
    end;


    fmg;
    hist(stns.(st).(fld).data(ix),50);
    if ( ~isempty(regexp(fld,'anom')) )
      axis([-0.3,+0.3,0,15]);
    else
      axis([0.1,0.8,0,11]);
    end;
    titlename([upper(st),' ',fldstr,' histogram']);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[lower(st),'_',fld,'_hist.tiff']));
    end;
  end;

end;

clear doPrint doRegression figspath fld fldstr csubfn subfn axlm fh cst st cs sts c
