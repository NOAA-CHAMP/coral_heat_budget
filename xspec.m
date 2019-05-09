function stn = xspec(stnm_or_stn,doSpecs)
%function stn = xspec(stnm_or_stn,doSpecs)

  figspath = get_thesis_path('../figs');

  if ( ~exist('doSpecs','var') || isempty(doSpecs) )
    doSpecs = true;
  end;

  stn = get_station_from_station_name(stnm_or_stn);

  if ( ~isfield(stn,'ndbc_air_t') )
    stn = load_all_ndbc_data(stn);
  end;

  if (doSpecs)
    plot_spec(stn,'ndbc_sea_t',[],[],[],[],'tiff');
    plot_spec(stn,'ndbc_air_t',[],[],[],[1e-5,1e7],'tiff');
    plot_spec(stn,'ndbc_wind1_speed',[],[],[],[1e-3,1e7],'tiff');
    plot_spec(stn,'ndbc_barom',[],[],[],[1e-8,1e10],'tiff');

    if ( isfield(stn,'ndbc_dew_t') )
      stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
      stn = station_relhumid_to_spechumid(stn,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
      plot_spec(stn,'ndbc_spechumid',[],[],[],[1e-10,1e1],'tiff');
    end;
    if ( isfield(stn,'ndbc_tide') )
      plot_spec(stn,'ndbc_tide',[],[],[],[1e-6,1e5],'tiff');
    end;
  end;

  if ( strcmpi(stn.station_name,'mlrf1') )
    if ( ~isfield(stn,'bic_surf_par') )
      rawstn = load_station_data(stn);
      flds = grepstruct(rawstn,'bic_');
      for fldix=1:numel(flds)
        stn.(flds{fldix}) = rawstn.(flds{fldix});
        stn.(flds{fldix}).date(stn.(flds{fldix}).data<0) = [];
        stn.(flds{fldix}).data(stn.(flds{fldix}).data<0) = [];
      end;
      rawstn=[]; clear rawstn;
    end;

    if (doSpecs)
      plot_spec(stn,'bic_surf_par',[],[],[],[1e1,1e10],'tiff');
    end;

    x = importdata('mlrf2-cleaning-dates.csv');
    stn.cleaning_date = datenum(x);
    x=[]; clear x


    % UV band of interest (305, 330, or 380nm narrowband)

    uvrng = '380nm';
    sfcuv = 'bic_surf_380nm'; shluv = 'bic_shallow_380nm';
    uvfld = 'kd_380nm'; uvlbl = 'K_d^3^8^0^n^m';
    % uvrng = '330nm';
    % sfcuv = 'bic_surf_330nm'; shluv = 'bic_shallow_330nm';
    % uvfld = 'kd_330nm'; uvlbl = 'K_d^3^3^0^n^m';
    % uvrng = '305nm';
    % sfcuv = 'bic_surf_305nm'; shluv = 'bic_shallow_305nm';
    % uvfld = 'kd_305nm'; uvlbl = 'K_d^3^0^5^n^m';


    fmg;
    plot_ts(stn.bic_surf_par,stn.bic_shallow_par);
    titlename('Molasses Reef MLRF2 PAR');
    xlim([min(stn.bic_surf_par.date),max(stn.bic_surf_par.date)]);
    datetick3('x',2,'keeplimits');
    arrow([stn.cleaning_date(1),2500],[stn.cleaning_date(1),2300]); 
    legend('Surface PAR','U/W PAR','Cleaning', 'Location','SouthWest');
    for ix=2:length(stn.cleaning_date);
      arrow([stn.cleaning_date(ix),2500],[stn.cleaning_date(ix),2300]);
    end;
    print('-dtiff',fullfile(figspath,'mlrf2-par-and-cleanings.tiff'));
    print('-dpng',fullfile(figspath,'mlrf2-par-and-cleanings.png'));

    fmg;
    plot_ts(stn.(sfcuv),stn.(shluv));
    titlename(['Molasses Reef MLRF2 ' uvrng ' UV']);
    xlim([min(stn.(sfcuv).date),max(stn.(sfcuv).date)]);
    datetick3('x',2,'keeplimits');
    arrtop = ceil(max(stn.(sfcuv).data))*1.20;
    arrbtm = ceil(max(stn.(sfcuv).data))*1.00;
    arrow([stn.cleaning_date(1),arrtop],[stn.cleaning_date(1),arrbtm]); 
    legend(['Surface ' uvrng ' UV'],['U/W ' uvrng ' UV'],'Cleaning', 'Location','SouthWest');
    for ix=2:length(stn.cleaning_date);
      arrow([stn.cleaning_date(ix),arrtop],[stn.cleaning_date(ix),arrbtm]);
    end;
    print('-dtiff',fullfile(figspath,['mlrf2-' uvrng '-and-cleanings.tiff']));
    print('-dpng',fullfile(figspath,['mlrf2-' uvrng '-and-cleanings.png']));


    % PAR (400-700nm) diffuse attenuation
    stn = calc_kd(stn,'bic_surf_par','bic_shallow_par',2);
    stn.kd = stn.kd_bic_surf_par_bic_shallow_par;

    % UV diffuse attenuation
    stn = calc_kd(stn,sfcuv,shluv,2);
    stn.(uvfld) = stn.(['kd_' sfcuv '_' shluv]);


    % Only use data taken < 3 hours from local mid-day (16.5-17.5 GMT)
    % Noontime solar elevation year-round lies between 48 and 88 deg Alt
    ix = find( get_daylight(stn.kd.date,stn.lat,stn.lon,35) & ...
               ismember(get_hour(stn.kd.date),15:19) & ...
               (stn.kd.data>0.00001) );
    stn.kd.date = stn.kd.date(ix);
    stn.kd.data = stn.kd.data(ix);
    stn.(uvfld).date = stn.(uvfld).date(ix);
    stn.(uvfld).data = stn.(uvfld).data(ix);
    keepix = [];
    for ix=1:length(stn.cleaning_date)
      dtdif = stn.kd.date - stn.cleaning_date(ix);
      % Skip the day of a cleaning, and stop a week afterward
      keepix = union(keepix,find( 1 <= dtdif & dtdif < 8 ));
    end;
    stn.kd.date = stn.kd.date(keepix);
    stn.kd.data = stn.kd.data(keepix);
    stn.(uvfld).date = stn.(uvfld).date(keepix);
    stn.(uvfld).data = stn.(uvfld).data(keepix);

    fmg;
    [cum,tid] = grp_ts(stn.kd.data,stn.kd.date,@get_week,@nanmedian,2);
    plot(7.*(tid-1),cum,'*', 'MarkerSize',10);
    [cum,tid] = grp_ts(stn.kd.data,stn.kd.date,@get_month,@nanmedian,2);
    plot(datenum(0,tid,15),cum,'rs', 'MarkerSize',12,'MarkerFaceColor','r');
    xlim([0,366]); datetick3('x',3,'keeplimits');
    ylabel('K_d [m^-^1]');
    ylim([0,0.4]);
    titlename('Median K_d^P^A^R');
    legend('Weekly','Monthly');
    print('-dtiff',fullfile(figspath,'mlrf2-clean-kd-par-climatology.tiff'));
    print('-dpng',fullfile(figspath,'mlrf2-clean-kd-par-climatology.png'));

    fmg;
    [cum,tid] = grp_ts(stn.(uvfld).data,stn.(uvfld).date,@get_week,@nanmedian,2);
    plot(7.*(tid-1),cum,'*', 'MarkerSize',10);
    [cum,tid] = grp_ts(stn.(uvfld).data,stn.(uvfld).date,@get_month,@nanmedian,2);
    plot(datenum(0,tid,15),cum,'rs', 'MarkerSize',12,'MarkerFaceColor','r');
    xlim([0,366]); datetick3('x',3,'keeplimits');
    ylabel('K_d [m^-^1]');
    ylim([0,0.4]);
    titlename(['Median ' uvlbl]);
    legend('Weekly','Monthly');
    print('-dtiff',fullfile(figspath,['mlrf2-clean-' uvfld '-climatology.tiff']));
    print('-dpng',fullfile(figspath,['mlrf2-clean-' uvfld '-climatology.png']));

    station_boxplots(stn,'kd','K_d^P^A^R [m^-^1]',[0,1],[],[],0,[0,1,1,0],[0,1,1,0]);
    station_boxplots(stn,uvfld,[uvlbl ' [m^-^1]'],[0,1],[],[],0,[0,1,1,0],[0,1,1,0]);

    % % [cum,tid] = grp_ts(stn.kd.data,stn.kd.date,[],'numel',2);
    % % fmg; plot(tid,cum,'*'); xlim([0,366]); datetick3('x',3,'keeplimits');
    % [cum,tid] = grp_ts(stn.kd.data,stn.kd.date,@get_week,'numel',2);
    % fmg; plot(7.*(tid-1),cum,'*'); xlim([0,366]); datetick3('x',3,'keeplimits');
    % titlename('Total # Data Points');
    % % print('-dtiff',fullfile(figspath,'mlrf2-clean-kd-par-climatology-N.tiff'));
    % % print('-dpng',fullfile(figspath,'mlrf2-clean-kd-par-climatology-N.png'));

    yrs=unique(get_year(stn.kd.date));
    cs = {'k*','b*','r*','m*'};
    lhs = [];
    fmg;
    for yrix=1:numel(yrs)
      ix = find(get_year(stn.kd.date)==yrs(yrix));
      lh=plot(get_yearday(stn.kd.date(ix)),stn.kd.data(ix),cs{yrix});
      lhs(end+1) = lh(1);
      xlim([0,366]);
      datetick3('x',3,'keeplimits');
    end;
    legend(lhs,num2str(yrs(:)));
    ylabel('K_d [m^-^1]');
    ylim([0,1]);
    titlename('K_d^P^A^R vs. year-day');
    print('-dtiff',fullfile(figspath,'mlrf2-clean-kd-par-by-month-year.tiff'));
    print('-dpng',fullfile(figspath,'mlrf2-clean-kd-par-by-month-year.png'));

    lhs = [];
    fmg;
    for yrix=1:numel(yrs)
      ix = find(get_year(stn.(uvfld).date)==yrs(yrix));
      lh=plot(get_yearday(stn.(uvfld).date(ix)),stn.(uvfld).data(ix),cs{yrix});
      lhs(end+1) = lh(1);
      xlim([0,366]);
      datetick3('x',3,'keeplimits');
    end;
    legend(lhs,num2str(yrs(:)));
    ylabel('K_d [m^-^1]');
    ylim([0,1]);
    titlename([uvlbl ' vs. year-day']);
    print('-dtiff',fullfile(figspath,['mlrf2-clean-' uvfld '-par-by-month-year.tiff']));
    print('-dpng',fullfile(figspath,['mlrf2-clean-' uvfld '-par-by-month-year.png']));

  end;

return;
