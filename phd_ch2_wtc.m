1;

more off

%dts=datenum([2009,2011],[11,12],[1,31]); 
%dts=datenum([2004,2008],[1,12],[1,31]); 
%dts=datenum([2004,2011],[1,12],[1,31]); 
dts=datenum([2001,2011],[1,12],[1,31]); 
%dts=datenum([1995,2011],[1,12],[1,31]); 

nDays = ceil(dts(2)-dts(1));

%maxScale = ceil(nDays./2);
maxScale = 350;

hlp=40;
rawsfld = 'ndbc_sea_t';
sfld = [rawsfld,'_',num2str(hlp),'_h_decimate'];

rawaflds = {'ndbc_air_t','ndbc_spechumid','ndbc_wind1_speed','ndbc_wind1_lshore','ndbc_wind1_xshore','hourly_avhrr_weekly_sst_lshore','hourly_avhrr_weekly_sst_xshore'};
rawaflds = {'ndbc_air_t'};
rawaflds = {'ndbc_wind1_speed','ndbc_wind1_lshore','ndbc_wind1_xshore','hourly_avhrr_weekly_sst_lshore','hourly_avhrr_weekly_sst_xshore'};

stnms = {'fwyf1','mlrf1','lonf1','smkf1','sanf1'};
stnms = {'mlrf1','lonf1'};
stnms = {'fwyf1','mlrf1'};

for cstnm=stnms
  stnm=cstnm{:};
  disp(upper(stnm));

  if ( ~exist('stn','var') || ~strcmpi(stn.station_name,stnm) || ~isfield(stn,rawsfld) )
    stn=[]; clear stn;
    stn = get_station_from_station_name(stnm);
    stn = load_all_ndbc_data(stn);
  end;

  stn = verify_variable(stn,{sfld,'ndbc_wind1_u','ndbc_wind1_v'});

  if ( isfield(stn,'ndbc_dew_t') && ~isfield(stn,'ndbc_spechumid') )
    stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
    stn = station_relhumid_to_spechumid(stn,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
  end;

  if ( ~isfield(stn,'hourly_avhrr_weekly_sst') )
    switch (stnm),
     case 'sanf1',
      stn = get_avhrr_weekly_field(stn,true,'nearest',3);
     case 'dryf1',
      stn = get_avhrr_weekly_field(stn,true,{@nanmean,2,2},3);
     otherwise,
      stn = get_avhrr_weekly_field(stn,true);
    end;
  end;

  if ( ~isfield(stn,'isobath_orientation') )
    stn = station_optimal_isobath_orientation(stn);
  end;
  if ( ~isfield(stn,'ndbc_wind1_xshore') )
    stn = station_reorient_vectors(stn,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');
  end;
  if ( ~isfield(stn,'hourly_avhrr_weekly_sst_xshore') )
    stn = station_reorient_vectors(stn,'isobath_orientation','hourly_avhrr_weekly_sst_x','hourly_avhrr_weekly_sst_y');
  end;


  for cfld=rawaflds

    rawafld = cfld{:};
    if ( ~isfield(stn,rawafld) )
      warning('No field %s.%s',stnm,rawafld);
      continue;
    end;

    afld = [rawafld,'_',num2str(hlp),'_h_decimate'];
    disp(afld);
    stn = verify_variable(stn,afld);

    aix=find(dts(1)<=stn.(afld).date&stn.(afld).date<dts(2));

    fh=fmg;
    try,
        wtc_ts(stn.(afld),stn.(sfld),aix,[],'maxPoints',(24*nDays),...
            'S0',(hlp/24),'MaxScale',maxScale,'MonteCarloCount',30,'BlackandWhite');
        datetick3('x',17,'keeplimits');
        titlename([upper(stn.station_name),' WTC of ',num2str(hlp),'hlp ',...
            strrep(rawafld,'_','\_'),' vs ',strrep(rawsfld,'_','\_'),': ',datestr(dts(1)),' to ',datestr(dts(2))]);
        %set(gca,'FontSize',7);
        set(gca,'FontWeight','bold');
        print('-dtiff',fullfile(get_thesis_path('../figs'),...
            [lower(stn.station_name),'-',afld,'-wtc-',sfld,'-',...
            datestr(floor(dts(1)),29),'-',datestr(floor(dts(2)),29),'.tif']));
    catch,
        catchwarn;
        close(fh); fh=[];
    end;

    %disp('Hit enter to analyze next field');
    %pause;
    pause(0.5);
    close(fh);

    clear rawafld afld aix fh ans

  end; %for cfld

  pause(1);

end; %for cstnm

clear dts nDays hlp rawsfld sfld rawaflds cfld ans

more on
