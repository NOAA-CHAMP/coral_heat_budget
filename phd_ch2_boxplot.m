1;

tic,

figspath = get_thesis_path('../figs');

for cstnm={'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1'};
%for cstnm={'smkf1'};

    stnm=cstnm{:};
    stn=get_station_from_station_name(stnm);
    stn=load_all_ndbc_data(stn);
    if (isfield(stn,'ndbc_dew_t'));
        stn=station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
        stn=station_relhumid_to_spechumid(stn,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
    end;

    stn=station_optimal_isobath_orientation(stn);
    stn=verify_variable(stn,{'ndbc_wind1_u','ndbc_wind1_v'});
    stn=station_reorient_vectors(stn,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');

    %for cfld={'ndbc_sea_t','ndbc_air_t','ndbc_barom','ndbc_wind1_speed','ndbc_wind1_xshore','ndbc_wind1_lshore','ndbc_tide','ndbc_spechumid'};
    for cfld={'ndbc_sea_t','ndbc_air_t','ndbc_wind1_speed','ndbc_wind1_xshore','ndbc_wind1_lshore','ndbc_spechumid'};
        fld=cfld{:};
        if (isfield(stn,fld));
            disp([stnm,'.',fld]);
            fh = figure;
            boxplot_ts(stn.(fld),[],'allcolors','k',...
                       'title',[upper(stnm),' ',strrep(fld,'_','\_')]);
            ylabel('');
            grid on;
            switch (fld),
             case 'ndbc_sea_t',         ylim([0,40]);
             case 'ndbc_air_t',         ylim([0,40]);
             case 'ndbc_wind1_speed',   ylim([0,40]);
             case 'ndbc_wind1_xshore',  ylim([-40,+40]);
             case 'ndbc_wind1_lshore',  ylim([-40,+40]);
             case 'ndbc_spechumid',     ylim([0,0.03]);
            end;
            set(gca,'FontSize',12);
            print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',fld,'.tif']));
            close(fh);
        end;
    end;
    stn=[]; clear stn;
    pause(1);
end;

clear figspath cstnm stnm cfld fld fh
toc,
