1;

for c={'fwyf1','mlrf1','lonf1','smkf1','looe1','sanf1','plsf1'};
    stnm=c{:}; disp(stnm); stn=get_station_from_station_name(stnm);
    if ( strcmpi(stnm,'sanf1') || strcmpi(stnm,'plsf1') || strcmpi(stnm,'dryf1') )
        stn=get_avhrr_weekly_field(stn,true,@nanmean);
    else
        stn=get_avhrr_weekly_field(stn,true);
    end;
    stn=station_optimal_isobath_orientation(stn);
    stn=station_reorient_vectors(stn,'isobath_orientation','avhrr_weekly_sst_x','avhrr_weekly_sst_y');
    station_boxplots(stn,'avhrr_weekly_sst_xshore','\partial_x_sSST',[-1e-3,1e-3],[],[],0,[1,1,0,0],[1,1,0,0]);
    station_boxplots(stn,'avhrr_weekly_sst_lshore','\partial_l_sSST',[-1e-3,1e-3],[],[],0,[1,1,0,0],[1,1,0,0]);
    stn=[]; clear stn; pack
    pause(1);
end;
clear c stnm
