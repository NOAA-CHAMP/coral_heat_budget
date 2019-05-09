1;

if ( ~exist('stn','var') )
    stnm = 'smkf1';
    stn = get_station_from_station_name(stnm);
    stn = load_all_ndbc_data(stn);
    stn = verify_variable(stn,'ndbc_sea_t_7_d_avg');

    dt=+3.5;
    stn = get_avhrr_weekly_field(stn);
    stn.raw_avhrr_weekly_sst.date=stn.raw_avhrr_weekly_sst.date + dt;
    stn.raw_avhrr_weekly_sst_field.date=stn.raw_avhrr_weekly_sst_field.date + dt;
end;

for seas=1:4;
    disp(seas);
    if ( 1 )
        seasix=find(get_season(stn.raw_avhrr_weekly_sst_field.date)==seas);
        for rix=1:17;
            for cix=1:17;
                stn.raw_avhrr_weekly_sst_field.corr2(seas,rix,cix)=...
                    corr2(stn.raw_avhrr_weekly_sst_field.field(seasix,rix,cix),...
                    stn.raw_avhrr_weekly_sst.data(seasix));

                [ix,fx] = intersect_dates(stn.ndbc_sea_t_7_d_avg.date,...
                    stn.raw_avhrr_weekly_sst_field.date(seasix));
                stn.raw_avhrr_weekly_sst_field.ndbc_corr2(seas,rix,cix)=...
                    corr2(stn.raw_avhrr_weekly_sst_field.field(seasix(fx),rix,cix),...
                    stn.ndbc_sea_t_7_d_avg.data(ix));
            end; %for cix
        end; %for rix
    end; %if

    cs = [0.50:0.05:1.00];

    if (0)
        fmg;
        contourf(stn.raw_avhrr_weekly_sst_field.lon,...
            stn.raw_avhrr_weekly_sst_field.lat,...
            squeeze(stn.raw_avhrr_weekly_sst_field.corr2(seas,:,:)),cs);
        caxis([min(cs),max(cs)]);
        colorbar;
        titlename([stn.station_name,' Weekly SST Corr. Coef. by Season ',num2str(seas)]);
        print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(mfilename),'-corr2-',lower(stn.station_name),'.tif']));
    end;

    if (1)
        fmg;
        contourf(stn.raw_avhrr_weekly_sst_field.lon,...
            stn.raw_avhrr_weekly_sst_field.lat,...
            squeeze(stn.raw_avhrr_weekly_sst_field.ndbc_corr2(seas,:,:)),cs);
        caxis([min(cs),max(cs)]);
        colorbar;
        hold on; plot(stn.lon,stn.lat,'wp');
        titlename([stn.station_name,' \mu_7_din situ Corr. Coef. by Season ',num2str(seas)]);
        print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(mfilename),'-ndbc_corr2-',lower(stn.station_name),'.tif']));
    end;

end; %for seas
