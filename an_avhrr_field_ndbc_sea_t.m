1;

if ( ~exist('stn','var') )
    stnm = 'smkf1';
    stn = get_station_from_station_name(stnm);
    stn = load_all_ndbc_data(stn);
    stn = get_avhrr_field(stn);
end;

for seas=1:4;
    disp(seas);
    if ( 1 )
        seasix=find(get_season(stn.avhrr_sst.date)==seas);
        for rix=1:size(stn.avhrr_sst_field.field,2);
            for cix=1:size(stn.avhrr_sst_field.field,3);
                stn.avhrr_sst_field.corr2(seas,rix,cix)=...
                    corr2(stn.avhrr_sst_field.field(seasix,rix,cix),...
                    stn.avhrr_sst.data(seasix));

                [ix,fx] = intersect_dates(stn.ndbc_sea_t.date,...
                    stn.avhrr_sst_field.date(seasix));
                stn.avhrr_sst_field.ndbc_corr2(seas,rix,cix)=...
                    corr2(stn.avhrr_sst_field.field(seasix(fx),rix,cix),...
                    stn.ndbc_sea_t.data(ix));
            end; %for cix
        end; %for rix
    end; %if

    cs = [0.50:0.05:1.00];

    if (1)
        fmg;
        contourf(stn.avhrr_sst_field.lon,...
            stn.avhrr_sst_field.lat,...
            squeeze(stn.avhrr_sst_field.corr2(seas,:,:)),cs);
        caxis([min(cs),max(cs)]);
        colorbar;
        hold on; plot(stn.lon,stn.lat,'wp');
        titlename([stn.station_name,' Synoptic SST Corr. Coef. by Season ',num2str(seas)]);
        print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(mfilename),'-synoptic-corr2-',lower(stn.station_name),'.tif']));
    end;

    if (1)
        fmg;
        contourf(stn.avhrr_sst_field.lon,...
            stn.avhrr_sst_field.lat,...
            squeeze(stn.avhrr_sst_field.ndbc_corr2(seas,:,:)),cs);
        caxis([min(cs),max(cs)]);
        colorbar;
        hold on; plot(stn.lon,stn.lat,'wp');
        titlename([stn.station_name,' \mu_7_din situ Corr. Coef. by Season ',num2str(seas)]);
        print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(mfilename),'-synoptic-ndbc_corr2-',lower(stn.station_name),'.tif']));
    end;

end; %for seas
