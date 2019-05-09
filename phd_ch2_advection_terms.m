1;

stn = get_station_from_station_name('looe1'); stn = get_looe1_adcp(stn);

fmg; hist(stn.adcp_x.data(ts_jfm(stn.adcp_x)),1000); xlim([-1,1]); titlename('LOOE1 cross-shore JFM');
fmg; hist(stn.adcp_x.data(ts_jas(stn.adcp_x)),1000); xlim([-1,1]); titlename('LOOE1 cross-shore JAS');

fmg; hist(stn.adcp_l.data(ts_jfm(stn.adcp_l)),1000); xlim([-1,1]); titlename('LOOE1 longshore JFM');
fmg; hist(stn.adcp_l.data(ts_jas(stn.adcp_l)),1000); xlim([-1,1]); titlename('LOOE1 longshore JAS');

stn=[]; ans=[]; clear stn ans


stn = get_station_from_station_name('mlrf1'); stn = get_avhrr_weekly_field(stn,true);

fmg; hist(stn.raw_avhrr_weekly_sst_x.data(ts_jfm(stn.raw_avhrr_weekly_sst_x)),1000); xlim([-8,+8].*1e-4); titlename('MLRF1 AVHRR \partial_x_sSST JFM');
fmg; hist(stn.raw_avhrr_weekly_sst_x.data(ts_jas(stn.raw_avhrr_weekly_sst_x)),1000); xlim([-8,+8].*1e-4); titlename('MLRF1 AVHRR \partial_x_sSST JAS');

fmg; hist(stn.raw_avhrr_weekly_sst_l.data(ts_jfm(stn.raw_avhrr_weekly_sst_l)),1000); xlim([-8,+8].*1e-4); titlename('MLRF1 AVHRR \partial_l_sSST JFM');
fmg; hist(stn.raw_avhrr_weekly_sst_l.data(ts_jas(stn.raw_avhrr_weekly_sst_l)),1000); xlim([-8,+8].*1e-4); titlename('MLRF1 AVHRR \partial_l_sSST JAS');
