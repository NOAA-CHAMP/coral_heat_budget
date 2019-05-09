1;

stn=get_station_from_station_name('looe1'); stn=get_looe1_adcp(stn);
find_date_ranges(stn.adcp_x.date,(3/24))
find_date_ranges(stn.adcp_sfc_x.date,(3/24))
datenum(2008,6,11)-datenum(2007,2,25),
ix=18777:30095;
stn.x.date=stn.adcp_x.date(ix); stn.x.data=stn.adcp_x.data(ix);
stn.l.date=stn.adcp_l.date(ix); stn.l.data=stn.adcp_l.data(ix);
stn.s.date=stn.adcp_sfc_speed.date(ix); stn.s.data=stn.adcp_sfc_speed.data(ix);
plot_spec(stn,'x',[],[],[],[1e-5,1e3]); titlename('LOOE1 Cross-shore'); print('-dtiff','../figs/looe1-xs-spec.tiff');
plot_spec(stn,'l',[],[],[],[1e-5,1e3]); titlename('LOOE1 Alongshore'); print('-dtiff','../figs/looe1-ls-spec.tiff');
plot_spec(stn,'s',[],[],[],[1e-5,1e3]); titlename('LOOE1 Near-surface'); print('-dtiff','../figs/looe1-sfc-spec.tiff');
