1;

if (~exist('stn','var') || ~strcmpi(stn.station_name,'looe1') )
  stn = []; clear stn;
  stn = get_station_from_station_name('looe1'); stn=get_looe1_microcat(stn); stn=get_looe1_adcp(stn);
end;

plot_spec(stn,'adcp_x'); grid off; xlim([0.1,10]);

for yr=2005:2009; stn.adcp_x_warm=subset_ts(stn.adcp_x,@(x)(datenum(yr,7,1)<=x.date&x.date<datenum(yr,10,1))); stn.adcp_x_cool=subset_ts(stn.adcp_x,@(x)(datenum(yr,12,1)<=x.date&x.date<datenum(yr+1,3,1))); plot_spec(stn,{'adcp_x_warm','adcp_x_cool'}); grid off; appendtitlename([' ',num2str(yr)]); xlim([0.1,10]); end;

for yr=2005:2009; stn.adcp_btm_x_warm=subset_ts(stn.adcp_btm_x,@(x)(datenum(yr,7,1)<=x.date&x.date<datenum(yr,10,1))); stn.adcp_btm_x_cool=subset_ts(stn.adcp_btm_x,@(x)(datenum(yr,12,1)<=x.date&x.date<datenum(yr+1,3,1))); plot_spec(stn,{'adcp_btm_x_warm','adcp_btm_x_cool'}); grid off; appendtitlename([' ',num2str(yr)]); xlim([0.1,10]); end;

for yr=2005:2009; stn.adcp_sfc_x_warm=subset_ts(stn.adcp_sfc_x,@(x)(datenum(yr,7,1)<=x.date&x.date<datenum(yr,10,1))); stn.adcp_sfc_x_cool=subset_ts(stn.adcp_sfc_x,@(x)(datenum(yr,12,1)<=x.date&x.date<datenum(yr+1,3,1))); plot_spec(stn,{'adcp_sfc_x_warm','adcp_sfc_x_cool'}); grid off; appendtitlename([' ',num2str(yr)]); xlim([0.1,10]); end;

% stn=[]; clear stn
% pack
