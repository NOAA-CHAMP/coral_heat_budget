    looe = []; clear looe;
    %looe = load('Looe_Ley_ADCP_25FEB2007_thru_19OCT2007.mat');
    %looe = load('Looe_Key_ADCP_19OCT2007_thru_11JUN2008.mat');
    looe = load('Looe_Key_ADCP_21JUN2008_thru_25JAN2009.mat');
    dts = datenum(2000+looe.SerYear,looe.SerMon,looe.SerDay,looe.SerHour,looe.SerMin,looe.SerSec);
    datestr(dts(1:8)),
    looe.adcp_u.date = dts;
    jumpix = find(diff(looe.adcp_u.date) <= 0);
    datestr(looe.adcp_u.date(jumpix(1)-2:jumpix(1)+5)),

    %jumpix = find(diff(looe.adcp_u.date) <= 0 | diff(looe.adcp_u.date) > (1.1/24));

  bn1=34;
  bn2=2;
  figure;
  plot(looe1.adcp_u.date,([looe1.adcp_dir.prof(:,bn1)-looe1.adcp_dir.prof(:,bn2)]));
  maxigraph; datetick3;
  titlename(['Dir (' num2str(bn1) '-' num2str(bn2) ')']);
  hold on;
  plot(stn.ndbc_ncep_30a_net_heat_flux.date,real(stn.ndbc_ncep_30a_net_heat_flux.data),'r-');
