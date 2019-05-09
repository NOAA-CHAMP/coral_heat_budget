1;

if ( ~isfield(stn,'adseatemp_ndbc_erai_erai_30a_net_flux_term') )
  stn = optimize_station_heat_budget('looe1','erai','avhrr_weekly','ndbc','tpxo_tide','erai','ad_seatemp');
end;
stn = verify_variable(stn,'adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_sum');


for dr=['x','l'];
  for fld={'adcp_baroclinic_btm_','adcp_btm_','adcp_baroclinic_sfc_','adcp_sfc_','adcp_'};
    % %plot_spec(looe1,[fld{:},dr],[],[],[],[],'tiff');
    % plot_spec(looe1,[fld{:},dr]);
  end;
end;



fmg;

strd=40;
fld = 'adcp_baroclinic_x_40hlp';
[cs,ch]=contourf(stn.(fld).date(1:strd:end),...
                 stn.adcp_bin_heights,...
                 stn.(fld).prof(1:strd:end,:)');
clabel(cs,ch); colorbar; datetick3;
ylim([0,34]);
caxis([-0.25,+0.25]);

lhs(1)=cs(1); legs{1}='ADCP u^x^s';

lhs(2)=...
    plot(stn.adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_sum.date,...
         stn.adseatemp_ndbc_erai_erai_30a_net_flux_term_40_h_sum.data+25,'b-');
legs{2}='\Sigma_4_0_hQ_0 + 25';

lhs(3)=plot_ts(stn.adcp_seatemp,'r');
legs{3}='T_s';

lhs(4)=...
    plot(stn.adseatemp_ndbc_erai_erai_30a_wind_stress.date,stn.adseatemp_ndbc_erai_erai_30a_wind_stress.data,'k:');
legs{4}='\tau';

lhs(5)=plot_ts(stn.ndbc_wind1_speed_72_h_lp,'k');
legs{5}='U_1_0';

legend(lhs,legs, 'Location','Best');

titlename([upper(stn.station_name),' Q_0 vs. ADCP ',strrep(fld,'_','\_')]);


% xlim(datenum([2008,2008],[06,09],[17,04]));
xlim(datenum([2006,2007],[10,02],[17,09]));
