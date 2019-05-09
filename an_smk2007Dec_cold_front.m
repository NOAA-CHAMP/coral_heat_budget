1;

dtlims = datenum(2007,[11,12],[11,31]);

if ( ~exist('smkf1','var') || ~isfield(smkf1,'ndbc_erai_erai_30a_net_flux') )
    smkf1 = optimize_station_heat_budget('smkf1','erai','avhrr_weekly','ndbc','tpxo_tide','erai');
end;
if ( ~isfield(smkf1,'ndbc_erai_erai_30a_net_flux') )
    error('Failed to get heat budget results??');
end;

smkf1 = verify_variable(smkf1,{'ndbc_erai_erai_30a_net_flux_term_1_d_sum','ndbc_erai_erai_30a_net_flux_term_5_d_sum',...
    'ndbc_erai_erai_30a_avhrr_hc_dTdt_1_d_sum','ndbc_erai_erai_30a_avhrr_hc_dTdt_5_d_sum'});
if ( ~isfield(smkf1,'surface_budget_1d_sum') )
    [smkf1.surface_budget_1d_sum.data,smkf1.surface_budget_1d_sum.date] = grp_ts(smkf1.ndbc_erai_erai_30a_net_flux_term.data,smkf1.ndbc_erai_erai_30a_net_flux_term.date,@floor,@nansum);
    [smkf1.total_budget_1d_sum.data,smkf1.total_budget_1d_sum.date] = grp_ts(smkf1.ndbc_erai_erai_30a_avhrr_hc_dTdt.data,smkf1.ndbc_erai_erai_30a_avhrr_hc_dTdt.date,@floor,@nansum);
end;

%    {'ndbc_air_t','ndbc_sea_t','ndbc_erai_erai_30a_wind_stress_xshore','ndbc_erai_erai_30a_net_flux_term_5_d_sum'},
%    [],[],{'T_a','T_s','\tau^x^,^y','\Sigma_1_dQ/\rhoC_ph, \Sigma_1_d\partial_tT_s'},datenum(2007,[11,12],[21,31]),
[smkf1,lhs,ahs,fh]=multiplot_station(smkf1,...
    {'ndbc_air_t','ndbc_sea_t','ndbc_erai_erai_30a_wind_stress_xshore','surface_budget_1d_sum'},...
    [],[],{'Air temp. ^oC','Sea temp. ^oC','Wind stress N/m^2','Ocean warming ^oC / day'},dtlims,...
    {[12,32],[12,32],[-0.5,+0.5],[-3.5,+3.5]},true,{'b-','b-','b-','b-'});
axes(ahs(3)); hold on; plot(smkf1.ndbc_erai_erai_30a_wind_stress_lshore.date,smkf1.ndbc_erai_erai_30a_wind_stress_lshore.data,'k-');
legend('Location','NorthEast', 'Cross-shore','Alongshore');
%axes(ahs(4)); hold on; plot(smkf1.ndbc_erai_erai_30a_avhrr_hc_dTdt.date,smkf1.ndbc_erai_erai_30a_avhrr_hc_dTdt.data,'k-');
%axes(ahs(4)); hold on; plot(smkf1.ndbc_erai_erai_30a_avhrr_hc_dTdt_5_d_sum.date,smkf1.ndbc_erai_erai_30a_avhrr_hc_dTdt_5_d_sum.data,'k-');
axes(ahs(4)); hold on; plot(smkf1.total_budget_1d_sum.date,smkf1.total_budget_1d_sum.data,'k-');
legend('Location','NorthEast', 'Surface budget','Total budget');

if (0)
    looe1 = get_station_from_station_name('looe1'); looe1 = get_looe1_microcat(looe1); looe1 = get_looe1_adcp(looe1);
    fmg;
end;
