1;

fmg; boxplot_tses({stn.ndbc_sea_t_dly_diff_flux,stn.b_ndbc_erai_erai_30a_net_flux_dly,stn.b_ndbc_erai_erai_30a_avhrr_dt_flux_dly,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly},[],'legend',{'Actual','Q_0+Q_b','Q_0+Q_b+\rhoC_ph(F_qu\bullet\nablaT_s+K\nabla^2T-s)','\partial_tT_s'},'widthandcolor',true);
