1;

stn.ndbc_sea_t_degapped = degap_days_ts(stn.ndbc_sea_t);

stn.ndbc_sea_t_diff.date=stn.ndbc_sea_t_degapped.date(1:end-1);
stn.ndbc_sea_t_diff.data=diff(stn.ndbc_sea_t_degapped.data);
stn=filter_gaps(stn,'ndbc_sea_t_degapped','ndbc_sea_t_diff',(6/24));

stn=verify_variable(stn,'ndbc_sea_t_diff_24_hour_sum');

stn = verify_variable(stn,'ndbc_sea_t_degapped_1_day_average');

stn.ndbc_sea_t_degapped_1_day_average_diff.date=stn.ndbc_sea_t_degapped.date(1:end-1);
stn.ndbc_sea_t_degapped_1_day_average_diff.data=diff(stn.ndbc_sea_t_degapped.data);
stn=filter_gaps(stn,'ndbc_sea_t_degapped','ndbc_sea_t_degapped_1_day_average_diff',(6/24));


stn.t.date=unique(floor(stn.ndbc_sea_t_degapped.date));
stn.t.data=grpstats(stn.ndbc_sea_t_degapped.data,floor(stn.ndbc_sea_t_degapped.date));
stn.q0.date=unique(floor(stn.benthic_ndbc_erai_30a_net_heat_flux_term.date));
stn.q0.data=grpstats(stn.benthic_ndbc_erai_30a_net_heat_flux_term.data,floor(stn.benthic_ndbc_erai_30a_net_heat_flux_term.date),'sum');
stn.hc.date=unique(floor(stn.ndbc_erai_30a_ww3_gom_qe_hc_dTdt.date));
stn.hc.data=grpstats(stn.ndbc_erai_30a_ww3_gom_qe_hc_dTdt.data,floor(stn.ndbc_erai_30a_ww3_gom_qe_hc_dTdt.date),'sum');
stn.sgs.date=unique(floor(stn.ndbc_erai_30a_ww3_gom_qe_hc_sgs_final_budget_24_hour_average.date));
stn.sgs.data=grpstats(stn.ndbc_erai_30a_ww3_gom_qe_hc_sgs_final_budget_24_hour_average.data,floor(stn.ndbc_erai_30a_ww3_gom_qe_hc_sgs_final_budget_24_hour_average.date),'sum');

stn.td.date=stn.t.date(1:end-1);
stn.td.data=diff(stn.t.data);

scatter_fit_ts(stn.td,stn.q0)
scatter_fit_ts(stn.td,stn.hc)
scatter_fit_ts(stn.td,stn.sgs)

scatter_fit_ts(stn.td,stn.q0,[],@(x)(find(x.data~=0)))
scatter_fit_ts(stn.td,stn.hc,[],@(x)(find(x.data~=0)))
scatter_fit_ts(stn.td,stn.sgs,[],@(x)(find(x.data~=0)))

scatter_fit_ts(stn.td,stn.hc,[],@(x)(find(x.data~=0)),[],[],[],[],true,1)
scatter_fit_ts(stn.td,stn.hc,[],@(x)(find(x.data~=0)),[],[],[],[],true,2)
scatter_fit_ts(stn.td,stn.hc,[],@(x)(find(x.data~=0)),[],[],[],[],true,3)
scatter_fit_ts(stn.td,stn.hc,[],@(x)(find(x.data~=0)),[],[],[],[],true,4)
scatter_fit_ts(stn.td,stn.hc,[],@(x)(find(x.data~=0)),[],[],[],[],true,-1)
scatter_fit_ts(stn.td,stn.hc,[],@(x)(find(x.data~=0)),[],[],[],[],true,-2)
scatter_fit_ts(stn.td,stn.hc,[],@(x)(find(x.data~=0)),[],[],[],[],true,-3)

scatter_fit_ts(stn.td,stn.sgs,[],@(x)(find(x.data~=0)),[],[],[],[],true,1)
scatter_fit_ts(stn.td,stn.hc,[],@(x)(find(x.data~=0)),[],[],[],[],true,1)


stn.wt.date=unique(get_yearweek(stn.ndbc_sea_t_degapped.date));
stn.wt.data=grpstats(stn.ndbc_sea_t_degapped.data,get_yearweek(stn.ndbc_sea_t_degapped.date));

stn.wq0.data=grpstats(stn.benthic_ndbc_erai_30a_net_heat_flux_term.data,get_yearweek(stn.benthic_ndbc_erai_30a_net_heat_flux_term.date),'sum');
stn.wq0.date=unique(get_yearweek(stn.benthic_ndbc_erai_30a_net_heat_flux_term.date));

stn.whc.date=unique(get_yearweek(stn.ndbc_erai_30a_ww3_gom_qe_hc_dTdt.date));
stn.whc.data=grpstats(stn.ndbc_erai_30a_ww3_gom_qe_hc_dTdt.data,get_yearweek(stn.ndbc_erai_30a_ww3_gom_qe_hc_dTdt.date),'sum');

stn.wsgs.date=unique(get_yearweek(stn.ndbc_erai_30a_ww3_gom_qe_hc_sgs_final_budget_24_hour_average.date));
stn.wsgs.data=grpstats(stn.ndbc_erai_30a_ww3_gom_qe_hc_sgs_final_budget_24_hour_average.data,get_yearweek(stn.ndbc_erai_30a_ww3_gom_qe_hc_sgs_final_budget_24_hour_average.date),'sum');

stn.wtd.date=stn.wt.date(1:end-1);
stn.wtd.data=diff(stn.wt.data);

scatter_fit_ts(stn.wtd,stn.wq0,[],@(x)(find(x.data~=0)),[],[],[],[],true,0)
scatter_fit_ts(stn.wtd,stn.wq0,[],@(x)(find(x.data~=0)),[],[],[],[],true,7)

scatter_fit_ts(stn.wtd,stn.whc,[],@(x)(find(x.data~=0)),[],[],[],[],true,0)
scatter_fit_ts(stn.wtd,stn.whc,[],@(x)(find(x.data~=0)),[],[],[],[],true,7)

scatter_fit_ts(stn.wtd,stn.wsgs,[],@(x)(find(x.data~=0)),[],[],[],[],true,0)
scatter_fit_ts(stn.wtd,stn.wsgs,[],@(x)(find(x.data~=0)),[],[],[],[],true,7)
