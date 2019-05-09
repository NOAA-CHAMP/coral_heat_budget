1;

if ( ~exist('stn','var') )
    stn = get_station_from_station_name('looe1');
end;
if( ~isfield(stn,'adcp_u') )
    stn = get_looe1_adcp(stn);
end;
%stn = rmfield(stn,grepstruct(stn,'adcp.*_(speed|dir)_'));
if( ~isfield(stn,'adcp_baroclinic_btm_u') )
  error('Did GET_LOOE1_ADCP fail to calculate adcp_baroclinic_btm_u?');
end;

if( ~isfield(stn,'b_microcat_seatemp_erai_erai_30a_net_heat_flux') )
  stn = station_heat_budget(stn,'erai','avhrr_weekly_sst','erai',[],'erai','microcat_seatemp');
end;

stn = verify_variable(stn,'adcp_baroclinic_btm_u_24_h_sum');
stn = verify_variable(stn,'microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum');
stn = verify_variable(stn,'b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum');


multiplot_station(stn,{...
    'microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum',...
    'adcp_baroclinic_btm_u_24_h_sum'});

multiplot_station(stn,{...
    'b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum',...
    'adcp_baroclinic_btm_u_24_h_sum'});



ix = 1:numel(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum.date);

scatter_fit_ts(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ix(1:24:end),[],'\Sigma_2_4_hQ_0/\rhoC_ph (NO warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')

scatter_fit_ts(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ix(1:24:end),[],'\Sigma_2_4_h(Q_0+Q_b)/\rhoC_ph (warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')

stn = verify_variable(stn,'erai_wind_speed_24_h_lp');

if (0)

ixix = find(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum.data>-0.6);
[ix,ig]=intersect_dates(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum.date(ixix),stn.erai_wind_speed_24_h_lp.date(stn.erai_wind_speed_24_h_lp.data<=9));
scatter_fit_ts(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ixix(ix),[],'\Sigma_2_4_hQ_0/\rhoC_ph (warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')

[ix,ig]=intersect_dates(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum.date(ixix),stn.erai_wind_speed_24_h_lp.date(stn.erai_wind_speed_24_h_lp.data<=9)); ixixix=ixix(ix);
scatter_fit_ts(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ixixix(1:24:end),[],'\Sigma_2_4_hQ_0/\rhoC_ph (warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')
[ix,ig]=intersect_dates(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum.date(ixix),stn.erai_wind_speed_24_h_lp.date(stn.erai_wind_speed_24_h_lp.data<=12)); ixixix=ixix(ix);
scatter_fit_ts(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ixixix(1:24:end),[],'\Sigma_2_4_hQ_0/\rhoC_ph (warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')
[ix,ig]=intersect_dates(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum.date(ixix),stn.erai_wind_speed_24_h_lp.date(stn.erai_wind_speed_24_h_lp.data<=15)); ixixix=ixix(ix);
scatter_fit_ts(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ixixix(1:24:end),[],'\Sigma_2_4_hQ_0/\rhoC_ph (warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')
[ix,ig]=intersect_dates(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum.date(ixix),stn.erai_wind_speed_24_h_lp.date(stn.erai_wind_speed_24_h_lp.data>12)); ixixix=ixix(ix);
scatter_fit_ts(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ixixix(1:24:end),[],'\Sigma_2_4_hQ_0/\rhoC_ph (warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')
scatter_fit_ts(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,@(x)(x(1:24:end)),[],'\Sigma_2_4_hQ_0/\rhoC_ph (warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')
scatter_fit_ts(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,@(x)(1:24:numel(x.date)),[],'\Sigma_2_4_hQ_0/\rhoC_ph (warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')
ixix = find(stn.b_microcat_seatemp_erai_erai_30a_net_heat_flux_term_24_h_sum.data>-0.6);
nansummary(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum.data)
ixix = find(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum.data>-0.5);
[ix,ig]=intersect_dates(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum.date(ixix),stn.erai_wind_speed_24_h_lp.date(stn.erai_wind_speed_24_h_lp.data<15)); ixixix=ixix(ix);
scatter_fit_ts(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ixixix,[],'\Sigma_2_4_hQ_0/\rhoC_ph (no warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')
ixix = find(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum.data>-0.9);
[ix,ig]=intersect_dates(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum.date(ixix),stn.erai_wind_speed_24_h_lp.date(stn.erai_wind_speed_24_h_lp.data<12)); ixixix=ixix(ix);
scatter_fit_ts(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ixixix,[],'\Sigma_2_4_hQ_0/\rhoC_ph (no warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')
scatter_fit_ts(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,stn.adcp_baroclinic_btm_u_24_h_sum,ixixix(1:24:end),[],'\Sigma_2_4_hQ_0/\rhoC_ph (no warm layer)','\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg')
scatter_fit_ts(stn.adcp_baroclinic_btm_u_24_h_sum,stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,[],ixixix(1:24:end),'\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg','\Sigma_2_4_hQ_0/\rhoC_ph (no warm layer)')
scatter_fit_ts(stn.adcp_baroclinic_btm_u_24_h_sum,stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,@(x)(find(-1<=x.data&x.data<=+1)),ixixix(1:24:end),'\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg','\Sigma_2_4_hQ_0/\rhoC_ph (no warm layer)')
ixix = find(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum.data>-20);
[ix,ig]=intersect_dates(stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum.date(ixix),stn.erai_wind_speed_24_h_lp.date(stn.erai_wind_speed_24_h_lp.data<15)); ixixix=ixix(ix);
scatter_fit_ts(stn.adcp_baroclinic_btm_u_24_h_sum,stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,@(x)(find(-1<=x.data&x.data<=+1)),ixixix(1:24:end),'\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg \epsilon [-1,+1]','\Sigma_2_4_hQ_0/\rhoC_ph (no warm layer, wind_2_4_h_l_p<15kts)')
scatter_fit_ts(stn.adcp_baroclinic_btm_u_24_h_sum,stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,@(x)(find(-5<=x.data&x.data<=+5)),ixixix(1:24:end),'\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg \epsilon [-5,+5]','\Sigma_2_4_hQ_0/\rhoC_ph (no warm layer, wind_2_4_h_l_p<15kts)')
scatter_fit_ts(stn.adcp_baroclinic_btm_u_24_h_sum,stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,[],ixixix(1:24:end),'\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg \epsilon [-5,+5]','\Sigma_2_4_hQ_0/\rhoC_ph (no warm layer, wind_2_4_h_l_p<15kts)'), %@(x)(find(-5<=x.data&x.data<=+5))
scatter_fit_ts(stn.adcp_baroclinic_btm_u_24_h_sum,stn.microcat_seatemp_erai_erai_30_net_heat_flux_term_24_h_sum,[],ixixix(1:24:end),'\Sigma_2_4_hu^.(\nablah=0) bottom - depth-avg','\Sigma_2_4_hQ_0/\rhoC_ph (no warm layer, wind_2_4_h_l_p<15kts)'), %@(x)(find(-5<=x.data&x.data<=+5))

end; %if(0)
