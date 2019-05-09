1;

if ( ~exist('stn','var') )
    stn = get_station_from_station_name('looe1');
    stn = station_optimal_isobath_orientation(stn);
end;
if ( ~isfield(stn,'adcp_x') )
    x=get_looe1_adcp; stn.adcp_u=x.adcp_u; stn.adcp_v=x.adcp_v; stn.adcp_x=x.adcp_x; stn.adcp_l=x.adcp_l; x=[]; clear x; pack;
    %stn = get_looe1_microcat(stn);
end;

if ( ~isfield(stn,'gom_hycom_xshore') )
    stn = get_gom_hycom(stn);
    stn = station_reorient_vectors(stn,'isobath_orientation','gom_hycom_u','gom_hycom_v');
end;

if ( ~isfield(stn,'fkeys_hycom_xshore') )
    stn = get_fkeys_hycom(stn);
    stn = station_reorient_vectors(stn,'isobath_orientation','fkeys_hycom_u','fkeys_hycom_v');
end;


stn = verify_variable(stn,{'adcp_x_1_d_avg','gom_hycom_xshore_1_d_avg','fkeys_hycom_xshore_1_d_avg','adcp_l_1_d_avg','gom_hycom_lshore_1_d_avg','fkeys_hycom_lshore_1_d_avg'});

fmg;
spt(2,3,1); hist(stn.adcp_x_1_d_avg.data,100);                xlim([-0.25,+0.25]);
spt(2,3,2); hist(stn.gom_hycom_xshore_1_d_avg.data,100);      xlim([-0.25,+0.25]);
spt(2,3,3); hist(stn.fkeys_hycom_xshore_1_d_avg.data,100);    xlim([-0.25,+0.25]);
spt(2,3,4); hist(stn.adcp_l_1_d_avg.data,100);                xlim([-1.00,+1.00]);
spt(2,3,5); hist(stn.gom_hycom_lshore_1_d_avg.data,100);      xlim([-1.00,+1.00]);
spt(2,3,6); hist(stn.fkeys_hycom_lshore_1_d_avg.data,100);    xlim([-1.00,+1.00]);
print('-dtiff',fullfile(get_thesis_path('../figs'),'ch2_fig2.tiff'));
