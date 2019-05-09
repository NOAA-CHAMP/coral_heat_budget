1;

if ( ~exist('stn','var') || isempty(stn) )
    stn = optimize_station_heat_budget('mlrf1','erai','avhrr_weekly','ndbc',[],'ww3');
end;

fmg; plot_ts(stn.ww3_avhrr_advected_heat); ylim([-.2,.2]); titlename('Advected Heat');
xlim(datenum([2001,2002],10,1)); datetick3('x',2,'keeplimits');
fmg; plot_ts(stn.ww3_ndbc_stokes_speed); titlename('Stokes current');
xlim(datenum([2001,2002],10,1)); datetick3('x',2,'keeplimits');
fmg; plot_ts(stn.hourly_avhrr_weekly_sst_xshore,stn.hourly_avhrr_weekly_sst_lshore); legend('Cross-shore','Longshore'); titlename('\nablaT_A_V_H_H_R');
xlim(datenum([2001,2002],10,1)); datetick3('x',2,'keeplimits');

fmg; plot_ts(stn.avhrr_weekly_diffused_heat); ylim([-.2,.2]); titlename('Diffused Heat');
xlim(datenum([2001,2002],10,1)); datetick3('x',2,'keeplimits');
fmg; plot_ts(stn.hourly_avhrr_weekly_sst_l,stn.hourly_avhrr_weekly_sst_xx,stn.hourly_avhrr_weekly_sst_yy,ts_op(stn.hourly_avhrr_weekly_sst_xx,stn.hourly_avhrr_weekly_sst_yy,'+')); legend('DEL2','XX','YY','XX+YY'); titlename('\nabla^2T_A_V_H_H_R');
xlim(datenum([2001,2002],10,1)); datetick3('x',2,'keeplimits');

station_plot_fields(stn,'avhrr_weekly_sst_field',datenum(2002,6,10):datenum(2002,6,25));
titlename(['AVHRR SST around ',upper(stn.station_name)]);
