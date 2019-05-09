1;
% Script to reproduce figures for FKEYS HYCOM vs. AVHRR SST/Ardhuin empirical currents estimates during a cold-front in 2007

stnm = 'mlrf1';

if ( ~exist('stn','var') )
    stn = get_station_from_station_name(stnm);
end;
if ( ~isfield(stn,'isobath_orientation') )
    stn = station_optimal_isobath_orientation(stn);
end;
if ( ~isfield(stn,'fkeys_hycom_u_field') )
    stn = get_fkeys_hycom(stn);
end;

if ( ~isfield(stn,'fkeys_hycom_x_field') )
    x=[]; clear x;
    x.x.lat=stn.fkeys_hycom_u_field.lat;
    x.x.lon=stn.fkeys_hycom_u_field.lon;
    x.x.u_field=stn.fkeys_hycom_u_field.field;
    x.x.v_field=stn.fkeys_hycom_v_field.field;

    x = station_reorient_field(x,stn.isobath_orientation,'x','u_field','v_field','x_field','l_field');

    stn.fkeys_hycom_x_field.date  = stn.fkeys_hycom_u_field.date;
    stn.fkeys_hycom_x_field.field = x.x.x_field;
    stn.fkeys_hycom_x_field.lat   = stn.fkeys_hycom_u_field.lat;
    stn.fkeys_hycom_x_field.lon   = stn.fkeys_hycom_u_field.lon;

    stn.fkeys_hycom_l_field.date  = stn.fkeys_hycom_u_field.date;
    stn.fkeys_hycom_l_field.field = x.x.l_field;
    stn.fkeys_hycom_l_field.lat   = stn.fkeys_hycom_u_field.lat;
    stn.fkeys_hycom_l_field.lon   = stn.fkeys_hycom_u_field.lon;
    x=[]; clear x;
end;

if ( ~isfield(stn.fkeys_hycom_seatemp_field,'gradient_xs') )
    x=[]; clear x;
    x.x.lat=stn.fkeys_hycom_seatemp_field.lat;
    x.x.lon=stn.fkeys_hycom_seatemp_field.lon;
    x.x.date=stn.fkeys_hycom_seatemp_field.date;
    x.x.gradient_x=stn.fkeys_hycom_seatemp_field.gradient_x;
    x.x.gradient_y=stn.fkeys_hycom_seatemp_field.gradient_y;

    x = station_reorient_field(x,stn.isobath_orientation,'x','gradient_x','gradient_y','gradient_xs','gradient_ls');

    stn.fkeys_hycom_seatemp_field.gradient_xs = x.x.gradient_xs;
    stn.fkeys_hycom_seatemp_field.gradient_ls = x.x.gradient_ls;
    x=[]; clear x;
end;

dt = datenum(2007,12,17);

[ig,ix] = min(abs(stn.fkeys_hycom_u_field.date-dt));
%fmg; contourf(stn.fkeys_hycom_u_field.lon,stn.fkeys_hycom_u_field.lat,squeeze(stn.fkeys_hycom_u_field.field(ix,:,:))); set(gca,'clim',[-1,+1]); colorbar; plot(stn.lon,stn.lat,'wp',stn.lon,stn.lat,'kp','MarkerSize',8); titlename([upper(stn.station_name),' FKEYS HYCOM u  ',datestr(dt)]);
%fmg; contourf(stn.fkeys_hycom_v_field.lon,stn.fkeys_hycom_v_field.lat,squeeze(stn.fkeys_hycom_v_field.field(ix,:,:))); set(gca,'clim',[-1,+1]); colorbar; plot(stn.lon,stn.lat,'wp',stn.lon,stn.lat,'kp','MarkerSize',8); titlename([upper(stn.station_name),' FKEYS HYCOM v  ',datestr(dt)]);
fmg; contourf(stn.fkeys_hycom_x_field.lon,stn.fkeys_hycom_x_field.lat,squeeze(stn.fkeys_hycom_x_field.field(ix,:,:))); set(gca,'clim',[-1,+1]); colorbar; plot(stn.lon,stn.lat,'wp',stn.lon,stn.lat,'kp','MarkerSize',8); titlename([upper(stn.station_name),' FKEYS HYCOM u\bullet\nablah ',datestr(dt)]);
%fmg; contourf(stn.fkeys_hycom_l_field.lon,stn.fkeys_hycom_l_field.lat,squeeze(stn.fkeys_hycom_l_field.field(ix,:,:))); set(gca,'clim',[-1,+1]); colorbar; plot(stn.lon,stn.lat,'wp',stn.lon,stn.lat,'kp','MarkerSize',8); titlename([upper(stn.station_name),' FKEYS HYCOM KS ',datestr(dt)]);

plot_ngdc_bathy_station(stn,[0:-2:-30,-40:-10:-300],[],[],@contour); colorbar off; axis([-80.43,-80.33,24.96,25.06]); daspect([cosd(25),1,1]);
quiver(stn.fkeys_hycom_u_field.lon,stn.fkeys_hycom_u_field.lat,squeeze(stn.fkeys_hycom_u_field.field(ix,:,:)),squeeze(stn.fkeys_hycom_v_field.field(ix,:,:))); plot(stn.lon,stn.lat,'wp',stn.lon,stn.lat,'kp','MarkerSize',8); titlename([upper(stn.station_name),' FKEYS HYCOM sfc. currents ',datestr(dt)]);
print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-fkeys-currents-2007Dec17.tif']));

[ig,ix] = min(abs(stn.fkeys_hycom_seatemp_field.date-dt));
%fmg; contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,1e3.*squeeze(stn.fkeys_hycom_seatemp_field.gradient_xs(ix,:,:))); set(gca,'clim',[-1,+1]); colorbar; plot(stn.lon,stn.lat,'wp',stn.lon,stn.lat,'kp','MarkerSize',8); titlename([upper(stn.station_name),' FKEYS HYCOM \nablaT_s\bullet\nablah ',datestr(dt)]);
plot_ngdc_bathy_station(stn,[0:-2:-30,-40:-10:-300],[],[],@contour); colorbar off; axis([-80.43,-80.33,24.96,25.06]); daspect([cosd(25),1,1]);
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,1e3.*squeeze(stn.fkeys_hycom_seatemp_field.gradient_xs(ix,:,:)),[-2:0.2:2]);
set(gca,'clim',[-.7,+.7]);
ch=colorbar;
ylabel(ch,'^oC\bulletkm^-^1');
plot(stn.lon,stn.lat,'wp',stn.lon,stn.lat,'kp','MarkerSize',12);
titlename([upper(stn.station_name),' FKEYS HYCOM \nablaT_s\bullet\nablah ',datestr(dt)]);
print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-fkeys-seatemp-gradient_xs-2007Dec17.tif']));
