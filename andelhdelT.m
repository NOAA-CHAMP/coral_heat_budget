1;

interpMethod = 'linear';

% stnm = 'fwyf1';
% stnm = 'mlrf1';
% stnm = 'tnrf1';
% stnm = 'lonf1'; interpMethod = 'triangle,linear';
% stnm = 'smkf1';
% stnm = 'looe1';
% stnm = 'sanf1'; interpMethod = 'triangle,linear';
stnm = 'dryf1';

if ( ~exist('stn','var') || ~isfield(stn,'station_name') || ~strcmpi(stn.station_name,stnm) )
  stn = []; clear stn;
  stn = get_station_from_station_name(stnm);
end;

if ( ~isfield(stn,'isobath_orientation') )
  stn = station_optimal_isobath_orientation(stn);
end;

if ( ~isfield(stn,'gom_hycom_seatemp_field') )
  stn = get_gom_hycom(stn);
end;
if ( ~isfield(stn,'gom_hycom_seatemp_x') )
  stn = calc_field_terms(stn,'gom_hycom_seatemp_field','gom_hycom_seatemp',interpMethod);
end;


stn = station_reorient_vectors(stn,'isobath_orientation','gom_hycom_seatemp_x','gom_hycom_seatemp_y');
stn = station_reorient_field(stn,'isobath_orientation','gom_hycom_seatemp_field');

stn = verify_variable(stn,'gom_hycom_seatemp_x_30_day_lowpass');
stn = verify_variable(stn,'gom_hycom_seatemp_y_30_day_lowpass');
stn = verify_variable(stn,'gom_hycom_seatemp_l_30_day_lowpass');
stn = verify_variable(stn,'gom_hycom_seatemp_xshore_30_day_lowpass');
stn = verify_variable(stn,'gom_hycom_seatemp_lshore_30_day_lowpass');

stn = verify_variable(stn,'gom_hycom_seatemp_x_90_day_lowpass');
stn = verify_variable(stn,'gom_hycom_seatemp_y_90_day_lowpass');
stn = verify_variable(stn,'gom_hycom_seatemp_l_90_day_lowpass');
stn = verify_variable(stn,'gom_hycom_seatemp_xshore_90_day_lowpass');
stn = verify_variable(stn,'gom_hycom_seatemp_lshore_90_day_lowpass');

stn = get_avhrr_weekly_field(stn,true);
stn = station_reorient_vectors(stn,'isobath_orientation','avhrr_weekly_sst_x','avhrr_weekly_sst_y');
stn = station_reorient_field(stn,'isobath_orientation','avhrr_weekly_sst_field');

stn = verify_variable(stn,'avhrr_weekly_sst_x_30_day_lowpass');
stn = verify_variable(stn,'avhrr_weekly_sst_y_30_day_lowpass');
stn = verify_variable(stn,'avhrr_weekly_sst_l_30_day_lowpass');
stn = verify_variable(stn,'avhrr_weekly_sst_xshore_30_day_lowpass');
stn = verify_variable(stn,'avhrr_weekly_sst_lshore_30_day_lowpass');

stn = verify_variable(stn,'avhrr_weekly_sst_x_90_day_lowpass');
stn = verify_variable(stn,'avhrr_weekly_sst_y_90_day_lowpass');
stn = verify_variable(stn,'avhrr_weekly_sst_l_90_day_lowpass');
stn = verify_variable(stn,'avhrr_weekly_sst_xshore_90_day_lowpass');
stn = verify_variable(stn,'avhrr_weekly_sst_lshore_90_day_lowpass');


if ( ~exist('res','var') )
  res = [];
end;

[stn,res] = get_tsg_station(stn,res,20);

seasname = {'JFM','AMJ','JAS','OND'};
figspath = get_thesis_path('../figs'); fbasename = [lower(stn.station_name) '-tsg-dsst-'];

for seas=1:4; rad=0.175; ix=find(get_season(stn.smooth_tsg.date)==seas); yrs = get_year(stn.smooth_tsg.date(ix)); fmg; set(gca,'color',[.75,.75,.75]); map_freef([stn.lon-rad,stn.lon+rad,stn.lat-rad,stn.lat+rad],'none'); scatter(stn.smooth_tsg.lon(ix),stn.smooth_tsg.lat(ix),(40.*(abs(stn.smooth_tsg.dsst(ix)).^0.5))+eps,stn.smooth_tsg.dsst(ix)); caxis([-.8,.8]); colormap(flipud(flag(3))); cbh=colorbar; ylabel(cbh,'Cross-shore \partial_xSST [^oC/km]'); set_pcolor_cursor; plot(stn.lon,stn.lat,'kp'); quiver(stn.lon,stn.lat,sind(stn.isobath_orientation+90),cosd(stn.isobath_orientation+90),rad/5,'k:','MarkerSize',10); titlename([upper(stn.station_name) ' Years ' num2str(min(yrs)) '-' num2str(max(yrs)) ' Season ' seasname{seas}]); print('-dtiff',fullfile(figspath,[fbasename,seasname{seas},'.tiff'])); end;


fmg; plot_ts(stn.gom_hycom_seatemp_xshore,stn.avhrr_weekly_sst_xshore,stn.smooth_tsg_ts,'ro');
legend('GOM','AVHRR 1km','SFP TSG'); titlename('\nablah^.\nablaT: Model, Satellite, Truth');

fmg; plot_ts(stn.gom_hycom_seatemp_xshore_30_day_lowpass,stn.avhrr_weekly_sst_xshore_30_day_lowpass,stn.smooth_tsg_ts,'ro');
legend('GOM','AVHRR 1km','SFP TSG (raw)'); titlename('\nablah^.\nablaT 30dLP: Model, Satellite, Truth');

fmg; plot_ts(stn.gom_hycom_seatemp_xshore_90_day_lowpass,stn.avhrr_weekly_sst_xshore_90_day_lowpass,stn.smooth_tsg_ts,'ro');
legend('GOM','AVHRR 1km','SFP TSG (raw)'); titlename('\nablah^.\nablaT 90dLP: Model, Satellite, Truth');


rad=0.075;
minlon=stn.lon-rad; maxlon=stn.lon+rad;
minlat=stn.lat-rad; maxlat=stn.lat+rad;

% dt = datenum(2003,11,29);
dt = datenum(2004,12,15);
% dt = datenum(2005,08,13);

gomix = find(stn.gom_hycom_seatemp_field.date==dt);
sstix = find(stn.avhrr_weekly_sst_field.date==dt);

ix=gomix; fld='gom_hycom_seatemp_field'; ixes=ix-4:ix+3; fmg; contourf(stn.(fld).lon,stn.(fld).lat,squeeze(nanmean(stn.(fld).gradient_xshore(ixes,:,:)))); colorbar; axis([minlon,maxlon,minlat,maxlat,0,1,-6e-4,+6e-4]); plot(stn.lon,stn.lat,'kp'); titlename([upper(stn.station_name) ' \nablah^.\nabla' strrep(fld,'_','\_') ' ' datestr(stn.(fld).date(min(ixes))) ' to ' datestr(stn.(fld).date(max(ixes)))]);

ix=sstix; fld='avhrr_weekly_sst_field'; ixes=ix-4:ix+3; fmg; contourf(stn.(fld).lon,stn.(fld).lat,squeeze(nanmean(stn.(fld).gradient_xshore(ixes,:,:)))); colorbar; axis([minlon,maxlon,minlat,maxlat,0,1,-6e-4,+6e-4]); plot(stn.lon,stn.lat,'kp'); titlename([upper(stn.station_name) ' \nablah^.\nabla' strrep(fld,'_','\_') ' ' datestr(stn.(fld).date(min(ixes))) ' to ' datestr(stn.(fld).date(max(ixes)))]);



mos=10:12;
fld='avhrr_weekly_sst_field'; ixes=find(ismember(get_month(stn.(fld).date),mos)); fmg; contourf(stn.(fld).lon,stn.(fld).lat,squeeze(nanmean(stn.(fld).gradient_x(ixes,:,:)))); colorbar; axis([minlon,maxlon,minlat,maxlat,0,1,-6e-4,+6e-4]); plot(stn.lon,stn.lat,'kp'); titlename([upper(stn.station_name) ' i^\^^.\nabla' strrep(fld,'_','\_') ' months ' num2str(mos)]);
fld='avhrr_weekly_sst_field'; ixes=find(ismember(get_month(stn.(fld).date),mos)); fmg; contourf(stn.(fld).lon,stn.(fld).lat,squeeze(nanmean(stn.(fld).gradient_xshore(ixes,:,:)))); colorbar; axis([minlon,maxlon,minlat,maxlat,0,1,-6e-4,+6e-4]); plot(stn.lon,stn.lat,'kp'); titlename([upper(stn.station_name) ' \nablah^.\nabla' strrep(fld,'_','\_') ' months ' num2str(mos)]);
