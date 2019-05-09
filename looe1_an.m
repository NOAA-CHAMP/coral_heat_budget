1;

%-- Unknown date --%
fmg; plot_ts(lonf1.ndbc_sea_t); titlename('LONF1 Sea T');
lonf1 = load_all_ndbc_data([],'lonf1');
fmg; plot_ts(lonf1.ndbc_sea_t); titlename('LONF1 Sea T');
datetick3('x',2,'keeplimits')
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stns = get_fknms_thermistors();
delete data/FKNMS_thermograph.mat
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stns = get_fknms_thermistors();
for ix=1:55; figure(ix); end;
stnms=grepstruct(stns,'^FKNMS_'); for ix=1:numel(stnms); stnm=stnms{ix}; if (isfield(stns.(stnm),'fknms_seatemp')); figure(100+ix); maxigraph; grid on; plot_ts(stns.(stnm).fknms_seatemp); datetick3('x',17,'keeplimits'); titlename(strrep(stnm,'_','\_')); else disp(['MISSING ' stnm]); end; end;
mlrf1 = get_station_from_station_name('mlrf1');
help sw_dist
sw_dist([mlrf1.lat,stns.FKNMS_HEN_CHIX.lat],[mlrf1.lon,stns.FKNMS_HEN_CHIX.lon],'km')
fmg; plot_ts(stns.FKNMS_HEN_CHIX.fknms_seatemp,stns.FKNMS_HEN_CHIX.hourly_fknms_seatemp_qc);
fmg; plot_ts(stns.FKNMS_HEN_CHIX.fknms_seatemp,stns.FKNMS_HEN_CHIX.hourly_fknms_seatemp_qc,'r');
fmg; plot_ts(stns.FKNMS_HEN_CHIX.fknms_seatemp,'b',stns.FKNMS_HEN_CHIX.hourly_fknms_seatemp_qc,'r');
fmg; plot_ts(stns.FKNMS_HEN_CHIX.fknms_seatemp,stns.FKNMS_HEN_CHIX.hourly_fknms_seatemp_qc);
reviewanim(100:139)
reviewanim(100:139,3)
help close
for ix=100:140; if ishandle(ix); close(ix); end; end;
fmg; plot_ts(stns.FKNMS_HEN_CHIX.fknms_seatemp,'b.',stns.FKNMS_HEN_CHIX.hourly_fknms_seatemp_qc,'r.');
help qa_ts
stn = qa_ts(stns,'FKNMS_HEN_CHIX.fknms_seatemp',7,true);
stn = qa_ts(stns.FKNMS_HEN_CHIX,'fknms_seatemp',7,true);
fmg;
txtcnt = ...
uicontrol('Style', 'text', 'Units', 'normalized', ...
'Position', [0.02 0.53 0.12 0.40], 'FontSize', 7, ...
'HorizontalAlignment', 'left', 'String', 'TESTOLA');
get(txtcnt)
fmg;
txtcnt = uicontrol('Style', 'text', 'Units', 'normalized','Position', [0.02 0.53 0.12 0.40], 'FontSize', 7,'HorizontalAlignment', 'left', 'String', 'Now is the time for all good men to come to the aid of their countries tis a far far better thing than I have ever done before');
txtcnt = uicontrol('Style', 'text', 'Units', 'normalized','Position', [0.02 0.53 0.12 0.40], 'FontSize', 7,'HorizontalAlignment', 'left', 'String', '\tNow is the time for all good men to come to the aid of their countries tis a far far better thing than I have ever done before');
fmg;
txtcnt = uicontrol('Style', 'text', 'Units', 'normalized','Position', [0.02 0.53 0.12 0.40], 'FontSize', 7,'HorizontalAlignment', 'left', 'String', sprintf('\tNow is the time for all good men to come to the aid of their countries tis a far far better thing than I have ever done before'));
fmg; txtcnt = uicontrol('Style', 'text', 'Units', 'normalized','Position', [0.02 0.53 0.12 0.40], 'FontSize', 7,'HorizontalAlignment', 'left', 'String', sprintf('    Now is the time for all good men to come to the aid of their countries tis a far far better thing than I have ever done before'));
stn = qa_ts(stns.FKNMS_HEN_CHIX,'fknms_seatemp',7,true);
edited
stn = qa_ts(stns.FKNMS_HEN_CHIX,'fknms_seatemp',7,true);
stn.hourly_fknms_seatemp_qc
stn = qa_ts(stns.FKNMS_KW_CHANL,'fknms_seatemp',7,true);
stn
help interp_ts
stn = qa_ts(stns.FKNMS_KW_CHANL,'fknms_seatemp',7,true);
dbstop verify_variable
stn = qa_ts(stns.FKNMS_KW_CHANL,'fknms_seatemp',7,true);
op
varname
dbquit
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
%-- 8/25/11  9:25 AM --%
fc
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
pd
help optim_q0
stn = optim_q0('mlrf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('mlrf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
stn = get_fkeys_hycom(stn);
quiver_ts(stn)
quiver_field(stn)
hold on;
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(nanmean(stn.fkeys_hycom_seatemp_field.field,1)));
nansummary(stn.fkeys_hycom_seatemp.data)
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(nanmean(stn.fkeys_hycom_seatemp_field.field,1)),[17:2:31]]);
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(nanmean(stn.fkeys_hycom_seatemp_field.field,1)),[17:2:31]);
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(nanmean(stn.fkeys_hycom_seatemp_field.field,1)),[17:2:31]);f
contour(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(nanmean(stn.fkeys_hycom_seatemp_field.field,1)),[17:2:31]);
quiver_field(stn,[],[],[],'fkeys_hycom_seatemp_field)
stn.ngdc_92m_bathy
help ndims
quiver_field(stn,[],[],[],'fkeys_hycom_seatemp_field',[17:2:31])
quiver_field(stn,'default',[],[],'fkeys_hycom_seatemp_field',[17:2:31])
dbstop quiver_field
quiver_field(stn,'default',[],[],'fkeys_hycom_seatemp_field',[17:2:31])
bathcntrs
size(squeeze(nanmean(stn.(bathfld).field,1))
size(squeeze(nanmean(stn.(bathfld).field,1)))
squeeze(nanmean(stn.(bathfld).field,1))
dbquit
quiver_field(stn)
dbquit
dbstatus
dbclear all
quiver_field(stn)
contour(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(nanmean(stn.fkeys_hycom_seatemp_field.field,1)),[25:0.2:27]);
quiver_field(stn)
[cs,ch]=contour(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(nanmean(stn.fkeys_hycom_seatemp_field.field,1)),[25:0.2:27]); clabel(cs,ch);
quiver_field(stn,'default',[],[],'none')
contour(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(nanmean(stn.fkeys_hycom_seatemp_field.field,1)),[25:0.2:27]);
quiver_field(stn,'default',[],[],'none')
[cs,ch]=contour(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(nanmean(stn.fkeys_hycom_seatemp_field.field,1)),[25:0.2:27]); clabel(cs,ch);
help map_freef
quiver_field(stn,'default',[],[],'default',[-10:-10:-50])
quiver_field(stn,'default',[],[],'default',[-5 -10:-10:-40 -100])
quiver_field(stn,'default',[],[],'default',[-5 -10:-10:-40 -80])
contour(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(iqr(stn.fkeys_hycom_seatemp_field.field,1)));
quiver_field(stn,'default',[],[],'default',[-5 -10:-10:-40 -80])
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(iqr(stn.fkeys_hycom_seatemp_field.field,1)));
nansummary(iqr(stn.fkeys_hycom_seatemp_field.field))
nansummary(std(stn.fkeys_hycom_seatemp_field.field))
quiver_field(stn,'default',[],[],'default',[-5 -10:-10:-40 -80])
quiver_field(stn)
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data<-0.5,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data<-0.5&stn.fkeys_hycom_v.data>0.5,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data<-0.3&stn.fkeys_hycom_v.data>0.3,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data<-0.2&stn.fkeys_hycom_v.data>0.2,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data<-0.1&stn.fkeys_hycom_v.data>0.2,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data>0.5&stn.fkeys_hycom_v.data>0.5,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data>0.5&stn.fkeys_hycom_v.data<0.5,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data>0.2&stn.fkeys_hycom_v.data<0.2,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data>0.1&stn.fkeys_hycom_v.data<0.1,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data>0.5&stn.fkeys_hycom_v.data<-0.5,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data>0.3&stn.fkeys_hycom_v.data<-0.3,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data>0.2&stn.fkeys_hycom_v.data<-0.3,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data>0.2&stn.fkeys_hycom_v.data<-0.2,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data>0.1&stn.fkeys_hycom_v.data<-0.2,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data<-0.2&stn.fkeys_hycom_v.data>0.2,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data<-0.1&stn.fkeys_hycom_v.data>0.2,1)))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data<-0.1&stn.fkeys_hycom_v.data<-0.2,1)))
print('-dpng','../figs/mlrf1-fkeys_hycom-eddy-example');
dir ../figs/mlrf1-fkeys_hycom-eddy-example.png
datestr(stn.fkeys_hycom_seatemp_field.date(54))
datestr(stn.fkeys_hycom_seatemp_field.date(50))
quiver_field(stn,stn.fkeys_hycom_u_field.date(find(stn.fkeys_hycom_u.data<-0.1&stn.fkeys_hycom_v.data<-0.2,1)))
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(stn.fkeys_hycom_seatemp_field.field(50,:,:)),[
adsf
])
nansummary(stn.fkeys_hycom_seatemp_field.field(50,:,:))
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(stn.fkeys_hycom_seatemp_field.field(50,:,:)),[20:0.2:24]);
quiver_field(stn,datenum(2004,1,13,12,0,0))
quiver_field(stn,datenum(2004,1,13,12,0,0),[],[],'default',[-5 -10:-40 -80))
quiver_field(stn,datenum(2004,1,13,12,0,0),[],[],'default',[-5 -10:-40 -80])
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(stn.fkeys_hycom_seatemp_field.field(50,:,:)),[20:0.2:24]);
minlon = min(stn.(ufld).lon(:));
maxlon = max(stn.(ufld).lon(:));
dlon = min(diff(unique(stn.(ufld).lon(:))));
minlat = min(stn.(ufld).lat(:));
maxlat = max(stn.(ufld).lat(:));
dlat = min(diff(unique(stn.(ufld).lat(:))));
fmg; map_freef([-80.49,-80.27,24.93,25.11], [-5 -10:-10:-40 -80]);
fmg; map_freef([-80.49,-80.27,24.93,25.11]);
contourf(stn.fkeys_hycom_seatemp_field.lon,stn.fkeys_hycom_seatemp_field.lat,squeeze(stn.fkeys_hycom_seatemp_field.field(50,:,:)),[20:0.2:24]);
print('-dpng','../figs/mlrf1-fkeys_hycom-seatemp-eddy-example');
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('mlrf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('mlrf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
print('-dtiff','../figs/mlrf1-ms1_clim-ndbc_sea_t_erai_erai_30a_erai_avhrr_hc_dTdt-with-even-better-Ktheta.tiff');
ylim([22,31])
ylim([22.5,30.5])
ylim([22,30.5])
ylim([22.3,30.3])
ylim([22,30.5])
ylim([22.3,30.5])
print('-dtiff','../figs/mlrf1-ms1_clim-ndbc_sea_t_erai_erai_30a_erai_avhrr_hc_dTdt-with-even-better-Ktheta.tiff');
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('mlrf1','erai','avhrr_weekly_sst','erai',[],'erai','hourly_misst_sst');
stn = optim_q0('mlrf1','erai','none','erai',[],'erai','hourly_misst_sst');
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('mlrf1','erai','none','erai',[],'erai','hourly_misst_sst');
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('mlrf1','erai','none','erai',[],'erai','hourly_misst_sst');
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('mlrf1','erai','none','erai',[],'erai','hourly_misst_sst');
help annsubs
annsubs(stn,'erai','none','erai',[],'erai','hourly_misst_sst',6:8,2002,2010,'ndbc_sea_t');
annsubs(stn,'erai','none','erai',[],'erai','hourly_misst_sst',[],6:8,2002,2010,'ndbc_sea_t');
plot_ts(stn.misst_sst)
plot_ts(stn.misst_sst,'r')
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('mlrf1','erai','gom_hycom','erai',[],'erai','ndbc_sea_t');
stn=[]; ans=[]; clear stn ans
stn = optim_q0('lonf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
ylim('default')
stn=[]; ans=[]; clear stn ans
stn = optim_q0('lonf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
ylim('default')
stn.station_name
91+45
182+45
274+45
datenum(2011,1,1)+319-1
datestr(datenum(2011,1,1)+319-1)
datestr(datenum(2011,1,1)+227-1)
datestr(datenum(2011,1,1)+136-1)
datestr(datenum(2011,1,1)+45-1)
stn=[]; ans=[]; clear stn ans
stn = optim_q0('lonf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
for ix=3:8; figure(ix); ylim([10,40]); end;
for ix=3:8; figure(ix); ylim([15,50]); end;
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
%-- 8/26/11 11:48 AM --%
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('lonf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
for ix=1:3; figure(ix); ylim([15,50]); end;
for ix=1:3; figure(ix); ylim([17,32]); end;
for ix=1:3; figure(ix); ylim([17-3,32+3]); end;
for ix=1:3; figure(ix); ylim([17-6,32+6]); end;
for ix=1:3; figure(ix); ylim([15,50]); end;
stn=[]; ans=[]; clear stn ans
stn = optim_q0('lonf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
stn=[]; ans=[]; clear stn ans
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
type get_ngdc_bathy_station.m
pd
rawxyz = load('data/coast/LGramer1-80.mat','lon','lat','depth');
rawxyz = load('coast/LGramer1-80.mat','lon','lat','depth');
stns = get_fknms_thermistors();
stns
stnms=grepstruct(stns,'^FKNMS_'); for ix=1:numel(stnms); stnm=stnms{ix}; stn=get_ngdc_bathy_station(stns.(stnm),[],rawxyz); stn=[]; clear stn; end;
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stns = get_fknms_thermistors();
stns
FKNMS_KW_CHANL: [1x1 struct]
rawxyz = load('coast/LGramer1-80.mat','lon','lat','depth');
grepstruct(stns,'CHE')
grepstruct(stns,'BROA')
stn = get_ngdc_bathy_station('FKNMS_BROAD_CRK');
stn
fmg; [cs,ch]=contour(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,[-2:-2:-20]); clabel(cs,ch);
stn=[]; ans=[]; clear stn ans
cs=[]; clear cs
stn = stns.FKNMS_BROAD_CRK;
stn = get_ngdc_bathy_station(stn);
fmg; [cs,ch]=contour(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,[-2:-2:-20]); clabel(cs,ch); plot(stn.lon,stn.lat,'kp');
stn
stn.ngdc_92m_bathy
rix=400:600; cix=350:550; fmg; [cs,ch]=contour(stn.ngdc_92m_bathy.lon(rix,cix),stn.ngdc_92m_bathy.lat(rix,cix),stn.ngdc_92m_bathy.field(rix,cix),[-2:-2:-20]); clabel(cs,ch); plot(stn.lon,stn.lat,'kp');
rix=425:525; cix=390:470; fmg; [cs,ch]=contour(stn.ngdc_92m_bathy.lon(rix,cix),stn.ngdc_92m_bathy.lat(rix,cix),stn.ngdc_92m_bathy.field(rix,cix),[0:-1:-20]); clabel(cs,ch); plot(stn.lon,stn.lat,'kp');
help clabel
rix=425:525; cix=390:470; fmg; [cs,ch]=contour(stn.ngdc_92m_bathy.lon(rix,cix),stn.ngdc_92m_bathy.lat(rix,cix),stn.ngdc_92m_bathy.field(rix,cix),[0:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp');
rix=400:525; cix=370:470; fmg; [cs,ch]=contour(stn.ngdc_92m_bathy.lon(rix,cix),stn.ngdc_92m_bathy.lat(rix,cix),stn.ngdc_92m_bathy.field(rix,cix),[0:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp');
rix=400:525; cix=390:500; fmg; [cs,ch]=contour(stn.ngdc_92m_bathy.lon(rix,cix),stn.ngdc_92m_bathy.lat(rix,cix),stn.ngdc_92m_bathy.field(rix,cix),[0:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp');
rix=400:525; cix=390:500; lons=stn.ngdc_92m_bathy.lon(rix,cix); lats=stn.ngdc_92m_bathy.lat(rix,cix); fmg; map_freef([min(lons),max(lons),min(lats),max(lats)],'none'); [cs,ch]=contour(lons,lats,stn.ngdc_92m_bathy.field(rix,cix),[0:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp');
help map_freef
[min(lons),max(lons),min(lats),max(lats)]
rix=400:525; cix=390:500; lons=stn.ngdc_92m_bathy.lon(rix,cix); lats=stn.ngdc_92m_bathy.lat(rix,cix); fmg; map_freef([min(lons(:)),max(lons(:)),min(lats(:)),max(lats(:))],'none'); [cs,ch]=contour(lons,lats,stn.ngdc_92m_bathy.field(rix,cix),[0:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp');
find_date_ranges(stn.fknms_seatemp.date,1)
find_date_ranges(stn.fknms_seatemp.date,10)
find_date_ranges(stn.fknms_seatemp.date,5)
rix=400:525; cix=390:500; lons=stn.ngdc_92m_bathy.lon(rix,cix); lats=stn.ngdc_92m_bathy.lat(rix,cix); fmg; map_freef([min(lons(:)),max(lons(:)),min(lats(:)),max(lats(:))],'none'); [cs,ch]=contour(lons,lats,stn.ngdc_92m_bathy.field(rix,cix),[-1:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp');
daspect([1,1,1])
rix=430:500; cix=400:470; lons=stn.ngdc_92m_bathy.lon(rix,cix); lats=stn.ngdc_92m_bathy.lat(rix,cix); fmg; map_freef([min(lons(:)),max(lons(:)),min(lats(:)),max(lats(:))],'none'); [cs,ch]=contour(lons,lats,stn.ngdc_92m_bathy.field(rix,cix),[-1:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp');
rix=430:500; cix=400:470; lons=stn.ngdc_92m_bathy.lon(rix,cix); lats=stn.ngdc_92m_bathy.lat(rix,cix); fmg; map_freef([min(lons(:)),max(lons(:)),min(lats(:)),max(lats(:))],'none'); [cs,ch]=contour(lons,lats,stn.ngdc_92m_bathy.field(rix,cix),[-2:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp');
stn=[]; ans=[]; clear stn ans
stn = stns.FKNMS_HEN_CHIX;
stn = get_ngdc_bathy_station(stn);
rix=430:500; cix=400:470; lons=stn.ngdc_92m_bathy.lon(rix,cix); lats=stn.ngdc_92m_bathy.lat(rix,cix); fmg; map_freef([min(lons(:)),max(lons(:)),min(lats(:)),max(lats(:))],'none'); [cs,ch]=contour(lons,lats,stn.ngdc_92m_bathy.field(rix,cix),[-2:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp'); titlename('Hens & Chickens bathymetry');
print('-dtiff','figs/fknms_hen_chix_bathy');
stn=[]; ans=[]; clear stn ans
stn = stns.FKNMS_BROAD_CRK;
stn = get_ngdc_bathy_station(stn);
rix=430:500; cix=400:470; lons=stn.ngdc_92m_bathy.lon(rix,cix); lats=stn.ngdc_92m_bathy.lat(rix,cix); fmg; map_freef([min(lons(:)),max(lons(:)),min(lats(:)),max(lats(:))],'none'); [cs,ch]=contour(lons,lats,stn.ngdc_92m_bathy.field(rix,cix),[-2:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp'); titlename('Broad Creek bathymetry');
print('-dtiff','figs/fknms_broad_crk_bathy');
stn=[]; ans=[]; clear stn ans
stnms=grepstruct(stns,'^FKNMS_[K-Z]'); for ix=1:numel(stnms); stnm=stnms{ix}; stn=[]; clear stn; end;
stnms=grepstruct(stns,'^FKNMS_[K-Z]'); for ix=1:numel(stnms); stnm=stnms{ix}; disp(stnm); stn=[]; clear stn; end;
stnms=grepstruct(stns,'^FKNMS_[K-Z]'); for ix=1:numel(stnms); stnm=stnms{ix}; disp(stnm); stn=get_ngdc_bathy_station(stns.(stnm),[],rawxyz); stn=[]; clear stn; end;
stns=[]; ans=[]; clear stns ans
stn.station_name = 'cheeca'; stn.lon=-80.61475; stn.lat=24.904083;
stn = get_ngdc_bathy_station(stn,[],rawxyz);
edit ancheeca.m
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
ancheeca
rix=430:500; cix=400:470; lons=stn.ngdc_92m_bathy.lon(rix,cix); lats=stn.ngdc_92m_bathy.lat(rix,cix); fmg; map_freef([min(lons(:)),max(lons(:)),min(lats(:)),max(lats(:))],'none'); [cs,ch]=contour(lons,lats,stn.ngdc_92m_bathy.field(rix,cix),[-2:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp'); titlename('Cheeca Rocks bathymetry');
edit ancheeca.m
print('-dtiff','figs/cheeca_bathy');
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
%-- 8/26/11  3:45 PM --%
pd
fclose('all'); close all; clear all; clear classes; clear java; clear mex; dbclear all; more on; pack
stn = optim_q0('lonf1','erai','avhrr_weekly_sst','erai',[],'erai','ndbc_sea_t');
9238/60
stn=[]; ans=[]; clear stn ans; pack
stn = get_looe1_adcp;
fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date,stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.profile); clabel(cs,ch);
fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date,stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.prof); clabel(cs,ch);
size(stn.adcp_baroclinic_u_40hlp.date),size(stn.adcp_bin_heights),size(stn.adcp_baroclinic_u_40hlp.prof)
fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date,stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.prof'); clabel(cs,ch);
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.prof(1:strd:end,:)'); clabel(cs,ch);
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.prof(1:strd:end,:)'); colorbar;
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.prof(1:strd:end,:)'); datetick3; colorbar;
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_u_40hlp.prof(1:strd:end,:)'); datetick3; colorbar;
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_u_40hlp.prof(1:strd:end,:)',[-.1:0.02:.1]); datetick3; colorbar;
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_u_40hlp.prof(1:strd:end,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([0,23]);
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_u_40hlp.prof(1:strd:end,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([2,22]);
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_baroclinic_u_40hlp.prof(1:strd:end,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([2,22]); titlename('Baroclinic u^.\nablah 40hlp');
x = load_all_ndbc_data([],'smkf1');
x = verify_variable(x,'ndbc_air_t_3_day_lowpass');
help plot_ts
lh = plot_ts(x.ndbc_air_t_3_day_lowpass,'k-','LineWidth',3)
titlename('u^.\nablah 40hlp'); grid on;
ylim([2,25])
lh = plot_ts(x.ndbc_air_t,'k-','LineWidth',3)
lh = plot_ts(stn.adcp_seatemp,'r');
ylim([2,28])
ylim([2,27])
stn = get_looe1_microcat(stn);
lh = plot_ts(stn.microcat_seatemp,'r');
ylim([2,33])
ylim([2,29])
ylim([2,33])
x = verify_variable(x,'ndbc_wind1_speed_3_day_lowpass');
x = verify_variable(x,'ndbc_wind1_u_3_day_lowpass');
x = verify_variable(x,'ndbc_wind1_v_3_day_lowpass');
lh = plot_ts(x.ndbc_wind1_v_3_day_lowpass,'m');
ylim([-15,33])
x = verify_variable(x,'ndbc_wind1_v_1_day_lowpass');
lh = plot_ts(x.ndbc_wind1_v_1_day_lowpass,'m');
ylim([-20,33])
strd=40; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u_40hlp.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_baroclinic_v_40hlp.prof(1:strd:end,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([2,22]); titlename('Baroclinic u^.(\nablah=0) 40hlp');
ylim([-20,33])
lh = plot_ts(x.ndbc_air_t,'k-','LineWidth',3)
lh = plot_ts(stn.microcat_seatemp,'r');
lh = plot_ts(x.ndbc_wind1_v_1_day_lowpass,'m');
lh = plot_ts(x.ndbc_wind1_u_1_day_lowpass,'g');
x = verify_variable(x,'ndbc_wind1_u_1_day_lowpass');
lh = plot_ts(x.ndbc_wind1_u_1_day_lowpass,'g');
quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:24:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:24:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:24:end),x.ndbc_wind1_v_1_day_lowpass.data(1:24:end));
help quiver
more on
help quiver
quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:24:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:24:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:24:end),x.ndbc_wind1_v_1_day_lowpass.data(1:24:end),0,'k');
quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:24:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:24:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:24:end),x.ndbc_wind1_v_1_day_lowpass.data(1:24:end),0,'kx');
strd=1; quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:strd:end),x.ndbc_wind1_v_1_day_lowpass.data(1:strd:end),0,'kx');
strd=8; quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:strd:end),x.ndbc_wind1_v_1_day_lowpass.data(1:strd:end),0,'kx');
strd=6; quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:strd:end),x.ndbc_wind1_v_1_day_lowpass.data(1:strd:end),0,'kx');
strd=6; quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:strd:end),x.ndbc_wind1_v_1_day_lowpass.data(1:strd:end),0.1,'kx');
strd=6; quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:strd:end),x.ndbc_wind1_v_1_day_lowpass.data(1:strd:end),0.2,'kx');
strd=3; quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:strd:end),x.ndbc_wind1_v_1_day_lowpass.data(1:strd:end),0.2,'kx');
lh = plot_ts(stn.adcp_seatemp,'b');
strd=3; quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:strd:end),x.ndbc_wind1_v_1_day_lowpass.data(1:strd:end),0.1,'kx');
strd=12; quiver(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u_1_day_lowpass.date(1:strd:end))),x.ndbc_wind1_u_1_day_lowpass.data(1:strd:end),x.ndbc_wind1_v_1_day_lowpass.data(1:strd:end),0.1,'kx');
strd=1; quiver(x.ndbc_wind1_u.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u.date(1:strd:end))),x.ndbc_wind1_u.data(1:strd:end),x.ndbc_wind1_v.data(1:strd:end),0.1,'kx');
strd=1; quiver(x.ndbc_wind1_u.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u.date(1:strd:end))),x.ndbc_wind1_u.data(1:strd:end),x.ndbc_wind1_v.data(1:strd:end),0.2,'kx');
ylim([-4,28])
colorbar
strd=6; quiver(x.ndbc_wind1_u.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u.date(1:strd:end))),x.ndbc_wind1_u.data(1:strd:end),x.ndbc_wind1_v.data(1:strd:end),0.3,'kx');
lh=plot_ts(x.ndbc_wind1_speed)
stn = get_erai_station(stn);
lh = plot_ts(stn.erai_air_t,'k:','LineWidth',2)
lh = plot_ts(stn.erai_air_t,'k:','LineWidth',4)
strd=3; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_baroclinic_u.prof(1:strd:end,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([2,22]); titlename('Baroclinic u^.\nablah');
strd=3; fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u.date(1:strd:end),stn.adcp_bin_heights,stn.adcp_baroclinic_u.prof(1:strd:end,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([2,22]); xlim(datenum(2008,12,[9,27])); titlename('Baroclinic u^.\nablah');
ix=find(datenum(2008,12,9)<=stn.adcp_baroclinic_u.date&stn.adcp_baroclinic_u.date<=datenum(2008,12,27));
strd=ix(1:3:end); fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u.date(strd),stn.adcp_bin_heights,stn.adcp_baroclinic_u.prof(strd,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([2,22]); titlename('Baroclinic u^.\nablah');
ylim([-4,28])
strd=ix(1:1:end); fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u.date(strd),stn.adcp_bin_heights,stn.adcp_baroclinic_u.prof(strd,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([2,22]); titlename('Baroclinic u^.\nablah');
ylim([-4,28])
strd=ix(1:1:end); fmg; [cs,ch]=contourf(stn.adcp_u.date(strd),stn.adcp_bin_heights,stn.adcp_u.prof(strd,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([-4,28]); titlename('Raw u^.\nablah');
set_hovmuller_cursor
lh = plot_ts(stn.erai_air_t,'k:','LineWidth',4)
lh=plot_ts(x.ndbc_wind1_speed)
lh = plot_ts(stn.microcat_seatemp,'r');
lh = plot_ts(stn.adcp_seatemp,'b');
strd=6; quiver(x.ndbc_wind1_u.date(1:strd:end),repmat(0,size(x.ndbc_wind1_u.date(1:strd:end))),x.ndbc_wind1_u.data(1:strd:end),x.ndbc_wind1_v.data(1:strd:end),0.3,'kx');
set_datetick_cursor
oi
pd
oi = load('data/mlrf1_oisst2.mat');
oi
oi.station
oi.station.oisst2_sst
find_date_ranges(oi.station.oisst2_sst.date)
oi.station.oisst2_sst.date([1 end])
datestr(oi.station.oisst2_sst.date(1))
numel(find(isfinite(oi.station.oisst2_sst.date)))
find_date_ranges(oi.station.oisst2_sst.date)
find_date_ranges(oi.station.oisst2_sst.date(isfinite(oi.station.oisst2_sst.date)))
min(diff(oi.station.oisst2_sst.date))
oi=[]; ans=[]; clear oi ans
oi = load('data/smkf1_oisst2.mat');
fmg; plot_ts(x.ndbc_sea_t,oi.station.oisst2_sst);
x = get_misst_station(x);
fmg; plot_ts(x.ndbc_sea_t,oi.station.oisst2_sst,x.misst_sst); legend('In situ','OISST v2','MISST');
oi=[]; ans=[]; clear oi ans
dir *oi*m
type read_oisst2
more on
type read_oisst2
dbstop read_oisst2
mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1); mlrf1 = get_misst_station(mlrf1); mlrf1 = read_oisst2(mlrf1);
stations
allfldsd
allflds
dbstatus
dbclear all
fmg; plot_ts(mlrf1.ndbc_sea_t,mlrf1.oisst2_sst,mlrf1.misst_sst); legend('In situ','OISST v2','MISST');
fmg; plot_ts(mlrf1.ndbc_sea_t,'-',mlrf1.mlrf1.misst_sst,'-',oisst2_sst,'-'); legend('In situ','OISST v2','MISST');
fmg; plot_ts(mlrf1.ndbc_sea_t,'-',mlrf1.misst_sst,'-',mlrf1.oisst2_sst,'-'); legend('In situ','OISST v2','MISST');
fmg; plot_ts(mlrf1.ndbc_sea_t,'-',mlrf1.misst_sst,'-',mlrf1.oisst2_sst,'-'); legend('In situ','MISST','OISST v2');
datetick3('x',2,'keeplimits')
mlrf1 = verify_variable(mlrf1,'ndbc_sea_t_1_day_average');
fmg; plot_ts(mlrf1.ndbc_sea_t_1_day_average,'-',mlrf1.misst_sst,'-',mlrf1.oisst2_sst,'-'); legend('In situ','MISST','OISST v2');
mlrf1 = verify_variable(mlrf1,'ndbc_sea_t_30_day_average');
mlrf1 = verify_variable(mlrf1,'oisst2_sst_30_day_average');
mlrf1 = verify_variable(mlrf1,'misst_sst_30_day_average');
fmg; plot_ts(mlrf1.ndbc_sea_t_30_day_average,'-',mlrf1.misst_sst_30_day_average,'-',mlrf1.oisst2_sst_30_day_average,'-'); legend('In situ','MISST','OISST v2');
annotline([],30.4,'MMM','red')
annotline([],31.4,'MMM+1','yellow')
delete(ans)
prctile(mlrf1.ndbc_sea_t_30_day_average.data,99.3)
prctile(mlrf1.ndbc_sea_t_30_day_average.data,99.2)
prctile(mlrf1.ndbc_sea_t_30_day_average.data,99)
prctile(mlrf1.misst_sst_30_day_average.data,99)
prctile(mlrf1.oisst2_sst_30_day_average.data,99)
mmmi=prctile(mlrf1.ndbc_sea_t_30_day_average.data,99)
mmmm=prctile(mlrf1.misst_sst_30_day_average.data,99)
mmmo=prctile(mlrf1.oisst2_sst_30_day_average.data,99)
annotline([],mmmi,'In situ','blue')
help  annotline
annotline([],mmmm,'MISST',[.2 .6 .2])
annotline([],mmmo,'OISSTv2','red')
datetick3('x',2,'keeplimits')
prctile(mlrf1.ndbc_sea_t_30_day_average.data,98.9)
prctile(mlrf1.ndbc_sea_t_30_day_average.data,98.7)
prctile(mlrf1.ndbc_sea_t_30_day_average.data,98.8)
mmmo=prctile(mlrf1.oisst2_sst_30_day_average.data,98.8)
mmmi=prctile(mlrf1.ndbc_sea_t_30_day_average.data,98.8)
mmmm=prctile(mlrf1.misst_sst_30_day_average.data,98.8)
annotline([],mmmi,'In situ','blue'); annotline([],mmmm,'MISST',[.2 .6 .2]); annotline([],mmmo,'OISSTv2','red');
fmg; plot_ts(mlrf1.ndbc_sea_t_30_day_average,'-',mlrf1.misst_sst_30_day_average,'-',mlrf1.oisst2_sst_30_day_average,'-'); legend('In situ','MISST','OISST v2');
mmmi=prctile(mlrf1.ndbc_sea_t_30_day_average.data,98.8); mmmm=prctile(mlrf1.misst_sst_30_day_average.data,98.8); mmmo=prctile(mlrf1.oisst2_sst_30_day_average.data,98.8); annotline([],mmmi,'In situ','blue'); annotline([],mmmm,'MISST',[.2 .6 .2]); annotline([],mmmo,'OISSTv2','red');
help commandline
lookfor command
lookfor history
edit anoi.m
mlrf1=[]; ans=[]; clear mlrf1 ans
smkf1 = anoi('smkf1')
smkf1=[]; ans=[]; clear smkf1 ans
sanf1 = anoi('sanf1')
datetick3('x',2,'keeplimits')
smkf1 = anoi('smkf1');
prctile(smkf1.ndbc_sea_t_30_day_average.data,98)
mlrf1 = anoi('mlrf1');
prctile(mlrkf1.ndbc_sea_t_30_day_average.data,98)
prctile(mlrf1.ndbc_sea_t_30_day_average.data,98)
datetick3('x',2,'keeplimits')
prctile(mlrf1.ndbc_sea_t_30_day_average.data(ismember(get_year(),[1992:2006])),98)
prctile(mlrf1.ndbc_sea_t_30_day_average.data(ismember(get_year(mlrf1.ndbc_sea_t_30_day_average.date),[1992:2006])),98)
prctile(mlrf1.ndbc_sea_t_30_day_average.data(ismember(get_year(mlrf1.ndbc_sea_t_30_day_average.date),[1992:2006])),99)
prctile(mlrf1.ndbc_sea_t_30_day_average.data(ismember(get_year(mlrf1.ndbc_sea_t_30_day_average.date),[1992:2006])),98.5)
prctile(smkf1.ndbc_sea_t_30_day_average.data(ismember(get_year(smkf1.ndbc_sea_t_30_day_average.date),[1992:2006])),98.5)
prctile(smkf1.ndbc_sea_t_30_day_average.data(ismember(get_year(smkf1.ndbc_sea_t_30_day_average.date),[1992:2006])),98.4)
prctile(smkf1.ndbc_sea_t_30_day_average.data(ismember(get_year(smkf1.ndbc_sea_t_30_day_average.date),[1992:2006])),98.3)
prctile(mlrf1.ndbc_sea_t_30_day_average.data(ismember(get_year(mlrf1.ndbc_sea_t_30_day_average.date),[1992:2006])),98.3)
mlrf1=[]; ans=[]; clear mlrf1 ans
smkf1=[]; ans=[]; clear smkf1 ans
sanf1=[]; ans=[]; clear sanf1 ans
mlrf1 = anoi('mlrf1');
smkf1 = anoi('smkf1');
edit anoi.m
mlrf1 = anoi('mlrf1');
smkf1 = anoi('smkf1');
sanf1 = anoi('sanf1')
set_hovmuller_cursor
set_datetick_cursor
set_hovmuller_cursor
a
help read_avhrr_subset
help read_avhrr_png
type read_avhrr_png
sst = read_avhrr_png('http://optics.marine.usf.edu/subscription/modis/FLKEYS/2011/daily/239/T20112390325.QKM.FLKEYS.PASS.L3D.SST4.png');
sstbytes = imread(('http://optics.marine.usf.edu/subscription/modis/FLKEYS/2011/daily/239/T20112390325.QKM.FLKEYS.PASS.L3D.SST4.png');
sstbytes = imread('http://optics.marine.usf.edu/subscription/modis/FLKEYS/2011/daily/239/T20112390325.QKM.FLKEYS.PASS.L3D.SST4.png');
sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;
sst(sstbytes >= 251) = nan;
minval = +3.0; maxval = +35.0;
sst(minval > sst | sst > maxval) = nan;
sstbytes = []; clear sstbytes
fmg; contourf(sst); colorbar;
sstbytes = imread('http://optics.marine.usf.edu/subscription/modis/FLKEYS/2011/daily/239/T20112390325.QKM.FLKEYS.PASS.L3D.SST4.png');
sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;
fmg; contourf(sst); colorbar;
fmg; contourf(sst,[29:.5:35]); colorbar;
sstbytes=[]; sst=[]; clear sstbytes sst
sstbytes = imread('http://optics.marine.usf.edu/subscription/modis/FLKEYS/2011/daily/239/T20112390325.QKM.FLKEYS.PASS.L3D.SST.png');
sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;
sst(sstbytes >= 251) = nan;
fmg; contourf(sst,[29:.5:35]); colorbar;
strd=ix(1:1:end); fmg; [cs,ch]=contourf(stn.adcp_baroclinic_u.date(strd),stn.adcp_bin_heights,stn.adcp_baroclinic_u.prof(strd,:)',[-.1:0.02:.1]); datetick3; colorbar; ylim([-4,28]); titlename('Baroclinic u^.\nablah');



sst = (32.0 * cast(sstbytes,'double') / 235.0) + 10.0;
