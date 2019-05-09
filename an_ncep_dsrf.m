1;

%-- Unknown date --%
help find
more on
help find
[tix,hix,q0ix]=intersect_all_dates(((30+(0.005/60))/(24*60)),stn.ndbc_sea_t.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date,stn.tmd_tide_i_depth.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates(((30+(0.005/60))/(24*60)),stn.tmd_tide_i_depth.date,stn.ndbc_sea_t.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates(((30+(1/60))/(24*60)),stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
datestr(stn.ndbc_sea_t.date(1)),datestr(stn.tmd_tide_i_depth.date(1))
datestr(stn.ndbc_sea_t.date(1)),datestr(stn.tmd_tide_i_depth.date(2))
(0.00/60.0)==0.0
[tix,hix,q0ix]=intersect_all_dates(((30+(1/60))/(24*60)),stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates(((30+(1/60))/(24*60)),stn.tmd_tide_i_depth.date,stn.ndbc_sea_t.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.tmd_tide_i_depth.date,stn.ndbc_sea_t.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date,stn.tmd_tide_i_depth.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date,stn.tmd_tide_i_depth.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.tmd_tide_i_depth.date,stn.ndbc_sea_t.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.tmd_tide_i_depth.date,stn.ndbc_sea_t.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
length(find(diff(stn.ndbc_sea_t.date<(1/24))))
length(find(diff(stn.ndbc_sea_t.date<(0/24))))
length(find(diff(stn.ndbc_sea_t.date<(2/24))))
length(find(diff(stn.ndbc_sea_t.date)<(1/24)))
length(find(diff(stn.ndbc_sea_t.date)<(5/24)))
length(find(diff(stn.ndbc_sea_t.date)<(24/24)))
length(find(diff(stn.ndbc_sea_t.date(:))<(5/24)))
length(find(diff(stn.ndbc_sea_t.date(:))<(0.9/24)))
length(find(diff(stn.ndbc_sea_t.date(:))<0))
length(find(diff(stn.ndbc_sea_t.date(:))<0.0000001))
length(find(diff(stn.ndbc_sea_t.date(:))<0.00001))
length(find(diff(stn.ndbc_sea_t.date(:))<0.0001))
length(find(diff(stn.ndbc_sea_t.date(:))<0.001))
length(find(diff(stn.ndbc_sea_t.date(:))<0.01))
length(find(diff(stn.ndbc_sea_t.date(:))<0.015))
length(find(diff(stn.ndbc_sea_t.date(:))<0.15))
length(find(diff(stn.ndbc_sea_t.date(:))<0.10))
length(find(diff(stn.ndbc_sea_t.date(:))<0.05))
length(find(diff(stn.ndbc_sea_t.date(:))<0.03))
length(find(diff(stn.ndbc_sea_t.date(:))<0.01))
length(find(diff(stn.ndbc_sea_t.date(:))<0.02))
1/24
length(find(diff(stn.ndbc_sea_t.date(:))<(2/24)))
length(find(diff(stn.ndbc_sea_t.date(:))<(0.5/24)))
length(find(diff(stn.ndbc_sea_t.date(:))<(0.4/24)))
length(find(diff(stn.ndbc_sea_t.date(:))<(0.49/24)))
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_dt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_qedt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
stn = verify_variable(stn,'gom_hycom_qedt_heat_flux_24_hour_average');
[tix,hix,q0ix]=intersect_all_dates([],stn.ndbc_sea_t.date,stn.tmd_tide_i_depth.date,stn.gom_hycom_qedt_heat_flux_24_hour_average.date); size(tix),size(hix),size(q0ix)
grepstruct(stn,'diag')
figure; maxigraph; plot(stn.ndbc_ncep_30a_cordiags.date,[stn.ndbc_ncep_30a_cordiags.ustar,stn.ndbc_ncep_30a_cordiags.tstar,stn.ndbc_ncep_30a_cordiags.qstar,]);
figure; maxigraph; plot(stn.ndbc_ncep_30a_cordiags.date,[stn.ndbc_ncep_30a_cordiags.ustar,stn.ndbc_ncep_30a_cordiags.tstar,stn.ndbc_ncep_30a_cordiags.qstar,]); datetick3; legend('u_*','\theta_*','q_*', 'Location','Best');
figure; hist(stn.ndbc_ncep_30a_cordiags.ustar,1000);
figure; hist(stn.ndbc_ncep_30a_cordiags.qstar,1000);
figure; hist((stn.ndbc_ncep_30a_cordiags.tstar.*stn.ndbc_ncep_30a_cordiags.qstar),1000);
figure; hist((1.2*9.8*stn.ndbc_ncep_30a_cordiags.ustar.*stn.ndbc_ncep_30a_cordiags.tstar),1000);
figure; hist((0.61*9.8*stn.ndbc_ncep_30a_cordiags.ustar.*stn.ndbc_ncep_30a_cordiags.qstar),1000);
figure; hist(((1.2*9.8*stn.ndbc_ncep_30a_cordiags.ustar.*stn.ndbc_ncep_30a_cordiags.tstar)./stn.ndbc_air_t.data),1000);
[tix,aix]=intersect_dates(stn.ndbc_ncep_30a_cordiags.date,stn.ndbc_air_t.date);
figure; hist(((1.2*9.8*stn.ndbc_ncep_30a_cordiags.ustar.*stn.ndbc_ncep_30a_cordiags.tstar)./stn.ndbc_air_t.data(aix)),1000);
figure; hist(((1.2*9.8*stn.ndbc_ncep_30a_cordiags.ustar(tix).*stn.ndbc_ncep_30a_cordiags.tstar(tix))./stn.ndbc_air_t.data(aix)),1000);
figure; hist(((1.2*9.8*stn.ndbc_ncep_30a_cordiags.ustar(tix).*stn.ndbc_ncep_30a_cordiags.tstar(tix))./1),1000);
figure; hist(((1.2*9.8*stn.ndbc_ncep_30a_cordiags.ustar(tix).*stn.ndbc_ncep_30a_cordiags.tstar(tix))./stn.ndbc_air_t.data(aix)),1000);
var(((1.2*9.8*stn.ndbc_ncep_30a_cordiags.ustar(tix).*stn.ndbc_ncep_30a_cordiags.tstar(tix))./stn.ndbc_air_t.data(aix))
var((1.2*9.8*stn.ndbc_ncep_30a_cordiags.ustar(tix).*stn.ndbc_ncep_30a_cordiags.tstar(tix))./stn.ndbc_air_t.data(aix))
std((1.2*9.8*stn.ndbc_ncep_30a_cordiags.ustar(tix).*stn.ndbc_ncep_30a_cordiags.tstar(tix))./stn.ndbc_air_t.data(aix))
figure; maxigraph; plot(stn.ndbc_ncep_30a_cordiags.date,[stn.ndbc_ncep_30a_cordiags.Cd,stn.ndbc_ncep_30a_cordiags.Ch,stn.ndbc_ncep_30a_cordiags.Ce,]); datetick3; legend('C_D','C_h','C_e', 'Location','Best');
%-- 12/27/10  8:33 PM --%
fclose('all'); close all; clear all; dbclear all; pack;
load('looe1_ms.mat')
load('../data/looe1_ms.mat')
station
stn = station; station=[]; clear station
fclose('all'); close all; clear all; dbclear all; pack;
stn = anlooe1([],[],1.20,1.20,0);
station=stn; save('../data/looe1_ms.dat','station'); station=[]; clear station
dbstop errbdgt
errbdgt
dts = stn.ndbc_ncep_30a_net_heat_flux.date(flix);
figure; maxigraph; plot(stn.ndbc_ncep_30a_sensible_heat_flux.date,stn.ndbc_ncep_30a_sensible_heat_flux.data,'k-'); plot(dts,sQsh2,'r.');
figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_sensible_heat_flux.date,stn.ndbc_ncep_30a_sensible_heat_flux.data,'k-'); plot(dts,sQsh2,'r.'); datetick3
station=stn; save('../data/looe1_ms.mat','station'); station=[]; clear station
max(diff(sort(stn.ndbc_air_t.data)))
min(diff(sort(stn.ndbc_air_t.data)))
min(diff(unique(stn.ndbc_air_t.data)))
min(diff(unique(stn.ndbc_wind1_speed.data)))
dbstatus
dbclear all
errbdgt
figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_sensible_heat_flux.date,stn.ndbc_ncep_30a_sensible_heat_flux.data,'k-'); plot(dts,sQsh2,'r.'); datetick3;
x = load_station_data('smkf1');
min(diff(unique(x.ndbc_air_t.data)))
min(diff(unique(x.air_t.data)))
min(diff(unique(x.sea_t.data)))
x
x=[]; clear x
stn
whos(stn)
whos stn
whos stn.adcp_u
help whos
s=whos('stn')
s=whos('stn.adcp_u')
help whos
help who
stn.adcp_u
35337*36
35337*36*2*4
grepstruct(stn,'adcp')
35337*36*2*6*2*3
35337*36*2*6*2*3*8
whos stn
help get_looe1_adcp
type get_looe1_adcp
min(diff(unique(stn.microcat_seatemp.data)))
max(diff(unique(stn.microcat_seatemp.data)))
type get_looe1_adcp
max(diff(unique(stn.microcat_seatemp.data)))
min(diff(unique(stn.adcp_seatemp.data)))
min(diff(unique(stn.microcat_seatemp.data)))
d=diff(unique(stn.microcat_seatemp.data));
min(d(d>eps))
min(d(d>1e-10))
min(d(d>1e-7))
x
x=load('../data/SFP/MC 143 Looe (Recov May 27 2009).asc',43,0);
x=csvread('../data/SFP/MC 143 Looe (Recov May 27 2009).asc',43,0);
help csvread
x=csvread('../data/SFP/MC 143 Looe (Recov May 27 2009).asc',0,0,[43 0 5895 0]);
x=csvread('../data/SFP/MC 143 Looe (Recov May 27 2009).asc',43,0,[43 0 5895 0]);
x=csvread('../data/SFP/foo.csv');
min(diff(unique(x)))
errbdgt
cov_ts(th,w)
[ig,r_th_w] = cov_ts(th,w);
r_th_w
scatter_fit_ts(stn.ncep_air_t,stn.ndbc_air_t)
scatter_fit_ts(stn.ncep_spechumid,stn.ndbc_spechumid)
x
clear x
x = load_all_ndbc_data([],'smkf1');
x = get_ncep_station(x,'narr');
scatter_fit_ts(stn.ncep_spechumid,x.ncep_spechumid)
scatter_fit_ts(x.ndbc_spechumid,x.ncep_spechumid)
scatter_fit_ts(stn.ndbc_spechumid,x.ncep_spechumid)
help station_dewp_to_relhumid
x=station_dewp_to_relhumid(x,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
help station_relhumid_to_spechumid
x=station_relhumid_to_spechumid(x,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
scatter_fit_ts(x.ndbc_spechumid,x.ncep_spechumid)
grepstruct(x,
grepstruct(x,'ncep')
scatter_fit_ts(x.ndbc_barom,x.ncep_barom)
scatter_fit_ts(stn.ncep_air_t,x.ncep_air_t)
figure; plot(stn.ncep_air_t.date,stn.ncep_air_t.data,x.ncep_air_t.date,x.ncep_air_t.data);
figure; plot(stn.ncep_air_t.date,stn.ncep_air_t.data,'b.',x.ncep_air_t.date,x.ncep_air_t.data,'g.'); datetick3;
figure; plot(stn.ncep_air_t.date,stn.ncep_air_t.data,'.',x.ncep_air_t.date,x.ncep_air_t.data,'.'); datetick3;
figure; maxigraph; plot(stn.ncep_air_t.date,stn.ncep_air_t.data,'.',x.ncep_air_t.date,x.ncep_air_t.data,'.'); datetick3; legend('LOOE1','SMKF1','Loc','Best'); titlename('T_a');
x
multiplot_station(x,{'ncep_dewp','ncep_dsrf','ncep_dlrf','ncep_barom','ncep_spechumid','ncep_cloud_cover','ncep_wind_speed'});
multiplot_station(x,{'ncep_dewp','ncep_dsrf','ncep_dlrf','ncep_barom','ncep_spechumid','ncep_cloud_cover','ncep_wind_speed','ncep_air_t','ncep_precip'});
scatter_fit_ts(x.ndbc_relhumid,x.ncep_relhumid)
scatter_fit_ts(x.ndbc_spechumid,x.ncep_spechumid)
x=[]; clear x
std(stn.ndbc_air_t.data)
[yr,mo,dy,hr,mn,sc]=datevec(stn.ndbc_air_t.date);
jd=datenum(yr,mo,dy)-datenum(yr,1,1)+1;
dh=datenum(1,mo,dy,hr);
dh=datenum(1,mo,dy,hr,0,0);
disp screeched
disp stretched
length(unique(dh))
udh=unique(dh);
for ix=1:length(udh); dhm(ix)=nanmean(stn.ndbc_air_t.data(dh==udh(ix))); end;
size(dhm)
for ix=1:length(udh); dhs(ix)=nanstd(stn.ndbc_air_t.data(dh==udh(ix))); end;
figure; maxigraph; plot(udh,[dhm,dhs]); datetick;
figure; maxigraph; plot(udh,[dhm;dhs]); datetick3;
datetick3('x',15
datetick3('x',15)
datetick3('x',16)
help datestr
datetick3('x',19)
help datestr
datetick3('x',6)
errbdgt
dts = stn.ndbc_ncep_30a_net_heat_flux.date(flix);
f.date=dts; f.data=sQsh2;
scatter_fit_ts(f,stn.ndbc_ncep_30a_sensible_heat_flux);
grepstruct(stn,'humid')
stn.ndbc_spechumid,stn.ndbc_sea_spechumid
vapor(nanmin(stn.ndbc_air_t.data))
vapor(nanmax(stn.ndbc_air_t.data))
vapor(nanmin(stn.ndbc_air_t.data))/vapor(nanmax(stn.ndbc_air_t.data))-
vapor(nanmin(stn.ndbc_air_t.data))/vapor(nanmax(stn.ndbc_air_t.data))
nanmedian(vapor(stn.ndbc_air_t.data))
nanmean(vapor(stn.ndbc_air_t.data))
var(stn.ndbc_sea_spechumid.data)
0.00001
std(stn.ndbc_sea_spechumid.data)
std(stn.ndbc_spechumid.data)
mean(stn.ndbc_spechumid.data)
% Simple estimates based on natural variance
sqs = 0.0037;
sqs2 = sqs.^2;
sqa = 0.0038;
sqa2 = sqa.^2;
q = qs - qa;					q2 = q.^2;
sq2 = sqs2 + sqa2 - (2.*qs.*qa.*COVqsqa);
sq = sqrt(abs(sq2));
[ig,r_q_w] = cov_ts(q,w);
%DEBUG:
r_q_w = 0;
sQlh2 = airdens2.*Le2.*Cdr.*Cqs.*( ...
(U2.*((q2.*sU2) + sq2)) + (V2.*((q2.*sV2) + sq2)) + (Ug2.*((q2.*sUg2) + sq2)) ...
+ ((2.*sw.*sq.*r_q_w)./(w.*q)) ...
);
%DEBUG:
disp({'sQsh2',nanmin(sQsh2),nanmean(sQsh2),nanmax(sQsh2),});
if ( doPlot )
figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_sensible_heat_flux.date,stn.ndbc_ncep_30a_sensible_heat_flux.data,'k-'); plot(dts,sQsh2,'r.'); datetick3;
end;
errbdgt
nanmean(q)
nanmean(q2)
nanmean(sq2)
nanmean(w.*q)
nanmean(sQlh2)
nanmean(sQlh)
nanmean(w.*th)
nanmean(w.*q)
1/nanmean(w.*q)
1/nanmax(w.*q)
1/nanmin(w.*q)
nanmean(1/(w.*q))
nanmean(1/(w.*th))
nanmean(1./(w.*q))
nanmean(1/(w.*q))
nanmean((1/(w.*q)))
nanmean(1./(w.*q))
nanmean(1./(w.*th))
nanmax(1./(w.*q))
dewp_to_relhumid(0.05,0.05)
dewp_to_relhumid(30.05,30.05)-dewp_to_relhumid(30.00,30.00)
dewp_to_relhumid(30.05,30.00)-dewp_to_relhumid(30.00,30.05)
dewp_to_relhumid(30.00,30.00)-dewp_to_relhumid(30.00,30.05)
dewp_to_relhumid(30.00,25.00)-dewp_to_relhumid(30.05,25.05)
nanmean(a-qa)
d
d=[]; clear d;
d = stn.ndbc_dew_t.data(dix);			d2 = d.^2;
nanmean(a-d)
nanmax(a-d)
nanmin(a-d)
dewp_to_relhumid(30.00,25.00)-dewp_to_relhumid(30.05,25.05)
dewp_to_relhumid(30.05,25.00)-dewp_to_relhumid(30.00,25.05)
relhumid_to_spechumid(30.05,90.5)-relhumid_to_spechumid(30.05,90)-
relhumid_to_spechumid(30.05,90.5)-relhumid_to_spechumid(30.05,90)
0.0001
relhumid_to_spechumid(30.05,50.5)-relhumid_to_spechumid(30.05,50)
help relhumid_to_spechumid
relhumid_to_spechumid(10.05,90.5)-relhumid_to_spechumid(10.05,90)
relhumid_to_spechumid(10.05,90.5)-relhumid_to_spechumid(10.00,90)
relhumid_to_spechumid(10.00,90.5)-relhumid_to_spechumid(10.05,90)
errbdgt
figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_latent_heat_flux.date,stn.ndbc_ncep_30a_latent_heat_flux.data,'k-'); plot(dts,sQlh,'r.'); datetick3; titlename('\sigmaQ_L_H');
figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_sensible_heat_flux.date,stn.ndbc_ncep_30a_sensible_heat_flux.data,'k-'); plot(dts,sQsh,'r.'); datetick3; titlename('\sigmaQ_S_H');
nanmean(sQsh)
nanmean(sQlh)
u,v
[u,v]=spddir_to_uv(14,130)
[u1,v1]=spddir_to_uv(14,130);[u2,v2]=spddir_to_uv(14.5,139);
[u1,v1]=spddir_to_uv(14,130);[u2,v2]=spddir_to_uv(14.5,139);u2-u1,v2-v2
[u1,v1]=spddir_to_uv(14.5,130);[u2,v2]=spddir_to_uv(14,139);u2-u1,v2-v2
mps2kts(0.55)
[u1,v1]=spddir_to_uv(14.9,130);[u2,v2]=spddir_to_uv(14,139);u2-u1,v2-v2
[u1,v1]=spddir_to_uv(14.9,130);[u2,v2]=spddir_to_uv(14,130);u2-u1,v2-v2
[u1,v1]=spddir_to_uv(15.1,130);[u2,v2]=spddir_to_uv(14,139.26);u2-u1,v2-v2
[u1,v1]=spddir_to_uv(15.1,230);[u2,v2]=spddir_to_uv(14,239.26);u2-u1,v2-v2
help spddir_to_uv
[u1,v1]=spddir_to_uv(15.1,230),[u2,v2]=spddir_to_uv(14,239.26),u2-u1,v2-v2
[u1,v1]=spddir_to_uv(15.1,230),[u2,v2]=spddir_to_uv(14,239.26),u2-u1,v2-v1
[u1,v1]=spddir_to_uv(15.1,230);[u2,v2]=spddir_to_uv(14,239.26);u2-u1,v2-v1
[u1,v1]=spddir_to_uv(15.1,130);[u2,v2]=spddir_to_uv(14,39.26);u2-u1,v2-v1
[u1,v1]=spddir_to_uv(15.1,130);[u2,v2]=spddir_to_uv(14,139.26);u2-u1,v2-v1
[u1,v1]=spddir_to_uv(14,130);[u2,v2]=spddir_to_uv(15.1,139.26);u2-u1,v2-v1
[u1,v1]=spddir_to_uv(14,130);[u2,v2]=spddir_to_uv(15.1,134.5);u2-u1,v2-v1
[u1,v1]=spddir_to_uv(15.1,130);[u2,v2]=spddir_to_uv(14,134.5);u2-u1,v2-v1
dewp_to_relhumid(30.09,25.00)-dewp_to_relhumid(30.00,25.31)
dewp_to_relhumid(30.09,25.31)-dewp_to_relhumid(30.00,25)
dewp_to_relhumid(30,25.00)-dewp_to_relhumid(30.09,25.31)
dewp_to_relhumid(30.09,25.00)-dewp_to_relhumid(30.00,25.31)
dewp_to_relhumid(30,25.31)-dewp_to_relhumid(30.09,25)
err
help grpstats
[ix1,ix2] = intersect_dates(stn.ndbc_ncep_30a_latent_heat_flux.date,dts);
[yr,mo,dy,hr,mn,sc] = datevec(stn.ndbc_ncep_30a_latent_heat_flux.date(ix1));
dy = datenum(yr,mo,dy,0,0,0);
fld.date = unique(dy);
fld.data = grpstats(stn.ndbc_ncep_30a_latent_heat_flux.data(ix1),dy);
err.date = unique(dy);
err.data = grpstats(sQlh(ix2),dy,@sum);
figure; maxigraph; hold on;
lh = [];
lh(end+1)=plot(fld.date,fld.data,'r.-');
plot(err.date,fld.data-err.data,'r+');
plot(err.date,fld.data+err.data,'r+');
lh(end+1)=plot(stn.daily_oaflux_latent_heat_flux.date,stn.daily_oaflux_latent_heat_flux.data,'k-');
[ix1,ix2]=intersect_dates(stn.daily_oaflux_latent_heat_flux.date,stn.daily_oaflux_latent_flux_err.date);
plot(stn.daily_oaflux_latent_heat_flux.date(ix1),[stn.daily_oaflux_latent_heat_flux.data(ix1)-stn.daily_oaflux_latent_flux_err.data(ix2),stn.daily_oaflux_latent_heat_flux.data(ix1)+stn.daily_oaflux_latent_flux_err.data(ix2)],'k+');
legend(lh,'Q_L_H \pmerr','OAFlux Q_0 \pmerr', 'Location','Best');
titlename('Daily heat fluxes [W/m^2]');
[ix1,ix2] = intersect_dates(stn.ndbc_ncep_30a_latent_heat_flux.date,dts);
[yr,mo,dy,hr,mn,sc] = datevec(stn.ndbc_ncep_30a_latent_heat_flux.date(ix1));
dy = datenum(yr,mo,dy,0,0,0);
fld.date = unique(dy);
fld.data = grpstats(stn.ndbc_ncep_30a_latent_heat_flux.data(ix1),dy);
err.date = unique(dy);
err.data = grpstats(sQlh(ix2),dy);
figure; maxigraph; hold on;
lh = [];
lh(end+1)=plot(fld.date,fld.data,'r.-');
plot(err.date,fld.data-err.data,'r+');
plot(err.date,fld.data+err.data,'r+');
lh(end+1)=plot(stn.daily_oaflux_latent_heat_flux.date,stn.daily_oaflux_latent_heat_flux.data,'k-');
[ix1,ix2]=intersect_dates(stn.daily_oaflux_latent_heat_flux.date,stn.daily_oaflux_latent_flux_err.date);
plot(stn.daily_oaflux_latent_heat_flux.date(ix1),[stn.daily_oaflux_latent_heat_flux.data(ix1)-stn.daily_oaflux_latent_flux_err.data(ix2),stn.daily_oaflux_latent_heat_flux.data(ix1)+stn.daily_oaflux_latent_flux_err.data(ix2)],'k+');
datetick3;
legend(lh,'Q_L_H \pmerr','OAFlux Q_0 \pmerr', 'Location','Best');
titlename('Daily heat fluxes [W/m^2]');
[ix1,ix2] = intersect_dates(stn.ndbc_ncep_30a_sensible_heat_flux.date,dts);
[yr,mo,dy,hr,mn,sc] = datevec(stn.ndbc_ncep_30a_sensible_heat_flux.date(ix1));
dy = datenum(yr,mo,dy,0,0,0);
fld.date = unique(dy);
fld.data = grpstats(stn.ndbc_ncep_30a_sensible_heat_flux.data(ix1),dy);
err.date = unique(dy);
err.data = grpstats(sQsh(ix2),dy);
figure; maxigraph; hold on;
lh = [];
lh(end+1)=plot(fld.date,fld.data,'r.-');
plot(err.date,fld.data-err.data,'r+');
plot(err.date,fld.data+err.data,'r+');
lh(end+1)=plot(stn.daily_oaflux_sensible_heat_flux.date,stn.daily_oaflux_sensible_heat_flux.data,'k-');
[ix1,ix2]=intersect_dates(stn.daily_oaflux_sensible_heat_flux.date,stn.daily_oaflux_sensible_flux_err.date);
plot(stn.daily_oaflux_sensible_heat_flux.date(ix1),[stn.daily_oaflux_sensible_heat_flux.data(ix1)-stn.daily_oaflux_sensible_flux_err.data(ix2),stn.daily_oaflux_sensible_heat_flux.data(ix1)+stn.daily_oaflux_sensible_flux_err.data(ix2)],'k+');
datetick3;
legend(lh,'Q_S_H \pmerr','OAFlux Q_0 \pmerr', 'Location','Best');
titlename('Daily Sensible Heat Flux [W/m^2]');
errbdgt
dewp_to_relhumid(30,25.31)-dewp_to_relhumid(30.09,25)
dewp_to_relhumid(30,25)-dewp_to_relhumid(30.09,25.31)
dewp_to_relhumid(30,25.31)-dewp_to_relhumid(30.09,25)
dewp_to_relhumid(30.09,25.31)-dewp_to_relhumid(30,25)
[u1,v1]=spddir_to_uv(14,130);[u2,v2]=spddir_to_uv(15.1,139.26);u2-u1,v2-v1
relhumid_to_spechumid(10.00,91.8)-relhumid_to_spechumid(10.09,90)
relhumid_to_spechumid(10,91.8)-relhumid_to_spechumid(10.09,90)
relhumid_to_spechumid(10,90)-relhumid_to_spechumid(10.09,91.8)
relhumid_to_spechumid(10.09,90)-relhumid_to_spechumid(10,91.8)
relhumid_to_spechumid(10.08,100)-relhumid_to_spechumid(10,100)
relhumid_to_spechumid(30.08,100)-relhumid_to_spechumid(30,100)
nanmean(stn.ndbc_sea_t.data)
relhumid_to_spechumid(26.08,100)-relhumid_to_spechumid(26,100)
errbdgt
print('-dtiff','../figs/looe1-latent-error-bars.tiff');
print('-dtiff','../figs/looe1-sensible-error-bars.tiff');
stn
fclose('all'); close all; clear all; dbclear all; pack;
%-- 12/30/10  1:27 PM --%
fclose('all'); close all; clear all; dbclear all; pack;
stn = tryfkeys('mlrf1',[],1.20,1.20,0);
errbdgt
stn=rmfield(stn,'ndbc_dew_t')
errbdgt
length(find(~isreal(th)))
length(find(~isreal(w)))
length(find(Ug2 < 0))
length(find(Ug2 == 0))
length(find(Ug < 0))
length(find(U2 == 0))
length(find((U.*U) == 0))
figure; plot(U2)
length(find(Ug2 < eps))
length(find(U2 < eps))
length(find(U2 < -eps))
length(find(U2 < -eps/2))
length(find(U2 < -eps./2))
length(find(U2 < 0))
length(find(U2 <= 0))
length(find(U2 < 0))
length(find(V2 < 0))
length(find(Ug2 < 0))
length(find(~isreal(w)))
isreal(w)
length(find(imag(w)>0))
length(find(imag(w)~=0))
ix=find(imag(w)~=0);
U2(ix)
[U2(ix),V2(ix),Ug2(ix)]
errbdgt
dbstop errbdgt
errbdgt
figure; plot(dts,[w,sw]);
figure; plot(dts,[w,sw]); datetick3;
figure; plot(dts,[w,sw],'.); datetick3;
figure; plot(dts,[w,sw],'.'); datetick3;
length(find(~isfinite((U2.*((th2.*sU2) + sth2)) + (V2.*((th2.*sV2) + sth2)) + (Ug2.*((th2.*sUg2) + sth2)))))
length(find(~isfinite(((2.*sw.*sth.*r_th_w)./(w.*th)))))
length(find(~isfinite(sw)))
length(find(~isfinite(sth)))
length(find(~isfinite( r_th_w )))
length(find(~isfinite( w )))
length(find(~isfinite( th )))
length(find(th==0))
length(find(w==0))
dbquit
dbclear all
errbdgt
lookfor auto
help corrmtx
help soundsc
soundsc(stn.ndbc_sea_t.data)
soundsc(stn.ndbc_wind1_speed.data)
soundsc(stn.ndbc_ncep_30a_absorbed_heat_flux.data);
soundsc(stn.ndbc_ncep_30a_absorbed_heat_flux_24_hour_average.data);
soundsc(real(stn.ndbc_ncep_30a_absorbed_heat_flux.data));
soundsc(stn.tmd_tide_i_depth.data);
help sound
load handel
sound(y,Fs)
clear y Fs
sound(y,Fs,16)
sound(y,Fs)
load handel
sound(y,Fs,16)
sound(y,Fs,8)
clear y Fs
help sound
load handel
sound(y,4096)
sound(y,1024)
sound(y,(8192*4))
sound(y,(8192*2))
clear y Fs
soundsc(stn.tmd_tide_i_depth.data,4096);
soundsc(stn.tmd_tide_i_depth.data,1024);
soundsc([U,V]);
soundsc([U,V],4096);
soundsc([s,a],4096);
soundsc([s,a],8192);
soundsc([s,q0],8192);
x=load_station_data('mlrf1');
scatter_fit_ts(x.wind1_speed,x.wind2_speed);
x.wind1_u
x = verify_variable(x,'wind1_u');
x = verify_variable(x,'wind1_v');
x = verify_variable(x,'wind2_u');
x = verify_variable(x,'wind2_v');
scatter_fit_ts(x.wind1_u,x.wind2_u);
scatter_fit_ts(x.wind1_u,x.wind2_u,@(x)(find(abs(x.data)>eps)));
scatter_fit_ts(x.wind1_u,x.wind2_u,@(x)(find(abs(x.data)>eps)),@(x)(find(abs(x.data)>eps)));
scatter_fit_ts(x.wind1_u,x.wind2_u,@(x)(find(abs(x.data)>0.1)),@(x)(find(abs(x.data)>0.1)));
scatter_fit_ts(x.wind1_v,x.wind2_v,@(x)(find(abs(x.data)>0.1)),@(x)(find(abs(x.data)>0.1)));
scatter_fit_ts(x.wind1_u,x.wind2_u,@(x)(find(abs(x.data)>=1)),@(x)(find(abs(x.data)>=1)));
scatter_fit_ts(x.wind1_u,x.wind1_v,@(x)(find(abs(x.data)>=1)),@(x)(find(abs(x.data)>=1)));
daspect([1 1 1])
x=[]; clear x
[a,b,c,d,e,f]=station_instrument_heights('mlrf1')
help station_instrument_heights
x=load_station_data('mlrf1');
x = verify_variable(x,'wind1_u');
x = verify_variable(x,'wind1_v');
x = verify_variable(x,'wind2_u');
x = verify_variable(x,'wind2_v');
scatter_fit_ts(x.wind1_u,x.wind2_u);
scatter_fit_ts(x.wind1_v,x.wind2_v);
3+0.1*0.1
x=load_station_data('fwyf1');
x = verify_variable(x,'wind1_u');
x = verify_variable(x,'wind1_v');
x = verify_variable(x,'wind2_u');
x = verify_variable(x,'wind2_v');
scatter_fit_ts(x.wind1_u,x.wind2_u);
scatter_fit_ts(x.wind1_v,x.wind2_v);
scatter_fit_ts(x.wind1_u,x.wind2_u);
x=[]; clear x
x=load_station_data('smkf1');
x = verify_variable(x,'wind1_u');
x = verify_variable(x,'wind1_v');
x = verify_variable(x,'wind2_u');
x = verify_variable(x,'wind2_v');
scatter_fit_ts(x.wind1_u,x.wind2_u);
scatter_fit_ts(x.wind1_v,x.wind2_v);
scatter_fit_ts(x.wind1_speed,x.wind2_speed);
x=[]; clear x
corerrbdgt
help matlabhome
help mhome
help path
help genpath
matlabroot
corerrbdgt
dwq
sQlh(1)
hl
[ix1,ix2]=intersect_dates(dts,stn.ndbc_ncep_30a_latent_heat_flux.date);
stn.ndbc_ncep_30a_latent_heat_flux.data(ix2(1))
stn.ndbc_ncep_30a_latent_heat_flux.data(qlhix(1))
corerrbdgt
dwq
sQlh(1)
stn.ndbc_ncep_30a_latent_heat_flux.data(qlhix(1))
hl
corerrbdgt
size(aix)
size(a)
clear qlhix
corerrbdgt
dwq
sQlh(ix)
hl
qlh(ix)
size(a)
nanmean(stn.ndbc_barom)
nanmean(stn.ndbc_barom.data)
nanmedian(stn.ndbc_barom.data)
nanstd(stn.ndbc_barom.data)
nanvar(stn.ndbc_barom.data)
stn.ncep_dlrf
clear qlhix
corerrbdgt
dwq
s(ix),a(ix)
dwth
dwq
hs
hl
qlh(ix)
qsh(ix)
clear qlhix
corerrbdgt
dwq
sQlh(ix)
hl
qlh(ix)
hs
qsh(ix)
dwth
sQsh(ix)
for ix=aix(:)'; foo(ix)=0; end;
size(foo)
size(aix)
for ix=aix(:); foo(ix)=0; end;
clear foo
for ix=aix(:); foo(ix)=0; end;
size(foo)
size(aix)
clear foo
for ix=aix(:); foo(ix)=0; end;
size(foo)
size(aix)*2
clear foo
numel(a)
clear qlhix
corerrbdgt
clear qlhix
corerrbdgt
dwth
dwth(1:10)
sQsh(1:10)
sQsh(1:10)'
hs(1:10)
qsh(1:10)'
dwq(1:10)
sQlh(1:10)'
hl(1:10)
qlh(1:10)'
qlh(1:10)'-hl(1:10)
qsh(1:10)'-hs(1:10)
clear qlhix
corerrbdgt
dwq(1:10)
sQlh(1:10)'
stn
dwq(1:10),sQlh(1:10)'
dwth(1:10),sQsh(1:10)'
errbdgt
dwth(1:10),sQsh(1:10)'
dwq(1:10),sQlh(1:10)'
errbdgt
dwq(1:10),sQlh(1:10)'
errbdgt
dwq(1:10),sQlh(1:10)'
mean(dwq)
nanmean(dwq)
nanmean(dwth)
mean(dwq(isfinite(dwq)))
mean(dwq(isfinite(dwth)))
mean(dwth(isfinite(dwth)))
which median
which iqr
mean(sQsh(isfinite(sQsh)))
mean(sQlh(isfinite(sQlh)))
84/12
.85/7.2
7.2/.85
Cdn
Cth
y(16:18)
Cdr
Cqs
stn.ncep_dsrf
scatter_fit_ts(stn.bic_surf_par,stn.ncep_dsrf);
multiplot_station(stn,{'bic_surf_par','ncep_dsrf'});
is
is=grpstats(stn.bic_surf_par.data,floor(stn.bic_surf_par.date));
ra=grpstats(stn.ncep_dsrf.data,floor(stn.ncep_dsrf.date));
is.data=is; ra.data=ra;
is
is.date=floor(stn.bic_surf_par.date);
ra.date=floor(stn.ncep_dsrf.date);
scatter_fit_ts(is,ra);
is.date=unique(floor(stn.bic_surf_par.date));
ra.date=unique(floor(stn.ncep_dsrf.date));
scatter_fit_ts(is,ra);
stn.sat_net_heat_flux
scatter_fit_ts(stn.ncep_net_heat_flux,stn.sat_net_heat_flux)
scatter_fit_ts(stn.ncep_net_heat_flux,stn.ndbc_ncep_30a_net_heat_flux)
scatter_fit_ts(stn.ncep_net_heat_flux,stn.sat_net_heat_flux)
x =
help get_satpar_insol
x
x=[]; clear x
x=get_satpar_insol('mlrf1');
scatter_fit_ts(stn.ncep_dsrf,x.sat_insol_in);
scatter_fit_ts(stn.bic_surf_par,x.sat_insol_in);
scatter_fit_ts(stn.bic_surf_par,stn.ncep_dsrf);
scatter_fit_ts(stn.bic_surf_par,stn.ncep_par);
help grpscatter
help scattergrp
help scatterhist
[bicix,parix]=intersect_dates(stn.bic_surf_par.date,stn.ncep_par.date);
scatterhist(stn.bic_surf_par.data(bicix),stn.ncep_par.data(parix));
scatter_fit_ts(is,ra);
is.data=grpstats(stn.bic_surf_par.data,floor(stn.bic_surf_par.date),@nanmax);
ra.data=grpstats(stn.ncep_dsrf.data,floor(stn.ncep_dsrf.date),@nanmax);
scatter_fit_ts(is,ra);
clear ra
ra.date=unique(floor(stn.ncep_dsrf.date));
ra.data=grpstats(stn.ncep_dsrf.data,floor(stn.ncep_dsrf.date),@max);
ra
length(find(isnan(stn.ncep_dsrf.data)))
rawra=stn.ncep_dsrf.data;
rawra.data=stn.ncep_dsrf.data;
rawra.date=stn.ncep_dsrf.date;
rawra.date(isnan(rawra.data))=[];
rawra.data(isnan(rawra.data))=[];
ra.date=unique(floor(rawra.date));
ra.data=grpstats(rawra.data,floor(rawra.date),@max);
scatter_fit_ts(is,ra);
clear is ra x
grepstruct(stn,'oafl')
scatter_fit_ts(stn.ncep_srf,stn.daily_oaflux_srf)
nanmean(stn.ncep_srf.data)
nanmean(stn.daily_oaflux_srf.data)
rawra.date=stn.ncep_dsrf.date;
rawra.data=stn.ncep_dsrf.data;
rawra.date(isnan(rawra.data))=[];
rawra.data(isnan(rawra.data))=[];
ra.date=unique(floor(rawra.date));
ra.data=grpstats(rawra.data,floor(rawra.date),@mean);
nanmean(ra.data)
scatter_fit_ts(ra,stn.daily_oaflux_srf)
scatter_fit_ts(stn.bic_surf_par,stn.daily_oaflux_srf)
scatter_fit_ts(ra,stn.daily_oaflux_srf)
