1;

ufld = 'adcp_baroclinic_x_3hlp';
bfld = 'b_ndbc_erai_erai_30a_avhrr_dt_1_d_sum';
wfld = 'ndbc_erai_erai_30a_wind_stress_xshore';

if ( exist('stn','var') )
  stn = verify_variable(stn,bfld);
end;
if ( ~exist('stn','var') || ~isfield(stn,bfld) )
  error('Please place a heat-budget result in struct STN first! STN.%s',bfld);
end;
stn = verify_variable(stn,wfld);

if ( ~exist('looe1','var') || ~isfield(looe1,'adcp_seatemp') )
  looe1 = []; clear looe1
  looe1 = get_station_from_station_name('looe1');
  looe1 = get_looe1_microcat(looe1);
  looe1 = get_looe1_adcp(looe1);
end;

if ( ~isfield(looe1,'dTdz') )
  looe1.dTdz = ts_op(looe1.microcat_seatemp,looe1.adcp_seatemp,'-');
end;
if ( ~isfield(looe1,ufld) )
  disp(['Trying to generate LOOE1.',ufld]);
  looe1 = verify_variable(looe1,ufld);
end;

% Series of cold fronts throughout which all necessary time series have data
%%%dts = datenum(2006,1,1):(1/24):datenum(2006,1,31);
%%dts = datenum(2008,10,23):(1/24):datenum(2008,11,23);
%dts = datenum(2008,10,15):(1/24):datenum(2008,11,30);
dts = datenum(2008,10,15):(1/24):datenum(2008,11,16);

[tix,uix,bix,wix,ig] = intersect_all_dates([],looe1.dTdz.date,looe1.(ufld).date,stn.(bfld).date,stn.(wfld).date,dts);

fmg;

subplot(5,1,1:4);
hold on;
contourf(looe1.(ufld).date(uix),looe1.adcp_bin_heights,looe1.(ufld).prof(uix,:)'),
xlim([dts(1),dts(end)+4.5]);
datetick('x',2,'keeplimits');
set(gca,'xticklabel',[]);
ylim([0,20.5]);
ylabel('Height above bottom [m]');
set(gca,'clim',[-0.45,+0.45]);
lh=colorbar('East');
set(lh,'FontSize',7);
text(dts(end)+0.8, 7.0,'Onshore','Rotation',270,'FontSize',9);
text(dts(end)+0.8,17.0,'Offshore','Rotation',270,'FontSize',9);
text(dts(end)+0.5,19.7,'m\bullets^-^1','FontSize',9);

text(datenum(2008,10,20,12,0,0),-0.8,'Cooling','FontSize',9);
text(datenum(2008,10,24,12,0,0),-0.8,sprintf('Daily\nwarming'),'FontSize',9);
text(datenum(2008,10,29, 0,0,0),-0.8,'Strong cold front','FontSize',9);
%text(datenum(2008,11,18, 0,0,0),-0.8,'Cold front #2','FontSize',9);
grid on;

subplot(5,1,5);
hold on;
plot(looe1.dTdz.date(tix),looe1.dTdz.data(tix),'b-','Color',[.3,.3,.3]);
plot(stn.(bfld).date(bix),stn.(bfld).data(bix),'m--','Color',[.7,.7,.7],'LineW',2.0);
plot(stn.(wfld).date(wix),10.*stn.(wfld).data(wix),'k:','Color',[0,0,0],'LineW',1.5)
xlim([dts(1),dts(end)+4.5]);
datetick('x',2,'keeplimits');
ylim([-3.0,+3.0]);
grid on;
lh=legend('\fontsize{7}LOOE1 \fontsize{10}\DeltaT_s/\Deltaz',...
          ['\fontsize{7}',upper(stn.station_name),' \fontsize{10}\Sigma_1_d\partial_tT_s'],...
          ['\fontsize{7}',upper(stn.station_name),' \fontsize{10}\tau\bullet\nablah [dyn/cm^2]'],...
          'Location','East');

print('-dtiff',fullfile(get_thesis_path('../figs'),[mfilename,'-',datestr(dts(1),'yyyy-mm-dd'),'-',datestr(dts(end),'yyyy-mm-dd'),'.tif']));

cm=colormap(gray);
print('-dtiff',fullfile(get_thesis_path('../figs'),[mfilename,'-',datestr(dts(1),'yyyy-mm-dd'),'-',datestr(dts(end),'yyyy-mm-dd'),'-GRAY.tif']));

%clear bfld bix ig lh t td tix ts ud ufld uix us wfld wix x cm dts ans
