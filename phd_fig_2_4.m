1;

disspath = get_thesis_path('../DISS');

doNorm = true;

xlm=[1e-1,3e3];

ylm=[1e-3,1.5e6];
%ylm_ts=[1e-5,1.5e6];
ylm_td=[1e-5,1.5e4];
ylm_cu=[1e-5,1.5e2];
ylm_qa=[1e-9,1.5e2];

lineWidth=2.5;

if ( ~exist('fwyf1','var') )
  fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
  lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1);
  sanf1 = get_station_from_station_name('sanf1'); sanf1 = load_all_ndbc_data(sanf1);
  mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1);
  smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1);
  looe1 = get_station_from_station_name('looe1'); looe1 = load_all_ndbc_data(looe1);
  looe1 = get_looe1_adcp(looe1);

  smkf1 = station_dewp_to_relhumid(smkf1,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
  smkf1 = station_relhumid_to_spechumid(smkf1,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
end;

pnl='a'-1;

% Row 1
pnl=pnl+1;
[P,W,fh,lh]=plot_spec(fwyf1,'ndbc_wind1_speed',[],[],xlm,ylm,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));

pnl=pnl+1;
[P,W,fh,lh]=plot_spec(lonf1,'ndbc_wind1_speed',[],[],xlm,ylm,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));

pnl=pnl+1;
[P,W,fh,lh]=plot_spec(sanf1,'ndbc_wind1_speed',[],[],xlm,ylm,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));

% Row 2
pnl=pnl+1;
[P,W,fh,lh]=plot_spec(mlrf1,'ndbc_air_t',[],[],xlm,ylm,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));

pnl=pnl+1;
[P,W,fh,lh]=plot_spec(smkf1,'ndbc_air_t',[],[],xlm,ylm,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));

pnl=pnl+1;
[P,W,fh,lh]=plot_spec(smkf1,'ndbc_spechumid',[],[],xlm,ylm_qa,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));


% Row 3
pnl=pnl+1;
[P,W,fh,lh]=plot_spec(smkf1,'ndbc_barom',[],[],xlm,ylm,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));

pnl=pnl+1;
[P,W,fh,lh]=plot_spec(smkf1,'ndbc_tide',[],[],xlm,ylm_td,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));

pnl=pnl+1;
[P,W,fh,lh]=plot_spec(lonf1,'ndbc_tide',[],[],xlm,ylm_td,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));


% Row 4
pnl=pnl+1;
[P,W,fh,lh]=plot_spec(looe1,'adcp_l',[],[],xlm,ylm_cu,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));

pnl=pnl+1;
[P,W,fh,lh]=plot_spec(looe1,'adcp_x',[],[],xlm,ylm_cu,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));

pnl=pnl+1;
[P,W,fh,lh]=plot_spec(looe1,'adcp_sfc_speed',[],[],xlm,ylm_cu,[],true);
y=ylim; text(2e-1,y(2)/10,['(',char(pnl),')'],'FontSize',30);
set(lh,'LineWidth',lineWidth,'Color','k');
print('-dtiff','-r300',fullfile(disspath,[mfilename,'_',char(pnl),'.tif']));



clear P W ans disspath doNorm fh
clear lh lineWidth pnl xlm y yl ylm ylm_*
