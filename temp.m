ww3 = seasonalize_ww3_region('uks');

dts=datenum(2005,2,1):(3/24):datenum(2017,6,1)-(1/24);

fld = get_ww3_multi_region([],'fknms',[],[]);


rh(1)=rectangle('Position',bbox2rect([ -80.60, -79.90, 24.70, 25.30]),'EdgeColor','r')
rh(2)=rectangle('Position',bbox2rect([ -81.30, -80.50, 24.40, 25.10]),'EdgeColor','r')
rh(3)=rectangle('Position',bbox2rect([ -82.00, -81.20, 24.30, 25.00]),'EdgeColor','r')


rh=rectangle('Position',bbox2rect([ -80.60, -80.00, 24.60, 25.30]),'EdgeColor','r')
rectangle(bbox2rect([ -80.60, -80.00, 24.60, 25.30])


          %  http://polar.ncep.noaa.gov/pub/history/waves/multi_1/201703/gribs/



for yrmo in 
for yr in range(begyr,endyr+1):

  for mo in range(begmo,endmo+1):







  begyr = 1987;
  endyr = 2016;
  begmo = 1
  endmo = 12
else
  begyr = 2013
  endyr = 2013
  begmo = 3
  endmo = 12



dayspermo = (31,28,31,30,31,30,31,31,30,31,30,31)








  begyr = 1987;
  endyr = 2016;
  begmo = 1
  endmo = 12

else:
  outpath = datapath;
  gridszs = "1.5/1.5";
  begyr = 2013
  endyr = 2013
  begmo = 3
  endmo = 12

dayspermo = (31,28,31,30,31,30,31,31,30,31,30,31)


client = ECMWFDataServer()

for yr in range(begyr,endyr+1):

  for mo in range(begmo,endmo+1):

    yrmo = "%04d%02d" %(yr,mo);

    enddy = dayspermo[mo-1];
    if ( mo == 2 and calendar.isleap(yr) ):
      enddy = 29;

    # E.g., "date": "2013-02-01/to/2013-02-28",
    datestr = "%04d-%02d-01/to/%04d-%02d-%02d" %(yr,mo,yr,mo,enddy);





        "target": "/cygdrive/e/heat_budget/data/ECMWF/ERA_Interim_%s_fc.grib"%yrmo,






  % % NOT helpful
  % scatter_fit_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_1_d_avg,stn.ndbc_sea_t_diff_flux,[],[],'\partial_tT','\Delta\mu_1_dT',[],[],true,0);








  stn.erai_air_t_interp.date = fld.dts;
  stn.erai_air_t_interp.data = interp_field(fld.lats,fld.lons,fld.Ta,stn.lat,stn.lon,'linear','warn');




FROM plot_daily_clim_heat_budget.m:
  if ( ~exist('begyr','var') || ~exist('endyr','var') )
    [t,sqt,bqt,dt,q,qe] = ...
        intersect_tses(stn.(sfld),stn.(sqtfld),stn.(bq0tfld),stn.(bdTfld),stn.(hcdTdt),stn.([hcdTdt,'_err']));
    if ( ~exist('begyr','var') );	begyr = get_year(t.date(1));		end;
    if ( ~exist('endyr','var') );	endyr = get_year(t.date(end)+1);	end;
  end;





%stnm = smkf1.station_name;
%td = smkf1.ndbc_tide_m;
%wd = smkf1.ndbc_wind1_speed;
%br = smkf1.ndbc_barom;
stnm = lonf1.station_name;
td = lonf1.ndbc_tide_m;
wd = lonf1.ndbc_wind1_speed;
br = lonf1.ndbc_barom;





1;
%%%% SCRIPT to create spectral analysis figures for JMR MS. from Ch. 2

if ~exist('fwyf1','var'); fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1); end;
if ~exist('lonf1','var'); lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1); end;
if ~exist('mlrf1','var'); mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1); end;
if ~exist('smkf1','var'); smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1); end;
if ~exist('sanf1','var'); sanf1 = get_station_from_station_name('sanf1'); sanf1 = load_all_ndbc_data(sanf1); end;
if ~exist('dryf1','var'); dryf1 = get_station_from_station_name('dryf1'); dryf1 = load_all_ndbc_data(dryf1); end;
if ~exist('looe1','var'); looe1 = get_station_from_station_name('looe1'); looe1 = get_looe1_adcp(looe1); end;
if ~isfield(smkf1,'ndbc_spechumid')
  smkf1 = station_dewp_to_relhumid(smkf1,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
  smkf1 = station_relhumid_to_spechumid(smkf1,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
end;
if ~isfield(smkf1,'ndbc_wind1_speed_10m')
  fwyf1 = station_wind_at_height(fwyf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  lonf1 = station_wind_at_height(lonf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  mlrf1 = station_wind_at_height(mlrf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  smkf1 = station_wind_at_height(smkf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  sanf1 = station_wind_at_height(sanf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
  dryf1 = station_wind_at_height(dryf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_air_t');
end;
if ~isfield(smkf1,'ndbc_tide_m')
  lonf1.ndbc_tide_m = lonf1.ndbc_tide; lonf1.ndbc_tide_m.data = unitsratio('m','ft')*(lonf1.ndbc_tide.data-nanmean(lonf1.ndbc_tide.data));
  smkf1.ndbc_tide_m = smkf1.ndbc_tide; smkf1.ndbc_tide_m.data = unitsratio('m','ft')*(smkf1.ndbc_tide.data-nanmean(smkf1.ndbc_tide.data));
end;
if ~isfield(looe1,'adcp_x_contig')
  %looe1.adcp_x_contig = subset_ts(looe1.adcp_x,@(x)(find(datenum(2005,3,25)<=x.date&x.date<=datenum(2008,6,11))));
  looe1.adcp_x_contig = subset_ts(looe1.adcp_x,@(x)(find(datenum(2005,3,25)<=x.date&x.date<=datenum(2009,1,25))));
  looe1.adcp_x_contig = subset_ts(looe1.adcp_x_contig,@ts_isfinite);
  %looe1.adcp_l_contig = subset_ts(looe1.adcp_l,@(x)(find(datenum(2005,3,25)<=x.date&x.date<=datenum(2008,6,11))));
  looe1.adcp_l_contig = subset_ts(looe1.adcp_l,@(x)(find(datenum(2005,3,25)<=x.date&x.date<=datenum(2009,1,25))));
  looe1.adcp_l_contig = subset_ts(looe1.adcp_l_contig,@ts_isfinite);
end;

%% MAKE FIGURES

%CI = 0.68;
%%CI = 0.80;
%%CI = 0.84;
CI = 0.85;
%%CI = 0.90;
%CI = 0.95;

doPrint = true;
%doPrint = false;

if 1
ch2_spec(fwyf1.ndbc_wind1_speed_10m,20,true,CI,'FWYF1 U10'); axis([3/24,3000,0,220]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'fwyf1-pmtm-U10.tif')); end;
ch2_spec(mlrf1.ndbc_air_t,20,true,CI,'MLRF1 T_a'); axis([3/24,3000,0,150]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'mlrf1-pmtm-Ta.tif')); end;
ch2_spec(smkf1.ndbc_spechumid,20,true,CI,'SMKF1 q_a'); axis([3/24,3000,0,3.2e-5]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-qa.tif')); end;
ch2_spec(smkf1.ndbc_barom,20,true,CI,'SMKF1 P_a'); axis([3/24,3000,0,400]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-Pa.tif')); end;
ch2_spec(smkf1.ndbc_tide_m,20,true,CI,'SMKF1 h_t_i_d_e'); axis([3/24,3000,0,27]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-h.tif')); end;
end;

if 1+1
%ch2_spec(looe1.adcp_x_contig,5,true,CI,'LOOE1 u_x_s'); axis([3/24,400,0,0.06]);
ch2_spec(looe1.adcp_x_contig,5,true,CI,'LOOE1 u_x_s'); axis([3/24,400,0,0.05]);
%ch2_spec(looe1.adcp_x_contig,20,true,CI,'LOOE1 u_x_s'); axis([3/24,400,0,0.05]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'looe1-pmtm-xs.tif')); end;
end;

if 1
%ch2_spec(lonf1.ndbc_sea_t,5,true,CI,'LONF1 T_s'); axis([3/24,3000,0,300]);
ch2_spec(lonf1.ndbc_sea_t,5,true,CI,'LONF1 T_s'); axis([3/24,3000,0,220]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'lonf1-pmtm-Ts.tif')); end;

%%ch2_spec(smkf1.ndbc_sea_t,2,true,CI,'SMKF1 T_s'); axis([3/24,3000,0,300]);
%ch2_spec(smkf1.ndbc_sea_t,5,true,CI,'SMKF1 T_s'); axis([3/24,3000,0,150]);
ch2_spec(smkf1.ndbc_sea_t,5,true,CI,'SMKF1 T_s'); axis([3/24,3000,0,220]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-Ts.tif')); end;

%ch2_spec(smkf1.ndbc_sea_t,2,false,CI,'SMKF1 T_s'); axis([19,31,0,10]);
ch2_spec(smkf1.ndbc_sea_t,20,false,CI,'SMKF1 T_s'); axis([23,31,0,1.8]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-Ts-diurnal.tif')); end;
end;

if 1
%ch2_spec(smkf1.ndbc_sea_t,2,false,CI,'SMKF1 T_s'); axis([19,31,0,10]);
ch2_spec(smkf1.ndbc_sea_t,20,false,CI,'SMKF1 T_s'); axis([26,30,0,0.11]);
set(gcf,'OuterPos',[.2,0,.6,1]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-Ts-inertial-INSET.tif')); end;
end;


%% CH2_SPEC_FIGS.m:
%% EXPLORATIONS

if 0
ch2_spec(fwyf1.ndbc_wind1_speed_10m,5,true,CI,'FWYF1 U10'); axis([3/24,3000,0,250]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'fwyf1-pmtm-U10.tif')); end;
ch2_spec(mlrf1.ndbc_air_t,5,true,CI,'MLRF1 T_a'); axis([3/24,3000,0,200]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'mlrf1-pmtm-Ta.tif')); end;
ch2_spec(smkf1.ndbc_spechumid,5,true,CI,'SMKF1 q_a'); axis([3/24,3000,0,1e-4]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-qa.tif')); end;
ch2_spec(smkf1.ndbc_barom,5,true,CI,'SMKF1 P_a'); axis([3/24,3000,0,1000]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-Pa.tif')); end;
ch2_spec(smkf1.ndbc_tide_m,5,true,CI,'SMKF1 h_t_i_d_e'); axis([3/24,3000,0,50]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'smkf1-pmtm-h.tif')); end;
end;



if 0
ch2_spec(fwyf1.ndbc_sea_t,2,true,CI,'FWYF1 T_s'); axis([3/24,3000,0,300]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'fwyf1-pmtm-Ts.tif')); end;

ch2_spec(fwyf1.ndbc_sea_t,2,false,CI,'FWYF1 T_s'); axis([19,30,0,10]);
if doPrint; pause(0.5); print('-dtiff',fullfile(get_thesis_path('../figs'),'fwyf1-pmtm-Ts-diurnal.tif')); end;
pause;
end;


if 0
% ch2_spec(fwyf1.ndbc_wind1_speed,5,false,CI,'FWYF1 U'); axis([23,31,0,10]);
% annotline(2*pi/sw_f(fwyf1.lat)/3600);
% ch2_spec(lonf1.ndbc_wind1_speed,5,false,CI,'LONF1 U'); axis([23,31,0,10]);
% annotline(2*pi/sw_f(lonf1.lat)/3600);
% ch2_spec(mlrf1.ndbc_wind1_speed,5,false,CI,'MLRF1 U'); axis([23,31,0,10]);
% annotline(2*pi/sw_f(mlrf1.lat)/3600);
% ch2_spec(smkf1.ndbc_wind1_speed,5,false,CI,'SMKF1 U'); axis([23,31,0,10]);
% annotline(2*pi/sw_f(smkf1.lat)/3600);
ch2_spec(sanf1.ndbc_wind1_speed,5,false,CI,'SANF1 U'); axis([23,31,0,10]);
annotline(2*pi/sw_f(sanf1.lat)/3600);
appendtitlename(' INERTIAL WIND VAR.');
ch2_spec(sanf1.ndbc_barom,5,false,CI,'SANF1 P_a'); axis([23,31,0,10]);
annotline(2*pi/sw_f(sanf1.lat)/3600);
% ch2_spec(dryf1.ndbc_wind1_speed,5,false,CI,'DRYF1 U'); axis([23,31,0,10]);
% annotline(2*pi/sw_f(dryf1.lat)/3600);
end;

if 0
ch2_spec(fwyf1.ndbc_wind1_speed,5,true,CI,'FWYF1 U'); axis([3/24,3000,0,500]);
ch2_spec(lonf1.ndbc_wind1_speed,5,true,CI,'LONF1 U'); axis([3/24,3000,0,500]);
ch2_spec(mlrf1.ndbc_wind1_speed,5,true,CI,'MLRF1 U'); axis([3/24,3000,0,500]);
ch2_spec(smkf1.ndbc_wind1_speed,5,true,CI,'SMKF1 U'); axis([3/24,3000,0,500]);
ch2_spec(sanf1.ndbc_wind1_speed,5,true,CI,'SANF1 U'); axis([3/24,3000,0,500]);
end;

if 0
ch2_spec(fwyf1.ndbc_wind1_speed_10m,5,true,CI,'FWYF1 U10'); axis([3/24,3000,0,500]);
ch2_spec(lonf1.ndbc_wind1_speed_10m,5,true,CI,'LONF1 U10'); axis([3/24,3000,0,500]);
ch2_spec(sanf1.ndbc_wind1_speed_10m,5,true,CI,'SANF1 U10'); axis([3/24,3000,0,500]);
ch2_spec(smkf1.ndbc_wind1_speed_10m,5,true,CI,'SMKF1 U10'); axis([3/24,3000,0,500]);
end;

if 0
ch2_spec(smkf1.ndbc_air_t,5,true,CI,'SMKF1 T_a');
ch2_spec(smkf1.ndbc_spechumid,5,true,CI,'SMKF1 q_a');
end;





  fmg;
  % %hist(P,W);
  % %plot(W,P,'k-',W,Pc(:,1),'k^',W,Pc(:,2),'kv');
  % %plot(W,Pc(:,1),'k^',W,Pc(:,2),'kv');
  % %%set(gca,'XScale','log','YScale','log');
  % plot(W,P,'k-',W,Pc(:,1),'k:',W,Pc(:,2),'k:');
  plot(W,P,'k-',W,Pc(:,1),'k.',W,Pc(:,2),'k.');
  if ( perDay )
    %set(gca,'XScale','log','YScale','linear');
    set(gca,'XScale','log','YScale','log');
    xlabel('Days');
  else
    %set(gca,'XScale','linear','YScale','linear');
    set(gca,'XScale','linear','YScale','log');
    xlabel('Hours');
  end;





alldat = [fwyf1.ndbc_sea_t.data;...
          mlrf1.ndbc_sea_t.data;...
          lonf1.ndbc_sea_t.data; ...
          smkf1.ndbc_sea_t.data;...
          sanf1.ndbc_sea_t.data;...
          dryf1.ndbc_sea_t.data];
% allgps = [repmat(1,[1,numel(fwyf1.ndbc_sea_t.data)]),...
%           repmat(2,[1,numel(mlrf1.ndbc_sea_t.data)]),...
%           repmat(3,[1,numel(lonf1.ndbc_sea_t.data)]),...
%           repmat(4,[1,numel(smkf1.ndbc_sea_t.data)]),...
%           repmat(5,[1,numel(sanf1.ndbc_sea_t.data)]),...
%           repmat(6,[1,numel(dryf1.ndbc_sea_t.data)])];
allgps = [repmat('FWY',[numel(fwyf1.ndbc_sea_t.data),1]);...
          repmat('MLR',[numel(mlrf1.ndbc_sea_t.data),1]);...
          repmat('LON',[numel(lonf1.ndbc_sea_t.data),1]);...
          repmat('SMK',[numel(smkf1.ndbc_sea_t.data),1]);...
          repmat('SAN',[numel(sanf1.ndbc_sea_t.data),1]);...
          repmat('DRY',[numel(dryf1.ndbc_sea_t.data),1])];





[ts.date,ts.data] = gap_expand(ts.date,ts.data);





[B,Stats,fh,lh]=scatter_fit(yearfrac(lonf1.ndbc_tide.date(ix)),lonf1.ndbc_tide.data(ix).*12*.0254,'Year','Tide [m]'); titlename(['LONF1 Tide Height ',char(fn)]);
set(lh(1),'Color',[.5,.5,.5],'LineStyle',':');
set(lh(2),'Color','k','LineWidth',2);
%legend hide;
axis([2001,2010,-0.6,1.6]);





1;

fld='ndbc_sea_t';

ys=1987:2013;
for cst={fwyf1,mlrf1,lonf1,smkf1,sanf1};
  t=cst{:}.(fld);
  disp(cst{:}.station_name);
  for yr=ys(:)';
    yrix = find(get_year(t.date)==yr);
    if (numel(yrix)<(305*24));
      disp(yr);
      ys(ys==yr)=[];
    else
      [ig1,ig2,ixes]=find_date_ranges(t.date(yrix),(1.1/24));
      clear ig1 ig2
      rgs=t.date(ixes);
      yrix=find(get_year(rgs)==yr);
      if (numel(yrix)<(305*24));
        disp(yr);
        ys(ys==yr)=[];
      end;
    end;
  end;
  t=[]; cst=[];
end;
clear st cst;
disp(ys);





[w,f,m,l,k,s,d] = intersect_tses(lkwf1.(fld),...
                                 fwyf1.(fld),...
                                 mlrf1.(fld),...
                                 lonf1.(fld),...
                                 smkf1.(fld),...
                                 sanf1.(fld),...
                                 dryf1.(fld) ...
                                 );

[f,m,l,k,s,d] = intersect_tses(fwyf1.(fld),...
                               mlrf1.(fld),...
                               lonf1.(fld),...
                               smkf1.(fld),...
                               sanf1.(fld),...
                               dryf1.(fld) ...
                               );





% % SUBSET_TS calls could be used to balance out annual cycle
%                                subset_ts(smkf1.ndbc_sea_t,@(x)(x.date<datenum(2005,8,1))),...
%                                subset_ts(dryf1.ndbc_sea_t,@(x)(get_year(x.date)~=2004)) ...

% [f,m,l,k,s,d] = intersect_tses(fwyf1.ndbc_sea_t,...
%                                dryf1.ndbc_sea_t ...

% [f,m,l,k,s] = intersect_tses(fwyf1.ndbc_sea_t,...
%                                smkf1.ndbc_sea_t,...
%                                sanf1.ndbc_sea_t ...

[f,m,l] = intersect_tses(fwyf1.ndbc_sea_t,...
                               mlrf1.ndbc_sea_t,...
                               lonf1.ndbc_sea_t ...
                               );

% % Check seasonal cycle and interannual bias
% %fmg; hist(get_week(m.date),52); xlim([1,52]); titlename('Data Density by Week');
% fmg; hist(get_month(m.date),12); xlim([1,12]); titlename('Data Density by Month');
% fmg; hist(get_season(m.date),4); xlim([1, 4]); titlename('Data Density by Season');
% fmg; hist(get_year(m.date),numel(unique(get_year(m.date)))); titlename('Data Density by Year');









   case 'LKWF1',
    qlh_adj = 0.90;
    doWarms = [false];
    kds = { ...
        %[.475,1.275, 45] ...
        %[.375,1.375,  0] ...
        %[.300,1.250,320] ...
        %[.375,1.375,320] ...
        %[0.150,0.500,354] ...
        %[0.150,0.500,  0] ...
        %[0.150,0.500, 45] ...
        %[0.100,0.400, 91] ...
        %[0.050,0.350, 45] ...		% No good!
        %[0.100,0.350, 91] ...
        %[0.050,0.400, 91] ...
        %[0.050,0.375, 91] ...
        %[0.050,0.350, 91] ...		% Quite good:   0.4,5.2,1.2
        %[0.050,0.350, 67] ...
        %[0.050,0.400, 79] ...		% Also not bad: 0.4,7.3,1.1
        %[0.050,0.375, 79] ...		% Even better:  0.4,4.3,1.1
        [0.050,0.400, 67] ...		% BEST SO FAR:  0.4,3.5,1.0
          };
    advfacs = { ...
        [0.00,1.0, 45] ...
              };
    kths = { ...
        [0,20, 45] ...
           };






   case 'LKWF1',
    qlh_adj = 0.90;
    doWarms = [false];
    kds = { ...
        %[.475,1.275, 45] ...
        %[.375,1.375,  0] ...
        %[.300,1.250,320] ...
        %[.375,1.375,320] ...
        %[0.150,0.500,354] ...
        %[0.150,0.500,  0] ...
        %[0.150,0.500, 45] ...
        %[0.100,0.400, 91] ...
        [0.050,0.350, 91] ...		% BEST SO FAR
          };
    advfacs = { ...
        [0.00,1.0, 45] ...
              };
    kths = { ...
        [0,20, 45] ...
           };




if(0)
end;



if (0)
  [lh,ax,fh]=multiplot_datetick(multidts,multidat,...
                                [upper(stn.station_name),' Sea Temperature and Variability Metrics'],...
                                [],fldabbr,[],ylms,10,{'k-','b-','g-','r-','c-'});
  xlim(datenum(2006,[5,9],[2,29]));
  datetick('x',2,'keeplimits');
  tah=annotation('textarrow',[0.2904,0.2904],[0.4140,0.5665]);
  set(tah,'String',{'Eddy Passage'},'LineWidth',2,'Color','r','FontName','TimesNewRoman','FontSize',12);
  %print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'_',mfilename,'_',hcsum,'-2006-eddy.tif']));
end;





fidx = 0;
clear multidts;
clear multidat;
clear fldabbr;
clear ylms;

fidx = fidx + 1;
multidts{fidx} = stn.(ts).date;
multidat{fidx} = stn.(ts).data;
fldabbr{fidx} = 'T_s';
ylms{fidx} = 'default';

fidx = fidx + 1;
multidts{fidx} = stn.(tsvar).date;
multidat{fidx} = stn.(tsvar).data;
fldabbr{fidx} = '\mu_3_d\sigma_1_d(T_s)';
ylms{fidx} = [0,1];

% fidx = fidx + 1;
% multidts{fidx} = stn.(tsavg).date;
% multidat{fidx} = stn.(tsavg).data;
% fldabbr{fidx} = '\mu_3_d(T_s)';
% ylms{fidx} = 'default';

% fidx = fidx + 1;
% multidts{fidx} = stn.(tsanom).date;
% multidat{fidx} = stn.(tsanom).data;
% fldabbr{fidx} = 'anom^m^i^n_3_d(T_s)';
% ylms{fidx} = [0,1];

% fidx = fidx + 1;
% multidts{fidx} = stn.(hcasum).date;
% multidat{fidx} = stn.(hcasum).data;
% fldabbr{fidx} = '\Sigma_3_d|\partial_tT_s|';
% ylms{fidx} = [0,1];

% fidx = fidx + 1;
% multidts{fidx} = stn.(qasum).date;
% multidat{fidx} = stn.(qasum).data;
% fldabbr{fidx} = '\Sigma_3_d|Q_0/\rhoC_ph|';
% ylms{fidx} = [0,1];

fidx = fidx + 1;
multidts{fidx} = stn.(ahcsum).date;
multidat{fidx} = stn.(ahcsum).data;
fldabbr{fidx} = '|\Sigma_3_d\partial_tT_s|';
ylms{fidx} = [0,1];

fidx = fidx + 1;
multidts{fidx} = stn.(aqsum).date;
multidat{fidx} = stn.(aqsum).data;
fldabbr{fidx} = '|\Sigma_3_dQ_0/\rhoC_ph|';
ylms{fidx} = [0,1];









[lh,ax,fh]=multiplot_datetick(multidts,multidat,...
                   [upper(stn.station_name),' Sea Temperature and Variability Metrics'],...
                   [],fldabbr,[],ylms,10,{'k-','b-','g-','r-','c-'});
set(ax(2),'ylim','default');
set(ax(3),'ylim','default');
set(ax(4),'ylim','default');
print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'_',mfilename,'_',hcsum,'.tif']));

xlim(datenum(2006,[5,9],[2,29]));
datetick('x',2,'keeplimits');
set(ax(2),'ylim',[0,1]);
set(ax(3),'ylim',[0,1]);
set(ax(4),'ylim',[0,1]);

tah=annotation('textarrow',[0.2904,0.2904],[0.4140,0.5665]);
set(tah,'String',{'Eddy Passage'},'LineWidth',2,'Color','r','FontName','TimesNewRoman','FontSize',12);

print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'_',mfilename,'_',hcsum,'-2006-eddy.tif']));





%% NOTE: Taking absolute value here


fidx = fidx + 1;
multidts{fidx} = stn.(hcasum).date;
multidat{fidx} = stn.(hcasum).data;
fldabbr{fidx} = '|\Sigma_3_d\partial_tT_s|';
ylms{fidx} = [0,1];

fidx = fidx + 1;
multidts{fidx} = stn.(qasum).date;
multidat{fidx} = stn.(qasum).data;
fldabbr{fidx} = '|\Sigma_3_dQ_0/\rhoC_ph|';
ylms{fidx} = [0,1];




1;

hcrmsec=[];
flds={};
for cfld=fieldnames(s)';
  fld=cfld{:};
  if ( isfield(s.(fld),'hcrmsec') )
    disp(sprintf('% 8.2f\t%s',s.(fld).hcrmsec,fld));
    hcrmsec(end+1)=s.(fld).hcrmsec;
    flds{end+1}=fld;
  end;
end

fmg;
plot(hcrmsec,1:numel(hcrmsec),'.-');
set(gca,'ytick',1:numel(hcrmsec),'yticklabel',fieldnames(s),'xscale','log');
xlim([0.1,200]);
titlename([upper(s.(fld).station_name),' Sensitivity Analysis']);
print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(s.(fld).station_name),'-anhcrmsec.tif']));





  [stn.(dly_srfld).data,stn.(dly_srfld).date] = grp_ts(stn.(srfld).data,stn.(srfld).date,@floor,@nansum,24);
  [stn.(dly_srtfld).data,stn.(dly_srtfld).date] = grp_ts(stn.(srtfld).data,stn.(srtfld).date,@floor,@nansum,24);
  [stn.(dly_srfld).data,stn.(dly_srfld).date] = grp_ts(stn.(srfld).data,stn.(srfld).date,@floor,@nansum,24);
  [stn.(dly_asrtfld).data,stn.(dly_asrtfld).date] = grp_ts(stn.(asrtfld).data,stn.(asrtfld).date,@floor,@nansum,24);
  [stn.(dly_lrfld).data,stn.(dly_lrfld).date] = grp_ts(stn.(lrfld).data,stn.(lrfld).date,@floor,@nansum,24);
  [stn.(dly_lrtfld).data,stn.(dly_lrtfld).date] = grp_ts(stn.(lrtfld).data,stn.(lrtfld).date,@floor,@nansum,24);
  [stn.(dly_qlhfld).data,stn.(dly_qlhfld).date] = grp_ts(stn.(qlhfld).data,stn.(qlhfld).date,@floor,@nansum,24);
  [stn.(dly_qlhtfld).data,stn.(dly_qlhtfld).date] = grp_ts(stn.(qlhtfld).data,stn.(qlhtfld).date,@floor,@nansum,24);
  [stn.(dly_qshfld).data,stn.(dly_qshfld).date] = grp_ts(stn.(qshfld).data,stn.(qshfld).date,@floor,@nansum,24);
  [stn.(dly_qshtfld).data,stn.(dly_qshtfld).date] = grp_ts(stn.(qshtfld).data,stn.(qshtfld).date,@floor,@nansum,24);
  [stn.(dly_q0fld).data,stn.(dly_q0fld).date] = grp_ts(stn.(q0fld).data,stn.(q0fld).date,@floor,@nansum,24);
  [stn.(dly_qtfld).data,stn.(dly_qtfld).date] = grp_ts(stn.(qtfld).data,stn.(qtfld).date,@floor,@nansum,24);
  [stn.(dly_sq0fld).data,stn.(dly_sq0fld).date] = grp_ts(stn.(sq0fld).data,stn.(sq0fld).date,@floor,@nansum,24);
  [stn.(dly_sqtfld).data,stn.(dly_sqtfld).date] = grp_ts(stn.(sqtfld).data,stn.(sqtfld).date,@floor,@nansum,24);
  [stn.(dly_bq0fld).data,stn.(dly_bq0fld).date] = grp_ts(stn.(bq0fld).data,stn.(bq0fld).date,@floor,@nansum,24);
  [stn.(dly_bq0tfld).data,stn.(dly_bq0tfld).date] = grp_ts(stn.(bq0tfld).data,stn.(bq0tfld).date,@floor,@nansum,24);
  [stn.(dly_bdTfld).data,stn.(dly_bdTfld).date] = grp_ts(stn.(bdTfld).data,stn.(bdTfld).date,@floor,@nansum,24);
  [stn.(dly_bdTffld).data,stn.(dly_bdTffld).date] = grp_ts(stn.(bdTffld).data,stn.(bdTffld).date,@floor,@nansum,24);
  [stn.(dly_hcdTdt).data,stn.(dly_hcdTdt).date] = grp_ts(stn.(hcdTdt).data,stn.(hcdTdt).date,@floor,@nansum,24);
  [stn.(dly_hcdTdtf).data,stn.(dly_hcdTdtf).date] = grp_ts(stn.(hcdTdtf).data,stn.(hcdTdtf).date,@floor,@nansum,24);








  %% Print a brief report (as for Table 4 of Gramer & Mariano "... Heat Budget") 
  dump_robust_fit(stn,dly_dsfld); % HEADER LINE
  dump_robust_fit(stn,dly_dsfld,dly_srtfld);
  dump_robust_fit(stn,dly_dsfld,dly_asrtfld);
  dump_robust_fit(stn,dly_dsfld,dly_lrtfld);
  dump_robust_fit(stn,dly_dsfld,dly_qlhtfld);
  dump_robust_fit(stn,dly_dsfld,dly_qshtfld);
  dump_robust_fit(stn,dly_dsfld,dly_sqtfld);
  dump_robust_fit(stn,dly_dsfld,dly_bq0tfld);
  dump_robust_fit(stn,dly_dsfld,dly_bdTfld);
  dump_robust_fit(stn,dly_dsfld,dly_hcdTdt);

ORIGINAL DUMP_ROBUST_FIT.m:
function [B,Stats] = dump_robust_fit(stn,fitsfld,fitqfld)
%function [B,Stats] = dump_robust_fit(stn,fitsfld,fitqfld)
%
% Dump to the Command Window (DISP(SPRINTF(...))), statistics (R2, SI, A, B,
% RMSE) from a robust linear least-squares fit between FITSFLD and FITQFLD.
% Scatter Index (SI) in this case is the RMSE / Std. Dev. of STN.(FITSFLD); B
% is "slope error" based on regression slope beta, i.e., 100% x abs(1-beta).
%
% Last Saved Time-stamp: <Fri 2013-05-10 19:09:00 Eastern Daylight Time gramer>

  if ( ~exist('fitqfld','var') || isempty(fitqfld) )
    disp(sprintf('%-50s %-4s %-6s %-8s %-8s %-7s',['VS ',fitsfld],'R2','SI%','bias','slope%','RMSE'));
  else
    [B,Stats] = scatter_fit_ts(stn.(fitqfld),stn.(fitsfld),[],[],[],[],'none');
    % disp({fitqfld,Stats.R2_2,num2str(100*Stats.s/25,'%g%%'),B(1),B(2),Stats.s,});
    disp(sprintf('%-50s %4.2f %5.1f%% %+6.3f %7.3f%% %7.3f',...
                 fitqfld,Stats.R2_2,abs(100*Stats.s/nanstd(stn.(fitsfld).data)),...
                 B(1),100*abs(1-B(2)),Stats.s));
  end;

  if ( nargout < 1 )
    B=[]; Stats=[]; clear B Stats;
  end;

return;







  %stn = station_heat_flux_term(stn,srfld,[srfld,'_term'],sfld,sal,mhfld);
  stn = station_heat_flux_term(stn,srfld,srtfld,sfld,sal,mhfld);

      %stn = station_heat_flux_term(stn,asrfld,[asrfld,'_term'],sfld,sal,mhfld);
      stn = station_heat_flux_term(stn,asrfld,asrtfld,sfld,sal,mhfld);

    %stn = station_heat_flux_term(stn,lrfld,[lrfld,'_term'],sfld,sal,mhfld);
    stn = station_heat_flux_term(stn,lrfld,lrtfld,sfld,sal,mhfld);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function optimize_station_heat_budget_disp_robust_fit(stn,fitsfld,fitqfld)
  [B,Stats] = scatter_fit_ts(stn.(fitqfld),stn.(fitsfld),[],[],[],[],'none');
  % disp({fitqfld,Stats.R2_2,num2str(100*Stats.s/25,'%g%%'),B(1),B(2),Stats.s,});
  disp(sprintf('%-50s %3.2g %4.0g%% %+7.3g %7.3g %7.3g',fitqfld,Stats.R2_2,100*Stats.s/nanmean(stn.(fitsfld).data),B(1),B(2),Stats.s));
return;







  [res.qbo.data,res.qbo.date] = grp_ts(qbo.data,qbo.date,per,sumfun,minN);

  [res.udT.data,res.udT.date] = grp_ts(udT.data,udT.date,per,sumfun,minN);

  [res.kd2T.data,res.kd2T.date] = grp_ts(kd2T.data,kd2T.date,per,sumfun,minN);




  % Special handling for huge advection swings???
  [res.raw_udT.data,res.raw_udT.date] = grp_ts(udT.data,udT.date,per,sumfun,minN);
  res.udT.date=res.raw_udT.date; res.udT.data=N.*cumsum(res.raw_udT.data);

  [res.raw_kd2T.data,res.raw_kd2T.date] = grp_ts(kd2T.data,kd2T.date,per,sumfun,minN);
  res.kd2T.date=res.raw_kd2T.date; res.kd2T.data=N.*cumsum(res.raw_kd2T.data);




  [res.t,res.dtf,res.q0,res.qt,res.sq0,res.sqt,res.bq0,res.bq0t,res.udT,res.kd2T,res.bdT,...
   res.hcdTdt,res.sr,res.asr,res.lr,res.qlh,res.qsh,res.qrh,res.qbo] = ...
      intersect_tses(res.t,res.dtf,res.q0,res.qt,res.sq0,res.sqt,res.bq0,res.bq0t,res.udT,res.kd2T,res.bdT,...
                     res.hcdTdt,res.sr,res.asr,res.lr,res.qlh,res.qsh,res.qrh,res.qbo);




  plot(res.udT.date,res.udT.data,'k-','LineWidth',1.5,'Color',[.5,.5,.5]);
  plot(res.kd2T.date,res.kd2T.data,'k:','LineWidth',1.5,'Color',[.5,.5,.5]);


  plot(res.udT.date(mrkix),res.udT.data(mrkix),'k+','LineWidth',1.5,'Color',[.5,.5,.5]);
  plot(res.kd2T.date(mrkix),res.kd2T.data(mrkix),'k+','LineWidth',1.5,'Color',[.5,.5,.5]);


  plh(end+1)=plot(res.udT.date(1),res.udT.data(1),'k-+','LineWidth',1.5,'Color',[.5,.5,.5]);
  legs(end+1) = {'G&M: F_qu\bullet\nablaT_s'};
  plh(end+1)=plot(res.kd2T.date(1),res.kd2T.data(1),'k:+','LineWidth',1.5,'Color',[.5,.5,.5]);
  legs(end+1) = {'G&M: K_H\nabla^2T_s'};










  plh(end+1:end+6) = ...
      plot(res.radif.date(1),res.radif.data(1),'rs-');
           res.turif.date(1),res.turif.data(1),'bs-');
           res.aradif.date(1),res.aradif.data(1),'m^-');
           res.qbo.date(1),res.qbo.data(1),'k-.+');
           res.climradif.date(1),res.climradif.data(1),'ro:');
           res.climturif.date(1),res.climturif.data(1),'bo:');





if ( ~isfield(stn.fkeys_hycom_seatemp_field,'gradient_xs') )
    x=[]; clear x;
    x.x.lat=stn.fkeys_hycom_seatemp_field.lat;
    x.x.lon=stn.fkeys_hycom_seatemp_field.lon;
    x.x.u_field=stn.fkeys_hycom_seatemp_field.gradient_x;
    x.x.v_field=stn.fkeys_hycom_seatemp_field.gradient_y;

    stn2 = station_reorient_field(stn,stn.isobath_orientation,'fkeys_hycom_seatemp_field','gradient_x','gradient_y','gradient_xs','gradient_ls');

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






set(get(lh,'xlabel'),'String','m\bullets^-^1','Rotation',0,'VerticalAlignment','middle');



set(lh,'FontSize',7,'YAxisLocation','right');
set(get(lh,'ylabel'),'String','m\bullets^-^1','Rotation',0,'HorizontalAlignment','center');



lh=colorbar('East');
set(lh,'FontSize',7);
set(get(lh,'ylabel'),'String','m\bullets^-^1','Rotation',0);
text(dts(end)+0.8, 7,'Onshore','Rotation',270,'FontSize',9);
text(dts(end)+0.8,17,'Offshore','Rotation',270,'FontSize',9);

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
%set(lh,'FontSize',7);
% text(datenum(2008,10,22),-2.8,'Daily warming','FontSize',10);
% text(datenum(2008,10,28),-2.8,'Cold front #1','FontSize',10);
% text(datenum(2008,11,17),-2.8,'Cold front #2','FontSize',10);




%print('-dtiff',fullfile(get_thesis_path('../figs'),[mfilename,'-',datestr(dts(1),'yyyy-mm-dd'),'-',datestr(dts(end),'yyyy-mm-dd')]));






          'Location','SouthEast', 'Orientation','horizontal');




text(datenum(2008,10,21,12,0,0),1,'Cooling','FontSize',14);
text(datenum(2008,10,24, 0,0,0),1,'Daily warming','FontSize',14);
text(datenum(2008,10,30, 0,0,0),1,'Strong cold front','FontSize',14);
%text(datenum(2008,11,18),1,'Cold front #2','FontSize',14);




      scatter_fit_ts(stn.b_ndbc_erai_erai_30a_net_flux_term_dly,stn.optim.daydt,[],[],'(Q_0+Q_b)/\rhoC_ph','\Delta_1_dT_s',[],[],true);
      scatter_fit_ts(stn.b_ndbc_erai_erai_30a_avhrr_dt_dly,stn.optim.daydt,[],[],'u^.\nabla_hT_s+K_H\nabla_h^2T_s+(Q_0+Q_b)/\rhoC_ph','\Delta_1_dT_s',[],[],true);
      scatter_fit_ts(stn.optim.dayq,stn.optim.daydt,[],[],'\partial_tT_s','\Delta_1_dT_s',[],[],true);





  %q0s  = [-500:50:500];
  %q0s  = [-500:50:0];
  q0s  = [-200:50:0];
  q0s = linspace
  %bets = [0.002,0.02,0.20];
  %bets = 2.*logspace(-3,-1,20);
  bets = 2.*logspace(-3,-0.8,20);




%{
  [ig,zx] = max(dTdthc_SS); [ix,jx] = ind2sub(size(dTdthc_SS),zx);
  text(q0s(ix),bets(jx),dTdthc_SS(jx,ix),'SS');
  [ig,zx] = max(dTdthc_US); [ix,jx] = ind2sub(size(dTdthc_US),zx);
  text(q0s(ix),bets(jx),dTdthc_US(jx,ix),'US');
  [ig,zx] = max(dTdthc_SU); [ix,jx] = ind2sub(size(dTdthc_SU),zx);
  text(q0s(ix),bets(jx),dTdthc_SU(jx,ix),'SU');
  [ig,zx] = max(dTdthc_UU); [ix,jx] = ind2sub(size(dTdthc_UU),zx);
  text(q0s(ix),bets(jx),dTdthc_UU(jx,ix),'UU');
%}




  clear dTdthc_SS_c dTdthc_US_c dTdthc_SU_c dTdthc_UU_c

          dTdthc_SS_c(bix,qix,1:3) = [0,0,1].*dTdthc_SS(bix,qix);
          dTdthc_US_c(bix,qix,1:3) = [0,1,0].*dTdthc_US(bix,qix);
          dTdthc_SU_c(bix,qix,1:3) = [0,1,1].*dTdthc_SU(bix,qix);
          dTdthc_UU_c(bix,qix,1:3) = [1,0,0].*dTdthc_UU(bix,qix);




(bix,qix)
(bix,qix)
(bix,qix)
(bix,qix)






  if ( ~exist('minN','var') || strcmpi(minN,'default') )
    switch (lower(char(cf))),
     case 'get_yearhour',	minN=[]; % We have no guess of the sampling frequency!
     case 'floor',		minN=23;
     case 'get_yeartriad',	minN=(24*2)+(23*1);
     case 'get_yearpentad',	minN=(24*3)+(23*2);
     case 'get_yearweek',	minN=(24*4)+(23*3);
     case 'get_yearmonth',	minN=24*((30*7)+(29*4)+(27));
     case 'get_yearseason',	minN=24*87;
     case 'get_year',		minN=24*340;
     otherwise,			minN=[];
    end;
  end;







%%%%DEBUG:  [s.ts,s.td,s.q,s.bq,s.dt,s.hc]=intersect_tses(s.raw_ts,s.raw_td,s.raw_q,s.raw_bq,s.raw_dt,s.raw_hc);
  [s.ts,s.q,s.bq,s.dt,s.hc]=intersect_tses(s.raw_ts,s.raw_q,s.raw_bq,s.raw_dt,s.raw_hc);
  [ig,s.td]=intersect_tses(s.raw_ts,s.raw_td);



'\partial_tT 37m', ...




  stn.opts.grid_interp_method = get_opt(stn.opts,'grid_interp_method','linear');
  disp(['** Using ',stn.opts.grid_interp_method,' grid interpolation **']);


  %%%% ??? DEBUG
  % disp('%% Forcing use of ERAI winds **');
  % Wfld='erai_wind_speed';
  % Dfld='erai_wind_dir';

  % disp('%% Forcing use of MLRF2 in situ PAR -> insolation **');
  % dsrfld='bic_surf_dsrf';


  if ( ~strcmp(sfld,[ISPFX '_sea_t']) )
    disp(['Using sea temperature ',sfld]);
  end;











   case 'SEASONAL_UNSTEADY9',
    res.u = res.u_US;
    R = repmat(R,size(res.u));
    seasix = find(ismember(get_month(dts),[1:4,9:12]));
    res.u(seasix) = res.u_UU(seasix);




  % Create fields for comparison of daily averages, sums, and changes
%%%%%%%%%%%%%%%%%%%%%%% TEMPORARILY COMMENT OUT???
%{
  [stn.(dly_sfld).data,stn.(dly_sfld).date] = grp_ts(stn.(sfld).data,stn.(sfld).date,@floor,@nanmean,23);
  stn.(dly_dsfld).date = stn.(dly_sfld).date(2:end);
  stn.(dly_dsfld).data = diff(stn.(dly_sfld).data);
  stn = filter_gaps(stn,dly_sfld,dly_dsfld,1.5);
  stn = station_heat_flux_term_inverse(stn,dly_dsffld,dly_dsfld,sfld,sal,mhfld);
  [stn.(dly_q0fld).data,stn.(dly_q0fld).date] = grp_ts(stn.(q0fld).data,stn.(q0fld).date,@floor,@nansum,24);
  [stn.(dly_qtfld).data,stn.(dly_qtfld).date] = grp_ts(stn.(qtfld).data,stn.(qtfld).date,@floor,@nansum,24);
  [stn.(dly_sq0fld).data,stn.(dly_sq0fld).date] = grp_ts(stn.(sq0fld).data,stn.(sq0fld).date,@floor,@nansum,24);
  [stn.(dly_sqtfld).data,stn.(dly_sqtfld).date] = grp_ts(stn.(sqtfld).data,stn.(sqtfld).date,@floor,@nansum,24);
  [stn.(dly_bq0fld).data,stn.(dly_bq0fld).date] = grp_ts(stn.(bq0fld).data,stn.(bq0fld).date,@floor,@nansum,24);
  [stn.(dly_bq0tfld).data,stn.(dly_bq0tfld).date] = grp_ts(stn.(bq0tfld).data,stn.(bq0tfld).date,@floor,@nansum,24);
  [stn.(dly_bdTfld).data,stn.(dly_bdTfld).date] = grp_ts(stn.(bdTfld).data,stn.(bdTfld).date,@floor,@nansum,24);
  [stn.(dly_bdTffld).data,stn.(dly_bdTffld).date] = grp_ts(stn.(bdTffld).data,stn.(bdTffld).date,@floor,@nansum,24);
  [stn.(dly_hcdTdt).data,stn.(dly_hcdTdt).date] = grp_ts(stn.(hcdTdt).data,stn.(hcdTdt).date,@floor,@nansum,24);
  [stn.(dly_hcdTdtf).data,stn.(dly_hcdTdtf).date] = grp_ts(stn.(hcdTdtf).data,stn.(hcdTdtf).date,@floor,@nansum,24);
%}









    {'_24_h_lp','SEASONAL_UNSTEADY',1.00,{[0.050,0.150, 45],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','SU',               1.00,{[0.050,0.150, 45],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','SEASONAL_UNSTEADY',1.00,{[0.050,0.150, 69],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','SU',               1.00,{[0.050,0.150, 69],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','SEASONAL_UNSTEADY',1.00,{[0.050,0.200, 45],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','SU',               1.00,{[0.050,0.200, 45],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','SEASONAL_UNSTEADY',1.00,{[0.050,0.200, 69],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','SU',               1.00,{[0.050,0.200, 69],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2




    %appendtitlename([' ',strrep(lp,'_','\_'),' ',sc,',wf=',num2str(wf)]);




  switch ( lower(stn.station_name) ),
   case 'looe1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
    end;
   case 'mlrf1',
    if ( ~isempty(regexp(sfld,'sea_t')) )
    end;
   case 'smkf1',
    if ( ~isempty(regexp(sfld,'sea_t')) )
      %%%% 2013 Mar 09
      %%%%bad_years=[2004];
    end;
  end;




reopt_wf_fwyf1.m:
    {'_24_h_lp','SU',0.66,{[.025,.225,102],},{[0,0.5,102]},{[   0, 200, 45]},},... %err  3.6,0.3
    {'_24_h_lp','SU',0.66,{[.025,.225,102],},{[0,0.5,102]},{[   0, 400, 45]},},... %err  3.6,0.3
    {'_24_h_lp','SU',0.66,{[.025,.225,102],},{[0,0.5,102]},{[ 200, 400, 45]},},... %err  3.6,0.3
    {'_24_h_lp','SU',0.66,{[.025,.225,102],},{[0,0.5,102]},{[   0,2000, 45]},},... %err  3.6,0.3
    {'_24_h_lp','SU',0.66,{[.025,.225,102],},{[0,0.5,102]},{[   0,4000, 45]},},... %err  3.6,0.3
    {'_24_h_lp','SU',0.66,{[.025,.225,102],},{[0,0.5,102]},{[2000,4000, 45]},},... %err  3.6,0.3



    %appendtitlename([' ',strrep(lp,'_','\_'),' ',sc,',wf=',num2str(wf)]);




      'FWYF1'
      stn.opts.hc_scaling = 'SU';
      stn.opts.warming_factor = 0.66;

      'MLRF1'
      stn.opts.hc_scaling = 'SU';
      stn.opts.warming_factor = 0.66;

      'SMKF1'
      stn.opts.hc_scaling = 'SU';
      stn.opts.warming_factor = 0.66;






  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	056	5	06/03/2000 1901	11/09/2000 1406	good
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	TL23	5	11/09/2000 1402	06/19/2001 1527	good
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	Tl46	5	06/19/2001 1527	10/23/2001 1813	good
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	TL92	5	10/23/2001 1813	07/10/2002 1725	
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	TL83	5	07/10/2002 1712	11/21/2002 1628	

  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	046	5	06/03/2000 1847	11/08/2000 2056	no data
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL92	5	11/08/2000 2049	06/19/2001 1539	good
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL49	5	06/19/2001 1539	10/23/2001 1821	good
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL23	5	10/23/2001 1821	07/10/2002 1744	
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL82	5	07/10/2002 1736	11/21/2002 1641	

  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	049	5	06/03/2000 1830	11/08/2000 2030	good
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL17	5	11/08/2000 2024	06/19/2001 1551	good
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL56	5	06/19/2001 1551	10/23/2001 1829	good
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL26	5	10/23/2001 1829	07/10/2002 1800	
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL57	5	07/10/2002 1755	11/20/2002 1923	


  25.06900
 -80.31950
  25.07383
 -80.32450
  25.07950
 -80.33433
  25.08783
 -80.34717



TL1
TL2
TL3
Tl4

25 04.14;	-80 19.17;
25 04.43;	-80 19.47;
25 04.77;	-80 20.06;
25 05.27;	-80 20.83;
  25.068999999999999
 -80.319500000000005
  25.073833333333333
 -80.324500000000000
  25.079499999999999
 -80.334333333333333
  25.087833333333332
 -80.347166666666666








  stns.b.lat1 =  25.0932;
  stns.b.lon1 = -80.3550;
  stns.b.lat2 =  25.1090;
  stns.b.lon2 = -80.3803;

  stns.csh.lat1 =  25.0740;
  stns.csh.lon1 = -80.3178;
  stns.csh.lat2 =  25.0733;
  stns.csh.lon2 = -80.3183;
  stns.csh.lat3 =  25.0673;
  stns.csh.lon3 = -80.3183;

  stns.cdp.lat1 =  25.0740;
  stns.cdp.lon1 = -80.3178;
  stns.cdp.lat2 =  25.0733;
  stns.cdp.lon2 = -80.3183;
  stns.cdp.lat3 =  25.0673;
  stns.cdp.lon3 = -80.3183;




  emix=find(arrayfun(@(x)(isempty(x.mag)|isempty(x.dir)|isempty(x.u)|isempty(x.v)|isnan(x.mag)|isnan(x.dir)|isnan(x.u)|isnan(x.v)),vbn));




  n = numel(vbn)/50;
  if ( n ~= floor(n) )
    error('Shape error for currents profiles');
  end;
  nemix=find(arrayfun(@(x)(~isempty(x.u)&~isempty(x.v)),vbn));

  stn.adcp_speed.prof = repmat(nan,[n,50]);
  dat = reshape([vbn.mag]',[50,n])';
  stn.adcp_speed.prof(nemix) = dat(nemix)./100;
  dat=[]; clear dat;

  stn.adcp_dir.prof = repmat(nan,[n,50]);
  dat = reshape([vbn.dir]',[50,n])';
  stn.adcp_dir.prof(nemix) = dat(nemix)./100;
  dat=[]; clear dat;

  stn.adcp_u.prof = repmat(nan,[n,50]);
  dat = reshape([vbn.u]',[50,n])';
  stn.adcp_u.prof(nemix) = dat(nemix)./100;
  dat=[]; clear dat;

  stn.adcp_v.prof = repmat(nan,[n,50]);
  dat = reshape([vbn.v]',[50,n])';
  stn.adcp_v.prof(nemix) = dat(nemix)./100;
  dat=[]; clear dat;

  vbn=[]; clear vbn;




  dr = repmat(nan,[n,50]);
  dr(nemix) = [vbn.dir];
  stn.adcp_dir.prof = dr;
  stn.adcp_u.prof = reshape([vbn.u],[50,n])';
  stn.adcp_v.prof = reshape([vbn.v],[50,n])';
  vbn=[]; clear vbn;





  n = numel(vbn)/50;
  if ( n ~= floor(n) )
    error('Shape error for currents profiles');
  end;
  spd = repmat(nan,[n,50]);
  spd(nmix) = reshape([vbn.mag],[50,n])';
  stn.adcp_speed.prof = reshape([vbn.mag],[50,n])';
  stn.adcp_dir.prof = reshape([vbn.dir],[50,n])';
  stn.adcp_u.prof = reshape([vbn.u],[50,n])';
  stn.adcp_v.prof = reshape([vbn.v],[50,n])';
  vbn=[]; clear vbn;



    p=((max(dat)-min(dat))*sin((dts-(13/24))*2*pi))+min(dat);
    p(p<0)=0;
    p = p + min(dat);





  for ix=1:2:numel(substitute_field_names)
    % s = [substitute_field_names{ix},'=',substitute_field_names{ix+1}];
    % disp(['Eval ',s]);
    % eval(s);
    assignin('caller',substitute_field_names{ix},substitute_field_names{ix+1});
  end;




  legend('T_3','T_7','Mrt T_7','u^x^s_7-u^x^s_2_2','T_2_2','Mrt T_2_2','Q_0/100','P_2_2','U^x^s', 'Location','SouthEast');




  fmg;
  yr=2008; mo=7; dys=1:20;
  for dix=1:numel(dys)
    dy = dys(dix);
    fs = @(x)(find(datenum(yr,mo,dy)==floor(x.date)));
    t0 = stn.ndbc_sea_t.data(fs(stn.ndbc_sea_t)); t0=t0(1);
    subplot_tight(floor(sqrt(numel(dys))),ceil(sqrt(numel(dys))),dix);
    hold on;
    ts=subset_ts(klgf1.cm_deep_udTdz,fs); plot(get_hour(ts.date),ts.data.*10,'b-');
    %ts=subset_ts(klgf1.cm_duTdz,fs); plot(get_hour(ts.date),ts.data.*10,'b-');
    %ts=subset_ts(klgf1.cm_deep_seatemp,fs); plot(get_hour(ts.date),ts.data-nanmean(ts.data),'r-');
    % 0.1K/hour warming at MLRF1 == approx. 1200 W/m^2 net surface heat flux
    ts=subset_ts(stn.simple_ndbc_erai_erai_30a_net_flux,fs); plot(get_hour(ts.date),ts.data./1000,'m-');
    ts=subset_ts(stn.ndbc_sea_t,fs); plot(get_hour(ts.date),ts.data-t0,'k-');
    ts=subset_ts(stn.bic_surf_par,fs); plot(get_hour(ts.date),ts.data./1500,'-','Color',[.9,.9,0]);
    %ts=subset_ts(stn.ncep_par,fs); plot(get_hour(ts.date)+2,ts.data./1500,'-','Color',[0,.5,0]);
    %ts=subset_ts(stn.ncep_dsrf,fs); plot(get_hour(ts.date)+2,ts.data./1500,'-','Color',[0,.5,0]);
    %ts=subset_ts(stn.ndbc_erai_erai_30a_wind_stress_xshore,fs); plot(get_hour(ts.date),ts.data,'c-');
    ts=subset_ts(stn.ndbc_erai_erai_30a_wind_stress,fs); plot(get_hour(ts.date),ts.data,'c-');
    %axis([0,23,-10,+10]);
    axis([0,23,-1.2,+1.2]);
    grid on;
    xlabel(num2str(dy));
  end;
  suptitlename(['NCORE vs. MLRF1 ',datestr(datenum(yr,mo,dys(1))),' - ',datestr(datenum(yr,mo,dys(end)))]);



  fmg;
  yr=2007; mo=7; dys=1:45;
  yr=2008; mo=7; dys=1:45;
  %yr=2007; mo=10; dys=15:60;
  for dix=1:numel(dys)
    dy = dys(dix);
    %  subplot_tight(floor(sqrt(numel(dys))),ceil(sqrt(numel(dys))),dix);
    %  hold on;
    fs = @(x)(find(datenum(yr,mo,dy)==floor(x.date)));
    %ts=subset_ts(klgf1.cm_deep_udTdz,fs); plot(get_hour(ts.date),ts.data.*10,'b-');
    %ts=subset_ts(klgf1.cm_duTdz,fs); plot(get_hour(ts.date),ts.data.*10,'b-');
    ts=subset_ts(klgf1.cm_deep_seatemp,fs); plot(get_hour(ts.date),ts.data-nanmean(ts.data),'r-');
    %%%%ts=subset_ts(klgf1.cm_dTdz,fs); plot(get_hour(ts.date),[0;diff(cumsum(ts.data))]./10,'r-');
    %%%ts=subset_ts(klgf1.cm_dTdz,fs); plot(get_hour(ts.date),[0;diff(ts.data)],'r-');
    % 0.1K/hour warming at MLRF1 == approx. 1200 W/m^2 net surface heat flux
    ts=subset_ts(stn.simple_ndbc_erai_erai_30a_net_flux,fs); plot(get_hour(ts.date),ts.data./1000,'m-');
    %ts=subset_ts(stn.ndbc_sea_t,fs); plot(get_hour(ts.date),ts.data-nanmean(ts.data),'k-');
    %ts=subset_ts(stn.bic_surf_par,fs); plot(get_hour(ts.date),ts.data./1500,'-','Color',[.9,.9,0]);
    %ts=subset_ts(stn.ncep_par,fs); plot(get_hour(ts.date)+2,ts.data./1500,'-','Color',[0,.5,0]);
    %ts=subset_ts(stn.ncep_dsrf,fs); plot(get_hour(ts.date)+2,ts.data./1500,'-','Color',[0,.5,0]);
    %ts=subset_ts(stn.ndbc_erai_erai_30a_wind_stress_xshore,fs); plot(get_hour(ts.date),ts.data,'c-');
    %ts=subset_ts(stn.ndbc_erai_erai_30a_wind_stress,fs); plot(get_hour(ts.date),ts.data,'c-');
    %axis([0,23,-10,+10]);
  end;
  axis([0,23,-2,+2]);
  grid on;
  suptitlename(['NCORE vs. MLRF1 ',datestr(datenum(yr,mo,dys(1))),' - ',datestr(datenum(yr,mo,dys(end)))]);







 set(gca,'FontSize',7);





[all_ht,all_hc]=intersect_tses(looe1.adcp_baroclinic_u1dTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
fmg;
boxplot_ts(all_ht,@(d)(datestr(d,3)),'grouporder',allmos,'allcol','k');
boxplot_ts(all_hc,@(d)(datestr(d,3)),'grouporder',allmos,'allcol','r');
ylim([-.1,+.1]); xlabel('Bottom Heat Transport'); grid on; ylabel('^oC/hour');
titlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport climatology']);





[all_ht,all_hc]=intersect_tses(looe1.adcp_udTdz,stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc);
fmg;
spt(2,1,1); boxplot_ts(subset_ts(all_ht,@(x)(find(ismember(get_hour(x.date),hrs)))),@(d)(datestr(d,3))); ylim([-.1,+.1]); xlabel('Heat Transport'); grid on; ylabel('^oC/hour');
spt(2,1,2); boxplot_ts(subset_ts(all_hc,@(x)(find(ismember(get_hour(x.date),hrs)))),@(d)(datestr(d,3))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
suptitlename(['Predicted ',upper(stn.station_name),' horizontal convection vs. observed ',upper(looe1.station_name),' heat transport climatology']);




spt(2,1,1); boxplot_ts(subset_ts(all_ht,@(x)(find(ismember(get_hour(x.date),hrs))))); ylim([-.1,+.1]); xlabel('Heat Transport'); grid on; ylabel('^oC/hour');
spt(2,1,2); boxplot_ts(subset_ts(all_hc,@(x)(find(ismember(get_hour(x.date),hrs))))); ylim([-.1,+.1]); xlabel('Horizontal Convection'); grid on; ylabel('^oC/hour');
p


t0 = [];
  if ( isempty(t0) ); t0 = stn.ndbc_sea_t.data(fs(stn.ndbc_sea_t)); t0=t0(1); end;





fmg; plot_ts(stn.simple_ndbc_erai_erai_30a_net_flux,stn.ndbc_erai_erai_30a_net_flux,stn.b_ndbc_erai_erai_30a_net_flux,stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux,stn.ndbc_sea_t_dly_diff_flux); legend('Simple','Sfc.','Sfc./Btm.','HC','Implied');





[stn.ndbc_sea_t_daily.data,stn.ndbc_sea_t_daily.date] = ...
    grp_ts(stn.ndbc_sea_t.data,stn.ndbc_sea_t.date,@floor,@nanmean,24);
stn.ndbc_sea_t_daily_diff.date = stn.ndbc_sea_t_daily.date(1:end-1);
stn.ndbc_sea_t_daily_diff.data = diff(stn.ndbc_sea_t_daily.data);
stn = filter_gaps(stn,'ndbc_sea_t_daily','ndbc_sea_t_daily_diff',1.5);







klgf1 = verify_variable(klgf1,'cm_deep_seatemp_1_d_var_3_d_avg');
  fmg;
  plot_ts(klgf1.cm_deep_seatemp_1_d_var_3_d_avg,...






  %[sp,dr] = ts_uv_to_spddir(u,v,doWind)

  klgf1.cm_shallow_speed_raw.date = klgf1.cm_shallow_u_raw.date;
  klgf1.cm_shallow_speed_raw.data = uv_to_spd(klgf1.cm_shallow_u_raw.data,klgf1.cm_shallow_v_raw.data);
  klgf1.cm_shallow_dir_raw.date = klgf1.cm_shallow_u_raw.date;
  klgf1.cm_shallow_dir_raw.data = uv_to_dir_curr(klgf1.cm_shallow_u_raw.data,klgf1.cm_shallow_v_raw.data);

  klgf1.cm_deep_speed_raw.date = klgf1.cm_deep_u_raw.date;
  klgf1.cm_deep_speed_raw.data = uv_to_spd(klgf1.cm_deep_u_raw.data,klgf1.cm_deep_v_raw.data);
  klgf1.cm_deep_dir_raw.date = klgf1.cm_deep_u_raw.date;
  klgf1.cm_deep_dir_raw.data = uv_to_dir_curr(klgf1.cm_deep_u_raw.data,klgf1.cm_deep_v_raw.data);

  klgf1.cm_shallow_speed.date = klgf1.cm_shallow_u.date;
  klgf1.cm_shallow_speed.data = uv_to_spd(klgf1.cm_shallow_u.data,klgf1.cm_shallow_v.data);
  klgf1.cm_shallow_dir.date = klgf1.cm_shallow_u.date;
  klgf1.cm_shallow_dir.data = uv_to_dir_curr(klgf1.cm_shallow_u.data,klgf1.cm_shallow_v.data);

  klgf1.cm_deep_speed.date = klgf1.cm_deep_u.date;
  klgf1.cm_deep_speed.data = uv_to_spd(klgf1.cm_deep_u.data,klgf1.cm_deep_v.data);
  klgf1.cm_deep_dir.date = klgf1.cm_deep_u.date;
  klgf1.cm_deep_dir.data = uv_to_dir_curr(klgf1.cm_deep_u.data,klgf1.cm_deep_v.data);


  marf1.cm_shallow_speed_raw.date = marf1.cm_shallow_u_raw.date;
  marf1.cm_shallow_speed_raw.data = uv_to_spd(marf1.cm_shallow_u_raw.data,marf1.cm_shallow_v_raw.data);
  marf1.cm_shallow_dir_raw.date = marf1.cm_shallow_u_raw.date;
  marf1.cm_shallow_dir_raw.data = uv_to_dir_curr(marf1.cm_shallow_u_raw.data,marf1.cm_shallow_v_raw.data);

  marf1.cm_deep_speed_raw.date = marf1.cm_deep_u_raw.date;
  marf1.cm_deep_speed_raw.data = uv_to_spd(marf1.cm_deep_u_raw.data,marf1.cm_deep_v_raw.data);
  marf1.cm_deep_dir_raw.date = marf1.cm_deep_u_raw.date;
  marf1.cm_deep_dir_raw.data = uv_to_dir_curr(marf1.cm_deep_u_raw.data,marf1.cm_deep_v_raw.data);

  marf1.cm_shallow_speed.date = marf1.cm_shallow_u.date;
  marf1.cm_shallow_speed.data = uv_to_spd(marf1.cm_shallow_u.data,marf1.cm_shallow_v.data);
  marf1.cm_shallow_dir.date = marf1.cm_shallow_u.date;
  marf1.cm_shallow_dir.data = uv_to_dir_curr(marf1.cm_shallow_u.data,marf1.cm_shallow_v.data);

  marf1.cm_deep_speed.date = marf1.cm_deep_u.date;
  marf1.cm_deep_speed.data = uv_to_spd(marf1.cm_deep_u.data,marf1.cm_deep_v.data);
  marf1.cm_deep_dir.date = marf1.cm_deep_u.date;
  marf1.cm_deep_dir.data = uv_to_dir_curr(marf1.cm_deep_u.data,marf1.cm_deep_v.data);







  % QA/QC
  baddt = unique(floor(k1lar101.cm_seapres_raw.date(abs(k1lar101.cm_seapres_raw.data-median(k1lar101.cm_seapres_raw.data))>5)));
  k1lar101.cm_u_qc = subset_ts(k1lar101.cm_u_raw,@(x)(find(~ismember(x.date,baddt))));
  k1lar101.cm_v_qc = subset_ts(k1lar101.cm_v_raw,@(x)(find(~ismember(x.date,baddt))));
  k1lar101.cm_w_qc = subset_ts(k1lar101.cm_w_raw,@(x)(find(~ismember(x.date,baddt))));
  k1lar101.cm_seatemp_qc = subset_ts(k1lar101.cm_seatemp_raw,@(x)(find(~ismember(x.date,baddt))));
  k1lar101.cm_seapres_qc = subset_ts(k1lar101.cm_seapres_raw,@(x)(find(~ismember(x.date,baddt))));




  % QA/QC
  baddt = unique(union(...
      floor(klgf1.cm_shallow_u_raw.date(abs(klgf1.cm_shallow_u_raw.data-median(klgf1.cm_shallow_u_raw.data))>5)),...
      floor(klgf1.cm_deep_u_raw.date(abs(klgf1.cm_deep_u_raw.data-median(klgf1.cm_deep_u_raw.data))>5))...
      ));




  %for fld={'cm_u_raw','cm_u_raw','cm_v_raw','cm_v_raw','cm_w_raw','cm_w_raw','cm_seapres_raw','cm_seapres_raw','cm_seatemp_raw','cm_seatemp_raw'};




1;

doCovs = false;

if ( ~exist('s','var') )
  s=[];
end;

% stnms={'fwyf1','mlrf1','lonf1','smkf1','looe1','sanf1','dryf1'};
stnms={'mlrf1','lonf1'};

%for perfun=[@floor,@get_yearmonth]; for seasfilt=[@ts_isfinite,@ts_boreal_warm,@ts_boreal_cool]; for filterBadDates = [false,true];
for perfun=[@floor]; for cseasfilt={@ts_isfinite,@ts_boreal_warm,@ts_boreal_cool}; for filterBadDates = [true];

seasfilt = cseasfilt{:};

diary off;
more off;

switch (char(perfun))
 case 'floor';               cumfun = @nanmean;      minN = 24;      maxgaps=(25/24);
 case 'get_yearweek';        cumfun = @nanmean;      minN = 24*4.5;  maxgaps=8;
 case 'get_yearmonth';       cumfun = @nanmean;      minN = 24*20;   maxgaps=32;
 otherwise,                  error('Unsupported perfun');
end;

%seasfilt = [];
%seasfilt = @ts_boreal_warm;
%seasfilt = @ts_boreal_cool;

switch (char(seasfilt))
 case 'ts_isfinite';         per = 'year';
 case 'ts_boreal_warm';      per = 'warm';
 case 'ts_boreal_cool';      per = 'cool';
 otherwise,                  error('Unsupported seasfilt');
end;

%filterBadDates = false;
%filterBadDates = true;

if ( filterBadDates ); filtstr = 'filtered'; else filtstr = 'unfiltered'; end;

fname = fullfile(get_thesis_path('../data'),[mfilename,'-',char(perfun),'-',char(seasfilt),'-',filtstr,'-results.log']);
diary(fname);

disp(mfilename);
timenow;
disp(char(perfun));
disp(char(seasfilt));
if ( filterBadDates )
  disp('WILL FILTER "BAD DATES"');
else
  disp('Dates unfiltered');
end;

sacs=0;
sars=0;
esacs=0;
esars=0;
swcs=0;
swrs=0;
eswcs=0;
eswrs=0;
sw2cs=0;
sw2rs=0;
esw2cs=0;
esw2rs=0;
cs=0;
rs=0;
windrs=0;
rhrs = 0;

for snmc=stnms
  snm=snmc{:};

  if ( ~isfield(s,snm) )
    s.(snm) = get_station_from_station_name(snm);
    s.(snm) = load_all_ndbc_data(s.(snm));
    s.(snm) = station_optimal_isobath_orientation(s.(snm));
    s.(snm) = verify_variable(s.(snm),'ndbc_wind1_u');
    s.(snm) = verify_variable(s.(snm),'ndbc_wind1_v');
    s.(snm) = station_reorient_vectors(s.(snm),s.(snm).isobath_orientation,...
                                       'ndbc_wind1_u','ndbc_wind1_v');

    [six,aix]=intersect_dates(s.(snm).ndbc_sea_t.date,s.(snm).ndbc_air_t.date);
    s.(snm).ndbc_sea_air_diff.date=s.(snm).ndbc_sea_t.date(six);
    s.(snm).ndbc_sea_air_diff.data=s.(snm).ndbc_sea_t.data(six)-s.(snm).ndbc_air_t.data(aix);

    s.(snm).ndbc_wind2.date = s.(snm).ndbc_wind1_xshore.date;
    s.(snm).ndbc_wind2.data = ((s.(snm).ndbc_wind1_xshore.data .^ 2) + (s.(snm).ndbc_wind1_lshore.data .^ 2));
    s.(snm).ndbc_wind.date = s.(snm).ndbc_wind2.date;
    s.(snm).ndbc_wind.data = sqrt(s.(snm).ndbc_wind2.data);
  end;
  if ( ~isfield(s.(snm),'erai_air_t') )
    x = get_erai_station(snm);
    disp('Adjusting ERAI radiation, waves');
    x = adjust_erai_station(x);
    x = adjust_erai_station_waves(x);
    s.(snm).erai_air_t = x.erai_air_t;
    s.(snm).erai_relhumid = x.erai_relhumid;
    s.(snm).erai_spechumid = x.erai_spechumid;
    s.(snm).erai_wind_speed = x.erai_wind_speed;
    s.(snm).erai_wind_dir = x.erai_wind_dir;
    s.(snm).erai_dsrf = x.erai_dsrf;
    s.(snm).erai_dlrf = x.erai_dlrf;
    s.(snm).erai_precip = x.erai_precip;
    s.(snm).erai_dsrf_adj = x.erai_dsrf_adj;
    s.(snm).erai_dlrf_adj = x.erai_dlrf_adj;
    s.(snm).erai_precip_adj = x.erai_precip_adj;
    s.(snm).erai_sigwavehgt = x.erai_sigwavehgt;
    s.(snm).erai_sigwavehgt_adj = x.erai_sigwavehgt_adj;
    x=[]; clear x;
    s.(snm) = verify_variable(s.(snm),'erai_wind_u');
    s.(snm) = verify_variable(s.(snm),'erai_wind_v');
    s.(snm) = station_reorient_vectors(s.(snm),s.(snm).isobath_orientation,...
                                       'erai_wind_u','erai_wind_v');

    s.(snm).erai_wind2.date = s.(snm).erai_wind_xshore.date;
    s.(snm).erai_wind2.data = ((s.(snm).erai_wind_xshore.data .^ 2) + (s.(snm).erai_wind_lshore.data .^ 2));
    s.(snm).erai_wind.date = s.(snm).erai_wind2.date;
    s.(snm).erai_wind.data = sqrt(s.(snm).erai_wind2.data);
  end;
  if ( ~isfield(s.(snm),'avhrr_weekly_sst') )
    switch (snm)
     case 'sanf1',
      s.(snm) = get_avhrr_weekly_field(s.(snm),true,'nearest',3,~filterBadDates);
     case 'dryf1',
      s.(snm) = get_avhrr_weekly_field(s.(snm),true,{'nanmean',2,2},3,~filterBadDates);
     otherwise,
      s.(snm) = get_avhrr_weekly_field(s.(snm),true,'linear',5,~filterBadDates);
    end;

    s.(snm) = rmfield(s.(snm), grepstruct(s.(snm),'_field'));
    s.(snm) = rmfield(s.(snm), grepstruct(s.(snm),'raw_'));

    s.(snm) = station_reorient_vectors(s.(snm),s.(snm).isobath_orientation,...
                                       'hourly_avhrr_weekly_sst_x','hourly_avhrr_weekly_sst_y');
  end;

  [s.(snm).ts_per.data,s.(snm).ts_per.date]=grp_ts(s.(snm).ndbc_sea_t.data,s.(snm).ndbc_sea_t.date,perfun,cumfun,minN);

  [s.(snm).sst_per.data,s.(snm).sst_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst.data,s.(snm).hourly_avhrr_weekly_sst.date,perfun,cumfun,minN);
  [s.(snm).sst_xs_per.data,s.(snm).sst_xs_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_xshore.data,s.(snm).hourly_avhrr_weekly_sst_xshore.date,perfun,cumfun,minN);
  [s.(snm).sst_ls_per.data,s.(snm).sst_ls_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_lshore.data,s.(snm).hourly_avhrr_weekly_sst_lshore.date,perfun,cumfun,minN);
  try,
    [s.(snm).sst_l_per.data,s.(snm).sst_l_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_l.data,s.(snm).hourly_avhrr_weekly_sst_l.date,perfun,cumfun,minN);
    [s.(snm).sst_dl_per.data,s.(snm).sst_dl_per.date]=grp_ts(s.(snm).hourly_avhrr_weekly_sst_dl.data,s.(snm).hourly_avhrr_weekly_sst_dl.date,perfun,cumfun,minN);
  catch,
  end;

  [s.(snm).ta_per.data,s.(snm).ta_per.date]=grp_ts(s.(snm).ndbc_air_t.data,s.(snm).ndbc_air_t.date,perfun,cumfun,minN);
  [s.(snm).W_kts_per.data,s.(snm).W_kts_per.date]=grp_ts(s.(snm).ndbc_wind1_speed.data,s.(snm).ndbc_wind1_speed.date,perfun,cumfun,minN);
  [s.(snm).U_per.data,s.(snm).U_per.date]=grp_ts(s.(snm).ndbc_wind1_xshore.data,s.(snm).ndbc_wind1_xshore.date,perfun,cumfun,minN);
  [s.(snm).V_per.data,s.(snm).V_per.date]=grp_ts(s.(snm).ndbc_wind1_lshore.data,s.(snm).ndbc_wind1_lshore.date,perfun,cumfun,minN);
  [s.(snm).W_per.data,s.(snm).W_per.date]=grp_ts(s.(snm).ndbc_wind.data,s.(snm).ndbc_wind.date,perfun,cumfun,minN);
  [s.(snm).W2_per.data,s.(snm).W2_per.date]=grp_ts(s.(snm).ndbc_wind2.data,s.(snm).ndbc_wind2.date,perfun,cumfun,minN);


  [s.(snm).e_ta_per.data,s.(snm).e_ta_per.date]=grp_ts(s.(snm).erai_air_t.data,s.(snm).erai_air_t.date,perfun,cumfun,minN);
  [s.(snm).e_qa_per.data,s.(snm).e_qa_per.date]=grp_ts(s.(snm).erai_spechumid.data,s.(snm).erai_spechumid.date,perfun,cumfun,minN);
  [s.(snm).e_W_kts_per.data,s.(snm).e_W_kts_per.date]=grp_ts(s.(snm).erai_wind_speed.data,s.(snm).erai_wind_speed.date,perfun,cumfun,minN);
  [s.(snm).e_U_per.data,s.(snm).e_U_per.date]=grp_ts(s.(snm).erai_wind_xshore.data,s.(snm).erai_wind_xshore.date,perfun,cumfun,minN);
  [s.(snm).e_V_per.data,s.(snm).e_V_per.date]=grp_ts(s.(snm).erai_wind_lshore.data,s.(snm).erai_wind_lshore.date,perfun,cumfun,minN);
  [s.(snm).e_W_per.data,s.(snm).e_W_per.date]=grp_ts(s.(snm).erai_wind.data,s.(snm).erai_wind.date,perfun,cumfun,minN);
  [s.(snm).e_W2_per.data,s.(snm).e_W2_per.date]=grp_ts(s.(snm).erai_wind2.data,s.(snm).erai_wind2.date,perfun,cumfun,minN);

  % NOTE USE OF "@nansum" FOR INSOLATION!
  [s.(snm).e_qswi_per.data,s.(snm).e_qswi_per.date]=grp_ts(s.(snm).erai_dsrf_adj.data,s.(snm).erai_dsrf_adj.date,perfun,@nansum,minN);
  [s.(snm).e_qlwi_per.data,s.(snm).e_qlwi_per.date]=grp_ts(s.(snm).erai_dlrf_adj.data,s.(snm).erai_dlrf_adj.date,perfun,cumfun,minN);

  [s.(snm).e_hs_per.data,s.(snm).e_hs_per.date]=grp_ts(s.(snm).erai_sigwavehgt_adj.data,s.(snm).erai_sigwavehgt_adj.date,perfun,cumfun,minN);

  % Calculate one-day differences with gaps removed
  s.(snm).ts_per_dif.date=s.(snm).ts_per.date(2:end); s.(snm).ts_per_dif.data=diff(s.(snm).ts_per.data);
  s.(snm).ta_per_dif.date=s.(snm).ta_per.date(2:end); s.(snm).ta_per_dif.data=diff(s.(snm).ta_per.data);
  s.(snm).e_ta_per_dif.date=s.(snm).e_ta_per.date(2:end); s.(snm).e_ta_per_dif.data=diff(s.(snm).e_ta_per.data);
  s.(snm).e_qa_per_dif.date=s.(snm).e_qa_per.date(2:end); s.(snm).e_qa_per_dif.data=diff(s.(snm).e_qa_per.data);
  s.(snm).sst_per_dif.date=s.(snm).sst_per.date(2:end); s.(snm).sst_per_dif.data=diff(s.(snm).sst_per.data);
  % REMOVE GAPS...
  s.(snm) = filter_gaps(s.(snm),'ts_per','ts_per_dif',maxgaps);
  s.(snm) = filter_gaps(s.(snm),'ta_per','ta_per_dif',maxgaps);
  s.(snm) = filter_gaps(s.(snm),'e_ta_per','e_ta_per_dif',maxgaps);
  s.(snm) = filter_gaps(s.(snm),'e_qa_per','e_qa_per_dif',maxgaps);
  s.(snm) = filter_gaps(s.(snm),'sst_per','sst_per_dif',maxgaps);


  %% REMOVE ANOMALY???
  % s.(snm).U_per = anomalize_daily_mean_ts(s.(snm).U_per);
  % s.(snm).V_per = anomalize_daily_mean_ts(s.(snm).V_per);
  % s.(snm).W_per = anomalize_daily_mean_ts(s.(snm).W_per);
  % s.(snm).W2_per = anomalize_daily_mean_ts(s.(snm).W2_per);
  % s.(snm).e_hs_per = anomalize_daily_mean_ts(s.(snm).e_hs_per);
  % s.(snm).sst_xs_per = anomalize_daily_mean_ts(s.(snm).sst_xs_per);
  % s.(snm).sst_ls_per = anomalize_daily_mean_ts(s.(snm).sst_ls_per);
  % try,
  %   s.(snm).sst_l_per = anomalize_daily_mean_ts(s.(snm).sst_l_per);
  %   s.(snm).sst_dl_per = anomalize_daily_mean_ts(s.(snm).sst_dl_per);
  % catch,
  % end;
  % s.(snm).ts_per_dif = anomalize_daily_mean_ts(s.(snm).ts_per_dif);
  % s.(snm).ta_per_dif = anomalize_daily_mean_ts(s.(snm).ta_per_dif);
  % s.(snm).e_ta_per_dif = anomalize_daily_mean_ts(s.(snm).e_ta_per_dif);
  % s.(snm).e_qa_per_dif = anomalize_daily_mean_ts(s.(snm).e_qa_per_dif);
  % s.(snm).sst_per_dif = anomalize_daily_mean_ts(s.(snm).sst_per_dif);

  if ( ismember(snm,{'lonf1','smkf1'}) )
    s.(snm) = station_dewp_to_relhumid(s.(snm),'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
    s.(snm) = station_relhumid_to_spechumid(s.(snm),'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
    [s.(snm).qa_per.data,s.(snm).qa_per.date]=grp_ts(s.(snm).ndbc_spechumid.data,s.(snm).ndbc_spechumid.date,perfun,cumfun,minN);
    s.(snm).qa_per_dif.date=s.(snm).qa_per.date(2:end); s.(snm).qa_per_dif.data=diff(s.(snm).qa_per.data);
    s.(snm) = filter_gaps(s.(snm),'qa_per','qa_per_dif',maxgaps);
  end;


  if ( doCovs )

    [c,r,p,ix]=cov_ts(s.(snm).ta_per_dif,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper Tair vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    sacs = sacs + c;
    sars = sars + r;

    [c,r,p,ix]=cov_ts(s.(snm).e_ta_per_dif,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' ERAI Dper Tair vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    esacs = esacs + c;
    esars = esars + r;


    if ( ismember(snm,{'lonf1','smkf1'}) )
      [c,r,p,ix]=cov_ts(s.(snm).qa_per_dif,s.(snm).ts_per_dif,seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' Dper qa vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    end;

    [c,r,p,ix]=cov_ts(s.(snm).e_qa_per_dif,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' ERAI Dper qa vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);


    [c,r,p,ix]=cov_ts(s.(snm).U_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Per U vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).V_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Per V vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).W_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Per Wind vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    swcs = swcs + c;
    swrs = swrs + r;

    [c,r,p,ix]=cov_ts(s.(snm).e_W_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' ERAI per Wind vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    eswcs = eswcs + c;
    eswrs = eswrs + r;

    [c,r,p,ix]=cov_ts(s.(snm).W2_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Per Wind2 vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    sw2cs = sw2cs + c;
    sw2rs = sw2rs + r;

    [c,r,p,ix]=cov_ts(s.(snm).e_W2_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' ERAI per Wind2 vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    esw2cs = esw2cs + c;
    esw2rs = esw2rs + r;


    if (0)
    end;

    [c,r,p,ix]=cov_ts(s.(snm).ndbc_air_t,s.(snm).ndbc_sea_t,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Tair vs. Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    % sacs = sacs + c;
    % sars = sars + r;

    [c,r,p,ix]=cov_ts(s.(snm).ndbc_wind1_speed,s.(snm).ndbc_sea_t,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Wind vs. Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    % swcs = swcs + c;
    % swrs = swrs + r;

    [c,r,p,ix]=cov_ts(s.(snm).ndbc_wind1_speed,s.(snm).ndbc_sea_air_diff,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Wind vs. Tsea-Tair: R=' num2str(r)]);
    cs = cs + c;
    rs = rs + r;

    ug = s.(snm).ndbc_air_t;
    ug.date(ug.data<1 | ~isfinite(ug.data)) = [];
    ug.data(ug.data<1 | ~isfinite(ug.data)) = [];
    ug.data = (6./ug.data);
    [c,r,p,ix]=cov_ts(s.(snm).ndbc_wind2,ug,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' U^2 + V^2 vs. Ug^2 ~ 6/Ta: C=' num2str(c)]);
    % cs = cs + c;


    if ( ismember(snm,{'lonf1','smkf1'}) )
      s.(snm).expcdew = s.(snm).ndbc_dew_t;
      s.(snm).expcdew.data = exp(0.07*s.(snm).expcdew.data);
      s.(snm).expmcta = s.(snm).ndbc_air_t;
      s.(snm).expmcta.data = exp(-0.07*s.(snm).expmcta.data);
      % [c,r,p,ix]=cov_ts(s.(snm).ndbc_sea_air_diff,s.(snm).ndbc_wind1_speed,seasfilt);
      [c,r,p,ix]=cov_ts(s.(snm).expcdew,s.(snm).ndbc_wind1_speed,seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' exp(C*Tdew) vs. exp(-C*Ta): R=' num2str(r) ' (',num2str(p),')']);
      windrs = windrs + r;

      expta = s.(snm).ndbc_air_t;
      expta.data = exp(0.06*expta.data);
      [c,r,p,ix]=cov_ts(s.(snm).ndbc_relhumid,expta,seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' RH vs. exp(0.06*Ta): R=' num2str(r) ' (',num2str(p),')']);
      rhrs = rhrs + r;
    end;

    [c,r,p,ix]=cov_ts(s.(snm).e_qswi_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Adj ERAI QSWI vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).e_qlwi_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Adj ERAI QLWI vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).e_hs_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Adj ERAI Hs vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).sst_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' AVHRR SST vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).sst_per_dif,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).sst_xs_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST XS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    [c,r,p,ix]=cov_ts(s.(snm).sst_ls_per,s.(snm).ts_per_dif,seasfilt);
    disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST LS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

    try,
      [c,r,p,ix]=cov_ts(s.(snm).sst_l_per,s.(snm).ts_per_dif,seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST L vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);

      [c,r,p,ix]=cov_ts(s.(snm).sst_dl_per,s.(snm).ts_per_dif,seasfilt);
      disp([s.(snm).station_name,' ',char(seasfilt),' Dper AVHRR SST DL vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c),' p=',num2str(p),' N=',num2str(numel(ix))]);
    catch,
    end;

  end; %if ( doCovs )

end;
sacs=sacs/6,
sars=sars/6,
esacs=esacs/6,
esars=esars/6,
swcs=swcs/6,
swrs=swrs/6,
eswcs=eswcs/6,
eswrs=eswrs/6,
sw2cs=sw2cs/6,
sw2rs=sw2rs/6,
esw2cs=esw2cs/6,
esw2rs=esw2rs/6,
cs=cs/6,
rs=rs/6,
windrs=windrs/2,
rhrs=rhrs/2,

clear c r six aix ug expcd expmct expta

timenow;

more on;

diary off;

end; end; end;









  % %stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld,tufld,tvfld);
  % stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
%%%%%%%%%%???DEBUG
warning('EXPERIMENT WITH REDUCED DEPTH');
  stn = station_mean_tide_height(stn,mhfld,stn.depth,hfld);
%%%%%%%%%%???DEBUG



   case 'MLRF1',
    if ( use_old_erai_only_options )
      stn.opts.Ppen = 0.473;
      cbds = {8e-4};
      doWarms = [false];
      kds = { ...
          [0.035,0.190, 55] ...
          % [0.050,0.350, 91] ... %Best for new options
          % [0.045,0.375, 45] ... %ERAI met, air_t, old Ppen
            };
      advfacs = { ...
          1 ...
                };
      kths = { ...
          0 ...
             };
      % Essentially TURN OFF outlier removal
      qlh_err_cutoff = 100000;
      bq0_err_cutoff = 100000;
      hc_err_cutoff = 100.00;
    else
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          % [0.075,0.275,110] ...     % Best for WW3 1999-2011?
          [0.050,0.350, 91] ... %Best so far
          [0.050,0.200,288] ... % Best match w/1m post-cleaning measurements
          [0.350,0.700,288] ... % Best(?) for mean depth 3.55m
          [0.350,0.700, 91] ... % Best(??) for mean depth 3.55m
            };
%%%%%%%%%%???DEBUG
warning('EXPERIMENT WITH REDUCED DEPTH');
%%%%%%%%%%???DEBUG
      advfacs = { ...
          [0.0,1.0, 45] ... %Best so far
                };
      kths = { ...
          [0,20, 45] ... %Best so far
             };
    end;









%%%%%%%%%%???DEBUG
stn.depth=1;
%%%%%%%%%%???DEBUG





    set(gca,'XAxis','top');
    set(gca,'XAxis','top');



    % titlename([upper(stn.station_name),' daily BIC clean btm. insol. [MJ/m^2/day]']);





  %%%
  %% Finally estimate tidal time series based on tide-cycle mean bathymetry

  % Assume tide height already includes our estimate of site depth: subtract
  % from mean before adding, and also add the error in bathymetry estimate of
  % site depth relative to ours, to compensate for overall error in mean.
  % site_bathy_depth = -stn.(bfld).field(ix,jx);
  if ( isfield(stn,hfld) )
    % Nothing to do
  elseif ( isnumeric(hfld) )
    val = hfld;
    hfld = 'fixed_tide_depth';
    stn.(hfld).date = now;
    stn.(hfld).data = val;
  else
    error('HFLD must either be a field name or numeric value');
  end;

  % if ( isscalar(bfld) && isnumeric(bfld) )
  %   stn.(mhfld).date = stn.(hfld).date;
  %   stn.(mhfld).data = stn.(hfld).data + mean_bathy_depth;
  % else
  %   site_bathy_depth = nanmean(stn.(hfld).data);
  %   stn.(mhfld).date = stn.(hfld).date;
  %   stn.(mhfld).data = stn.(hfld).data + mean_bathy_depth - site_bathy_depth;
  % end;
  site_bathy_depth = nanmean(stn.(hfld).data);
  stn.(mhfld).date = stn.(hfld).date;
  stn.(mhfld).data = stn.(hfld).data + mean_bathy_depth - site_bathy_depth;







    % Experiment - compare with underwater PAR measurements at ~1m depth
    stn = station_mean_tide_height(stn,mhfld,-2.55,hfld);



%%%%??? DEBUG
  stn.(srfld) = ts_op(stn.erai_dsrf,stn.(usrfld),'+');




    stn.(diagfld).tau_Ppen.date = dts;
    stn.(diagfld).tau_Ppen.data = Ppen .* tau;
    stn.(diagfld).tau_Ppen.data(badix) = 0;



  stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,stn.opts,asrdiagfld);
  bdsrfld = ['b_' RAPFX '_' ISPFX '_dsrf'];
  stn.(bdsrfld) = ts_op(stn.(srfld),

  stn.([bdsrfld,'_daily_dose']) = par_dose(stn.(bdsrfld));




if (0)
end;

ylm=def_ylm;
lineWidth=def_lineWidth;

stn = get_station_from_station_name('looe1');
stn = get_looe1_microcat(stn);
ylm=def_ylm_ts;
[Pxx,Wd,fh,lhs]=plot_spec(stn,'microcat_seatemp',[],[],xlm,ylm,printFmt,doNorm);
set(lhs,'LineWidth',lineWidth);

stn = get_looe1_adcp(stn);
for cfld={'adcp_seatemp','adcp_speed','adcp_x','adcp_l'};
    fld=cfld{:};
    switch (lower(fld)),
     case 'adcp_seatemp',
      ylm=def_ylm_ts;
    end;
    [Pxx,Wd,fh,lhs]=plot_spec(stn,fld,[],[],xlm,ylm,printFmt,doNorm);
    set(lhs,'LineWidth',lineWidth);
end;
%%%%???DEBUG
close all;
stn=[]; clear stn Pxx Wd fh lhs;






for kth = {'ndbc_wind1_speed'};




optimize...m:
  if ( adjust_waves )
    switch (WAVEPFX),
     case 'erai',	stn = adjust_erai_station_waves(stn);
     case 'ww3',	stn = adjust_ww3_station_waves(stn);
     case 'ndbc',	
      stn.(wpfld).data = ((stn.(wpfld).data.*1.00) + 0.0);
      stn.(whfld).data = ((stn.(whfld).data.*0.25) + 0.3);
      stn.(wdfld).data = ((stn.(wdfld).data.*1.00) + 0.0);
     otherwise,		error('Unknown wave source "%s"',WAVEPFX);
    end;
  end;


fh=1; figure(fh);
fh=fh+1; figure(fh);





%%%%???DEBUG
figspath = get_thesis_path('../figs/save');




disp('Current=1% of LP wind speed (NOT YET IMPLEMENTED)');
disp('Current=2% of LP wind speed (NOT YET IMPLEMENTED)');


for kd = {[0.05,0.35,137],[0.05,0.30, 91],0.2};
  kd{:},




function stn = optimize_station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,comparisonsfld)
%function stn = optimize_station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,comparisonsfld)
%
% Optimize all parameters of reef ocean heat budget by graphing results for the caller
%
% Last Saved Time-stamp: <Sat 2012-08-04 10:34:50  lew.gramer>



function opts = station_heat_budget_options(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld)
%function opts = station_heat_budget_options(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld)
%
% Set initial site-specific options for ocean heat budget. Usually called
% from, e.g., STATION_HEAT_BUDGET (see).
%
% Last Saved Time-stamp: <Tue 2012-07-31 14:37:14  lew.gramer>



function stn = station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%function stn = station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%
% Calculate errors for each term in reef ocean heat budget. Uses TOGA/COARE
% error algorithm for turbulent fluxes, and published relative and absolute
% errors, confirmed by local in situ regressions where ever available, to
% calculate error budget for all other terms. See STATION_HEAT_BUDGET or
% OPTIM_Q0 for args. Adds new fields to STN with '_err' appended to name.
%
% Last Saved Time-stamp: <Tue 2012-07-31 16:08:02  lew.gramer>



function dump_station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%function dump_station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%
% Output a report with DISP of the median and "variance" (max-min/6) of error
% in terms of heat budget as calculated by, e.g., STATION_HEAT_BUDGET_ERRORS.
%
% Last Saved Time-stamp: <Sun 2012-04-15 17:07:34  lew.gramer>


function annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr,mos,begyr,endyr,comparisonsfld)
%function annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr,mos,begyr,endyr,comparisonsfld)
%
% Plot yearly cumulative comparison (one subplot per year of sea temperature)
% for various forms (intermediate results) of the reef ocean heat budget.
% First eight args (STN to AFLD) are documented in STATION_HEAT_BUDGET (v.).
% COMMENTSTR (DEFAULT: STN.commentstr or '') is a string appended to plot
% title; MOS, BEGYR, ENDYR constrain which months of which years are plotted;
% finally STN.(COMPARISONSFLD) plots a sea temperature other than STN.(SFLD)
% for comparison with the heat budget results. NOTE: Resets on any gap>=7d.
%
% Last Saved Time-stamp: <Sun 2012-04-15 17:07:52  lew.gramer>


function res = compare_flux_climatologies(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr)
%function res = compare_flux_climatologies(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr)
%
% Plot annual or interannual heat-budget term climatologies from Gramer and
% Mariano (2012), OAFlux (Yu and Weller 2007) and ISCCP (Zhang et al. 2004).
%
% Last Saved Time-stamp: <Tue 2012-04-03 13:15:58  Lew.Gramer>


function plot_budget_years(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr,mos,begyr,endyr,comparisonsfld)
%function plot_budget_years(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr,mos,begyr,endyr,comparisonsfld)
%
% Plot yearly cumulative comparison (one subplot per year of sea temperature)
% for various forms (intermediate results) of the reef ocean heat budget.
% First eight args (STN to AFLD) are documented in STATION_HEAT_BUDGET (v.).
% COMMENTSTR (DEFAULT: STN.commentstr or '') is a string appended to plot
% title; MOS, BEGYR, ENDYR constrain which months of which years are plotted;
% finally STN.(COMPARISONSFLD) plots a sea temperature other than STN.(SFLD)
% to compare with heat budget result. NOTE: Resets comparison on any gap>=7d.
%
% Lew.Gramer@noaa.gov: ADAPTED from ANNSUBS.m (v.) on 2012 Mar 27
%
% Last Saved Time-stamp: <Tue 2012-07-31 16:47:20  lew.gramer>


function stn = ms1_clim(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr,tempfld,budgetfld,mos,begyr,endyr)
%function stn = ms1_clim(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr,tempfld,budgetfld,mos,begyr,endyr)
%
% Plot comparison of ocean heat budget climatology (expressed as an implied
% sea teperature time series) with climatology of actual sea temperature.
%
% Last Saved Time-stamp: <Wed 2011-11-23 21:04:51  Lew.Gramer>


function compare_heat_budget_cumsums(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%function compare_heat_budget_cumsums(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld)
%
% Plot comparison of simple cumulative sums of daily climatology vs. heat budget
% estimates. This compares both annual amplitude and interannual variability.
%
% USES: dsffld,climq0fld,sq0fld,bq0fld,qtAdvffld,bdTffld,hcdTdtf
%
% Last Saved Time-stamp: <Sun 2012-01-15 14:17:17  Lew.Gramer>


function res = compare_flux_climatologies(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr)
%function res = compare_flux_climatologies(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr)
%
% Plot annual or interannual heat-budget term climatologies from Gramer and
% Mariano (2012), OAFlux (Yu and Weller 2007) and ISCCP (Zhang et al. 2004).
%
% Last Saved Time-stamp: <Tue 2012-04-03 13:15:58  Lew.Gramer>


function res = regress_temp_budget_terms(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,doPlot,saveFile)
%function res = regress_temp_budget_terms(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,doPlot,saveFile)
%
% Regress heat budget terms vs. variability in DAILY MEAN sea temperature.
% If DOPLOT (DEFAULT: false), plots ROBUSTFIT linear regressions; displays a
% report using DISP of simple R^2, RMSE, Scatter Index, and of R^2, RMSE, SI,
% bias, and slope from regression; if SAVEFILE (DEFAULT: false), also output
% to the CSV file FULLFILE(DATAPATH,[STNM,'-',MNM,'-',hcdTdt,'.csv']).
%
% See STATION_HEAT_BUDGET for the meaning of the other arguments.
%
% Last Saved Time-stamp: <Tue 2012-07-31 17:01:01  lew.gramer>




  %
  if ( exist('substitute_field_names','var') && ...
       iscell(substitute_field_names) && ...
       mod(numel(substitute_field_names),2) == 0 )
    for ix=1:2:numel(substitute_field_names)
      assignin('caller',substitute_field_names{ix},substitute_field_names{ix+1});
    end;
  end;






  if ( isfield(stn_or_stnm,'station_name') )
    stnm = stn_or_stnm.station_name;
  elseif ( ischar(stn_or_stnm) )
    stnm = stn_or_stnm;
  else
    error('STN_OR_STNM must be name string or struct with .station_name field!');
  end;

  station_heat_budget_field_names;

  opts = [];

  switch (stnm),

   case {'fwyf1'},
    opts.default_salinity = get_opt(opts,'default_salinity',35.5);
    % Results of running OPTIM_Q0 on ERAI met data
    % *NOTE*: With *in situ* met data, Warm Layer and different/higher Kd are better!
    opts.kd = get_opt(opts,'kd',[0.10,0.30,0]);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    % opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',8.0e-4);
    % opts.hc_scaling = get_opt(opts,'hc_scaling','US');






  scatter_fit_ts_seasons(stn.gom_hycom_seatemp,stn.ndbc_sea_t_7_d_avg,[],[],'GOM HYCOM',[STNM,' \mu_7_dT_s'],[],[],true);
  subplots_set('xlim',[18,32],'ylim',[18,32]);
  print('-dtiff',fullfile(figspath,[stnm,'-gom_hycom_seatemp-scatter-ndbc_sea_t_7_d_avg.tif']));

  scatter_fit_ts_seasons(stn.fkeys_hycom_seatemp,stn.ndbc_sea_t_7_d_avg,[],[],'FKEYS HYCOM',[STNM,' \mu_7_dT_s'],[],[],true);
  subplots_set('xlim',[18,32],'ylim',[18,32]);
  print('-dtiff',fullfile(figspath,[stnm,'-fkeys_hycom_seatemp-scatter-ndbc_sea_t_7_d_avg.tif']));





rawaflds = {'ndbc_air_t','ndbc_spechumid','ndbc_wind1_speed','ndbc_wind1_lshore','ndbc_wind1_xshore','hourly_avhrr_weekly_sst','hourly_avhrr_weekly_sst_lshore','hourly_avhrr_weekly_sst_xshore'};
%rawaflds = {'hourly_avhrr_weekly_sst','hourly_avhrr_weekly_sst_lshore','hourly_avhrr_weekly_sst_xshore'};
%rawaflds = {'ndbc_air_t'};

stnms = {'fwyf1','mlrf1','lonf1','smkf1','sanf1'};
%stnms = {'smkf1'};
%stnms = {'mlrf1'};





FROM phd_proposal_fig_7.m:
1;

if ( ~exist('stn','var') )
    if ( ~exist('stnm','var') )
        %stnm='fwyf1';
        %stnm='mlrf1';
        %stnm='lonf1';
        stnm='dryf1';
    end;
    stn = get_station_from_station_name(stnm);
    stn = load_all_ndbc_data(stn);
    stn = qa_ts(stn,'ndbc_sea_t');
    stn = qa_ts(stn,'ndbc_air_t');
    stn = qa_ts(stn,'ndbc_wind1_speed');
end;

%stn = verify_variable(stn,{'ndbc_wind1_u_7_d_var','ndbc_wind1_v_7_d_var'});
%    'ndbc_wind1_u_3_d_avg_7_d_var_0_d_asof_sum_ndbc_wind1_v_3_d_avg_7_d_var',...

flds = {...
    'ndbc_sea_t_qc',...
    'ndbc_sea_t_qc_1_d_var_3_d_avg',...
    'ndbc_air_t_qc',...
    'ndbc_air_t_qc_1_d_var_3_d_avg',...
    'ndbc_wind1_u_qc_3_d_avg',...
    'ndbc_wind1_v_qc_3_d_avg',...
    'ndbc_wind1_speed_qc_3_d_avg',...
    'ndbc_wind1_u_qc_7_d_var_7_d_avg_sum_ndbc_wind1_v_qc_7_d_var',...
    };
ylbls = {...
    'T_s',...
    '\mu_3_d(\sigma_1_dT_s)',...
    'T_a',...
    '\mu_3_d(\sigma_1_dT_a)',...
    'W_U',...
    'W_V',...
    '\mu_3_d(W)',...
    '\sigma_7_d(W_U)+\sigma_7_d(W_V)',...
    };

stn = verify_variable(stn,flds);

multiplot_station(stn,flds,[],[],ylbls,...
    datenum([1987,2008],[6,6],[1,1]),...
    {[0,40],[0,2],[0,40],[0,4],[-50,+50],[-50,+50],[0,20],[0,20]});
datetick3('x',10,'keeplimits');





  disp('SCALED UP BY 1km');
  s.(snm).sst_xs_per.data=s.(snm).sst_xs_per.data.*1e3;
  s.(snm).sst_ls_per.data=s.(snm).sst_ls_per.data.*1e3;
  s.(snm).sst_l_per.data=s.(snm).sst_l_per.data.*1e3*1e3;
  s.(snm).sst_dl_per.data=s.(snm).sst_dl_per.data.*1e3*1e3;





  [c,r,p]=cov_ts(s.(snm).hourly_avhrr_weekly_sst_xshore,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' AVHRR SST XS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c)]);

  [c,r,p]=cov_ts(s.(snm).hourly_avhrr_weekly_sst_lshore,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' AVHRR SST LS vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c)]);


  [c,r,p]=cov_ts(s.(snm).hourly_avhrr_weekly_sst_l,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' AVHRR SST L vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c)]);
  % sacs = sacs + c;
  % sars = sars + r;

  [c,r,p]=cov_ts(s.(snm).hourly_avhrr_weekly_sst_dl,s.(snm).ts_per_dif,seasfilt);
  disp([s.(snm).station_name,' ',char(seasfilt),' AVHRR SST DL vs. Dper Tsea: R2=',num2str(r.^2),' C=',num2str(c)]);
  % sacs = sacs + c;
  % sars = sars + r;







  if ( isfield(stn,'raw_erai_dsrf') )
  end;
  if ( isfield(stn,'erai_dsrf') )
  end;




for perfun=[@floor,@get_yearmonth]; for seasfilt=[@finite_ts,@ts_boreal_warm,@ts_boreal_cool];
for filterBadDates = [false,true];

disp(char(perfun));
disp(char(seasfilt));
if ( filterBadDates )
  disp('WILL FILTER "BAD DATES"');
else
  disp('Dates unfiltered');
end;

end; end; end;

FOOBAR







perfun = @floor;                cumfun = @nanmean;      minN = 24;      maxgaps=(25/24);
%perfun = @get_yearweek;         cumfun = @nanmean;      minN = 24*4.5;  maxgaps=8;
%perfun = @get_yearmonth;        cumfun = @nanmean;      minN = 24*20;   maxgaps=32;

seasfilt = [];
%seasfilt = @ts_boreal_warm;
%seasfilt = @ts_boreal_cool;

%filterBadDates = false;
filterBadDates = true;

if ( filterBadDates ); filtstr = 'filtered'; else filtstr = 'unfiltered'; end;

diary off;
fname = fullfile(get_thesis_path('../data'),[mfilename,'-',char(perfun),'-',char(seasfilt),'-',filtstr,'-results.log']);
diary(fname);







if ~exist('s.fwyf1','var'); s.fwyf1 = get_station_from_station_name('fwyf1'); s.fwyf1 = load_all_ndbc_data(s.fwyf1); end;
s.fwyf1 = verify_variable(s.fwyf1,'ndbc_wind1_u');
s.fwyf1 = verify_variable(s.fwyf1,'ndbc_wind1_v');
[six,aix]=intersect_dates(s.fwyf1.ndbc_sea_t.date,s.fwyf1.ndbc_air_t.date);
s.fwyf1.ndbc_sea_air_diff.date=s.fwyf1.ndbc_sea_t.date(six);
s.fwyf1.ndbc_sea_air_diff.data=s.fwyf1.ndbc_sea_t.data(six)-s.fwyf1.ndbc_air_t.data(aix); 

if ~exist('s.mlrf1','var'); s.mlrf1 = get_station_from_station_name('mlrf1'); s.mlrf1 = load_all_ndbc_data(s.mlrf1); end;
s.mlrf1 = verify_variable(s.mlrf1,'ndbc_wind1_u');
s.mlrf1 = verify_variable(s.mlrf1,'ndbc_wind1_v');
[six,aix]=intersect_dates(s.mlrf1.ndbc_sea_t.date,s.mlrf1.ndbc_air_t.date);
s.mlrf1.ndbc_sea_air_diff.date=s.mlrf1.ndbc_sea_t.date(six);
s.mlrf1.ndbc_sea_air_diff.data=s.mlrf1.ndbc_sea_t.data(six)-s.mlrf1.ndbc_air_t.data(aix); 

if ~exist('s.lonf1','var'); s.lonf1 = get_station_from_station_name('lonf1'); s.lonf1 = load_all_ndbc_data(s.lonf1); end;
s.lonf1 = verify_variable(s.lonf1,'ndbc_wind1_u');
s.lonf1 = verify_variable(s.lonf1,'ndbc_wind1_v');
s.lonf1 = station_dewp_to_relhumid(s.lonf1,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
s.lonf1 = station_relhumid_to_spechumid(s.lonf1,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
[six,aix]=intersect_dates(s.lonf1.ndbc_sea_t.date,s.lonf1.ndbc_air_t.date);
s.lonf1.ndbc_sea_air_diff.date=s.lonf1.ndbc_sea_t.date(six);
s.lonf1.ndbc_sea_air_diff.data=s.lonf1.ndbc_sea_t.data(six)-s.lonf1.ndbc_air_t.data(aix); 

if ~exist('s.smkf1','var'); s.smkf1 = get_station_from_station_name('smkf1'); s.smkf1 = load_all_ndbc_data(s.smkf1); end;
s.smkf1 = verify_variable(s.smkf1,'ndbc_wind1_u');
s.smkf1 = verify_variable(s.smkf1,'ndbc_wind1_v');
s.smkf1 = station_dewp_to_relhumid(s.smkf1,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
s.smkf1 = station_relhumid_to_spechumid(s.smkf1,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
[six,aix]=intersect_dates(s.smkf1.ndbc_sea_t.date,s.smkf1.ndbc_air_t.date);
s.smkf1.ndbc_sea_air_diff.date=s.smkf1.ndbc_sea_t.date(six);
s.smkf1.ndbc_sea_air_diff.data=s.smkf1.ndbc_sea_t.data(six)-s.smkf1.ndbc_air_t.data(aix); 

if ~exist('s.sanf1','var'); s.sanf1 = get_station_from_station_name('sanf1'); s.sanf1 = load_all_ndbc_data(s.sanf1); end;
s.sanf1 = verify_variable(s.sanf1,'ndbc_wind1_u');
s.sanf1 = verify_variable(s.sanf1,'ndbc_wind1_v');
[six,aix]=intersect_dates(s.sanf1.ndbc_sea_t.date,s.sanf1.ndbc_air_t.date);
s.sanf1.ndbc_sea_air_diff.date=s.sanf1.ndbc_sea_t.date(six);
s.sanf1.ndbc_sea_air_diff.data=s.sanf1.ndbc_sea_t.data(six)-s.sanf1.ndbc_air_t.data(aix); 

if ~exist('s.dryf1','var'); s.dryf1 = get_station_from_station_name('dryf1'); s.dryf1 = load_all_ndbc_data(s.dryf1); end;
s.dryf1 = verify_variable(s.dryf1,'ndbc_wind1_u');
s.dryf1 = verify_variable(s.dryf1,'ndbc_wind1_v');
[six,aix]=intersect_dates(s.dryf1.ndbc_sea_t.date,s.dryf1.ndbc_air_t.date);
s.dryf1.ndbc_sea_air_diff.date=s.dryf1.ndbc_sea_t.date(six);
s.dryf1.ndbc_sea_air_diff.data=s.dryf1.ndbc_sea_t.data(six)-s.dryf1.ndbc_air_t.data(aix); 






          fmg;
          climsq = squeeze(stn.optim.climsq(:,doWarmix,advix,kthix));
          climsq_err = stn.optim.climsq_err(end,doWarmix,advix,kthix);
          climsq_minus_err = ts_op(climsq(end),climsq_err,'-');
          climsq_plus_err = ts_op(climsq(end),climsq_err,'+');
          % climsq_minus_err = ts_op(climsq(end),ts_op(climsq_err,24,'*'),'-');
          % climsq_plus_err = ts_op(climsq(end),ts_op(climsq_err,24,'*'),'+');
          % % climsq_minus_err.date = climsq_err.date;
          % % % climsq_minus_err.data = climsq(end).data - cumsum(climsq_err.data);
          % % climsq_minus_err.data = climsq(end).data - (climsq_err.data.*24);
          % % climsq_plus_err.date = climsq_err.date;
          % % % climsq_plus_err.data = climsq(end).data + cumsum(climsq_err.data);
          % % % plot_ts(stn.optim.climt,climsq,climsq_minus_err,'k^',climsq_plus_err,'kv');
          % % climsq_plus_err.data = climsq(end).data + (climsq_err.data.*24);

          % lhs=plot_ts(stn.optim.climt,'LineWidth',3,climsq,climsq_minus_err,'k:',climsq_plus_err,'k:');

          climssqt = squeeze(stn.optim.climssqt(end,doWarmix,advix,kthix));
          climsdt = squeeze(stn.optim.climsdt(end,doWarmix,advix,kthix));

          lhs = plot_ts(stn.optim.climt,'LineWidth',3,climsq);
          lhs(end+1) = plot_ts(climsq_minus_err,'k:');
          plot_ts(climsq_plus_err,'k:');
          lhs(end+1) = plot_ts(climssqt,'r.');
          lhs(end+1) = plot_ts(climsdt,'mo');

          datetick3('x',3);
          titlename([STNM ' Daily Clim: ' doWarmStr ' ' stn.optim.cbdstrs{default_cbdix}]);
          % % legend({'T_s',stn.optim.kdstrs{:},[stn.optim.kdstrs{end},' \pm error']});
          kdstrs = strcat( stn.optim.kdstrs',{' Errs '},...
                           cellstr(num2str(stn.optim.dayerror(:,doWarmix,advix,kthix),'%.1f')),{', '},...
                           cellstr(num2str(stn.optim.climerror(:,doWarmix,advix,kthix),'%.1f')) );
          % legend(lhs,{'T_s',kdstrs{:},[stn.optim.kdstrs{end},' \pm error']}, 'Location','South');
          legend(lhs,...
                 {'T_s',kdstrs{:},[stn.optim.kdstrs{end},' \pm error'],'Q_0','Q_0(\gamma)+Q_b',}, ...
                 'Location','South');
          xlim(stn.optim.climt.date([1 end]));
          ylim([minmin([stn.optim.climsq.data]),maxmax([stn.optim.climsq.data])]);
          ylim([16,34]);
          ylabel('^oC');
          % appendtitlename([' (' strrep(QEPFX,'_','\_') ' Adv=' stn.optim.advfacstrs{advix} ' K_\theta=' stn.optim.kthstrs{kthix} ')' stn.commentstr]);
          appendtitlename([' (' strrep(hcdTdt,'_','\_') ' Adv=' stn.optim.advfacstrs{advix} ' K_\theta=' stn.optim.kthstrs{kthix} ')' stn.commentstr]);
          appendtitlename([' (' num2str(begyr) '-' num2str(endyr) ')']);








          plot_ts(stn.optim.climt,'LineWidth',3,climsq,climsq_minus_err,'k:',climsq_plus_err,'k:',climssqt,'r.',climsdt,'mo');




            if ( doPlot )


              %% Calculate and plot time series and errors from this option
              clear t q qe
              [tix,qix,qeix] = intersect_all_dates([],stn.(sfld).date,stn.(hcdTdt).date,stn.([hcdTdt,'_err']).date);
              t.date = stn.(sfld).date(tix);			t.data = stn.(sfld).data(tix);
              q.date = stn.(hcdTdt).date(qix);			q.data = stn.(hcdTdt).data(qix);
              qe.date = stn.([hcdTdt,'_err']).date(qeix);	qe.data = stn.([hcdTdt,'_err']).data(qeix);







            end; %if ( ~isempty(baddays) )

            %DEBUG:
            disp(numel(stn.(hcdTdt).data));

            stn = station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld);
            %DEBUG:
            disp('POST-QC'); dump_station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld);


            if ( doPlot )








  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Optimization parameters

  % Number of points to use in finite-difference templates for gradients
  npts = 5;

  % "Fudge factor" for Latent Heat Flux
  qlh_adj = 1.0;

  % Apply TOGA-COARE 3.0a Warm-layer adjustment for sea temp. sensor depth?
  doWarms = [false, true];
  % Attentuation coefficient/climatology for "penetrative" (PAR+NUV) insolation
  kds = {...
      0.1,[0.05,0.15,45],[0.05,0.15,137],[0.05,0.15,228],[0.05,0.15,320],...
      0.2,[0.10,0.30,45],[0.10,0.30,137],[0.10,0.30,228],[0.10,0.30,320],...
      0.3,[0.20,0.40,45],[0.20,0.40,137],[0.20,0.40,228],[0.20,0.40,320],...
        };

  % Bulk benthic heat coefficient for flow over sea-floor
  % % %cbds = {0,3.8e-5,2.9e-4,3.8e-4,8.0e-4,16.0e-4,24.0e-4};
  % % Cbd~CD^2, based on estimate CD~0.017 in Davis and Monismith (2011)
  % cbds = {2.9e-4};
  % Davis and Monismith (2011) CD times Cbh, where Cbh~1e-3
  cbds = {0.017*1e-2};

  % Heat advection fudge factor (0=no advection)
  advfacs = {...
      0,...
      1,...
            };
  % Heat diffusion fudge factor (0=no diffusion)
  kths = {...
      0 ,...
      5, [0 ,10, 45],...
      10, [0 ,20, 45],...
         };


  switch ( STNM ),

   case 'FWYF1',
    if ( use_old_erai_only_options )
      stn.opts.Ppen = 0.473;
      cbds = {8e-4};
      doWarms = [false];
      kds = { ...
          [0.035,0.220, 76] ...
            };
      advfacs = { ...
          1 ...
                };
      kths = { ...
          0 ...
             };
      % Essentially TURN OFF outlier removal
      qlh_err_cutoff = 100000;
      bq0_err_cutoff = 100000;
      hc_err_cutoff = 100.00;
    else
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          %[0.025,0.300,110] ...     % Best for no waves 1987-2011 and(?) WW3 1999-2011
          [0.025,0.300,137] ...     % Best for ERAI waves 1993-2011
            };
      advfacs = { ...
          [0.00,1.00, 91] ...
                };
      kths = { ...
          [0,10, 91] ...
             };
    end;


   case 'MLRF1',
    if ( use_old_erai_only_options )
      stn.opts.Ppen = 0.473;
      cbds = {8e-4};
      doWarms = [false];
      kds = { ...
          [0.035,0.190, 55] ...
          [0.050,0.350, 91] ... %Best for new options
          [0.045,0.375, 45] ... %ERAI met, air_t, old Ppen
            };
      advfacs = { ...
          1 ...
                };
      kths = { ...
          0 ...
             };
      % Essentially TURN OFF outlier removal
      qlh_err_cutoff = 100000;
      bq0_err_cutoff = 100000;
      hc_err_cutoff = 100.00;
    else
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          % [0.075,0.275,110] ...     % Best for WW3 1999-2011?
          [0.050,0.350, 91] ... %Best so far
            };
      advfacs = { ...
          [0.0,1.0, 45] ... %Best so far
                };
      kths = { ...
          [0,20, 45] ... %Best so far
             };
    end;


   case 'LONF1',
    if ( use_old_erai_only_options )
      stn.opts.Ppen = 0.473;
      cbds = {8e-4};
      doWarms = [false];
      kds = { ...
          [0.200,0.370,228] ...
            };
      advfacs = { ...
          1 ...
                };
      kths = { ...
          [0,5,183] ...
             };
      % Essentially TURN OFF outlier removal
      qlh_err_cutoff = 100000;
      bq0_err_cutoff = 100000;
      hc_err_cutoff = 100.00;
    else
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      %stn.opts.hc_scaling = 'SU';
      %stn.commentstr = [stn.commentstr,' (',stn.opts.hc_scaling,') '];
      doWarms = [false];
      kds = { ...
          % [0.700,1.200, 20] ... % Ended closest to balance
          [0.675,1.250, 10] ... % Best for ERAI waves adv=[0,1,91],K=[0,2,91]
            };
      advfacs = { ...
          [0.0,1.0, 91] ... % Best so far
                };
      kths = { ...
          [0, 2, 91] ... % We want SOME diffusion if we can have it
             };
    end;


   case 'SMKF1',
    if ( use_old_erai_only_options )
      stn.opts.Ppen = 0.473;
      cbds = {8e-4};
      doWarms = [false];
      kds = { ...
          [0.045,0.375, 45] ...
            };
      advfacs = { ...
          1 ...
                };
      kths = { ...
          0 ...
             };
      % Essentially TURN OFF outlier removal
      qlh_err_cutoff = 100000;
      bq0_err_cutoff = 100000;
      hc_err_cutoff = 100.00;
    else
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      %stn.opts.hc_scaling = 'SU';
      %stn.commentstr = [stn.commentstr,' (',stn.opts.hc_scaling,') '];
      doWarms = [true];
      kds = { ...
          % [0.100,0.400,110] ...    % "Best" for *all* HC scalings [US][SU]
          [0.100,0.600, 70] ...    % Best for NO GRADIENT case - Best so far?
            };
      advfacs = { ...
          [0.0,0.25, 91] ...
                };
      kths = { ...
          [0, 2, 91] ...
             };
    end;

   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    % So-so RMSE 5.5 (best was 2.1), but best *climatological* RMSE (0.8)
    kds = { ...
        [0.010,0.300, 81] ... %Best so far
          };
    advfacs = { ...
        0 ...
              };
    kths = { ...
        [0, 2, 91] ... %Best so far
           };

   case 'DRYF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
%%%%???DEBUG
    stn.opts.hc_scaling = 'SU';
    stn.commentstr = [stn.commentstr,' (',stn.opts.hc_scaling,') '];
%%%%???DEBUG
    qlh_adj = 1.1;
    doWarms = [false];
%%%%???DEBUG
    doWarms = [true];
    kds = { ...
        % [0.300,0.700,342] ...
        % [0.300,0.700,  0] ... %Best so far?
        % [0.250,0.750,  0] ...
        [0.200,0.800,  0] ...
        [0.250,0.750,342] ... %Best so far?
        [0.200,0.800,342] ...
          };
    advfacs = { ...
        [0.00,1.00, 45] ... %Best so far
              };
    kths = { ...
        [0, 2, 45] ...
           };

   case 'LOOE1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          [0.050,0.400,137] ... %BEST FOR 0, 0, all years
            };
      advfacs = { ...
          [0.00,0.50, 45] ...
                };
      kths = { ...
          0 ...                 %BEST
             };

    elseif ( ~isempty(regexp(sfld,'adcp')) || ~isempty(regexp(sfld,'ad_')) )

      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          [0.010,0.400,113] ... %BEST SO FAR
            };
      advfacs = { ...
          [0.00,0.50, 45] ...
                };
      kths = { ...
          [0,20, 45] ...
             };

    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;

   case 'HAWK1',
    error('Hawk Channel SFP station HAWK1 *not yet implemented*');

   case 'NCR31',
    error('NCORE site "3" *not yet implemented*');

   case 'NCRC1',
    error('NCORE site "C" *not yet implemented*');

   case {...
       'BNPIN',...
       'BNPMI',...
       'BNPON',...
       'BNPNN',...
       'BNPPA',...
       'TAVRK',...
       'CONSH',...
       'CONDP',...
        },
    % stn.opts.hc_scaling = 'SS';
    % stn.opts.hc_scaling = 'SU';
    % stn.opts.hc_scaling = 'UU';
    stn.commentstr = [stn.commentstr,' (',stn.opts.hc_scaling,') '];
    doWarms = [false];
    kds = { ...
        % [0.700,1.200, 20] ... % Ended closest to balance
        [0.675,1.250, 10] ... % Best for ERAI waves adv=[0,1,91],K=[0,2,91]
          };
    advfacs = { ...
        [0.0,1.0, 91] ... % Best so far
              };
    kths = { ...
        [0, 2, 91] ... % We want SOME diffusion if we can have it
           };

   otherwise,
    error('Station %s options not implemented yet!',STNM);
  end;







  if ( strcmpi(ISPFX,'erai') )
    warning('Special ERAI options');
    begyr = 1996;
    adjust_waves = false;
    adjust_reanalysis = false;
    reanalysis_shortwave = true;
    reanalysis_longwave = true;
  end;




  s1dfld = [sfld '_1_d_avg'];
  stn = verify_variable(stn,s1dfld);
  a1dfld = [afld '_1_d_avg'];
  stn = verify_variable(stn,a1dfld);
  sq01dfld = [sq0fld '_1_d_avg'];
  stn = verify_variable(stn,sq01dfld);
  sr1dfld = [srfld '_1_d_avg'];
  stn = verify_variable(stn,sr1dfld);
  asr1dfld = [asrfld '_1_d_avg'];
  stn = verify_variable(stn,asr1dfld);
  lr1dfld = [lrfld '_1_d_avg'];
  stn = verify_variable(stn,lr1dfld);
  qlh1dfld = [qlhfld '_1_d_avg'];
  stn = verify_variable(stn,qlh1dfld);
  qsh1dfld = [qshfld '_1_d_avg'];
  stn = verify_variable(stn,qsh1dfld);
  qcool1dfld = [qcoolfld '_1_d_avg'];
  stn = verify_variable(stn,qcool1dfld);
  fqudTffld1d = [fqudTffld '_1_d_avg'];
  stn = verify_variable(stn,fqudTffld1d);
  qbofld1d = [qbofld '_1_d_avg'];
  stn = verify_variable(stn,qbofld1d);
  kd2Tffld1d = [kd2Tffld '_1_d_avg'];
  stn = verify_variable(stn,kd2Tffld1d);
  hcdTdtf1d = [hcdTdtf '_1_d_avg'];
  stn = verify_variable(stn,hcdTdtf1d);
  hcdTdthcf1d = [hcdTdthcf '_1_d_avg'];






  % % Net transport ex mixing efficiency (DEFAULT: 8% mixing efficiency)
  % R = get_opt(opts,'hc_R',(1.00-0.08));



  % %DEBUG:
  % res.dTdt_hour=res.dTdtq0+res.dTdthc;
  % if (doDebug)
  %     disp([descStr,' Iterating HC']); tic,
  % end;
  % for ix=2:length(res.dTdt)
  %   % % res.dTdtq0(ix) = res.dTdtq0(ix) + res.dTdtq0(ix-1) + res.dTdthc(ix-1);
  %   % % res.dTdtx(ix) = res.dTdtx(ix) + res.dTdtx(ix-1) - res.dTdthc(ix-1);
  %   % res.dTdtq0(ix) = res.dTdtq0(ix) + res.dTdthc(ix-1);
  %   % res.dTdtx(ix) = res.dTdtx(ix) - res.dTdthc(ix-1);
  %   res.dTdtq0(ix) = res.dTdtq0(ix) - res.dTdthc(ix-1);
  %   res.dTdtx(ix) = res.dTdtx(ix) + res.dTdtx(ix-1) + res.dTdthc(ix-1);

  %   res.dTdx(ix)=(res.dTdtq0(ix)-res.dTdtx(ix))./res.dx(ix);
  %   res.dTdx(~isfinite(res.dTdx)) = 0;
  %   %%%% Should we actually recalc res.B0/uf/u at each new hour also? NO!!
  %   res.dTdthc(ix)=-R.*dt.*res.u(ix).*res.dTdx(ix);
  %   res.dTdt(ix)=res.dTdtq0(ix)+res.dTdthc(ix);
  % end;
  % %DEBUG:
  % toc,








              %[cum,tid] = grp_ts(qe.data,qe.date,'daily',@nanmean,23);
              %[cum,tid] = grp_ts(qe.data,qe.date,'daily',@nansum,23);
              [cum,tid] = grp_ts(qe.data,qe.date,'daily',@(x)(nanmean(abs(x))),23);
              %[cum,tid] = grp_ts(qe.data,qe.date,'daily',@(x)(nanmax(abs(x))),23);
              cum = 24*cum;
              %[cum,tid] = grp_ts(qe.data,qe.date,'daily',@(x)(nansum(abs(x))),23);
              stn.optim.climsq_err(kdix,doWarmix,advix,kthix).date = tid;
              stn.optim.climsq_err(kdix,doWarmix,advix,kthix).data = cum;






            baddts = [];
            % qlh_err_cutoff = 1000;
            qlh_err_cutoff = 500;
            stn = verify_variable(stn,[qlhfld,'_err_1_d_sum']);
            baddts = [baddts;...
                      stn.([qlhfld,'_err_1_d_sum']).date(abs(stn.([qlhfld,'_err_1_d_sum']).data)>qlh_err_cutoff)];
            % bq0_err_cutoff = 1500;
            bq0_err_cutoff = 1000;
            stn = verify_variable(stn,[bq0fld,'_err_1_d_sum']);
            baddts = [baddts;...
                      stn.([bq0fld,'_err_1_d_sum']).date(abs(stn.([bq0fld,'_err_1_d_sum']).data)>bq0_err_cutoff)];
            hc_err_cutoff = 10.00;
            % hc_err_cutoff = 13.00;
            % hc_err_cutoff = 1.00;
            stn = verify_variable(stn,[hcdTdt,'_err_1_d_sum']);
            baddts = [baddts;...
                      stn.([hcdTdt,'_err_1_d_sum']).date(abs(stn.([hcdTdt,'_err_1_d_sum']).data)>hc_err_cutoff)];


            baddays = unique(floor(baddts));
            disp(['** Days with anomalous errors removed: ',num2str(numel(baddays)),' **']);
            %DEBUG:
            find_date_ranges(baddays,1);
            %DEBUG:            disp('** NOT REMOVED **'); if (0);
            if ( ~isempty(baddays) )











            %DEBUG:            dump_station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld);





            %DEBUG:            disp(['** Days with anomalous errors NOT removed: ',num2str(numel(baddays)),' **']);




  skd2T = 4e-4.*abs(kd2T);					s2kd2T = signedSQR(skd2T);



              %DEBUG:
              if ( numel(baddays) > 100 )
                disp(datestr(baddays([1,end])));
              else
                disp(datestr(baddays));
              end;



  sKperHr = 1.5e-4;




  % Adjust heat storage errors for uncertainty of T (0.5K), S (1.0PSU), and h (18m) in K/hr calculation
  sKperHour = 2e-4;






  [dix,six,aix,qaix,qsix,Wix,hix,lix,Hsix,qlhix,qshix,qswiix,qswix,qlwiix,...
   sq0ix,sqtix,qbix,qbtix,bq0ix,bq0tix,fqudTix,kd2Tix,qtix,hchcix,hcix] = ...
      intersect_all_dates([],...


                          stn.(qtAdvfld).date,...


  qt = stn.(qtAdvfld).data(qtix);





  ix=ix+1; nm{ix}='Q_S_W/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qsw,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  SSE(ix)=sum((qsw.data-dt.data).^2); R2(ix)=1-(SSE(ix)./(N.*nrm)); RMSE(ix)=sqrt(SSE(ix)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='\gammaQ_S_W/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(aqsw,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  SSE(ix)=sum((aqsw.data-dt.data).^2); R2(ix)=1-(SSE(ix)./(N.*nrm)); RMSE(ix)=sqrt(SSE(ix)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='Q_L_W/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qlw,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  SSE(ix)=sum((qlw.data-dt.data).^2); R2(ix)=1-(SSE(ix)./(N.*nrm)); RMSE(ix)=sqrt(SSE(ix)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='Q_L_H/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qlh,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  SSE(ix)=sum((qlh.data-dt.data).^2); R2(ix)=1-(SSE(ix)./(N.*nrm)); RMSE(ix)=sqrt(SSE(ix)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='Q_S_H/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qsh,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  SSE(ix)=sum((qsh.data-dt.data).^2); R2(ix)=1-(SSE(ix)./(N.*nrm)); RMSE(ix)=sqrt(SSE(ix)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='Q_0/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(sq0,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  SSE(ix)=sum((sq0.data-dt.data).^2); R2(ix)=1-(SSE(ix)./(N.*nrm)); RMSE(ix)=sqrt(SSE(ix)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='(Q_0(\gamma)+Q_b)/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(bq0,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  SSE(ix)=sum((bq0.data-dt.data).^2); R2(ix)=1-(SSE(ix)./(N.*nrm)); RMSE(ix)=sqrt(SSE(ix)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='(Q_0(\gamma)+Q_b)/\rhoC_ph+u_s_f_c^.\nablaT_k_m+K_h\nabla^2T_k_m';
  [B,Stats]=scatter_fit_ts(bdt,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  SSE(ix)=sum((bdt.data-dt.data).^2); R2(ix)=1-(SSE(ix)./(N.*nrm)); RMSE(ix)=sqrt(SSE(ix)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='\partial_tT_s';
  [B,Stats]=scatter_fit_ts(hcdt,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  SSE(ix)=sum((hcdt.data-dt.data).^2); R2(ix)=1-(SSE(ix)./(N.*nrm)); RMSE(ix)=sqrt(SSE(ix)./N); SI(ix)=RMSE(ix)./nrm;










  RMSE(ix)=sqrt(sum((aqsw.data-dt.data)^2)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='Q_L_W/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qlw,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  RMSE(ix)=sqrt(sum((qlw.data-dt.data)^2)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='Q_L_H/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qlh,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  RMSE(ix)=sqrt(sum((qlh.data-dt.data)^2)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='Q_S_H/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(qsh,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  RMSE(ix)=sqrt(sum((qsh.data-dt.data)^2)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='Q_0/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(sq0,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  RMSE(ix)=sqrt(sum((sq0.data-dt.data)^2)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='(Q_0(\gamma)+Q_b)/\rhoC_ph';
  [B,Stats]=scatter_fit_ts(bq0,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  RMSE(ix)=sqrt(sum((bq0.data-dt.data)^2)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='(Q_0(\gamma)+Q_b)/\rhoC_ph+u_s_f_c^.\nablaT_k_m+K_h\nabla^2T_k_m';
  [B,Stats]=scatter_fit_ts(bdt,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  RMSE(ix)=sqrt(sum((bdt.data-dt.data)^2)./N); SI(ix)=RMSE(ix)./nrm;

  ix=ix+1; nm{ix}='\partial_tT_s';
  [B,Stats]=scatter_fit_ts(hcdt,dt,[],[],[STNM,' \Sigma',cumstr,nm{ix}],['\DeltaT_s/\Delta',cumstr],fh,[],true);
  rR2(ix)=Stats.regress_stats(1); rRMSE(ix)=Stats.s; rSI(ix)=Stats.s./nrm; rbias(ix)=B(1); rslope(ix)=B(2);
  RMSE(ix)=sqrt(sum((hcdt.data-dt.data)^2)./N); SI(ix)=RMSE(ix)./nrm;






  % Make sure calculated and climatology vectors are the same length!
  [ig,ix]=intersect(res.t.date,res.climsr.date);
  res.climdate = res.climdate(ix);
  res.climsr = res.climsr(ix);
  res.climlr = res.climlr(ix);
  res.climqlh = res.climqlh(ix);
  res.climqsh = res.climqsh(ix);





%%%%%%%%%%%%%%%%%%%%%
%% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%
function res = c_f_c_pare_result(res,dts)
  flds = fieldnames(res);
  for fld=flds{:}';
    dtix = 
    res.
  end;
return;

function res = c_f_c_grp_ts(res,fld,ts,per,sumfun,minN)
  [res.(fld),dt] = grp_ts(ts.data,ts.date,per,sumfun,minN);
return;




  if ( strcmp(per,'yearly') )
    xdat = res.date;
  else
    xdat = res.date;
  end;
  xmrkdat = xdat([1:mrk:end-1,end]);





  if ( strcmp(per,'yearly') )
    xdat = res.date;
  else
    xdat = res.date;
  end;
  mrkix = [1:mrk:length(xdat)-1,length(xdat)];
  xmrkdat = xdat([1:mrk:end-1,end]);


  fh = fmg;

if (0)
  ax(1)=subplot_tight(3,1,[1]);
  hold on; grid on;
  % plot(xdat,[res.t,res.t(1)+res.sqt-res.sqt(1),res.t(1)+res.bq0t-res.bq0t(1),res.t(1)+res.hc_dTdt-res.hc_dTdt(1)]);
  % legend(ax(1),'T ','T(0)+\SigmaQ_0/\rhoC_ph ','T(0)+\Sigma(Q_0(\gamma)+Q_b)/\rhoC_ph ','T(0)+\Sigma\partial_tT_H_C ',...
  %        'Location','SouthEast', 'Orientation','horizontal');
  plot(xdat,res.t,'k-',...
       xdat,res.t(1)+res.bq0t-res.bq0t(1),'b-',...
       xdat,res.t(1)+res.hc_dTdt-res.hc_dTdt(1),'r-');
  plot(xmrkdat,res.t(1:mrk:end),'k.',...
       xmrkdat,res.t(1)+res.bq0t(1:mrk:end)-res.bq0t(1),'bs',...
       xmrkdat,res.t(1)+res.hc_dTdt(1:mrk:end)-res.hc_dTdt(1),'ro');
  plh=plot(1,res.t(1),'k.-',...
           1,res.t(1)+res.bq0t(1)-res.bq0t(1),'bs-',...
           1,res.t(1)+res.hc_dTdt(1)-res.hc_dTdt(1),'ro-');
  legend(plh,'T','T(0)+\Sigma(Q_0(\gamma)+Q_b)/\rhoC_ph','T(0)+\Sigma\partial_tT_H_C',...
         'Location','SouthEast', 'Orientation','horizontal');
  axis([min(xdat),max(xdat), 5,45]);

  titlename([upper(stn.station_name) ': ' ...
             strrep(bdTfld,'_','\_') ' ' upper(per) ' ' upper(char(sumfun)) ' climatology ' ...
             stn.commentstr]);

  ax(2)=subplot_tight(3,1,[2,3]);
  hold on; grid on;
end;

  % PLOT (v.) does not accept property-value pairs as non-final arguments!
  if ( ~strcmp(per,'daily') )
    % Daily implied actual flux is way too noisy - blots out other plots
    plot(xdat,res.dtf,'k-','LineWidth',2);
  end;
  plot(xdat,res.radif,'r-',...
       xdat,res.turif,'b-',...
       xdat,res.aradif,'m-',...
       xdat,res.qbo,'k-.',...
       xdat,res.climradif,'r:',...
       xdat,res.climturif,'b:');
  if ( ~strcmp(per,'daily') )
    % Daily implied actual flux is way too noisy - blots out other plots
    plot(xmrkdat,res.dtf(1:mrk:end),'k.','LineWidth',2);
  end;
  plot(xmrkdat,res.radif(1:mrk:end),'rs',...
       xmrkdat,res.turif(1:mrk:end),'bs',...
       xmrkdat,res.aradif(1:mrk:end),'m^',...
       xmrkdat,res.qbo(1:mrk:end),'k+',...
       xmrkdat,res.climradif(1:mrk:end),'ro',...
       xmrkdat,res.climturif(1:mrk:end),'bo');
  plh=[];
  legs={};
  if ( ~strcmp(per,'daily') )
    % Daily implied actual flux is way too noisy - blots out other plots
    plh(end+1)=plot(1,res.dtf(1),'k.-','LineWidth',2);
    legs(end+1) = {'Actual'};
  end;
  plh(end+1:end+6) = ...
      plot(1,res.radif(1),'rs-',...
           1,res.turif(1),'bs-',...
           1,res.aradif(1),'m^-',...
           1,res.qbo(1),'k-.+',...
           1,res.climradif(1),'ro:',...
           1,res.climturif(1),'bo:');
  legs(end+1:end+6) = { ...
      'G&M: Q_S_W+Q_L_W',...
      'G&M: Q_L_H+Q_S_H',...
      'G&M: \gammaQ_S_W+Q_L_W',...
      'G&M: Q_b',...
      'ISCCP: Q_S_W+Q_L_W',...
      'OAFlux: Q_L_H+Q_S_H',...
                      };
  legh=legend(plh,legs);
  set(legh,'FontSize',7);
  axis([min(xdat),max(xdat), -400,400]);









  % Allowing leap-days might cause problems when comparing cum stats on
  % *different time series*, e.g., a time series which includes data for
  % one or more 29th's of Feb, and one which does not. If caller does not
  % like this behavior, they can simply specify custom CUMFUN and MINN.
  switch ( per ),
   case 'hourly',   n=365*24; minN=1;      mrk=24; N=1;
   case 'daily',    n=365;    minN=23;     mrk=30; N=24;
   case 'pentad',   n=73;     minN=23*4;   mrk= 5; N=24*5;
   case 'weekly',   n=52;     minN=23*6;   mrk= 4; N=24*7;
   case 'monthly',  n=12;     minN=23*21;  mrk= 1; N=24*[31,28,31,30,31,30,31,31,30,31,30,31]';
   case 'seasonal', n=4;      minN=23*68;  mrk= 1; N=24*[31+28+31,30+31+30,31+31+30,31+30+31]';
   % With interannual comparison, everything is different...
   case 'yearly',   n=nyrs;   minN=23*273; mrk= 1; N=1;
   otherwise,       error('Do not know how to do a "%s"-period climatology!',char(per));
  end;






%%%%???DEBUG
  % % Gross quality control on the relative errors
  % sQsh(~isfinite(sQsh)|0>sQsh|sQsh>1) = 1;
  % sQlh(~isfinite(sQlh)|0>sQlh|sQlh>1) = 1;
%%%%???DEBUG

  stn.([qshfld '_relerr']).date = dts(:);
  stn.([qshfld '_relerr']).data = sQsh(:);
  stn.([qlhfld '_relerr']).date = dts(:);
  stn.([qlhfld '_relerr']).data = sQlh(:);

  % Above calculates relative (normalized) errors - convert to absolute errors!
  sQsh = sQsh(:) .* abs(hs(:));
  s2Qsh = signedSQR(sQsh);
  sQlh = sQlh(:) .* abs(hl(:));
  s2Qlh = signedSQR(sQlh);


  stn.([qshfld '_err']).date = dts(:);
  stn.([qshfld '_err']).data = sQsh(:);
  stn.([qlhfld '_err']).date = dts(:);
  stn.([qlhfld '_err']).data = sQlh(:);






      % Plot histograms of estimated representation error for terms of Q0
      fmg;
      spt(2,2,1); hist(stn.([qlhfld,'_err']).data(stn.([qlhfld,'_err']).data<500),1000); xlim([0,80]);
      annotation('textbox',[0.40,0.90,.01,.01],'String','(a)','LineStyle','none','FontSize',14,'FontWeight','bold');
      spt(2,2,2); hist(stn.([qshfld,'_err']).data(stn.([qshfld,'_err']).data<500),1000); xlim([0,80]);
      annotation('textbox',[0.40,0.44,.01,.01],'String','(b)','LineStyle','none','FontSize',14,'FontWeight','bold');
      spt(2,2,3); hist(stn.([srfld,'_err']).data(0<stn.([srfld,'_err']).data&stn.([srfld,'_err']).data<500),1000); xlim([0,80]);
      annotation('textbox',[0.84,0.90,.01,.01],'String','(c)','LineStyle','none','FontSize',14,'FontWeight','bold');
      spt(2,2,4); hist(stn.([lrfld,'_err']).data(stn.([lrfld,'_err']).data<500),1000); xlim([0,80]);
      annotation('textbox',[0.84,0.44,.01,.01],'String','(d)','LineStyle','none','FontSize',14,'FontWeight','bold');
      suptitlename([STNM,' ',strrep(sq0fld,'_','\_'),' error histograms']);
      if ( doPrint )
        print('-dtiff',fullfile(figspath,[stnm,'-error-histograms-',sq0fld,'.tiff']));
      end;






  % % Best of Markovic et al (2009) and regression: 1% relative error, 7 W/m^2 bias
  % sQswi = 7 + (0.01.*qswi);					s2Qswi = signedSQR(sQswi);






%%%%???DEBUG
stn=safe_rmfield(stn,'ndbc_hfbulk_net_heat_flux');
%%%%???DEBUG
  if ( ~isfield(stn,'ndbc_hfbulk_net_heat_flux') )




  % Climatology has one value per period!
  raw_climsr = subset_ts(stn.(climsrfld),@(x)(find(ismember(get_year(x.date),unique(get_year(t.date))))));
  [climsr,climlr,climqlh,climqsh] = ...
      intersect_tses(raw_climsr,stn.(climlrfld),stn.(climqlhfld),stn.(climqshfld));




  if ( strcmp(per,'yearly') )
    % xdat = unique(get_year(t.date));
    xdat = unique(get_year(res.date));
  else
    xdat = 1:n;
  end;
  xmrkdat = min(xdat):mrk:max(xdat);






  % Climatology has one value per period!
  deldt = 0.75*min(diff(stn.(climsrfld).date));
  [ig,climsr,climlr,climqlh,climqsh] = ...
      intersect_tses(deldt,t,stn.(climsrfld),stn.(climlrfld),stn.(climqlhfld),stn.(climqshfld));





   case 'DRYF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    % doWarms = [true,false];
    doWarms = [false];
    kds = { ...
        [0.300,0.700,273] ...
        [0.300,0.700,319] ...
        [0.300,0.700, 45] ... %Best so far?
        [0.300,0.700,  0] ...
        [0.200,0.800,  0] ...
        [0.100,0.900,  0] ...
        [0.050,0.950,  0] ...
          };
    advfacs = { ...
        [0.00,1.00, 45] ... %Best so far
              };
    kths = { ...
        [0, 2, 45] ...
           };




   case 'LOOE1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          % [0.050,0.400,113] ...
          [0.050,0.400,137] ... %BEST FOR 0, 0, all years
          % [0.025,0.400,113] ... %BEST FOR 0, [0,20,45], sans 2008
          % [0.025,0.400,137] ...
            };
      advfacs = { ...
          % 0 ...                 % BEST FOR BOTH
          % [0.00,0.25, 45] ...
          [0.00,0.50, 45] ...
          % [0.00,0.50, 91] ...
          % [0.00,0.50,137] ...
          % [0.00,0.50,182] ...
          % [0.00,1.00, 91] ...
          % [0.00,1.00,137] ...
          % [0.00,1.00,182] ...
                };
      kths = { ...
          0 ...                 %BEST FOR 0, 0, all years
          % [0,20, 45] ...        %BEST FOR 0, [0,20,45], sans 2008
             };

    elseif ( ~isempty(regexp(sfld,'adcp')) || ~isempty(regexp(sfld,'ad_')) )

      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          [0.010,0.400,113] ... %BEST SO FAR
            };
      advfacs = { ...
          [0.00,0.50, 45] ...
                };
      kths = { ...
          [0,20, 45] ...
             };

    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;





   case 'DRYF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    % doWarms = [true,false];
    doWarms = [false];
    kds = { ...
        0.5 ...
        [0.300,0.700,  0] ...
        [0.300,0.700, 91] ...
        [0.300,0.700, 45] ... %Best so far
        [0.200,0.800, 45] ...
        [0.100,0.900, 45] ...
        [0.050,0.950, 45] ...
          };
    advfacs = { ...
        [0.00,1.00, 45] ... %Best so far
              };
    kths = { ...
        [0, 2, 45] ...
           };





   case 'DRYF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    % doWarms = [true,false];
    doWarms = [false];
    kds = { ...
        [0.050,0.950, 45] ...
        [0.100,0.900, 45] ...
        0.4 ...
        0.5 ...
        0.6 ...
        [0.300,0.500, 45] ...
        [0.300,0.700, 45] ... %Best so far
        [0.300,0.900, 45] ...
          };
    advfacs = { ...
        [0.00,1.00, 45] ... %Best so far
              };
    kths = { ...
        [0, 2, 45] ...
           };




  % 2012 Mar 31: TRIAL
  if ( strcmpi(stn.station_name,'dryf1') )
    stn.(mhfld).data = stn.(mhfld).data - 1.0;
    stn.(tufld).data = stn.(tufld).data .* 2.5;
    stn.(tvfld).data = stn.(tvfld).data .* 2.5;
  end;



  stn = station_tmd_tide(stn);
  % 2012 Mar 31: TRIAL
  %stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
  stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld,tufld,tvfld);




   case 'DRYF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.300,0.500, 45] ... %Best so far
        [0.200,0.600,  0] ...
        [0.200,0.800,  0] ...
        [0.200,1.000,  0] ...
        [0.200,1.200,  0] ...
        [0.100,0.600,  0] ...
        [0.100,0.800,  0] ...
        [0.100,1.000,  0] ...
          };
    advfacs = { ...
        [0.00,1.00, 45] ... %Best so far
        [0.00,1.00, 91] ...
        [0.00,0.50, 45] ...
        [0.00,0.50, 91] ...
        0 ...
              };
    kths = { ...
        [0,20, 45] ...
        [0,20, 91] ...
        [0,10, 45] ...
        [0,10, 91] ...
        0 ... %Best so far
           };





   case 'SMKF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.100,0.400,110] ...    % "Best" for *all* HC scalings [US][SU]
        [0.100,0.600, 60] ...    % Best for NO GRADIENT case - Best so far?
        [0.100,0.500, 60] ...    % Best for NO GRADIENT case - Best so far?
        [0.200,0.500, 60] ...    % Best for NO GRADIENT case - Best so far?
        [0.100,0.600, 70] ...    % Best for NO GRADIENT case - Best so far?
          };
    advfacs = { ...
        [0.0,0.25, 91] ...
              };
    kths = { ...
        [0, 2, 91] ...
           };





  switch ( lower(stn.station_name) ),
   case 'looe1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
%%%%DEBUG??? TRIAL
      % bad_years=[2009];
    end;
   case 'mlrf1',
    if ( ~isempty(regexp(sfld,'sea_t')) )
    end;
   case 'smkf1',
    if ( ~isempty(regexp(sfld,'sea_t')) )
      bad_years=[2004];
    end;
  end;




   case 'looe1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
%%%%DEBUG??? TRIAL
      % bad_years=[2009];
    end;



    if ( ~isempty(regexp(sfld,'sea_t')) )
      % bad_years=[2004,2007,2009];
      bad_years=[2004];
    end;



   case 'SMKF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    % stn.opts.hc_scaling = 'SS';
    % stn.opts.hc_scaling = 'SU';
    % stn.opts.hc_scaling = 'UU';
    doWarms = [true];
    kds = { ...
        [0.100,0.400,110] ...    % "Best" for *all* HC scalings [US][SU]
        [0.100,0.600, 60] ...    % Best for NO GRADIENT case - Best so far?
          };
    advfacs = { ...
        0 ...
        % % [0.0,0.5, 45] ...
        % [0.0,0.5, 91] ...
              };
    kths = { ...
        0 ...
        % % [0,10, 45] ...
        % [0,10, 91] ...
           };





  %% Filter out periods when the heat budget does NOT work well...

  bad_years=[];

  switch ( lower(stn.station_name) ),
   case 'looe1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      bad_years=[2009];
    end;
   case 'mlrf1',
    if ( ~isempty(regexp(sfld,'sea_t')) )
      % bad_years=[1998,2000,2005,2010];
      % bad_years=[2005,2010];
    end;
   case 'smkf1',
    if ( ~isempty(regexp(sfld,'sea_t')) )
%%%%DEBUG??? TRIAL
      % bad_years=[2004,2007,2009];
    end;
  end;

  if ( ~isempty(find(ismember(get_year(stn.(sfld).date),bad_years))) )
    warning('** Removing %s year(s) %s',sfld,num2str(bad_years(:)'));
    stn.commentstr = [stn.commentstr,' (sans ',num2str(bad_years(:)'),') '];
    for badyr = bad_years(:)';
      stn.(sfld).data(get_year(stn.(sfld).date)==badyr)=[];
      stn.(sfld).date(get_year(stn.(sfld).date)==badyr)=[];
    end;
  end;






   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.075,0.275,110] ...     % Best for WW3 1999-2011?
        [0.050,0.350, 45] ...
        [0.050,0.350, 91] ...
        [0.050,0.400, 45] ...
        [0.050,0.400, 91] ... %Best so far
          };
    advfacs = { ...
        [0.0,1.0, 45] ...
        [0.0,1.0, 91] ... %Best so far
              };
    kths = { ...
        [0,20, 45] ...
        [0,20, 91] ... %Best so far
           };






   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        % [0.075,0.275,110] ...     % Best for WW3 1999-2011?
        [0.050,0.400, 91] ... %Best so far
        [0.050,0.400,110] ...
          };
    advfacs = { ...
        [0.0,1.0, 45] ... %Best so far
        [0.0,1.0, 91] ...
              };
    kths = { ...
        [0,20, 91] ... %Best so far
           };




   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.010,0.300, 81] ... %Best so far
          };
    advfacs = { ...
        0 ...
        % [0.00,0.25, 45] ...
        % [0.00,0.25, 91] ... %Best so far
              };
    kths = { ...
        % [0, 5, 91] ...
        % [0, 5, 45] ...
        [0, 2, 91] ... %Best so far
        % [0, 2, 45] ...
           };

    %RMSE 2.2?
    kds = { ...
        [0.010,0.300, 81] ... %Best so far
          };
    advfacs = { ...
        [0.00,0.25, 91] ... %Best so far
              };
    kths = { ...
        [0, 2, 45] ...
           };

    %RMSE 2.5??
    kds = { ...
        [0.010,0.300, 45] ...
          };
    advfacs = { ...
        0 ...
              };
    kths = { ...
        [0, 2, 45] ...
           };








   case 'LOOE1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          % [0.050,0.400,113] ...
          [0.050,0.400,137] ... %BEST FOR 0, 0, all years
          % [0.025,0.400,113] ... %BEST FOR 0, [0,20,45], sans 2008
          % [0.025,0.400,137] ...
            };
      advfacs = { ...
          % 0 ...                 % BEST FOR BOTH
          % [0.00,0.25, 45] ...
          [0.00,0.50, 45] ...
          % [0.00,0.50, 91] ...
          % [0.00,0.50,137] ...
          % [0.00,0.50,182] ...
          % [0.00,1.00, 91] ...
          % [0.00,1.00,137] ...
          % [0.00,1.00,182] ...
                };
      kths = { ...
          0 ...                 %BEST FOR 0, 0, all years
          % [0,20, 45] ...        %BEST FOR 0, [0,20,45], sans 2008
             };

    elseif ( ~isempty(regexp(sfld,'adcp')) || ~isempty(regexp(sfld,'ad_')) )

      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          [0.010,0.400,113] ... %BEST SO FAR
            };
      advfacs = { ...
          [0.00,0.50, 45] ...
                };
      kths = { ...
          [0,20, 45] ...
             };

    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;









   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.010,0.300, 81] ... %Best so far
        [0.050,0.300, 81] ...
        [0.010,0.300, 45] ...
        [0.050,0.300, 45] ...
          };
    advfacs = { ...
        [0.00,0.25,120] ... %Best so far
        [0.00,0.25, 91] ...
        [0.00,0.25, 45] ...
        0 ...
              };
    kths = { ...
        [0,20, 91] ...
        [0,20, 45] ...
        [0,10, 91] ...
        [0,10, 45] ...
        0 ...
           };







   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.075,0.275,110] ...     % Best for ERAI waves 1993-2011 and WW3 1999-2011
        [0.050,0.400, 91] ... %Best so far
        [0.050,0.400,110] ...
        [0.050,0.400,137] ...
          };
    advfacs = { ...
        [0.0,1.0, 45] ... %Best so far
        [0.0,0.5, 91] ...
        [0.0,1.0, 91] ...
              };
    kths = { ...
        [0,20, 45] ... 
        [0,10, 45] ...
        [0,20, 91] ... %Best so far
        [0,10, 91] ...
           };





    % % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
    % %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
    % %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);




   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.700,1.200, 20] ... % Best so far
        [0.650,1.250, 20] ...
        [0.650,1.300, 20] ...
        [0.600,1.300, 20] ...
          };
    advfacs = { ...
        [0.0,1.0, 68] ...
        [0.0,1.0, 91] ... % Best so far
              };
    kths = { ...
        [0, 2, 91] ... % We want SOME diffusion if we can have it
           };






   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.700,1.200, 20] ... % Best so far
        [0.650,1.250, 20] ...
          };
    advfacs = { ...
        [0.0,1.0, 45] ... % Best so far
        [0.0,1.0, 91] ...
              };
    kths = { ...
        [0, 2, 91] ... % We want SOME diffusion if we can have it
           };





   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.700,1.200, 20] ... % Best so far
        [0.700,1.100, 20] ...
          };
    advfacs = { ...
        [0.0,1.0, 45] ... % Best so far
              };
    kths = { ...
        0 ... % Best so far
        [0,10, 45] ...
        [0,10, 91] ...
        [0, 5, 45] ...
        [0, 5, 91] ...
           };





   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.800,1.100, 20] ...
        [0.800,1.100, 81] ...
        [0.700,1.200, 20] ... % Best so far
        [0.700,1.200, 81] ...
        [0.700,1.100, 20] ...
        [0.700,1.100, 81] ...
          };
    advfacs = { ...
        [0.0,0.5, 45] ...
        [0.0,1.0, 45] ...
              };
    kths = { ...
        0 ...
        [0,10, 45] ...
        [0,10, 91] ...
        [0, 5, 45] ...
           };








   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.075,0.275,110] ...     % Best for ERAI waves 1993-2011 and WW3 1999-2011
        [0.050,0.400, 91] ...
        [0.050,0.400,110] ...
        [0.050,0.400,137] ...
          };
    advfacs = { ...
        [0.0,1.0, 45] ...
        [0.0,0.5, 45] ...
        0 ...
              };
    kths = { ...
        [0,20, 45] ...
        0 ...
           };





   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.800,1.100, 20] ... % Best so far
        [0.800,1.100, 81] ...
        [0.800,1.100,142] ...
        [0.800,1.100,203] ...
        [0.800,1.100,274] ...
        [0.800,1.100,335] ...
          };
    advfacs = { ...
        0 ...
        [0.0,0.5, 45] ...
        [0.0,1.0, 45] ...
              };
    kths = { ...
        0 ...
        [0,10, 45] ...
        [0,20, 45] ...
           };





   case 'FWYF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.025,0.300,110] ...     % Good for WW3 1999-2011 *and* ERAI waves 1987-2011
          };
    advfacs = { ...
        [0.00,1.00, 45] ...
        [0.00,1.00, 91] ...
              };
    kths = { ...
        [0,10, 45] ...
        [0,20, 91] ...
           };



   case 'FWYF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.025,0.300,110] ...     % Good for WW3 1999-2011 *and* ERAI waves 1987-2011
          };
    advfacs = { ...
        [0.00,1.00, 45] ...
        [0.00,0.50, 45] ...
        [0.00,0.25, 45] ...
        0 ...
              };
    kths = { ...
        [0,20, 45] ...
        [0,10, 45] ...
        0 ...
           };










  switch ( STNM ),

   case 'FWYF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.025,0.300,110] ...     % Good for WW3 1999-2011 *and* ERAI waves 1987-2011
          };
    advfacs = { ...
        [0.0,1.0, 45] ...
        [0.00,1.00, 91] ...
              };
    kths = { ...
        [0,20, 45] ...
           };

   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.075,0.275,110] ...     % For AVHRR, and ERAI waves 1993-2011 or WW3 waves 1999-2011
          };
    advfacs = { ...
        [0.0,1.0, 45] ...
              };
    kths = { ...
        [0,20, 45] ...
           };

   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.700,1.000, 20] ... % Best so far
        [0.700,1.000, 81] ...
        [0.700,1.000,142] ...
        [0.700,1.000,203] ...
        [0.700,1.000,274] ...
        [0.700,1.000,335] ...
        [0.800,1.100, 20] ...
          };
    advfacs = { ...
        1 ...
        [0.0,1.0, 45] ...
              };
    kths = { ...
        [0,20, 45] ...
           };


   case 'SMKF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    % stn.opts.hc_scaling = 'SS';
    % stn.opts.hc_scaling = 'SU';
    % stn.opts.hc_scaling = 'UU';
    doWarms = [true];
    kds = { ...
        % [0.100,0.400,110] ...    % "Best" for *all* HC scalings [US][SU]
        [0.100,0.600, 60] ...    % Best for NO GRADIENT case
          };
    advfacs = { ...
        0 ...
              };
    kths = { ...
        0 ...
           };


   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.0100,0.3000, 81] ...
          };
    advfacs = { ...
        [0.00,1.00, 45] ...
        [0.00,0.25,120] ...
              };
    kths = { ...
        [0,20, 45] ...
           };


   case 'DRYF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.300,0.500, 20] ...
        [0.300,0.500, 45] ...
        [0.300,0.500,335] ...
        [0.300,0.500,360] ...
        ...
        [0.200,0.600, 20] ...
        [0.200,0.600, 45] ...
        [0.200,0.600,335] ...
        [0.200,0.600,360] ...
          };
    advfacs = { ...
        [0.00,1.00, 45] ...
              };
    kths = { ...
        [0,20, 45] ...
           };


   case 'LOOE1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          % [0.050,0.400,113] ...
          [0.050,0.400,137] ... %BEST FOR 0, 0, all years
          % [0.025,0.400,113] ... %BEST FOR 0, [0,20,45], sans 2008
          % [0.025,0.400,137] ...
            };
      advfacs = { ...
          % 0 ...                 % BEST FOR BOTH
          % [0.00,0.25, 45] ...
          [0.00,0.50, 45] ...
          % [0.00,0.50, 91] ...
          % [0.00,0.50,137] ...
          % [0.00,0.50,182] ...
          % [0.00,1.00, 91] ...
          % [0.00,1.00,137] ...
          % [0.00,1.00,182] ...
                };
      kths = { ...
          0 ...                 %BEST FOR 0, 0, all years
          % [0,20, 45] ...        %BEST FOR 0, [0,20,45], sans 2008
             };

    elseif ( ~isempty(regexp(sfld,'adcp')) || ~isempty(regexp(sfld,'ad_')) )

      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          [0.010,0.400,113] ... %BEST SO FAR
            };
      advfacs = { ...
          [0.00,0.50, 45] ...
                };
      kths = { ...
          [0,20, 45] ...
             };

    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;


   case 'HAWK1',
    error('Hawk Channel SFP station HAWK1 *not yet implemented*');

   otherwise,
    error('Station %s options not implemented yet!',STNM);
  end;











   case 'SMKF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    % stn.opts.hc_scaling = 'SS';
    % stn.opts.hc_scaling = 'SU';
    % stn.opts.hc_scaling = 'UU';
    doWarms = [true];
    kds = { ...
        [0.100,0.400,110] ...    %"Best" for *all* HC scalings [US][SU]
        [0.100,0.550, 70] ...    %"Best" for NO GRADIENT case
        [0.050,0.550, 70] ...
        [0.100,0.600, 60] ...
        [0.050,0.600, 60] ...
          };
    advfacs = { ...
        % [0.0,0.5, 45] ...
        0 ...
              };
    kths = { ...
        % [0,20, 45] ...
        0 ...
           };








        [0.100,0.550, 70] ...    %"Best" for NO GRADIENT case
        [0.050,0.550, 70] ...


        [0.100,0.450, 90] ...
        [0.100,0.500, 80] ...



        [0.100,0.550, 70] ...    %"Best" for NO GRADIENT case



  %%%% ??? DEBUG
  if ( doADCP )
    if ( ~isfield(stn,'adcp_u') || ~isfield(stn,'adcp_v') )
      disp(['Processing ADCP currents ** from LOOE1 ** for advection']);
      x = get_looe1_adcp;
      stn.adcp_u = x.adcp_u;
      stn.adcp_v = x.adcp_v;
      x=[]; clear x;
    else
      disp(['Processing ADCP currents for advection']);
    end;
  end;
  %%%% ??? DEBUG




      %%%% DEBUG???
      % If we have in situ currents, TRY them for km-scale advection!
      if ( doADCP )
        ufld = 'adcp_u_40_h_lp';
        vfld = 'adcp_v_40_h_lp';
        stn = verify_variable(stn,{ufld,vfld});
        hufld = [ufld];
        hvfld = [vfld];
        stn.opts.km_scale_advection = true;
        stn.opts.calculate_advection = get_opt(stn.opts,'calculate_advection',true);
        stn.opts.calculate_diffusion = true;
        disp('** (INSITU+STOKES) ADVECTION **');
      else
      %%%% DEBUG???
        stn.opts.km_scale_advection = false;
        stn.opts.calculate_advection = get_opt(stn.opts,'calculate_advection',true);
        stn.opts.calculate_diffusion = true;
        disp('** ONLY STOKES ADVECTION **');
      end;







%%%%???DEBUG
ignore_benthos=true;





        ...
        [0.050,0.250,137] ...
        [0.050,0.250,110] ...
        [0.050,0.250, 73] ...





%%%%???DEBUG: REMOVE "PROBLEM YEARS" FOR DEBUGGING
% SMKF1
% stn.(sfld).data(stn.(sfld).date<datenum(1999,1,1))=[];
% stn.(sfld).date(stn.(sfld).date<datenum(1999,1,1))=[];
% stn.(sfld).data(stn.(sfld).date<datenum(1998,1,1))=[];
% stn.(sfld).date(stn.(sfld).date<datenum(1998,1,1))=[];
% stn.(sfld).data(stn.(sfld).date<datenum(1997,1,1))=[];
% stn.(sfld).date(stn.(sfld).date<datenum(1997,1,1))=[];
% stn.(sfld).data(stn.(sfld).date<datenum(1995,1,1))=[];
% stn.(sfld).date(stn.(sfld).date<datenum(1995,1,1))=[];
% FWYF1
% stn.(sfld).data(get_year(stn.(sfld).date)==1996)=[];
% stn.(sfld).date(get_year(stn.(sfld).date)==1996)=[];
% LOOE1 - ADCP Ts
% stn.(sfld).data(get_year(stn.(sfld).date)==2009)=[];
% stn.(sfld).date(get_year(stn.(sfld).date)==2009)=[];
% stn.(sfld).data(datenum(2007,4,1)<=stn.(sfld).date&stn.(sfld).date<datenum(2007,5,1))=[];
% stn.(sfld).date(datenum(2007,4,1)<=stn.(sfld).date&stn.(sfld).date<datenum(2007,5,1))=[];
%%%%???DEBUG






   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        % [0.065,0.275,130] ...
        [0.065,0.275,110] ...     % For AVHRR, and ERAI waves 1993-2011 or WW3 waves 1999-2011
        [0.075,0.275,110] ...
        [0.065,0.275, 91] ...
        [0.075,0.275, 91] ...
        % [0.065,0.275, 70] ...
        % ...
        % [0.100,0.300,130] ...
        % [0.100,0.300, 50] ...
        % [0.100,0.300, 91] ...
          };
    advfacs = { ...
        [0.0,1.0, 45] ...
              };
    kths = { ...
        [0,20, 45] ...
           };








%%%%???DEBUG
stn=safe_rmfield(stn,'ndbc_hfbulk_net_heat_flux');
%%%%???DEBUG





  %%%%???DEBUG
  stn.bogus_bulk_srf = stn.erai_srf;
  stn.bogus_bulk_lrf = stn.erai_lrf;
  stn = station_ndbc_hfbulk(stn,'bogus_bulk_srf','bogus_bulk_lrf','erai_relhumid');
  %%%%???DEBUG






  % [six,qix,gix,eix,nix,oix,lix,xix] = ...
  %     intersect_all_dates([],...
  %                         stn.(dsffld).date,stn.(fluxfld).date,...
  %                         stn.gom_hycom_net_heat_flux.date,...
  %                         stn.erai_actual_net_heat_flux.date,...
  %                         stn.ncep_net_heat_flux.date,...
  %                         stn.daily_oaflux_net_heat_flux.date,...
  %                         stn.landy_net_heat_flux.date,...
  %                         stn.nocs_net_heat_flux.date);




        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...




   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, *NO* gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.065,0.265,110] ...     % Good for WW3 waves 1999-2011
        ...
        [0.065,0.275,130] ...
        [0.065,0.275,110] ...     % For ERAI waves and AVHRR gradients 1993-2011
        [0.065,0.275, 50] ...
        [0.065,0.275, 70] ...
        [0.065,0.275, 91] ...
        ...
        [0.100,0.300,130] ...
        [0.100,0.300, 50] ...
        [0.100,0.300, 91] ...
          };
    advfacs = { ...
        0 ...
        [0,1, 45] ...
        1 ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
        [0,20, 45] ...
           };











            baddays = unique(floor(baddts));
%%%%DEBUG???
baddays=[];
%%%%DEBUG???
            if ( ~isempty(baddays) )






  [t,q0,qt,sq0,sqt,bq0,bq0t,bdT,hc_dTdt,sr,asr,lr,qlh,qsh,qrh,qbo,climsr,climlr,climqlh,climqsh] = ...
      intersect_tses(stn.(sfld),stn.(q0fld),stn.(qtfld),stn.(sqtfld),stn.(bq0tfld),stn.(sq0fld),stn.(bq0fld),stn.(bdTfld),stn.(hcdTdt),stn.(srfld),stn.(asrfld),stn.(lrfld),stn.(qlhfld),stn.(qshfld),stn.(qrhfld),stn.(qbofld),stn.(climsrfld),stn.(climlrfld),stn.(climqlhfld),stn.(climqshfld));








function x = chkann(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr)
%function x = chkann(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,commentstr)
%
% Plot annual or interannual climatologies for key heat budget terms
%
% Last Saved Time-stamp: <Tue 2012-03-27 17:14:30  Lew.Gramer>

  if ( ~exist('per','var') || isempty(per) )
    per = 'daily';
  end;
  if ( exist('commentstr','var') && ~isempty(commentstr) )
    stn.commentstr = commentstr;
  elseif ( ~isfield(stn,'commentstr') )
    stn.commentstr = '';
  end;

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

  % Would prefer MEDIAN, but published climatologies generally use MEAN (and
  % for, e.g., daily insolation, these give very different results!)
  sumfun = @nanmean;

  yrs = get_year(stn.(sfld).date);
  uyrs = unique(yrs);
  nyrs = numel(uyrs);

  % Allowing leap-days might cause problems when comparing cum stats on
  % *different time series*, e.g., a time series which includes data for
  % one or more 29th's of Feb, and one which does not. If caller does not
  % like this behavior, they can simply specify custom CUMFUN and MINN.
  switch ( per ),
   case 'hourly',   n=365*24; minN=1;      mrk=24; N=1;
   case 'daily',    n=365;    minN=23;     mrk=30; N=24;
   case 'pentad',   n=73;     minN=23*5;   mrk= 5; N=24*5;
   case 'weekly',   n=52;     minN=23*7;   mrk= 4; N=24*7;
   case 'monthly',  n=12;     minN=23*28;  mrk= 1; N=24*[31,28,31,30,31,30,31,31,30,31,30,31]';
   case 'seasonal', n=4;      minN=23*90;  mrk= 1; N=24*[31+28+31,30+31+30,31+31+30,31+30+31]';
   % With interannual comparison, everything is different...
   case 'yearly',   n=nyrs;   minN=23*200; mrk= 1; N=1;
   otherwise,       error('Do not know how to do a "%s"-period climatology!',char(per));
  end;

  % [t,q0,qt,sq0,sqt,bq0,bq0t,bdT,hc_dTdt,sr,asr,lr,qlh,qsh,qrh,qbo] = ...
  %     intersect_tses();


  %SEE function [cum,tid] = grp_ts(dat,dts,per_or_cumfun,sumfun,minN)

  x.t = grp_ts(stn.(sfld).data,stn.(sfld).date,per,sumfun,minN);

  x.q0 = grp_ts(stn.(q0fld).data,stn.(q0fld).date,per,sumfun,minN);
  x.qt = N.*cumsum(grp_ts(stn.(qtfld).data,stn.(qtfld).date,per,sumfun,minN));
  x.sq0 = grp_ts(stn.(sq0fld).data,stn.(sq0fld).date,per,sumfun,minN);
  x.sqt = N*cumsum(grp_ts(stn.(sqtfld).data,stn.(sqtfld).date,per,sumfun,minN));
  x.bq0 = grp_ts(stn.(bq0fld).data,stn.(bq0fld).date,per,sumfun,minN);
  x.bq0t = N.*cumsum(grp_ts(stn.(bq0tfld).data,stn.(bq0tfld).date,per,sumfun,minN));

  x.bdT = N.*cumsum(grp_ts(stn.(bdTfld).data,stn.(bdTfld).date,per,sumfun,minN));
  x.hc_dTdt = N.*cumsum(grp_ts(stn.(hcdTdt).data,stn.(hcdTdt).date,per,sumfun,minN));

  x.sr = grp_ts(stn.(srfld).data,stn.(srfld).date,per,sumfun,minN);
  x.asr = grp_ts(stn.(asrfld).data,stn.(asrfld).date,per,sumfun,minN);
  x.lr = grp_ts(stn.(lrfld).data,stn.(lrfld).date,per,sumfun,minN);
  x.qlh = grp_ts(stn.(qlhfld).data,stn.(qlhfld).date,per,sumfun,minN);
  x.qsh = grp_ts(stn.(qshfld).data,stn.(qshfld).date,per,sumfun,minN);
  x.qrh = grp_ts(stn.(qrhfld).data,stn.(qrhfld).date,per,sumfun,minN);

  x.qbo = grp_ts(stn.(qbofld).data,stn.(qbofld).date,per,sumfun,minN);

  x.radif = x.sr + x.lr;
  x.aradif = x.asr + x.lr;
  x.turif = x.qlh + x.qsh + x.qrh;
  x.coolif = x.lr + x.turif;

  % Climatology has one value per period!
  climkeepix = find(stn.(hcdTdt).date(1) <= stn.(climsrfld).date & stn.(climsrfld).date <= stn.(hcdTdt).date(end));
  x.climsr = grp_ts(stn.(climsrfld).data(climkeepix),stn.(climsrfld).date(climkeepix),per,sumfun,1);
  x.climlr = grp_ts(stn.(climlrfld).data(climkeepix),stn.(climlrfld).date(climkeepix),per,sumfun,1);
  x.climqlh = grp_ts(stn.(climqlhfld).data(climkeepix),stn.(climqlhfld).date(climkeepix),per,sumfun,1);
  x.climqsh = grp_ts(stn.(climqshfld).data(climkeepix),stn.(climqshfld).date(climkeepix),per,sumfun,1);

  x.climq0 = x.climsr + x.climlr + x.climqlh + x.climqsh;

  x.climradif = x.climsr + x.climlr;
  x.climturif = x.climqlh + x.climqsh;


  fh = fmg;

if (0)
  ax(1)=subplot_tight(3,1,[1]);
  hold on; grid on;
  % plot(1:n,[x.t,x.t(1)+x.sqt-x.sqt(1),x.t(1)+x.bq0t-x.bq0t(1),x.t(1)+x.hc_dTdt-x.hc_dTdt(1)]);
  % legend(ax(1),'T ','T(0)+\SigmaQ_0/\rhoC_ph ','T(0)+\Sigma(Q_0(\gamma)+Q_b)/\rhoC_ph ','T(0)+\Sigma\partial_tT_H_C ',...
  %        'Location','SouthEast', 'Orientation','horizontal');
  plot(1:n,x.t,'k-',...
       1:n,x.t(1)+x.bq0t-x.bq0t(1),'b-',...
       1:n,x.t(1)+x.hc_dTdt-x.hc_dTdt(1),'r-');
  plot(1:mrk:n,x.t(1:mrk:end),'k.',...
       1:mrk:n,x.t(1)+x.bq0t(1:mrk:end)-x.bq0t(1),'bs',...
       1:mrk:n,x.t(1)+x.hc_dTdt(1:mrk:end)-x.hc_dTdt(1),'ro');
  plh=plot(1,x.t(1),'k.-',...
           1,x.t(1)+x.bq0t(1)-x.bq0t(1),'bs-',...
           1,x.t(1)+x.hc_dTdt(1)-x.hc_dTdt(1),'ro-');
  legend(plh,'T','T(0)+\Sigma(Q_0(\gamma)+Q_b)/\rhoC_ph','T(0)+\Sigma\partial_tT_H_C',...
         'Location','SouthEast', 'Orientation','horizontal');
  xlim([1,n]); ylim([5,45]);
  titlename([upper(stn.station_name) ': ' ...
             strrep(bdTfld,'_','\_') ' ' upper(per) ' ' upper(char(sumfun)) ' climatology ' ...
             stn.commentstr]);

  ax(2)=subplot_tight(3,1,[2,3]);
  hold on; grid on;
end;

  plot(1:n,x.radif,'r-',...
       1:n,x.turif,'b-',...
       1:n,x.aradif,'m-',...
       1:n,x.qbo,'k-.',...
       1:n,x.climradif,'r:',...
       1:n,x.climturif,'b:');
  plot(1:mrk:n,x.radif(1:mrk:end),'rs',...
       1:mrk:n,x.turif(1:mrk:end),'bs',...
       1:mrk:n,x.aradif(1:mrk:end),'m^',...
       1:mrk:n,x.qbo(1:mrk:end),'k+',...
       1:mrk:n,x.climradif(1:mrk:end),'ro',...
       1:mrk:n,x.climturif(1:mrk:end),'bo');
  plh=plot(1,x.radif(1),'rs-',...
           1,x.turif(1),'bs-',...
           1,x.aradif(1),'m^-',...
           1,x.qbo(1),'k-.+',...
           1,x.climradif(1),'ro:',...
           1,x.climturif(1),'bo:');
  legh=legend(plh,...
              'G&M: Q_S_W+Q_L_W',...
              'G&M: Q_L_H+Q_S_H',...
              'G&M: \gammaQ_S_W+Q_L_W',...
              'G&M: Q_b',...
              'ISCCP: Q_S_W+Q_L_W',...
              'OAFlux: Q_L_H+Q_S_H');
  axis([1,n, -350,350]);


  titlename([upper(stn.station_name) ': ' ...
             strrep(bdTfld,'_','\_') ' ' upper(per) ' ' upper(char(sumfun)) ' climatology ' ...
             stn.commentstr]);

  if ( nargout < 1 )
    x = []; clear x;
  end;

return;











  %SEE function [cum,tid] = grp_ts(dat,dts,per_or_cumfun,sumfun,minN)

  x.t = grp_ts(stn.(sfld).data,stn.(sfld).date,per,sumfun,minN);

  x.q0 = grp_ts(stn.(q0fld).data,stn.(q0fld).date,per,sumfun,minN);
  x.qt = N.*cumsum(grp_ts(stn.(qtfld).data,stn.(qtfld).date,per,sumfun,minN));
  x.sq0 = grp_ts(stn.(sq0fld).data,stn.(sq0fld).date,per,sumfun,minN);
  x.sqt = N*cumsum(grp_ts(stn.(sqtfld).data,stn.(sqtfld).date,per,sumfun,minN));
  x.bq0 = grp_ts(stn.(bq0fld).data,stn.(bq0fld).date,per,sumfun,minN);
  x.bq0t = N.*cumsum(grp_ts(stn.(bq0tfld).data,stn.(bq0tfld).date,per,sumfun,minN));

  x.bdT = N.*cumsum(grp_ts(stn.(bdTfld).data,stn.(bdTfld).date,per,sumfun,minN));
  x.hc_dTdt = N.*cumsum(grp_ts(stn.(hcdTdt).data,stn.(hcdTdt).date,per,sumfun,minN));

  x.sr = grp_ts(stn.(srfld).data,stn.(srfld).date,per,sumfun,minN);
  x.asr = grp_ts(stn.(asrfld).data,stn.(asrfld).date,per,sumfun,minN);
  x.lr = grp_ts(stn.(lrfld).data,stn.(lrfld).date,per,sumfun,minN);
  x.qlh = grp_ts(stn.(qlhfld).data,stn.(qlhfld).date,per,sumfun,minN);
  x.qsh = grp_ts(stn.(qshfld).data,stn.(qshfld).date,per,sumfun,minN);
  x.qrh = grp_ts(stn.(qrhfld).data,stn.(qrhfld).date,per,sumfun,minN);

  x.qbo = grp_ts(stn.(qbofld).data,stn.(qbofld).date,per,sumfun,minN);







  %SEE function [cum,tid] = grp_ts(dat,dts,per_or_cumfun,sumfun,minN)

  x.t = grp_ts(stn.(sfld).data,stn.(sfld).date,per,sumfun,minN);

  x.q0 = grp_ts(stn.(q0fld).data,stn.(q0fld).date,per,sumfun,minN);
  x.qt = N.*cumsum(grp_ts(stn.(qtfld).data,stn.(qtfld).date,per,sumfun,minN));
  x.sq0 = grp_ts(stn.(sq0fld).data,stn.(sq0fld).date,per,sumfun,minN);
  x.sqt = N*cumsum(grp_ts(stn.(sqtfld).data,stn.(sqtfld).date,per,sumfun,minN));
  x.bq0 = grp_ts(stn.(bq0fld).data,stn.(bq0fld).date,per,sumfun,minN);
  x.bq0t = N.*cumsum(grp_ts(stn.(bq0tfld).data,stn.(bq0tfld).date,per,sumfun,minN));

  x.bdT = N.*cumsum(grp_ts(stn.(bdTfld).data,stn.(bdTfld).date,per,sumfun,minN));
  x.hc_dTdt = N.*cumsum(grp_ts(stn.(hcdTdt).data,stn.(hcdTdt).date,per,sumfun,minN));

  x.sr = grp_ts(stn.(srfld).data,stn.(srfld).date,per,sumfun,minN);
  x.asr = grp_ts(stn.(asrfld).data,stn.(asrfld).date,per,sumfun,minN);
  x.lr = grp_ts(stn.(lrfld).data,stn.(lrfld).date,per,sumfun,minN);
  x.qlh = grp_ts(stn.(qlhfld).data,stn.(qlhfld).date,per,sumfun,minN);
  x.qsh = grp_ts(stn.(qshfld).data,stn.(qshfld).date,per,sumfun,minN);
  x.qrh = grp_ts(stn.(qrhfld).data,stn.(qrhfld).date,per,sumfun,minN);

  x.qbo = grp_ts(stn.(qbofld).data,stn.(qbofld).date,per,sumfun,minN);









  % Allowing leap-days might cause problems when comparing cum stats on
  % *different time series*, e.g., a time series which includes data for
  % one or more 29th's of Feb, and one which does not. If caller does not
  % like this behavior, they can simply specify custom CUMFUN and MINN.
  switch ( per ),
   case 'hourly',   n=365*24; minN=1;      mrk=24; N=1;
   case 'daily',    n=365;    minN=23;     mrk=30; N=24;
   case 'pentad',   n=73;     minN=23*5;   mrk= 5; N=24*5;
   case 'weekly',   n=52;     minN=23*7;   mrk= 4; N=24*7;
   case 'monthly',  n=12;     minN=23*28;  mrk= 1; N=24*[31,28,31,30,31,30,31,31,30,31,30,31]';
   case 'seasonal', n=4;      minN=23*90;  mrk= 1; N=24*[31+28+31,30+31+30,31+31+30,31+30+31]';
   case 'yearly',   n=nyrs;   minN=23*350; mrk= 1; N(mod(uyrs,4)~=0)=365*24; N(mod(uyrs,4)==0)=366*24;
   otherwise,       error('Do not know how to do a "%s"-period climatology!',char(per));
  end;





   case 'yearly',   n = nyrs;   minN=23*350; mrk =  1; N(mod(uyrs,4)~=0)=365*24; N(mod(uyrs,4)==0)=366*24;


   case 'yearly',   n = nyrs;   minN=23*350; mrk =  1; N = repmat((365*24),[numel(uyrs),1]); N(mod(uyrs,4)==0) = 366*24;






% LOOE1 - MC Ts
% stn.(sfld).data(get_year(stn.(sfld).date)==2008)=[];
% stn.(sfld).date(get_year(stn.(sfld).date)==2008)=[];



   case 'LOOE1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          % [0.050,0.350,137] ...
          % [0.050,0.350,113] ...
          % [0.050,0.400,113] ...



       [0.025,0.400,113] ...
       [0.020,0.400,113] ...
   


   case 'LOOE1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          [0.050,0.350,137] ...
          [0.050,0.400,137] ...
          [0.050,0.350,113] ...
          [0.050,0.400,113] ...
          [0.025,0.400,113] ... %BEST FOR ADCP TEMP
            };
      advfacs = { ...
          0 ...
          % [0.00,0.50, 45] ...
          % [0.00,0.50, 91] ...
          % [0.00,0.50,137] ...
          % [0.00,0.50,182] ...
          % [0.00,1.00, 91] ...
          % [0.00,1.00,137] ...
          % [0.00,1.00,182] ...
                };
      kths = { ...
          0 ...
          % { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };

    elseif ( ~isempty(regexp(sfld,'adcp')) || ~isempty(regexp(sfld,'ad_')) )

      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
       [0.025,0.400,113] ... %BEST SO FAR
            };
      advfacs = { ...
          [0.00,0.50, 45] ... %BEST SO FAR
          [0.00,1.00, 45] ...
                };
      kths = { ...
          0 ...
          [0,20, 45] ...
          % { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ... %BEST SO FAR
             };

    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;






   case 'LOOE1',
    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          [0.050,0.350,137] ...
          [0.050,0.400,137] ...
          [0.050,0.350,113] ...
          [0.050,0.400,113] ...
          [0.025,0.400,113] ... %BEST FOR ADCP TEMP
            };
      advfacs = { ...
          0 ...
          % [0.00,0.50, 45] ...
          % [0.00,0.50, 91] ...
          % [0.00,0.50,137] ...
          % [0.00,0.50,182] ...
          % [0.00,1.00, 91] ...
          % [0.00,1.00,137] ...
          % [0.00,1.00,182] ...
                };
      kths = { ...
          0 ...
          % { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };

    elseif ( ~isempty(regexp(sfld,'adcp')) || ~isempty(regexp(sfld,'ad_')) )

      % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
      doWarms = [true];
      kds = { ...
          [0.025,0.400,113] ... %BEST SO FAR
            };
      advfacs = { ...
          [0.00,0.50,319] ...
          [0.00,0.50,  0] ...
          [0.00,0.50, 91] ...
          [0.00,0.50,137] ...
          [0.00,0.50, 45] ...
          ...
          [0.00,1.00,319] ...
          [0.00,1.00,  0] ...
          [0.00,1.00, 91] ...
          [0.00,1.00,137] ...
          [0.00,1.00, 45] ... %BEST SO FAR
                };
      kths = { ...
          0 ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ... %BEST SO FAR
             };

    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;







            %% 
            % Final QC: remove points with anomalous error est. FROM THIS ESTIMATE

            % %NO! hc_err_cutoff = prctile(abs(stn.([hcdTdt,'_err']).data),99.99);
            % % PRCTILE would continually shrink STN.(SRFLD) with each loop!
            % % We still do that below, but usually only by a few points...
            % % hc_err_cutoff = 0.08;
            % hc_err_cutoff = 0.10;
            % baddts = stn.([hcdTdt,'_err']).date(abs(stn.([hcdTdt,'_err']).data)>hc_err_cutoff);





Fq=cos(2 pi * (year-day - 45)/366)


    if ( ~isempty(regexp(sfld,'microcat')) || ~isempty(regexp(sfld,'mc_')) )
      % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [true];
      kds = { ...
          [0.500,0.800,137] ...
          [0.500,0.800,113] ...
          [0.400,0.700,137] ...
          [0.400,0.700,113] ...
          [0.300,0.600,137] ...
          [0.300,0.600,113] ...
          [0.200,0.500,137] ...
          [0.200,0.500,113] ...
          [0.100,0.400,137] ...
          [0.100,0.400,113] ...
          [0.050,0.350,137] ...
          [0.025,0.400,113] ... %BEST FOR ADCP TEMP
            };
      advfacs = { ...
          0 ...
          [0.00,1.00,319] ...
          [0.00,1.00,  0] ...
          [0.00,1.00, 45] ...
          [0.00,1.00, 91] ...
                };
      kths = { ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };




    if ( strcmp(sfld,'microcat_seatemp') )
    elseif ( strcmp(sfld,'adcp_seatemp') )





   case 'LOOE1',
    if ( strcmp(sfld,'microcat_seatemp') )
      % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [false,true];
      kds = { ...
          [0.050,0.350,137] ...
          [0.025,0.400,105] ... %BEST FOR ADCP TEMP
            };
      advfacs = { ...
          [0.00,1.00, 91] ...
          1 ...
                };
      kths = { ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };
    elseif ( strcmp(sfld,'adcp_seatemp') )
      % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [true];
      kds = { ...
          [0.025,0.400,137] ...
          [0.025,0.400,129] ...
          [0.025,0.400,121] ...
          [0.025,0.400,113] ... %BEST SO FAR
          [0.025,0.400,105] ...
          [0.025,0.400, 97] ...
            };
      advfacs = { ...
          [0.00,1.00, 45] ...
                };
      kths = { ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };
    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;







    elseif ( strcmp(sfld,'adcp_seatemp') )
      % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [true];
      kds = { ...
          [0.025,0.400,137] ... %BEST SO FAR
          [0.050,0.350,137] ...
            };
      advfacs = { ...
          [0.00,1.00, 45] ...
          [0.00,1.00, 91] ...
          [0.00,1.00,137] ...
                };
      kths = { ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };
    else







          [0.100,0.300, 45] ...
          [0.100,0.300,137] ...
          [0.100,0.300,183] ...
          ....
          [0.200,0.400, 45] ...
          [0.200,0.400,137] ...
          [0.200,0.400,183] ...



          [0.050,0.350,137] ... %BEST SO FAR
          [0.050,0.350, 91] ...
          [0.050,0.350,182] ...
          [0.025,0.400,137] ...
          [0.025,0.400, 91] ...
          [0.025,0.400,182] ...



   case 'LOOE1',
    if ( strcmp(sfld,'microcat_seatemp') )
      % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [false,true];
      kds = { ...
          [0.050,0.350,137] ... %BEST SO FAR
          [0.050,0.350, 91] ...
          [0.050,0.350,182] ...
          [0.025,0.400,137] ...
          [0.025,0.400, 91] ...
          [0.025,0.400,182] ...
            };
      advfacs = { ...
          [0.00,1.00, 90] ...
                };
      kths = { ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };
    elseif ( strcmp(sfld,'adcp_seatemp') )
      % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [true];
      kds = { ...
          [0.050,0.350,137] ...
          [0.100,0.300, 45] ...
          [0.100,0.300,137] ...
          [0.100,0.300,183] ...
          ....
          [0.200,0.400, 45] ...
          [0.200,0.400,137] ...
          [0.200,0.400,183] ...
            };
      advfacs = { ...
          [0.00,1.00, 90] ...
          1 ...
                };
      kths = { ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };
    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;







   case 'LOOE1',
    if ( strcmp(sfld,'microcat_seatemp') )
      % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [false,true];
      kds = { ...
          [0.050,0.350,137] ...
          [0.100,0.300, 45] ...
          [0.100,0.300,137] ...
          [0.100,0.300,183] ...
          ....
          [0.200,0.400, 45] ...
          [0.200,0.400,137] ...
          [0.200,0.400,183] ...
            };
      advfacs = { ...
          [0.00,1.00, 90] ...
          1 ...
                };
      kths = { ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };
    elseif ( strcmp(sfld,'adcp_seatemp') )
      % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [true];
      kds = { ...
          [0.050,0.350,137] ...
          [0.100,0.300, 45] ...
          [0.100,0.300,137] ...
          [0.100,0.300,183] ...
          ....
          [0.200,0.400, 45] ...
          [0.200,0.400,137] ...
          [0.200,0.400,183] ...
            };
      advfacs = { ...
          [0.00,1.00, 90] ...
          1 ...
                };
      kths = { ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };
    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;






   case 'SMKF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    stn.opts.hc_scaling = get_opt(stn.opts,'hc_scaling','SS');
    doWarms = [true];
    kds = { ...
        [0.150,0.350,110] ...
        [0.100,0.400,110] ...
        [0.050,0.450,110] ...
        [0.050,0.400,110] ...
          };
    advfacs = { ...
        0 ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
        20 ...
           };




   case 'SMKF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.0100,0.2850, 81] ...
        ...
        0.3 ...
        ...
        [0.200,0.450,  0] ...
        [0.200,0.450, 30] ...
        [0.200,0.450, 60] ...
        [0.200,0.450, 90] ...
        [0.200,0.450,120] ...
        [0.200,0.450,137] ...
        [0.200,0.450,150] ...
        [0.200,0.450,180] ...
        [0.200,0.450,210] ...
        [0.200,0.450,240] ...
        [0.200,0.450,270] ...
        [0.200,0.450,300] ...
        [0.200,0.450,330] ...
          };
    advfacs = { ...
        0 ...
        [0.00,1.00, 90] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };






            %%%%DEBUG???
            %annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,[],[1:3]);







              %%%%???DEBUG
              % [cum,tid] = grp_ts(qe.data,qe.date,'daily',@nanmean,23);
              [cum,tid] = grp_ts(qe.data,qe.date,'daily',@(x)(nanmean(abs(x))),23);
              % [cum,tid] = grp_ts(qe.data,qe.date,'daily',@(x)(nanmax(abs(x))),23);
              cum = 24*cum;



%%%%???DEBUG
              % [cum,tid] = grp_ts(qe.data,qe.date,@floor,@nansum,24);
              [cum,tid] = grp_ts(qe.data,qe.date,@floor,@(x)(nansum(abs(x))),24);
%%%%???DEBUG






    % %%%% DEBUG???
    % maxGrad = 2e-4;
    % disp(['** Limiting SST gradients to ',num2str(maxGrad*1e3),'K/km **']);
    % stn.(hkmtxsfld).data(stn.(hkmtxsfld).data>maxGrad)=maxGrad;
    % stn.(hkmtxsfld).data(stn.(hkmtxsfld).data<-maxGrad)=-maxGrad;
    % stn.(hkmtlsfld).data(stn.(hkmtlsfld).data>maxGrad)=maxGrad;
    % stn.(hkmtlsfld).data(stn.(hkmtlsfld).data<-maxGrad)=-maxGrad;
    % %%%% DEBUG???






    % %%%% DEBUG???
    % maxAdv = 0.045;
    maxAdv = Inf;
    if ( maxAdv < Inf)
      disp(['** Limiting QC''d advection to ',num2str(maxAdv),'K/hr **']);
      stn = qa_ts(stn,udTfld);
      stn.(udTfld).data(stn.(udTfld).data>maxAdv) = maxAdv;
      stn.(udTfld).data(stn.(udTfld).data<-maxAdv) = -maxAdv;
    end;
    % %%%% DEBUG???








   case 'FWYF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.050,0.275,110] ...
        [0.035,0.265,110] ...     % Good for WW3 1999-2011 *and* ERAI waves 1987-2011
        [0.035,0.275,110] ...
        [0.025,0.275,110] ...
          };
    advfacs = { ...
        [0.00,1.00, 90] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };







   case 'SMKF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false,true];
    kds = { ...
        [0.0100,0.2850, 81] ...
        ...
        0.3 ...
        ...
        [0.200,0.450,  0] ...
        [0.200,0.450, 30] ...
        [0.200,0.450, 60] ...
        [0.200,0.450, 90] ...
        [0.200,0.450,120] ...
        [0.200,0.450,137] ...
        [0.200,0.450,150] ...
        [0.200,0.450,180] ...
        [0.200,0.450,210] ...
        [0.200,0.450,240] ...
        [0.200,0.450,270] ...
        [0.200,0.450,300] ...
        [0.200,0.450,330] ...
          };
    advfacs = { ...
        0 ...
        [0.00,1.00, 90] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };









   case 'LOOE1',
    if ( strcmp(sfld,'microcat_seatemp') )
      % SMKF1 NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [true];
      kds = { ...
          [0.050,0.350,137] ...
          [0.100,0.300, 45] ...
          [0.100,0.300,137] ...
          [0.100,0.300,183] ...
            };
      advfacs = { ...
          0,...
          1,...
                };
      kths = { ...
          0 ...
          5 ...
          [0,10, 45] ...
          [0,20, 45] ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };
    elseif ( strcmp(sfld,'adcp_seatemp') )
      % SMKF1 NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
      doWarms = [true];
      kds = { ...
          [0.050,0.350,137] ...
          [0.100,0.300, 45] ...
          [0.100,0.300,137] ...
          [0.100,0.300,183] ...
            };
      advfacs = { ...
          0,...
          1,...
                };
      kths = { ...
          0 ...
          5 ...
          [0,10, 45] ...
          [0,20, 45] ...
          { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
             };
    else
      error('No optimization parameters defined for LOOE1 "%s"',sfld);
    end;





   case 'FWYF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.035,0.265,110] ...     % Good for WW3 waves 1999-2011
        [0.075,0.275, 90] ...     % More options 1993-2011
        [0.075,0.275,130] ...
        [0.075,0.275,110] ...
        [0.050,0.275, 90] ...     % For ERAI waves and AVHRR gradients 1993-2011
        [0.050,0.275,130] ...
        [0.050,0.275,110] ...
          };
    advfacs = { ...
        [0.0,0.50, 90] ...
        [0.0,0.75, 90] ...
        [0.0,1.00, 90] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };







    stn.bulk_srf.date = stn.erai_dsrf.date;
    stn.bulk_srf.data = (1-stn.erai_ndbc_albedo.data).*stn.erai_dsrf.data;





DRYF1
        0.4 ...
        [0.300,0.500, 20] ...
        [0.300,0.500, 81] ...
        [0.300,0.500,142] ...
        [0.300,0.500,203] ...
        [0.300,0.500,274] ...
        [0.300,0.500,335] ...



   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.650,0.900, 20] ...
        [0.700,0.900, 20] ...
        [0.650,0.950, 20] ...
        [0.700,0.950, 20] ...
          };
    advfacs = { ...
        % 0 ...
        % [0.0,0.5, 90] ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        % 0 ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };






   case 'DRYF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.050,0.250,137] ...
        [0.500,0.700, 20] ...
        [0.300,0.500, 20] ...
        [0.100,0.300, 20] ...
        [0.010,0.300, 20] ...
        [0.500,0.700, 81] ...
        [0.300,0.500, 81] ...
        [0.100,0.300, 81] ...
        [0.010,0.300, 81] ...
          };
    advfacs = { ...
        0 ...
        [0.0,1.00, 90] ...
        1 ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };







   case 'DRYF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % stn.opts.kd = get_opt(stn.opts,'kd',[0.045,0.375,45]);
    % stn.opts.do_warm_layer = get_opt(stn.opts,'do_warm_layer',false);
    % stn.opts.convective_drag_coefficient = get_opt(stn.opts,'convective_drag_coefficient',8.0e-4);
    % stn.opts.K_theta = get_opt(stn.opts,'K_theta',0);

    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.20,0.40,320], };
    % kths = { 0 };

    % *MLRF* NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.050,0.250,137] ...
        [0.700,0.900, 20] ...
        [0.500,0.700, 20] ...
          };
    advfacs = { ...
        0 ...
        [0.0,1.00, 90] ...
        1 ...
              };
    kths = { ...
        0 ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% ??? DEBUG DEBUG DEBUG






MLRF1:
advfacs:
        [0.0,0.25, 90] ...
        [0.0,0.50, 90] ...
        [0.0,0.75, 90] ...
        [0.0,1.00, 90] ...
        [0.0,0.25,120] ...
        [0.0,0.50,120] ...
        [0.0,0.75,120] ...
        [0.0,1.00,120] ...




    figfname= fullfile(figspath,[lower(stn.station_name),'-chkann-',hcdTdt,'.']);
      print('-dpng',[figfname 'png']);





chkann.m:
  legend(plh,...
         [upper(RAPFX) ' Q_S_W + ' upper(RAPFX) '/insitu Q_L_W'],...
         'Q_L_H+Q_S_H',...
         ['\gamma' upper(RAPFX) ' Q_S_W + ' upper(RAPFX) '/insitu Q_L_W'],...
         'Modeled: Q_b^O',...
         'ISCCP: Q_S_W+Q_L_W',...
         'OAFlux: Q_L_H+Q_S_H');






        { stn.(tspdfld),@(T)(min(20,((T./2.0).^2).*20)) } ...







   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    % kds = { ...
    %     [0.0100,0.2850, 81] ...
    %     [0.0050,0.2850, 81] ...
    %     [0.0050,0.2900, 81] ...  %BEST SO FAR
    %       };
    kds = { ...
        [0.0050,0.3000, 76] ...
        [0.0100,0.3000, 81] ...
          };
    advfacs = { ...
        [0.00,0.25,120] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };








   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.0100,0.2850, 81] ...  %BEST SO FAR
        [0.0050,0.2850, 81] ...
          };
    advfacs = { ...
        [0.00,0.25,110] ...
        [0.00,0.30,110] ...
        [0.00,0.25,120] ...
        [0.00,0.30,120] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };








   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.0150,0.2850, 71] ...
          };
    advfacs = { ...
        % 0 ...
        [0.00,0.25, 90] ...
        % [0.00,0.50, 90] ...
        % [0.00,0.75, 90] ...
        % [0.00,0.10, 90] ...
        % [0.00,0.25,181] ...
        % [0.00,0.50,181] ...
        % [0.00,0.75,181] ...
        % [0.00,0.10,181] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };




%%%%???DEBUG
stn.(scoolfld) = stn.(sfld);
%%%%???DEBUG






    % % Plot relative contributions to total error
    % fmg; grpplot_ts(ts_op(stn.([asrfld,'_err']),stn.([q0fld,'_err']),'/'),[],[],[],'m'); grpplot_ts(ts_op(stn.([lrfld,'_err']),stn.([q0fld,'_err']),'/'),[],[],[],'r'); grpplot_ts(ts_op(stn.([qlhfld,'_err']),stn.([q0fld,'_err']),'/'),[],[],[],'g'); grpplot_ts(ts_op(stn.([qshfld,'_err']),stn.([q0fld,'_err']),'/'),[],[],[],'b'); legend('\sigma\gammaQ_S_W/\Sigma','\sigmaQ_L_W/\Sigma','\sigmaQ_L_H/\Sigma','\sigmaQ_S_H/\Sigma', 'Location','Best'); titlename([STNM,' relative error contribution (',strrep(q0fld,'_','\_'),')']);







   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    % %%%%???DEBUG
    % cbds = { ...
    %     8.0e-4 ...
    %        };
    % % kds = { ...
    % %     [0.600,1.100,355] ...
    % %     [0.600,1.200,355] ...
    % %     [0.600,1.300,355] ...
    % %     [0.600,1.350,355] ...
    % %     [0.600,1.400,355] ...
    % %     [0.500,1.400,355] ...
    % %     [0.500,1.500,355] ...
    % %     ...
    % %     [0.600,1.100,325] ...
    % %     [0.600,1.200,325] ...
    % %     [0.600,1.300,325] ...
    % %     [0.600,1.350,325] ...
    % %     [0.600,1.400,325] ...
    % %     [0.500,1.400,325] ...
    % %     [0.500,1.500,325] ...
    % %       };
    % % %%%%??? DEBUG Kds for 2000-2010 - huge!
    % % kds = { ...
    % %     [0.600,1.100,355] ...
    % %     [0.600,1.200,355] ...
    % %     [0.600,1.200,325] ...
    % %     [0.600,1.100,325] ...
    % %     [0.500,1.400,325] ...
    % %     [0.500,1.500,325] ...
    %%%%??? DEBUG Kds for 2000-2007
    kds = { ...
        [0.500,1.500,325] ...
        [0.500,0.900, 20,182.6225] ...
        [0.550,0.900, 20,182.6225] ...
        [0.600,0.900, 20,182.6225] ...
        ...
        [0.600,0.900, 20] ...
        [0.650,0.900, 20] ...
        [0.700,0.900, 20] ...
          };
    advfacs = { ...
        % 0 ...
        % [0.0,0.5, 90] ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        % 0 ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };









    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        % 0 ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
        { stn.(tspdfld),@(T)(min(20,((T./2.0).^2).*20)) } ...
           };





  if ( ~isfield(stn,sfld) )
    if ( regexp(sfld,'misst_sst') )
      stn = get_misst_station(stn);
      stn.hourly_misst_sst = interp_ts(stn.misst_sst);
    elseif ( regexp(sfld,'avhrr_weekly_sst') )
      stn = get_avhrr_weekly_field(stn,true,stn.opts.grid_interp_method,npts);
      stn.hourly_avhrr_weekly_sst = interp_ts(stn.avhrr_weekly_sst);
    elseif ( regexp(sfld,'erai_') )
      if ( ~strcmpi(ISPFX,'erai') && ~strcmpi(RAPFX,'erai') )
        stn = get_erai_station(stn);
        if ( adjust_reanalysis )
          stn = adjust_erai_station(stn);
        end;
      end;
    elseif ( regexp(sfld,'microcat') )
      stn = get_looe1_microcat(stn);
    elseif ( regexp(sfld,'adcp') )
      stn = get_looe1_adcp(stn);
    elseif ( regexp(sfld,'ndbc') )
      if ( ~strcmpi(ISPFX,'ndbc') )
        stn = load_all_ndbc_data(stn);
      end;
    elseif ( regexp(sfld,'^(sea|ct_|ctd_)') )
      if ( ~strcmpi(ISPFX,'icon') )
        stn = load_station_data(stn);
      end;
    end;
  end;









  % % ULRF error estimated by regressing QlwO vs. in situ blackbody radiation (say what??)
  % sQlwo = 4;							s2Qlwo = signedSQR(sQlwo);






   case 'FWYF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.035,0.265,110] ...
          };
    advfacs = { ...
        [0.0,0.5, 90] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };










   case 'MLRF1',
    % % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    % doWarms = [true];
    % kds = { ...
    %     [0.050,0.250,110] ...
    %       };
    % advfacs = { ...
    %     1 ...
    %           };
    % kths = { ...
    %     { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
    %        };
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, *NO* gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.065,0.265,110] ...
        [0.070,0.270,110] ...
          };
    advfacs = { ...
        1 ...
        [0.0,0.5,  0] ...
        [0.0,0.5, 90] ...
        [0.0,0.5,182] ...
        [0.0,0.5,273] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };






   case 'MLRF1',
    % % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    % doWarms = [true];
    % kds = { ...
    %     [0.050,0.250,110] ...
    %       };
    % advfacs = { ...
    %     1 ...
    %           };
    % kths = { ...
    %     { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
    %        };
    % NDBC met, ERAI barom/spechumid, new Ppen, ERAI waves, *NO* gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.065,0.265,110] ...
        [0.070,0.270,110] ...
          };
    advfacs = { ...
        1 ...
        [0.0,0.5,  0] ...
        [0.0,0.5, 90] ...
        [0.0,0.5,182] ...
        [0.0,0.5,273] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };









%%%%???DEBUG
% %% LONF1
% badix = find(datenum(2004,5,16)<=stn.(afld).date&stn.(afld).date<datenum(2004,5,27));
% stn.(afld).date(badix)=[]; stn.(afld).data(badix)=[];
% badix = find(datenum(2008,6,1)<=stn.(afld).date&stn.(afld).date<datenum(2008,7,1));
% stn.(afld).date(badix)=[]; stn.(afld).data(badix)=[];
%
% %% SANF1
% badix = find(datenum(2001,1,23)<=stn.(sfld).date&stn.(sfld).date<datenum(2001,1,26));
% stn.(sfld).date(badix)=[]; stn.(sfld).data(badix)=[];
% badix = find(datenum(2003,1,23)<=stn.(sfld).date&stn.(sfld).date<datenum(2003,1,26));
% stn.(sfld).date(badix)=[]; stn.(sfld).data(badix)=[];
%%%%???DEBUG




   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    stn = verify_variable(stn,[Ufld,'_3_d_lp']);
    stn = verify_variable(stn,[Ufld,'_30_d_lp']);
    stn = verify_variable(stn,[Ufld,'_90_d_lp']);
    kds = { ...
        [0.050,0.250,110] ...
        { stn.([Ufld,'_3_d_lp']),@(U3d)(sin((pi/2).*min(U3d,35)./35)+1) } ...
        { stn.([Ufld,'_30_d_lp']),@(U30d)(sin((pi/2).*min(U30d,35)./35)+1) } ...
        { stn.([Ufld,'_90_d_lp']),@(U90d)(sin((pi/2).*min(U90d,35)./35)+1) } ...
          };
    advfacs = { ...
        1 ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };





    kds = { ...
        [0.050,0.250,110] ...
        { stn.(Wfld),@(W)( min(0.25,((((W./35).^1).*0.20)+0.05)) ) } ...
        { stn.(Wlpfld),@(Wlp)( min(0.25,((((Wlp./35).^1).*0.20)+0.05)) ) } ...
          };



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% RADIATIVE FLUXES

  % Estimate short-wave radiative flux errors

% ??? VERIFY against Lewis, Wei, van Dommelen, Voss (2011)
% ??? (see same comment in STATION_ABSORBED_INSOLATION.m)

  % % Consistent with Markovic et al (2009), 5% relative error, 7 W/m^2 bias in ERA40 Qswi
  % sQswi = 7 + (0.05.*qswi);					s2Qswi = signedSQR(sQswi);
  % % Based on regressing FWYF1 ERAI vs. Lew QC'd RSMAS in situ DAILY AVERAGES: bias +10, RMSE 40
  % sQswi = 50 + (0.08.*qswi);					s2Qswi = signedSQR(sQswi);

  % Based on regressing FWYF1 ERAI vs. Lew QC'd RSMAS in situ DAILY AVERAGES: bias -3, RMSE 40
  %sQswi = 43 + (0.01.*qswi);					s2Qswi = signedSQR(sQswi);
  sQswi = 7 + (0.01.*qswi);					s2Qswi = signedSQR(sQswi);

  %%%%???DEBUG: Remove "error" from night-time values
  sQswi(qswi<1) = 0;						s2Qswi = signedSQR(sQswi);


  % Assuming Albedo has constant relative uncertainty of 4% *of Qswi*
  %sQsw = 43 + (0.05.*qswi);					s2Qsw = signedSQR(sQsw);
  sQsw = 7 + (0.05.*qswi);					s2Qsw = signedSQR(sQsw);

  %%%%???DEBUG: Remove "error" from night-time values
  sQsw(qswi<1) = 0;						s2Qsw = signedSQR(sQsw);







   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.050,0.250,110] ...
        { stn.(Wfld),@(W)( min(0.25,((((W./35).^2).*0.20)+0.05)) ) } ...
        { stn.(Wlpfld),@(Wlp)( min(0.25,((((Wlp./35).^2).*0.20)+0.05)) ) } ...
          };
    advfacs = { ...
        1 ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };





  [ig,r_sgamma_sQsw] = cov_ts(sgamma,sQsw,@isfinite,@isfinite);
%%%%???DEBUG
  r_sgamma_sQsw = 0;



ADJUST_WW3_STATION_WAVES.m:
if(0)
end;
  %%%%
  %% Significant wave height

  % Piecewise breakpoint in Peak Wave Period [s]
  cutoff_per = 8.0;

  % Regression against all seasons - piecewise by wave PERIOD
  [chopix,ig] = intersect_dates(stn.ww3_sigwavehgt.date,stn.ww3_peakwaveper.date(stn.ww3_peakwaveper.data<=cutoff_per));
  a = +0.104; 
  b = +0.884;  % RMSE ~ 0.16[m], r^2 ~ 0.73
  stn.ww3_sigwavehgt_adj.data(chopix) = ( (stn.ww3_sigwavehgt.data(chopix) .* b) + a );

  [rollix,ig] = intersect_dates(stn.ww3_sigwavehgt.date,stn.ww3_peakwaveper.date(stn.ww3_peakwaveper.data>cutoff_per));
  % Power2 CFIT
  a = +0.483; 
  b = +0.493;
  c = +0.219;  % RMSE ~ 0.13[m], r^2 ~ 0.58
  stn.ww3_sigwavehgt_adj.data(rollix) = ( ( a .* (stn.ww3_sigwavehgt.data(rollix) .^ b) ) + c );

  % LINEAR ADJ/INSITU FIT:
  %  Wv3: -0.012+1.03x r2~0.70 RMSE~0.15
  %  Wv4: -0.078+0.85x r2~0.43 RMSE~0.22
  %  Wv7: -0.062+1.04x r2~0.71 RMSE~0.15






  %%%%
  %% Significant wave height

  % Regression against all seasons - piecewise by wave PERIOD
  [chopix,ig] = intersect_dates(stn.ww3_sigwavehgt.date,stn.ww3_peakwaveper.date(stn.ww3_peakwaveper.data<=6.5));
  a = +0.115; 
  b = +0.892;  % RMSE ~ 0.15[m], r^2 ~ 0.77
  stn.ww3_sigwavehgt_adj.data(chopix) = ( (stn.ww3_sigwavehgt.data(chopix) .* b) + a );

  [rollix,ig] = intersect_dates(stn.ww3_sigwavehgt.date,stn.ww3_peakwaveper.date(stn.ww3_peakwaveper.data>6.5));
  % % Linear ROBUST
  % a = +0.330; 
  % b = +0.336;  % RMSE ~ 0.15[m], r^2 ~ 0.54
  % stn.ww3_sigwavehgt_adj.data(rollix) = ( (stn.ww3_sigwavehgt.data(rollix) .* b) + a );

  % Power2 CFIT
  a = +0.432; 
  b = +0.723;
  c = +0.243;  % RMSE ~ 0.15[m], r^2 ~ 0.54
  stn.ww3_sigwavehgt_adj.data(rollix) = ( ( a .* (stn.ww3_sigwavehgt.data(rollix) .^ b) ) + c );









        warning('No empirical adjustment defined for ERA-Interim Waves yet??');



  stn = station_tmd_tide(stn);
  stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);

%%%%???DEBUG
stn.(mhfld) = stn.(hfld);
%%%%???DEBUG




ANNSUBS.m:
  if ( ~exist('commentstr','var') || isempty(commentstr) )
    if ( isfield(stn,'commentstr') )
      commentstr = stn.commentstr;
    else
      commentstr = '';
    end;
  end;

  x.commentstr = commentstr;
  station_heat_budget_field_names;
  %% HACK ALERT: string COMMENTSTR is *not* a variable name...
  commentstr = x.commentstr;
  clear x






  % Some basic quality control (SHOULD have been done during initial extraction!)
  if ( isfield(newdat,'ww3_peakwaveper') )
    % Peak Period below 2s is NOT real - just means a model was spinning up!
    baddts = newdat.ww3_peakwaveper.date(newdat.ww3_peakwaveper.data<2);
    if ( ~isempty(baddts) )
      [badix,ig] = intersect_dates(newdat.ww3_sigwavehgt.date,baddts);
      newdat.ww3_sigwavehgt.data(badix) = nan;
      [badix,ig] = intersect_dates(newdat.ww3_peakwavedir.date,baddts);
      newdat.ww3_peakwavedir.data(badix) = nan;
      [badix,ig] = intersect_dates(newdat.ww3_peakwaveper.date,baddts);
      newdat.ww3_peakwaveper.data(badix) = nan;

      if ( isfield(newdat,'ww3_peakwave_u') && isfield(newdat,'ww3_peakwave_v') )
        [badix,ig] = intersect_dates(newdat.ww3_peakwave_u.date,baddts);
        newdat.ww3_peakwave_u.data(badix) = nan;
        [badix,ig] = intersect_dates(newdat.ww3_peakwave_v.date,baddts);
        newdat.ww3_peakwave_v.data(badix) = nan;
      end;
    end;
  end;








   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false,true];
    kds = { ...
        [0.0150,0.2850, 83] ...  % Both best for "No Warm Layer"
        [0.0150,0.2850, 71] ...
        ...
        [0.0150,0.2850, 60] ...  % Both best for "Warm Layer"
        [0.0125,0.2875, 71] ...
          };
    advfacs = { ...
        0 ...
        [0.0,0.25,  0] ...
        [0.0,0.25, 90] ...
        [0.0,0.25,181] ...
        [0.0,0.25,273] ...
              };
    kths = { ...
        0 ...
           };



   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false,true];
    kds = { ...
        [0.015,0.285, 83] ...    % Best for "No Warm Layer"
        ...
        [0.015,0.285, 60] ...    % Both best for "Warm Layer"
        [0.010,0.290, 83] ...
        ...
        [0.015,0.285, 71] ...    % More options to try
        [0.0125,0.2875, 71] ...
        [0.0125,0.2875, 83] ...
          };
    advfacs = { ...
        0 ...
              };
    kths = { ...
        0 ...
           };






   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false,true];
    kds = { ...
        [0.035,0.265, 45] ...    % All good for "No Warm Layer"
        [0.025,0.275, 45] ...
        ...
        [0.035,0.265, 60] ...
        [0.025,0.275, 60] ...
        [0.015,0.285, 60] ...
        ...
        [0.025,0.275, 83] ...
        [0.015,0.285, 83] ...
        [0.010,0.290, 83] ...
        ...
        [0.015,0.285, 45] ...    % Also good for "Warm Layer"
        [0.010,0.290, 60] ...
          };
    advfacs = { ...
        0 ...
              };
    kths = { ...
        0 ...
           };






   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false,true];
    kds = { ...
        [0.050,0.250, 45] ...
        [0.035,0.265, 45] ...
        [0.025,0.275, 45] ...
        [0.015,0.285, 45] ...
        [0.010,0.290, 45] ...
        ...
        [0.050,0.250, 60] ...
        [0.035,0.265, 60] ...
        [0.025,0.275, 60] ...
        [0.015,0.285, 60] ...
        [0.010,0.290, 60] ...
        ...
        [0.050,0.250, 83] ...
        [0.035,0.265, 83] ...
        [0.025,0.275, 83] ...
        [0.015,0.285, 83] ...
        [0.010,0.290, 83] ...
          };
    advfacs = { ...
        0 ...
        [0.0,0.5, 90] ...
              };
    kths = { ...
        0 ...
           };





   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false,true];
    kds = { ...
        [0.050,0.250, 83] ...
        [0.035,0.265, 83] ...
        [0.025,0.275, 83] ...
        [0.015,0.285, 83] ...
        ...
        [0.050,0.250,110] ...  %MLRF1
        [0.035,0.265,110] ...  %FWYF1
        [0.025,0.275,110] ...
        [0.015,0.285,110] ...
          };
    advfacs = { ...
        0 ...
        [0.0,0.5, 90] ...
        1 ...
              };
    kths = { ...
        0 ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
        { stn.(tspdfld),@(T)(min(20,((T./2.0).^2).*20)) } ...
           };








%%%%??? DEBUG
begyr = 1970;
% USF AVHRR 1km SST data (for Florida) before 1996 appears highly suspect
if ( regexp(KMPFX,'avhrr') ); begyr = 1996; end;
if ( regexp(sfld,'avhrr') ); begyr = 1996; end;
%%%%??? DEBUG - make periods for WAVEPFX='ww3','erai','ndbc' all match
begyr = 2000;
endyr = 2012;
% %%%%??? DEBUG - make periods consistent with SMKF1
% endyr = 2007;
begdt = datenum(begyr,1,2);
enddt = datenum(endyr,1,1);
badyrs = [];
if (~isempty(badyrs)); stn.commentstr = [stn.commentstr,' (ex ',num2str(badyrs,'%g,'),') ']; end;
t.data(begdt>t.date|t.date>enddt|ismember(get_year(t.date),badyrs))=[];
t.date(begdt>t.date|t.date>enddt|ismember(get_year(t.date),badyrs))=[];
q.data(begdt>q.date|q.date>enddt|ismember(get_year(q.date),badyrs))=[];
q.date(begdt>q.date|q.date>enddt|ismember(get_year(q.date),badyrs))=[];
begyr = min(get_year(t.date));
endyr = max(get_year(t.date));
%%%%??? DEBUG




'LONF1'
    % % Include tidal currents in advection estimate
    % stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);


    % % % % kds = { ...
    % % % %     0.7 ...
    % % % %     0.65 ...
    % % % %     [0.500,0.900,201,182.6225] ...
    % % % %     [0.500,0.850,201,182.6225] ...
    % % % %     [0.500,0.800,201,182.6225] ...
    % % % %     [0.550,0.900,201,182.6225] ...
    % % % %     [0.550,0.850,201,182.6225] ...
    % % % %     [0.550,0.800,201,182.6225] ...
    % % % %     [0.550,0.750,201,182.6225] ...
    % % % %       };
    % % % kds = { ...
    % % %     0.7 ...
    % % %     [0.550,0.900,201,182.6225] ...
    % % %     [0.600,0.900,201,182.6225] ...
    % % %     [0.600,0.850,201,182.6225] ...
    % % %     [0.700,0.800,201,182.6225] ...
    % % %     [0.700,0.850,201,182.6225] ...
    % % %     ...
    % % %     [0.600,0.800,201-91,182.6225] ...
    % % %     [0.600,0.850,201-91,182.6225] ...
    % % %     [0.650,0.850,201-91,182.6225] ...
    % % %     [0.700,0.850,201-91,182.6225] ...
    % % %     [0.700,0.900,201-91,182.6225] ...
    % % %     [0.600,0.900,201-91,182.6225] ...
    % % %     [0.650,0.900,201-91,182.6225] ...
    % % %     ...
    % % %     [0.600,0.800,201-182,182.6225] ...
    % % %     [0.600,0.850,201-182,182.6225] ...
    % % %     [0.650,0.850,201-182,182.6225] ...
    % % %     [0.700,0.850,201-182,182.6225] ...
    % % %     [0.700,0.900,201-182,182.6225] ...
    % % %     [0.600,0.900,201-182,182.6225] ...
    % % %     [0.650,0.900,201-182,182.6225] ...
    % % %       };
    % % kds = { ...
    % %     0.7 ...
    % %     [0.650,0.850,201-182,182.6225] ...
    % %     [0.750,0.850,201-182,182.6225] ...
    % %     [0.700,0.850,201-182,182.6225] ...
    % %     [0.700,0.800,201-182,182.6225] ...
    % %     [0.650,0.850,350,182.6225] ...
    % %     [0.750,0.850,350,182.6225] ...
    % %     [0.700,0.850,350,182.6225] ...
    % %     [0.700,0.800,350,182.6225] ...
    % %       };
    % kds = { ...
    %     0.65 ...
    %     0.70 ...
    %     0.75 ...
    %     ...
    %     [0.550,0.850,325] ...
    %     [0.550,0.850,355] ...
    %     [0.650,0.750,325] ...
    %     [0.650,0.750,355] ...
    %     ...
    %     [0.600,0.900,325] ...
    %     [0.600,0.900,355] ...
    %     [0.700,0.800,325] ...
    %     [0.700,0.800,355] ...
    %     ...
    %     [0.600,1.000,325] ...
    %     [0.600,1.000,355] ...
    %     [0.600,1.100,325] ...
    %     [0.600,1.100,355] ...
    %     [0.700,1.000,325] ...
    %     [0.700,1.000,355] ...
    %     [0.700,1.100,325] ...
    %     [0.700,1.100,355] ...
    %       };



%%%%???DEBUG
cum=cum(60:end);
tid=tid(60:end);
%%%%???DEBUG
%%%%???DEBUG
cum=cum(60:end);
tid=tid(60:end);
%%%%???DEBUG
%%%%???DEBUG
cum=cum(60:end);
tid=tid(60:end);
%%%%???DEBUG







From BIGMAP.m:
% % NO GOOD - even my PC TITAN cannot handle this size NGDC region!

% % For Gramer & Mariano
% lons=[sanf1.lon,looe1.lon,smkf1.lon,lonf1.lon,mlrf1.lon,fwyf1.lon]; lats=[sanf1.lat,looe1.lat,smkf1.lat,lonf1.lat,mlrf1.lat,fwyf1.lat];
% % % For ICRS paper
% % lons=[smkf1.lon,lonf1.lon,mlrf1.lon,fwyf1.lon]; lats=[smkf1.lat,lonf1.lat,mlrf1.lat,fwyf1.lat];
% 
% stn.station_name='mid_freef_ngdc';
% stn.lon=mean([min(lons),max(lons)]);
% stn.lat=mean([min(lats),max(lats)]);
% rawmatfname = fullfile(get_ecoforecasts_path('coast'),'LGramer1-80.mat');
% rawxyz = load(rawmatfname,'lon','lat','depth');
% stn = get_ngdc_bathy_station(stn,100,rawxyz);


fmg;

% % For ICRS paper
% %map_freef([mlrf1.lon-0.97,mlrf1.lon+0.97,mlrf1.lat-1.00,mlrf1.lat+1.00],isodepths);
% [ig,ig,cs,ch,cobj,clhs]=map_freef([stn.lon-0.66,stn.lon+0.66,stn.lat-0.60,stn.lat+0.60],isodepths);

% For Gramer & Mariano
% % %map_freef([smkf1.lon-1.10,smkf1.lon+1.10,smkf1.lat-1.00,smkf1.lat+1.00],isodepths);
% % %[ig,ig,cs,ch,cobj,clhs]=map_freef_ngdc([smkf1.lon-0.20,fwyf1.lon+0.10,smkf1.lat-0.15,fwyf1.lat+0.15],isodepths);
% % [ig,ig,cs,ch,cobj,clhs]=map_freef_ngdc([stn.lon-0.20,stn.lon+0.10,stn.lat-0.15,stn.lat+0.15],isodepths);
% % set(cobj,'LineWidth',1.5);
% % set(ch,'Color',isocolor);
% % set(clhs,'Color',isocolor);
% [ig,ig,cs,ch,cobj,clhs]=map_freef([stn.lon-1.10,stn.lon+1.10,stn.lat-1.00,stn.lat+1.00],'none');
[ig,ig,cs,ch,cobj,clhs]=map_freef([-82.1,-79.9,24.2,25.8],'none');
set(cobj,'LineWidth',1.5);







text(fwyf1.lon,fwyf1.lat,'. \leftarrow FWYF1','Rotation',+000, 'FontWeight','demi');
text(mlrf1.lon,mlrf1.lat,'. \leftarrow MLRF1','Rotation',-045, 'FontWeight','demi');
text(lonf1.lon,lonf1.lat,'. \leftarrow LONF1','Rotation',+030, 'FontWeight','demi');
text(smkf1.lon,smkf1.lat,'. \leftarrow SMKF1','Rotation',-060, 'FontWeight','demi');






   case 'SMKF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    % doWarms = [true];
    doWarms = [false];
    kds = { ...
        0.2 ...
        0.3 ...
        0.4 ...
        ...
        [0.050,0.250,110] ...
        ...
        [0.200,0.400,137] ...
        [0.300,0.500,137] ...
        [0.500,0.700,137] ...
        [0.200,0.400,  0] ...
        [0.300,0.500,  0] ...
        [0.500,0.700,  0] ...
        ...
        [0.250,0.550,330] ...
        [0.250,0.550,  0] ...
        [0.250,0.550, 30] ...
        [0.250,0.550, 60] ...
        [0.250,0.550, 90] ...
        [0.250,0.550,120] ...
          };
        % ...
        % [0.250,0.325,  0] ...
        % [0.250,0.325, 30] ...
        % [0.250,0.325, 45] ...
    % % Ignore unphysical SST gradients - see SFP TSG data analysis
    % stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    % % Ignore alongshore component of heat advection?
    % stn.opts.add_alongshore_advection = get_opt(stn.opts,'add_alongshore_advection',false);
    % % Include tidal currents in advection estimate
    % stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    advfacs = { ...
        0 ...
              };
        % [0.0,0.5, 90] ...
        % 1 ...
    kths = { ...
        0 ...
           };
        % { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...








    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    % doWarms = [true];
    doWarms = [false,true];
    kds = { ...
        [0.200,0.400,137] ...
        [0.300,0.500,137] ...
        [0.500,0.700,137] ...
        [0.200,0.400,  0] ...
        [0.300,0.500,  0] ...
        [0.500,0.700,  0] ...
        ...
        [0.250,0.325,  0] ...
        [0.250,0.325, 30] ...
        [0.250,0.325, 45] ...
          };
    % % Ignore unphysical SST gradients - see SFP TSG data analysis
    % stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    % % Ignore alongshore component of heat advection?
    % stn.opts.add_alongshore_advection = get_opt(stn.opts,'add_alongshore_advection',false);
    % Include tidal currents in advection estimate
    stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    advfacs = { ...
        0 ...
        [0.0,0.5, 90] ...
        1 ...
              };
        % [0, 1, 183] ...
    kths = { ...
        0 ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };
        % 2 ...
        % [0, 5, 45] ...







  % Assuming Albedo has constant relative uncertainty of 4% *of Qswi*
  sQsw = 1.04 .* sQswi;						s2Qsw = signedSQR(sQsw);



  % Assuming Albedo has constant relative uncertainty of 4% *of Qswi*
  sQsw = 0.04 .* sQswi;						s2Qsw = signedSQR(sQsw);






  % Time series of sea surface reflectances
  [albix,srfix] = intersect_dates(stn.(albfld).date,stn.(srfld).date(qswix));
  Albedo = stn.(albfld).data(albix);

  % Assuming Albedo is time-varying and with constant 4% relative uncertainty
  sAlbedo = 0.04 .* Albedo;					s2Albedo = signedSQR(sAlbedo);

  [ig,r_sAlbedo_sQswi] = cov_ts(sAlbedo,sQswi);

  s2Qsw = ( (Albedo2.*s2Qswi) + (qswi2.*s2Albedo) + (2.*Albedo.*qswi.*sAlbedo.*sQswi.*r_sAlbedo_sQswi)  );
  s2Qsw(qswi<0.1) = 0;
  sQsw = signedsqrt(s2Qsw);









  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% RADIATIVE FLUXES

  % Estimate short-wave radiative flux errors

% ??? VERIFY against Lewis, Wei, van Dommelen, Voss (2011)
% ??? (see same comment in STATION_ABSORBED_INSOLATION.m)

  % % Consistent with Markovic et al (2009), 5% relative error, 7 W/m^2 bias in ERA40 Qswi
  % sQswi = 7 + (0.05.*qswi);					s2Qswi = signedSQR(sQswi);
  % % Based on regressing FWYF1 ERAI vs. Lew QC'd RSMAS in situ DAILY AVERAGES: bias +10, RMSE 40
  % sQswi = 50 + (0.08.*qswi);					s2Qswi = signedSQR(sQswi);
  % Based on regressing FWYF1 ERAI vs. Lew QC'd RSMAS in situ DAILY AVERAGES: bias -3, RMSE 40
  sQswi = 43 + (0.01.*qswi);					s2Qswi = signedSQR(sQswi);


  % Time series of sea surface reflectances according to Reanalysis
  [albix,srfix] = intersect_dates(stn.(albfld).date,stn.(srfld).date(qswix));
  Albedo = stn.(albfld).data(albix);
  % %%%% DEBUG??? NEED qswoix!
  % Albedo = stn.(usrfld).data(qswoix) ./ qswi;
  % badix = (Albedo>0.5 | ~isfinite(Albedo));
  % Albedo(badix) = mean(Albedo(~badix));

  % HACK: Assuming Albedo is time-varying but with ZERO uncertainty
  s2Qsw = signedSQR(1 - Albedo) .* s2Qswi;
  s2Qsw(qswi<0.1) = 0;
  sQsw = signedsqrt(s2Qsw);

  % Assuming gamma is time-varying and with constant relative uncertainty
  % (100% - gamma) - normally ~4% relative uncertainty in net insolation.
  % See ESTIM_GAMMA_ERROR.m for the basis of this relative error estimate
  % %%%% ??? HACK
  % %[asrfix,srfix] = intersect_dates(stn.(asrfld).date,stn.(srfld).date(qswix));
  % %gamma = stn.(asrfld).data(asrfix) ./ stn.(srfld).data(qswix(srfix));
  % %gamma(~isfinite(gamma)) = mean(gamma(isfinite(gamma)));
  % %gamma(stn.(asrfld).data(asrfix) < 0.1) = 0;
  % %gamma2 = signedSQR(gamma);
  % %%%% ??? HACK
  [gamix,srfix] = intersect_dates(stn.(gamfld).date,stn.(srfld).date(qswix));
  gamma = stn.(gamfld).data(gamix);
  gamma(~isfinite(gamma)) = mean(gamma(isfinite(gamma)));
  gamma2 = signedSQR(gamma);

  sgamma = 1.00 .* (1 - gamma);
  s2gamma = signedSQR(sgamma);

  [ig,r_sgamma_sQsw] = cov_ts(sgamma,sQsw);

  s2gammaQsw = ( (gamma2.*s2Qsw) + (qsw2.*s2gamma) + (2.*gamma.*qswi.*sgamma.*sQsw.*r_sgamma_sQsw)  );
  sgammaQsw = signedsqrt(s2gammaQsw);







  
OLDER VERSION...

  % Assuming gamma is time-varying and with constant 100% relative
  % uncertainty (should be ~4% relative uncertainty in net insolation).
  % See ESTIM_GAMMA_ERROR.m for the basis of this relative error estimate
  % %%%% ??? HACK
  % [asrfix,srfix] = intersect_dates(stn.(asrfld).date,stn.(srfld).date(qswix));
  % gamma = stn.(asrfld).data(asrfix) ./ stn.(srfld).data(qswix(srfix));
  % gamma(~isfinite(gamma)) = mean(gamma(isfinite(gamma)));
  % gamma(stn.(asrfld).data(asrfix) < 0.1) = 0;
  % gamma2 = signedSQR(gamma);
  % %%%% ??? HACK
  [gamix,srfix] = intersect_dates(stn.(gamfld).date,stn.(srfld).date(qswix));
  gamma = stn.(gamfld).data(gamix);
  gamma(~isfinite(gamma)) = mean(gamma(isfinite(gamma)));
  gamma2 = signedSQR(gamma);

  s2gamma = 1.00 .* (1 - gamma2);
  sgamma = signedsqrt(s2gamma);





  if ( strcmp(sfld,'ndbc_sea_t_1da') )
    stn = verify_variable(stn, 'ndbc_sea_t_1_d_avg');
    stn.(sfld) = stn.ndbc_sea_t_1_d_avg;
  end;




%%
%% VARIABLES WHICH ARE SET BY THIS SCRIPT:
%%   bathyfld = 'ngdc_92m_bathy';
%%   slopefld = 'ngdc_offshore_slope';
%%   bathorifld = 'isobath_orientation';
%%   hfld = [TIDEPFX '_i_depth'];
%%   mhfld = ['mean_' hfld];
%%   tufld = [TIDEPFX '_u'];
%%   tvfld = [TIDEPFX '_v'];
%%   tspdfld = [TIDEPFX '_speed'];
%%   tdirfld = [TIDEPFX '_dir'];
%%
%%   sfld = [ISPFX '_sea_t'];
%%   afld = [ISPFX '_air_t'];
%%   pfld = [RAPFX '_barom'];
%%   dfld = [RAPFX '_dew_t'];
%%   rhfld = [RAPFX '_relhumid'];
%%   qafld = [RAPFX '_spechumid'];
%%   qsfld = [SPREFIX ISPFX '_sea_spechumid'];
%%
%%   WINDINFIX = '_wind';
%%   Wfld = [ISPFX WINDINFIX '_speed'];
%%   Dfld = [ISPFX WINDINFIX '_dir'];
%%   Ufld = [ISPFX WINDINFIX '_u'];
%%   Vfld = [ISPFX WINDINFIX '_v'];
%%   Ulpfld = [Ufld QE_LOWPASS];
%%   Vlpfld = [Vfld QE_LOWPASS];
%%   Wlpfld = [Wfld QE_LOWPASS];
%%   Dlpfld = [Dfld QE_LOWPASS];
%%
%%   rfld = [RAPFX '_precip'];
%%   cfld = [RAPFX '_cloud_cover'];
%%   pblzfld = [RAPFX '_pblz'];
%%   dsrfld = [RAPFX '_dsrf'];
%%   usrfld = [RAPFX '_usrf'];
%%   srfld = [RAPFX '_srf'];
%%   asrfld = ['absorbed_' RAPFX '_srf'];
%%   gamfld = ['absorbed_' RAPFX '_gamma'];
%%   dlrfld = [RAPFX '_' ISPFX '_dlrf'];
%%   ulrfld = [SPREFIX RAPFX '_' ISPFX '_ulrf'];
%%   lrfld = [SPREFIX RAPFX '_' ISPFX '_lrf'];
%% % Water-benthos fluxes
%%   qbfld = ['b_' RAPFX '_srf'];
%%   btfld = ['b_' RAPFX '_t'];
%%   qbofld = ['b_' RAPFX '_qbo'];
%% % Turbulent and net surface fluxes
%%   qlhfld = [TURPFX '_latent_flux'];
%%   qshfld = [TURPFX '_sensible_flux'];
%%   qrhfld = [TURPFX '_rain_flux'];
%%   adj_qlhfld = ['adj_' qlhfld];
%%   qradfld = [SPREFIX RAPFX '_' ISPFX '_arf'];
%%   qradtfld = [qradfld '_term'];
%%   qturfld = [TURPFX '_turbulent_flux'];
%%   qturtfld = [qturfld '_term'];
%%
%%   q0fld = [TURPFX '_net_flux'];
%%   q0lpfld = [q0fld Q0_LOWPASS];
%%   qtfld = [q0fld '_term'];
%%   qlh30fld = [TUR30PFX '_latent_flux'];
%%   qsh30fld = [TUR30PFX '_sensible_flux'];
%%   qrh30fld = [TUR30PFX '_rain_flux'];
%%   qtur30fld = [TUR30PFX '_turbulent_flux'];
%%   q030fld = [TUR30PFX '_net_flux'];
%%   qt30fld = [q030fld '_term'];
%%   q030lpfld = [q030fld Q0_LOWPASS];
%%
%%   sqradfld = [SPREFIX RAPFX '_' ISPFX '_rf'];
%%   sq0fld = ['simple_' q0fld];
%%   sqtfld = ['simple_' qtfld];
%%   sq030fld = ['simple_' q030fld];
%%   sqt30fld = ['simple_' qt30fld];
%%
%%   climsrfld = [CLIMPFX '_srf'];
%%   climasrfld = ['absorbed_' CLIMPFX '_srf'];
%%   climlrfld = [CLIMPFX '_lrf'];
%%   climevapfld = [CLIMPFX '_evap'];
%%   climqlhfld = [CLIMPFX '_latent_flux'];
%%   climqshfld = [CLIMPFX '_sensible_flux'];
%%   climq0fld = [CLIMPFX '_net_flux'];
%%   climqtfld = [climq0fld '_term'];
%%
%%   raevapfld = [RAPFX '_evap'];
%%   raqlhfld = [RAPFX '_latent_flux'];
%%   raqshfld = [RAPFX '_sensible_flux'];
%%   raq0fld = [RAPFX '_net_flux'];
%%   raqtfld = [raq0fld '_term'];
%%
%%   whfld = [WAVEPFX '_sigwavehgt'];
%%   wpfld = [WAVEPFX '_peakwaveper'];
%%   wdfld = [WAVEPFX '_peakwavedir'];
%%   ufld = [KMPFX '_u'];
%%   vfld = [KMPFX '_v'];
%%   Tfld = [KMPFX '_seatemp_field'];
%%   kmtfld = [KMPFX '_seatemp'];
%%   kmtxfld = [kmtfld '_x'];
%%   kmtyfld = [kmtfld '_y'];
%%   kmtlfld = [kmtfld '_l'];
%%   kmtxsfld = [kmtfld '_xshore'];
%%   kmtlsfld = [kmtfld '_lshore'];
%%   hufld = ['hourly_' KMPFX '_u'];
%%   hvfld = ['hourly_' KMPFX '_v'];
%%   hkmtfld = ['hourly_' kmtfld];
%%   hkmtxfld = ['hourly_' kmtxfld];
%%   hkmtyfld = ['hourly_' kmtyfld];
%%   hkmtlfld = ['hourly_' kmtlfld];
%%   hkmtxsfld = ['hourly_' kmtxsfld];
%%   hkmtlsfld = ['hourly_' kmtlsfld];
%%   netufld = [TIDEPFX '_' KMPFX '_u'];
%%   netvfld = [TIDEPFX '_' KMPFX '_v'];
%%   ssufld = [STOKESPFX '_u'];
%%   ssvfld = [STOKESPFX '_v'];
%%   sssfld = [STOKESPFX '_speed'];
%%   ssdfld = [STOKESPFX '_dir'];
%%   ssxsfld = [STOKESPFX '_xshore'];
%%   sslsfld = [STOKESPFX '_lshore'];
%%   qeufld = [QEPFX '_u'];
%%   qevfld = [QEPFX '_v'];
%%   qesfld = [QEPFX '_speed'];
%%   qedfld = [QEPFX '_dir'];
%%   qexsfld = [QEPFX '_xshore'];
%%   qelsfld = [QEPFX '_lshore'];
%%   udTfld = [QEPFX '_advected_heat'];
%%   rawudTfld = ['raw_' udTfld];
%%   qtAdvfld = [TURPFX '_' QEPFX '_qtadv'];
%%   model_K_theta = 20;
%%   kd2Tfld = [KMPFX '_diffused_heat'];
%%   rawkd2Tfld = ['raw_' kd2Tfld];
%%   dTfld = [TURPFX '_' QEPFX '_dt'];
%%   dTffld = [dTfld '_flux'];
%%   dTflpfld = [dTffld Q0_LOWPASS];
%%   qbotfld = [qbofld '_term'];
%%   bq0fld = ['b_' q0fld];
%%   bq0tfld = [bq0fld '_term'];
%%   bq0lpfld = [bq0fld Q0_LOWPASS];
%%   bdTfld = ['b_' dTfld];
%%   bdTffld = [bdTfld '_flux'];
%%   bdTflpfld = [bdTffld Q0_LOWPASS];
%%   hcfactor = [HCPFX '_termFactor'];
%%   hcdTdt = [HCPFX '_dTdt'];
%%   hcdTdtf = [hcdTdt '_flux'];
%%   hcdTdthc = [hcdTdt 'hc'];
%%   hcdTdthcf = [hcdTdthc '_flux'];
%%
%% % Sub-grid-scale heat diffusion (calculated as residual)
%%   sgskd2Tfld = [HCPFX '_sgs_diffused_heat'];
%%   sgskfld = [HCPFX '_sgs_K_theta'];
%%   sgsdTdt = [HCPFX '_sgs_final_budget'];
%%
%% % Climatological total heat budget terms - for comparison 
%%   climqtAdvfld = [CLIMPFX '_' QEPFX '_qtadv'];
%%   climdTfld = [CLIMPFX '_' QEPFX '_dt'];
%%   climdTffld = [climdTfld '_flux'];
%%   climbq0fld = ['b_' climq0fld];
%%   climbq0tfld = [climbq0fld '_term'];
%%   climbdTfld = ['b_' climdTfld];
%%   climbdTffld = [climbdTfld '_flux'];
%%   climhcfactor = [CLIMHCPFX '_termFactor'];
%%   climhcdTdt = [CLIMHCPFX '_dTdt'];








     case 'sanf1',	stn.(orifld)=90;  % Arbitrarily chosen 2011 Apr 22



  % %%%%???DEBUG
  % ignore_benthos = true;




  % Bulk benthic heat transfer coefficient (Cbd^2) for flow over sea-floor
  % %cbds = {0,3.8e-5,2.9e-4,3.8e-4,8.0e-4,16.0e-4,24.0e-4};
  % cbds = {8.0e-4};
  % To coincide with estimate CD~0.017 in Davis and Monismith (2011)
  cbds = {2.9e-4};




%%%%???DEBUG
  ignore_benthos = true;





%%%%??? DEBUG
    % frac = 0.10;
%%%%??? DEBUG



        % stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'-');
%%%%??? DEBUG
        stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'+');
%%%%??? DEBUG



  % If caller is interested in insolation absorbed by sea bed. (Assumes all
  % NIR/UVR radiation is absorbed in water column: only PAR reaches bottom.)
  if ( exist('qbfld','var') && ~isempty(qbfld) )
    stn.(qbfld).date = dts;
    stn.(qbfld).data = qsw .* PAR_PER_INSOL .* tau .* (1 - Ab);
    stn.(qbfld).data(badix) = 0;
  end;





CANONICAL Kd (from MLRF1 and FWYF1)
    %    [0.050,0.250,137] ...





   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true,false];

    kds = { ...
        0.7 ...
        0.65 ...
        [0.500,0.900,201,182.6225] ...
        [0.500,0.850,201,182.6225] ...
        [0.500,0.800,201,182.6225] ...
        [0.550,0.850,201,182.6225] ...
        [0.550,0.800,201,182.6225] ...
        [0.550,0.750,201,182.6225] ...
          };
    % stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    advfacs = { ...
        0 ...
        [0.0,0.5, 90] ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        0 ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
        [ 0,20,183] ...
           };





    kds = { ...
        [0.600,0.850,320] ...
        0.7 ...
        [0.550,0.850,274] ...
        [0.550,0.850,228] ...
        [0.550,0.900,320] ...
        [0.550,0.900,274] ...
        ...
        [0.500,0.900,201,182.6225] ...
        [0.550,0.850,201,182.6225] ...
        [0.650,0.750,201,182.6225] ...
        [0.600,0.800,201,182.6225] ...
          };

   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];

    % stn = read_usf_kd_par(stn);
    % stn.hourly_amodis_kd_par = interp_ts(stn.amodis_kd_par);
    % [dsrix,ig] = intersect_dates(stn.(dsrfld).date,stn.hourly_amodis_kd_par.date);
    % stn.(dsrfld) = subset_ts(stn.(dsrfld),dsrix);
    %     stn.hourly_amodis_kd_par ...

    % U30dfld = [Ufld '_30_d_lp']
    % stn = verify_variable(stn,U30dfld);
    %     { stn.(U30dfld),@(U30d)(0.69+(U30d.*0.009)) } ...
    %    [0.600,0.750,320] ...
    kds = { ...
        [0.600,0.850,320] ...
        [0.600,0.850,274] ...
        [0.600,0.850,228] ...
        [0.600,0.850,183] ...
        0.7 ...
        [0.550,0.850,274] ...
        [0.550,0.850,228] ...
        [0.550,0.850,183] ...
        0.725 ...
        [0.550,0.900,320] ...
        [0.550,0.900,274] ...
        [0.550,0.900,228] ...
        [0.550,0.900,183] ...
        [0.550,0.900, 31,182.6225] ...
          };
    % stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    advfacs = { ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };








LONF1
    kds = { ...
        [0.600,0.750,274] ...
        [0.600,0.750,  0] ...
        [0.600,0.750,320] ...
        ...
        [0.600,0.800,274] ...
        [0.600,0.800,  0] ...
        [0.600,0.800,320] ...
        [0.600,0.850,274] ...
        [0.600,0.850,  0] ...
        [0.600,0.850,320] ...
          };




   case 'FWYF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    %    [0.050,0.250,137] ...
    kds = { ...
        [0.035,0.265,110] ...
          };
    advfacs = { ...
        [0.0,0.5, 90] ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };

   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
%        [0.050,0.250,137] ...
%        [0.070,0.250, 91] ...
    kds = { ...
        [0.050,0.250,137] ...
        ...
        [0.050,0.250, 91] ...
        [0.050,0.250,100] ...
        [0.050,0.250,110] ...
        ...
        [0.040,0.270,137] ...
        ...
        [0.040,0.270, 91] ...
        [0.040,0.270,100] ...
        [0.040,0.270,110] ...
          };
    advfacs = { ...
        1 ...
              };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };





   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    % stn = read_usf_kd_par(stn);
    % stn.hourly_amodis_kd_par = interp_ts(stn.amodis_kd_par);
    % [dsrix,ig] = intersect_dates(stn.(dsrfld).date,stn.hourly_amodis_kd_par.date);
    % stn.(dsrfld) = subset_ts(stn.(dsrfld),dsrix);
    %     stn.hourly_amodis_kd_par ...
    kds = { ...
        [0.550,0.800,320] ...
        [0.550,0.800,  0] ...
        [0.550,0.800, 45] ...
        [0.550,0.800, 91] ...
        [0.550,0.800,137] ...
        [0.550,0.800,183] ...
        [0.550,0.800,228] ...
        [0.550,0.800,274] ...
          };
    % stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    advfacs = { ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };






    kds = { ...
        [0.500,0.800,274,183] ...
        [0.550,0.800,274,183] ...
        [0.600,0.750,274,183] ...
        [0.550,0.750,274] ...
        [0.600,0.800,274] ...
        [0.550,0.800,274] ...
        [0.600,0.750,274] ...
          };




   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true,false];
    kds = { ...
        [0.450,0.650,274,183] ...
        [0.550,0.750,274,183] ...
        [0.600,0.800,274,183] ...
        [0.450,0.650,274] ...
        [0.550,0.750,274] ...
        [0.600,0.800,274] ...
        [0.550,0.800,274] ...
        [0.600,0.750,274] ...
          };
    % stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    advfacs = { ...
        0 ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };




  %% Brian Haus Wave experiment site #7 (bhwv7) - on reef slope

  % ix = 1:length(stns.bhwv7.triaxys_sigwaveper.data);
  % ix = find(0.1<stns.bhwv7.triaxys_sigwavehgt.data & stns.bhwv7.triaxys_peakwaveper.data<=15);



%%%%???DEBUG
ix=[];
%%%%???DEBUG





   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false,true];
    kds = { ...
        [0.550,0.750,  0,183] ...
        [0.550,0.750, 45,183] ...
        [0.550,0.750, 91,183] ...
        [0.550,0.750,137,183] ...
        [0.550,0.750,183,183] ...
        [0.550,0.750,228,183] ...
        [0.550,0.750,274,183] ...
        [0.550,0.750,320,183] ...
        [0.600,0.800,274] ...
          };
    % stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    advfacs = { ...
        0 ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };




   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
%%%%???DEBUG
%%%%???DEBUG
    kds = { ...
        [0.550,0.750,  0] ...
        [0.550,0.750, 45] ...
        [0.550,0.750,228] ...
        [0.550,0.750,274] ...
        [0.600,0.800,  0] ...
        [0.600,0.800, 45] ...
        [0.600,0.800,228] ...
        [0.600,0.800,274] ...
        [0.550,0.750,  0,183] ...
        [0.550,0.750, 45,183] ...
        [0.550,0.750, 91,183] ...
        [0.550,0.750,137,183] ...
        [0.550,0.750,183,183] ...
        [0.550,0.750,228,183] ...
        [0.550,0.750,274,183] ...
        [0.550,0.750,320,183] ...
          };
    % stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    advfacs = { ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };






        %%%%???DEBUG - TRY ELIMINATING BENTHIC EXCHANGE SUB-MODEL        
        disp('** Ignoring absorption and benthic exchange **');
        stn.(bq0fld) = stn.(sq0fld);
        stn.commentstr = [stn.commentstr,' (NO Q_b) '];
        %%%%???DEBUG







   case 'LONF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
%%%%???DEBUG
%%%%???DEBUG
    kds = { ...
        [0.500,0.700,  0] ...
        [0.500,0.700, 45] ...
        [0.500,0.700, 91] ...
        [0.500,0.700,137] ...
        [0.500,0.700,183] ...
        [0.500,0.700,228] ...
        [0.500,0.700,274] ...
        [0.500,0.700,320] ...
        [0.600,0.800,  0] ...
        [0.600,0.800, 45] ...
        [0.600,0.800, 91] ...
        [0.600,0.800,137] ...
        [0.600,0.800,183] ...
        [0.600,0.800,228] ...
        [0.600,0.800,274] ...
        [0.600,0.800,320] ...
        [0.500,0.700,  0,183] 
          };
    % stn.opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    advfacs = { ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    kths = { ...
        0 ...
        10 ...
        20 ...
        { stn.(Wfld),@(W)(min(40,((W./35).^2).*40)) } ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };





   case 'LONF1',
    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.20,0.37,228], };
    % kths = { [0 ,5 ,183], };

    % *MLRF* NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];
    kds = { ...
        [0.050,0.250,137] ...
        [0.100,0.300,137] ...
        [0.100,0.300,183] ...
        [0.150,0.375,137] ...
        [0.150,0.375,183] ...
          };
    advfacs = { ...
        0 ...
        1 ...
              };
    % LONF1 diffusion - related to tidal currents?
    %[0 ,5 ,183],...
    kths = { ...
        [0, 5, 45] ...
        [0, 5, 183] ...
           };





   case 'SANF1',
    % Number of points to use in finite-difference templates for gradients
    npts = 3;

    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.20,0.37,228], };
    % kths = { 0 };

    % *MLRF* NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false,true];
    kds = { ...
        [0.050,0.250,137] ...
          };
    advfacs = { ...
        1 ...
              };
    kths = { ...
        0 ...
        [0, 5, 45] ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };





   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.050,0.250,137] ...
          };
    advfacs = { ...
        1 ...
              };
    % % Ignore unphysical SST gradients - see SFP TSG data analysis
    % stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    % advfacs = { ...
    %     1 ...
    %     [0.0,0.5, 90] ...
    %           };
    % kths = { ...
    %     0 ...
    %     [0,20, 45] ...
    %     { stn.(Wfld),@(W)(min(40,((W./35).^2).*40)) } ...
    %     { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
    %     { stn.(Wlpfld),@(Wlp)(min(40,((Wlp./35).^2).*40)) } ...
    %     { stn.(Wlpfld),@(Wlp)(min(20,((Wlp./35).^2).*20)) } ...
    %        };
    kths = { ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
           };





FWYF1:

    % % Enhanced albedo - justification?
    % stn.opts.albedo_increment = get_opt(stn.opts,'albedo_increment',1.5);


    % % Ignore unphysical SST gradients - see SFP TSG data analysis
    % stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    % % Ignore alongshore component of heat advection
    % stn.opts.add_alongshore_advection = get_opt(stn.opts,'add_alongshore_advection',false);
    advfacs = { ...
        [0.0,0.5, 90] ...
        { stn.(Wlpfld),@(Wlp)(min(0.5,((Wlp./35)).*0.5)) } ...
        { stn.(Wlpfld),@(Wlp)(min(0.5,((Wlp./35).^2).*0.5)) } ...
              };






        stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'-');
%%%%???DEBUG - TRY ELIMINATING BENTHIC EXCHANGE
        % stn.(bq0fld) = stn.(sq0fld);
        % stn.commentstr = [stn.commentstr,' (NO Q_b) '];
%%%%???DEBUG





%%%%??? DEBUG
% if ( strcmp(stnm,'fwyf1') )
%   stn.commentstr = [stn.commentstr,' (No Fall Adv)'];
%   stn.(hkmtxfld).data(get_month(stn.(hkmtxfld).date)>9) = 0;
%   stn.(hkmtyfld).data(get_month(stn.(hkmtyfld).date)>9) = 0;
% end;
%%%%??? DEBUG







%%%%???DEBUG
    % % Ignore unphysical SST gradients - see SFP TSG data analysis
    % stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
%%%%???DEBUG



   case 'FWYF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true,false];
    % % Enhanced albedo - justification?
    % stn.opts.albedo_increment = get_opt(stn.opts,'albedo_increment',1.5);
    kds = { ...
        [0.050,0.250,137] ...
          };
    % Ignore unphysical SST gradients - see SFP TSG data analysis
    stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    % % Ignore alongshore component of heat advection
    % stn.opts.add_alongshore_advection = get_opt(stn.opts,'add_alongshore_advection',false);
    advfacs = { ...
        0 ...
              };
    kths = { ...
        0 ...
        [0,20, 45] ...
        { stn.(Wfld),@(W)(min(20,(W./35).*20)) } ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
        { stn.(Wlpfld),@(Wlp)(min(40,((Wlp./35).^2).*40)) } ...
        { stn.(Wlpfld),@(Wlp)(min(20,((Wlp./35).^2).*20)) } ...
           };









   case 'FWYF1',
...
    kths = { ...
        [0, 4, 45] ...


   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];
    kds = { ...
        [0.050,0.250,137] ...
          };
    % Ignore unphysical SST gradients - see SFP TSG data analysis
    stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    advfacs = { ...
        1 ...
              };
    kths = { ...
        [0,20, 45] ...
        { stn.(Wfld),@(W)(min(20,(W./35).*20)) } ...
        { stn.(Wfld),@(W)(min(20,((W./35).^2).*20)) } ...
        { stn.(Wfld),@(W)(min(10,(W./35).*10)) } ...
        { stn.(Wfld),@(W)(min(10,((W./35).^2).*10)) } ...
           };





%{
  if ( ~isfield(stn,bathyfld) )
    stn = get_ngdc_bathy_station(stn);
  end;
  if ( ~isfield(stn,slopefld) )
    stn = station_ngdc_offshore_slope(stn);
  end;
  if ( ~isfield(stn,bathorifld) )
    stn = station_optimal_isobath_orientation(stn);
  end;

  warning('off','Ecoforecasts:mergedNonTS');
  stn = load_all_ndbc_data(stn);
%}






FWYF1:
    % kths = { ...
    %     0 ...
    %     [0, 4, 45] ...
    %     2 ...





   case 'SMKF1',
    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.045,0.375,45], };
    % kths = { 0 };

    % % ERAI met, NDBC air_t, higher Ppen - amplitude problem but good annual...
    % kds = { [0.045,0.375, 91], }

    % % ERAI met, NDBC air_t, new and improved higher Ppen
    % kds = { [0.060,0.400, 45], };
    % kths = { 0 };

    % % ERAI met, NDBC air_t, new Ppen - Kd optimized to close interannual budget
    % % Q0 w/Warm Layer also matches OAFlux annual amplitude! But 2005 TOO HOT
    % doWarms = [true];
    % % Satellite chlor_a from USF suggests a SEMI-annual Kd cycle
    % % kds = { [0.060,0.250, 45, 183], };
    % % kds = { [0.100,0.250, 45, 183], };
    % kds = { ...
    %     % [0.060,0.250, 45, 183],...
    %     % [0.100,0.250, 45, 183],...
    %     % [0.060,0.250, 45],...
    %     [0.100,0.250, 45],...
    %     [0.100,0.300, 45],...
    %       };
    % kths = { 0 };

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true,false];
    kds = { ...
        [0.050,0.250,137] ...
        [0.100,0.300,137] ...
        [0.100,0.300,183] ...
        [0.150,0.375,137] ...
        [0.150,0.375,183] ...
          };

    % Ignore unphysical SST gradients - see SFP TSG data analysis
    stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    % % Ignore alongshore component of heat advection?
    % stn.opts.add_alongshore_advection = get_opt(stn.opts,'add_alongshore_advection',false);
    advfacs = { ...
        0 ...
        [0, 1, 183] ...
        1 ...
              };
    kths = { ...
        0 ...
        [0, 5, 45] ...
        [0,10, 45] ...
           };

%%%%??? DEBUG
    doWarms = [true];
    kds = { ...
        [0.200,0.400,137] ...
        [0.300,0.500,137] ...
        [0.500,0.700,137] ...
          };
    advfacs = { ...
        0 ...
        1 ...
              };
    kths = { ...
        0 ...
        2 ...
        [0, 5, 45] ...
           };
%%%%??? DEBUG











   case 'FWYF1',
    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.035,0.220, 76], };
    % kths = { 0 };

    % % ERAI met, ERAI air_t, new Ppen - warm layer adjustment matches
    % % climatology closely, but grossly overestimates interannual budget
    % doWarms = [true];
    % kds = { ...
    %     [0.035,0.220, 91], ...
    %       };
    % kths = { ...
    %     [0,20, 45] ...
    %        };

    % % ERAI met, NDBC air_t, new Ppen - warm layer adjustment matches
    % % climatology closely, but grossly overestimates interannual budget
    % % doWarms = [true];
    % doWarms = [false];
    % kds = { ...
    %     [0.025,0.220, 91], ...
    %       };
    % kths = { ...
    %     [0,10, 45] ...
    %        };

    % % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves
    % doWarms = [true];
    % kds = { ...
    %     [0.050,0.350,137] ...
    %     [0.050,0.350,183] ...
    %     [0.050,0.350,228] ...
    %     0.2 ...
    %     [0.100,0.350,137] ...
    %     [0.100,0.350,183] ...
    %     [0.100,0.350,228] ...
    %     0.225 ...
    %       };
    % kths = { ...
    %     0 ...
    %     [0,10, 45] ...
    %     [0,20, 45] ...
    %        };

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];

    % % Enhanced albedo - justification?
    % stn.opts.albedo_increment = get_opt(stn.opts,'albedo_increment',1.5);
    kds = { ...
        [0.050,0.250,137] ...
          };

    % Ignore unphysical SST gradients - see SFP TSG data analysis
    stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    % Ignore alongshore component of heat advection
    stn.opts.add_alongshore_advection = get_opt(stn.opts,'add_alongshore_advection',false);
    advfacs = { ...
        1 ...
              };
    kths = { ...
        0 ...
        [0, 4, 45] ...
        2 ...
           };










  if ( isfield(stn,'opts') )
    if ( isfield(stn.opts,'kd') )
      appendtitlename([' K_d=' num2str(stn.opts.kd) ]);
    end;
    if ( isfield(stn.opts,'K_theta') )
      appendtitlename([' K\theta=' num2str(stn.opts.K_theta)]);
    end;
  end;





num2str(kth,'%g,')



   case 'MLRF1',
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];

    % % Enhanced albedo - justification?
    % stn.opts.albedo_increment = get_opt(stn.opts,'albedo_increment',1.5);
    kds = { ...
        [0.050,0.250,137] ...
          };

    % Ignore unphysical SST gradients - see SFP TSG data analysis
    stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    advfacs = { ...
        1 ...
              };
    kths = { ...
        % 0 ...
        % [0,10,  0] ...
        % [0,10, 45] ...
        % [0,10, 91] ...
        % [0,20,  0] ...
        [0,20, 45] ...
        % [0,20, 91] ...
        {stn.(Wfld),@(dat)(min(20,dat.*(20/35)))}
           };



'SMKF1'
%%%%??? DEBUG
    doWarms = [false,true];
    kds = { ...
        [0.100,0.300, 46] ...
        [0.100,0.300,137] ...
        [0.100,0.300,183] ...
        [0.150,0.375, 46] ...
        [0.150,0.375, 92] ...
        [0.150,0.375,137] ...
        [0.150,0.375,183] ...
        [0.100,0.300, 92] ...
          };
    advfacs = { ...
        0 ...
        1 ...
              };
    kths = { ...
        0 ...
        [0, 5, 45] ...
           };
%%%%??? DEBUG





    kds = { ...
        [0.100,0.300,137] ...
        [0.150,0.375,137] ...
        [0.150,0.375,183] ...
        [0.100,0.300, 92] ...
          };




          fmg;
          climsq = squeeze(stn.optim.climsq(:,doWarmix,advix,kthix));
          climsq_err = squeeze(stn.optim.climsq_err(:,doWarmix,advix,kthix));
          % plot_ts(stn.optim.climt,climsq,ts_op(climsq,climsq_err,'-'),ts_op(climsq,climsq_err,'+'));
          climsq_minus_err.date = climsq_err.date;
          climsq_minus_err.data = climsq.data - cumsum(climsq_err.data);
          climsq_plus_err.date = climsq_err.date;
          climsq_plus_err.data = climsq.data + cumsum(climsq_err.data);
          plot_ts(stn.optim.climt,climsq,climsq_minus_err,climsq_plus_err);






   case 'MLRF1',
    % % ERAI met, air_t, old Ppen
    % kds = { [0.045,0.375,45], };
    % kths = { 0 };

    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.035,0.190, 55], };
    % kths = { 0 };

    % % ERAI met, NDBC air_t, new Ppen - Kd optimized to close interannual budget
    % doWarms = [true];
    % kds = { ...
    %     [0.050,0.250,137], ...
    %       };
    % kths = { ...
    %     0, ...
    %     % [-10, -5,137] ...
    %     % [-12, -4,137] ...
    %     % [-15,  0,137] ...
    %        };

    % % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves
    % doWarms = [true];
    % kds = { ...
    %     [0.050,0.350,137] ...
    %     [0.050,0.350,183] ...
    %     [0.050,0.350,228] ...
    %     0.2 ...
    %     [0.100,0.350,137] ...
    %     [0.100,0.350,183] ...
    %     [0.100,0.350,228] ...
    %     0.225 ...
    %       };
    % kths = { ...
    %     0 ...
    %     [0,10, 45] ...
    %     [0,20, 45] ...
    %        };

    % % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients
    % doWarms = [true];
    % kds = { ...
    %     [0.050,0.350,137] ...
    %     [0.100,0.300, 45] ...
    %     [0.100,0.300,137] ...
    %     [0.100,0.300,183] ...
    %     [0.050,0.250,137] ...
    %       };
    % advfacs = { ...
    %     0,...
    %     1,...
    %           };
    % kths = { ...
    %     0 ...
    %     [0,10, 45] ...
    %     [0,20, 45] ...
    %        };

    % % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk ULRF
    % doWarms = [true];
    % kds = { ...
    %     [0.050,0.350,137] ...
    %     [0.050,0.350,160] ...
    %     [0.100,0.300,137] ...
    %     [0.100,0.300,160] ...
    %     [0.050,0.250,137] ...
    %     [0.050,0.250,160] ...
    %     [0.050,0.200,137] ...
    %     [0.050,0.200,160] ...
    %       };
    % advfacs = { ...
    %     0,...
    %     1,...
    %           };
    % kths = { ...
    %     0 ...
    %     [0,10, 45] ...
    %     [0,20, 45] ...
    %        };

    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [false];

    % % Enhanced albedo - justification?
    % stn.opts.albedo_increment = get_opt(stn.opts,'albedo_increment',1.5);
    kds = { ...
        [0.050,0.250,137] ...
          };

    % Ignore unphysical SST gradients - see SFP TSG data analysis
    stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    advfacs = { ...
        1 ...
              };

    kths = { ...
        [0, 5, 45] ...
           };

%%%%DEBUG???
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves, AVHRR gradients, bulk USRF, bulk ULRF
    doWarms = [true];

    % % Enhanced albedo - justification?
    % stn.opts.albedo_increment = get_opt(stn.opts,'albedo_increment',1.5);
    kds = { ...
        [0.050,0.250,137] ...
        [0.100,0.300,137] ...
        [0.050,0.300,137] ...
          };

    % Ignore unphysical SST gradients - see SFP TSG data analysis
    stn.opts.max_sst_gradient = get_opt(stn.opts,'max_sst_gradient',7e-4);
    advfacs = { ...
        1 ...
              };
    kths = { ...
        0 ...
        [0, 5, 45] ...
        [0,10, 45] ...
           };
%%%%DEBUG???

















%%%%???DEBUG
%a=0; b=1.0;
%%%%???DEBUG




%%%%DEBUG???
    stn.opts.albedo_increment = 0;
%%%%DEBUG???





  % stn.erai_dsrf_adj.date = stn.erai_dsrf.date;
  %   stn.erai_dsrf_adj.date(ix,1) = stn.erai_dsrf.date(ix);
  %   stn.erai_dsrf_adj.date(ix,1) = stn.erai_dsrf.date(ix);
  %   stn.erai_dsrf_adj.date(ix,1) = stn.erai_dsrf.date(ix);
  stn.erai_dsrf_adj.date = stn.erai_dsrf.date;



    stn.erai_dlrf_adj.date(ix,1) = stn.erai_dlrf.date(ix);





  stn.erai_dsrf_adj = struct('date',[],'data',[]);


  stn.erai_dlrf_adj = struct('date',[],'data',[]);





%%%%???DEBUG
    adjust_reanalysis = false;






    % Compare simple cumulative sums of daily climatology vs. our (last) estimate
    err_dt = 30*24;

    [dix,cix,six,bix,aix,kix,hix] = ...
        intersect_all_dates([],stn.(dsffld).date,stn.(climq0fld).date,stn.(sq0fld).date,stn.(bq0fld).date,stn.(qtAdvffld).date,stn.(bdTffld).date,stn.(hcdTdtf).date);

    fmg;
    lh=[];
    lh(end+1)=plot(stn.(dsffld).date(dix(1):end),nancumsum(stn.(dsffld).data(dix(1):end)),'k');
    lh(end+1)=plot(stn.(climq0fld).date(cix(1):end),nancumsum(stn.(climq0fld).data(cix(1):end).*24),'c');

    cs.date = stn.(sq0fld).date(six(1):end);
    cs.data = nancumsum(stn.(sq0fld).data(six(1):end));
    lh(end+1)=plot(cs.date,cs.data,'r-.');
    [csix,serrix] = intersect_dates(cs.date,stn.([sq0fld,'_err']).date);
    csp = cs.data(csix) + nancumsum(stn.([sq0fld,'_err']).data(serrix));
    plot(stn.([sq0fld,'_err']).date(serrix(1:err_dt:end)),csp(1:err_dt:end),'r.');
    csn = cs.data(csix) - nancumsum(stn.([sq0fld,'_err']).data(serrix));
    plot(stn.([sq0fld,'_err']).date(serrix(1:err_dt:end)),csn(1:err_dt:end),'r.');

    bcs.date = stn.(bq0fld).date(bix(1):end);
    bcs.data = nancumsum(stn.(bq0fld).data(bix(1):end));
    lh(end+1)=plot(stn.(bq0fld).date(bix(1):end),nancumsum(stn.(bq0fld).data(bix(1):end)),'m:');
    [bcsix,qerrix] = intersect_dates(bcs.date,stn.([q0fld,'_err']).date);
    bcsp = bcs.data(bcsix) + nancumsum(stn.([q0fld,'_err']).data(qerrix));
    plot(stn.([q0fld,'_err']).date(qerrix(1:err_dt:end)),bcsp(1:err_dt:end),'m.');
    bcsn = bcs.data(bcsix) - nancumsum(stn.([q0fld,'_err']).data(qerrix));
    plot(stn.([q0fld,'_err']).date(qerrix(1:err_dt:end)),bcsn(1:err_dt:end),'m.');

    lh(end+1)=plot(stn.(qtAdvffld).date(aix(1):end),nancumsum(stn.(qtAdvffld).data(aix(1):end)),'y--');
    lh(end+1)=plot(stn.(bdTffld).date(kix(1):end),nancumsum(stn.(bdTffld).data(kix(1):end)),'y:');
    lh(end+1)=plot(stn.(hcdTdtf).date(hix(1):end),nancumsum(stn.(hcdTdtf).data(hix(1):end)),'b--');
    legend(lh,'Actual \rhoC_ph\partial_tT_s','OAFlux/ISCCP Q_0',...
           'G&M Q_0 \pm\sigmaQ_0','G&M Q_0(\gamma)+Q_b','G&M Q_0(\gamma)+Q_b + \rhoC_phu^.\nablaT',...
           'G&M Q_0(\gamma)+Q_b + \rhoC_ph(u^.\nablaT + K\nabla^2T)',...
           'G&M \rhoC_ph\partial_tT_s',...
           'Location','Best');
    datetick3('x',20,'keeplimits');
    ylim([-3e6,+3e6]);    ylabel('W/m^2');
    titlename([STNM ' simple cumulative sums: Hourly fluxes K_d=' num2str(opts.kd) ' K\theta=' num2str(opts.K_theta) stn.commentstr ]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-simple-sum-',hcdTdtf,'.tiff']));
    end;










    cs.date = stn.(sq0fld).date(six(1):end);
    cs.data = nancumsum(stn.(sq0fld).data(six(1):end));
    lh(end+1)=plot(cs.date,cs.data,'r-.');
    [csix,serrix] = intersect_dates(cs.date,stn.([sq0fld,'_err']).date);
    csp = nancumsum(stn.(sq0fld).data(six)) + stn.([sq0fld,'_err']).data(serrix);
    plot(stn.([sq0fld,'_err']).date(serrix(1:24:end)),csp(1:24:end),'r^');
    csn = nancumsum(stn.(sq0fld).data(six)) - stn.([sq0fld,'_err']).data(serrix);
    plot(stn.([sq0fld,'_err']).date(serrix(1:24:end)),csn(1:24:end),'rv');

qerrix,
stn.([q0fld,'_err']).date,
    lh(end+1)=plot(stn.(bq0fld).date(bix(1):end),nancumsum(stn.(bq0fld).data(bix(1):end)),'m:');
    csp = nancumsum(stn.(bq0fld).data(six)) + stn.([q0fld,'_err']).data(qerrix);
    plot(stn.([q0fld,'_err']).date(qerrix(1:24:end)),csp(1:24:end),'m^');
    csn = nancumsum(stn.(bq0fld).data(six)) - stn.([q0fld,'_err']).data(qerrix);
    plot(stn.([q0fld,'_err']).date(qerrix(1:24:end)),csn(1:24:end),'mv');





%%%% ??? DEBUG
    advfacs = { ...
        0 ...
              };
    kths = { ...
        0 ...
           };
%%%% ??? DEBUG





      %%%%??? DEBUG
      maxGrad = 1e-3;
      disp(['** Zeroing unphysical gradients (|dSST/di|>',num2str(maxGrad*1e3),'K/km) **']);
      badix = find(abs(stn.(hkmtxfld).data)>maxGrad | abs(stn.(hkmtyfld).data)>maxGrad);
      zeroix = [];
      for ix=badix(:)'; zeroix = [zeroix ix-(24*7)+1:ix+(24*7)-1]; end;
      badix = unique(zeroix);
      for ix=badix(:)';
        allix = 
        stn.(hkmtxfld).data(allix) = 0;
        stn.(hkmtyfld).data(allix) = 0;
      end;
      %%%%??? DEBUG






    % stn = station_bulk_ulr(stn,sfld,ulrfld);
    stn = station_bulk_ulr(stn,sfld,ulrfld,dlrfld);



  % %%%% ??? DEBUG
  %if ( ~reanalysis_shortwave )
  if ( ~isfield(stn,usrfld) )
    disp('** Using reanalysis downward, bulk upward shortwave **');
    stn.commentstr = [stn.commentstr,' (bulk USRF) '];
    stn = station_bulk_albedo(stn,albfld,Wfld,cfld);
    %%%% ??? DEBUG
    albinc = 1.5;
    disp(['** Albedo + ',num2str(albinc),'%! **']);
    stn.commentstr = [stn.commentstr,' (+',num2str(albinc),'%) '];
    stn.(albfld).data = stn.(albfld).data + (albinc/100);
    %%%% ??? DEBUG
    stn.(usrfld) = ts_op(stn.(dsrfld),stn.(albfld),'.*');
    stn.(usrfld).data = -stn.(usrfld).data;
    stn.(srfld) = ts_op(stn.(dsrfld),stn.(usrfld),'+');
  end;
  % %%%% ??? DEBUG




  % %%%% ??? DEBUG

  disp('** Using reanalysis longwave fluxes **');
  stn.(dlrfld) = stn.([RAPFX '_dlrf']);
  stn.(ulrfld) = stn.([RAPFX '_ulrf']);
  stn.commentstr = [stn.commentstr,' (',upper(RAPFX),' LRF) '];

  % disp('** Using reanalysis downward, bulk upward longwave **');
  % stn.(dlrfld) = stn.([RAPFX '_dlrf']);
  % stn = station_bulk_ulr(stn,sfld,ulrfld);
  % stn.commentstr = [stn.commentstr,' (',upper(RAPFX),' DLRF) '];


  % disp('** Using bulk longwave fluxes **');
  % stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);

  % disp('** Using (simple) bulk longwave fluxes **');
  % stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);
  % stn = station_bulk_ulr(stn,sfld,ulrfld);
  % stn.commentstr = [stn.commentstr,' (simple ULRF) '];

  % disp('** Using reanalysis downward, SST-bulk upward longwave **');
  % stn.(dlrfld) = stn.([RAPFX '_dlrf']);
  % if ( ~isfield(stn,'avhrr_weekly_sst') ); stn = get_avhrr_weekly_field(stn,true); end;
  % if ( ~isfield(stn,'hourly_avhrr_weekly_sst') ); stn.hourly_avhrr_weekly_sst = interp_ts(stn.avhrr_weekly_sst); end;
  % stn = station_bulk_ulr(stn,'hourly_avhrr_weekly_sst',ulrfld);
  % stn.commentstr = [stn.commentstr,' (',upper(RAPFX),' DLRF, AVHRR ULRF) '];

  %%%% ??? DEBUG








%ix=find(ismember(get_jday(tsg.date),[196:227]) & ds<10 & (abs(sind(mlrf1.isobath_orientation-tsg.hdg-90))<0.25));

%%%%DEBUG???
ix=find(ismember(get_jday(tsg.date),[196:227]) & ds<10);
fmg;
[cs,ch]=contour(mlrf1.ngdc_92m_bathy.lon,mlrf1.ngdc_92m_bathy.lat,mlrf1.ngdc_92m_bathy.field,[0:-2:-10,-20,-30]);
clabel(cs,ch,[-2 -6 -10 -20 -30]);
axis([mlrf1.lon-0.11,mlrf1.lon+0.11,mlrf1.lat-0.11,mlrf1.lat+0.11]); axis square;
axes; hold on; set(gca,'Color','none'); scatter(tsg.lon(ix),tsg.lat(ix),[],get_year(tsg.date(ix)));
set_pcolor_cursor; colorbar('East'); plot(mlrf1.lon,mlrf1.lat,'kp'); linkaxes;
axis([mlrf1.lon-0.11,mlrf1.lon+0.11,mlrf1.lat-0.11,mlrf1.lat+0.11]); axis square;
titlename('TSG Datapoints by Year');
%%%%DEBUG???






    grpplot_ts(stn.(climq0fld),@get_month,@nanmean,[],'rs');
    [cum,tid]=grp_ts(stn.(climq0fld).data,stn.(climq0fld).date,@get_month,@nanmean,[]); plot(1:12,cum,'ks','MarkerSize',6); 





    if ( ~isfield(stn,hkmtxsfld) )
      % Cross- and long-shore components of sea temperature gradient
      stn = station_reorient_vectors(stn,bathorifld,hkmtxfld,hkmtyfld);
    end;

    hkmtxslpfld = [hkmtxsfld,'_90_d_lp']; stn=verify_variable(stn,hkmtxslpfld);
    hkmtlslpfld = [hkmtlsfld,'_90_d_lp']; stn=verify_variable(stn,hkmtlslpfld);

    % %%%% DEBUG???
    % maxGrad = 6e-4;
    % disp(['** Limiting SST gradients to ',num2str(maxGrad*1e3),'K/km **']);
    % stn.(hkmtxsfld).data(stn.(hkmtxsfld).data>maxGrad)=maxGrad;
    % stn.(hkmtxsfld).data(stn.(hkmtxsfld).data<-maxGrad)=-maxGrad;
    % stn.(hkmtlsfld).data(stn.(hkmtlsfld).data>maxGrad)=maxGrad;
    % stn.(hkmtlsfld).data(stn.(hkmtlsfld).data<-maxGrad)=-maxGrad;
    % %%%% DEBUG???

    % % % Include tidal currents in advection estimate?
    % % stn.(udTxsfld) = ts_op(stn.(hkmtxsfld),stn.(netxsfld),'.*');
    % % Exclude tidal currents from advection estimate!
    % stn.(udTxsfld) = ts_op(stn.(hkmtxsfld),stn.(qexsfld),'.*');
    stn.(udTxsfld) = ts_op(stn.(hkmtxslpfld),stn.(qexsfld),'.*');
    % Convert to units of [K/hr]
    stn.(udTxsfld).data = -(3600*stn.(udTxsfld).data);
    stn = station_heat_flux_term_inverse(stn,udTfxsfld,udTxsfld,sfld,[],mhfld);

    % % % Include tidal currents in advection estimate?
    % % stn.(udTlsfld) = ts_op(stn.(hkmtlsfld),stn.(netlsfld),'.*');
    % % Exclude tidal currents from advection estimate!
    % stn.(udTlsfld) = ts_op(stn.(hkmtlsfld),stn.(qelsfld),'.*');
    stn.(udTlsfld) = ts_op(stn.(hkmtlslpfld),stn.(qelsfld),'.*');
    % Convert to units of [K/hr]
    stn.(udTlsfld).data = -(3600*stn.(udTlsfld).data);
    stn = station_heat_flux_term_inverse(stn,udTflsfld,udTlsfld,sfld,[],mhfld);

    % Include both vector components of model heat advection
    stn.(udTfld) = ts_op(stn.(udTxsfld),stn.(udTlsfld),'+');
    stn = station_heat_flux_term_inverse(stn,udTffld,udTfld,sfld,[],mhfld);








    cs = nancumsum(stn.(sq0fld).data(six));
    plot(stn.([sq0fld,'_err']).date(errix(1:24:end)),cs(1:24:end)+stn.([sq0fld,'_err']).data(errix(1:24:end)),'r^');
    plot(stn.([sq0fld,'_err']).date(errix(1:24:end)),cs(1:24:end)-stn.([sq0fld,'_err']).data(errix(1:24:end)),'rv');




    % Compare one-day averages (sub-sampled) with OAFlux climatology

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% ??? DEBUG DEBUG DEBUG
    %{

if(0)
    [cix,six] = intersect_dates(stn.(climq0fld).date,stn.(sq01dfld).date);

    fmg;
    boxplot_ts(stn.(climq0fld),'month','index',cix,...
               'title',[STNM,' OAFlux/ISCCP Q_0',stn.commentstr]);
    ylim([-1000,1000]); ylabel('W/m^2');

    fmg;
    boxplot_ts(stn.(sq01dfld),'month','index',six,...
               'title',[STNM,' Gramer&Mariano Sea-surface Q_0',stn.commentstr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',sq01dfld,'.tiff']));
    end;


    [cix,six] = intersect_dates(stn.(climsrfld).date,stn.(sr1dfld).date);

    fmg;
    boxplot_ts(stn.(climsrfld),'month','index',cix,...
               'title',[STNM,' OAFlux/ISCCP Q_S_W',stn.commentstr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',climsrfld,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(sr1dfld),'month','index',six,...
               'title',[STNM,' Gramer&Mariano Q_S_W',stn.commentstr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',sr1dfld,'.tiff']));
    end;


    [cix,six] = intersect_dates(stn.(climqlhfld).date,stn.(qlh1dfld).date);

    fmg;
    boxplot_ts(stn.(climqlhfld),'month','index',cix,...
               'title',[STNM,' OAFlux Q_L_H',stn.commentstr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',climqlhfld,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(qlh1dfld),'month','index',six,...
               'title',[STNM,' Gramer&Mariano Q_L_H',stn.commentstr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',qlh1dfld,'.tiff']));
    end;


    fmg;
    sh=boxplot_ts(stn.(sr1dfld),'month','allcolors','b');
    lh=boxplot_ts(stn.(qcool1dfld),'month','allcolors','r');
    ylim([-1000,1000]); ylabel('W/m^2');
    legend([sh(1),lh(1)], 'Q_S_W','Q_L_W+Q_L_H+Q_S_H', 'Location','South');
    titlename([STNM,' Gramer&Mariano Air-Sea fluxes (1d avg)',stn.commentstr]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',sr1dfld,'-vs-',qcool1dfld,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(udTffld1d),'month','allcolors','k',...
               'title',[STNM,' Gramer&Mariano Km-scale Heat advection (1d avg)',stn.commentstr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',udTffld1d,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(kd2Tffld1d),'month','allcolors','k',...
               'title',[STNM,' Gramer&Mariano Km-scale Heat diffusion (1d avg)',stn.commentstr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',kd2Tffld1d,'.tiff']));
    end;

    fmg;
    boxplot_ts(stn.(hcdTdthcf1d),'month','allcolors','k',...
               'title',[STNM,' Gramer&Mariano Horizontal Convection (1d avg)',stn.commentstr]);
    ylim([-1000,1000]); ylabel('W/m^2');
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',hcdTdthcf1d,'.tiff']));
    end;
end;

    fmg;
    lh=boxplot_ts(stn.(sr1dfld),'month','allcol','b');
    rh=boxplot_ts(stn.(climsrfld),'month','allcol','r');
    legend([lh(1),rh(1)],'ERA-Interim','OAFlux/ISCCP');
    titlename([STNM,' Net insolation Q_S_W',stn.commentstr]);
    ylabel('Wm^-^2'); ylim([-1000,1000]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',sr1dfld,'-vs-',climsrfld,'.tiff']));
    end;

    fmg;
    lh=boxplot_ts(stn.(lr1dfld),'month','allcol','b');
    rh=boxplot_ts(stn.(climlrfld),'month','allcol','r');
    legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux/ISCCP');
    titlename([STNM,' Net longwave flux Q_L_W',stn.commentstr]);
    ylabel('Wm^-^2'); ylim([-1000,1000]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',lr1dfld,'-vs-',climlrfld,'.tiff']));
    end;

    fmg;
    lh=boxplot_ts(stn.(qlh1dfld),'month','allcol','b');
    rh=boxplot_ts(stn.(climqlhfld),'month','allcol','r');
    legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux');
    titlename([STNM,' Latent heat flux Q_L_H',stn.commentstr]);
    ylabel('Wm^-^2'); ylim([-1000,1000]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',qlh1dfld,'-vs-',climqlhfld,'.tiff']));
    end;

    fmg;
    lh=boxplot_ts(stn.(qsh1dfld),'month','allcol','b');
    rh=boxplot_ts(stn.(climqshfld),'month','allcol','r');
    legend([lh(1),rh(1)],'Gramer&Mariano','OAFlux');
    titlename([STNM,' Sensible heat flux Q_S_H',stn.commentstr]);
    ylabel('Wm^-^2'); ylim([-1000,1000]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-boxplot-',qsh1dfld,'-vs-',climqshfld,'.tiff']));
    end;

    nd = ts_op(stn.(a1dfld),stn.(s1dfld),'-',@(x)(intersect_dates(x.date,stn.(climafld).date)));
    od = ts_op(stn.(climafld),stn.(climsfld),'-',@(x)(intersect_dates(x.date,stn.(a1dfld).date)));
    fmg; grpplot_ts(nd,@get_week,@nanmedian,0,'b'); grpplot_ts(od,@get_week,@nanmedian,0,'r');
    legend('In situ','OAFlux');
    titlename([STNM,' Air-Sea Temperature Difference']);
    ylabel('K'); ylim([-5,2]);
    if ( doPrint )
      print('-dtiff',fullfile(figspath,[stnm,'-median-air-sea-temperature-',ISPFX,'-vs-',CLIMPFX,'.tiff']));
    end;

    %}
    %%%% ??? DEBUG DEBUG DEBUG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












        legend('T_s','Q_0','Q_0(\gamma)+Q_b','Q_0(\gamma)+Q_b+u^.\nablaT+K\nabla^2T','\partial_tT_s',...




  if ( doADCP && (~isfield(stn,'adcp_u') || ~isfield(stn,'adcp_v')) )
    warning('DOADCP flag set with no available ADCP data!');
    doADCP = false;
  end;




,ufld,vfld
  if ( ~strcmp(ufld,[KMPFX '_u']) || ~strcmp(vfld,[KMPFX '_v']) )
    disp(['Using ocean currents ',ufld,',',vfld]);
  end;





            if ( doPlot )

              % annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,[],[1:3]);
              annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,afld,stn.commentstr);





%%%% DEBUG???
    advfacs = { ...
        1,...
              };
    kths = { ...
        0 ...
           };
%%%% DEBUG???



%%%% DEBUG???
    advfacs = { ...
        1,...
              };
%%%% DEBUG???


%%%% DEBUG???
    advfacs = { ...
        0,...
              };
    kths = { ...
        0 ...
           };
%%%% DEBUG???





    % Fluxes WITH or WITHOUT warm-layer adjustment
    % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
    %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
    %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);
    stn.all_zeros = stn.(sfld); stn.all_zeros.data(:) = 0;
    stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                            pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
                            Dfld,'all_zeros','all_zeros',wpfld,whfld,pblzfld,doWarm,max_wl);
                            % Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);
                            % % Dfld,ssufld,ssvfld,wpfld,whfld,pblzfld,doWarm,max_wl);
                            % % % Dfld,ssufld,ssvfld,wpfld,whfld,600,doWarm,max_wl);
                            % % % % Dfld,netufld,netvfld,wpfld,whfld,pblzfld,doWarm,max_wl);







    % Fluxes WITH or WITHOUT warm-layer adjustment
    % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
    %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
    %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);
    stn.all_zeros = stn.(sfld); stn.all_zeros.data(:) = 0;
    stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                            pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
                            Dfld,'all_zeros','all_zeros',wpfld,whfld,pblzfld,doWarm,max_wl);
                            % Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);
                            % % Dfld,ssufld,ssvfld,wpfld,whfld,pblzfld,doWarm,max_wl);
                            % % % Dfld,ssufld,ssvfld,wpfld,whfld,600,doWarm,max_wl);
                            % % % % Dfld,netufld,netvfld,wpfld,whfld,pblzfld,doWarm,max_wl);



%%%% DEBUG???
  % % Basic sanity check
  % gam(0>=gam | gam>1) = 1;
%%%% DEBUG???





          % ylim([21,35]);
          % ylim([15,50]);
          % ylim([15,35]);
          ylim([16,34]);
          ylabel('^oC');
          % appendtitlename([' (' strrep(KMPFX,'_','\_') ' K_\theta=' num2str(kth,'%g,') ')']);




%%%%??? DEBUG
begyr = 1996;
if ( strcmpi(KMPFX,'none') ); begyr = 1970; end;




    % kds = { ...
    %     0.20 ...
    %     0.25 ...
    %     0.30 ...
    %     [0.050,0.350,137] ...
    %       };


%%%%    stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld,opts);





  doWarms = [false, true];
  % doWarms = [true];
  % doWarms = [false];

%%%%%% REALLY GOOD!!! MLRF1
%      [0.05,0.35,45],...
%%%%%% THE BEST SO FAR!!! MLRF1
%      [0.045,0.375,45],...

  kds = {...
      % 0.1,[0.05,0.15,45],[0.05,0.15,137],[0.05,0.15,228],[0.05,0.15,320],...
      % 0.2,[0.10,0.30,45],[0.10,0.30,137],[0.10,0.30,228],[0.10,0.30,320],...
      % 0.3,[0.20,0.40,45],[0.20,0.40,137],[0.20,0.40,228],[0.20,0.40,320],...
      % 0.05, 0.1, 0.2, ...
      % 0.2, 0.3, 0.4,...
      % 0.4, 0.5, 0.6,...
%      [0.045,0.175,45],...
%      [0.045,0.375,45],...
      % [0.20,0.37,228],...
      % [0.25,0.30,320],...
      % [0.20,0.38,228],...
%      [0.20,0.40,228],...
% LONF1
%      [0.20,0.37,228],...
% SANF1
      [0.035,0.175,45],...
      % 0.21,[0.045,0.375,45],[0.045,0.375,137],[0.045,0.375,228],[0.045,0.375,320],...
      % 0.4, [0.300,0.500,45],[0.300,0.500,137],[0.300,0.500,228],[0.300,0.500,320],...
%      0.50,...
%      [0.50,0.70,45],...
      % 0.3,[0.20,0.40,45],[0.20,0.40,137],[0.20,0.40,228],[0.20,0.40,320],...
      % 0.4,[0.30,0.50,45],[0.30,0.50,137],[0.30,0.50,228],[0.30,0.50,320],...
      % 0.6,[0.50,0.70,45],[0.50,0.70,137],[0.50,0.70,228],[0.50,0.70,320],...
      % 0.2,[0.10,0.30,0],[0.10,0.30,91],[0.10,0.30,182],[0.10,0.30,274],...
      % 0.4,[0.30,0.50,0],[0.30,0.50,91],[0.30,0.50,182],[0.30,0.50,274],...
      % 0.6,[0.50,0.70,0],[0.50,0.70,91],[0.50,0.70,182],[0.50,0.70,274],...
      % 0.8,[0.70,0.90,0],[0.70,0.90,91],[0.70,0.90,182],[0.70,0.90,274],...
        };

  % % cbds = {0,3.8e-5,3.8e-4,8.0e-4,16.0e-4,24.0e-4};
  % % cbds = {3.8e-4,8.0e-4};
  cbds = {8.0e-4};
  % cbds = {2.0e-4};

  advfacs = { ...
      0,...
            };

  kths = {...
% SANF1
      0 ,...
      % 2 ,...
      % [0 ,3 ,183],...
      % [0 ,4 ,160],...
      % [0 ,5 ,183],...
      % [0 ,6 ,183],...
      % [0 ,10,183],...
% LONF1
%      [0 ,5 ,183],...
      % [1 ,6 ,183],...
      % [2 ,6 ,183],...
      % 4,...
      % .5,[0 ,1 ,0],[0 ,1 ,91],[0 ,1 ,183],[0 ,1 ,274],...
      % 2 ,[1 ,3 ,0],[1 ,3 ,91],[1 ,3 ,183],[1 ,3 ,274],...
      % 3 ,[2 ,4 ,0],[2 ,4 ,91],[2 ,4 ,183],[2 ,4 ,274],...
      % 4 ,[3 ,5 ,0],[3 ,5 ,91],[3 ,5 ,183],[3 ,5 ,274],...
      % 5 ,[4 ,6 ,0],[4 ,6 ,91],[4 ,6 ,183],[4 ,6 ,274],...
      % [0,2,0],[1,3,0],[3,5,0],[5,7,0],[7,9,0],...
      % 20,[10,30,0],[10,30,91],[10,30,183],[10,30,274],...
      % 0 ,...
      % [4 ,6 ,274],...
      % 10,[ 5,15,0],[ 5,15,91],[ 5,15,183],[ 5,15,274],...
      % 20,[10,30,0],[10,30,91],[10,30,183],[10,30,274],...
      %
      % 0,[ 5,15,45],[ 5,15,137],[ 5,15,228],[ 5,15,320],...
      %   [10,30,45],[10,30,137],[10,30,228],[10,30,320],...
      % 0,...
      % 10,[ 5,15,228],[ 5,15,320],...
      % 20,[10,30,320],...
         };









    end; %for yr = firstyear:lastyear

    disp(['Saving ',matfname]);
    save(matfname,'result');

  end; %if ( exist(matfname,'file') ) else


  stn = merge_station_data(stn,result);
  result = []; clear result
  %DEBUG:  disp('Merged');

  % Calculate hourly averages of all fields
  for fldix=1:numel(flds)
    fld = flds{fldix};
    if ( ~isempty(fld) )
      fld = ['rsmas_',fld];
      rawfld = ['raw_',fld];
      %DEBUG:      disp(fld);
      stn.(fld) = interp_ts(stn.(rawfld));
    end;
  end;







          %% Km-scale Heat Diffusion
          if ( opts.calculate_diffusion )
            stn = station_calc_kdel2t(stn,opts.K_theta,Tfld,...
                                      rawkd2Tfld,kd2Tfld,...
                                      bq0tfld,bdTfld,opts.grid_interp_method,false);
            stn = station_heat_flux_term_inverse(stn,kd2Tffld,kd2Tfld,sfld,[],mhfld);
          else
            stn.(bdTfld) = stn.(bq0tfld);
          end; %if ( isempty(opts.laplacian_climatology) ) else
          stn = station_heat_flux_term_inverse(stn,bdTffld,bdTfld,sfld,[],mhfld);
          if ( isfield(stn,bdTflpfld) ); stn = rmfield(stn,bdTflpfld); end;
          stn = verify_variable(stn,bdTflpfld);






%%%% DEBUG???
% max_wl = nanmax(stn.(mhfld).data) - 1.5;


    % Fluxes WITH or WITHOUT warm-layer adjustment
    % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
    %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
    %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);
    stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                            pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
                            Dfld,ssufld,ssvfld,wpfld,whfld,pblzfld,doWarm,max_wl);
                            % Dfld,ssufld,ssvfld,wpfld,whfld,600,doWarm,max_wl);
                            % % Dfld,netufld,netvfld,wpfld,whfld,pblzfld,doWarm,max_wl);




    doWarms = [true];
    kds = { ...
        [0.015,0.350,137] ...
        [0.015,0.350,183] ...
        [0.015,0.350,228] ...
        0.2 ...
        [0.100,0.350,137] ...
        [0.100,0.350,183] ...
        [0.100,0.350,228] ...
        0.2 ...
          };
    kths = { ...
        0 ...
        [0,20, 45] ...
           };




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG???
            cum = cum(30:end);
            tid = tid(30:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG???






%%%% ??? DEBUG:
stn.erai_sigwavehgt.data(get_year(stn.erai_sigwavehgt.date)<1999)=[];
stn.erai_sigwavehgt.date(get_year(stn.erai_sigwavehgt.date)<1999)=[];
stn.erai_peakwaveper.data(get_year(stn.erai_peakwaveper.date)<1999)=[];
stn.erai_peakwaveper.date(get_year(stn.erai_peakwaveper.date)<1999)=[];
stn.erai_peakwavedir.data(get_year(stn.erai_peakwavedir.date)<1999)=[];
stn.erai_peakwavedir.date(get_year(stn.erai_peakwavedir.date)<1999)=[];



FWYF1:
    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves
    doWarms = [true];
    kds = { ...
        [0.015,0.350,137] ...
        [0.015,0.350,183] ...
        [0.015,0.350,228] ...
        [0.100,0.350,137] ...
        [0.100,0.350,183] ...
        [0.100,0.350,228] ...
          };
    kths = { ...
        0 ...
        [0,20, 45] ...
           };



    % NDBC met, ERAI barom/spechumid, new Ppen, WW3 waves
    % doWarms = [true];
    doWarms = [false];
    kds = { ...
        [0.010,0.350, 45] ...
        [0.010,0.350, 91] ...
        [0.010,0.350,137] ...
        [0.010,0.350,183] ...
        [0.015,0.350, 45] ...
        [0.015,0.350, 91] ...
        [0.015,0.350,137] ...
        [0.015,0.350,183] ...
          };
    kths = { ...
        0 ...
        [5,15,  0] ...
        [0,20,  0] ...
        [5,15, 45] ...
        [0,20, 45] ...
        [5,15, 90] ...
        [0,20, 90] ...
           };



    % NDBC met, ERAI barom/spechumid, new Ppen
    doWarms = [true];
    % doWarms = [false];
    kds = { ...
        [0.050,0.250, 45] ...
        [0.050,0.250, 91] ...
        [0.050,0.250,137] ...
        [0.050,0.350, 45] ...
        [0.050,0.350, 91] ...
        [0.050,0.350,137] ...
          };
    kths = { ...
        0 ...
        10 ...
        [0,20, 45] ...
           };





  legs={};
  legs{end+1}=['Actual'];
  legs{end+1}=['Gramer&Mariano'];
  %legs{end+1}=['GoM HYCOM'];
  legs{end+1}=['ERA-Interim'];
  legs{end+1}=['NCEP NARR'];
  legs{end+1}=['OAFlux 1^o'];
  legs{end+1}=['Large&Yeager CORE.2'];
  legs{end+1}=['NOCS v2'];
  % legs{end+1}=['Actual (',num2str(get_year(s.date(1))),'-',num2str(get_year(s.date(end))),')'];
  % legs{end+1}=['Gramer&Mariano (',num2str(get_year(q.date(1))),'-',num2str(get_year(q.date(end))),')'];
  % %legs{end+1}=['GoM HYCOM (',num2str(get_year(g.date(1))),'-',num2str(get_year(g.date(end))),')'];
  % legs{end+1}=['ERA-Interim (',num2str(get_year(e.date(1))),'-',num2str(get_year(e.date(end))),')'];
  % legs{end+1}=['NCEP NARR (',num2str(get_year(n.date(1))),'-',num2str(get_year(n.date(end))),')'];
  % legs{end+1}=['OAFlux 1^o (',num2str(get_year(o.date(1))),'-',num2str(get_year(o.date(end))),')'];
  % legs{end+1}=['Large&Yeager CORE.2 (',num2str(get_year(l.date(1))),'-',num2str(get_year(l.date(end))),')'];
  % legs{end+1}=['NOCS v2 (',num2str(get_year(x.date(1))),'-',num2str(get_year(x.date(end))),')'];



  legend('Actual','Gramer&Mariano','GoM HYCOM',...
         'ERA-Interim','NCEP NARR','OAFlux','Large&Yeager CORE.2','NOCS v2');



  % fmg;
  % plot(unique(yrmo),mo30a,'m-');
  % plot(unique(hfyrmo),mohf,'kd');
  % plot(unique(ncepyrmo),moncep,'cs');
  % plot(unique(gomyrmo),mogom,'r--');
  % plot(stn.landy_net_heat_flux.date,stn.landy_net_heat_flux.data,'bo:');
  % plot(stn.monthly_nocs_net_heat_flux.date,stn.monthly_nocs_net_heat_flux.data,'g^--');
  % ylabel('Mean Monthly Net Heat Flux [ W / m^2 ]'); 
  % xlim([datenum(2002,12,1),datenum(2007,2,1)]);
  % datetick3;
  % % set(gca,'color',[.95 .95 .95]);
  % title(['Estimates of Q_0 at ' stn.station_name]);
  % legend('Gramer & Mariano','Smith (1988) bulk formulae','NCEP NARR 32km reanalysis','GoM 4km HYCOM + NCODA',...
  %        'Large & Yeager (2009)','NOC Southampton (2009)', 'Location','Best');
  % % print('-dtiff','../figs/mlrf1-cmp_monthly_fluxes.tiff');



%%%% ??? HACK HACK HACK
disp('%%%% ??? HACK HACK HACK - forcing ERAI air_t');
afld='erai_air_t';
%%%% ??? HACK HACK HACK





FWYF1 - WARM
    doWarms = [false];
    kds = { ...
        [0.035,0.220,137], ...
        [0.035,0.220, 91], ...
        [0.025,0.220,137], ...
        [0.025,0.220, 91], ...
        [0.025,0.220, 45], ...
          };
    kths = { ...
        0, ...
        [0,10,320] ...
        [0,10,  0] ...
        [0,10, 45] ...
           };






        [0,20, 45] ...



    kths = { ...
        0, ...
        [0,10, 45] ...
        [0,20, 45] ...
        [5,35,  0] ...
           };


    kds = { ...
        [0.035,0.220, 76], ...
        [0.050,0.250, 76], ...
        [0.035,0.220, 91], ...
        [0.050,0.250, 91], ...
          };





    kds = { ...
        [0.035,0.220,137], ...
        [0.050,0.250,137], ...
        [0.035,0.220, 76], ...
        [0.050,0.250, 76], ...
        [0.050,0.250, 45], ...
          };



    % ERAI met, NDBC air_t, new Ppen - Kd optimized to OAFlux
    kds = { ...
        [0.050,0.250,137], ...
          };
    kths = { ...
        -5, ...
        [ -8, -6, 45] ...
        [-10, -5, 45] ...
        [ -8, -6, 91] ...
        [-10, -5, 91] ...
        [ -8, -6,137] ...
        [-10, -5,137] ...
           };









    sq01dfld = [sq0fld '_1_d_avg'];
    stn = verify_variable(stn,sq01dfld);
    [cix,six] = intersect_dates(stn.(climq0fld).date,stn.(sq01dfld).date);
    x.station_name = stn.station_name;
    x.(climq0fld).date=stn.(climq0fld).date(cix); x.(climq0fld).data=stn.(climq0fld).data(cix);
    x.(sq01dfld).date=stn.(sq01dfld).date(six); x.(sq01dfld).data=stn.(sq01dfld).data(six);
    % station_boxplots(x,climq0fld,'OAFlux Q_0',[-1000,1000],[],[],0,[0,1,0,0],true);
    fmg;
    boxplot_ts(stn.(climq0fld),'month','index',cix,...
               'title',[upper(stn.station_name),' OAFlux Q_0']);
    ylim([-1000,1000]); ylabel('W/m^2');
    % station_boxplots(x,sq01dfld,'Gramer&Mariano Sea-surface Q_0',[-1000,1000],[],[],0,[0,1,0,0],true);
    fmg;
    boxplot_ts(stn.(sq01dfld),'month','index',six,...
               'title',[upper(stn.station_name),' Gramer&Mariano Sea-surface Q_0']);
    ylim([-1000,1000]); ylabel('W/m^2');

    sr1dfld = [srfld '_1_d_avg'];
    stn = verify_variable(stn,sr1dfld);
    [cix,six] = intersect_dates(stn.(climsrfld).date,stn.(sr1dfld).date);
    x.station_name = stn.station_name;
    x.(climsrfld).date=stn.(climsrfld).date(cix); x.(climsrfld).data=stn.(climsrfld).data(cix);
    x.(sr1dfld).date=stn.(sr1dfld).date(six); x.(sr1dfld).data=stn.(sr1dfld).data(six);
    % station_boxplots(x,climsrfld,'OAFlux Q_S_W',[-1000,1000],[],[],0,[0,1,0,0],true);
    fmg;
    boxplot_ts(stn.(climsrfld),'month','index',cix,...
               'title',[upper(stn.station_name),' OAFlux Q_S_W']);
    ylim([-1000,1000]); ylabel('W/m^2');
    % station_boxplots(x,sr1dfld,'Gramer&Mariano Q_S_W',[-1000,1000],[],[],0,[0,1,0,0],true);
    fmg;
    boxplot_ts(stn.(sr1dfld),'month','index',six,...
               'title',[upper(stn.station_name),' Gramer&Mariano Q_S_W']);
    ylim([-1000,1000]); ylabel('W/m^2');

    qlh1dfld = [qlhfld '_1_d_avg'];
    stn = verify_variable(stn,qlh1dfld);
    [cix,six] = intersect_dates(stn.(climqlhfld).date,stn.(qlh1dfld).date);
    x.station_name = stn.station_name;
    x.(climqlhfld).date=stn.(climqlhfld).date(cix); x.(climqlhfld).data=stn.(climqlhfld).data(cix);
    x.(qlh1dfld).date=stn.(qlh1dfld).date(six); x.(qlh1dfld).data=stn.(qlh1dfld).data(six);
    % station_boxplots(x,climqlhfld,'OAFlux Q_L_H',[-1000,1000],[],[],0,[0,1,0,0],true);
    fmg;
    boxplot_ts(stn.(climqlhfld),'month','index',cix,...
               'title',[upper(stn.station_name),' OAFlux Q_L_H']);
    ylim([-1000,1000]); ylabel('W/m^2');
    % station_boxplots(x,qlh1dfld,'Gramer&Mariano Q_L_H',[-1000,1000],[],[],0,[0,1,0,0],true);
    fmg;
    boxplot_ts(stn.(qlh1dfld),'month','index',six,...
               'title',[upper(stn.station_name),' Gramer&Mariano Q_L_H']);
    ylim([-1000,1000]); ylabel('W/m^2');

    asr1dfld = [asrfld '_1_d_avg'];
    stn = verify_variable(stn,asr1dfld);
    qcool1dfld = [qcoolfld '_1_d_avg'];
    stn = verify_variable(stn,qcool1dfld);
    fmg;
    % sh=boxplot(stn.(sr1dfld).data,get_month(stn.(sr1dfld).date), 'notch','on', 'whisker',2, 'symbol','rx','colors','rrr');
    sh=boxplot_ts(stn.(sr1dfld),'month','allcolors','r');
    % lh=boxplot(stn.(qcool1dfld).data,get_month(stn.(qcool1dfld).date), 'notch','on', 'whisker',2, 'symbol','b+','colors','bbb');
    lh=boxplot_ts(stn.(qcool1dfld),'month','allcolors','b');
    ylim([-1000,1000]); ylabel('W/m^2');
    legend([sh(1),lh(1)], 'Q_S_W','Q_L_W+Q_L_H+Q_S_H', 'Location','South');
    titlename([upper(stn.station_name),' Gramer&Mariano Air-Sea fluxes (1d avg)']);

    hcdTdthcf1d = [hcdTdthcf '_1_d_avg'];
    stn = verify_variable(stn,hcdTdthcf1d);
    fmg;
    % boxplot(stn.(hcdTdthcf1d).data,get_month(stn.(hcdTdthcf1d).date), 'notch','on', 'whisker',2, 'symbol','k+','colors','kkk');
    boxplot_ts(stn.(hcdTdthcf1d),'month','allcolors','k',...
               'title',[upper(stn.station_name),' Gramer&Mariano Horizontal Convection (1d avg)']);
    ylim([-1000,1000]); ylabel('W/m^2');









    % ERAI met, NDBC air_t, new Ppen - Kd optimized to OAFlux
    kds = { ...
        [0.050,0.250,137], ...
          };
    kths = { ...
        0, ...
        -5, ...
        -6, ...
        -7, ...
        -8, ...
           };





    % ERAI met, NDBC air_t, new Ppen - Kd optimized to OAFlux
    kds = { ...
        [0.050,0.250,137], ...
          };
    kths = { ...
        0, ...
        [ -5, 5, 45], ...
        -5, ...
           };







    % kds = { ...
    %     [0.050,0.200, 55], ...
    %       };
    % kths = { ...
    %     5, ...
    %        };
    kds = { ...
        [0.050,0.200, 55], ...
        [0.030,0.250,165], ...
        [0.050,0.300,165], ...
        [0.050,0.300,150], ...
        [0.100,0.250,150], ...
          };
    kths = { ...
        5, ...
        [5,10,  0], ...
        [5,10, 45], ...
        [5,10, 91], ...
        [5,10,137], ...
        [0,15,  0], ...
        [0,15, 45], ...
        [0,15, 91], ...
        [0,15,137], ...
           };
    kds = { ...
        [0.100,0.250,150], ...
           };
    kths = { ...
        [0,15, 45], ...
           };

    kds = { ...
        [0.050,0.250,137], ...
        [0.050,0.250,160], ...
        [0.050,0.250,183], ...
          };
    kths = { ...
        0, ...
        7, ...
        [ -5,10, 45], ...
        [-15,15, 45], ...
           };

    kds = { ...
        [0.050,0.250,183], ...
        [0.050,0.250,160], ...
        [0.050,0.250,137], ...
          };
    kths = { ...
        0, ...
        [ -5,10, 45], ...
        [ -5,10,137], ...
        [ -5, 5, 45], ...
           };
    kds = { ...
        [0.050,0.250,183], ...
        [0.050,0.250,160], ...
        [0.050,0.250,137], ...
          };
    kths = { ...
        0, ...
        [ -5,10, 45], ...
        [ -5,10,137], ...
        [ -5, 5, 45], ...
           };












    % ERAI met, NDBC air_t, new Ppen - Kd optimized to OAFlux
    kds = { ...
        [0.050,0.200, 55], ...
        [0.100,0.200, 55], ...
        [0.100,0.200,137], ...
        [0.100,0.200,228], ...
        [0.100,0.200,320], ...
        [0.100,0.300, 55], ...
        [0.100,0.300,137], ...
        [0.100,0.300,228], ...
        [0.100,0.300,320], ...
          };
    kths = { ...
        0, ...
        2, ...
        5, ...
           };




From XSPEC.m:
    % ix = find( get_daylight(stn.kd.date,stn.lat,stn.lon,35) & ...
    %            ismember(get_hour(stn.kd.date),15:17) & ...
    %            (stn.kd.data>0.00001) );





      % keepix = union(keepix,find( (18/24) < dtdif & dtdif < 10 ));


    % [cum,tid] = grp_ts(stn.kd.data,stn.kd.date,@get_jday,@nanmedian,2);
    % fmg; plot(tid,cum,'*'); xlim([0,366]); datetick3('x',3,'keeplimits');






               ismember(get_hour(stn.kd.date),14:18) & ...

               ismember(get_hour(kd.date),14:20) & ...

               ismember(get_hour(kd.date),14:18) & ...



From OPTIM_Q0.m:
   case 'mlrf1',
    % % ERAI met, air_t, old Ppen
    % kds = { [0.045,0.375,45], };
    % kths = { 0 };

    % % ERAI met, NDBC air_t, old Ppen
    % kds = { [0.035,0.190, 55], };
    % kths = { 0 };

    % ERAI met, NDBC air_t, new Ppen - Kd optimized to OAFlux
    kds = { ...
        [0.050,0.200, 55], ...
        [0.100,0.200, 55], ...
        [0.100,0.200,137], ...
        [0.100,0.200,228], ...
        [0.100,0.200,320], ...
        [0.100,0.300, 55], ...
        [0.100,0.300,137], ...
        [0.100,0.300,228], ...
        [0.100,0.300,320], ...
          };
    kths = { ...
        0, ...
        2, ...
        5, ...
           };








    data = stn.ndbc_sea_t_erai_erai_30a_wind_stress.data;
    save(matfname,'date','data');
    clear date data;
  end;





  % doPlot = [0,1,1];
  % % doPrint = [0,1,1];
  % doPrint = [0,0,0];




    % ERAI met, NDBC air_t, new and improved higher Ppen
    kds = { ...
        [0.060,0.380, 30], ...
        [0.060,0.380, 45], ...
        [0.060,0.380, 60], ...
        [0.060,0.400, 30], ...
        [0.060,0.400, 45], ...
        [0.060,0.400, 60], ...
        [0.060,0.420, 30], ...
        [0.060,0.420, 45], ...
        [0.060,0.420, 60], ...
          };



        [0.060,0.420, 45], ...
        [0.060,0.500, 45], ...



    % ERAI met, NDBC air_t, new and improved higher Ppen
    kds = { ...
        [0.060,0.300,  0], ...
        [0.060,0.300, 45], ...
        ...
        [0.060,0.400,  0], ...
        [0.060,0.400, 45], ...
        ...
        [0.060,0.500,  0], ...
        [0.060,0.500, 45], ...
        ...
        [0.060,0.600,  0], ...
        [0.060,0.600, 45], ...
          };




    % ERAI met, NDBC air_t, new and improved higher Ppen
    kds = { ...
        [0.045,0.375, 45], ...
        [0.050,0.375, 45], ...
        [0.055,0.375, 45], ...
        [0.055,0.380, 45], ...
        [0.035,0.800, 45], ...
        [0.035,0.900, 45], ...
          };



    % % ERAI met, NDBC air_t, higher Ppen - amplitude problem but good annual...
    % kds = { [0.045,0.375, 91], }

    % ERAI met, NDBC air_t, new and improved higher Ppen
    kds = { ...
        [0.045,0.375, 45], ...
        [0.045,0.375, 91], ...
        ...
        [0.050,0.375, 45], ...
        [0.040,0.375, 91], ...
          };




    kds = { ...
        [0.045,0.375, 91], ...
        [0.045,0.375,228], ...
        ...
        [0.045,0.400, 45], ...
        [0.045,0.400,228], ...
        ...
        [0.100,0.200,228], ...
          };



    kds = { ...
        [0.050,0.350,  0], ...
        [0.050,0.350, 45], ...
        [0.050,0.350, 91], ...
        [0.050,0.350,137], ...
        ...
        [0.045,0.375,  0], ...
        [0.045,0.375, 45], ...
        [0.045,0.375, 91], ...
        [0.045,0.375,137], ...
        ...
        [0.045,0.400,  0], ...
        [0.045,0.400, 45], ...
        [0.045,0.400, 91], ...
        [0.045,0.400,137], ...
          };




   case 'smkf1',
    % ERAI met, NDBC air_t, old Ppen
    kds = { [0.045,0.375,45], };
    kths = { 0 };

    % ERAI met, NDBC air_t, new and improved higher Ppen
    kds = { ...
        [0.045,0.700, 45], ...
        [0.045,0.700, 91], ...
        [0.045,0.700,137], ...
        [0.035,0.800, 45], ...
        [0.035,0.800, 91], ...
        [0.035,0.800,137], ...
          };
    kths = { 0 };




        0.2, ...
        [0.200,0.370,228], ...
        [0.050,0.600,  0], ...
        [0.050,0.600, 45], ...





          % pos(1) = 0;



  % Use only months with at least 25 days worth of data
  for varix=1:nvars
    for stix=1:nstns
      for yr=yrs(:)';
        for mo=1:12
          ix = find(get_year(stns{stix}.(vars{varix}).date)==yr&get_month(stns{stix}.(vars{varix}).date)==mo);
          if ( numel(ix) < (25*24) )
            stns{stix}.(vars{varix}).date(ix) = [];
            stns{stix}.(vars{varix}).data(ix) = [];
          end;
        end;
      end;
    end;
  end;







    res = intersect_all_dates([],dts{:});
    ixen(1:nstns,varix) = res{:};



    for stix=1:nstns
      ixen(stix,varix) = res{stix};
    end;




  %stnms = { 'fwyf1','mlrf1','lonf1','smkf1','sanf1' };




    stn.h.date = [];
    stn.h.data = [];

        ix = find(get_year(stn.ndbc_tide.date) == yr);
        stn.h.date(end+1:end+numel(ix),1) = stn.ndbc_tide.date(ix);
        stn.h.data(end+1:end+numel(ix),1) = stn.ndbc_tide.data(ix);



    stn.h.data(get_year(stn.h.date)>=2010) = [];
    stn.h.date(get_year(stn.h.date)>=2010) = [];
    [B,Stats,fh,lh] = scatter_fit(stn.h.date,stn.h.data);
    ylim([-1.5,+3.5]);
    {numel(Stats.w),Stats.se(2),Stats.t(2),Stats.p(2),},
    set(lh(1),'Marker','.','MarkerSize',1,'LineStyle','none','LineWidth',1,'Color',[.3,.3,.3]);
    set(lh(2),'Marker','none','LineStyle','-','LineWidth',3,'Color','k');









    stn.t.date = [datenum(1987,1,1):(1/24):datenum(2012,1,1)]';
    stn.t.data = repmat(0,size(stn.t.date));

    [tix,nix] = intersect_dates(stn.t.date,stn.ndbc_sea_t.date);
    stn.t.date(tix) = stn.ndbc_sea_t.date(nix);
    stn.t.data(tix) = stn.ndbc_sea_t.data(nix);





    stn.t.date = [];
    stn.t.data = [];

    yrs = unique(get_year(stn.ndbc_sea_t.date));
    for yr=1987:2011
      ix = find(get_year(stn.ndbc_sea_t.date)==yr);
      if ( numel(ix) >= 8000 )
        stn.t.date(end+1:end+numel(ix),1) = stn.ndbc_sea_t.date(ix);
        stn.t.data(end+1:end+numel(ix),1) = stn.ndbc_sea_t.data(ix);
      else
        stn.t.date(end+1:end+8760,1) = datenum(yr,1,1):(1/24):(datenum(yr,1,1)+364.99);
        stn.t.data(end+1:end+8760,1) = repmat(0,[8670,1]);
      end;
    end;








    x.station_name = stn.station_name;
    x.t = stn.ndbc_sea_t;
    [n,yrs] = grp_ts(x.t.data,x.t.date,@get_year,@numel,0);

    fullyrs = yrs(n>=8000);
    gapyrs = yrs(n<8000);

    x.t.data = x.t.data(ismember(get_year(x.t.date),fullyrs));
    x.t.date = x.t.date(ismember(get_year(x.t.date),fullyrs));
    gapstuffing.t.date = datenum(gapyrs,7,1);
    gapstuffing.t.data = repmat(0,[numel(gapyrs) 1]);
    x = merge_station_data(x,gapstuffing);
    clear gapstuffing;

    station_anova_multcompare(x,'t',[],'^oC',[],[1,0,0,0]);
    view(270,90);
    xlim([25.5,28]);






    station_anova_multcompare(stn,'ndbc_sea_t',@(x)(find(ismember(get_year(x.date),[1993:1996,1998:2001,2003,2005:2007,2009]))),'^oC',[],[1,0,0,0]);

    station_anova_multcompare(x,'t',[],'^oC',[],[1,0,0,0]);

    @(x)(find(ismember(get_year(x.date),[1993:1996,1998:2001,2003,2005:2007,2009])))

    gapstuffing.date = datenum(gapyrs,1,1);



DRYF1
        0.2, ...
        0.3, ...
        0.4, ...
        [0.20,0.40,183], ...
        [0.20,0.40,228], ...
        [0.20,0.40,274], ...
        ...
        [0.25,0.35,320], ...
        [0.20,0.35,320], ...
        [0.25,0.30,320], ...
        ...
        [0.30,0.50,183], ...
        [0.30,0.50,228], ...
        [0.30,0.50,274], ...
        [0.30,0.50,320], ...




if ( strcmpi(KMPFX,'none') ); begyr = 1970; end;






  vars = {'ndbc_air_t', 'ndbc_sea_t', 'tau_xshore', 'ndbc_wind1_speed'};
  ylbs = {'T_a',        'T_s',        '\tau^x^s',   'U'};
  ylms = {[2,36],       [2,36],       [-1,+1],      [2,36]};



        stns{stix} = verify_variable(stns{stix},'ndbc_wind1_u');
        stns{stix}.ndbc_wind1_u.data = kts2mps(stns{stix}.ndbc_wind1_u.data);
        stns{stix} = verify_variable(stns{stix},'ndbc_wind1_v');
        stns{stix}.ndbc_wind1_v.data = kts2mps(stns{stix}.ndbc_wind1_v.data);




  fmg;
  for varix=1:nvars
    for stix=1:nstns
      %plotix = (stix-1)*nvars+varix;
      plotix = (varix-1)*nstns+stix;
      subplot_tight(nvars,nstns,plotix);
      mos = get_month(stns{stix}.(vars{varix}).date(ixen{stix,varix}));
      dat = stns{stix}.(vars{varix}).data(ixen{stix,varix});
      boxplot(dat, mos, 'notch','on', 'whisker',2);
      xlim([0.5 12.5]);
      ylim(ylms{varix});
      if ( plotix<=nstns ); titlename(upper(stnms{stix})); end;
      %if ( mod(plotix,nstns) == 1 ); ylabel(ylbs{varix}); end;
      ylabel(ylbs{varix}); 
      %xlabel('Year-Month');
      %ylabel([upper(stnms{stix}) ' ' ylbs{varix}]);
    end;
  end;










   case 'mlrf1',
    % % ERAI met, air_t
    % kds = { [0.045,0.375,45], };
    % kths = { 0 };

    % ERAI met, NDBC air_t
    kds = { [0.035,0.200, 45], };
    kds = { [0.035,0.195, 45], };
    kds = { [0.035,0.190, 55], };
    kths = { 0 };







OPTIM_Q0:

if (1)
      fmg; plot_ts(stn.optim_q0.climt,squeeze(stn.optim_q0.climsq(:,2,kthix)));
      datetick3('x',3);
      titlename([stnm ' Daily Clim: Warm Layer ' stn.optim_q0.cbdstrs{default_cbdix}]);
      legend({'T',stn.optim_q0.kdstrs{:}});
      xlim(stn.optim_q0.climt.date([1 end]));
      ylim([minmin([stn.optim_q0.climsq.data]),maxmax([stn.optim_q0.climsq.data])]);
      % ylim([21,35]);
      % ylim([15,50]);
      % ylim([15,35]);
      ylim([21,33]);
      appendtitlename([' (' strrep(KMPFX,'_','\_') ' K_\theta=' num2str(kth,'%g,') ')']);
      appendtitlename([' (' num2str(begyr) '-' num2str(endyr) ')']);
end;









      fmg; plot(squeeze(stn.optim_q0.error(:,doWarmix,:))); titlename([stnm ' Total Error: ' doWarmStr]);
      legend(stn.optim_q0.kthstrs);
      set(gca,'XTick',[1:numel(kds)],'XTickLabel',stn.optim_q0.kdstrs);
      ylim([nanmin(stn.optim_q0.error(:)),nanmax(stn.optim_q0.error(:))]);

    fmg; plot(squeeze(stn.optim_q0.error(:,2,:))); titlename([stnm ' Total Error: Warm Layer']);
    legend(stn.optim_q0.kthstrs);
    set(gca,'XTick',[1:numel(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.error(:)),nanmax(stn.optim_q0.error(:))]);

    fmg; plot(squeeze(nanmedian(stn.optim_q0.seaserror(:,1,:,:),4))); titlename([stnm ' Mdn Seas Err: No Warm Layer']);
    legend(stn.optim_q0.kthstrs);
    set(gca,'XTick',[1:numel(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.seaserror(:)),nanmax(stn.optim_q0.seaserror(:))]);

    fmg; plot(squeeze(nanmedian(stn.optim_q0.seaserror(:,2,:,:),4))); titlename([stnm ' Mdn Seas Err: Warm Layer']);
    legend(stn.optim_q0.kthstrs);
    set(gca,'XTick',[1:numel(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.seaserror(:)),nanmax(stn.optim_q0.seaserror(:))]);








      fmg; plot_ts(stn.(sfld),squeeze(stn.optim_q0.sq(:,1,kthix))); titlename([stnm ' Time Series: No Warm Layer']);
      legend({'T',stn.optim_q0.kdstrs{:}});
      ylim([minmin([stn.optim_q0.sq.data]),maxmax([stn.optim_q0.sq.data])]);
      fmg; plot_ts(stn.(sfld),squeeze(stn.optim_q0.sq(:,2,kthix))); titlename([stnm ' Time Series: Warm Layer']);
      legend({'T',stn.optim_q0.kdstrs{:}});
      ylim([minmin([stn.optim_q0.sq.data]),maxmax([stn.optim_q0.sq.data])]);














%%%%DEBUG???
    opts.convective_drag_coefficient = get_opt(opts,'convective_drag_coefficient',2.0e-4);







   case {'looe1'},
    opts.kd = get_opt(opts,'kd',[0.045,0.375,45]);

    opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    opts.convective_drag_coefficient = get_opt(opts,'convective_drag_coefficient',8.0e-4);







  if ( ~exist('stn_or_stnm','var') || isempty(stn_or_stnm) )
    stn.station_name = 'aoat_broad_key_2';
    disp(stn.station_name);
    % Midway between BRDCRK2 transect start and end points
    stn.lon = -80.201641666666660;
    stn.lat = 25.309874999999998;
    stn.depth = 3.9;
  else
    stn = get_station_from_station_name(stn_or_stnm);
  end;










    prevWarns(end+1) = warning('off','get_wera_station:BadFileShape');
    prevWarns(end+1) = warning('off','get_wera_station:SizeMismatch');
    prevWarns(end+1) = warning('off','anwera:NoURL');
    prevWarns(end+1) = warning('off','anwera:BadCoords');
    prevWarns(end+1) = warning('off','anwera:NoLocalAver');
    prevWarns(end+1) = warning('off','anwera:NoLocalAcc');
    prevWarns(end+1) = warning('off','anwera:NoData');
    prevWarns(end+1) = warning('off','anwera:BadData');





    warning('off','anwera:NoURL');
    warning('off','anwera:BadCoords');
    warning('off','anwera:NoLocalAver');
    warning('off','anwera:NoLocalAcc');
    warning('off','anwera:NoData');
    warning('off','anwera:BadData');

    warning('on','anwera:NoURL');
    warning('on','anwera:BadCoords');
    warning('on','anwera:NoLocalAver');
    warning('on','anwera:NoLocalAcc');
    warning('on','anwera:NoData');
    warning('on','anwera:BadData');






%%%% DEBUG???
    for yr=2005:2005
%%%% DEBUG???
jds = 1:2;




            result.wera_u.date(end+1,1) = dt;
            result.wera_u.data(end+1,1) = U(result.wera_ix);
            result.wera_v.date(end+1,1) = dt;
            result.wera_v.data(end+1,1) = V(result.wera_ix);






   case 2006,
    %INTERRUPTED DOWNLOAD
    jds = 54:365;





for yr=2006:2011
for yr=2005:2005
for yr=2006:2011
  switch (yr),
   case {2004,2008,2012,2016},
    jds = 1:366;
   case 2005,
    % ONLY DUE TO INTERRUPTED DOWNLOAD!
    jds = 308:365;
   case 2011,
    jds = 1:33;
   otherwise,
    jds = 1:365;
  end;
  disp({yr,jds});
  for jd=jds(:)'
    for hr=0:23
      ds=sprintf('%04d%03d%02d00',yr,jd,hr);
      anwera(ds,0);
    end;
  end;
end;




%for yr=2005:2011
for yr=2005:2005
%for yr=2006:2011





if( ~isfield(stn,'adcp_baroclinic_btm_u') )
  stn.adcp_baroclinic_btm_u.date = stn.adcp_baroclinic_u.date;
  stn.adcp_baroclinic_btm_u.data = nanmean(stn.adcp_baroclinic_u.prof(:,1:6),2);
  stn.adcp_baroclinic_btm_v.date = stn.adcp_baroclinic_v.date;
  stn.adcp_baroclinic_btm_v.data = nanmean(stn.adcp_baroclinic_v.prof(:,1:6),2);
end;



     %case 'looe1',	stn.(orifld)=77.00; % Straight isobaths in GoogleEarth: Changed GET_LOOE1_ADCP






    stn.adcp_sfc_u.data = nanmean(stn.adcp_u.prof(:,24:29),2);
    stn.adcp_sfc_v.data = nanmean(stn.adcp_v.prof(:,24:29),2);
    stn.adcp_sfc_w.data = nanmean(stn.adcp_w.prof(:,24:29),2);






  % QEPFX - Quasi-Eulerian currents and heat advection
  if ( ~exist('QEPFX','var') || isempty(QEPFX) )
    switch (KMPFX),
     case 'fkeys_hycom',	QEPFX = [WAVEPFX '_fkeys_qe'];
     case 'gom_hycom',		QEPFX = [WAVEPFX '_gom_qe'];
     case 'avhrr_weekly_sst',	QEPFX = [WAVEPFX '_avhrr_qe'];
     case 'none',		QEPFX = [WAVEPFX '_none_qe'];
     otherwise,			error('Unknown km-scale model "%s"',KMPFX);
    end;
  end;






    opts.hc_debug = get_opt(opts,'hc_debug',true);




  kths = {  0,...
            2,[ 1, 3,0],[ 1, 3,91],[ 1, 3,183],[ 1, 3,274],...
           20,[10,30,0],[10,30,91],[10,30,183],[10,30,274],...
         };




      appendtitlename([' (K_\theta=' num2str(kth,'%g,') ')']);



      titlename([stnm ' Daily Clim: No Warm Layer K_\theta=' kthstrs{kthix}]);



function stn = optim_q0(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld)
%function stn = optim_q0(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld)

  doPlot = true;

  set_more off
  timenow,
  tic,

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

  datapath = get_thesis_path('../data');


  stn = get_station_from_station_name(stn_or_stnm);

  opts = station_heat_budget_options(stn);

  if ( ~isfield(stn,bathyfld) )
    stn = get_ngdc_bathy_station(stn);
  end;
  if ( ~isfield(stn,slopefld) )
    stn = station_ngdc_offshore_slope(stn);
  end;
  if ( ~isfield(stn,bathorifld) )
    stn = station_optimal_isobath_orientation(stn);
  end;

  warning('off','Ecoforecasts:mergedNonTS');
  stn = load_all_ndbc_data(stn);

  %%%
  %% Sea temperature (in situ NDBC, site-specific, remote sensing, or reanalysis)

  if ( ~isfield(stn,sfld) )
    if ( regexp(sfld,'misst_sst') )
      stn = get_misst_station(stn);
      stn.hourly_misst_sst = interp_ts(stn.misst_sst);
    elseif ( regexp(sfld,'avhrr_weekly_sst') )
      stn = get_avhrr_weekly_field(stn,true);
      stn.hourly_avhrr_weekly_sst = interp_ts(stn.avhrr_weekly_sst);
    elseif ( regexp(sfld,'erai_') )
      if ( ~strcmpi(ISPFX,'erai') && ~strcmpi(RAPFX,'erai') )
        stn = get_erai_station(stn);
      end;
    elseif ( regexp(sfld,'microcat') )
      stn = get_looe1_microcat(stn);
    elseif ( regexp(sfld,'adcp') )
      stn = get_looe1_adcp(stn);
    elseif ( regexp(sfld,'ndbc') )
      if ( ~strcmpi(ISPFX,'ndbc') )
        stn = load_all_ndbc_data(stn);
      end;
    elseif ( regexp(sfld,'^(sea|ct_|ctd_)') )
      if ( ~strcmpi(ISPFX,'icon') )
        stn = load_station_data(stn);
      end;
    else
      error(['Unknown sea temperature field ',sfld]);
    end;
  end;


  if ( ~isfield(stn,afld) || ~isfield(stn,sfld) )
    switch (ISPFX),
     case 'ndbc',	%stn = load_all_ndbc_data(stn);
     case 'erai',	stn = get_erai_station(stn);
     case 'icon',	stn = load_station_data(stn);
     otherwise,		error('Unknown in situ dataset "%s"',ISPFX);
    end;
  end;
  warning('on','Ecoforecasts:mergedNonTS');

  stn = verify_variable(stn, afld);
  stn = verify_variable(stn, pfld);
  stn = verify_variable(stn, Wfld);
  stn = verify_variable(stn, sfld);


  stn = station_tmd_tide(stn);
  stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
  % Assume max warm-layer depth somewhere near bottom boundary layer top
  max_wl = nanmax(stn.(mhfld).data) - 0.5;


  if ( ~isfield(stn,dsrfld) || ~isfield(stn,cfld) )
    switch (RAPFX),
     case 'erai',		stn = get_erai_station(stn);
     case 'ncep',		stn = get_ncep_station(stn,'narr');
     %%% Not Yet Implemented
     % case 'era40',		stn = get_era40_station(stn);
     % case 'cfsr',		stn = get_ncep_station(stn,'cfsr');
     otherwise,		error('Unavailable gridded/reanalysis dataset "%s"',RAPFX);
    end;
    stn = station_heat_flux_term(stn,raq0fld,raqtfld,sfld,[],nanmean(stn.(mhfld).data));
  end;

  if ( ~isfield(stn,whfld) )
    % If waves not from reanalysis, user must want model (WaveWatch III) or wind estimate
    switch (WAVEPFX),
     case 'ww3',	stn = get_ww3_station(stn);
     case 'ndbc',	stn = station_wind_to_wave(stn,Wfld,Dfld,wpfld,whfld,wdfld);
     case 'erai',	stn = get_erai_station(stn); %IF NOT ALSO OUR REANALYSIS DATASET
     otherwise,		error('Unknown wave source "%s"',WAVEPFX);
    end;
  end;

  %% Low-pass filter winds for quasi-Eulerian currents
  stn = verify_variable(stn,Ulpfld);
  stn = verify_variable(stn,Vlpfld);
  stn.(Wlpfld).date = stn.(Ulpfld).date;
  stn.(Wlpfld).data = uv_to_spd(stn.(Ulpfld).data,stn.(Vlpfld).data);
  stn.(Dlpfld).date = stn.(Ulpfld).date;
  stn.(Dlpfld).data = uv_to_dir(stn.(Ulpfld).data,stn.(Vlpfld).data);

  stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wlpfld,Dlpfld,whfld,wpfld,wdfld);
  stn.(netufld) = ts_op(stn.(tufld),stn.(ssufld),'+');
  stn.(netvfld) = ts_op(stn.(tvfld),stn.(ssvfld),'+');


  %% Kilometer-scale Ocean Data

  if ( ~isfield(stn,ufld) )
    switch (KMPFX),
     case 'fkeys_hycom',
      stn = get_fkeys_hycom(stn,[],[],[],[],'linear');
     case 'gom_hycom',
      stn = get_gom_hycom(stn);
     case 'avhrr_weekly_sst',
      if ( ~isfield(stn,Tfld) )
        disp('Loading AVHRR_WEEKLY_SST instead of hydrodynamic model data...');
        stn = get_avhrr_weekly_field(stn,true);
      end;
      opts.km_scale_advection = false;
      disp('ONLY STOKES');
     case 'none',
      disp('Loading NO hydrodynamic model data...');
      opts.km_scale_advection = false;
      disp('ONLY STOKES');
     otherwise,
      error('Unknown km-scale data source "%s"',KMPFX);
    end;
    more off;
  end;


  %% Surface fluxes

  stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);


  kds = {...
      0.2,[0.10,0.30,0],[0.10,0.30,91],[0.10,0.30,182],[0.10,0.30,274],...
      % 0.4,[0.30,0.50,0],[0.30,0.50,91],[0.30,0.50,182],[0.30,0.50,274],...
      % 0.6,[0.50,0.70,0],[0.50,0.70,91],[0.50,0.70,182],[0.50,0.70,274],...
      % 0.8,[0.70,0.90,0],[0.70,0.90,91],[0.70,0.90,182],[0.70,0.90,274],...
        };

  % cbds = {0,3.8e-5,3.8e-4,8.0e-4,16.0e-4,24.0e-4};
  % cbds = {3.8e-4,8.0e-4};
  cbds = {8.0e-4};

  kths = { 0,...
            2,[ 1, 3,0],[ 1, 3,91],[ 1, 3,183],[ 1, 3,274],...
           20,[10,30,0],[10,30,91],[10,30,183],[10,30,274],...
         };

  % default_cbd = 3.8e-4; %From literature
  default_cbd = 8.0e-4; %Most Keys sites
  [ig,default_cbdix] = min(abs([cbds{:}] - default_cbd));

  for ix=1:length(kds)
    stn.optim_q0.kdstrs{ix} = num2str(kds{ix},'%g,');
  end;
  for ix=1:length(cbds)
    stn.optim_q0.cbdstrs{ix} = num2str(cbds{ix},'C_b_d=%g');
  end;
  for ix=1:length(cbds)
    stn.optim_q0.kthstrs{ix} = num2str(kths{ix},'%g,');
  end;

  for doWarmix=1:2
    doWarm = logical(doWarmix-1);

    % tic,

    % Fluxes WITH or WITHOUT warm-layer adjustment
    % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
    %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
    %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);
    stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                            pfld,sfld,[],[],TURPFX,dsrfld,dlrfld,rfld,...
                            Dfld,netufld,netvfld,wpfld,whfld,pblzfld,doWarm,max_wl);
    % Algorithm sometimes returns complex numbers!
    stn.(qlhfld).data = real(stn.(qlhfld).data);
    stn.(qshfld).data = real(stn.(qshfld).data);
    stn.(qrhfld).data = real(stn.(qrhfld).data);
    stn.(qturfld) = ts_op(stn.(qlhfld),stn.(qshfld),'+');
    if ( isfield(stn,qrhfld) && is_valid_ts(stn.(qrhfld)) )
      stn.(qturfld) = ts_op(stn.(qturfld),stn.(qrhfld),'+');
    end;
    badix = find(~isfinite(stn.(qturfld).data));
    stn.(qturfld).date(badix) = [];
    stn.(qturfld).data(badix) = [];
    stn = station_heat_flux_term(stn,qturfld,qturtfld,sfld,[],mhfld);

    % toc,

    % Just in case something above reset it!
    more off;

    for kdix=1:length(kds)
      kd = kds{kdix};

      opts.kd = kd;
      stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,opts);

      stn.(qradfld) = ts_op(stn.(asrfld),stn.(lrfld),'+');
      badix = find(~isfinite(stn.(qradfld).data));
      stn.(qradfld).date(badix) = [];
      stn.(qradfld).data(badix) = [];
      stn = station_heat_flux_term(stn,qradfld,qradtfld,sfld,[],mhfld);

      stn.(q0fld) = ts_op(stn.(qradfld),stn.(qturfld),'+');
      % badix = find(~isfinite(stn.(q0fld).data));
      % stn.(q0fld).date(badix) = [];
      % stn.(q0fld).data(badix) = [];
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],mhfld);

      % Net flux without absorption calculation or benthic flux - for comparison
      stn.(sqradfld) = ts_op(stn.(srfld),stn.(lrfld),'+');
      stn.(sq0fld) = ts_op(stn.(sqradfld),stn.(qturfld),'+');
      stn = station_heat_flux_term(stn,sq0fld,sqtfld,sfld,[],mhfld);

      for cbdix=1:length(cbds)
        cbd = cbds{cbdix};
        opts.benthic_debug = false;
        opts.convective_drag_coefficient = cbd;

        %% Benthic Heat Exchanges
        stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld,opts);
        % badix = find(~isfinite(stn.(qbofld).data) | abs(stn.(qbofld).data)>2e3);
        badix = find(~isfinite(stn.(qbofld).data));
        stn.(qbofld).date(badix) = [];
        stn.(qbofld).data(badix) = [];
        stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],mhfld);

        stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'-');
        stn = station_heat_flux_term(stn,bq0fld,bq0tfld,sfld,[],mhfld);
        if ( isfield(stn,bq0lpfld) ); stn = rmfield(stn,bq0lpfld); end;
        stn = verify_variable(stn,bq0lpfld);

        %% Km-scale Heat Diffusion
        opts.grid_interp_method = get_opt(opts,'grid_interp_method','linear');
        % kth = 20;
        kth = 2;
        stn = station_calc_kdel2t(stn,kth,Tfld,...
                                  rawkd2Tfld,kd2Tfld,...
                                  bq0tfld,bdTfld,opts.grid_interp_method);
        stn = station_heat_flux_term_inverse(stn,bdTffld,bdTfld,sfld,[],mhfld);
        if ( isfield(stn,bdTflpfld) ); stn = rmfield(stn,bdTflpfld); end;
        stn = verify_variable(stn,bdTflpfld);


        %% Apply Horizontal Convection to total fluxes
        bet = stn.(slopefld);
        [tix,hix,tspdix,aix,Wix,sq0ix,bq0ix,dTix] = ...
            intersect_all_dates([],stn.(sfld).date,stn.(mhfld).date,stn.(tspdfld).date,stn.(afld).date,stn.(Wfld).date,stn.(sq0fld).date,stn.(bq0fld).date,stn.(bdTflpfld).date);
        dts = stn.(sfld).date(tix);
        t = stn.(sfld).data(tix);
        s = repmat(36,size(stn.(sfld).data(tix)));
        h = stn.(mhfld).data(hix);
        tspd = stn.(tspdfld).data(tspdix);
        at = stn.(afld).data(aix);
        W = stn.(Wfld).data(Wix);
        sq0 = stn.(sq0fld).data(sq0ix);
        bq0 = stn.(bq0fld).data(bq0ix);
        dT = stn.(bdTflpfld).data(dTix);

        opts.hc_debug = get_opt(opts,'hc_debug',false);
        % opts.hc_R = get_opt(opts,'hc_R',(1.00-0.08));
        opts.hc_scaling = get_opt(opts,'hc_scaling','US');
        opts.hc_max_onset_secs = get_opt(opts,'hc_max_onset_secs',12*3600);
        res = horizontal_convection(t,s,h,dT,bet,opts,dts,dT,W);
        stn.(hcdTdt).date = dts;
        stn.(hcdTdt).data = res.dTdt;
        stn = station_heat_flux_term_inverse(stn,hcdTdtf,hcdTdt,sfld,[],mhfld);

%%%%??? DEBUG
% stn.(hcdTdt) = stn.(bq0tfld);
%%%%??? DEBUG

        % annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,[],[1:3]);

        if ( doPlot )

         %% Calculate and plot time series and errors from this option
         clear t q
         [tix,qix] = intersect_dates(stn.(sfld).date,stn.(hcdTdt).date);
         t.date = stn.(sfld).date(tix);	t.data = stn.(sfld).data(tix);
         q.date = stn.(hcdTdt).date(qix);	q.data = stn.(hcdTdt).data(qix);

         % Evaluate total error for this option
         t0 = t.data(1);
         sq.date = q.date;
         sq.data = t0 + cumsum(q.data) - q.data(1);
         stn.optim_q0.error(kdix,doWarmix,cbdix) = sqrt(sum((t.data-sq.data).^2));

         stn.optim_q0.doWarm(kdix,doWarmix,cbdix) = doWarm;
         stn.optim_q0.kd{kdix,doWarmix,cbdix} = kd;
         stn.optim_q0.cbd(kdix,doWarmix,cbdix) = cbd;
         stn.optim_q0.q(kdix,doWarmix,cbdix) = q;
         stn.optim_q0.sq(kdix,doWarmix,cbdix) = sq;

         % Evaluate climatological error for this option
         if ( ~isfield(stn.optim_q0,'climt') )
           [cum,tid] = grp_ts(t.data,t.date,'daily');
           stn.optim_q0.climt.data = cum;
           stn.optim_q0.climt.date = tid;
         end;
         % [cum,tid] = grp_ts(q.data,q.date,'daily',@nansum);
         [cum,tid] = grp_ts(q.data,q.date,'daily',@nanmean);
         cum = 24*cum;
         stn.optim_q0.climq(kdix,doWarmix,cbdix).data = cum;
         stn.optim_q0.climq(kdix,doWarmix,cbdix).date = tid;
         t0 = stn.optim_q0.climt.data(1);
         sq.date = tid;
         sq.data = t0 + cumsum(cum) - cum(1);
         stn.optim_q0.climsq(kdix,doWarmix,cbdix).data = sq.data;
         stn.optim_q0.climsq(kdix,doWarmix,cbdix).date = sq.date;
         stn.optim_q0.climerror(kdix,doWarmix,cbdix) = sqrt(sum((stn.optim_q0.climt.data-sq.data).^2));

         % Evaluate seasonal amplitude error for this option
         stn.optim_q0.minseassq = +Inf;
         stn.optim_q0.maxseassq = -Inf;
         yrs = unique(get_year(t.date));
         for yrix = 1:length(yrs)
           yr = yrs(yrix);
           ix = find(get_year(t.date)==yr);
           if ( isempty(ix) )
             warning('NO DATA for year %g??',yr);
           else
             t0 = t.data(ix(1));
             sq.date = q.date(ix);
             sq.data = t0 + cumsum(q.data(ix)) - q.data(ix(1));
             stn.optim_q0.seasyear(kdix,doWarmix,cbdix,yrix) = yr;
             stn.optim_q0.seassq(kdix,doWarmix,cbdix,yrix) = sq;
             stn.optim_q0.minseassq = nanmin(stn.optim_q0.minseassq,nanmin(sq.data(:)));
             stn.optim_q0.maxseassq = nanmax(stn.optim_q0.maxseassq,nanmax(sq.data(:)));
             stn.optim_q0.seaserror(kdix,doWarmix,cbdix,yrix) = sqrt(sum((t.data(ix)-sq.data).^2));
           end; %if isempty(ix) else
         end; %for yrix

        end; %if doPlot

      end; %for cbdix

    end; %for kdix
  end; %for doWarmix

  if ( doPlot )
    stnm = upper(stn.station_name);

if (0)
    fmg; plot(squeeze(stn.optim_q0.error(:,1,:))); titlename([stnm ' Total Error: No Warm Layer']);
    legend(stn.optim_q0.cbdstrs);
    set(gca,'XTick',[1:length(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.error(:)),nanmax(stn.optim_q0.error(:))]);

    fmg; plot(squeeze(stn.optim_q0.error(:,2,:))); titlename([stnm ' Total Error: Warm Layer']);
    legend(stn.optim_q0.cbdstrs);
    set(gca,'XTick',[1:length(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.error(:)),nanmax(stn.optim_q0.error(:))]);

    fmg; plot(squeeze(nanmedian(stn.optim_q0.seaserror(:,1,:,:),4))); titlename([stnm ' Mdn Seas Err: No Warm Layer']);
    legend(stn.optim_q0.cbdstrs);
    set(gca,'XTick',[1:length(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.seaserror(:)),nanmax(stn.optim_q0.seaserror(:))]);

    fmg; plot(squeeze(nanmedian(stn.optim_q0.seaserror(:,2,:,:),4))); titlename([stnm ' Mdn Seas Err: Warm Layer']);
    legend(stn.optim_q0.cbdstrs);
    set(gca,'XTick',[1:length(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.seaserror(:)),nanmax(stn.optim_q0.seaserror(:))]);
end;

    fmg; plot_ts(stn.optim_q0.climt,squeeze(stn.optim_q0.climsq(:,1,default_cbdix))); titlename([stnm ' Daily Clim: No Warm Layer']);
    legend({'T',stn.optim_q0.kdstrs{:}});
    ylim([minmin([stn.optim_q0.climsq.data]),maxmax([stn.optim_q0.climsq.data])]);
    ylim([21,35]); appendtitlename([' (K_\theta=' num2str(kth,'%g') ')']);
    fmg; plot_ts(stn.optim_q0.climt,squeeze(stn.optim_q0.climsq(:,2,default_cbdix))); titlename([stnm ' Daily Clim: Warm Layer']);
    legend({'T',stn.optim_q0.kdstrs{:}});
    ylim([minmin([stn.optim_q0.climsq.data]),maxmax([stn.optim_q0.climsq.data])]);
    ylim([21,35]); appendtitlename([' (K_\theta=' num2str(kth,'%g') ')']);

if (0)
    fmg; plot_ts(stn.(sfld),squeeze(stn.optim_q0.sq(:,1,default_cbdix))); titlename([stnm ' Time Series: No Warm Layer']);
    legend({'T',stn.optim_q0.kdstrs{:}});
    ylim([minmin([stn.optim_q0.sq.data]),maxmax([stn.optim_q0.sq.data])]);
    fmg; plot_ts(stn.(sfld),squeeze(stn.optim_q0.sq(:,2,default_cbdix))); titlename([stnm ' Time Series: Warm Layer']);
    legend({'T',stn.optim_q0.kdstrs{:}});
    ylim([minmin([stn.optim_q0.sq.data]),maxmax([stn.optim_q0.sq.data])]);

    yrs = unique(get_year(stn.optim_q0.sq(1,1,default_cbdix).date));
    nrows = floor(sqrt(length(yrs)));
    ncols = ceil(length(yrs)/nrows);
    for doWarmix=1:2
      doWarm = logical(doWarmix-1);
      fmg;
      if (~doWarm)	suptitle([stnm ' Annual TS: No Warm Layer']);
      else		suptitle([stnm ' Annual TS: Warm Layer']);
      end;
      for yrix=1:length(yrs)
        yr = yrs(yrix);
        if ( is_valid_ts(stn.optim_q0.seassq(1,doWarmix,cbdix,yrix)) )
          ix = find(get_year(stn.(sfld).date)==yr);
          t.date = stn.(sfld).date(ix);
          t.data = stn.(sfld).data(ix);

          subplot_tight(nrows,ncols,yrix);
          plot_ts(t,stn.optim_q0.seassq(:,doWarmix,default_cbdix,yrix));
          datetick('x',17,'keeplimits');
          % legend({'T',stn.optim_q0.kdstrs{:}});
          xlabel(num2str(yr));
          xlim([min(stn.optim_q0.seassq(kdix,doWarmix,cbdix,yrix).date),...
                max(stn.optim_q0.seassq(kdix,doWarmix,cbdix,yrix).date)]);
          ylim([stn.optim_q0.minseassq,stn.optim_q0.maxseassq]);
          % ylim([minmin([stn.optim_q0.seassq.data]),maxmax([stn.optim_q0.seassq.data])]);
        end; %if isempty(ix) else
      end;
    end; %for doWarm
end;
  end;

  toc,
  timenow,
  set_more;

return;











%%%% DEBUG???
  % stn.(rawkl).data(:) = 0;




%%%% DEBUG???
firstix = find(stn.(bq0lpfld).date>=datenum(1993,9,3));
stn.(bdTflpfld).date = stn.(bq0lpfld).date(firstix:end);
stn.(bdTflpfld).data = stn.(bq0lpfld).data(firstix:end);
%%%% DEBUG???






  bdTlpfld = [bdTfld '_36_hour_lowpass'];



        %% Calculate and plot time series and errors from this option
        [tix,qix] = intersect_dates(stn.(sfld).date,stn.(bq0tfld).date);
        t.date = stn.(sfld).date(tix);	t.data = stn.(sfld).data(tix);
        q.date = stn.(bq0tfld).date(qix);	q.data = stn.(bq0tfld).data(qix);





   case {'sanf1'},
    % PRELIMINARY Results of running OPTIM_Q0
    opts.kd = get_opt(opts,'kd',[0.10,0.30,274]);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    opts.convective_drag_coefficient = get_opt(opts,'convective_drag_coefficient',3.8e-4);




   case {'smkf1'},
    % opts.kd = get_opt(opts,'kd',[0.10,0.25,274]);



    % *NOTE* *NOTE* *NOTE* : Warm Layer is better with *in situ* met data!





        % Evaluate climatological error for this option
        if ( ~isfield(stn.optim_q0,'climt') )
          [cum,tid] = grp_ts(t.data,t.date,'daily');
          stn.optim_q0.climt.data = cum;
          stn.optim_q0.climt.date = tid;
        end;
        [cum,tid] = grp_ts(q.data,q.date,'daily');
        stn.optim_q0.climq(kdix,doWarmix,cbdix).data = cum;
        stn.optim_q0.climq(kdix,doWarmix,cbdix).date = tid;
        t0 = stn.optim_q0.climt.data(1);
        sq.date = tid;
        sq.data = t0 + 24*(cumsum(cum) - cum(1));
        stn.optim_q0.climsq(kdix,doWarmix,cbdix).data = sq.data;
        stn.optim_q0.climsq(kdix,doWarmix,cbdix).date = sq.date;
        stn.optim_q0.climerror(kdix,doWarmix,cbdix) = sqrt(sum((stn.optim_q0.climt.data-sq.data).^2));









    fmg; plot_ts(stn.(sfld),squeeze(stn.optim_q0.sq(:,1,default_cbdix+1))); titlename([stnm ' Time Series: No Warm Layer Enhanced C_b_d']);
    legend({'T',stn.optim_q0.kdstrs{:}});
    fmg; plot_ts(stn.(sfld),squeeze(stn.optim_q0.sq(:,2,default_cbdix+1))); titlename([stnm ' Time Series: Warm Layer Enhanced C_b_d']);
    legend({'T',stn.optim_q0.kdstrs{:}});







  % kds = {0.1,0.2,0.4,0.6,[0.05,0.15,274],[0.10,0.30,274],[0.30,0.50,274],[0.50,0.70,274]};





   case {'mlrf1'},
    % Results of running OPTIM_Q0
    opts.kd = get_opt(opts,'kd',[0.10,0.30,274]);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    opts.convective_drag_coefficient = get_opt(opts,'convective_drag_coefficient',3.8e-4);





    fmg; plot_ts(stn.(sfld),squeeze(stn.optim_q0.sq(:,1,default_cbdix+2))); titlename([stnm ' Time Series: No Warm Layer EnHANced C_b_d']);
    legend({'T',stn.optim_q0.kdstrs{:}});
    fmg; plot_ts(stn.(sfld),squeeze(stn.optim_q0.sq(:,2,default_cbdix+2))); titlename([stnm ' Time Series: Warm Layer EnHANced C_b_d']);
    legend({'T',stn.optim_q0.kdstrs{:}});





        %% Benthic Heat Exchanges
        % stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld,opts);
        stn = station_benthic_exchange(stn,sfld,ssufld,ssvfld,qbfld,btfld,qbofld,opts);
        % badix = find(~isfinite(stn.(qbofld).data) | abs(stn.(qbofld).data)>2e3);
        badix = find(~isfinite(stn.(qbofld).data));
        stn.(qbofld).date(badix) = [];
        stn.(qbofld).data(badix) = [];
        stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],mhfld);





    stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                            pfld,sfld,[],[],TURPFX,dsrfld,dlrfld,rfld,...
                            Dfld,ssufld,ssvfld,wpfld,whfld,pblzfld,doWarm,max_wl);



  default_cbd = 3.8e-3; %LONF1



  % (BULK or REANALYSIS for long-wave component of net radiation?)
  %qradfld = [RAPFX '_arf'];
  qradfld = [SPREFIX RAPFX '_' ISPFX '_arf'];
  qradtfld = [qradfld '_term'];
  qturfld = [TURPFX '_turbulent_heat_flux'];
  qturtfld = [qturfld '_term'];

  q0fld = [TURPFX '_net_heat_flux'];
  qtfld = [q0fld '_term'];

  % Fluxes without warm-layer adjustment
  qlh30fld = [TUR30PFX '_latent_heat_flux'];
  qsh30fld = [TUR30PFX '_sensible_heat_flux'];
  qrh30fld = [TUR30PFX '_rain_heat_flux'];
  qtur30fld = [TUR30PFX '_turbulent_heat_flux'];
  q030fld = [TUR30PFX '_net_heat_flux'];
  qt30fld = [q030fld '_term'];

  % Fluxes without absorption correction
  % (BULK or REANALYSIS for long-wave component of net radiation?)
  sqradfld = [SPREFIX RAPFX '_' ISPFX '_rf'];
  %sqradfld = [RAPFX '_rf'];
  sq0fld = ['simple_' q0fld];
  sqtfld = ['simple_' qtfld];
  % ... AND without warm-layer adjustment
  sq030fld = ['simple_' q030fld];
  sqt30fld = ['simple_' qt30fld];










   case {'lonf1'},
    opts.kd = get_opt(opts,'kd',0.40);
    % opts.kd = get_opt(opts,'kd',[0.10,0.25,274]);
    % opts.kd = get_opt(opts,'kd',[0.20,0.40,91]);
    % opts.kd = get_opt(opts,'kd',[0.30,0.50,274]);
    % opts.kd = get_opt(opts,'kd',[0.40,0.70,274]);
    opts.grid_interp_method = get_opt(opts,'grid_interp_method','triangular,linear');





function stn = optim_q0(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld)
%function stn = optim_q0(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld)

  set_more off
  timenow,
  tic,

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

  datapath = get_thesis_path('../data');


  stn = get_station_from_station_name(stn_or_stnm);

  opts = station_heat_budget_options(stn);

  if ( ~isfield(stn,bathyfld) )
    stn = get_ngdc_bathy_station(stn);
  end;
  if ( ~isfield(stn,slopefld) )
    stn = station_ngdc_offshore_slope(stn);
  end;
  if ( ~isfield(stn,bathorifld) )
    stn = station_optimal_isobath_orientation(stn);
  end;

  stn = station_tmd_tide(stn);
  stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
  % Assume max warm-layer depth somewhere near bottom boundary layer top
  max_wl = nanmax(stn.(mhfld).data) - 0.5;

  warning('off','Ecoforecasts:mergedNonTS');
  stn = load_all_ndbc_data(stn);
  if ( ~isfield(stn,afld) || ~isfield(stn,sfld) )
    switch (ISPFX),
     case 'ndbc',	%stn = load_all_ndbc_data(stn);
     case 'erai',	stn = get_erai_station(stn);
     case 'icon',	stn = load_station_data(stn);
     otherwise,		error('Unknown in situ dataset "%s"',ISPFX);
    end;
  end;
  warning('on','Ecoforecasts:mergedNonTS');

  stn = verify_variable(stn, afld);
  stn = verify_variable(stn, pfld);
  stn = verify_variable(stn, Wfld);
  stn = verify_variable(stn, sfld);

  if ( ~isfield(stn,dsrfld) || ~isfield(stn,cfld) )
    switch (RAPFX),
     case 'erai',		stn = get_erai_station(stn);
     case 'ncep',		stn = get_ncep_station(stn,'narr');
     %%% Not Yet Implemented
     % case 'era40',		stn = get_era40_station(stn);
     % case 'cfsr',		stn = get_ncep_station(stn,'cfsr');
     otherwise,		error('Unavailable gridded/reanalysis dataset "%s"',RAPFX);
    end;
    stn = station_heat_flux_term(stn,raq0fld,raqtfld,sfld,[],nanmean(stn.(mhfld).data));
  end;

  if ( ~isfield(stn,whfld) )
    % If waves not from reanalysis, user must want model (WaveWatch III) or wind estimate
    switch (WAVEPFX),
     case 'ww3',	stn = get_ww3_station(stn);
     case 'ndbc',	stn = station_wind_to_wave(stn,Wfld,Dfld,wpfld,whfld,wdfld);
     case 'erai',	stn = get_erai_station(stn); %IF NOT ALSO OUR REANALYSIS DATASET
     otherwise,		error('Unknown wave source "%s"',WAVEPFX);
    end;
  end;

  %% Low-pass filter winds for quasi-Eulerian currents
  stn = verify_variable(stn,Ulpfld);
  stn = verify_variable(stn,Vlpfld);
  stn.(Wlpfld).date = stn.(Ulpfld).date;
  stn.(Wlpfld).data = uv_to_spd(stn.(Ulpfld).data,stn.(Vlpfld).data);
  stn.(Dlpfld).date = stn.(Ulpfld).date;
  stn.(Dlpfld).data = uv_to_dir(stn.(Ulpfld).data,stn.(Vlpfld).data);

  stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wlpfld,Dlpfld,whfld,wpfld,wdfld);
  stn.(netufld) = ts_op(stn.(tufld),stn.(ssufld),'+');
  stn.(netvfld) = ts_op(stn.(tvfld),stn.(ssvfld),'+');


  stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);


  % kds = {0.1,0.2,0.4,0.6,[0.05,0.15,274],[0.10,0.30,274],[0.30,0.50,274],[0.50,0.70,274]};
  kds = {0.2,0.4,[0.10,0.30,91],[0.30,0.50,91],[0.10,0.30,182],[0.30,0.50,182],[0.10,0.30,274],[0.30,0.50,274],[0.10,0.30,0],[0.30,0.50,0]};
  cbds = {3.8e-4};
  kths = {2,20,200,[1,3,1],[10,30,1],[100,300,1]};

  [ig,default_cbdix] = min(abs(cbds - 3.8e-4));

  for ix=1:length(kds)
    stn.optim_q0.kdstrs{ix} = num2str(kds{ix},'%g,');
  end;
  for ix=1:length(cbds)
    stn.optim_q0.cbdstrs{ix} = num2str(cbds{ix},'%g,');
  end;

  for doWarmix=1:2
    doWarm = logical(doWarmix-1);

    % tic,

    % Fluxes WITH or WITHOUT warm-layer adjustment
    % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
    %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
    %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);
    stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                            pfld,sfld,[],[],TURPFX,dsrfld,dlrfld,rfld,...
                            Dfld,netufld,netvfld,wpfld,whfld,pblzfld,doWarm,max_wl);
    % Algorithm sometimes returns complex numbers!
    stn.(qlhfld).data = real(stn.(qlhfld).data);
    stn.(qshfld).data = real(stn.(qshfld).data);
    stn.(qrhfld).data = real(stn.(qrhfld).data);
    stn.(qturfld) = ts_op(stn.(qlhfld),stn.(qshfld),'+');
    if ( isfield(stn,qrhfld) && is_valid_ts(stn.(qrhfld)) )
      stn.(qturfld) = ts_op(stn.(qturfld),stn.(qrhfld),'+');
    end;
    stn = station_heat_flux_term(stn,qturfld,qturtfld,sfld,[],mhfld);

    % toc,

    % Just in case something above reset it!
    more off;

    for kdix=1:length(kds)
      kd = kds{kdix};

      opts.kd = kd;
      stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,opts);

      stn.(qradfld) = ts_op(stn.(asrfld),stn.(lrfld),'+');
      stn = station_heat_flux_term(stn,qradfld,qradtfld,sfld,[],mhfld);

      stn.(q0fld) = ts_op(stn.(qradfld),stn.(qturfld),'+');
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],mhfld);

      % Net flux without absorption calculation or benthic flux - for comparison
      stn.(sqradfld) = ts_op(stn.(srfld),stn.(lrfld),'+');
      stn.(sq0fld) = ts_op(stn.(sqradfld),stn.(qturfld),'+');
      stn = station_heat_flux_term(stn,sq0fld,sqtfld,sfld,[],mhfld);

      for cbdix=1:length(cbds)
        cbd = cbds{cbdix};
        opts.benthic_debug = false;
        opts.convective_drag_coefficient = cbd;

        %% Benthic Heat Exchanges
        stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld,opts);
        stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],mhfld);
        stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'-');
        stn = station_heat_flux_term(stn,bq0fld,bq0tfld,sfld,[],mhfld);

        for kthix=1:length(kths)
          kth = kths{kthix};
          opts.K_theta = get_opt(opts,'K_theta',kth);

          

          %% Calculate and plot time series and errors from this option
          [tix,qix] = intersect_dates(stn.(sfld).date,stn.(bq0tfld).date);
          t.date = stn.(sfld).date(tix);	t.data = stn.(sfld).data(tix);
          q.date = stn.(bq0tfld).date(qix);	q.data = stn.(bq0tfld).data(qix);

          % Evaluate total error for this option
          t0 = t.data(1);
          sq.date = q.date;
          sq.data = t0 + cumsum(q.data) - q.data(1);
          stn.optim_q0.error(kdix,doWarmix,cbdix,kthix) = sqrt(sum((t.data-sq.data).^2));

          stn.optim_q0.doWarm(kdix,doWarmix,cbdix,kthix) = doWarm;
          stn.optim_q0.kd{kdix,doWarmix,cbdix,kthix} = kd;
          stn.optim_q0.cbd(kdix,doWarmix,cbdix,kthix) = cbd;
          stn.optim_q0.kth(kdix,doWarmix,cbdix,kthix) = kth;
          stn.optim_q0.q(kdix,doWarmix,cbdix,kthix) = q;
          stn.optim_q0.sq(kdix,doWarmix,cbdix,kthix) = sq;

          % Evaluate seasonal amplitude error for this option
          yrs = unique(get_year(t.date));
          for yrix = 1:length(yrs)
            yr = yrs(yrix);
            ix = find(get_year(t.date)==yr);
            if ( isempty(ix) )
              warning('NO DATA for year %g??',yr);
            else
              t0 = t.data(ix(1));
              sq.date = q.date(ix);
              sq.data = t0 + cumsum(q.data(ix)) - q.data(ix(1));
              stn.optim_q0.seasyear(kdix,doWarmix,cbdix,kthix,yrix) = yr;
              stn.optim_q0.seaserror(kdix,doWarmix,cbdix,kthix,yrix) = sqrt(sum((t.data(ix)-sq.data).^2));
            end; %if isempty(ix) else
          end; %for yrix
        end; %for kthix

      end; %for cbdix

    end; %for kdix
  end; %for doWarmix

  doPlot = true;
  if ( doPlot )
    stnm = upper(stn.station_name);

    fmg; plot(squeeze(stn.optim_q0.error(:,1,:))); titlename([stnm ' Total Error: No Warm Layer']);
    legend(stn.optim_q0.cbdstrs);
    set(gca,'XTick',[1:length(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.error(:)),nanmax(stn.optim_q0.error(:))]);

    fmg; plot(squeeze(stn.optim_q0.error(:,2,:))); titlename([stnm ' Total Error: Warm Layer']);
    legend(stn.optim_q0.cbdstrs);
    set(gca,'XTick',[1:length(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.error(:)),nanmax(stn.optim_q0.error(:))]);

    fmg; plot(squeeze(nanmedian(stn.optim_q0.seaserror(:,1,:,:),4))); titlename([stnm ' Mdn Seas Err: No Warm Layer']);
    legend(stn.optim_q0.cbdstrs);
    set(gca,'XTick',[1:length(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.seaserror(:)),nanmax(stn.optim_q0.seaserror(:))]);

    fmg; plot(squeeze(nanmedian(stn.optim_q0.seaserror(:,2,:,:),4))); titlename([stnm ' Mdn Seas Err: Warm Layer']);
    legend(stn.optim_q0.cbdstrs);
    set(gca,'XTick',[1:length(kds)],'XTickLabel',stn.optim_q0.kdstrs);
    ylim([nanmin(stn.optim_q0.seaserror(:)),nanmax(stn.optim_q0.seaserror(:))]);

    fmg; plot_ts(stn.ndbc_sea_t,squeeze(stn.optim_q0.sq(:,1,default_cbdix))); titlename([stnm ' Time Series: No Warm Layer']);
    legend(stn.optim_q0.kdstrs);
    fmg; plot_ts(stn.ndbc_sea_t,squeeze(stn.optim_q0.sq(:,2,default_cbdix))); titlename([stnm ' Time Series: Warm Layer']);
    legend(stn.optim_q0.kdstrs);
  end;

  toc,
  timenow,
  set_more;

return;





              stn.optim_q0.seast(kdix,doWarmix,cbdix,kthix,yrix).date = t.date(ix);
              stn.optim_q0.seast(kdix,doWarmix,cbdix,kthix,yrix).data = t.data(ix);







      for cbdix=1:length(cbds)
        cbd = cbds{cbdix};
        opts.benthic_debug = false;
        opts.convective_drag_coefficient = cbd;

        %% Benthic Heat Exchanges
        stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld,opts);
        stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],mhfld);
        stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'-');
        stn = station_heat_flux_term(stn,bq0fld,bq0tfld,sfld,[],mhfld);

        %% Calculate and plot time series and errors from this option
        [tix,qix] = intersect_dates(stn.(sfld).date,stn.(bq0tfld).date);
        t.date = stn.(sfld).date(tix);	t.data = stn.(sfld).data(tix);
        q.date = stn.(bq0tfld).date(qix);	q.data = stn.(bq0tfld).data(qix);

        % Evaluate total error for this option
        t0 = t.data(1);
        sq.date = q.date;
        sq.data = t0 + cumsum(q.data) - q.data(1);
        stn.optim_q0.error(kdix,doWarmix,cbdix) = sqrt(sum((t.data-sq.data).^2));

        stn.optim_q0.doWarm(kdix,doWarmix,cbdix) = doWarm;
        stn.optim_q0.kd{kdix,doWarmix,cbdix} = kd;
        stn.optim_q0.cbd(kdix,doWarmix,cbdix) = cbd;
        stn.optim_q0.q(kdix,doWarmix,cbdix) = q;
        stn.optim_q0.sq(kdix,doWarmix,cbdix) = sq;

        % Evaluate seasonal amplitude error for this option
        yrs = unique(get_year(t.date));
        for yrix = 1:length(yrs)
          yr = yrs(yrix);
          ix = find(get_year(t.date)==yr);
          if ( isempty(ix) )
            warning('NO DATA for year %g??',yr);
          else
            t0 = t.data(ix(1));
            sq.date = q.date(ix);
            sq.data = t0 + cumsum(q.data(ix)) - q.data(ix(1));
            stn.optim_q0.seasyear(kdix,doWarmix,cbdix,yrix) = yr;
            stn.optim_q0.seast(kdix,doWarmix,cbdix,yrix).date = t.date(ix);
            stn.optim_q0.seast(kdix,doWarmix,cbdix,yrix).data = t.data(ix);
            stn.optim_q0.seaserror(kdix,doWarmix,cbdix,yrix) = sqrt(sum((t.data(ix)-sq.data).^2));
          end; %if isempty(ix) else
        end; %for yrix

      end; %for cbdix





  cbds = {0,3.8e-5,3.8e-4,3.8e-3};




   case {'sanf1'},
    % opts.kd = get_opt(opts,'kd',[0.10,0.25,274]);




   case {'fwyf1'},
    opts.kd = get_opt(opts,'kd',[0.30,0.70,274]);
    % opts.hc_scaling = get_opt(opts,'hc_scaling','US');





  cbds = {0,0.5e-4,1.0e-4,2.0e-4,3.8e-4,4.0e-4,8.0e-4};




   case {'mlrf1'},
    % % opts.kd = get_opt(opts,'kd',[0.10,0.25,274]);
    % opts.kd = get_opt(opts,'kd',[0.10,0.40,274]);
    % % opts.kd = get_opt(opts,'kd',[0.10,0.4,92]);
    % % opts.kd = get_opt(opts,'kd',0.25);








%function stn = optim_q0(stn_or_stnm)





  kds = {0.2,0.4,[0.10,0.30,274],[0.10,0.30,91],[0.30,0.50,274]};
  cbds = {0,0.5e-4,1.0e-4,2.0e-4,3.8e-4,8.0e-4,16.0e-4};







  doPlot = true;
  if ( doPlot )
    fmg; plot(squeeze(stn.optim_q0.error(:,1,:))); titlename('No warm layer'); legend(num2str([squeeze(stn.optim_q0.cbd(3,1,:))]));
    set(gca,'XTick',[1:length(kds)],'XTickLabel',kdstrs);
    ylim([min(stn.optim_q0.error(:)),max(stn.optim_q0.error(:))]);

    fmg; plot(squeeze(stn.optim_q0.error(:,2,:))); titlename('Warm layer'); legend(num2str([squeeze(stn.optim_q0.cbd(3,1,:))]));
    set(gca,'XTick',[1:length(kds)],'XTickLabel',kdstrs);
    ylim([min(stn.optim_q0.error(:)),max(stn.optim_q0.error(:))]);

    fmg; plot(squeeze(mean(stn.optim_q0.seaserror(:,1,:,:),4))); titlename('No warm layer SEASONAL MEAN'); legend(num2str([squeeze(stn.optim_q0.cbd(3,1,:))]));
    set(gca,'XTick',[1:length(kds)],'XTickLabel',kdstrs);
    ylim([min(stn.optim_q0.seaserror(:)),max(stn.optim_q0.seaserror(:))]);

    fmg; plot(squeeze(mean(stn.optim_q0.seaserror(:,2,:,:),4))); titlename('Warm layer SEASONAL MEAN'); legend(num2str([squeeze(stn.optim_q0.cbd(3,1,:))]));
    set(gca,'XTick',[1:length(kds)],'XTickLabel',kdstrs);
    ylim([min(stn.optim_q0.seaserror(:)),max(stn.optim_q0.seaserror(:))]);
  end;






            stn.optim_q0.seasbq0t(kdix,doWarmix,cbdix,yrix).date = q.date(ix);
            stn.optim_q0.seasbq0t(kdix,doWarmix,cbdix,yrix).data = sq;





    % We want rain rate in mm/hr, not m/hr!
    result(stix).raw_erai_conv_precip.data = result(stix).raw_erai_conv_precip.data*1e3;
    result(stix).raw_erai_precip.data = result(stix).raw_erai_precip.data*1e3;
    result(stix).erai_conv_precip.data = result(stix).erai_conv_precip.data*1e3;
    result(stix).erai_precip.data = result(stix).erai_precip.data*1e3;





function stn = optim_q0(stn_or_stnm)
%function stn = optim_q0(stn_or_stnm)

  set_more off
  timenow,
  % tic,

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;

  datapath = get_thesis_path('../data');


  stn = get_station_from_station_name(stn_or_stnm);

  opts = station_heat_budget_options(stn);

  stn = load_all_ndbc_data(stn);
  stn = get_erai_station(stn);

  stn = verify_variable(stn, afld);
  stn = verify_variable(stn, pfld);
  stn = verify_variable(stn, Wfld);
  stn = verify_variable(stn, sfld);

  if ( ~isfield(stn,bathyfld) )
    stn = get_ngdc_bathy_station(stn);
  end;
  if ( ~isfield(stn,slopefld) )
    stn = station_ngdc_offshore_slope(stn);
  end;
  if ( ~isfield(stn,bathorifld) )
    stn = station_optimal_isobath_orientation(stn);
  end;

  stn = station_tmd_tide(stn);
  stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
  % Assume max warm-layer depth somewhere near bottom boundary layer top
  max_wl = nanmax(stn.(mhfld).data) - 0.5;

  %% Low-pass filter winds for quasi-Eulerian currents
  stn = verify_variable(stn,Ulpfld);
  stn = verify_variable(stn,Vlpfld);
  stn.(Wlpfld).date = stn.(Ulpfld).date;
  stn.(Wlpfld).data = uv_to_spd(stn.(Ulpfld).data,stn.(Vlpfld).data);
  stn.(Dlpfld).date = stn.(Ulpfld).date;
  stn.(Dlpfld).data = uv_to_dir(stn.(Ulpfld).data,stn.(Vlpfld).data);

  stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wlpfld,Dlpfld,whfld,wpfld,wdfld);
  stn.(netufld) = ts_op(stn.(tufld),stn.(ssufld),'+');
  stn.(netvfld) = ts_op(stn.(tvfld),stn.(ssvfld),'+');


  stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);


  % kds = {0.2};
  kds = {0.2,0.4,[0.10,0.30,274]};
  % cbds = {3.8e-4};
  cbds = {2.0e-4,3.8e-4,8.0e-3};

  for kdix=1:length(kds)
    opts.kd = kds{kdix};

    stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,opts);

    for doWarmix=1:2
      doWarm = logical(doWarmix-1);

      tic,

      % Fluxes WITH or WITHOUT warm-layer adjustment
      % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
      %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
      %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,doWarm,max_wl);
      stn = station_heat_flux(stn,Wfld,afld,rhfld,...
                              pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
                              Dfld,netufld,netvfld,wpfld,whfld,pblzfld,doWarm,max_wl);
      % Algorithm sometimes returns complex numbers!
      stn.(qlhfld).data = real(stn.(qlhfld).data);
      stn.(qshfld).data = real(stn.(qshfld).data);
      stn.(qrhfld).data = real(stn.(qrhfld).data);
      stn.(q0fld).data = real(stn.(q0fld).data);
      stn.(qradfld) = ts_op(stn.(asrfld),stn.(lrfld),'+');
      stn.(qturfld) = ts_op(stn.(qlhfld),stn.(qshfld),'+');
      stn = station_heat_flux_term(stn,qradfld,qradtfld,sfld,[],mhfld);
      stn = station_heat_flux_term(stn,qturfld,qturtfld,sfld,[],mhfld);
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],mhfld);

      % Net flux without absorption calculation or benthic flux - for comparison
      stn.(sqradfld) = ts_op(stn.(srfld),stn.(lrfld),'+');
      stn.(sq0fld) = ts_op(stn.(sqradfld),stn.(qturfld),'+');
      stn = station_heat_flux_term(stn,sq0fld,sqtfld,sfld,[],mhfld);

      toc,


      for cbdix=1:length(cbds)
        opts.convective_drag_coefficient = cbds{cbdix};

        %% Benthic Heat Exchanges
        stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld,opts);
        stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],mhfld);
        stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'-');
        stn = station_heat_flux_term(stn,bq0fld,bq0tfld,sfld,[],mhfld);

        [tix,qix] = intersect_dates(stn.(sfld).date,stn.(bq0tfld).date);
        t.date = stn.(sfld).date(tix);		t.data = stn.(sfld).data(tix);
        q.date = stn.(bq0tfld).date(qix);	q.data = stn.(bq0tfld).data(qix);

        % Evaluate total error for this option
        t0 = t.data(1);
        sq = t0 + cumsum(q.data) - q.data(1);
        stn.optim_q0.error(kdix,doWarmix,cbdix) = sqrt(sum((t.data-sq).^2));

        % Evaluate seasonal amplitude error for this option
        yrs = unique(get_year(t.date));
        for yrix = 1:length(yrs)
          yr = yrs(yrix);
          ix = find(get_year(t.date)==yr);
          if ( isempty(ix) )
            warning('NO DATA for year %g??',yr);
          else
            t0 = t.data(ix(1));
            sq = t0 + cumsum(q.data(ix)) - q.data(ix(1));
            stn.optim_q0.seasyear(kdix,doWarmix,cbdix,yrix) = yr;
            stn.optim_q0.seast(kdix,doWarmix,cbdix,yrix).date = t.date(ix);
            stn.optim_q0.seast(kdix,doWarmix,cbdix,yrix).data = t.data(ix);
            stn.optim_q0.seasbq0t(kdix,doWarmix,cbdix,yrix).date = q.date(ix);
            stn.optim_q0.seasbq0t(kdix,doWarmix,cbdix,yrix).data = sq;
            stn.optim_q0.seaserror(kdix,doWarmix,cbdix,yrix) = sqrt(sum((t.data(ix)-sq).^2));
          end; %if isempty(ix) else
        end; %for yrix
      end; %for cbdix

    end; %for doWarmix
  end; %for kdix

  % toc,
  timenow,
  set_more;

return;









    % stn.q0 = stn.(q0fld);
    % qstr=q0fld;
    % stn.q0 = stn.(sq0fld);
    % qstr=sq0fld;
    % stn.q0 = stn.(dTffld);
    % qstr=dTffld;
    stn.q0 = stn.(bdTffld);
    qstr=bdTffld;





   case {13,14,15},
    nrows = 5;
    ncols = 3;
   case {16},
    nrows = 4;
    ncols = 4;






  % qf = bdTfld;







    % qstr=bdTffld;

    %%%% ??? DEBUG: Base horizontal convection on surface heating only!
    % q=bq0; stn.commentstr = [stn.commentstr ' HC(Q0+Qb) '];
    %%%% ??? DEBUG: Base horizontal convection on total (km-scale) budget
    q=dT;








  % QEPFX - Quasi-Eulerian currents and heat advection
  if ( ~exist('QEPFX','var') || isempty(QEPFX) )
    switch (KMPFX),
     case 'fkeys_hycom',	QEPFX = [WAVEPFX '_fkeys_qe'];
     case 'gom_hycom',		QEPFX = [WAVEPFX '_gom_qe'];
     case 'avhrr_weekly_sst',	QEPFX = [WAVEPFX '_avhrr_qe'];
     case 'none',		QEPFX = [WAVEPFX '_none_qe'];
     otherwise,			error('Unknown km-scale model "%s"',KMPFX);
    end;
  end;


  % In case any of our PREFIX argument strings are too long!
  fix_varnamelengths;


  %%%
  %% All station struct fieldnames used to produce heat budget 






    if ( ~isfield(stn,mhfld) )
      % Calculate the mean tide depth experienced by a watermass moving over
      % an M2 tidal ellipse centered on the coordinates of our station
      % stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
      stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld,tufld,tvfld);
%%%% DEBUG???
stn.mean_tmd_tide_i_depth.data = stn.mean_tmd_tide_i_depth.data + 5;
%%%% DEBUG???
    end;






  bhcdTdt = ['benthic_' HCPFX '_dTdt'];



    switch (KMPFX),
     case 'avhrr_weekly_sst',
      opts.km_scale_advection = false;
      opts.gradient_climatology = [];
      opts.laplacian_climatology = [];
      disp('ONLY STOKES, AVHRR_WEEKLY_SST GRADIENTS, and FIELD LAPLACIAN');
     case 'none',
      opts.km_scale_advection = false;
      disp('ONLY STOKES');
    end;







    if ( ~isfield(stn,ufld) )
      switch (KMPFX),
        %function stn = get_fkeys_hycom(stn_or_stnm,mindt,maxdt,vars,flds,interpMethod,fkeyspath)
       case 'fkeys_hycom',		stn = get_fkeys_hycom(stn,[],[],[],[],'linear');
        %function stn = get_gom_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
       case 'gom_hycom',		stn = get_gom_hycom(stn);
       case 'avhrr_weekly_sst',
        if ( ~isfield(stn,Tfld) )
          disp('Loading AVHRR_WEEKLY_SST instead of hydrodynamic model data...');
          stn = get_avhrr_weekly_field(stn,true);
          opts.km_scale_advection = false;
          opts.gradient_climatology = [];
          opts.laplacian_climatology = [];
          disp('ONLY STOKES, AVHRR_WEEKLY_SST GRADIENTS, and FIELD LAPLACIAN');
        end;
       case 'none',
        disp('Loading NO hydrodynamic model data...');
        opts.km_scale_advection = false;
        disp('ONLY STOKES');
       otherwise,			error('Unknown km-scale data source "%s"',KMPFX);
      end;
      more off;
    end;







    %%%
    %% Sea temperature (in situ NDBC, site-specific, or remote sensing)

    if ( ~isfield(stn,sfld) )
      if ( regexp(sfld,'hourly_misst') )
        stn = get_misst_station(stn);
        stn.hourly_misst_sst = interp_ts(stn.misst_sst);
      elseif ( regexp(sfld,'avhrr_weekly') )
        stn = get_avhrr_weekly_field(stn,true);
        stn.avhrr_weekly_sst.date = stn.avhrr_weekly_sst_field.date;
        stn.avhrr_weekly_sst.data = ...
            interp_field(stn.avhrr_weekly_sst_field.lat,stn.avhrr_weekly_sst_field.lon,...
                         stn.avhrr_weekly_sst_field.field,stn.lat,stn.lon);
        stn.hourly_avhrr_weekly_sst = interp_ts(stn.avhrr_weekly_sst);







    %%%
    %% Sea temperature (in situ NDBC, site-specific, or remote sensing)

    if ( ~isfield(stn,sfld) )
      switch (sfld),
       case 'hourly_misst_sst',
        stn = get_misst_station(stn);
        stn.hourly_misst_sst = interp_ts(stn.misst_sst);

       case 'hourly_avhrr_weekly_sst',
        stn = get_avhrr_weekly_field(stn,true);
        stn.avhrr_weekly_sst.date = stn.avhrr_weekly_sst_field.date;
        stn.avhrr_weekly_sst.data = ...
            interp_field(stn.avhrr_weekly_sst_field.lat,stn.avhrr_weekly_sst_field.lon,...
                         stn.avhrr_weekly_sst_field.field,stn.lat,stn.lon);
        stn.hourly_avhrr_weekly_sst = interp_ts(stn.avhrr_weekly_sst);

       case 'microcat_seatemp',
        stn = get_looe1_microcat(stn);
       case 'adcp_seatemp',
        stn = get_looe1_adcp(stn);

       case 'ndbc_sea_t',
        if ( ~strcmpi(ISPFX,'ndbc') )
          stn = load_all_ndbc_data(stn);
        end;
       case {'sea_t','ct_shallow_seatemp','ctd_shallow_seatemp','ctd_deep_seatemp',},
        if ( ~strcmpi(ISPFX,'icon') )
          stn = load_station_data(stn);
        end;
       case 'erai_sea_t',
        if ( ~strcmpi(ISPFX,'erai') && ~strcmpi(RAPFX,'erai') )
          stn = get_erai_station(stn);
        end;

       otherwise,
        error(['Unknown sea temperature field ',sfld]);
      end;
    end;






%%%% ??? DEBUG
% mindt = datenum(yr,mos(1),7,0,0,0);
% maxdt = datenum(yr,mos(end)+1,21,0,0,0);
%%%% ??? DEBUG



%%%% ??? DEBUG
% plot(stn.hourly_misst_sst.date,stn.hourly_misst_sst.data,'r')
% ylim([16,32]);
%%%% ??? DEBUG





sssst = 'stokes_avhrr_weekly_heat_advection';
if (1)
    begyr=1998; endyr=1998; mos=1:5;
    annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,sfld,stn.commentstr,mos,begyr,endyr);
  if (0)
    if ( ~isfield(stn,'hourly_avhrr_weekly_sst_xshore') )
      stn = get_avhrr_weekly_field(stn,true);
      stn = station_reorient_vectors(stn,bathorifld,'hourly_avhrr_weekly_sst_x','hourly_avhrr_weekly_sst_y');
    end;
    if ( ~isfield(stn,sssst) )
      %% Use Stokes drift o del AVHRR field for heat advection
      % ax = ts_op(stn.(ssufld),stn.hourly_avhrr_weekly_sst_x,'*');
      % ay = ts_op(stn.(ssvfld),stn.hourly_avhrr_weekly_sst_y,'*');
      % stn.(sssst) = ts_op(ax,ay,'+');
      % ax = []; ay = []; clear ax ay

      % Use cross shore components only of Stokes drift and del AVHRR field for heat advection
      if ( ~isfield(stn,ssxsfld) )
        stn = station_reorient_vectors(stn,bathorifld,ssufld,ssvfld,ssxsfld,sslsfld);
      end;

      % % % % Allow heat advection both on- and offshore
      % % % stn.(sssst) = ts_op(stn.(ssxsfld),stn.hourly_avhrr_weekly_sst_xshore,'*');
      % % % Allow heat advection during *onshore Stokes drift only*
      % % stn.(sssst) = ts_op(stn.(ssxsfld),stn.hourly_avhrr_weekly_sst_xshore,'*',@(x)(find(x.data<=0)));
      % % Significantly moderate heat advection during *offshore Stokes drift only*
      % stn.(sssst) = ts_op(stn.(ssxsfld),stn.hourly_avhrr_weekly_sst_xshore,'*',@(x)(find(x.data<=0.05)));
      % Significantly moderate heat advection
      stn.(ssxsfld).data = stn.(ssxsfld).data .* 0.5;
      stn.(sssst) = ts_op(stn.(ssxsfld),stn.hourly_avhrr_weekly_sst_xshore,'*',@(x)(find(x.data<=0.05)));

      stn.(sssst).data = -stn.(sssst).data.*3600;

      stn.dTdt = ts_op(stn.(dTfld),stn.(sssst),'+');
    end;
    T0 = stn.(sfld).data(find(stn.(sfld).date>=datenum(begyr,mos(1),1),1));
    begix = find(stn.dTdt.date>=datenum(begyr,mos(1),1),1);
    plot(stn.dTdt.date(begix:end),...
         T0+cumsum(stn.dTdt.data(begix:end)),'r-.');
  end;
end;





  %rgnlat = [minlat:dlat:maxlat]';
  rgnlat = [maxlat:-dlat:minlat]';



    station.avhrr_sst_midy = interp1(rgnlat,[length(rgnlat):-1:1]',stn.lat,'nearest');




%DEBUG:
    station.avhrr_sst_field.lat = rgnlat(station.avhrr_sst_yix);
    disp(['RE-Saving ' matfname]);
    save(matfname,'station');
%DEBUG:





%%%% DEBUG??? REALLY BELONGS UP ABOVE BEFORE "save(...)"!
    % Field .date should be a column vector (Nx1)
    stn.avhrr_sst_field.date = stn.avhrr_sst_field.date(:);
    % Make sure our dates increase monotonically
    [stn.avhrr_sst_field.date,sortix] = sort(stn.avhrr_sst_field.date);
    stn.avhrr_sst_field.field = stn.avhrr_sst_field.field(sortix,:,:);








    station = safe_rmfield(station,{'avhrr_sst_field','avhrr_sst'});






          sstbytes_all = imread(fpath);

          [ig,ig,YYYY,MM,DD,hh,mm] = parseusfurl(url);

          sstbytes = sstbytes_all(station.avhrr_sst_xix,station.avhrr_sst_yix);

          sst = read_avhrr_subset(sstbytes);

          % Only take images with some cloud-free data surrounding our site!
          ctr = round(size(sst)./2);
          ctrsst = sst((ctr(1)-8+1):(ctr(1)+8),(ctr(2)-8+1):(ctr(2)+8));
          % if ( length(find(isnan(ctrsst))) < (0.50*numel(ctrsst)) && ...
          %      ~isnan(sst(ctr(1), ctr(2))) )
          if ( length(find(isnan(ctrsst))) < (0.50*numel(ctrsst)) )
            dts{boxix}{sstix{boxix}} = datenum(YYYY,MM,DD,hh,mm,0);
            sstix{boxix} = sstix{boxix} + 1;
          else
            skipped{boxix} = skipped{boxix} + 1;
            sstbytes = cast(0, 'uint8');
          end;

          % Try to save STATION-specific image subset for future reference.
          % NOTE: If data in image was unuseable, 'sstbytes' will just be [0].
          if ( numel(sstbytes) <= 1 || ~exist(pname, 'file') )
            imwrite(sstbytes, pname);
          end;

          clear sst;
          clear ctrsst;
          clear sstbytes;









  for yrix=1:length(yrs)
    yr = yrs(yrix);
    mindt = datenum(yr,mos(1),1);
    maxdt = datenum(yr,mos(end)+1,1);
%%%% ??? DEBUG
maxdt = datenum(yr,mos(end)+2,1);

    % dtix = find(get_year(stn.(comparisonsfld).date)==yr & ismember(get_month(stn.(comparisonsfld).date),[mos]));
    dtix = find(mindt<=stn.(comparisonsfld).date & stn.(comparisonsfld).date<=maxdt);
    t.date = stn.(comparisonsfld).date(dtix);
    t.data = stn.(comparisonsfld).data(dtix);

    % dtix = find(get_year(stn.(qfld).date)==yr & ismember(get_month(stn.(qfld).date),[mos]));
    dtix = find(mindt<=stn.(qfld).date & stn.(qfld).date<=maxdt);
    q.date = stn.(qfld).date(dtix);
    q.data = stn.(qfld).data(dtix);

    if ( ~isempty(t.data(isfinite(t.data))) && ~isempty(q.data(isfinite(q.data))) )
      qs.date = q.date;
      qs.data = t.data(1) + cumsum(q.data) - q.data(1);

      subplot_tight(nrows,ncols,yrix);
      hold on;
      plot(t.date,t.data,'k',qs.date,qs.data,'b:');
      % xlim([datenum(yr,min(mos),1),datenum(yr,max(mos),31,23,59,59)]);
      xlim([mindt,maxdt-(1/1000)]);
      % datetick('x',2,'keeplimits');
      datetick('x',17,'keeplimits');
      set_datetick_cursor;
      % ylim([10,40]);
      ylim([18,33]);
      grid on;
    end;
  end;








  nrows = 3;
  ncols = 2;

  if ( numel(begyr:endyr) > 6 )
    nrows = 3;
    ncols = 3;

    if ( numel(begyr:endyr) > 9 )
      nrows = 4;
      ncols = 3;

      if ( numel(begyr:endyr) > 12 )
        nrows = 5;
        ncols = 3;

        if ( numel(begyr:endyr) > 15 )
          nrows = 5;
          ncols = 4;

          if ( numel(begyr:endyr) > 20 )
            begyr = endyr - 19;
          end;

        end;

      end;

    end;

  end;




  if ( ~exist('comparisonsfld','var') || isempty(comparisonsfld) )
    comparisonsfld = 'ndbc_sea_t';
  else
    comparisonsfld = sfld;
  end;









%%function annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,SPREFIX,commentstr,mos,begyr,endyr,comparisonsfld)




  qsfld = [SPREFIX ISPFX '_sea_spechumid'];

  switch (ISPFX),
   case 'ndbc',		WINDINFIX = '_wind1';
   otherwise,		WINDINFIX = '_wind';
  end;
  Wfld = [ISPFX '_wind1_speed'];
  Dfld = [ISPFX '_wind1_dir'];
  Ufld = [ISPFX '_wind1_u'];
  Vfld = [ISPFX '_wind1_v'];

  %%%% ??? DEBUG: Low-pass filter winds for quasi-Eulerian currents
  Ulpfld = [Ufld '_72_hour_lowpass'];
  Vlpfld = [Vfld '_72_hour_lowpass'];
  Wlpfld = [Wfld '_72_hour_lowpass'];
  Dlpfld = [Dfld '_72_hour_lowpass'];
  %%%% ??? DEBUG: Low-pass filter winds for quasi-Eulerian currents











  nrows = 4;
  ncols = 2;

  if ( numel(begyr:endyr) > 8 )
    ncols = 3;
    if ( numel(begyr:endyr) > 12 )
      nrows = 5;
      if ( numel(begyr:endyr) > 15 )
        ncols = 4;
        if ( numel(begyr:endyr) > 20 )
          begyr = endyr - 19;
        end;
      end;
    end;
  end;









      %%%% ??? DEBUG: Low-pass filter winds for quasi-Eulerian currents
      stn.commentstr = [stn.commentstr ' (QE~72hlp) '];
      stn = station_stokes_drift(stn,sssfld,ssdfld,ssufld,ssvfld,Wlpfld,Dlpfld,whfld,wpfld,wdfld);









  sfldfname = '';
  if ( exist('sfld','var') && ~isempty(sfld) )
    sfldfname = ['_' sfld];
  end;


  matfname = fullfile(datapath,...
                      [lower(stn.station_name) '_' TURPFX '_' QEPFX ...
                      sfldfname '_heat_budget.mat']);








  %%%
  %% Call SCRIPT to set:
  %% Variable-name prefixes ("PFX") for various input and output datasets; AND,
  %% All station struct fieldnames used to produce heat budget 
  station_heat_budget_field_names;


  % Special handling for non-SEAKEYS stations
  switch (lower(stn.station_name))
   case 'looe1',
    sfld = 'microcat_seatemp';
    error('Have not finished implementing LOOE1 yet!');
  end;

  %%%
  %% Heat Budget Options structure
  if ( isfield(stn,'station_heat_budget_options') )
    opts = stn.station_heat_budget_options;
  else
    opts = station_heat_budget_options(stn);
  end;








  endjd = datenum(1,mos(end)+1,1) - datenum(1,1,1) + (23.5/24);





stn.avhrr_weekly_sst.date = stn.avhrr_weekly_sst_field.date;
stn.avhrr_weekly_sst.data = ...
    interp_field(stn.avhrr_weekly_sst_field.lat,stn.avhrr_weekly_sst_field.lon,...
                 stn.avhrr_weekly_sst_field.field,stn.lat,stn.lon,'linear');






stn = rmfield(stn,'hourly_fkeys_seatemp_xshore');
stn = rmfield(stn,'hourly_fkeys_seatemp_lshore');




    %%%% ??? DEBUG
    % q0fld = q030fld;
    % qtfld = qt30fld;
    % stn.commentstr = [stn.commentstr ' (NO WARM LAYER) '];
    %%%% ??? DEBUG

    % Adjust latent heat flux by a fixed factor, as in Monismith et al (2006)




  axis([1,n,18,31]);






  % Peak Period below 2s is NOT real - just means a model was spinning up!
  if ( isfield(ww3,'ww3_peakwaveper') )
    badix = find(ww3.ww3_peakwaveper.data<2);
    if ( ~isempty(badix) )
    end;
  end;






    % newdat.(newfld) = interp_ts(newdat.(rawfld));
    newdat.(newfld).date = ...
        [ww3.(fld).date(1):(1/24):ww3.(fld).date(end)]';
    goodix = find(isfinite(ww3.(fld).data));
    newdat.(newfld).data = ...
        interp1(ww3.(fld).date(goodix),ww3.(fld).data(goodix),newdat.(newfld).date);








    if ( nargout > 1 )
      ww3 = dat.ww3;
    end;




    disp(['Saving to ' ww3matfname '...']);
    station = stn;
    save(ww3matfname,'station','ww3');
    station = []; clear station;

    if ( nargout < 2 )
      ww3 = [];
      clear ww3;
    end;












    % startdt = datenum(2005,1,1);
    % startdt = datenum(2005,1,23);



    opts.K_theta = get_opt(opts,'K_theta',model_K_theta);
    if ( isnumeric(opts.K_theta) && ~isscalar(opts.K_theta) )
      stn.commentstr = [stn.commentstr ' K_\theta seas. clim. (' num2str(opts.K_theta,'%g ') ') '];
    elseif ( all(isfield(opts.K_theta,{'func','arg'})) )
      stn.commentstr = [stn.commentstr ' K_\theta~' strrep(opts.K_theta.arg,'_','\_') ' '];
      opts.K_theta = opts.K_theta.func(stn.(opts.K_theta.arg));
    end;





      stn.(grclim).date = stn.(qtfld).date;
      stn.(grclim).data = build_clim_opt(opts.gradient_climatology,grclim,stn.(qtfld).date);







  if ( ischar(ltfld) && ischar(dtfld) && isfield(stn,dtfld) )
    % stn.(ltfld) = ts_op(stn.(dtfld),stn.(kl),'+');
    stn.(ltfld) = ts_op(stn.(dtfld),stn.(kllp),'+');
  elseif ( ~isempty(ltfld) || ~isempty(dtfld) )
    warning('Ecoforecasts:StationKDel2T:BadFieldNames',...
            'Invalid fieldname(s) LTFLD or DTFLD for heat budget terms!');
  end;




  for facix = 1:length(facs)
    fac = facs(facix);
    % dat(:,facix) = x.sfc_btm.data(qix) + (fac.*x.dif.data(aix));
    dat(:,facix) = x.sfc_btm_adv.data(qix) + (fac.*x.dif.data(aix));
  end;





  if (~exist('K_theta','var') || isempty(K_theta))
    % From HYCOM v2.2 documentation: TEMDF2 * DX
    %    For FKEYS HYCOM ~ 2.5 [m^2/s]
    %    For GoM HYCOM ~ 20 [m^2/s]
    K_theta = 20.0e-3 * 1e3;
  end;
  if (~exist('fldnm','var') || isempty(fldnm))
    fldnm='fkeys_hycom_seatemp_field';
  end;
  if (~exist('rawkl','var') || isempty(rawkl))
    rawkl='native_fkeys_hycom_diffused_heat';
  end;
  if (~exist('kl','var') || isempty(kl))
    kl='fkeys_hycom_diffused_heat';
  end;
  if (~exist('dtfld','var') || isempty(dtfld))
    dtfld='fkeys_hycom_dt';
  end;
  if (~exist('ltfld','var') || isempty(ltfld))
    ltfld='fkeys_hycom_lt';
  end;
  if (~exist('interpMethod','var') || isempty(interpMethod))
    interpMethod='linear';
  end;







  % Make sure we have a Laplacian
  if ( ~isfield(stn.(fldnm),'laplacian') )
    stn = calc_field_terms(stn,fldnm);
  end;

  k = build_clim_opt(K_theta,'K_theta',stn.(fldnm).date);

  l = interp_field(stn.(fldnm).lat,stn.(fldnm).lon,stn.(fldnm).laplacian,stn.lat,stn.lon,interpMethod);
  stn.(rawkl).date = stn.(fldnm).date;
  % stn.(rawkl).data = -k.*3600.*l;
  stn.(rawkl).data = k.*3600.*l;

%%%% DEBUG???
  % stn.(rawkl).data(:) = 0;

  stn.(kl) = interp_ts(stn.(rawkl));




  if ( ischar(ltfld) && ischar(dtfld) && isfield(stn,dtfld) )
    % [ix1,ix2] = intersect_dates(stn.(dtfld).date,stn.(kl).date);
    % stn.(ltfld).date = stn.(dtfld).date(ix1);
    % stn.(ltfld).data = stn.(dtfld).data(ix1) + stn.(kl).data(ix2);
%%%% ??? DEBUG
    [ix1,ix2] = intersect_dates(stn.(dtfld).date,stn.(kllp).date);
    stn.(ltfld).date = stn.(dtfld).date(ix1);
    stn.(ltfld).data = stn.(dtfld).data(ix1) + stn.(kllp).data(ix2);
  elseif ( ~isempty(ltfld) || ~isempty(dtfld) )
    warning('Ecoforecasts:StationKDel2T:BadFieldNames',...
            'Invalid fieldname(s) LTFLD or DTFLD for heat budget terms!');
  end;





  stn.(kl) = interp_ts(stn.(rawkl));
  % stn.(kl).date = [stn.(rawkl).date(1):(1/24):stn.(rawkl).date(end)]';
  % stn.(kl).data = spline(stn.(rawkl).date,stn.(rawkl).data,stn.(kl).date);
  % stn = filter_gaps(stn,rawkl,kl);





  oris = minori:5:maxori;
  udT = repmat(nan,[length(dts),length(oris)]);
  qtAdv = repmat(nan,[length(qdts),length(oris)]);
  dT = repmat(nan,[length(ddts),length(oris)]);
  bdT = repmat(nan,[length(bdts),length(oris)]);
  for ix=1:length(oris)
    ori = oris(ix);

    % % stn = station_calc_udotdelt(stn,qeufld,qevfld,Tfld,kmtfld,...
    % stn = station_calc_udotdelt(stn,ssufld,ssvfld,Tfld,kmtfld,...
    %                             ['raw_' udTfld],udTfld,...
    %                             qtfld,qtAdvfld,grdInterpMethod);
    % commentstr = [commentstr ' U_S_S+V_S_S '];






From OPTIM_ORI.m:
  if ( ~isfield(stn,hufld) )
    stn.(hufld) = interp_ts(stn.(ufld));
    stn.(hvfld) = interp_ts(stn.(vfld));
  end;


    %%%% ??? DEBUG: Assume no km-scale advection at all
    stn.(ufld).data(:) = 0;
    stn.(vfld).data(:) = 0;
    stn.(hufld).data(:) = 0;
    stn.(hvfld).data(:) = 0;
    %%%% ??? DEBUG: Assume no km-scale advection at all

  if ( ~isfield(stn,netufld) )
    stn.(netufld) = ts_op(stn.(tufld),stn.(hufld),'+');
    stn.(netvfld) = ts_op(stn.(tvfld),stn.(hvfld),'+');
    stn = calc_quasi_eulerian(stn,STOKESPFX,KMPFX,QEPFX);
  end;
  %%%% DEBUG???


%%%% MLRF1 Stokes only ERAI waves
  oris = 75:5:95;










    % % Simplest approach - constant Kd
    % opts.kd = get_opt(opts,'kd',0.3);
    % % Follows seasonal pattern in Kd, but with a mean of 0.3 (to match constant)
    % opts.kd = get_opt(opts,'kd',[0.15,0.45,274]);



  % opts.K_theta.func = @(ts)(struct('date',ts.date,'data',2.0+(0.5.*ts.data)));
  opts.K_theta.func = @(ts)(struct('date',ts.date,'data',20.*exp(ts.data./10)));





    elseif ( isa(opts.K_theta,'function_handle') )
      stn.commentstr = [stn.commentstr ' K_\theta~' strrep(Wlpfld,'_','\_') ' '];
      opts.K_theta = opts.K_theta(stn.(Wlpfld));
    end;










  if ( is_valid_ts(op) && isnumeric(op.data) )

    disp( [ 'Using time series ' opnm ] );
    [opix,dtix] = intersect_dates(op.date,dts);
    if ( numel(opix) ~= numel(dts) )
      error('Date mistmatch between %s time series and DTS',opnm);
    end;
    val = op.data(opix);

  elseif ( isnumeric(op) )








    else
      commentstr = [commentstr ' \partial_x_sT seas. clim. (' num2str(opts.gradient_climatology) ') '];

      %DEBUG:
      grclim = [hkmtxsfld '_clim'];

      % stn.(grclim).date = stn.(qtfld).date;
      % stn.(grclim).data = build_clim_opt(opts.gradient_climatology,grclim,stn.(qtfld).date);

      grdclim = opts.gradient_climatology;
      mingr = grdclim(1); maxgr = grdclim(2); peakjd = grdclim(3);

      meangr = mean([mingr,maxgr]);
      grrng = maxgr - meangr;
      peakrads = (2*pi*peakjd/366);

      %DEBUG:
      grclim = [hkmtxsfld '_clim'];
      disp( [ 'Using seasonal sine-modulated cross-shore gradient = ' num2str([mingr,maxgr,peakjd]) ] );
      stn.(grclim).date = stn.(qtfld).date;
      stn.(grclim).data = meangr + grrng.*cos((2.*pi.*get_yearday(stn.(qtfld).date)/366)-peakrads);

      if ( ~isfield(stn,qexsfld) )
        stn = station_reorient_vectors(stn,bathorifld,qeufld,qevfld,qexsfld,qelsfld);
      end;

      stn.(udTfld) = ts_op(stn.(qexsfld),stn.(grclim),'*');
      stn.(qtAdvfld) = ts_op(stn.(udTfld),stn.(qtfld),'+');

    end;








    else

      grdclim = opts.gradient_climatology;
      commentstr = [commentstr ' \partial_x_sT seas. clim. (' num2str(grdclim) ') '];
      mingr = grdclim(1); maxgr = grdclim(2); peakjd = grdclim(3);

      meangr = mean([mingr,maxgr]);
      grrng = maxgr - meangr;
      peakrads = (2*pi*peakjd/366);

      %DEBUG:
      grclim = [hkmtxsfld '_clim'];
      disp( [ 'Using seasonal sine-modulated cross-shore gradient = ' num2str([mingr,maxgr,peakjd]) ] );
      stn.(grclim).date = stn.(qtfld).date;
      stn.(grclim).data = meangr + grrng.*cos((2.*pi.*get_yearday(stn.(qtfld).date)/366)-peakrads);

      if ( ~isfield(stn,qexsfld) )
        stn = station_reorient_vectors(stn,bathorifld,qeufld,qevfld,qexsfld,qelsfld);
      end;

      stn.(udTfld) = ts_op(stn.(qexsfld),stn.(grclim),'*');
      stn.(qtAdvfld) = ts_op(stn.(udTfld),stn.(qtfld),'+');
    end;










      lplclim = [hkmtlfld '_clim'];





      lplclim = opts.laplacian_climatology;
      commentstr = [commentstr ' \nabla^2T seas. clim. (' num2str(lplclim) ') '];
      minlpl = lplclim(1); maxlpl = lplclim(2); peakjd = lplclim(3);

      meanlpl = mean([minlpl,maxlpl]);
      lplrng = maxlpl - meanlpl;
      peakrads = (2*pi*peakjd/366);

      %DEBUG:
      lplclim = [hkmtlfld '_clim'];
      disp( [ 'Using seasonal sine-modulated Laplacian = ' num2str([minlpl,maxlpl,peakjd]) ] );
      stn.(lplclim).date = stn.(qtfld).date;
      stn.(lplclim).data = meanlpl + lplrng.*cos((2.*pi.*get_yearday(stn.(qtfld).date)/366)-peakrads);










    else
      % lplclim = [hkmtlfld '_clim'];
      % stn.(lplclim).date = stn.(qtfld).date;
      % stn.(lplclim).data = build_clim_opt(opts.laplacian_climatology,lplclim,stn.(qtfld).date);

      lplclim = opts.laplacian_climatology;
      commentstr = [commentstr ' \nabla^2T seas. clim. (' num2str(lplclim) ') '];
      minlpl = lplclim(1); maxlpl = lplclim(2); peakjd = lplclim(3);

      meanlpl = mean([minlpl,maxlpl]);
      lplrng = maxlpl - meanlpl;
      peakrads = (2*pi*peakjd/366);

      %DEBUG:
      lplclim = [hkmtlfld '_clim'];
      disp( [ 'Using seasonal sine-modulated Laplacian = ' num2str([minlpl,maxlpl,peakjd]) ] );
      stn.(lplclim).date = stn.(qtfld).date;
      stn.(lplclim).data = meanlpl + lplrng.*cos((2.*pi.*get_yearday(stn.(qtfld).date)/366)-peakrads);

      k = build_clim_opt(K_theta,'K_theta',stn.(lplclim).date);

      stn.(kd2Tfld).date = stn.(lplclim).date;
      stn.(kd2Tfld).data = 3600 .* k .* stn.(lplclim).data;
      stn.(dTfld) = ts_op(stn.(kd2Tfld),stn.(qtAdvfld),'+');

    end;







      if ( isscalar(Kd) )
        %DEBUG:
        disp(['Using constant Kd = ' num2str(Kd)]);

      else
        % Among other papers, Shi and Wang (2010) present clear evidence of
        % seasonality in, e.g., coastal Kd(490): allow caller to use that!
        switch ( numel(Kd) )
         case 2,	minkd = Kd(1); maxkd = Kd(2); peakjd = 91.5;
         case 3,	minkd = Kd(1); maxkd = Kd(2); peakjd = Kd(3);
         otherwise,	error('Invalid OPTS.kd vector length!');
        end;

        meankd = mean([minkd,maxkd]);
        kdrng = maxkd - meankd;
        peakrads = (2*pi*peakjd/366);

        %DEBUG:
        disp( [ 'Using seasonal sine-modulated Kd = ' num2str([minkd,maxkd,peakjd]) ] );
        Kd = meankd + kdrng.*cos((2.*pi.*get_yearday(stn.(qswfld).date(qswix))/366)-peakrads);

        %DEBUG:        [cum,tid]=grp_ts(Kd,stn.(qswfld).date(qswix),[],[],1);
        %DEBUG:        figure; plot(tid,cum); titlename(num2str(peakjd));
      end; %if ( isscalar(Kd) ) else







    elseif ( isnumeric(Kd) )
      if ( numel(Kd) == 2 )
        % Default peak year-day for seasonal climatology is in MARCH
        Kd(3) = 91.5;
      end;
      Kd = build_clim_opt(Kd,'Kd',stn.(qswfld).date(qswix));

      if ( isscalar(Kd) )
        %DEBUG:
        disp(['Using constant Kd = ' num2str(Kd)]);

      else
        % Among other papers, Shi and Wang (2010) present clear evidence of
        % seasonality in, e.g., coastal Kd(490): allow caller to use that!
        switch ( numel(Kd) )
         case 2,	minkd = Kd(1); maxkd = Kd(2); peakjd = 91.5;
         case 3,	minkd = Kd(1); maxkd = Kd(2); peakjd = Kd(3);
         otherwise,	error('Invalid OPTS.kd vector length!');
        end;

        meankd = mean([minkd,maxkd]);
        kdrng = maxkd - meankd;
        peakrads = (2*pi*peakjd/366);

        %DEBUG:
        disp( [ 'Using seasonal sine-modulated Kd = ' num2str([minkd,maxkd,peakjd]) ] );
        Kd = meankd + kdrng.*cos((2.*pi.*get_yearday(stn.(qswfld).date(qswix))/366)-peakrads);

        %DEBUG:        [cum,tid]=grp_ts(Kd,stn.(qswfld).date(qswix),[],[],1);
        %DEBUG:        figure; plot(tid,cum); titlename(num2str(peakjd));
      end; %if ( isscalar(Kd) ) else

    end; %if ( ischar(Kd) ) elseif ( isnumeric(Kd) )









  %%%% DEBUG???
  % opts.gradient_climatology = [-2e-4, +2e-4, 355];





    %%%% ??? DEBUG: Assume no km-scale advection at all
    stn.(ufld).data(:) = 0;
    stn.(vfld).data(:) = 0;
    stn.(hufld).data(:) = 0;
    stn.(hvfld).data(:) = 0;
    commentstr = [commentstr ' (NO KM-SCALE ADVEC) '];
    %%%% ??? DEBUG: Assume no km-scale advection at all







    disp(['Loading from ' matfname]);
    load(matfname,'station');
    flds = fieldnames(station);
    for fldix=1:length(flds)
      fld = flds{fldix};
      if ( ~isfield(stn,fld) )
        stn.(fld) = station.(fld);
      else
        disp(['Field already in STN: ' fld]);
      end;
    end;
    station = []; clear station;







if ( ~exist('stn','var') || ~isfield(stn,'station_name') || ~strcmpi(stn.station_name,stnm) )
  stn = []; clear stn;
  stn = station_heat_budget(stnm,'erai','gom_hycom',[],[],'ww3');
end;





%colormap(flipud(flag(3))); 

colormap(flipud(logflag(5))); 



fmg; plot_ts(stn.gom_hycom_seatemp_xshore,stn.avhrr_weekly_sst_xshore);
legend('GOM','AVHRR 1km'); titlename('\nablah^.\nablaT: Model vs. Satellite');






          d = [ 0 ; sw_dist(lat,lon,'km') ];
          % First value will be +Inf
          dsst = [ 0 ; diff(sst) ] ./ d;







  end; %if ( exist(matfname,'file') ) else

  [ig,hgs] = sw_dist(res.lat,res.lon,'km');
  % Convert squirrelly SW_DIST heading -> degrees True
  % SW_DIST outputs:
  %   phaseangle  = angle of line between stations with x axis (East).
  %                 Range of values are -180..+180. (E=0, N=90, S=-90)
  hgs = (-hgs + 90);
  hgs(hgs < 0) = hgs(hgs < 0) + 360;
  res.hdg = [ nan ; hgs ];

  set_more;

return;







  if ( ~isfield(res,'date') )
    for yrs = unique(res.yrs)'
    end;
  end;





  stn.smooth_tsg.lon  = mean([stn.tsg.lon(1:end-2),stn.tsg.lon(2:end-1),stn.tsg.lon(3:end)],2);
  stn.smooth_tsg.lat  = mean([stn.tsg.lat(1:end-2),stn.tsg.lat(2:end-1),stn.tsg.lat(3:end)],2);
  stn.smooth_tsg.d    = sum([stn.tsg.d(1:end-2),stn.tsg.d(2:end-1),stn.tsg.d(3:end)],2);
  stn.smooth_tsg.dsst = mean([stn.tsg.dsst(1:end-2),stn.tsg.dsst(2:end-1),stn.tsg.dsst(3:end)],2);







    % Ignore the along-shore component of model heat advection
    % % And fix isobath orientation 35o further to the right
    % stn.(bathorifld) = stn.(bathorifld) + 35;
    stn = station_cross_shore_advection(stn,bathorifld,...
                                        qeufld,qevfld,Tfld,kmtfld,...
                                        ['raw_' udTfld],udTfld,...
                                        qtfld,qtAdvfld,opts.grid_interp_method);
    commentstr = [commentstr ' cross-shore only (' num2str(stn.(bathorifld)) 'o) '];







      % fld = fld(:,:,end:-1:1);
      fld = flipdim(fld,3);
      fliplon = true;






  stn.tsg.lon = res.lon(goodix);
  stn.tsg.lat = res.lat(goodix);
  stn.tsg.dsst = res.crestdsst(goodix);






  if ( ~isfield(stn,'ngdc_92m_bathy') )
    stn = get_ngdc_bathy_station(stn);
  end;

  z = interp_field(stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.lon,...
                   stn.ngdc_92m_bathy.field,stn.lat,stn.lon);





disp(['INTERPMETHOD is ' interpMethod]);




    % Ignore the along-shore component of model heat advection
    % stn = station_cross_shore_advection(stn,90,...
    stn.(bathorifld) = 90;
    stn = station_cross_shore_advection(stn,bathorifld,...
                                        qeufld,qevfld,Tfld,kmtfld,...
                                        ['raw_' udTfld],udTfld,...
                                        qtfld,qtAdvfld,opts.grid_interp_method);
    commentstr = [commentstr ' cross-shore only (' num2str(stn.(bathorifld)) 'o) '];






%%%% DEBUG???
  end; %if ( exist(matfname,'file') ) else

%%%% DEBUG???
  % end; %if ( exist(matfname,'file') ) else






    % Use CONTOUR to find "reef crest" (z btw -20 and -50) coordinates
    fh=fmg; [ig,ig,cs,ch] = map_freef([-81.9,-80.0,24.2,26.0],[-50 -20]); close(fh);
    pause(2);
    coordix=find(~ismember(cs(1,:),[-50 -20]));

    res.crestlon=cs(1,coordix);
    res.crestlat=cs(2,coordix);



cresthi=-11;
crestlo=-60;
fh=fmg; [ig,ig,cshi,ch] = map_freef([-81.9,-80.0,24.2,26.0],cresthi); close(fh);
fh=fmg; [ig,ig,cslo,ch] = map_freef([-81.9,-80.0,24.2,26.0],crestlo); close(fh);
bkix=find(cshi(1,:)==cresthi); cs = [cslo(:,end:-1:2) , cshi(:,bkix(1)+1:1:bkix(2)-1)];
crestlon=cs(1,:); crestlat=cs(2,:);
save('frt_reef_crest.mat','crestlon','crestlat','cresthi','crestlo');




  % Trim extraneous trailing elements - if any
  flds = fieldnames(res);
  for fldix=1:length(flds)
    fld=flds{fldix};
    lastix = find(isfinite(res.(fld)),1,'last');
    res.(fld) = res.(fld)(1:lastix);
  end;







          badix = find(~isfinite(lat) | ~isfinite(lon) | ~isfinite(sst));
          lon(badix) = [];
          lat(badix) = [];
          sst(badix) = [];
          sss(badix) = [];
          chl(badix) = [];








        res.lon(ngot+1:ngot+npts)=x.data(:,1);
        res.lat(ngot+1:ngot+npts)=x.data(:,2);
        res.yr(ngot+1:ngot+npts)=x.data(:,3);
        res.jd(ngot+1:ngot+npts)=x.data(:,4);
        res.sst(ngot+1:ngot+npts)=x.data(:,5);
        res.sss(ngot+1:ngot+npts)=x.data(:,6);
        if ( size(x.data,2) == 6 )
          warning('No chl data in "%s"',fname);
        else
          res.chl(ngot+1:ngot+npts)=x.data(:,7);
          if ( size(x.data,2) > 7 )
            warning('Extra column(s) in "%s"',fname);
          end;
        end;






        res.lon(ngot+1:ngot+npts)=x.data(:,1);
        res.lat(ngot+1:ngot+npts)=x.data(:,2);
        res.yr(ngot+1:ngot+npts)=x.data(:,3);
        res.jd(ngot+1:ngot+npts)=x.data(:,4);
        res.sst(ngot+1:ngot+npts)=x.data(:,5);
        res.sss(ngot+1:ngot+npts)=x.data(:,6);
        if ( size(x.data,2) == 6 )
          warning('No chl data in "%s"',fname);
        else
          res.chl(ngot+1:ngot+npts)=x.data(:,7);
          if ( size(x.data,2) > 7 )
            warning('Extra column(s) in "%s"',fname);
          end;
        end;








  if ( ~isfield(stn,gxfld) )
    stn.(gxfld).date = stn.(fldnm).date;
    stn.(gxfld).data = ...
        interp_field(stn.(fldnm).lat,stn.(fldnm).lon,stn.(fldnm).gradient_x,stn.lat,stn.lon,interpMethod);
  end;
  if ( ~isfield(stn,gyfld) )
    stn.(gyfld).date = stn.(fldnm).date;
    stn.(gyfld).data = ...
        interp_field(stn.(fldnm).lat,stn.(fldnm).lon,stn.(fldnm).gradient_y,stn.lat,stn.lon,interpMethod);
  end;







  % Calculate along- and cross-shore field gradients
  if ( ~isfield(stn,xgfld) || ~isfield(stn,lgfld) )
    [stn,xgfld,lgfld] = station_reorient_vectors(stn,ori,gxfld,gyfld);
    xg = stn.(xgfld).data;
    lg = stn.(lgfld).data;
  end;






      % stn.(hfld).date = actualdts;




dsst = [ 0 ; diff(sst)./d ];





  if ( exist(matfname,'file') )

    disp(['Loading from ' matfname]);
    load(matfname,'station');
    flds = fieldnames(station);
    for fldix=1:length(flds)
      fld = flds{fldix};
      stn.(fld) = station.(fld);
    end;
    station = []; clear station;

%%%% DEBUG??? Fixups: stuff I may have added after the last recalc-from-scratch
    if ( ~isfield(stn,bathorifld) )
      stn = station_optimal_isobath_orientation(stn);
    end;
    if ( ~isfield(stn.(Tfld),'gradient_x') )
      % Calculate gradients and field Laplacians
      stn = calc_field_terms(stn,Tfld);
    end;
    if ( ~isfield(stn,tspdfld) || ~isfield(stn,tdirfld) )
      switch (TIDEPFX),
       case 'tmd_tide',		stn = station_tmd_tide(stn);
       case 'tpxo_tide',	stn = station_tmd_tide(stn);
      end;
    end;
    if ( ~isfield(stn,mhfld) )
      % Calculate the mean tide depth experienced by a watermass moving over
      % an M2 tidal ellipse centered on the coordinates of our station
      stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
    end;
%%%% DEBUG???











     % case 2010,	mos = 1:10;




stnms={'cmrc3','srvi2','lppr1','dbjm1','lciy2','41140','lkwf1','fwyf1','cryf1','mlrf1','lonf1','tnrf1','mose1','smkf1','looe1','amsf1','sanf1','plsf1','dryf1','42003',};



42003_erai.mat
cmrc3_erai.mat
fwyf1_erai.mat
lciy2_erai.mat
lkwf1_erai.mat
lonf1_erai.mat
looe1_erai.mat
lppr1_erai.mat
mlrf1_erai.mat
sanf1_erai.mat
smkf1_erai.mat
srvi2_erai.mat







     %DEBUG:     case 2010,		mos = 3;





  % Calculate cross-shore velocity component
  xfld = regexprep(unm,'_u$','_xshore');
  lfld = regexprep(unm,'_v$','_lshore');
  [uix,vix,fldix] = intersect_all_dates([],stn.(unm).date,stn.(vnm).date,...
                                        stn.(fldnm).date);
  dts = stn.(unm).date(uix);
  u = stn.(unm).data(uix);
  v = stn.(vnm).data(vix);
  % Calculate cross-shore velocity component
  [x,l] = reorient_vectors(ori,u,v);



  % Calc. advected field at native time resolution ([m/s]*[K/m]->[K/s] for heat)

  if ( ~isfield(stn.(fldnm),'gradient_x') || ~isfield(stn.(fldnm),'gradient_y') )
    % Calculate gradients and field Laplacians
    stn = calc_field_terms(stn,fldnm);
  end;
  gx = ...
      interp_field(stn.(fldnm).lat,stn.(fldnm).lon,stn.(fldnm).gradient_x,stn.lat,stn.lon,interpMethod);
  gy = ...
      interp_field(stn.(fldnm).lat,stn.(fldnm).lon,stn.(fldnm).gradient_y,stn.lat,stn.lon,interpMethod);
  % Calculate cross-shore field gradient
  [xg,lg] = reorient_vectors(ori,gx(fldix),gy(fldix));




    %%%
    %% Radiative Fluxes

    opts = [];

    % % Follows seasonal pattern from various regional coastal studies
    % opts.kd = [0.15,0.45];


    %%%
    %% Benthic Heat Exchanges

    opts = [];



    opts.scaling = 'SS';
    % opts.scaling = 'US';
    % opts.scaling = 'SU';
    % opts.scaling = 'UU';





%%%% DEBUG???
    sai_opts.kd = [0.1,0.3,274];


%%%% DEBUG???
    sai_opts.kd = [0.1,0.25,91];






    % ALSO fix any units we got wrong, e.g., rain should be mm/hr, not m/hr!
    result(stix).erai_conv_precip.data = result(stix).erai_conv_precip.data*1e3;
    result(stix).erai_precip.data = result(stix).erai_precip.data*1e3;



%%%% ??? DEBUG
% K_theta = 100;
% K_theta = 10;
      commentstr = [commentstr ' K_\theta=' num2str(K_theta) ' '];
%%%% ??? DEBUG



  x.non_qlh = x.q0 - x.qlh;
  x.non_qlh_qb = x.bq0 - x.qlh;

  x.climnon_qlh = x.climq0 - x.climqlh;

  % plot(1:n,[x.non_qlh,x.qlh,x.climnon_qlh,x.climqlh,]);
  % legend(ax(2),'Q_0-Q_L_H','Q_L_H','OAFlux: Q_0-Q_L_H','OAFlux: Q_L_H');
  % plot(1:n,[x.non_qlh_qb,x.qlh]);
  % legend(ax(2),'Q_0-Q_L_H','Q_L_H');
  % plot(1:n,[x.asr,x.lr,x.qlh,x.qsh,x.coolf]);







  ax(2)=subplot_tight(3,1,[2,3]); hold on; grid on;
  plot(1:n,x.radif,'r',...
       1:n,x.aradif,'m',...
       1:n,x.turif,'b',...
       1:n,x.climradif,'r',...
       1:n,x.climturif,'b');
  plh=plot(1:30:n,x.radif(1:30:end),'rs',...
           1:30:n,x.aradif(1:30:end),'ms',...
           1:30:n,x.turif(1:30:end),'bs',...
           1:30:n,x.climradif(1:30:end),'ro',...
           1:30:n,x.climturif(1:30:end),'bo');
  legend(plh,'ERAI Q_S_W + ERAI/insitu Q_L_W','\gammaQ_S_W+Q_L_W','Q_L_H+Q_S_H','ISCCP: Q_S_W+Q_L_W','OAFlux: Q_L_H+Q_S_H');





function x = chkann(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%function x = chkann(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)

  station_heat_budget_field_names;

  if ( ~exist('per','var') || isempty(per) )
    per = 'daily';
  end;

  sumfun = @nanmedian;

  switch ( per ),
   case 'hourly',   cumfun = @get_jhour_no_leap; n = 365*24; N = 1;
   case 'daily',    cumfun = @get_jday_no_leap;  n = 365;    N = 24;
   case 'pentad',   cumfun = @get_pentad;        n = 73;     N = 24*5;
   case 'weekly',   cumfun = @get_week;          n = 52;     N = 24*7;
   case 'monthly',  cumfun = @get_month;         n = 12;     N = 24*[31,28,31,30,31,30,31,31,30,31,30,31]';
   case 'seasonal', cumfun = @get_season;        n = 4;      N = 24*[31+28+31,30+31+30,31+31+30,31+30+31]';
   otherwise,      error('Do not know how to do a "%s" climatology!',char(per));
  end;

  x.t = grpstats(stn.(sfld).data,cumfun(stn.(sfld).date),sumfun);

  x.q0 = grpstats(stn.(q0fld).data,cumfun(stn.(q0fld).date),sumfun);
  x.qt = N.*cumsum(grpstats(stn.(qtfld).data,cumfun(stn.(qtfld).date),sumfun));
  x.sq0 = grpstats(stn.(sq0fld).data,cumfun(stn.(sq0fld).date),sumfun);
  x.sqt = N*cumsum(grpstats(stn.(sqtfld).data,cumfun(stn.(sqtfld).date),sumfun));
  x.bq0 = grpstats(stn.(bq0fld).data,cumfun(stn.(bq0fld).date),sumfun);
  x.bq0t = N.*cumsum(grpstats(stn.(bq0tfld).data,cumfun(stn.(bq0tfld).date),sumfun));

  x.bdT = N.*cumsum(grpstats(stn.(bdTfld).data,cumfun(stn.(bdTfld).date),sumfun));
  x.hc_dTdt = N.*cumsum(grpstats(stn.(hcdTdt).data,cumfun(stn.(hcdTdt).date),sumfun));

  x.sr = grpstats(stn.(srfld).data,cumfun(stn.(srfld).date),sumfun);
  x.asr = grpstats(stn.(asrfld).data,cumfun(stn.(asrfld).date),sumfun);
  x.lr = grpstats(stn.(lrfld).data,cumfun(stn.(lrfld).date),sumfun);
  x.qlh = grpstats(stn.(qlhfld).data,cumfun(stn.(qlhfld).date),sumfun);
  x.qsh = grpstats(stn.(qshfld).data,cumfun(stn.(qshfld).date),sumfun);
  x.qrh = grpstats(stn.(qrhfld).data,cumfun(stn.(qrhfld).date),sumfun);

  x.radif = x.qlh + x.qsh;
  x.coolif = x.lr + x.radif;

  x.non_qlh = x.q0 - x.qlh;
  x.non_qlh_qb = x.bq0 - x.qlh;

  x.climsr = grpstats(stn.(climsrfld).data,cumfun(stn.(climsrfld).date),sumfun);
  x.climlr = grpstats(stn.(climlrfld).data,cumfun(stn.(climlrfld).date),sumfun);
  x.climqlh = grpstats(stn.(climqlhfld).data,cumfun(stn.(climqlhfld).date),sumfun);
  x.climqsh = grpstats(stn.(climqshfld).data,cumfun(stn.(climqshfld).date),sumfun);
  x.climq0 = x.climsr + x.climlr + x.climqlh + x.climqsh;
  x.climnon_qlh = x.climq0 - x.climqlh;


  fh = figure,
  maxigraph;

  ax(1)=subplot_tight(3,1,[1]); hold on; grid on;
  % plot(1:n,[x.t,x.t(1)+x.sqt-x.sqt(1),x.t(1)+x.bq0t-x.bq0t(1),x.t(1)+x.hc_dTdt-x.hc_dTdt(1)]);
  % legend(ax(1),'T ','T(0)+\SigmaQ_0/\rhoC_ph ','T(0)+\Sigma(Q_0(\gamma)+Q_b)/\rhoC_ph ','T(0)+\Sigma\partial_tT_H_C ',...
  %        'Location','SouthEast', 'Orientation','horizontal');
  plot(1:n,[x.t,x.t(1)+x.bq0t-x.bq0t(1),x.t(1)+x.hc_dTdt-x.hc_dTdt(1)]);
  legend(ax(1),'T','T(0)+\Sigma(Q_0(\gamma)+Q_b)/\rhoC_ph','T(0)+\Sigma\partial_tT_H_C',...
         'Location','SouthEast', 'Orientation','horizontal');
  xlim([1,n]); ylim([10,40]);
  titlename([upper(stn.station_name) ': ' strrep(bdTfld,'_','\_') ' ' upper(per) ' ' upper(char(sumfun)) ' climatology']);

  ax(2)=subplot_tight(3,1,[2,3]); hold on; grid on;
  plot(1:n,[x.non_qlh,x.qlh,x.climnon_qlh,x.climqlh,]);
  legend(ax(2),'Q_0-Q_L_H','Q_L_H','OAFlux: Q_0-Q_L_H','OAFlux: Q_L_H');
  % plot(1:n,[x.non_qlh_qb,x.qlh]);
  % legend(ax(2),'Q_0-Q_L_H','Q_L_H');
  % plot(1:n,[x.asr,x.lr,x.qlh,x.qsh,x.coolf]);
  xlim([1,n]); ylim([-250,250]);

  if ( nargout < 1 )
    x = []; clear x;
  end;

return;









% MLRF1 GOM new gradient
oris=47:56;




    %%%% DEBUG???
    curInterpMethod = 'linear';
    grdInterpMethod = 'linear';
    commentstr = [commentstr ' ' strrep(KMPFX,'_','\_') ' ' curInterpMethod '/' grdInterpMethod ' '];
    stn = get_fkeys_hycom(stn,[],[],[],[],curInterpMethod);
    stn.(hufld) = interp_ts(stn.(ufld));
    stn.(hvfld) = interp_ts(stn.(vfld));
    stn.(netufld) = ts_op(stn.(tufld),stn.(hufld),'+');
    stn.(netvfld) = ts_op(stn.(tvfld),stn.(hvfld),'+');
    stn = calc_quasi_eulerian(stn,STOKESPFX,KMPFX,QEPFX);
    %%%% DEBUG???












  % Whether we loaded data, extracted data from netCDF files, or both, limit
  % our return to only dates the user requested - that have reasonable data
  for vix = 1:length(vars)
    fld = ['gom_hycom_' flds{vix}];

    for ix = 1:numel(stations)
      % Gross quality control
      if ( ~isfield(stations(ix),fld) || isempty(stations(ix).(fld)) )
        warning('Ecoforecasts:GomHycom:MissingField',...
                'No field "%s" found after load!',fld);
      elseif ( ~isempty(strfind(fld,'_field')) )
        % Make sure lon and lat are row vectors as for other models
        stations(ix).(fld).lon = stations(ix).(fld).lon(:);
        stations(ix).(fld).lat = stations(ix).(fld).lat(:);

        % With qtot in the mix, we are not enforcing tight bounds on any var
        ixes = find(-3000 >= stations(ix).(fld).field | stations(ix).(fld).field >= 3000);
        stations(ix).(fld).field(ixes) = nan;

        ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).field = stations(ix).(fld).field(ixes,:,:);
      else
        % With qtot in the mix, we are not enforcing tight bounds on any var
        ixes = find(-3000 <= stations(ix).(fld).data & stations(ix).(fld).data <= 3000);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).data = stations(ix).(fld).data(ixes);

        ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).data = stations(ix).(fld).data(ixes);
      end;
    end; %for ix = 1:numel(stations)
  end; %for vix = 1:length(vars)

  % If caller wants different interpolation method than we saved in MAT file,
  % re-call INTERP_FIELD on all the appropriate "_field" fields (if present)
  for ix = 1:numel(stations)
    if ( ~isfield(stations(ix),'gom_hycom_interp_method') || ...
         ~strcmpi(stations(ix).gom_hycom_interp_method,interpMethod) )
      disp(['Reinterpolating time series for ' stations(ix).station_name ...
            ' using ' upper(interpMethod)]);
      for vix = 1:length(vars)
        fld = ['gom_hycom_' flds{vix}];
        if ( isempty(strfind(fld,'_field')) )
          fnm = [fld '_field'];
          if ( isfield(stations(ix),fld) )
            if ( ~isfield(stations(ix),fnm) )
              warning('Ecoforecasts:GomHycom:NoFieldToReinterp',...
                      '%s: No field available to reinterpolate "%s"',...
                      stations(ix).station_name,fld);
            else
              f = stations(ix).(fnm);
              stations(ix).(fld).data = ...
                  interp_field(f.lat,f.lon,f.field,stations(ix).lat,stations(ix).lon,interpMethod,'warn');
            end;
          end;
        end;
      end;
    end;
  end;








  for ix = 1:numel(stns)
    matfname = fullfile(datapath,[lower(stns(ix).station_name) '_gom_hycom.mat']);
    % If we already did this before - do not need to load again.
    if ( exist(matfname,'file') )
      disp(['Reloading MAT file ' matfname]);
      load(matfname,'station');
      stations(ix) = station;
      station = []; clear station;
      needix(needix == ix) = [];
    else
      %%%% DEBUG ???
      %%%% HACK??? FIX GAPS DUE TO SERVER TIME-OUTS!
      oldmatfname = fullfile(datapath,[lower(stns(ix).station_name) '_gom_hycom_MINUS_2004_2009.mat']);
      if ( ~exist(oldmatfname,'file') )
        error('Missing GAP MAT file %s!',oldmatfname);
      else
        load(oldmatfname,'station');
        stations(ix) = station;
        station=[]; clear station;
      end;
      %%%% HACK??? FIX GAPS DUE TO SERVER TIME-OUTS!

      %%%% HACK??? FIX GAPS DUE TO SERVER TIME-OUTS!
      % if ( isfield(stns(ix),'station_name') )
      %   stations(ix).station_name = stns(ix).station_name;
      % end;
      % stations(ix).lon = stns(ix).lon; stations(ix).lat = stns(ix).lat;
      % stations(ix).gom_hycom_interp_method = lower(interpMethod);

    end;
  end;



    %%%% HACK??? FIX GAPS DUE TO SERVER TIME-OUTS!    allyrs = [2004 2009];



%%%% ??? DEBUG
      catch
keyboard;
        if ( exist('nc','var') && ~isempty(nc) )
          close(nc);
        end;
        clear nc;
        rethrow(lasterror);
        % lerr = lasterror;
        % msg = 'CAUGHT ERROR!';
        % if ( isfield(lerr,'identifier') )
        %   msg = [ msg ' ' lerr.identifier ];
        % end;
        % if ( isfield(lerr,'message') )
        %   msg = [ msg ' : ' lerr.message ];
        % end;
        % warning(msg);
      end;
%%%% ??? DEBUG









FROM get_pathfinder_sst.m:
  % baseurl = 'http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily/night/04km';
  baseurl = 'http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/avhrr/pathfinder/data_v5/daily';
  nighturl = [baseurl '/night/04km'];
  dayurl = [baseurl '/day/04km'];







function stn = get_gom_hycom(stn_or_stnm,mindt,maxdt,vars,flds,interpMethod,baseurl,doDownload)
%function stn = get_gom_hycom(stn_or_stnm,mindt,maxdt,vars,flds,interpMethod,baseurl,doDownload)

% Tries to download all netCDF files for local processing first, *unless*
% optional arg DODOWNLOAD is given and evaluates to False.






    %%%% ??? DEBUG:
    allyrs = [endyr20 (begyr30:endyr)+0.3];







    %%%
    %% Horizontal Convection

    bet = stn.(slopefld);
    % stn = station_horizontal_convection(stn,sfld,[],mhfld,tspdfld,bq0fld,bet,hc_opts);

    [tix,hix,tspdix,aix,Wix,q0ix,bq0ix,dTix] = ...
      intersect_all_dates([],stn.(sfld).date,stn.(mhfld).date,stn.(tspdfld).date,stn.(afld).date,stn.(Wfld).date,stn.(q0fld).date,stn.(bq0fld).date,stn.(bdTffld).date);
    t = stn.(sfld).data(tix);
    s = repmat(36,size(stn.(sfld).data(tix)));
    h = stn.(mhfld).data(hix);
    tspd = stn.(tspdfld).data(tspdix);
    at = stn.(afld).data(aix);
    W = stn.(Wfld).data(Wix);
    q0 = stn.(q0fld).data(q0ix);
    bq0 = stn.(bq0fld).data(bq0ix);
    dT = stn.(bdTffld).data(dTix);

    qstr=bdTffld;
    %%%% ??? DEBUG: Base horizontal convection on surface heating only!
    q=bq0; commentstr = [commentstr ' HC(Q0+Qb) '];
    %%%% ??? DEBUG: Base horizontal convection on total (km-scale) budget
    % q=dT;

    hc_opts.R = (1.00-0.08);

    hc_opts.scaling = 'SS';
    % hc_opts.scaling = 'US';
    % hc_opts.scaling = 'SU';
    % hc_opts.scaling = 'UU';
    hc_opts.maximum_onset_time = 6*3600;
    %%%% ??? DEBUG
    commentstr = [commentstr ' (HC ' hc_opts.scaling ' maxT:9h) '];
    res = horizontal_convection(t,s,h,q,bet,hc_opts);

    dts = stn.(sfld).date(tix);
    flds = fieldnames(res);
    for fldix = 1:length(flds)
      fld = flds{fldix};
      dat = res.(fld);
      res.(fld) = [];
      res.(fld).date = dts;
      res.(fld).data = dat;

      stnfld = ['hc_' fld];
      stn.(stnfld).date = dts;
      stn.(stnfld).data = dat;
    end;

    %%%% ??? DEBUG
    % Truncate start of all time series for display purposes
    startdt = datenum(2004,1,1);
    % startdt = datenum(2005,1,1);
    % startdt = datenum(2005,1,23);

    startix = find(dts>=startdt,1);

    %%%% DEBUG???
    ix = find(~isfinite(q0(startix:end)),1,'last');
    if ( ~isempty(ix) )
      startix = startix + ix;
    end;

    dts = dts(startix:end);
    t = t(startix:end);
    h = h(startix:end);
    tspd = tspd(startix:end);
    at = at(startix:end);
    W = W(startix:end);
    q0 = q0(startix:end);
    bq0 = bq0(startix:end);
    dT = dT(startix:end);
    fac = stn.hc_termFactor.data(startix:end);
    dTdt.date = stn.hc_dTdt.date(startix:end);
    dTdt.data = stn.hc_dTdt.data(startix:end);

    bt = stn.(btfld);
    [ig,startix] = min(abs(bt.date-dts(1)));
    [ig,endix] = min(abs(bt.date-dts(end)));
    bt.date = bt.date(startix:endix);
    bt.data = bt.data(startix:endix);
    %%%% ??? DEBUG

    dsr_lpfld = [asrfld '_24_hour_sum'];
    stn = verify_variable(stn,dsr_lpfld);
    [ig,dsr_lpix] = intersect_dates(dts,stn.(dsr_lpfld).date);

    % lhf_lpfld = [qlhfld '_24_hour_sum'];
    % stn = verify_variable(stn,lhf_lpfld);
    % [ig,lhf_lpix] = intersect_dates(dts,stn.(lhf_lpfld).date);

    udT_lpfld = [udTfld];
    stn = verify_variable(stn,udT_lpfld);
    [ig,udT_lpix] = intersect_dates(dts,stn.(udT_lpfld).date);

    %%%% ??? DEBUG
    T0 = t(1);
    bigfh=figure, maxigraph; hold on;
    plot(dts,t,dts,at,'k:',bt.date,bt.data,'m--',dts,T0+cumsum(q0.*fac),dts,T0+cumsum(bq0.*fac),dts,T0+cumsum(dT.*fac),dTdt.date,T0+cumsum(dTdt.data));
    % plot(stn.(dsr_lpfld).date(dsr_lpix),(stn.(dsr_lpfld).data(dsr_lpix)./1000),stn.(lhf_lpfld).date(lhf_lpix),-(stn.(lhf_lpfld).data(lhf_lpix)./1000),'y:','Color',[.8,.8,0]);
    % plot(stn.(dsr_lpfld).date(dsr_lpix),(stn.(dsr_lpfld).data(dsr_lpix)./1000),stn.(udT_lpfld).date(udT_lpix),T0+cumsum(stn.(udT_lpfld).data(udT_lpix)),'y:','Color',[.8,.8,0]);
    plot(stn.(dsr_lpfld).date(dsr_lpix),T0+(stn.(dsr_lpfld).data(dsr_lpix).*fac),stn.(udT_lpfld).date(udT_lpix),T0+cumsum(stn.(udT_lpfld).data(udT_lpix)),'y:','Color',[.8,.8,0]);
    % plot(dts,h,dts,tspd,'k:',dts,W,'o','Color',[0,.8,.2]);
    plot(dts,h,dts,tspd,'k:','Color',[0,.8,.2]);
    datetick3('x',2,'keeplimits');
    % legend('T_s','T_a','T_b','T_0+Q_0/\rhoC_ph','T_0+(Q_0+Q_b)/\rhoC_ph','T_0+\partial_tT_k_m',['T_0+\partial_tT ' hc_opts.scaling],'T_0+\Sigma_1_d Q_S_W\times1h/\rho^.C_p^.h','T_0+\Sigmau^.\nablaT_k_m','h_t_i_d_e','SPD_t_i_d_e','W', 'Location','SouthWest'); %'Best');
    legend('T_s','T_a','T_b','T_0+Q_0/\rhoC_ph','T_0+(Q_0+Q_b)/\rhoC_ph','T_0+\partial_tT_k_m',['T_0+\partial_tT ' hc_opts.scaling],'T_0+\Sigma_1_d Q_S_W\times1h/\rho^.C_p^.h','T_0+\Sigmau^.\nablaT_k_m','h_t_i_d_e','SPD_t_i_d_e', 'Location','SouthWest'); %'Best');
    titlename([ commentstr stn.station_name ' ' strrep(qstr,'_','\_') ]);
    disp(commentstr);
    grid on;
    % axis('tight');
    % axis('tight'); ylim([-50 150]);
    % axis([startdt,startdt+365,0,60]);
    % axis([startdt,startdt+93,0,26]);
    % axis([startdt,startdt+5,0,26]);
    axis([startdt,datenum(2008,12,31),-100,100]);
    %%%% ??? DEBUG
    datetick3('x',2,'keeplimits');


    % Plot daily climatology comparisons
    chkann(stn,'daily',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX);

    % [stn,hls,hxs]=multiplot_station(stn,{'absorbed_erai_srf','ndbc_erai_30a_latent_heat_flux','ndbc_erai_30a_net_heat_flux','benthic_erai_qbo','benthic_ndbc_erai_30a_net_heat_flux','ndbc_air_t','ndbc_sea_t'},[],[],{'SW','LH','Q_0','BH','Q_0+Q_b','T_a','T_s'},[datenum(2005,4,1),datenum(2005,5,1)],{[-1000,1000],[-1000,1000],[-1000,1000],[-1000,1000],[-1000,1000],[16,28],[16,28]});
    % axes(hxs(1));
    % hold on; plot(stn.erai_srf.date,stn.erai_srf.data,'ro','MarkerSize',1.5);












   % case 'mlrf1', minori= 40; maxori= 75;
   case 'mlrf1', minori= 75; maxori=110;


  % stn.(qtAdvfld) = stn.(qtfld);
    % stn.(qtAdvfld) = stn.(qtfld);


  stn.([bdTfld '_minus_adv']) = ts_op(stn.(bdTfld),stn.(qtAdvfld),'-');
  mdts = stn.([bdTfld '_minus_adv']).date;
  bdTma = repmat(nan,[length(mdts),length(oris)]);
    stn.([bdTfld '_minus_adv']) = ts_op(stn.(bdTfld),stn.(qtAdvfld),'-');
    bdTma(:,ix) = stn.([bdTfld '_minus_adv']).data;
  % res = bdTma; resnm='\Sigma [ K_\theta_H\nabla^2T_k_m + (Q_0+Q_b)/\rhoC_ph ]';
  for ix=1:length(oris)
    goodix = find( isfinite(bdTma(:,ix)) );
    bdTma(goodix,ix) = t(1) + cumsum(bdTma(goodix,ix));
  end;
  plot(mdts,bdTma,'r');
  % lh = plot(mdts,res);








stn.(Tfld) = rmfield(stn.(Tfld),{'gradient_x','gradient_y','laplacian'});






    % Inter-gridpoint distances in [m] for gradient calculation
    midlat = repmat(tstr.lat(midy),size(tstr.lon));
    dy = [0 ; cumsum(sw_dist(midlat,tstr.lon,'km'))] .* 1e3;
    midlon = repmat(tstr.lon(midy),size(tstr.lat));
    dx = [0 ; cumsum(sw_dist(tstr.lat,midlon,'km'))] .* 1e3;
    dt = tstr.date(fldix);

    %DEBUG:
    disp('Calculating GRADIENT');
    % [dTdx,dTdy,dTdt] = gradient(permute(fld,[2 3 1]),dx,dy,dt);
    [dTdx,dTdy,dTdt] = gradient(permute(fld,[3 2 1]),dx,dy,dt);
    dTdx = permute(dTdx,[3 2 1]);
    dTdy = permute(dTdy,[3 2 1]);







  x1 = floor(length(lon) / 2); x1 = max(1,x1); x2 = min(length(lon),x1+1);
  y1 = floor(length(lat) / 2); y1 = max(1,y1); y2 = min(length(lat),y1+1);

  % Use MKS units [but time in index units] for gradient and Laplacian
  dt = 1;
  % dx = sw_dist(lat([y1 y1]),lon([x1 x2]),'km')*1e3;
  % dy = sw_dist(lat([y1 y2]),lon([x1 x1]),'km')*1e3;
  dy = [0 ; cumsum(sw_dist(midlat,tstr.lon,'km'))] .* 1e3;
  dx = [0 ; cumsum(sw_dist(tstr.lat,midlon,'km'))] .* 1e3;

  dx = [0 ; cumsum(sw_dist(lat(repmat(y1,size(lon))),lon,'km'))]*1e3;
  dy = [0 ; cumsum(sw_dist(lat,lon(repmat(x1,size(lat))),'km'))]*1e3;






  % We are only interested in the central point??
  delT = squeeze([ dTdx(midx,midy,:) dTdy(midx,midy,:) ]);


  %DEBUG:  v(:) = 0;








    % stn = station_calc_udotdelt(stn,qeufld,qevfld,Tfld,...
    %                             ['raw_' udTfld],udTfld,...
    %                             qtfld,qtAdvfld);









   case 'mlrf1', minori= 10; maxori= 45;




%%%%% EXPERIMENTS WITH INTERPOLATION (v. INTERP2)

% linear/linear GOM
% oris=50:55;

% % linear/linear FKEYS
% oris=50:55;
% oris=52:.2:54;
% oris=52.6;

% % linear/nearest FKEYS
% % oris=50:55;
% % oris=52:.2:54;
% oris=53.2;

% % nearest/linear FKEYS
% % oris=50:60;
% % oris=53:.2:55;
% oris=53.9;

% % nearest/nearest FKEYS
% oris=54:.1:55;
% oris=54.9;






  %DEBUG:
  res.dTdt_hour=res.dTdtq0+res.dTdthc;
  disp('Iterating HC'); tic,
  for ix=2:length(res.dTdt)
    res.dTdtq0(ix) = res.dTdtq0(ix) + res.dTdtq0(ix-1) + res.dTdthc(ix-1);
    res.dTdtx(ix) = res.dTdtx(ix) + res.dTdtx(ix-1) - res.dTdthc(ix-1);
    res.dTdx(ix)=(res.dTdtq0(ix)-res.dTdtx(ix))./res.dx(ix);
    res.dTdx(~isfinite(res.dTdx)) = 0;
%%%% DO WE ACTUALLY NEED TO RECALC res.B/uf/u AT EACH NEW HOUR ALSO?
    res.dTdthc(ix)=-R.*dt.*res.u(ix).*res.dTdx(ix);
    res.dTdt(ix)=res.dTdtq0(ix)+res.dTdthc(ix);
  end;
  %DEBUG:
  toc,










  res.dTdtq0_sum = cumsum(res.dTdtq0);
  res.dTdtx_sum = cumsum(res.dTdtx);
  res.dTdthc_sum = cumsum(res.dTdthc);

  res.dTdtq0_sum(2:end) = res.dTdtq0_sum(2:end) + res.dTdthc_sum(1:end-1);
  res.dTdtx_sum(2:end) = res.dTdtx_sum(2:end) - res.dTdthc_sum(1:end-1);

  res.dTdx_sum=(res.dTdtq0_sum-res.dTdtx_sum)./res.dx;
  res.dTdx_sum(~isfinite(res.dTdx_sum)) = 0;
  res.dTdt_sum=res.dTdtq0_sum+res.dTdthc_sum;

  res.dTdt = res.dTdt_sum;









%%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%

function val = sai__get_opt(opts,nm,dflt)
  if ( isfield(opts,nm) )
    val = opts.(nm);
  else
    val = dflt;
  end;
return;




%%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%

function val = sbe__get_opt(opts,nm,dflt)
  if ( isfield(opts,nm) )
    val = opts.(nm);
  else
    val = dflt;
  end;
return;




%%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%

function val = hc__get_opt(opts,nm,dflt)
  if ( isfield(opts,nm) )
    val = opts.(nm);
  elseif ( nargin > 2 )
    val = dflt;
  else
    val = [];
  end;
return;






': Accumulated heat advection '




  [nt,nr,nc] = size(stn.(fldnm).laplacian);
  midr = ceil(nr/2); midc = ceil(nc/2);
  stn.(rawkl).date = stn.(fldnm).date;
  stn.(rawkl).data = -k*3600*squeeze(stn.(fldnm).laplacian(:,midr,midc));








    switch (KMPFX),
     case 'fkeys_hycom',	stn = get_fkeys_hycom(stn,[],[],[],[],interpMethod);
      more off;
    end;







      switch (yr),
       case {1988,1992,1996,2000,2004,2008,2012,2016,2020}, 


    %DEBUG:
    for yr=2004:2004
      %DEBUG:
      for jd=1:2






  matfname = fullfile(datapath,[lower(stn.station_name) '_fkeys_hycom_' upper(interpMethod) '.mat']);




%%%% MLRF1 FKEYS w/o Benthic
% oris=52:58;
% oris=56.2:0.2:57.8;
% oris=56.6;

%%%% MLRF1 FKEYS w/ Benthic
oris=52:58;
oris=53.2:0.2:54.8;
oris=54.2;




%%%% DEBUG???
%    if ( ~isfield(stn.(Tfld),'gradient_x') )
      % Calculate gradients and field Laplacians
      stn = calc_field_terms(stn,Tfld);
%    end;




    % Calculate gradient and Laplacian on model fields
    for vix = 1:length(vars)
      var = vars{vix};
      varl = var(1);
      fld = ['fkeys_hycom_' flds{vix}];
      if ( ~isempty(strfind(fld,'_field')) )
        station = calc_field_terms(station,fld);
      end;
    end;








    % Modify STATION.fkeys_hycom_*_field structs before saving STATION struct
    % to MAT file: add Nx1 fields .lon and .lat based on actual model grid.
    station = query_fkeys_hycom_station_field_coords(station,'fkeys_hycom_u_field');
    station = query_fkeys_hycom_station_field_coords(station,'fkeys_hycom_v_field');
    station = query_fkeys_hycom_station_field_coords(station,'fkeys_hycom_seatemp_field');








  plot(1:n,[x.non_qlh,x.qlh,]);






    for yr=begyr:endyr
      %DEBUG:
      disp(yr);

      if ( yr>=2008 )
        expt = '902';
      else
        expt = '303';
      end;

      for jd=1:366
        for hr=0:6:18

          dt = datenum(yr,1,1,hr,0,0) + jd - 1;
          [ig,dtix] = min(abs(alldts - dt));

          % First make sure we have ALL variables for this date/time
          fpatt = fullfile(fkeyspath,sprintf('%s_archv.%04d_%03d_%02d_3z*.nc',expt,yr,jd,hr));
          flist = dir(fpatt);
          if ( isempty(flist) )
            % warning('Missing all variables for date/time "%s"!',fpatt);
            continue;
          elseif ( length(flist) < length(unique(vars)) )
            warning('Missing some variables for date/time "%s"!',fpatt);
            continue;
          end;

          for vix = 1:length(vars)
            var = vars{vix};
            varl = var(1);
            fld = ['fkeys_hycom_' flds{vix}];
            %DEBUG:            disp(fld);

            %303_archv.2004_001_06_3zs.nc
            fbasename = sprintf('%s_archv.%04d_%03d_%02d_3z%s.nc',expt,yr,jd,hr,varl);
            fname = fullfile(fkeyspath,fbasename);
            if ( ~exist(fname,'file') )
              warning('Missing file "%s"!',fname);
              continue;
            end;
            %DEBUG:            disp(fname);
            nc = mDataset(fname);
            if ( isempty(nc) )
              warning('Unable to open "%s"!',fname);
              continue;
            end;

            station.(fld).date(dtix,1) = dt;
            if ( ~isempty(strfind(flds{vix},'_field')) )
              station.(fld).field(dtix,1:(2*yrad)+1,1:(2*xrad)+1) = ...
                  squeeze(cast(nc{var}(1,1,yix-yrad:yix+yrad,xix-xrad:xix+xrad),'double'));
            else
              datsz = getShape(nc{var});
              % 2-D data elements
              if ( length(datsz) == 3 )
                station.(fld).data(dtix,1) = ...
                    squeeze(cast(nc{var}(1,yix,xix),'double'));
              % 3-D data elements
              else
                station.(fld).data(dtix,1) = ...
                    squeeze(cast(nc{var}(1,1,yix,xix),'double'));
              end;
            end;
            close(nc); clear nc;

          end;
        end;
      end;
    end;











  % n = 366;
  n = 52;

  % N = 24;
  N = 24*7;










     % case 'looe1',	stn.(orifld)=80;  %? Or 70 as in GET_LOOE1_ADCP???





    %scalings = {'SS','US','SU','UU'}





    %%%% ??? DEBUG: Base horizontal convection on surface heating only!
    q=bq0; qstr=bq0fld; commentstr = [commentstr ' HC(Q0+Qb) '];
    %%%% ??? DEBUG: Base horizontal convection on total (km-scale) budget
    % q=dT; qstr=bdTffld;





From STATION_BENTHIC_EXCHANGE.m:
  % Benthic heating from incident insolation. ASSUMES benthic "active layer"
  % under sea bed is ~10m deep, and contains "a" mix of marine sediment and
  % coral rock: rock density 2600 [kg/m3; Hughes 1987] specific heat capacity
  % ~0.2 times that of seawater; marine sediment is 70% sandy clay, 30% water
  % mix by mass; dry sandy clay has density ~2600, specific heat ~0.3 times
  % seawater. Our rho_b, Cp_b are volume averages of the above properties.
  global hb
  if ( isempty(hb) )
    disp([mfilename ': Using default thermal active depth 10m']);
    hb = 10;
  end;

  rho_rock = 2.6.*rhow;
  rho_sand = ((2.6*0.7) + (1.0*0.3)).*rhow;
  Cp_rock = 0.2.*Cpw;
  Cp_sand = ((0.3*0.7) + (1.0*0.3)).*Cpw;

  rhob = (frac*rho_sand) + ((1-frac)*rho_rock);
  Cpb = (frac*Cp_sand) + ((1-frac)*Cp_rock);

  % Boundary layer depth - substrate in thermal equilibrium with water
  global hbl
  if ( isempty(hbl) )
    disp([mfilename ': Using default benthic boundary layer depth 3cm']);
    hbl = 0.03; %[m]
  end;


  % Conversion factor for benthic temperature change: [W/m2] -> [K/hr]
  Kperhour = 3600./(rhob*Cpb*hb);

  % Parameters for benthic long-wave heat flux
  sigma = 5.67e-8;  % Stefan-Boltzman constant
  % Benthos emissivity: Oke (1978), Evans et al. (1998)
  epsb = 0.98;
  epsw = 0.96; % Seawater emissivity

  % Parameters for benthos-water convective heat flux
  % Drag coefficient (e.g., Thibodeaux & Boyle 1987)
  global Cbd
  if ( isempty(Cbd) )
    disp([mfilename ': Using default convective coefficient 3.8e-4']);
    Cbd = 3.8e-4;
  end;





From RECALCAB.m:
    %%%
    %% Radiative Fluxes

%%%% ??? DEBUG
clear global
% global Kd;
% global Ab;
% % for Kd=[0.05 0.10 0.20 0.30];
% % for Ab=[0.1 0.2:0.2:1.0];
% %% for Kd=[0.10];
% Kd=[0.30];
% Ab=[];
% % for Ab=[0.1350 0.2350 0.3350]; %0.2350=value calculated for FWYF1,MLRF1,SMKF1
% % for Ab=[0.0690 0.1690 0.2690]; %0.1690=value calculated for LONF1
% commentstr = [' (Kd:' num2str(Kd) ' Ab:' num2str(Ab) ') '];
%%%% ??? DEBUG

      stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld);







  plot(x.asr.date,cumsum(x.asr.data),x.coolf.date,cumsum(x.coolf.data));




  %% Commented out: Also loads .amodis_chlor_a, which we do not want!
  parfld = 'bic_surf_par';
  if ( ~isfield(stn,parfld) )
    if ( ~isfield(stn,'sea_t') )
      stn = load_station_data(stn);
    else
      disp('No BIC data!');
    end;
  end;





      chl = stn.amodis_chlor_a_1_hour_spline.data.*100;



      % Kd = 0.20 + 0.10.*cos((2.*pi.*get_yearday(stn.(qswfld).date(qswix))/366)-(pi/4));
      Kd = 0.25 + 0.20.*sin((2.*pi.*get_yearday(stn.(qswfld).date(qswix))/366));

      %DEBUG:
      disp('Using seasonal sine-modulated Kd = 0.05-0.45');
      Kd = 0.25 + 0.20.*sin((2.*pi.*get_yearday(stn.(qswfld).date(qswix))/366));



  figure;
  maxigraph;

  ax(1)=subplot(2,1,1); hold on;
  plot_ts(x.t); ylabel(ax(1),'T_s'); ylim([10,40]);

  ax(2)=subplot(3,1,2); plot_ts(x.qt); ylabel(ax(2),'Q_0'); ylim([10,40]);
  % ax(2)=subplot(3,1,2); plot_ts(x.bq0t); ylabel(ax(2),'Q_0+Q_b'); ylim([10,40]);

  ax(3)=subplot(3,1,3); hold on;
  % plot_ts(x.sr);  ylabel(ax(3),'Q_S_W');
  % plot_ts(x.lr);  ylabel(ax(3),'Q_L_W');
  % plot_ts(x.qlh); ylabel(ax(3),'Q_L_H');
  % plot_ts(x.qlh_24_hour_sum); ylabel(ax(3),'\Sigma_1_dQ_L_H');
  % plot_ts(x.qsh); ylabel(ax(3),'Q_S_H');
  % plot_ts(x.qrh); ylabel(ax(3),'Q_R_H');

  plot_ts(x.asr_24_hour_sum,x.lr_24_hour_sum,x.qlh_24_hour_sum,x.qsh_24_hour_sum);






  % plot_ts(x.sr);  ylabel(ax(3),'Q_S_W');
  % plot_ts(x.lr);  ylabel(ax(3),'Q_L_W');
  % plot_ts(x.qlh); ylabel(ax(3),'Q_L_H');
  % plot_ts(x.qlh_24_hour_sum); ylabel(ax(3),'\Sigma_1_dQ_L_H');
  % plot_ts(x.qsh); ylabel(ax(3),'Q_S_H');
  % plot_ts(x.qrh); ylabel(ax(3),'Q_R_H');









%%%% DEBUG???
figure; maxigraph; plot(stn.(qswfld).date(qswix),Kd); datetick3;





From ANLOOE1.m:
%function stn = anlooe1(stn,R,cfac,wfac,lagoff)
%
% Load in situ data for Looe Key spar (and meteorology data from nearby SMKF1
% SEAKEYS station), and appropriate climatology, reanalysis, and model data,
% then calculate a full heat budget based on LOOE1 MicroCAT sea temperature.
%
% Last Saved Time-stamp: <Wed 2011-03-16 15:50:49  lew.gramer>
.
.
.
  if ( ~isfield(stn,'ncep_lrf') )
    stn = get_ncep_station(stn,'narr');
  end;

  rhfld = 'ndbc_relhumid';
  shfld = 'ndbc_spechumid';
  insfld = 'ncep_dsrf';
  uswfld = 'ncep_usrf';
  swfld = 'ncep_srf';

  if ( ~isfield(stn,rhfld) )
    if ( isfield(stn,'ndbc_dew_t') && (length(find(isfinite(stn.ndbc_dew_t.data))) > 0) )
      stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t',rhfld);
    else
      error('SMKF1 data should have contained field .ndbc_dew_t?!');
    end;
  end;
  if ( ~isfield(stn,shfld) )
    stn = station_relhumid_to_spechumid(stn,'ndbc_air_t',rhfld,shfld);
  end;

  stn = station_bulk_longwave(stn,'ndbc_air_t',shfld, ...
                              'ndbc_barom','ncep_cloud_cover','ndbc_sea_t', ...
                              'ncep_cloud_cover','ndbc_bulk_rh_dlrf', ...
                              'ndbc_bulk_rh_ulrf','ndbc_bulk_rh_lrf');
  stn = station_ndbc_hfbulk(stn,'ncep_srf','ndbc_bulk_rh_lrf');

  [stn,Q0_factor] = ...
      station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
                        'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lrf',...
                        'ndbc_ncep_30a','ncep_dsrf','ncep_dlrf','ncep_precip',...
                        'ndbc_wind1_dir','default','default','ww3_peakwaveper','ww3_sigwavehgt',...
                        'default',true);

  % Try Florida Keys HYCOM (and Gulf of Mexico HYCOM too)
  stn = tryfkeys(stn,R,cfac,wfac,lagoff);






    % axis([datenum(2005,4,1),datenum(2005,5,1),16,28]);
    % axis('tight');
    % axis([datenum(yr,1,1),datenum(yr,12,31),0,60]);
    axis([datenum(yr,1,1),datenum(yr,3,31),0,26]);
    % ylim([-50 150]);
    % ylim([0 30]);
    % ylim([0 50]);
    %%%% ??? DEBUG







    % axis([datenum(yr,1,1),datenum(yr,3,31),14,25]);





    [ig,asrix] = intersect_dates(dts,stn.(asrfld).date);


    plot(stn.(asrfld).date(asrix),T0+cumsum(stn.(asr).data(asrix).*2.5e-8),stn.(udT_lpfld).date(udT_lpix),T0+cumsum(stn.(udT_lpfld).data(udT_lpix)),'y:','Color',[.8,.8,0]);






    qlhadj = 0.90; % FWYF1
    qlhadj = 1.00; % MLRF1
    qlhadj = 1.00; % LONF1
    qlhadj = 0.90; % SMKF1
    if ( qlhadj ~= 1.00 )
      adj_qlhfld = ['adj_' qlhfld];
      %%%% ??? DEBUG
      commentstr = [commentstr ' Q_L_H*' num2str(qlhadj) ' '],
      stn.(adj_qlhfld).date = stn.(qlhfld).date;
      stn.(adj_qlhfld).data = stn.(qlhfld).data .* qlhadj;
      [swix,lwix,lhix,shix,rhix] = ...
          intersect_all_dates([],stn.(asrfld).date,stn.(lrfld).date,...
                              stn.(adj_qlhfld).date,stn.(qshfld).date,stn.(qrhfld).date);
      stn.(q0fld).date = stn.(qlhfld).date(lhix);
      stn.(q0fld).data = stn.(asrfld).data(swix) + stn.(lrfld).data(lwix) ...
          + stn.(adj_qlhfld).data(lhix) + stn.(qshfld).data(shix) + stn.(qrhfld).data(rhix);
      % Algorithm sometimes returns complex numbers!
      stn.(q0fld).data = real(stn.(q0fld).data);
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],mhfld);
    end;








% %%%% MLRF1 SSU
% oris=90:95;
oris=93;

% MLRF1 FKEYS fixed gradient
% oris=50:55;
% oris=52:.1:53;
oris=52.3;

%%%% DEBUG??? MLRF1 FKEYS
  oris = 50:55;
  oris = 51.2:.2:52.8;
  oris = 52.0:.05:52.4;
  oris = 52.2:.015:52.35;
oris=52.25;

%%%% DEBUG??? MLRF1 GOM_HYCOM
  oris = 10:5:40;
  oris = 32:1:38;
oris=35;

%%%% DEBUG??? FWYF1 FKEYS
  oris=0:0.75:5;
  oris=1.5:0.125:2.25;
  oris=1.78:0.015:1.87;
oris=1.81;

%%%% DEBUG??? FWYF1 GOM_HYCOM
  oris=50:5:80;
oris=60;


%%%% DEBUG??? SMKF1 FKEYS_HYCOM
  oris=60:1:67;
  oris=64:0.4:66;
oris=65;

%%%% DEBUG??? SMKF1 GOM_HYCOM
  oris=75:5:90;
  oris=70:1:77;
oris=75;


%%%% DEBUG??? LONF1 FKEYS_HYCOM
  oris=60:1.5:70;
  oris=69:0.5:72;
  oris=68.1:0.1:69;
  oris=67.8:0.05:68.15;
  oris=67.5:0.05:67.8;
oris=67.5;

%%%% DEBUG??? LONF1 GOM_HYCOM
oris=67.5;














%  oris=60:0.5:63;
%  oris=55:0.5:59;





  %%%
  %% Variable-name prefixes ("PFX") for various input and output datasets

  % ISPFX - In situ (station) data
  if ( ~exist('ISPFX','var') || isempty(ISPFX) )
    ISPFX = 'ndbc';
  end;

  % TIDEPFX - Tidal heights and currents
  if ( ~exist('TIDEPFX','var') || isempty(TIDEPFX) )
    TIDEPFX = 'tmd_tide';
  end;

  % RAPFX - Atmospheric reanalysis
  if ( ~exist('RAPFX','var') || isempty(RAPFX) )
    RAPFX = 'erai';	%Or 'ncep'
  end;

  % TURPFX - Turbulent surface fluxes
  % With warm-layer adjustment
  TURPFX = [ISPFX '_' RAPFX '_30a'];
  % WITHOUT warm-layer adjustment (just cold-skin and rain)
  TUR30PFX = [ISPFX '_' RAPFX '_30'];

  % WAVEPFX - Surface waves
  if ( ~exist('WAVEPFX','var') || isempty(WAVEPFX) )
    WAVEPFX = 'ww3';
  end;

  % KMPFX - Kilometer-scale ocean currents and sea temperature
  if ( ~exist('KMPFX','var') || isempty(KMPFX) )
    KMPFX = 'fkeys_hycom';	%Or 'gom_hycom'
  end;

  % STOKESPFX - Surface and near-surface currents (mean, Stokes drift)
  STOKESPFX = [WAVEPFX '_' ISPFX '_stokes'];

  % QEPFX - Quasi-Eulerian currents and heat advection
  switch (KMPFX),
   case 'fkeys_hycom',	QEPFX = [WAVEPFX '_fkeys_qe'];
   case 'gom_hycom',	QEPFX = [WAVEPFX '_gom_qe'];
   otherwise,		error('Unknown km-scale model "%s"',KMPFX);
  end;



  %%%
  %% All station struct fieldnames used to produce heat budget 

  bathyfld = 'ngdc_92m_bathy';
  slopefld = 'ngdc_offshore_slope';

  % Tide data (or model)
  hfld = [TIDEPFX '_i_depth'];
  mhfld = ['mean_' hfld];
  tufld = [TIDEPFX '_u'];
  tvfld = [TIDEPFX '_v'];

  % Meteorology (in situ)
  afld = [ISPFX '_air_t'];
  sfld = [ISPFX '_sea_t'];

  %%%% Air pressure and humidity
  %%%% In situ records often missing or incomplete. Air-sea flux algorithms
  %%%% are quite insensitive to normal ranges in these, so use reanalysis.
  % pfld = [ISPFX '_barom'];
  % dfld = [ISPFX '_dew_t'];
  % rhfld = [ISPFX '_relhumid'];
  % qafld = [ISPFX '_spechumid'];
  pfld = [RAPFX '_barom'];
  switch (RAPFX),
   case 'ncep',		dfld = [RAPFX '_dewp'];
   otherwise,		dfld = [RAPFX '_dew_t'];
  end;
  rhfld = [RAPFX '_relhumid'];
  qafld = [RAPFX '_spechumid'];

  qsfld = [ISPFX '_sea_spechumid'];

  Wfld = [ISPFX '_wind1_speed'];
  Dfld = [ISPFX '_wind1_dir'];
  Ufld = [ISPFX '_wind1_u'];
  Vfld = [ISPFX '_wind1_v'];


  % Meteorology (gridded/reanalysis)
  rfld = [RAPFX '_precip'];
  cfld = [RAPFX '_cloud_cover'];
  pblzfld = [RAPFX '_pblz'];


  % Radiative surface fluxes
  dsrfld = [RAPFX '_dsrf'];
  usrfld = [RAPFX '_usrf'];
  srfld = [RAPFX '_srf'];
  asrfld = ['absorbed_' RAPFX '_srf'];
  gamfld = ['absorbed_' RAPFX '_gamma'];

  % BULK or REANALYSIS for long-wave radiation?
  dlrfld = [RAPFX '_' ISPFX '_dlrf'];
  ulrfld = [RAPFX '_' ISPFX '_ulrf'];
  lrfld = [RAPFX '_' ISPFX '_lrf'];
  % dlrfld = [RAPFX '_dlrf'];
  % ulrfld = [RAPFX '_ulrf'];
  % lrfld = [RAPFX '_lrf'];

  % Water-benthos fluxes
  qbfld = ['benthic_' RAPFX '_srf'];
  btfld = ['benthic_' RAPFX '_t'];
  qbofld = ['benthic_' RAPFX '_qbo'];


  % Turbulent and net surface fluxes
  qlhfld = [TURPFX '_latent_heat_flux'];
  qshfld = [TURPFX '_sensible_heat_flux'];
  qrhfld = [TURPFX '_rain_heat_flux'];

  q0fld = [TURPFX '_net_heat_flux'];
  qtfld = [q0fld '_term'];

  % Fluxes without warm-layer adjustment
  q030fld = [TUR30PFX '_net_heat_flux'];
  qt30fld = [q030fld '_term'];


  % Ocean processes
  whfld = [WAVEPFX '_sigwavehgt'];
  wpfld = [WAVEPFX '_peakwaveper'];
  wdfld = [WAVEPFX '_peakwavedir'];

  ufld = [KMPFX '_u'];
  vfld = [KMPFX '_v'];
  Tfld = [KMPFX '_seatemp_field'];

  % Hourly fit to native data
  hufld = ['hourly_' KMPFX '_u'];
  hvfld = ['hourly_' KMPFX '_v'];


  % Near-bottom currents (mean + tide): used both for benthic heat exchange
  % and (friction velocity-dependent) horizontal convection calculations
  netufld = [TIDEPFX '_' KMPFX '_u'];
  netvfld = [TIDEPFX '_' KMPFX '_v'];


  % Surface and near-surface currents (mean, Stokes drift)
  ssufld = [STOKESPFX '_u'];
  ssvfld = [STOKESPFX '_v'];
  sssfld = [STOKESPFX '_speed'];
  ssdfld = [STOKESPFX '_dir'];


  % Quasi-Eulerian currents and heat advection
  % Careful not to double-add advection by km-scale currents!
  qeufld = [QEPFX '_u'];
  qevfld = [QEPFX '_v'];
  qesfld = [QEPFX '_speed'];
  qedfld = [QEPFX '_dir'];

  udTfld = [QEPFX '_advected_heat'];

  qtAdvfld = [TURPFX '_' QEPFX '_qtadv'];

  % Km-scale heat diffusion
  switch (KMPFX),
   case 'fkeys_hycom',	K_theta = 2.5;
   case 'gom_hycom',	K_theta = 20;
   otherwise,		error('Unknown km-scale model "%s"',KMPFX);
  end;
  kd2Tfld = [KMPFX '_diffused_heat'];


  % Total budget
  dTfld = [TURPFX '_' QEPFX '_dt'];
  dTffld = [dTfld '_heat_flux'];

  % Total budget with benthic flux
  qbotfld = [qbofld '_term'];
  bq0fld = ['benthic_' q0fld];
  bq0tfld = [bq0fld '_term'];
  bdTfld = ['benthic_' dTfld];
  bdTffld = [bdTfld '_heat_flux'];

  % Total budget with benthic flux and horizontal convection
  HCPFX = 'hc_';










%%%% ??? DEBUG
global Kd;
if ( isempty(Kd) )
    Kd=0.10;
end;







    % startix = find(dts>=datenum(2005,4,1),1);
    % startix = find(dts>=datenum(2008,11,27),1);




% for Ab=[0 0.2350 1]; %0.2350=value calculated for FWYF1





  mean_bathy_depth = -mean(mean(stn.(bfld).field(ix-radx:ix+radx,jx-radx:jx+radx)));
%%%% DEBUG???





  matfname = fullfile(datapath,[lower(stn.station_name) '_heat_budget.mat']);






    dsr_lpfld = [dsrfld '_30_hour_lowpass'];
    par_lpfld = [parfld '_30_hour_lowpass'];



    dsr_lpfld = [dsrfld];
    par_lpfld = [parfld];



    %%%% ??? DEBUG
    if ( ~isfield(stn,slopefld) )
      stn = station_ngdc_offshore_slope(stn);
    end;





%%%% ??? DEBUG
stn.(hfld).data = stn.(hfld).data + 20;
commentstr = [commentstr ' (h+20)'],
%%%% ??? DEBUG


%%%% ??? DEBUG
    disp(['NOT Saving to ' matfname]);
%%%% ??? DEBUG









%%%% DEBUG:  tbase = 15;






    % Estimate heat flux into water column
    qbo(contigix) = qblwo(contigix) - qblwi(contigix) - qbsh(contigix) + qbcdu(contigix);
    % At each transition (gap, dawn, year-end), rebalance energy budget with
    % an UNPHYSICAL radiative "pulse" from the benthos into the water column!
    if ( contigix(end) < length(dts) )
%%%% ??? DEBUG:
      % qbo(contigix(end)) = (tb(contigix(end))-t(contigix(end)+1))./Kperhour;
%%%% ??? DEBUG:
      % %DEBUG:      qbo(contigix(end)) = nan;
    end;








  % hb = 2; %[m]
%%%% DEBUG:  hb = 6; %[m]


%%%% ??? DEBUG:
hb = 10;




  % hbl = 0.03; %[m]
%%%% DEBUG:  hbl = 0.2; %[m]
%%%% DEBUG:  hbl = 2; %[m]

%%%% ??? DEBUG:
  hbl = 0.3; %[m]




%%%% ??? DEBUG:
Cbd = 2e-4;
Cbd = 0.5e-4;




  % Cbd = 1.4e-3;		% Drag coefficient (e.g., Thibodeaux & Boyle 1987)
%%%% DEBUG:  Cbd = 0;		% Drag coefficient (???)
%%%% DEBUG:

%%%% ??? DEBUG:
global Cbd
% Cbd = 3.8e-4;
%%%% ??? DEBUG:
Cbd = 7e-4;








    % q=q0; qstr=q0fld;





    %%%
    %% Save results to MAT file for future reference
    disp(['DEBUG??? ***NOT SAVING*** to ' matfname]);
    % disp(['Saving to ' matfname]);
    % station = stn;
    % save(matfname,'station');
    % station = []; clear station;







  figure; maxigraph; plot_ts(stn.erai_dsrf,'go',stn.erai_srf,'ro',stn.absorbed_erai_srf,'b.');
  axis([datenum(2005,4,1),datenum(2005,4,12),-200,1000]); datetick3('x',2,'keeplimits');
  hold on;
  % xlim([datenum(2005,4,1),datenum(2005,4,12)]); datetick3('x',2,'keeplimits');
  % plot(stn.absorbed_erai_gamma.date,stn.absorbed_erai_gamma.data*100,'k');
  plot(stn.aidiag.theta.date(stn.aidiag.theta.data>-20),stn.aidiag.theta.data(stn.aidiag.theta.data>-20)*10,'yo');
  plot(stn.aidiag.sun_angle_correction.date,stn.aidiag.sun_angle_correction.data*100,'ko');
  plot(stn.aidiag.tau.date,stn.aidiag.tau.data*100,'ko','Color',[.5,.5,.5]);
  plot(stn.ndbc_sea_t.date,(stn.ndbc_sea_t.data-25).*600,'k-','color',[.5,.5,.5]);
  plot(stn.bic_surf_par.date+(0.5/24),stn.bic_surf_par.data.*0.473,'k-');
  titlename('After After After');





    %%%% ??? DEBUG    Kd=0.5;
  %%%% ??? DEBUG  Ab = 0.40;


  % NIR and UVB component (fully absorbed)
  gam = (1 - PAR_PER_INSOL);
  % Downward PAR absorbed
  gam = gam + (PAR_PER_INSOL.*(1 - tau));
  % PAR reflected from benthos and absorbed
  gam = gam + (PAR_PER_INSOL.*tau.*Ab.*(1-tau));

%%%% ??? DEBUG
gam(0>gam | gam>1) = 1;

  stn.(aifld).date = dts;
  stn.(aifld).data = gam .* qsw;

  % When the sun is too low for reasonable estimates, ignore insolation
  badix = find(theta <= 0.5);
%%%% ??? DEBUG
  badix = find(theta < -12);
  stn.(aifld).data(badix) = 0;











%%%% ??? DEBUG
  badix = find(theta < 0);





  xlim([datenum(2005,4,1),datenum(2005,4,12)]); datetick3('x',2,'keeplimits');
  figure; maxigraph; plot_ts(stn.tau); titlename('\tau');
  xlim([datenum(2005,4,1),datenum(2005,4,12)]); datetick3('x',2,'keeplimits');



%%%% ??? DEBUG
yds=yds-(0.1/24);



  % When the sun is too low for reasonable estimates, ignore insolation
  badix = find(theta < 2);
%%%% DEBUG???
  badix = find(theta <= 0.5);
  stn.(aifld).data(badix) = 0;








    hold on; plot(stn.erai_srf.date,stn.erai_srf.data,'r:','LineWidth',2);




    % axis([datenum(2004,1,1),datenum(2004,12,31),15,34]); datetick3('x',2,'keeplimits');





STATION_HEAT_BUDGET.m:
  alt_dfld = [RAPFX '_dew_t'];
  alt_rhfld = [RAPFX '_relhumid'];
  alt_qafld = [RAPFX '_spechumid'];
  alt_qsfld = [RAPFX '_sea_spechumid'];

  % DEBUG: ??? For now, ALWAYS use reanalysis humidities
  rhfld = alt_rhfld;
  qafld = alt_qafld;
  qsfld = alt_qsfld;




    if ( ~isfield(stn,dfld) || ~isfield(stn.(dfld),'date') || numel(stn.(dfld).date)<(numel(stn.(afld).date)/3) )
      if ( isfield(stn,rhfld) )
        disp([rhfld '->' dfld]);
        stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
      elseif ( isfield(stn,qafld) )
        disp([qafld '->' dfld]);
        stn = station_spechumid_to_relhumid(stn,afld,qafld,rhfld);
        stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
      elseif ( isfield(stn,alt_rhfld) )
        dfld = alt_dfld;
        rhfld = alt_rhfld;
        qafld = alt_qafld;
        qsfld = alt_qsfld;
        % Some reanalysis have dew_t, some have relhumid!
        if ( ~isfield(stn,dfld) )
          disp([rhfld '->' dfld]);
          stn = station_relhumid_to_dewp(stn,afld,rhfld,dfld);
        end;
      else
        error('Found no humidity/dew-point data (%s,%s,%s,%s)!',dfld,alt_dfld,rhfld,alt_rhfld);
      end;
    end;
    if ( ~isfield(stn,rhfld) )
      stn = station_dewp_to_relhumid(stn,afld,dfld,rhfld);
    end;
    if ( ~isfield(stn,qafld) )
      stn = station_relhumid_to_spechumid(stn,afld,rhfld,qafld);
    end;



      %station_stokes_drift(stn,ssfld,sdfld,sufld,svfld,wsfld,wdfld,hsfld,tpfld,tdfld)





    dTdt.date = stn.(sfld).date(1:end-1);
    dTdt.data = diff(stn.(sfld).data);
    gapix = find(diff(stn.(sfld).date) > (2/24));
    dTdt.date(gapix) = [];
    dTdt.data(gapix) = [];








    disp(['RE-saving to ' matfname]);
    station = stn;
    save(matfname,'station');
    station = []; clear station;




















HORIZONTAL_CONVECTION.m:
  % Assuming 1-to-1 correspondence
  res.u_SS=res.Qv_SS./h;
  % Assuming fit in Fig. 10, panel (c) of Monismith et al (2006)
  res.u_SS=(5.0*res.u_SS) - 0.05;
  res.u_SS(res.u_SS<0) = 0;

  res.u = res.u_SS;




  % Assuming fit in Fig. 10, panel (c) of Monismith et al (2006)
  res.u=(5.0*res.u) - (5.0*0.012) + 0.01;




  % % Assuming 1-to-1 correspondence
  % res.u=res.Qv_SS./h;

  % Assuming fit in Fig. 10, panel (c) of Monismith et al (2006)
  res.u=(3.75*res.Qv_SS./h) - 0.01;
  res.u(res.u<0) = 0;






res.dTdtq0.date,T0+cumsum(res.dTdtq0.data),




% [cf,hf]=contourf(q0,bet,dTdthc.*3600./dt);
% titlename('\partial_tT_H_C = u_H_C [m/s] ^. \partial_xT_H_C [K/m] vs. Q_0 and \beta');
[cf,hf]=contourf(q0,bet,dTdt.*3600./dt);
titlename('\partial_tT = (Q_0/\rhoC_ph) + u_H_C [m/s] ^. \partial_xT_H_C [K/m] vs. Q_0 and \beta');





% [cf,hf]=contourf(q0,bet,dTdt);



% dTdt=dTdthc;




%q0=-200; [dt,dx]=meshgrid([1:24]*3600,1:1000); dTdtq0=dt.*q0./(rhoCp.*h); dTdtx=dt.*q0./(rhoCp.*(h+(bet.*dx))); dTdx=(dTdtq0-dTdtx)./dx; figure; contourf(dt/3600,dx,dTdx,[-.001:-.001:-.02]); colorbar;

%[bet,q0]=meshgrid(0.005:.005:0.20,0:10:1000); B=(g.*alph.*abs(q0))./(rhoCp); u=(B.*h./bet).^(1/3); figure; contourf(q0,bet,u); colorbar;


% B = (g .* alph .* abs(q0)) ./ (rhoCp);
% uf = (B .* h) .^ (1/3);
% bets=[0.006:0.002:0.20];
% for dx = [100 1000];
% u = repmat(nan,size(bets));
% dTdt = repmat(nan,size(bets));
% for dt = [1 12 24]*3600;
% dTdtq0 = dt*q0/(rhoCp*h);
% for betix=1:length(bets)
%   bet = bets(betix);
%   % Qv = 1.00 .* ((uf .* (dt) ./ h).^(3/2));
%   Qv = 0.15 .* (bet.^(-1/3)) .* uf .* h;
%   % Qv = 0.15 .* bet .* (uf.^3) .* ((dt)^2) ./ (h.^2);
%   u(betix) = Qv/h;
%   %%%% ??? DEBUG  u(betix) = 0.01;
%   dTdtx = dt*q0/(rhoCp*(h+(bet*dx)));
%   dTdx = ( dTdtq0 - dTdtx ) / dx;
%   dTdthc = -dt*u(betix)*dTdx;
%   dTdt(betix) = dTdtq0 + dTdthc;
% end;
% figure; maxigraph; plot(bets,(dTdt*24*3600/dt)); xlabel('\beta'); ylabel('K/day'); titlename(['\partial_tT ' num2str(dx) ' ' num2str(dt)]);
% % figure; maxigraph; plot(bets,u); xlabel('\beta'); ylabel('m/s'); titlename('u_H_C');
% % disp(dTdtq0);
% end;
% end;






figure;
maxigraph;
hold on;
ax(1) = axes;
[cf,hf]=contourf(q0,bet,dTdt*24);
ax(2) = axes; % STUPID colorbar behavior!
cbh = colorbar('peer',ax(2)); ylabel(cbh,'\partial_tT [K/hr]');
cbh = colorbar('peer',ax(1)); ylabel(cbh,'\partial_tT [K/hr]');
% linkaxes(ax);
[cf,hf]=contourf(q0,bet,dTdt*24);
[c,h]=contour(q0,bet,u); %clabel(c,h,'Color','r');
[c,h]=contour(q0,bet,dTdx); %clabel(c,h,'Color','w');
xlabel('Q_0 [Wm^-^2]'); ylabel('\beta = \nablah');
titlename('\partial_tT = (Q_0/\rhoC_ph) + u_H_C [m/s] ^. \partial_xT_H_C [K/m] vs. Q_0 and \beta');










1;
% Sensitivity analysis for Horizontal Convection

%[dt,dx]=meshgrid([1:24]*3600,1:1000); dTdtq0=dt.*q0./(rhoCp.*h); dTdtx=dt.*q0./(rhoCp.*(h+(bet.*dx))); dTdx=(dTdtq0-dTdtx)./dx; figure; contourf(dt/3600,dx,dTdx,[-.001:-.001:-.02]); colorbar;
%[bet,q0]=meshgrid(0.005:.005:0.20,0:10:1000); B=(g.*alph.*abs(q0))./(rhoCp); u=(B.*h./bet).^(1/3); figure; contourf(q0,bet,u); colorbar;

%dt=3600; [bet,q0]=meshgrid(0.005:.005:0.20,0:10:1000); B=(g.*alph.*abs(q0))./(rhoCp); u=(B.*h./bet).^(1/3); dx=u*dt; dTdtq0=dt.*q0./(rhoCp.*h); dTdtx=dt.*q0./(rhoCp.*(h+(bet.*dx))); dTdx=(dTdtq0-dTdtx)./dx; dTdthc=-dt.*u.*dTdx; dTdt=dTdtq0+dTdthc; figure; contourf(q0,bet,dTdt); colorbar;



q0 = -200;
rhoCp = 4.09e6;
h = 3.55;
g = 9.8;
alph = 2.9e-4;
B = (g .* alph .* abs(q0)) ./ (rhoCp);
uf = (B .* h) .^ (1/3);
bets=[0.006:0.002:0.20];
for dx = [100 1000];
u = repmat(nan,size(bets));
dTdt = repmat(nan,size(bets));
for dt = [1 12 24]*3600;
dTdtq0 = dt*q0/(rhoCp*h);
for betix=1:length(bets)
  bet = bets(betix);
  % Qv = 1.00 .* ((uf .* (dt) ./ h).^(3/2));
  Qv = 0.15 .* (bet.^(-1/3)) .* uf .* h;
  % Qv = 0.15 .* bet .* (uf.^3) .* ((dt)^2) ./ (h.^2);
  u(betix) = Qv/h;
  %%%% ??? DEBUG  u(betix) = 0.01;
  dTdtx = dt*q0/(rhoCp*(h+(bet*dx)));
  dTdx = ( dTdtq0 - dTdtx ) / dx;
  dTdthc = -dt*u(betix)*dTdx;
  dTdt(betix) = dTdtq0 + dTdthc;
end;
figure; maxigraph; plot(bets,(dTdt*24*3600/dt)); xlabel('\beta'); ylabel('K/day'); titlename(['\partial_tT ' num2str(dx) ' ' num2str(dt)]);
% figure; maxigraph; plot(bets,u); xlabel('\beta'); ylabel('m/s'); titlename('u_H_C');
% disp(dTdtq0);
end;
end;









  hfld = ['BOGUS_' TIDEPFX '_i_depth'];
stn.(hfld) = stn.([TIDEPFX '_i_depth']);
stn.(hfld).data = stn.(hfld).data + 30;



%%%%% DEBUG???
Kd=0.50;





goodix = find(ismember(get_year(stn.(qtfld).date),[1993 1997 1998 2001 2004 2005]));



% Remove any year with ANY NaNs
nanix = find(isnan(newq.data));
badyrs = unique(get_year(newq.date(nanix)));
badix = find(ismember(get_year(newq.date),badyrs));
newq.date(badix) = [];
newq.data(badix) = [];








plot(1:(1/24):365,nanmean(acc));



if ( ~doAcc )
  % acc = dat - repmat(dat(:,1),[1 24]);
  acc(:,1:end-1) = dat(:,2:end) - repmat(dat(:,1),[1 23]);
  acc(:,end) = 0;
  acc(2:end,end) = dat(2:end,1) - dat(1:end-1,23);
end;
% if ( ~doAcc )
%   acc(:,1) = 0;
%   acc(2:end,1) = dat(2:end,1) - dat(1:end-1,23);
%   acc(:,2:24) = diff(dat,[],2);
% end;















  acc(2:end,1) = dat(2:end,1) - dat(1:end-1,23);





if ( ~doAcc )
  acc(:,1) = 0;
  acc(2:end,1) = dat(2:end,1) - dat(1:end-1,23);
  acc(:,2:24) = diff(dat,[],2);
end;





acc = cumsum(dat,2);
if ( ~doAcc )
  acc(:,1) = 0;
  % acc(2:end,1) = dat(2:end,1) - dat(1:end-1,23);
  acc(:,2:24) = diff(dat,1,2);
end;





if ( ~doAcc )
  acc = dat - repmat(dat(:,1),[1 24]);
end;




% if ( ~doAcc )
%   acc(1,1) = 0;
%   acc(2:end,1) = dat(2:end,1) - dat(1:end-1,23);
%   acc(:,2:24) = diff(dat,1,2);
% end;





1;

doAcc = true;

qtfld = 'ndbc_sea_t'; doAcc = false;

% qtfld = 'ndbc_hfbulk_heat_flux_term';
% qtfld = 'ndbc_erai_30a_heat_flux_term';
% qtfld = 'ndbc_erai_30a_ww3_fkeys_qe_dt';
% qtfld = 'benthic_ndbc_erai_30a_ww3_fkeys_qe_dt';
% qtfld = 'benthic_ndbc_erai_30a_ww3_fkeys_qe_dt_netqf';

% goodix = ts_boreal_cool(stn.(qtfld));
goodix = 1:length(stn.(qtfld).data);

[newq.date,newq.data] = gap_expand(stn.(qtfld).date(goodix),stn.(qtfld).data(goodix));


begix = find(isfinite(newq.data),1);
endix = find(isfinite(newq.data),1,'last');
newq.date = newq.date(begix:endix);
newq.data = newq.data(begix:endix);

begix = find(get_hour(newq.date)==0,1);
endix = find(get_hour(newq.date)==23,1,'last');
newq.date = newq.date(begix:endix);
newq.data = newq.data(begix:endix);

ndys = length(unique(floor(newq.date)));

dat = reshape(newq.data,[24 ndys])';

acc = 0;
acc(1:ndys,2:25) = cumsum(dat,2);
if ( ~doAcc )
  acc(1,1) = 0;
  acc(2:end,1) = dat(2:end,1) - dat(1:end-1,23);
  acc(:,2:24) = diff(dat,1,2);
  acc(:,25) = acc(:,1);
end;

figure;
maxigraph;
plot(0:24,nanmean(acc));
titlename([stn.station_name ' ' strrep(qtfld,'_','\_')]);












dat = reshape(newq.data,[ndys 24]);




                                        % [bdTfld '_heat_flux_24_hour_average'],...






  if ( isfield(stn,bdTffld) )
    plot_fluxes(stn,firstyr,1,dys,{sfld,btfld},[],{qtfld,dTfld,bdTfld,[bdTffld '_netqf_heat_flux']},[],...
                {'NDBC sea temperature','Modeled substrate temperature',...
                 '(Q_0 == \gammaQ_S_W + Q_L_W + Q_L_H + Q_S_H)/\rhoC_ph',...
                 'K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + Q_0/\rhoC_ph',...
                 'K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + (Q_0+Q_b)/\rhoC_ph',...
                 'HC( K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + (Q_0+Q_b)/\rhoC_ph )',...
                });
  else
    plot_fluxes(stn,firstyr,1,dys,{sfld},[],{qtfld,dTfld},[],...
                {'NDBC sea temperature',...
                 '(Q_0 == \gammaQ_S_W + Q_L_W + Q_L_H + Q_S_H)/\rhoC_ph',...
                 'K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + Q_0/\rhoC_ph',...
                });
  end;



% % Cbds=3.6e-5;
% Cbds=6e-5;



  % % Assume >50% of insolation is in NIR or UV, totally absorbed by seawater.
  % % Calculate (1 - Kd_PAR) absorption of the other <50%, in both directions.
  % % Assumes specular reflection of insolation by the sea-floor; also assumes
  % % any insolation not reflected by the sea-floor is immediately absorbed and
  % % conducted into the sea-bed FOREVER TO DISAPPEAR: therefore, a separate
  % % calculation of benthic-water column heat exchange *must* be done.
  % % gam = (1-PAR_PER_INSOL) + PAR_PER_INSOL*(1 - tau + tau*Ab(1-tau))
  % %  = 1 + PAR_PER_INSOL*(- tau + tau*Ab(1 - tau))
  % %  = 1 - PAR_PER_INSOL*tau*(1 - Ab*(1 - tau))
  % gam = 1 - (PAR_PER_INSOL .* tau .* (1 - (Ab.*(1 - tau))));
  gam = ( 1 - PAR_PER_INSOL ) + ( PAR_PER_INSOL.*(1 - tau + (tau.*Ab.*(1-tau))) );

  % % This assumed water absorbed *ALL* radiation not reflected to space, i.e.,
  % % gam = (1-PAR_PER_INSOL) + PAR_PER_INSOL*(1-tau+tau*Ab(1-tau) + tau(1-Ab))
  % %  = 1 - PPI*(1-tau+tau + tau*Ab-tau*Ab - tau*tau*Ab) = 1 - PPI*(1-tau^2Ab)
  % gam = 1 - PAR_PER_INSOL*(1 - (tau.*tau.*Ab));

  % % Original calculation which was in error, assumed PAR_PER_INSOL == 0.5
  % gam = ( 1 - (Ab.*tau.*(1 - tau)) );










  gam = 1 - (PAR_PER_INSOL .* tau .* (1 - (Ab.*(1 - tau)) - (tau*(1-Ab))));




  % If caller is interested in insolation absorbed by sea bed
  if ( exist('qbfld','var') && ~isempty(qbfld) )
    stn.(qbfld).date = dts;
    % NOTE: Sand ~0.20 times, sandy clay Cp ~0.3 times that of seawater Cp;
    % mean 60% sediment/water mix implies ~0.6 x sw Cp; mean sediment/sw mix
    % density ~2.5 x sw rho: heat capacity of sea floor ~1.5 x sw. Conduc.?
    stn.(qbfld).data = qsw .* tau .* (1 - Ab);
    stn.(qbfld).data(badix) = 0;
  end;




      % resid(hblix,hix,ix) = error_ts(stn.(bdTffld),dTdt,'rmse');
      resid(hblix,hix,ix) = error_ts(stn.(bdTffld),dTdt,'mem');




      % Estimate convective heat flux from benthos into water column
      qbsh(contigix(ix)) = rhow.*Cpw.*Cbd.*(spd(contigix(ix)).^2).*(t(contigix(ix)) - tb(contigix(ix)));






  rawufld = 
  rawvfld = ts_op(stn.(tvfld),stn.(vfld),'+');



      %DEBUG:      qbsh(contigix(ix)) = 0;	% No citations found for this process, or my "drag coefficient"




%%%% DEBUG:
  % epsb = 0.96;
  % epsw = 0.98; % Seawater emissivity



  % Benthic thermal conductivity, 50/50 mix: marine sediment (1; Nobes et al
  % 1986) and porous water-sorbed carbonate (2.8; Thomas, Frost, Harvey 1973)
%  Kb = 1.9;	% [W/m/K]
  Kb = 1.0;	% [W/m/K]



  Cbd = 1.4e-5;		% Drag coefficient (???)




  Cpb = mean([0.2 , ((0.3*0.7) + (1.0*0.3))]).* Cpw;



%%%% ???
  stn.(bdTfld) = ts_op(stn.(qtfld),stn.(qbotfld),'+');





%{
%}



    %DEBUG:    qbswi(contigix(49:end)) = 0;	% Try "turning off" insolation after first 48 hours





  hb = 2; %[m]
  rho_rock = 2.6.*rhow;
  rho_sand = ((2.6*0.7) + (1.0*0.3)).*rhow;
  Cp_rock = 0.2.*Cpw;
  Cp_sand = ((0.3*0.7) + (1.0*0.3)).*Cpw;

  rhob = mean([2.6 , ((2.6*0.7) + (1.0*0.3))]).* rhow;
  Cpb = mean([0.2 , ((0.3*0.7) + (1.0*0.3))]).* Cpw;









  QEPFX = 'ww3_fkeys_qe';





  QEPFX = 'ww3_fkeys_qe';


QEPFX = [WAVEPFX '_fkeys_qe'];

stn = get_gom_hycom(stn);

fkeys






  %DEBUG:
  disp(mfilename); disp(length(gapix)-1);




%%%%%%%%
%%%%%%%% Following code was post-loop-fudging (still reef-boiling) code


function stn = station_benthic_exchange(stn,tfld,ufld,vfld,qbfld,tbfld,qbofld)
%function stn = station_benthic_exchange(stn,tfld,ufld,vfld,qbfld,tbfld,qbofld)
%
% Calculate benthic temperature STN.(TBFLD) and heat flux into water column
% STN.(QBOFLD), from sea temperature STN.(TFLD), hourly currents UFLD,VFLD,
% and benthic absorbed short-wave radiation, STN.(QBFLD), assuming: initial
% equilibrium between benthic and water temperature; constant (mean) values
% for seabed density, heat capacity, long-wave emissivity, and conductivity.
%

  %DEBUG:
  tic,

  stn = verify_variable(stn,tfld);
  stn = verify_variable(stn,ufld);
  stn = verify_variable(stn,vfld);
  stn = verify_variable(stn,qbfld);

  [tix,uix,vix,qbix] = ...
      intersect_all_dates([],stn.(tfld).date,stn.(ufld).date,stn.(vfld).date,stn.(qbfld).date);

  dts = stn.(tfld).date(tix);
  t = stn.(tfld).data(tix);
  u = stn.(ufld).data(uix);
  v = stn.(vfld).data(vix);
  qb = stn.(qbfld).data(qbix);


  % Approximate seawater density [kg/m^3] and specific heat capacity [J/kg/K]
  rhow = 1.02e3;
  Cpw = 4e3;


  % NOTE: Many of the following assumptions were validated with literature
  % (as cited below), and with an inter-tidal mud study (Guarini et al 1997)

  % Benthic heating from incident insolation. ASSUMES benthic "active layer"
  % under sea bed is ~2m deep, and contains 50/50 mix of marine sediment and
  % coral rock: rock density 2600 [kg/m3; Hughes 1987] specific heat capacity
  % ~0.2 times that of seawater; marine sediment is 70% sandy clay, 30% water
  % mix by mass; dry sandy clay has density ~2600, specific heat ~0.3 times
  % seawater. Our rho_b and Cp_b are mass averages of the above properties.
  hb = 2; %[m]
  rhob = mean([2.6 , ((2.6*0.7) + (1.0*0.3))]).* rhow;
  Cpb = mean([0.2 , ((0.3*0.7) + (1.0*0.3))]).* Cpw;

  % Conversion factor: [W/m2] -> [K/hr]
  Kperhour = 3600./(rhob*Cpb*hb);

  % Parameters for benthic long-wave heat flux
  sigma = 5.67e-8;  % Stefan-Boltzman constant
  % Benthos emissivity: Oke (1978), Evans et al. (1998)
  epsb = 0.98;
  eps = 0.96; % Seawater emissivity

  % Parameters for benthos-water sensible heat flux
  Cbd = 1.4e-3;		% Drag coefficient (???)
  spd = uv_to_spd(u,v);	% Tidal/mean current above log layer


  % Parameters for benthic heat conduction

  % Benthic thermal conductivity, 50/50 mix: marine sediment (1; Nobes et al
  % 1986) and porous water-sorbed carbonate (2.8; Thomas, Frost, Harvey 1973)
  Kb = 1.9;	% [W/m/K]

  % Temperature at base of thermally active sediment/rock layer: based on
  % annual mean sea temperatures at LONF1,MLRF1,SMKF1,FWYF1 for 1992-2010,
  % according to the method of Golosov and Kirillin (2010)
  tbase = 26.5;

  % Benthos temperature
  tb = repmat(nan,size(t));
  %DEBUG:  disp(size(tb));


  % We have no sea-bed temperature data, so we assume an initial temperature
  % equilibrium between the sea-floor and the water column

  % Assume only an initial equilibrium (free-running heat-exchange model)
  gapix = [1 ; length(dts)];

  % % Assume initial equilibrium at the start of any >1 hour gap in data
  % gapix = [1 ; (find(diff(dts) > (1.1/24))+1) ; length(dts)];

  % % Assume initial equilibrium at the start of any multi-day gap in data
  % gapix = [1 ; (find(diff(dts) > 1)+1) ; length(dts)];

  % % Assume on January 01 each year, sea and benthos at equilibrium
  % gapix = [1 ; find(get_yearday(dts)==0 & get_hour(dts)==0) ; length(dts)];

  % % Assume at local dawn each day [=.10 GMT], sea and benthos at equilibrium
  % gapix = [1 ; find(get_hour(dts)==10) ; length(dts)];

  %DEBUG:
  disp(mfilename); disp(length(gapix)-1);
  for gapixix = 1:length(gapix)-1
    contigix = gapix(gapixix):gapix(gapixix+1);
    tb(contigix(1),1) = t(contigix(1));
    qbswi = qb(contigix(1:end-1));

    qblwi = eps.*sigma.*( (t(contigix(1:end-1))+273.14) .^4 );

    clear qblwo qbsh qbcdu qbcdd
    for ix = 1:length(contigix)-1
      % Estimate outgoing long-wave flux for each hour
      qblwo(ix,1) = epsb.*sigma.*( (tb(contigix(ix))+273.14) .^4 );
      % Adjust estimate for sensible heat flux from benthos into water column
      qbsh(ix,1) = Cbd.*(spd(contigix(ix)).^2).*(t(contigix(ix)) - tb(contigix(ix)));

      % Adjust estimate for heat conduction across sea-bed/water interface
      qbcdu(ix,1) = -Kb * ((t(contigix(ix)) - tb(contigix(ix)))/hb);

      % Adjust estimate for heat conduction through sediment/coral rock layer
      %   Qcd = K(dT/dz), see definitions of KB,TBASE above for assumptions
      qbcdd(ix,1) = -Kb * ((tb(contigix(ix)) - tbase)/hb);

      % Estimate next hour's benthic temperature based on total fluxes
      tb(contigix(ix+1),1) = tb(contigix(ix)) + ...
          (Kperhour.*( qbswi(ix) - qblwo(ix) + qblwi(ix) + qbsh(ix) - qbcdu(ix) + qbcdd(ix) ));
    end;

    % Heat flux into water column
    qbo(contigix(1:end-1),1) = qblwo - qblwi - qbsh + qbcdu;
    % At each transition (gap, dawn, year-end), rebalance energy budget with
    % an unphysical radiative "pulse" from the benthos into the water column!
    if ( contigix(end) < length(dts) )
      qbo(contigix(end),1) = (tb(contigix(end))-t(contigix(end)+1))./Kperhour;
      % %DEBUG:      qbo(contigix(end),1) = nan;
    end;
  end;

  stn.(tbfld).date = dts;
  stn.(tbfld).data = tb;
  stn.(qbofld).date = dts(1:end-1);
  stn.(qbofld).data = qbo;

  %DEBUG:
  toc,

return;










%%%%%%%%
%%%%%%%% Following code was loop-fudging, "working" (but reef-boiling) code

function stn = station_benthic_exchange(stn,tfld,ufld,vfld,qbfld,tbfld,qbofld)
%function stn = station_benthic_exchange(stn,tfld,ufld,vfld,qbfld,tbfld,qbofld)
%
% Calculate benthic temperature STN.(TBFLD) and heat flux into water column
% STN.(QBOFLD), from sea temperature STN.(TFLD), hourly currents UFLD,VFLD,
% and benthic absorbed short-wave radiation, STN.(QBFLD), assuming: initial
% equilibrium between benthic and water temperature; constant (mean) values
% for seabed density, heat capacity, long-wave emissivity, and conductivity.
%

  %DEBUG:
  tic,

  stn = verify_variable(stn,tfld);
  stn = verify_variable(stn,ufld);
  stn = verify_variable(stn,vfld);
  stn = verify_variable(stn,qbfld);

  [tix,uix,vix,qbix] = ...
      intersect_all_dates([],stn.(tfld).date,stn.(ufld).date,stn.(vfld).date,stn.(qbfld).date);

  dts = stn.(tfld).date(tix);
  t = stn.(tfld).data(tix);
  u = stn.(ufld).data(uix);
  v = stn.(vfld).data(vix);
  qb = stn.(qbfld).data(qbix);


  % Approximate seawater density [kg/m^3] and specific heat capacity [J/kg/K]
  rhow = 1.02e3;
  Cpw = 4e3;


  % NOTE: Many of the following assumptions were validated with literature
  % (as cited below), and with an inter-tidal mud study (Guarini et al 1997)

  % Benthic heating from incident insolation. ASSUMES benthic "active layer"
  % under sea bed is ~2m deep, and contains 50/50 mix of marine sediment and
  % coral rock: rock density 2600 [kg/m3; Hughes 1987] specific heat capacity
  % ~0.2 times that of seawater; marine sediment is 70% sandy clay, 30% water
  % mix by mass; dry sandy clay has density ~2600, specific heat ~0.3 times
  % seawater. Our rho_b and Cp_b are mass averages of the above properties.
  hb = 2; %[m]
  rhob = mean([2.6 , ((2.6*0.7) + (1.0*0.3))]).* rhow;
  Cpb = mean([0.2 , ((0.3*0.7) + (1.0*0.3))]).* Cpw;

  % Conversion factor: [W/m2] -> [K/hr]
  Kperhour = 3600./(rhob*Cpb*hb);

  % Parameters for benthic long-wave heat flux
  sigma = 5.67e-8;  % Stefan-Boltzman constant
  % Benthos emissivity: Oke (1978), Evans et al. (1998)
  epsb = 0.98;
  eps = 0.96; % Seawater emissivity

  % Parameters for benthos-water sensible heat flux
  Cbd = 1.4e-3;		% Drag coefficient (???)
  spd = uv_to_spd(u,v);	% Tidal/mean current above log layer


  % Parameters for benthic heat conduction

  % Benthic thermal conductivity, 50/50 mix: marine sediment (1; Nobes et al
  % 1986) and porous water-sorbed carbonate (2.8; Thomas, Frost, Harvey 1973)
  Kb = 1.9;	% [W/m/K]

  % Temperature at base of thermally active sediment/rock layer: based on
  % annual mean sea temperatures at LONF1,MLRF1,SMKF1,FWYF1 for 1992-2010,
  % according to the method of Golosov and Kirillin (2010)
  tbase = 26.5;

  % Benthos temperature
  tb = repmat(nan,size(t));
  %DEBUG:  disp(size(tb));


  % We have no sea-bed temperature data, so we assume an initial temperature
  % equilibrium between the sea-floor and the water column

  % % Assume only an initial equilibrium (free-running heat-exchange model)
  % gapix = [1 ; length(dts)];

  % % Assume initial equilibrium at the start of any >1 hour gap in data
  % gapix = [1 ; (find(diff(dts) > (1.1/24))+1) ; length(dts)];

  % Assume initial equilibrium at the start of any multi-day gap in data
  gapix = [1 ; (find(diff(dts) > 1)+1) ; length(dts)];

  % % Assume on January 01 each year, sea and benthos at equilibrium
  % gapix = [1 ; find(get_yearday(dts)==0 & get_hour(dts)==0) ; length(dts)];

  % % Assume at local dawn each day [=.10 GMT], sea and benthos at equilibrium
  % gapix = [1 ; find(get_hour(dts)==10) ; length(dts)];

  %DEBUG:
  tb(1,1) = t(1);
  disp(mfilename); disp(length(gapix)-1);
  for gapixix = 1:length(gapix)-1
    contigix = gapix(gapixix):gapix(gapixix+1);
    % tb(contigix(1),1) = t(contigix(1));
    qbswi = qb(contigix(1:end-1));

    qblwi = eps.*sigma.*( (t(contigix(1:end-1))+273.14) .^4 );

    clear qblwo qbsh qbcdu qbcdd
    for ix = 1:length(contigix)-1
      % Estimate outgoing long-wave flux for each hour
      qblwo(ix,1) = epsb.*sigma.*( (tb(contigix(ix))+273.14) .^4 );
      % Adjust estimate for sensible heat flux from benthos into water column
      qbsh(ix,1) = Cbd.*(spd(contigix(ix)).^2).*(t(contigix(ix)) - tb(contigix(ix)));

      % Adjust estimate for heat conduction across sea-bed/water interface
      qbcdu(ix,1) = -Kb * ((t(contigix(ix)) - tb(contigix(ix)))/hb);

      % Adjust estimate for heat conduction through sediment/coral rock layer
      %   Qcd = K(dT/dz), see definitions of KB,TBASE above for assumptions
      qbcdd(ix,1) = -Kb * ((tb(contigix(ix)) - tbase)/hb);

      % Estimate next hour's benthic temperature based on total fluxes
      tb(contigix(ix+1),1) = tb(contigix(ix)) + ...
          (Kperhour.*( qbswi(ix) - qblwo(ix) + qblwi(ix) + qbsh(ix) - qbcdu(ix) + qbcdd(ix) ));
    end;

    % Heat flux into water column
    qbo(contigix(1:end-1),1) = qblwo - qblwi - qbsh + qbcdu;
    % At each transition (gap, dawn, year-end), rebalance energy budget with
    % an unphysical radiative "pulse" from the benthos into the water column!
    if ( contigix(end) < length(dts) )
      qbo(contigix(end),1) = qbo(contigix(end-1),1);
      % qbo(contigix(end),1) = (tb(contigix(end))-t(contigix(end)+1))./Kperhour;
      % %DEBUG:      qbo(contigix(end),1) = nan;
    end;
  end;

  stn.(tbfld).date = dts;
  stn.(tbfld).data = tb;
  stn.(qbofld).date = dts(1:end-1);
  stn.(qbofld).data = qbo;

  %DEBUG:
  toc,

return;











  stn = compare_heat_budgets(stn,sfld,RAPFX,KMPFX); appendtitlename(' (sans benthic)');





  % %%%
  % %% Benthic heat exchanges
  % stn = station_benthic_exchange(stn,sfld,tufld,tvfld,qbfld,btfld,qbofld);
  % stn.(dTfld) = ts_op(stn.(dTfld),stn.(qbofld),'+');
  % stn = station_heat_flux_term_inverse(stn,[dTfld '_heat_flux'],...
  %                                      dTfld,sfld,[],hfld);
  % stn = compare_heat_budgets(stn,sfld,RAPFX,KMPFX); appendtitlename(' (with benthic)');





if ( ~isfield(stn,'erai_sea_spechumid') )
  stn.erai_sea_spechumid.date = stn.ndbc_sea_t.date;
  % 0.98 factor from Stommel - accounts for salinity
  stn.erai_sea_spechumid.data = 0.98 .* relhumid_to_spechumid(stn.ndbc_sea_t.data,100);
end;






  tbase = 25;
  tbase = 20;



      %DEBUG:    qbsh = 0;





  % Benthic heating due to incident insolation. ASSUMES: benthic "warm layer"
  % under sea bed contains 70% sandy clay, 30% seawater mix by mass and is
  % ~2m deep; dry sandy clay has density ~2.6 times seawater, specific heat
  % capacity ~0.3 times seawater. Our rho_b and Cp_b are thus mass averages.
  rhow = 1e3;
  Cpw = 4e3;

  hb = 2; %[m]
  rhob = ((2.6*0.7) + (1.0*0.3)).* rhow;
  Cpb = ((0.3*0.7) + (1.0*0.3)).* Cpw;

  hourlyDel = 3600./(rhob*Cpb*hb);







  Cbd = 1.4e2;





  %% For now, require ALL input arguments except DOPLOT

  % if (~exist('unm','var') || isempty(unm))
  %   unm='gom_hycom_u';
  % end;
  % if (~exist('vnm','var') || isempty(vnm))
  %   vnm='gom_hycom_v';
  % end;
  % if (~exist('fldnm','var') || isempty(fldnm))
  %   fldnm='gom_hycom_seatemp_field';
  % end;
  % if (~exist('rawud','var') || isempty(rawud))
  %   rawud='daily_gom_hycom_advected_heat';
  % end;
  % if (~exist('ud','var') || isempty(ud))
  %   ud='gom_hycom_advected_heat';
  % end;
  % if (~exist('htfld','var') || isempty(htfld))
  %   % htfld='netqf';
  %   htfld='ndbc_ncep_30a_heat_flux_term';
  % end;
  % if (~exist('dt','var') || isempty(dt))
  %   dt='gom_hycom_dt';
  % end;

  if (~exist('doPlot','var') || isempty(doPlot))
    doPlot = false;
  end;




  if ( doPlot )
    [firstyr,ig,ig] = datevec(stn.(dt).date(1));
    dys = ceil(stn.(dt).date(end)) - datenum(firstyr,1,1) + 1;
    plot_fluxes(stn,firstyr,1,dys);
  end;










function stn = station_calc_udotdelt(stn,unm,vnm,fldnm,rawud,ud,htfld,dt,doPlot)
%function stn = station_calc_udotdelt(stn,unm,vnm,fldnm,rawud,ud,htfld,dt,doPlot)
%
% Calculate dT/dt as the sum of net ocean surface heating (Q0/rho*h*Cpe) and
% advective heat flux terms, per eqn. (1) of Gramer, 2010.
%
% All arguments after STN struct are FIELD NAME strings:
% Inputs: UNM=ocean current U-component time series; VNM=current V; FLDNM=sea
%  temperature field struct; HTFLD=net surface heating term (Q0/rho*h*Cp).
% Outputs: RAWUD=low time-resolution advected heat term (e.g., for Global or
%  Gulf of Mexico HYCOM currents, this time series has one value per day);
%  UD=hourly advected heat term; DT=sum of hourly advected heat UD and HTFLD.
%
% Last Saved Time-stamp: <Tue 2011-01-25 17:47:08  lew.gramer>

  %% For now, require ALL input arguments except DOPLOT

  % if (~exist('unm','var') || isempty(unm))
  %   unm='gom_hycom_u';
  % end;
  % if (~exist('vnm','var') || isempty(vnm))
  %   vnm='gom_hycom_v';
  % end;
  % if (~exist('fldnm','var') || isempty(fldnm))
  %   fldnm='gom_hycom_seatemp_field';
  % end;
  % if (~exist('rawud','var') || isempty(rawud))
  %   rawud='daily_gom_hycom_advected_heat';
  % end;
  % if (~exist('ud','var') || isempty(ud))
  %   ud='gom_hycom_advected_heat';
  % end;
  % if (~exist('htfld','var') || isempty(htfld))
  %   % htfld='netqf';
  %   htfld='ndbc_ncep_30a_heat_flux_term';
  % end;
  % if (~exist('dt','var') || isempty(dt))
  %   dt='gom_hycom_dt';
  % end;

  if (~exist('doPlot','var') || isempty(doPlot))
    doPlot = false;
  end;

  [uix,vix,fldix] = intersect_all_dates([],stn.(unm).date,stn.(vnm).date,stn.(fldnm).date);

  % Calculate advected heat at native time resolution of data ([m/s]*[K/m] == [K/s])
  [udotdelT,dTdx,dTdy] = calc_udotdelt(stn.(unm).data(uix),stn.(vnm).data(vix),stn.(fldnm),fldix);

  stn.(rawud).date = stn.(unm).date;
  % Convert to units of [K/hr]
  stn.(rawud).data = -(3600*udotdelT)';

  % Spline-fit an hourly timer series
  stn.(ud).date = [stn.(unm).date(1):(1/24):stn.(unm).date(end)]';
  stn.(ud).data = spline(stn.(rawud).date,stn.(rawud).data,stn.(ud).date);
  stn = filter_gaps(stn,rawud,ud);

  [ix1,ix2] = intersect_dates(stn.(htfld).date,stn.(ud).date);
  stn.(dt).date = stn.(htfld).date(ix1);
  stn.(dt).data = stn.(htfld).data(ix1) + stn.(ud).data(ix2);

  if ( doPlot )
    [firstyr,ig,ig] = datevec(stn.(dt).date(1));
    dys = ceil(stn.(dt).date(end)) - datenum(firstyr,1,1) + 1;
    plot_fluxes(stn,firstyr,1,dys);
  end;

return;








function stn = station_benthic_exchange(stn,tfld,ufld,vfld,qbfld,tbfld,qbofld)
%function stn = station_benthic_exchange(stn,tfld,ufld,vfld,qbfld,tbfld,qbofld)

  stn = verify_variable(stn,tfld);
  stn = verify_variable(stn,ufld);
  stn = verify_variable(stn,vfld);
  stn = verify_variable(stn,qbfld);

  [tix,uix,vix,qbix] = ...
      intersect_all_dates([],stn.(tfld).date,stn.(ufld).date,stn.(vfld).date,stn.(qbfld).date);

  dts = stn.(tfld).date(tix);
  t = stn.(tfld).data(tix);
  u = stn.(ufld).data(uix);
  v = stn.(vfld).data(vix);
  qb = stn.(qbfld).data(qbix);


  % Benthic heating due to incident insolation. ASSUMES: benthic "warm layer"
  % under sea bed contains 70% sandy clay, 30% seawater mix by mass and is
  % ~1m deep; dry sandy clay has density ~2.6 times seawater, specific heat
  % capacity ~0.3 times seawater. Our rho_b and Cp_b are thus mass averages.
  rhow = 1e3;
  Cpw = 4e3;

  hb = 1; %[m]
  rhob = ((2.6*0.7) + (1.0*0.3)).* rhow;
  Cpb = ((0.3*0.7) + (1.0*0.3)).* Cpw;

  hourlyDel = 3600./(rhob*Cpb*hb);

  % Parameters for benthic long-wave heat flux
  sigma = 5.67e-8;  % Stefan-Boltzman constant
  epsb = 0.98; % Oke (1978), Evans et al. (1998)
  eps = 0.97;

  % Parameters for benthic sensible heat flux
  Cbd = 1.4e-3;
  spd = uv_to_spd(u,v);

  tb = repmat(nan,size(t));
  %DEBUG:  disp(size(tb));

  %% At local dawn each day [=.10 GMT], assume sea and benthos at equilibrium
  %gapix = [1 ; find(get_hour(dts) == 10) ; length(dts)];
  % Or use free-running heat exchange model assuming only initial equilibrium
  gapix = [1 ; (find(diff(dts) > (1.1/24))+1) ; length(dts)];

  %DEBUG:
  disp(length(gapix)-1);
  for ix = 1:length(gapix)-1
    contigix = gapix(ix):gapix(ix+1);
    tb(contigix(1)) = t(contigix(1));
    qbswi = qb(contigix(1:end-1));
    qblwi = eps.*sigma.*( (t(contigix(1:end-1))+273.14) .^4 );
    % First guess is based solely on incoming radiative fluxes
    tb(contigix(2:end)) = tb(contigix(1)) + (hourlyDel.*cumsum(qbswi+qblwi));
    qblwo = epsb.*sigma.*( (tb(contigix(1:end-1))+273.14) .^4 );
    % Adjust estimate based on all radiative fluxes
    tb(contigix(2:end)) = tb(contigix(1)) + (hourlyDel.*cumsum(qbswi-qblwo+qblwi));

    % Adjust estimate for sensible heat flux from benthos into water column
    qsh = Cbd.*(spd(contigix(2:end)).^2).*(tb(contigix(2:end))-t(contigix(2:end)));
    %DEBUG:    qsh = 0;
    tb(contigix(2:end)) = tb(contigix(2:end)) - (hourlyDel.*cumsum(qsh));

    % Also calculate heat flux into water column!
    qbo(contigix(1:end-1)) = qblwo - qblwi + qsh;
    % At each transition (gap, or dawn), rebalance energy budget with a
    % radiative "pulse" from benthos into water column
    if ( contigix(end) < length(dts) )
      qbo(contigix(end)) = (tb(contigix(end))-t(contigix(end)+1))./hourlyDel;
      %DEBUG:
      qbo(contigix(end)) = nan;
    end;
  end;

  stn.(tbfld).date = dts;
  stn.(tbfld).data = tb;
  stn.(qbofld).date = dts(1:end-1);
  stn.(qbofld).data = qbo;

return;











  [tix,uix,vix,qbix] = ...
      intersect_all_dates([],stn.(tfld).date,stn.(qbfld).date);
  dts = stn.(tfld).date(tix);
  t = stn.(tfld).data(tix);
  qb = stn.(qbfld).data(qbix);







  [tix,qbix] = intersect_dates(stn.(tfld).date,stn.(qbfld).date);
  dts = stn.(tfld).date(tix);
  t = stn.(tfld).data(tix);
  qb = stn.(qbfld).data(qbix);

  [ig,uix] = intersect_dates(dts,stn.(ufld).date);
  u = interp1(stn.(ufld).date(uix(1):uix(end)),stn.(ufld).data(uix(1):uix(end)),dts,'spline','extrap',nan);
  [ig,vix] = intersect_dates(dts,stn.(vfld).date);
  v = interp1(stn.(vfld).date(vix(1):vix(end)),stn.(vfld).data(vix(1):vix(end)),dts,'spline','extrap',nan);










  [tix,uix,vix,qbix] = ...
      intersect_all_dates([],stn.(tfld).date,stn.(ufld).date,stn.(vfld).date,stn.(qbfld).date);

  dts = stn.(tfld).date(tix);
  t = stn.(tfld).data(tix);
  u = stn.(ufld).data(uix);
  v = stn.(vfld).data(vix);
  qb = stn.(qbfld).data(qbix);






  hrs = get_hour(dts);
  for hr=[1:9 11:23]
    thishr = find(hrs == hr);
    lasthr = 
  stn.(tbfld).data(intraday) = 






  % stn.(qtAdvfld) = ts_op(stn.(qtfld),stn.(udTfld),'plus');

  % stn.(dTfld) = ts_op(stn.(qtAdvfld),stn.(kd2Tfld),'plus');






  % Assume 50% of insolation is in NIR, totally absorbed in water column.
  % Calculate (1 - Kd_PAR) absorption of the other 50% - in both directions.
  % Assumes specular reflection of insolation by the sea-floor; also assumes
  % any insolation not reflected by the sea-floor is immediately reradiated
  % as NIR - and therefore completely absorbed by the water column.
  tau = exp(-Kd .* z .* sun_angle_correction);








  gam = 0.5 + (0.5.*(1 - tau + (tau.*(Ab.*tau))));

%   stn.(aifld).date = dts;
%   stn.(aifld).data = qsw .* ( 1 - (Ab.*tau.*(1 - tau)) );

                        gam = ( 1 - (Ab.*tau.*(1 - tau)) );









  %%%
  %% 
  if ( ~isfield(stn,dTfld) )
    stn = station_calc_kdel2t(stn,K_theta,'fkeys_hycom_seatemp_field',...
                              'native_fkeys_hycom_diffused_heat',...
                              'fkeys_hycom_diffused_heat',...
                              'fkeys_hycom_qedt','fkeys_hycom_qelt');
    stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_qelt_heat_flux',...
                                         'fkeys_hycom_qelt',...
                                         'ndbc_sea_t',[],hfld);
  end;








  l = del2(rotfld,dx,dy,inf);
  l = permute(l,[3 1 2]);
  stn.(fnm).laplacian = l;




    % Stokes ocean drift estimated from wind and surface wave data
    result(stix) = station_stokes_drift(result(stix),...
                                        'erai_stokes_speed','erai_stokes_dir',...
                                        'erai_stokes_u','erai_stokes_v',...
                                        'erai_wind_speed','erai_wind_dir',...
                                        'erai_sigwavehgt','erai_peakwaveper','erai_peakwavedir');






  % Stokes ocean drift estimated from wind and surface wave data
  station = station_stokes_drift(station,...
                                 'erai_stokes_speed','erai_stokes_dir',...
                                 'erai_stokes_u','erai_stokes_v',...
                                 'erai_wind_speed','erai_wind_dir',...
                                 'erai_sigwavehgt','erai_peakwaveper','erai_peakwavedir');







  % Now spline-interpolate raw 3-hourly data into hourly fields
  for ddix=1:length(dds)
    dd = dds(ddix);
    dd = dd{:};
    flds = dd{4};
    for fldix = 1:length(flds)
      fld = ['erai_' flds{fldix}];
      rawfld = ['raw_' fld];

      rawdts = result(1).(rawfld).date;
      dts = [rawdts(1):(1/24):rawdts(end)]';
      ndts = length(dts);
      for stix=1:length(stnms)
        rawdat = result(stix).(rawfld).data;
        result(stix).(fld).date(1:ndts,1) = dts(:);
        result(stix).(fld).data(1:ndts,1) = interp1(rawdts,rawdat,dts(:),'spline');
        result(stix) = filter_gaps(result(stix),rawfld,fld,[],(6/24),[],nan);
      end;
    end;
  end;








function dat = get_erai_station_apply_cvt(cvt,dts,dat)
  [op,valstr] = strtok(cvt);
  val = str2num(valstr);
  switch ( op ),
   case 'add',
    dat = dat + val;
   case 'mult',
    dat = dat .* val;
   case 'deaccum',
    resetix = find(get_hour(dts) == 3 | get_hour(dts) == 15);
    resets = dat(resetix);
    dat=[0;diff(dat)]./(3*3600);
dat(2:end+1,:,:) = dat;
    dat(resetix) = resets./(3*3600);
   case 'none',
   otherwise,
    error('Unrecognized conversion string "%s"',cvt);
  end;
return;





            % result(stix).(fld).data(end+1:end+ndts,1) ...
            %     = dat(:,result(stix).erai_latix,result(stix).erai_lonix);
            %% Stupid INTERP2/INTERP3 is just too hard to use
            xix = result(stix).erai_lonix;
            yix = result(stix).erai_latix;
            xerr = result(stix).lon - lons(xix);
            yerr = result(stix).lat - lats(yix);
            result(stix).(fld).data(end+1:end+ndts,1) ...
                = interp_field(lats,lons,dat,




            result(stix).(fld).data(end+1:end+ndts,1) = ...
                interp3(DT,LAT,LON,dat,dts,result(stix).lat,result(stix).lon,'linear');



        hrs = cast( nc{'time'}(:,:,:), 'double');
        dts = datenum(yr,mo,1) + (hrs./24);
        ndts = numel(dts);
        [DT,LAT,LON] = ndgrid(dts,lats,lons);
        for varix = 1:length(vars)
          var = vars{varix};
          fld = ['raw_erai_' flds{varix}];
          cvt = cvts{varix};
          dat = cast( nc{var}(:,:,:), 'double' );
          dat = get_erai_station_apply_cvt(cvt,dts,dat);

          for stix=1:length(stnms)
            result(stix).(fld).date(end+1:end+ndts,1) = dts(:);
            % result(stix).(fld).data(end+1:end+ndts,1) = dat(:,result(stix).erai_latix,result(stix).erai_lonix);
            result(stix).(fld).data(end+1:end+ndts,1) = ...
                interp3(DT,LAT,LON,dat,dts,result(stix).lat,result(stix).lon,'linear');
          end; %for stix
        end; %for varix





    %DEBUG:    flds = dd{2};
        %DEBUG:        flds = dd{2};
    %DEBUG:    flds = dd{2};




nc = mDataset('\\cygnus\gramer\home\RSMAS\Coastal\thesis\ERA_Interim_201005_fc.grib');
nj_info(nc),
ehrs=cast(nc{'time'}(:,:,:),'double');
edts=datenum(2010,05,01)+(ehrs./24);
edat=cast(nc{'Surface_solar_radiation_downwards'}(:,yix,xix),'double');
edif=[0;diff(edat)]; edif(edif<0) = 0;
x.date=edts; x.data=edif; x.data(x.data<0)=0;
scatter_fit_ts(x,mlrf2.bic_surf_par)
edif=[0;diff(edat)]./(3*3600); edif(edif<0) = 0;
x.date=edts; x.data=edif;
scatter_fit_ts(x,mlrf2.bic_surf_par)
help station_par_to_insol
mlrf2 = station_par_to_insol(mlrf2,'bic_surf_par','bic_surf_dsrf','bic_surf_usrf','bic_surf_srf');
[mlrf2.lon,mlrf2.lat,mlrf2.depth]=get_station_coords('mlrf1')
[mlrf2.lon,mlrf2.lat,mlrf2.depth]=get_station_coords('mlrf1');
mlrf2
mlrf2 = station_par_to_insol(mlrf2,'bic_surf_par','bic_surf_dsrf','bic_surf_usrf','bic_surf_srf');
scatter_fit_ts(x,mlrf2.bic_surf_dsrf)




  dtstr = atrs.CoordinateModelRunDate;




From 'fc' (already in 'an'):
            'Charnock', ... %                      true      Charnock @ surface
            'Mean_sea_level_pressure', ... % Pa    true      Mean_sea_level_pressure @ surface
            'N10_metre_U_wind_component', ... % m s-1 true   N10_metre_U_wind_component @ surface
            'N10_metre_V_wind_component', ... % m s-1 true   N10_metre_V_wind_component @ surface
            'N2_metre_dewpoint_temperature', ... % K true    N2_metre_dewpoint_temperature @ surface
            'N2_metre_temperature', ... % K        true      N2_metre_temperature @ surface
            'Sea_surface_temperature', ... % K     true      Sea_surface_temperature @ surface
            'Skin_temperature', ... % K            true      Skin_temperature @ surface
            'Surface_pressure', ... % Pa           true      Surface_pressure @ surface
            'Total_cloud_cover', ... %             true      Total_cloud_cover @ surface



% 'fc':
%            'East-West_surface_stress', ... % N m-2 s true   East-West_surface_stress @ surface
%            'North-South_surface_stress', ... % N m-2 s true North-South_surface_stress @ surface




  if ( ~isfield(stn,'fkeys_hycom_qedt_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_qedt_heat_flux','fkeys_hycom_qedt',...
                                         'ndbc_sea_t',[],hfld);
  end;









  if ( ~isfield(stn,'gom_hycom_dt_netqf') )
    stn = station_calc_udotdelt(stn,[],[],[],[],[],'netqf','gom_hycom_dt_netqf');
  end;


  if ( ~isfield(stn,'netqf_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'netqf_heat_flux','netqf',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;
  if ( ~isfield(stn,'gom_hycom_netqf_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'gom_hycom_netqf_heat_flux','gom_hycom_netqf',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;



  if ( ~isfield(stn,'gom_hycom_qedt_netqf') )
    stn = station_calc_udotdelt(stn,'gom_hycom_quasi_eulerian_u','gom_hycom_quasi_eulerian_v',...
                                'gom_hycom_seatemp_field',...
                                'daily_gom_hycom_quasi_eulerian_advected_heat',...
                                'gom_hycom_quasi_eulerian_advected_heat',...
                                'netqf','gom_hycom_qedt_netqf');
  end;



  % Now play with adding heat diffusion term
  if ( ~isfield(stn,'gom_hycom_qelt_heat_flux') )
    stn = station_calc_kdel2t(stn,20,'gom_hycom_seatemp_field',...
                              'daily_gom_hycom_diffused_heat','gom_hycom_diffused_heat',...
                              'gom_hycom_qedt','gom_hycom_qelt');
    stn = station_heat_flux_term_inverse(stn,'gom_hycom_qelt_heat_flux','gom_hycom_qelt',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');

    stn = station_calc_kdel2t(stn,20,'gom_hycom_seatemp_field',...
                              'daily_gom_hycom_diffused_heat','gom_hycom_diffused_heat',...
                              'gom_hycom_qedt_netqf','gom_hycom_qelt_netqf');
    stn = station_heat_flux_term_inverse(stn,'gom_hycom_qelt_netqf_heat_flux',...
                                         'gom_hycom_qelt_netqf',...
                                         'ndbc_sea_t',[],'tmd_tide_i_depth');
  end;



  stn = station_calc_udotdelt(stn,[],[],[],[],[],'netqf','gom_hycom_dt_netqf');










  if ( ~isfield(stn,'fkeys_hycom_dt_netqf') )
    stn = station_calc_udotdelt(stn,'fkeys_hycom_u','fkeys_hycom_v',...
                                'fkeys_hycom_seatemp_field',...
                                'daily_fkeys_hycom_advected_heat',...
                                'fkeys_hycom_advected_heat',...
                                'netqf','fkeys_hycom_dt_netqf');
  end;

  if ( ~isfield(stn,'fkeys_hycom_dt_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_dt_heat_flux','fkeys_hycom_dt',...
                                         'ndbc_sea_t',[],hfld);
  end;






  if ( ~isfield(stn,'netqf_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'netqf_heat_flux','netqf',...
                                         'ndbc_sea_t',[],hfld);
  end;
  if ( ~isfield(stn,'fkeys_hycom_netqf_heat_flux') )
    stn = station_heat_flux_term_inverse(stn,'fkeys_hycom_netqf_heat_flux','fkeys_hycom_netqf',...
                                         'ndbc_sea_t',[],hfld);
  end;




  if ( ~isfield(stn,'fkeys_hycom_qedt_netqf') )
    stn = station_calc_udotdelt(stn,'fkeys_hycom_quasi_eulerian_u','fkeys_hycom_quasi_eulerian_v',...
                                'fkeys_hycom_seatemp_field',...
                                'daily_fkeys_hycom_quasi_eulerian_advected_heat',...
                                'fkeys_hycom_quasi_eulerian_advected_heat',...
                                'netqf','fkeys_hycom_qedt_netqf');
  end;










  if (0)
    multiplot_station(stn,{'ndbc_sea_t_diff','fkeys_hycom_qelnetqf','fkeys_hycom_seatemp_laplacian','fkeys_hycom_sgs_diffused_heat','fkeys_hycom_sgs_thermal_diffusivity'});
    xlim([datenum(2004,12,9,6,0,0),datenum(2004,12,19,18,0,0)]);
  end;







  for varix = 1:length(vars)





overs = (newq > sdiff);
unders = (newq < sdiff);
newq(overs) = 






  if ( nargout < 2 )
    result = []; clear result;
  end;





  matfname = fullfile(datapath,[stnm '_erai.mat']);
  if ( exist(matfname,'file') )
      if ( stix < length(stnms) )
        disp(['Station ' stnm ' MAT file already exists:' matfname]);
      end;





  matfname = fullfile(datapath,[stnm '_erai.mat']);
  if ( exist(matfname,'file') )
    disp(['Loading from ' matfname]);
    load(matfname,'station');
    station = []; clear station;
  end;






          rdat = permute(dat,[2 3 1]);





%DEBUG:disp({'s2w',nanmin(s2w),nanmean(s2w),nanmax(s2w),});


%DEBUG:
disp({'sQsh',nanmin(sQsh),nanmean(sQsh),nanmax(sQsh),});
%DEBUG:disp({'s2Qsh',nanmin(s2Qsh),nanmean(s2Qsh),nanmax(s2Qsh),});

%DEBUG:
disp({'sQlh',nanmin(sQlh),nanmean(sQlh),nanmax(sQlh),});
%DEBUG:disp({'s2Qlh',nanmin(s2Qlh),nanmean(s2Qlh),nanmax(s2Qlh),});


%DEBUG:
disp({'sQsw',nanmin(sQsw),nanmean(sQsw),nanmax(sQsw),});
%DEBUG:disp({'s2Qsw',nanmin(s2Qsw),nanmean(s2Qsw),nanmax(s2Qsw),});


%DEBUG:
disp({'sgammaQsw',nanmin(sgammaQsw),nanmean(sgammaQsw),nanmax(sgammaQsw),});
%DEBUG:disp({'s2gammaQsw',nanmin(s2gammaQsw),nanmean(s2gammaQsw),nanmax(s2gammaQsw),});


%DEBUG:
disp({'sQlw',nanmin(sQlw),nanmean(sQlw),nanmax(sQlw),});
%DEBUG:disp({'s2Qlw',nanmin(s2Qlw),nanmean(s2Qlw),nanmax(s2Qlw),});


%DEBUG:
disp({'sQ0',nanmin(sQ0),nanmean(sQ0),nanmax(sQ0),});






% r_th_w = 0.1239; % Representative value for QA data from SMKF1 for 1992-2010






W = sqrt((U.*U) + (V.*V));					W2 = W.^2;






  figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_sensible_heat_flux.date,stn.ndbc_ncep_30a_sensible_heat_flux.data,'k-'); plot(dts,sQsh,'r.'); datetick3; titlename('\sigmaQ_S_H');



  figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_latent_heat_flux.date,stn.ndbc_ncep_30a_latent_heat_flux.data,'k-'); plot(dts,sQlh,'r.'); datetick3; titlename('\sigmaQ_L_H');











th = s - a;							th2 = th.^2;
s2th = s2s + s2a - (2.*ss.*sa.*COVsa);				sth = sqrt(abs(s2th));

% r_th_w = 0.1239; % Representative value for QA data from SMKF1 for 1992-2010
[ig,r_th_w] = cov_ts(th,w);
%DEBUG:r_th_w = 0;

s2Qsh = airdens2.*Cpa2.*Cdr.*Cth.*( ...
    (U2.*((th2.*s2U) + s2th)) + (V2.*((th2.*s2V) + s2th)) + (Ug2.*((th2.*s2Ug) + s2th)) ...
    + ((2.*sw.*sth.*r_th_w)./(w.*th)) ...
    );





q = qs - qa;							q2 = q.^2;
s2q = s2qs + s2qa - (2.*qs.*qa.*COVqsqa);			sq = sqrt(abs(s2q));

[ig,r_q_w] = cov_ts(q,w);
%DEBUG:r_q_w = 0;

s2Qlh = airdens2.*Le2.*Cd.*Ce.*( ...
    (U2.*((q2.*s2U) + s2q)) + (V2.*((q2.*s2V) + s2q)) + (Ug2.*((q2.*s2Ug) + s2q)) ...
    + ((2.*sw.*sq.*r_q_w)./(w.*q)) ...
    );
s2Qlh(~isfinite(s2Qlh)) = 0;
sQlh = sqrt(abs(s2Qlh));















s2th = s2s + s2a - (2.*ss.*sa.*COVsa);				sth = sqrt(abs(s2th));




Cdew.*exp((+Cdew).*d).*sd

Cair.*exp((-Cair).*a).*sa





C = ((c1.*d) ./ (d + c2)) - ((c1.*a) ./ (a + c2));






    eu=sqrt(ue(1)^2+(ue(2)*u(i))^2);
    et=sqrt(te(1)^2+(te(2)*dt(j))^2);
    eq=sqrt(qe(1)^2+(qe(2)*dq)^2);
    xq1=ecq/sq/cqnh;
    xq2=eq/dq;
    fq=cqnh/k*dcq/sq;
    fs=cunh/k*dcu/su+a/3;
    xt1=ect/st/ctnh;
    xt2=et/dt(j);
    ft=ctnh/k*dct/st;
    fq=cqnh/k*dcq/sq;

   dwq(j,i)=sqrt(xq1^2+xq2^2+((fq*(1-a)+fs)^2*(xt1^2+xt2^2)+((1-ft)-2*fq)^2*(xs1^2+xs2^2))/d(j,i).^2);







  plot(dts(ix2),stn.ndbc_ncep_30a_latent_heat_flux.data(ix1)+sQlh(ix2),'r+');





  figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_sensible_heat_flux.date,stn.ndbc_ncep_30a_sensible_heat_flux.data,'k-'); plot(dts,sQsh2,'r.'); datetick3; titlename('\sigma^2Q_S_H');



  figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_latent_heat_flux.date,stn.ndbc_ncep_30a_latent_heat_flux.data,'k-'); plot(dts,sQlh2,'r.'); datetick3; titlename('\sigma^2Q_L_H');





% Maximal error estimates based on natural variance
% sqs = 0.0037;
sqs = 7e-5;
sqs2 = sqs.^2;
% sqa = 0.0038;
sqa = 1.4e-4;
sqa2 = sqa.^2;







  flds = grepstruct(s,'nocs_');


  s = annocs(stnm);
  for fld = {'nocs_latent_heat_flux','nocs_lrf','nocs_sensible_heat_flux', ...
             'nocs_srf','nocs_net_heat_flux','nocs_heat_flux_term',...
             'nocs_bulk_latent_heat_flux','nocs_bulk_sensible_heat_flux', ...
             'nocs_bulk_net_heat_flux','nocs_bulk_heat_flux_term',...
            }
    stn.(fld{:}) = s.(fld{:});
  end;
  s = []; clear s;



sQsh2 = 0;
sQsh2 = sQsh2 + (U2.*((th2.*sU2) + sth2));
sQsh2 = sQsh2 + (V2.*((th2.*sV2) + sth2));
sQsh2 = sQsh2 + (Ug2.*((th2.*sUg2) + sth2));
sQsh2 = sQsh2 + ((2.*sw.*sth.*r_th_w)./(w.*th));
sQsh2 = airdens2.*Cpa2.*Cdr.*Cth.*sQsh2;





sQsh2 = airdens2.*Cpa2.*Cdr.*Cth.*( ...
    (U2.*((th2.*sU2) + sth2)) + (V2.*((th2.*sV2) + sth2)) + (Ug2.*((th2.*sUg2) + sth2)) ...
    + ((2.*sw.*sth.*r_th_w)./(w.*th)) ...
    );






sw = (w.^(-.5)) .* sqrt( (2*U2.*sU2) + (2*U2.*sU2) + (2*U2.*sU2) + COVU2V2 + COVU2Ug2 + COVV2Ug2 );




%%%% DEBUG
return;
%%%% DEBUG





min(diff(stn.(tf).date(tfix)))/(1/24),
disp('DEBUG');
    [hfix2,tfix2,q0fix2] = intersect_all_dates( [], ...
        stn.(hf).date, stn.(tf).date, stn.(q0f).date );
min(diff(stn.(tf).date(tfix2)))/(1/24),
%%%% DEBUG
keyboard;
return;
%%%% DEBUG








    [hfix,tfix,q0fix] = intersect_all_dates( [], ...
        stn.(hf).date, stn.(tf).date, stn.(q0f).date );
%%%% DEBUG
return;
%%%% DEBUG







    [tfix,hfix,q0fix] = intersect_all_dates( [], ...
        stn.(hf).date, stn.(tf).date, stn.(q0f).date );
%%%% DEBUG
return;
%%%% DEBUG





  % First try Gulf of Mexico HYCOM
  stn = trygom(stn,R,cfac,wfac,lagoff,false);





%function stn = query_gom_hycom_station_field_coords(stn,fldnm)
  [tsz,ysz,xsz] = size(stn.(fldnm).field);
  begx = xix-floor(xsz/2);	begx = max(begx,1);
  endx = xsz-begx+1;		endx = min(endx,xsz);
  begy = yix-floor(ysz/2);	begy = max(begy,1);
  endy = ysz-begy+1;		endy = min(endy,ysz);







            [boguslonix,boguslatix] = gridnbhd_km(lons, lats,...
                                        dat.(stnm).lon, dat.(stnm).lat, 0);
            boguslonix = boguslonix - dat.minlonix + 1,
            boguslatix = boguslatix - dat.minlatix + 1,





          fname = fullfile(datapath, 'ww3', ...
                           sprintf('%s.%s.%04d%02d.grb', dataset, var, yr, mo));






      matfname = fullfile(datapath, 'ww3', ...
                          sprintf('ww3.%s.%04d%02d.mat', datafld, yr, mo));
      disp(['Saving to ' matfname]);
      save(matfname,'dat');  
      %DEBUG:
      toc,






            [boguslonix,boguslatix] = gridnbhd_km(lons, lats,...
                                        dat.(stnm).lon, dat.(stnm).lat, 0);
            boguslonix = boguslonix - dat.minlonix + 1,
            boguslatix = boguslatix - dat.minlatix + 1,






      %DEBUG:      minlonix,      minlatix,
      if ( ~isfield(dat,'minlonix') )
        dat.minlonix = 58;
        dat.minlatix = 94;
      end;

      for stix = 1:length(stns.lons)

        stnm = lower(stns.codes{stix});

        % NOTE: If we ever add any variables to the list above, change the
        % following line temporarily to simply read "if (1)"!
        if ( ~isfield(dat,stnm) )
          %DEBUG:          disp(['Adding ' stnm ' to ' matfname]);
          dat.(stnm).lon = stns.lons(stix);
          dat.(stnm).lat = stns.lats(stix);
          ididx = [];
          if ( ~isempty(stncfg) )
            ididx = find(strcmpi(stnm,stncfg{1}));
          end;
          % [lonix,latix] = gridnbhd_km(dat.lon, dat.lat,...
          %                             dat.(stnm).lon, dat.(stnm).lat, 0);
          if ( ~isempty(ididx) )
            latix = double(stncfg{2}(ididx));
            lonix = double(stncfg{3}(ididx));
          else
            [lonix,latix] = gridnbhd_km(lons, lats,...
                                        dat.(stnm).lon, dat.(stnm).lat, 0);
            %DEBUG:          fprintf('%s,%d,%d\n',stnm,latix,lonix);
          end;

          dat.(stnm).lonix = lonix - minlonix + 1;
          dat.(stnm).latix = latix - minlatix + 1;

%%%%DEBUG:
%             disp(['Station ' stnm]);
%             [lonix,latix] = gridnbhd_km(dat.lon, dat.lat,...
%                                         dat.(stnm).lon, dat.(stnm).lat, 0),
%             disp(['Vs. ']);
%             [dat.(stnm).lonix,dat.(stnm).latix],











    for ix=1:length(result.microcat_fnames)
      fname = fullfile(looepath,result.microcat_fnames{ix});
      [t,c,dt] = textread(fname, ' %f, %f, %[^\n]\n', 'headerlines',48);
      result.microcat_seatemp.date = datenum(dt);
      result.microcat_seatemp.data = t;
      result.microcat_cond.date = result.microcat_seatemp.date;
      result.microcat_cond.data = c;
      result.microcat_salin.date = result.microcat_seatemp.date;
      result.microcat_salin.data = sw_salt(t,c);
    end;









  if ( ~isfield(stn,'fkeys_hycom_u') )
    stn = get_fkeys_hycom(stn);
  end;



  if ( ~isfield(stn,'ww3_stokes_speed') )
    stn = station_stokes_drift(stn,'ww3_stokes_speed','ww3_stokes_dir','ww3_stokes_u',...
                               'ww3_stokes_v','ndbc_wind1_speed','ndbc_wind1_dir',...
                               'ww3_sigwavehgt','ww3_peakwaveper','ww3_peakwavedir');
  end;
  if ( ~isfield(stn,'fkeys_hycom_quasi_eulerian_speed') )
    stn = calc_quasi_eulerian(stn,'ww3_stokes','fkeys_hycom','fkeys_hycom_quasi_eulerian');
  end;










%   if ( ~isfield(stn,'fkeys_hycom_qelt_heat_flux') )
%   end;






  stn = filter_gaps(stn,'ndbc_sea_t','fkeys_hycom_sgs_diffused_heat');
  stn = filter_gaps(stn,'fkeys_hycom_qelnetqf','fkeys_hycom_sgs_diffused_heat');







  % pcp = pchip(rawdts,rawdat,dts);

  % figure; maxigraph; hold on;
  % plot(rawdts,rawdat,'k.');
  % plot(dts,dat,'r--');
  % plot(dts,pcp,'b:');
  % datetick3;
  % legend('FKEYS HYCOM \nabla^2T','Spline','PCHIp', 'Location','Best');







  nffld = 'ndbc_ncep_30a_net_heat_flux';
  badix = find(abs(stn.(nffld).data) > 2000);
  if ( ~isempty(badix) )
    warning('Deleting %d bad points from %s',length(badix),nffld);
    stn.(nffld).date(badix) = [];
    stn.(nffld).data(badix) = [];
  end;

  htfld = 'ndbc_ncep_30a_heat_flux_term';
  badix = find(abs(stn.(htfld).data) > 10);
  if ( ~isempty(badix) )
    warning('Deleting %d bad points from %s',length(badix),htfld);
    stn.(htfld).date(badix) = [];
    stn.(htfld).data(badix) = [];
  end;












  if ( ~isfield(stn,'gom_hycom_dt') )
    disp('Calculating GoM HYCOM results also...');
    stn = trygom(stn,R,cfac,wfac,lagoff);
    close all;
  end;







                  'gom_hycom_qelt',...
                  'fkeys_hycom_qelt',...
                  ...



                  ...
                  'GoM 4km HYCOM u^.\nablaT + \nabla^2T + Stokes + Q_0',...
                  'FKEYS 1km HYCOM u^.\nablaT + \nabla^2T + Stokes + Q_0',...










  % Now play with adding heat diffusion term
  stn = station_calc_kdel2t(stn,2.5,'gom_hycom_seatemp_field',...
                            'daily_gom_hycom_diffused_heat','gom_hycom_diffused_heat',...
                            'gom_hycom_qedt','gom_hycom_qelt');
  stn = station_heat_flux_term_inverse(stn,'gom_hycom_qelt_heat_flux','gom_hycom_qelt',...
                                       'ndbc_sea_t',[],'tmd_tide_i_depth');
  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],'tmd_tide_i_depth',...
                                      'gom_hycom_qelt_heat_flux_24_hour_average',...
                                      'gom_hycom_qelqvf','gom_hycom_qelqf',R,cfac,wfac,...
                                      'gom_hycom_qelt','gom_hycom_qelnetqf',lagoff);








  if ( doPlot )
    [firstyr,ig,ig] = datevec(stn.(lt).date(1));
    dys = ceil(stn.(lt).date(end)) - datenum(firstyr,1,1) + 1;
    plot_fluxes(stn,firstyr,1,dys);
  end;







  if ( ~all(isfield(stn.(fld),{'gradient_x','gradient_y','gradient_t','laplacian'})) )

    stn.(fld)=rmfield(stn.(fld),{'gradient_x','gradient_y','gradient_t','laplacian'});






  if ( ~exist('dtix','var') || isempty(dtix) )
    dtix = 3;
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = true;
  end;
  if ( ~exist('fld','var') || isempty(fld) )
    fld = 'gom_hycom_seatemp_field';
  end;






%function stn = query_fkeys_hycom_station_field_coords(stn,fldnm)
  stn.(fldnm).lon = lons(begx:endx)';
  stn.(fldnm).lat = lats(begy:endy)';







  [tsz,ysz,xsz] = size(stn.(fldnm).data);
  xrad = floor(xsz/2);
  yrad = floor(ysz/2);

  stn.(fldnm).lon = lons(xix-xrad:xix+xrad)';
  stn.(fldnm).lat = gom_hycom_lats(yix-yrad:yix+yrad)';

  stn.(fldnm).field = stn.(fldnm).data;
  stn.(fldnm) = rmfield(stn.(fldnm),'data');




  % yix = yix - 1; ????





    % elseif ( strcmp(flds{vix},'seatemp_field') )






function stn = get_gom_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
%function stn = get_gom_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
%
% Add fields for ocean surface U and V current components, sea temperature,
% salinity, mixed-layer depth, and 9x9 U, V, and sea temperature fields from
% assimilative NRL HYCOM+NCODA Gulf of Mexico 1/25 Degree Analysis (Prasad
% and Hogan 2007), to struct STN. If a station name string STNM (five chars.,
% e.g., 'mlrf1') is given instead of a struct, or if struct has no valid .lon
% and .lat fields (e.g., -80.38, 25.01), but has a .station_name field, then
% GET_STATION_COORDS is called, to try to retrieve the station's location.
%
% Optional VARS is a string or cellstr specifying which variables to extract.
% No online catalog seems to be available, but see the URL templates below
% for a hard-won list of all the variables available from this dataset. If
% VARS is specified, optional FLDS may also be specified as the corresp.
% list of field names to be added to struct STN, e.g., 'seatemp'.
%
% DEFAULT VARS and corresponding FLDS cell arrays are as follows:
%   vars = { 'u', 'v', 'mld', 'qtot',          'temperature', 'salinity', 'u',       'v',       'temperature',  };
%   flds = { 'u', 'v', 'mld', 'net_heat_flux', 'seatemp',     'salinity', 'u_field', 'v_field', 'seatemp_field' };
%
% NOTE: The string 'gom_hycom_' is prepended to each of the names in FLDS.
%
% DEFAULT for BASEURL:
%  'http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1'
% If you have another ocean model with a THREDDS interface, specify it here.
%
% SAMPLE URL (files are stored by full year):
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1/2003
%
% CALLS: MDATASET (netCDF-Java), GET_STATION_COORDS, QUERY_GOM_HYCOM_INDICES,
%        QUERY_GOM_HYCOM_STATION_FIELD_COORDS.
%
% Last Saved Time-stamp: <Wed 2010-11-24 13:25:07 Eastern Standard Time gramer>

  set_more off;

  datapath = get_thesis_path('../data');

  if ( ischar(stn_or_stnm) )
    stn.station_name = stn_or_stnm;
  elseif ( isstruct(stn_or_stnm) )
    stn = stn_or_stnm;
    if ( ~isfield(stn,'station_name') )
      error('Station STRUCT with no station_name field!');
    end;
  else
    error('First arg must either be station STRUCT or station name string!');
  end;

  if ( ~isfield(stn,'lon') )
    if ( isfield(stn,'station_name') )
      [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
    else
      error('Station specified without coordinates or name!');
    end;
  end;

    if ( ~exist('vars','var') || isempty(vars) )
      vars = { 'u', 'v', 'mld', 'qtot',          'temperature', 'salinity', 'u',       'v',       'temperature',  };
      flds = { 'u', 'v', 'mld', 'net_heat_flux', 'seatemp',     'salinity', 'u_field', 'v_field', 'seatemp_field' };
    end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = lower(vars);
  end;


  % What are the first and last days with available model data?
  %%%% ??? NOTE: These are currently hard-wired: will need modifying later!
  zerodt = datenum(2003,1,1) + 002 - 1;
  lastdt = datenum(2010,1,1) + 197 - 1;
  % Model outputs are once-daily time series - assume midnight GMT
  alldts = zerodt:lastdt;


  matfname = fullfile(datapath,[lower(stn.station_name) '_gom_hycom.mat']);

  % If we already did this before, just load MAT file and subset
  if ( exist(matfname,'file') )
    disp(['Reloading from MAT file ' matfname]);
    x = load(matfname,'stn');
    allflds = grepstruct(x.stn,'gom_hycom');
    for fldix = 1:length(allflds)
      fld = allflds{fldix};
      stn.(fld) = x.stn.(fld);
    end;
    x = []; clear x;

  else

    disp('Loading original data from online netCDF...');

    % Grid-point radii for time-series fields (e.g., seatemp_field)
    yrad = 4;
    xrad = 4;

    if ( ~exist('mindt','var') || isempty(mindt) )
      mindt = zerodt;
    end;
    if ( ~exist('maxdt','var') || isempty(maxdt) )
      maxdt = lastdt;
    end;

    if ( ~exist('baseurl','var') || isempty(baseurl) )
      %% "Legacy" THREDDS interface
      %% baseurl = 'http://tds.hycom.org/opendap/nph-dods/datasets/hycom/GOMl0.04/expt_20.1';

      % New THREDDS interface
      baseurl = 'http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1';
    end;

    ixes = find(mindt <= alldts & alldts <= maxdt);
    alldts = alldts(ixes);

    %[begyr,ig,ig] = datevec(alldts(1));
    %[endyr,ig,ig] = datevec(alldts(end));
    % NEW IMPLEMENTATION: First time through, get ALL data
    [begyr,ig,ig] = datevec(zerodt);
    [endyr,ig,ig] = datevec(lastdt);


    [yix,xix] = query_gom_hycom_indices(stn.lon,stn.lat);

    for yr=begyr:endyr

      dts = [];

      url = sprintf('%s/%04d',baseurl,yr);
      %DEBUG:
      disp(url);
      nc = mDataset(url);
      if ( isempty(nc) )
        warning('Invalid URL "%s"!',url);
        continue;
      end;

      for vix = 1:length(vars)
        var = vars{vix};
        fld = ['gom_hycom_' flds{vix}];
        %DEBUG:
        disp(fld);

        if ( isempty(dts) )
          dts = cast(getTimes(nc{var}),'double');
          ndts = numel(dts);
        end;
        if ( yr == begyr )
          if ( ~isempty(strfind(flds{vix},'_field')) )
            stn.(fld) = struct('date',[],'field',[]);
          else
            stn.(fld) = struct('date',[],'data',[]);
          end;
        end;

        stn.(fld).date(end+1:end+ndts,1) = dts';
        if ( ~isempty(strfind(flds{vix},'_field')) )
          stn.(fld).data(end+1:end+ndts,1:(2*yrad)+1,1:(2*xrad)+1) = ...
              squeeze(cast(nc{var}(1:end,1,yix-yrad:yix+yrad,xix-xrad:xix+xrad),'double'));
        else
          datsz = getShape(nc{var});
          % 2-D data elements
          if ( length(datsz) == 3 )
            stn.(fld).data(end+1:end+ndts,1) = ...
                squeeze(cast(nc{var}(1:end,yix,xix),'double'));
          % 3-D data elements - get the first (highest) Z index
          else
            stn.(fld).data(end+1:end+ndts,1) = ...
                squeeze(cast(nc{var}(1:end,1,yix,xix),'double'));
          end;
        end;
      end;

      close(nc); clear nc;

    end;

    % Calculate speed and direction from model U and V currents
    if ( isfield(stn,'gom_hycom_u') && isfield(stn,'gom_hycom_v') )
      stn.gom_hycom_speed.date = stn.gom_hycom_u.date;
      stn.gom_hycom_speed.data = uv_to_spd(stn.gom_hycom_u.data,stn.gom_hycom_v.data);
      stn.gom_hycom_dir.date = stn.gom_hycom_u.date;
      stn.gom_hycom_dir.data = uv_to_dir_curr(stn.gom_hycom_u.data,stn.gom_hycom_v.data);
    end;

    % Modify STN.gom_hycom_*_field structs: change field .data to
    % .field, and add Nx1 fields .lon and .lat, before saving to MAT.
    stn = query_gom_hycom_station_field_coords(stn,'gom_hycom_u_field');
    stn = query_gom_hycom_station_field_coords(stn,'gom_hycom_v_field');
    stn = query_gom_hycom_station_field_coords(stn,'gom_hycom_seatemp_field');

    disp(['Saving to MAT file ' matfname]);
    save(matfname,'stn');

  end;

  % Limit our result fields to only those dates we requested
  for vix = 1:length(vars)
    fld = ['gom_hycom_' flds{vix}];

    % Gross quality control
    if ( ~isfield(stn,fld) )
      warning('No field "%s" found after load!',fld);
    % elseif ( strcmp(flds{vix},'seatemp_field') )
    elseif ( ~isempty(strfind(flds{vix},'_field')) )
      ixes = find(-4 >= stn.(fld).field | stn.(fld).field >= 40);
      stn.(fld).field(ixes) = nan;

      ixes = find(alldts(1) <= stn.(fld).date & stn.(fld).date <= alldts(end));
      stn.(fld).date = stn.(fld).date(ixes);
      stn.(fld).field = stn.(fld).field(ixes,:,:);
    else
      ixes = find(-3000 <= stn.(fld).data & stn.(fld).data <= 3000);
      stn.(fld).date = stn.(fld).date(ixes);
      stn.(fld).data = stn.(fld).data(ixes);

      ixes = find(alldts(1) <= stn.(fld).date & stn.(fld).date <= alldts(end));
      stn.(fld).date = stn.(fld).date(ixes);
      stn.(fld).data = stn.(fld).data(ixes);
    end;
  end;

  set_more;

return;










  minlon = -98.00;
  dlon = 0.0400;
  maxlon = -76.40;
  if ( minlon-dlon > stn.lon || stn.lon > maxlon+dlon )
    error('Ecoforecasts:gomHYCOM:BadLon',...
          'Longitude %f outside value range [%f,%f]',stn.lon,minlon,maxlon);
  end;

  xix = round( (stn.lon - minlon) ./ dlon );

  lons = minlon:dlon:maxlon;







        % if ( strcmp(flds{vix},'seatemp_field') )



    l = permute(l,[3 1 2]);




function calc_gom_hycom_terms(stn,dtix,doPlot)

  if ( ~exist('dtix','var') || isempty(dtix) )
    dtix = 3;
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = true;
  end;

  dts = stn.gom_hycom_seatemp_field.date;
  lon = stn.gom_hycom_seatemp_field.lon;
  lat = stn.gom_hycom_seatemp_field.lat;
  fld = stn.gom_hycom_seatemp_field.field;

  % Use MKS units for gradient and Laplacian
  dx = sw_dist(lat([1 1]),lon([1 2]),'km')*1e3;
  dy = sw_dist(lat([1 2]),lon([1 1]),'km')*1e3;
  [gx,gy] = gradient(squeeze(fld(dtix,:,:)),dx,dy);
  l = del2(squeeze(fld(dtix,:,:)),dx,dy);

  if ( doPlot )
    figure; maxigraph;
    ax1=subplot(2,2,1); contourf(lon,lat,squeeze(fld(dtix,:,:))); title('T');
    ax2=subplot(2,2,2); contourf(lon,lat,gx); title('\nablaT^x');
    ax3=subplot(2,2,4); contourf(lon,lat,gy); title('\nablaT^y');
    ax4=subplot(2,2,3); contourf(lon,lat,l); title('\nabla^2T');
    suptitle(datestr(dts(dtix)));
    colorbar('peer',ax1); colorbar('peer',ax2); colorbar('peer',ax3); colorbar('peer',ax4);
  end;

return;







%xlim([datenum(2007,11,22),datenum(2007,12,31)]); datetick3('x',2,'keeplimits');




    % If we get a WEIRD FILE with HYPERSAMPLING - compress to half-hourly means
    if ( min(diff(dts)) < (0.4/24) )
      halfhour = round(dts/(0.5/24))*(0.5/24);
      newdts = dts(1):(0.5/24.0):dts(end);
      t = interp1(dts,t,newdts);
      c = interp1(dts,t,newdts);
      dts = newdts;
    end;






    % stn.microcat_salin.data(end+1:end+ndts) = sw_salt(t(:),c(:),4.9);




    % [t,c,dt] = textread(fname, ' %f, %f, %[^\n]\n', 'headerlines',48);



    rawlines = importdata(fname);
    lasthdr = strmatch('start sample number =',rawlines);
    if ( isempty(lasthdr) )
      warning('Header could not be skipped in "%s"',fname);
      continue;
    end;
    % rawlines = char(rawlines(lasthdr+1:end));
    rawlines = rawlines(lasthdr+1:end);
    C = textscan(rawlines,' %f, %f, %s\n');




    fid = fopen(fname,'r');
    if ( fid < 0 )
      warning('Skipping missing file "%s"',fname);
      continue;
    end;






    [t,c,dt] = textread(fname, ' %f, %f, %[^\n]\n', 'headerlines',48);
    dts = datenum(dt);
    badix = find(datenum(1987,1,1)>dts | dts>datenum(2010,1,1) ...
                 | 1>t | t>36 | 0.001>c | c>1.000);
    dts(badix) = [];
    t(badix) = [];
    c(badix) = [];
    ndts = length(dts);

    stn.microcat_seatemp.date(end+1:end+ndts) = datenum(dt);
    stn.microcat_seatemp.data(end+1:end+ndts) = t;
    stn.microcat_cond.date(end+1:end+ndts) = stn.microcat_seatemp.date;
    stn.microcat_cond.data(end+1:end+ndts) = c;
    stn.microcat_salin.date(end+1:end+ndts) = stn.microcat_seatemp.date;
    % stn.microcat_salin.data(end+1:end+ndts) = sw_salt(t,c);
    stn.microcat_salin.data(end+1:end+ndts) = c;







  %DEBUG:  disp('tic'); tic,
  % stn.adcp_pct_good = squeeze(mean(adcp_pct_good,1));
  stn.adcp_avg_eacnt = squeeze(mean(adcp_eacnt,1));
  max_eacnt = max(stn.adcp_avg_eacnt(:,20:36),[],2);
  %DEBUG:  toc,

  %DEBUG:  disp('tic'); tic,
  % If I were a MATLAB geek, I could probably figure out how to remove this loop!






    % Also do baroclinic averages - AFTER removal of spurious returns
    stn.adcp_speed.data(ix) = nanmean(stn.adcp_speed.prof(ix,goodix));
    stn.adcp_u.data(ix) = nanmean(stn.adcp_u.prof(ix,goodix));
    stn.adcp_v.data(ix) = nanmean(stn.adcp_v.prof(ix,goodix));
    stn.adcp_w.data(ix) = nanmean(stn.adcp_w.prof(ix,goodix));
    stn.adcp_err.data(ix) = nanmean(stn.adcp_err.prof(ix,goodix));
    % Don't forget to use vectorial averaging
    stn.adcp_dir.data(ix) = uv_to_dir_curr(stn.adcp_u.data(ix),stn.adcp_v.data(ix));









  % Also do baroclinic averages - AFTER removal of spurious returns
  stn.adcp_speed.data = nanmean(stn.adcp_speed.prof,2);
  stn.adcp_u.data = nanmean(stn.adcp_u.prof,2);
  stn.adcp_v.data = nanmean(stn.adcp_v.prof,2);
  stn.adcp_w.data = nanmean(stn.adcp_w.prof,2);
  stn.adcp_err.data = nanmean(stn.adcp_err.prof,2);
  % Don't forget to use vectorial averaging...
  stn.adcp_dir.data = uv_to_dir_curr(stn.adcp_u.data,stn.adcp_v.data);






  minlon = -80.18;
  maxlon = -80.00;







%   plot_fluxes(stn,2004,109,6.5*365);
%   plot_fluxes(stn,2006,45,4.7*365);
%   figure(fh);




% %%%% DEBUG
% firstyr = 2008;
% dys = 365;
% %%%% DEBUG
%   fh = plot_fluxes(stn,firstyr,1,dys);
%   appendtitlename(sprintf(' R:%g C:%g W:%g', R, cfac, wfac));
%   if ( lagoff ~= 0 )
%     appendtitlename(sprintf(' lag:%d', lagoff));
%   end;








function [yix,xix] = query_gom_hycom_indices(lon,lat)
%function [yix,xix] = query_gom_hycom_indices(lon,lat)
%
% Return y- and x-indices (zero-based) for Gulf of Mexico 1/25-degree HYCOM

  datapath = get_thesis_path('../data');

  minlon = -98;
  dlon = 0.0400;
  xix = round( (lon - minlon) ./ dlon );

  % minlat = 18.0916;
  % dlat = 0.0339;  ... dlat = 0.0381;
  % yix = round( (lat - minlat) ./ dlat );
  load(fullfile(datapath, 'gom_hycom_lats.mat'));
  [ig,yix] = min( abs(gom_hycom_lats - lat) );
  yix = yix - 1;

return;




function [yix,xix] = query_flkeys_hycom_indices(lon,lat)
%function [yix,xix] = query_flkeys_hycom_indices(lon,lat)
%
% Return y- and x-indices (zero-based) for Florida Keys 1/100-degree HYCOM

  datapath = get_thesis_path('../data');

  minlon = -83.36;
  dlon = 0.0100;
  xix = round( (lon - minlon) ./ dlon );

  % minlat = 22.775372;
  % dlat = ??? 0.0339;  ... dlat = 0.0381;
  % yix = round( (lat - minlat) ./ dlat );
  load(fullfile(datapath, 'flkeys_hycom_lats.mat'));
  [ig,yix] = min( abs(flkeys_hycom_lats - lat) );
  yix = yix - 1;

return;










  if ( exist(matfname,'file') )
    if (doOverwrite)
      warning(['Deleting old file ' matfname]);
      delete(matfname);
    else
    end;
    
  else
    disp(['Subsetting raw XYZ data to ' matfname]);

    % Find all points "inside" our bounding rectangle
    inix = find( (abs(stn.lon-lon)<=(xrad/(111.12*cosd(stn.lat)))) ...
                 & (abs(stn.lat-lat)<=(yrad/(111.12))) );

    lon = lon(inix);
    lat = lat(inix);
    depth = depth(inix);

    [LON,LAT] = meshgrid(unique(lon),unique(lat));
    result.ngdc_92m_bathy.lon = LON';
    result.ngdc_92m_bathy.lat = LAT';
    result.ngdc_92m_bathy.field = griddata(lon,lat,depth,LON',LAT');
    save(matfname,'result');
  end;

  stn.ngdc_92m_bathy = result.ngdc_92m_bathy;
  result = []; clear result;

  if ( doPlot )
    figure;
    contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field);
    % set(gca,'ydir','rev');
    maxigraph;
    colorbar;
    titlename(['NDGC bathymetry surrounding ' stnm]);
  end;

return;






    % Find all points "inside" our bounding ellipse: ugh, geography sucks
    smaecc(1) = xrad;
    smaecc(2) = axes2ecc(xrad,yrad);
    % WGS 84 reference ellipsoid for the Earth, in [km]
    wgs84_ell = [6356.752 (1/298.25722356)];
    [elllat,elllon] = ellipse1(stn.lat,stn.lon,smaecc,90,[],wgs84_ell,[],100);
    inix = find(inside(lon,lat, elllon,elllat) > 0);






    ddeg = sqrt(((stn.lon - lon).^2) + ((stn.lat - lat).^2));
    inix = find(ddeg < (rad
    dkm = sw_dist([
    inix = find(






  % nlon=length(unique(LON));
  % nlat=length(unique(LAT));
  % stn.ngdc_92m_bathy.field = reshape(dat',[nlon nlat]);
  stn.ngdc_92m_bathy.field = griddata(lon,lat,dat,LON',LAT');





      clim.landy_lh = nc{'Q_lat'}(:,latix,lonix);
      clim.landy_lr = nc{'Q_lwdn'}(:,latix,lonix) - nc{'Q_lwup'}(:,latix,lonix);
      clim.landy_sh = nc{'Q_sen'}(:,latix,lonix);
      clim.landy_sr = nc{'Q_swnet'}(:,latix,lonix);
      clim.landy_ev = nc{'F_evap'}(:,latix,lonix);
      clim.landy_pr = nc{'F_prec'}(:,latix,lonix);
      clim.landy_ro = nc{'F_roff'}(:,latix,lonix);
      clim.landy_taux = nc{'taux'}(:,latix,lonix);
      clim.landy_tauy = nc{'tauy'}(:,latix,lonix);






  matfname = fullfile(datapath,sprintf('%s_nocs.mat',stnm));

  if ( exist(matfname,'file') )

    disp(['Loading NOCS data from ' matfname]);
    x = load(matfname,'stn');
    result = x.stn;
    x = []; clear x;

  else







  flds = {}; mns = [];
  fldbasenms = {'ncep_srf','ncep_lrf', ...
                'monthly_nocs_sensible_heat_flux','monthly_nocs_latent_heat_flux', ...
                'landy_sh','landy_lh_flux'};
  for ix = 1:length(fldbasenms)
    flds{ix} = [fldbasenms{ix}];
    mns(ix,1:52) = grpstats(real(stn.(flds{ix}).data),get_week(stn.(flds{ix}).date));
    mns(ix,1:52) = grpstats(real(stn.(flds{ix}).data),get_week(stn.(flds{ix}).date));
  end;
  figure; maxigraph; hold on;
  plot(1:52,mns); xlim([1 52]); legend(strrep(flds,'_','\_'));
  title('NCEP NARR Radiative Fluxes');

  flds = {}; mns = [];
  fldbasenms = { 'ndbc_ncep_30a_sensible_heat_flux','ndbc_ncep_30a_latent_heat_flux', ...
                 'monthly_nocs_sensible_heat_flux','monthly_nocs_latent_heat_flux', ...
                 'landy_sh','landy_lh_flux' };
  for ix = 1:length(fldbasenms)
    flds{ix} = [fldbasenms{ix}];
    mns(ix,1:52) = grpstats(real(stn.(flds{ix}).data),get_week(stn.(flds{ix}).date));
  end;
  figure; maxigraph; hold on;
  plot(1:52,mns); xlim([1 52]); legend(strrep(flds,'_','\_'));
  title('TOGA-COARE 3.0a Turbulent Fluxes');

  flds = {}; mns = [];
  fldbasenms = {'gom_hycom_advected_heat','netqf','gom_hycom_advected_heat'};
  for ix = 1:length(fldbasenms)
    flds{ix} = [fldbasenms{ix}];
    mns(ix,1:52) = grpstats(real(stn.(flds{ix}).data),get_week(stn.(flds{ix}).date));
  end;
  figure; maxigraph; hold on;
  plot(1:52,mns); xlim([1 52]); legend(strrep(flds,'_','\_'));
  title('GoM HYCOM, thermal siphon Advective Fluxes');








    matfname = fullfile(datapath,sprintf('%s_ms_PROBLEM.mat',stnm));




  dys = ceil(stn.(dt).date(end)) - datenum(firstyr,1,1) + 1;





  % fld = 'ndbc_ncep_30a_heat_flux_term';
  fld = 'ndbc_ncep_30a_net_heat_flux';





    plot(dts,real(dat),cncspec{mod(cix-1,length(cncspec))+1});
    plot(dts(1:96:end),real(dat(1:96:end)),mrkspec{mod(cix-1,length(mrkspec))+1});

    plot(dts,real(dat),cncspec{mod(cix-1,length(cncspec))+1});
    plot(dts(1:96:end),real(dat(1:96:end)),mrkspec{mod(cix-1,length(mrkspec))+1});

    plot(dts,real(dat),cncspec{mod(cix-1,length(cncspec))+1});
    plot(dts(1:96:end),real(dat(1:96:end)),mrkspec{mod(cix-1,length(mrkspec))+1});






  dt = repmat(1,size(tstr.field,1));






  midx = round(length(tstr.lon) / 2);
  midy = round(length(tstr.lat) / 2);

  dx = sw_dist(tstr.lat([midy midy]),tstr.lon([midx midx+1]),'km');
  dy = sw_dist(tstr.lat([midy midy+1]),tstr.lon([midx midx]),'km');






  dxdeg = mean(diff(tstr.lon(:)));
  dydeg = mean(diff(tstr.lat(:)));





%function [sst, LONS, LATS] = ansst(bbox_or_stanm, yr, wk, dataset, region)
  % Store/retrieve data in this M-file's local directory
  [pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = fullfile(pathroot, '../data', '');
  end;


  % Save images for future reference
  fpath = fullfile(avhrrpath, fname);








%function [sst, LONS, LATS] = query_avhrr_weekly_subset(bbox_or_stanm, yr, wk, dataset, region)


  % Store/retrieve data in this M-file's local directory
  [pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = fullfile(pathroot, '../data', '');
  end;




  % Save images for future reference
  fpath = fullfile(datapath, fname);
  if ( exist(fpath, 'file') )
    try
      sstbytes = imread(fpath);
    catch
      % Last ditch effort - delete corrupt file and load from URL anyway
      fprintf(2, '\n');
      warning('Corrupt cached file? Trying URL "%s"...', url);
      delete(fpath);
    end;
  end;
  if ( ~exist(fpath, 'file') )
    %DEBUG:    url,
    [fpath, fstatus] = urlwrite(url, fpath);
    if ( fstatus == 0 )
      fprintf(2, '\n');
      warning('SKIPPED! Failed to download "%s"...', url);
      return;
    end;
    sstbytes = imread(fpath);
  end;









%function dts = query_avhrr(indts, bbox_or_stanm_list, region)
  % Store/retrieve data in this M-file's local directory
  [pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = fullfile(pathroot, '../data', '');
  end;


    yrmo_ext = fullfile(datapath, ...
                        sprintf( '%06d-%03dx%03d.mat', yrmo, ...
                                 ((boxradius*2)+1), ((boxradius*2)+1) ));



          pname = fullfile(datapath, [stanm{boxix} '-' fname]);
%           if ( exist(pname, 'file') )
%             dts{boxix}{sstix{boxix}} = datenum(YYYY,MM,DD,hh,mm,0);
%             sstix{boxix} = sstix{boxix} + 1;
%             continue;
%           end;









%function [dts,ssts,urls] = anavhrr(indts, bbox_or_stanm, region, acceptable_bad_pct)
  % Store/retrieve data in this M-file's local directory
  [pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = fullfile(pathroot, '../data', '');
  end;



        pname = fullfile(avhrrpath, [stanm '-' fname]);










      %if ( yr == begyr )
      if ( ~isfield(stn,fld) )
        stn.(fld) = struct('date',[],'data',[]);
      end;


    vars = { 'u', 'v', 'mld', 'temperature', 'salinity', 'temperature',  };
    flds = { 'u', 'v', 'mld', 'seatemp',     'salinity', 'seatemp_field' };



% DEFAULT for BASEURL:
%  'http://tds.hycom.org/opendap/nph-dods/datasets/hycom/GOMl0.04/expt_20.1'
% If you have another ocean model with a THREDDS interface, specify it here.
%
% SAMPLE URL:
%  http://tds.hycom.org/opendap/nph-dods/datasets/hycom/GOMl0.04/expt_20.1/2003/2d/archv.2003_002_00_2d.nc?mld[0:1:0][0:1:384][0:1:540],qtot[0:1:0][0:1:384][0:1:540],ssh[0:1:0][0:1:384][0:1:540]
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1
%
% NetCDF filename and query templates for each variable:
%  2d/archv.2003_002_00_2d.nc?mld[0:1:0][0:1:384][0:1:540],qtot[0:1:0][0:1:384][0:1:540],ssh[0:1:0][0:1:384][0:1:540]
%  salt/archv.2003_002_00_3zs.nc?salinity[0:1:0][0:1:39][0:1:384][0:1:540]
%  temp/archv.2003_002_00_3zt.nc?temperature[0:1:0][0:1:39][0:1:384][0:1:540]
%  uvel/archv.2003_002_00_3zu.nc?u[0:1:0][0:1:39][0:1:384][0:1:540]
%  vvel/archv.2003_002_00_3zv.nc?v[0:1:0][0:1:39][0:1:384][0:1:540]
%  wvel/archv.2003_002_00_3zw.nc?w_velocity[0:1:0][0:1:39][0:1:384][0:1:540]
%
% CALLS: MDATASET (netCDF-Java), QUERY_GOM_HYCOM_INDICES, GET_STATION_COORDS




GET_PATHFINDER.m:
      dat = nc{'bsst'}(latix,lonix);
      close(nc); clear nc;

      % Pixel-to-SST conversion. Pathfinder documentation: "for version 5.0 pixel
      % value with slope of 0.075 and y-intercept of -3.0 to convert to oC".
      sst = (cast(dat,'double') .* 0.075) - 3.0;




    fnames = regexpi(s,'[> ]([12][90][^ ]*[.]hdf)[^a-z]','tokens');


      url = sprintf('%s/%04d/bsst/%04d%03d.s04d1pfrt-bsst.hdf',baseurl,yr,yr,jd);


  wasWarned = false;
      if ( isempty(nc) )
        if ( ~wasWarned )
          warning('HDF files for some days in 1995-2009 were not available!');
          wasWarned = true;
        end;
        continue;
      end;







ANWERA.m:
  % Export figures to a location relative to this M-file's local directory
  [pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  if ( ~exist('figspath', 'var') || isempty(figspath) )
    figspath = fullfile(pathroot, '../figs', '');
  end;





% %%%%  %DEBUG:
%   hafld = [hffld];
% %%%%  %DEBUG:



./(30*24*3600)


%   [stn.monthly_sea_t,clim.sea_t,stn.monthly_sea_t_anom] = ...
%       monthly_clim_ts(stn.ndbc_sea_t);
  [yr,mo,dy] = datevec(stn.ndbc_sea_t.date);
  dts = datenum(yr,mo,1);

  % Monthly mean sea temperature
  stn.monthly_sea_t.date = unique(dts);
  stn.monthly_sea_t.data = grpstats(stn.ndbc_sea_t.data,dts);

  % Monthly sea temperature anomaly
  clim.sea_t = grpstats(stn.ndbc_sea_t.data,mo);
  [yr,mo,dy] = datevec(stn.monthly_sea_t.date);
  moix = grp2idx(mo);
  stn.monthly_sea_t_anom.date = stn.monthly_sea_t.date;
  stn.monthly_sea_t_anom.data = stn.monthly_sea_t.data - clim.sea_t(moix);






    result.wxt_wdir.date = dts;
    result.wxt_wdir.data = grpstats(dat{6},bin);






    %DEBUG:
    return;
    %DEBUG:




function [t,f] = scatter_hc(stn,tfld,ffld,tfix,ffix)
%function [t,f] = scatter_hc(stn,tfld,ffld,tfix,ffix)
%
% Linearly regress the rate of change in STN.(TFLD) (determined using an
% N-point cenetered finite difference), against the values in STN.(FFLD).
% Useful for comparing sea temperature change to heat flux, for example!
%
% Last Saved Time-stamp: <Thu 2010-08-05 15:16:55  Lew.Gramer>




  % badix = find(abs(t.data) > 0.3);
  % t.date(badix) = [];
  % t.data(badix) = [];







  for ix = 1:length(accflds)
    stn = verify_variable(stn,accflds{ix});
    if (~isfield(stn,accflds{ix}))
      warning('No field "%s" could be found or calculated!', accflds{ix});
    else
      legs = { legs{:} accflds{ix} };
      fld = stn.(accflds{ix});
      dtix = find(begdt <= fld.date & fld.date < enddt);
      dts = fld.date(dtix);
      dat = fld.data(dtix);
      goodix = find(isfinite(dat));
      dat(goodix) = cumsum(dat(goodix));
      % plot(dts,real(dat),linspec{mod(cix-1,length(linspec))+1});
%DEBUG:
      [rix,aix] = intersect_dates(firstdts,fld.date);
      val0 = firstdat(rix(1));
      plot(dts,real(dat+val0),linspec{mod(cix-1,length(linspec))+1});
      cix = cix + 1;
    end;
  end;









  % % Include precipitation flux if we got it
  % if ( isfield(stn,rffld) )
  %   stn = station_heat_flux_term( MORE CODE NEEDED HERE );
  %   rffld = [rffld '_term'];
  %   [rfix,nqix] = intersect_dates(stn.(rffld).date,stn.(hcfld).date);
  %   stn.(hcfld).date = stn.(hcfld).date(nqix);
  %   stn.(hcfld).data = stn.(hcfld).data(nqix) + stn.(rffld).data(rfix);
  % end;

  % multiplot_station(stn,{'ndbc_sea_t','bic_ncep_26_heat_flux_term',...
  %                     'bic_ncep_26_net_heat_flux_24_hour_lowpass',...
  %                     'bic_ncep_26_net_heat_flux_24_hour_average',hcfld});
  % xlim([datenum(2009,7,6.5) datenum(2009,7,21.5)]);
  % datetick('x',2,'keeplimits');

  % plot_fluxes(stn,2009,187,202);
  % plot_fluxes(stn,2009,187,365);
  % plot_fluxes(stn,2008,1,365);
  % plot_fluxes(stn,2008,120,280);
  % plot_fluxes(stn,2008,205,365);

%   %LONF1
%   plot_fluxes(stn,2005,90,365*2.5);

%   %FWYF1
%   plot_fluxes(stn,2001,220,365*8.5);

%   %MLRF1
%   plot_fluxes(stn,1997,1,365*12.5);

%   %SMKF1
%   % plot_fluxes(stn,1998,91,365);
%   % plot_fluxes(stn,1999,1,365);
%   % plot_fluxes(stn,2002,191,365*6.5);
%   plot_fluxes(stn,1998,91,365*12);

  yrs = [];
  for yr = 1987:2010
    yrix = find(get_year(stn.(hcfld).date) == yr);
    if ( (length(yrix) < 3000) || (max(diff(stn.(hcfld).date(yrix))) > 14) )
%       warning('Too little data to plot %04d',yr);
    else
      yrs(end+1) = yr;
%       plot_fluxes(stn,yr,1,365);
%       ylim([-15 +15]);
%       appendtitlename(sprintf(' (Abs:%g R:%g c:%g w:%g)',doabs,R,cfac,wfac));
%       figstub = fullfile(figspath,sprintf('%s_tryhc_%04d',stn.station_name,yr));
%       print('-dpng',[figstub '001-365.png']);
% %       plot_fluxes(stn,yr,1,120);
% %       %print('-dpng',[figstub '001-120.png']);
% %       plot_fluxes(stn,yr,120,240);
% %       %print('-dpng',[figstub '120-240.png']);
% %       plot_fluxes(stn,yr,240,365);
% %       %print('-dpng',[figstub '240-365.png']);
    end;
  end;






    % % TOO SLOW! Comment out for now...
    % s([1 3],wk) = bootci(nboot, @nanstd, dat);





x = load(['../data/' stn.station_name '_thesis.mat']);
stnm = lower(stn.station_name);
stn.ncep_precip = x.(stnm).ncep_precip;
stn.ndbc_ncep_30_rain_heat_flux = x.(stnm).ndbc_ncep_30_rain_heat_flux;
stn.ndbc_ncep_30a_rain_heat_flux = x.(stnm).ndbc_ncep_30a_rain_heat_flux;
x = []; clear x; clear stnm;

tic,
  stn = station_bulk_longwave(stn,'ndbc_air_t','ncep_spechumid','ndbc_barom','ncep_cloud_cover','ndbc_sea_t','ncep_cloud_cover','ndbc_ncep_rh_dlrf','ndbc_ncep_rh_ulrf','ndbc_ncep_rh_lrf');
 stn = station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ncep_relhumid','ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_ncep_rh_lrf','ndbc_ncep_rh_30a','ncep_dsrf','ndbc_ncep_rh_dlrf','ncep_precip','ndbc_wind1_dir','default','default','default','default',true);
toc,







  % Mao, Lei and Patterson 2009 (WARMING paper)
  Ra = g.*alpha.*(abs(q0)./(cp.*rho)).*(h.^4)./(nu*(k^2));
  %DEBUG:  disp('MLP09 Ra:'); [nanmin(Ra), nanmean(Ra), nanmax(Ra),],
  x = h ./ b;
  eta = 0.1;

  % Steady-state
  %u = ((x./h).^(1/3)) .* exp(-1./(3.*b.*x.*eta)) .* (Ra.^(1/3)) .* (k./h);













    [yr, length(yrix), max(diff(stn.netqf.date(yrix))), ],







%         'ndbc_air_t', ...
%         'ndbc_sea_t_1_day_decimate', ...
%         'ndbc_sea_t_1_day_lowpass', ...
%         'ndbc_sea_t_7_day_lowpass', ...
%         'ndbc_wind1_speed', ...
%         'ndbc_tide', ...
%         'ndbc_air_sea_t', ...


%   accflds = { ...
%       'ncep_heat_flux_term', ...
%       'ndbc_hfbulk_heat_flux_term',...
%       'ndbc_hfbulk_dt',...
%       'ndbc_bulk_30a_heat_flux_term', ...
%       'ndbc_bulk_30a_dt', ...
%       'ndbc_ncep_26_heat_flux_term', ...
%       'ndbc_ncep_26_dt', ...
%       'ndbc_ncep_30_heat_flux_term', ...
%       'ndbc_ncep_30_dt', ...
%       'ndbc_ncep_30a_heat_flux_term', ...
%       'ndbc_ncep_30a_dt', ...
%       'ndbc_ncep_30a_heat_flux_term_1_day_lowpass', ...
%       'ndbc_ncep_30a_dt_1_day_lowpass', ...
%       'ndbc_ncep_30a_heat_flux_term_7_day_lowpass', ...
%       'ndbc_ncep_30a_dt_7_day_lowpass', ...
%             };

%         'ndbc_hfbulk_heat_flux_term',...
%         'ndbc_hfbulk_dt',...


%         'ndbc_ncep_30a_total_heat_flux_term', ...

%   accflds = { ...
%         'bic_ncep_heat_flux_term', ...
%         'bic_ncep_dt', ...
%         'bic_ncep_rh_heat_flux_term', ...
%         'bic_ncep_rh_dt', ...
%         'bic_ncep_26_heat_flux_term', ...
%         'bic_ncep_26_dt', ...
%             };

%   accflds = { ...
%         'ndbc_ncep_30a_latent_flux_term', ...
%         'ndbc_ncep_30a_sensible_flux_term', ...
%         'ndbc_ncep_30a_shortwave_flux_term', ...
%         'ndbc_ncep_30a_longwave_flux_term', ...
%         'ndbc_ncep_30a_heat_flux_term', ...
%         'ndbc_ncep_30a_dt', ...
%             };


%   accflds = { ...
%       'ndbc_hfbulk_heat_flux_term',...
%       'ndbc_ncep_26_heat_flux_term', ...
%       'ndbc_ncep_30_heat_flux_term', ...
%       'ndbc_ncep_30a_heat_flux_term', ...
%             };

%   accflds = { ...
%       'ndbc_ncep_30a_heat_flux_term', ...
%       'ndbc_ncep_30a_dt', ...
%             };

%       'ndbc_ncep_30a_heat_flux_term_1_day_decimate', ...
%       'ndbc_ncep_30a_dt_1_day_decimate', ...

%       'ndbc_ncep_26_heat_flux_term_1_day_decimate', ...
%       'ndbc_ncep_26_dt_1_day_decimate', ...
%       'ndbc_ncep_30_heat_flux_term_1_day_decimate', ...
%       'ndbc_ncep_30_dt_1_day_decimate', ...
%       'ndbc_ncep_30a_heat_flux_term_1_day_decimate', ...
%       'ndbc_ncep_30a_dt_1_day_decimate', ...


%       'ndbc_bulk_26_dt', ...
%       'ndbc_bulk_30_dt', ...
%       'ndbc_bulk_30a_dt', ...

%       'ndbc_ncep_26_dt', ...
%       'ndbc_ncep_30_dt', ...
%       'ndbc_ncep_30a_dt', ...


%       'ndbc_hfbulk_heat_flux_term', ...
%       'ndbc_bulk_30a_heat_flux_term', ...
%       'ndbc_ncep_30a_heat_flux_term', ...

%       'ndbc_bulk_26_heat_flux_term', ...
%       'ndbc_bulk_30_heat_flux_term', ...
%       'ndbc_bulk_30a_heat_flux_term', ...

%       'ndbc_ncep_26_heat_flux_term', ...
%       'ndbc_ncep_30_heat_flux_term', ...
%       'ndbc_ncep_30a_heat_flux_term', ...

  else
    % Diagnosing
    relflds = { ...
        'ndbc_sea_t', ...
        'ndbc_tide', ...
              };
%     accflds = { ...
%         'ncep_swhf_term', ...
% %         'ncep_srf_term', ...
%         'ncep_lwhf_term', ...
% %         'ncep_lrf_term', ...
%         'ndbc_ncep_26_latent_flux_term', ...
%         'ndbc_ncep_26_sensible_flux_term', ...
%         'ndbc_ncep_26_heat_flux_term', ...
%               };

    absflds = { ...
        'ndbc_ncep_30a_wind_stress', ...
              };

    accflds = { ...
        'ncep_swhf_term', ...
        'ncep_srf_term', ...
        'ncep_lwhf_term', ...
        'ncep_lrf_term', ...
        'ndbc_ncep_30a_latent_flux_term', ...
        'ndbc_ncep_30a_sensible_flux_term', ...
        'ndbc_ncep_30a_heat_flux_term', ...
              };

  end;








  if ( nargin < 2 )
    yr=2005;
  end;
  if ( nargin < 3 )
    begdy = 1;
  end;
  if ( nargin < 4 )
    enddy = begdy+3;
  end;
























  stn = station_tmd_tide(stn,stn.ndbc_sea_t.date);

  hfld = 'tide_i_depth';
  if ( ~isfield(stn,hfld) )
    hfld = 'tmd_tide_i_depth';
  end;
  %DEBUG:
  hfld = 'tmd_tide_i_depth';
  %DEBUG:  hfld = 'tpx_tide_i_depth';
  %DEBUG:
  disp(hfld);











%   % Unbalanced thermal forcing, stress divergence
%   Qv(warmix) = 1.00 .* sqrt(b .* (uf(warmix).^2) .* (24*3600) ./ h(warmix));		% OK for warming, tarrible for cooling





%   coolix = find(q0 > +50);
%   warmix = find(q0 < -50);
%   Qv = repmat(0,size(uf));
%   % Unbalanced thermal forcing, viscous/unsteady inertia
%   Qv(coolix) = 1.00 .* b .* (uf(coolix).^3) .* ((24*3600)^2) ./ (h(coolix).^2);		% VERY GOOD for cooling
%   % Unbalanced thermal forcing, stress divergence
%   Qv(warmix) = 1.00 .* sqrt(b .* (uf(warmix).^2) .* (24*3600) ./ h(warmix));		% Good for warming










% Tried to use l=b/h: No, result was *TERRIBLE*!
  u = Ra .* (((b.^2)./h).^4) .* (x.^3) .* k;
  %DEBUG:  disp('MLP10 u:'); [nanmin(u), nanmean(u), nanmax(u),],





%function stn = station_horizontal_convection(stn,tf,sf,hf,q0f,qvf,qf,R)

%   q24f = [q0f '_24_hour_average'];
%   stn = verify_variable(stn,q24f);
%   [q24fix,qfix] = intersect_dates(stn.(q24f).date,stn.(qf).date);
%   noopix = find( abs( stn.(q24f).data ) < 50 );

%   q40f = [qf '_40_hour_average'];
%   stn = verify_variable(stn, q40f);
%   coolix = find(diff(sign(stn.(q40).data)) < 0);
%   warmix = find(diff(sign(stn.(q40).data)) > 0);











%   % Unbalanced thermal forcing, viscous/unsteady inertia
%   Qv = 1.00 .* b .* (uf.^3) .* ((24*3600)^2) ./ (h.^2);		% GOOD for cooling

%   warmix = find(sign(q0) < 0);
%   % Unbalanced thermal forcing, stress divergence
%   Qv(warmix) = 1.00 .* sqrt(b .* (uf(warmix).^2) .* (24*3600) ./ h(warmix));		% GOOD for warmin!











function [T,Dt,B,Qv,Q] = thermal_exchange(t,s,h,q0,b,R)
%function [T,Dt,B,Qv,Q] = thermal_exchange(t,s,h,q0,b,R)
%
% Thermally-induced exchange flow (thermal syphon, horizontal advection).
%
% From physical data (any combination of same-sized vectors and scalars):
%   t = sea temperature [oC]
%   s = salinity [psu]
%   h = site depth [m] or pressure [db]
%   q0 = net surface heat flux [W m^-2]
%   b = bottom slope (rise/run)
%   R = mixing(Ri) + entrainment rate (DEFAULT: 0.30)
%
% ... use formulae in literature to calculate thermal exchange parameters:
%   T = onset time for steady-state thermal exchange flow [s]
%   Dt = critical (maximum) depth for viscosity-dominated exchange [m]
%   B = net surface buoyance flux [m^2 s^-3]
%   Qv = volumetric discharge rate [m^2 s^-1]
%   Q = heat exchange [K/hr]
%
% SEE: Farrow and Patterson 1993, Sturman et al. 1999, Monismith et al. 2006,
%  Hughes and Griffiths 2008, Chubarenko 2010, Mao et al. 2010a and 2010b
%
% Last Saved Time-stamp: <Wed 2010-06-23 16:20:25 Eastern Daylight Time gramer>

  if ( ~exist('R','var') || isempty(R) )
    R = 0.30;
  end;

  g = 9.79;
  alpha = sw_alpha(s,t,h);
  rho = sw_dens(s,t,h);
  cp = sw_cp(s,t,h);

  % Net surface buoyancy flux (ignoring runoff, for now...)
  B = (g .* alpha .* abs(q0)) ./ (rho .* cp);

  % Characteristic convective velocity scale
  uf = (B .* h) .^ (1/3);

  % Critical depth for viscous balance (assuming diurnal forcing)
  Dt = 0.3 .* uf .* (24*3600);


  % Volumetric discharge rate Qv [m2/s]

  % Sturman et al. 1999
  % l = 1000; %[m]
  % Qv = ((l .* B).^(1/3)) .* h;
  % Qv = 0.24 .* (B.^(1/3)) .* ( (l.*b./(1+b)) .^ (4/3) );
  % % h = l*b
  % Qv = 0.24 .* (B.^(1/3)) .* ( (h./(1+b)) .^ (4/3) );

  % Monismith et al. 2006
  % % Qv = 0.36 .* (1./(1+b)) .* uf .* h;
  % Qv = 0.3 .* uf .* h;
%   % Steady thermal balance, advective inertia
%   Qv = 0.15 .* (b.^(-1/3)) .* uf .* h;			% Bad

%   % Unbalanced thermal forcing, advective inertia
%   Qv = 1.00 .* sqrt((uf.^3) .* (24*3600) ./ h);		% Not great

%   % Balanced thermal forcing, viscous/unsteady inertia
%   Qv = 1.00 .* sqrt((uf.^3) .* (24*3600) ./ h);		% Also not great

%   % Balanced thermal forcing, stress divergence
%   Qv = 1.00 .* uf;						% The worst

%   % Unbalanced thermal forcing, viscous/unsteady inertia
%   Qv = 1.00 .* b .* (uf.^3) .* ((24*3600)^2) ./ (h.^2);	% GOOD for cooling

%   % Unbalanced thermal forcing, stress divergence
%   Qv = 1.00 .* sqrt(b .* (uf.^2) .* (24*3600) ./ h);		% GOOD for heating!

  % Heat exchange rate Q [K/hr]
  l = 3600 .* Qv ./ h;
%%%% ??? DEBUG
%l = 1e3;
%%%% ??? DEBUG
  D = h + (b.*l);
  % dT = 3600 .* ( (q0./(cp.*rho.*h)) - (q0./(cp.*rho.*(h+D))) );
  dT = 3600 .* ( (q0./(cp.*rho.*(h+D))) - (q0./(cp.*rho.*h)) );
  Q = R .* dT;


  % Onset time T for steady convection [s]

%   % Chubarenko 2010
%   % reducedg = dT .* alpha,
%   rho0 = sw_dens(s,(t+dT),h);
%   reducedg = (rho - rho0) ./ rho0;
%   T = sqrt(l./(g.*reducedg.*b));

  % Lei and Patterson 2005
  nu = 1.05e-6;		% Molecular kinematic viscosity of seawater at 35psu, 20oC [m2/s]
  % nu = 1e-3;		% Eddy kinematic viscosity estimated over a reef [m2/s]
  k = 1.46e-7;		% Thermal diffusivity of seawater at 35psu, 20oC [m2/s]
  % Rac = 657.5;		% Critical Rayleigh number
  % Ra = g.*alpha.*(q0./(cp.*rho.*h)).*(h.^4)./(nu*(k^2));
  % T = sqrt(Rac/Ra).*(h^.2)./k;

  % Mao, Lei and Patterson 2010
  Ra = g.*alpha.*(abs(q0)./(cp.*rho.*h)).*(l.^4)./(nu*(k^2));
  x = h / b;
  T = (x.^(2/3)) .* (Ra.^(-1/3)) .* (l.^(4/3)) ./ k;
  %DEBUG:
  [nanmin(T), nanmean(T), nanmax(T),]

  lc = (Ra.^(-1/4)) .* (b.^(-3/2)) .* l;
  %DEBUG:  [nanmin(lc), nanmean(lc), nanmax(lc),]
  u = Ra .* (b.^4) .* (l.^-4) .* (x.^3) .* k;
  %DEBUG:  [nanmin(u), nanmean(u), nanmax(u),]

  QMLP = Ra .* (b.^5) .* ((x./l).^4) .* k;
  %DEBUG:  [nanmin(QMLP), nanmean(QMLP), nanmax(QMLP),]

return;














  % Heat exchange rate [K/day]
  l = 24 .* 3600 .* Qv ./ p;
  D = p + (b.*l);
  dT = 3600 .* ( (q0./(cp.*rho.*p)) - (q0./(cp.*rho.*(p+D))) ),
  Q = R .* dT;

  % Chubarenko 2010
  % % T = (1)*3600*24;
  % reducedg = dT .* alpha,
  rho0 = sw_dens(s,(t+dT),p);
  reducedg = (rho - rho0) ./ rho0;
  T = sqrt(l./(g.*reducedg.*b));








  badix = find(theta < 2);
  z(badix) = [];
  q0(badix) = [];
  dts(badix) = [];






      for termfld = {'latent','sensible','shortwave','longwave'}







%%%%%%%%%% ??? DEBUG
          idepthfld = [];
%%%%%%%%%% ??? DEBUG




%       'ndbc_ncep_30_heat_flux_term', ...
%       'ndbc_ncep_30_dt', ...




  [wz,az,pz,stz] = station_instrument_heights(station.station_name);
  Q0_factor = Q0factor(station.ndbc_sea_t.data,[],stz);




      station.(htfld).date = station.(nffld).date;
      station.(htfld).data = (station.(nffld).data ./ Q0_factor) .* (60*60);









  [wz,az,pz,stz] = station_instrument_heights(station.station_name);
  Q0_factor = Q0factor(station.ndbc_sea_t,[],stz);

  for q0fld = { ...
      'ndbc_hfbulk',...
      'ndbc_bulk_rh','ndbc_bulk_26','ndbc_bulk_30','ndbc_bulk_30a',...
      'ndbc_ncep_26','ndbc_ncep_30','ndbc_ncep_30a',...
              }
    nffld = [q0fld{:} '_net_heat_flux'];
    htfld = [q0fld{:} '_heat_flux_term'];
    dtfld = [q0fld{:} '_dt'];

    if ( isfield(station,nffld) && ~isfield(station,htfld) )
      station.(htfld).date = station.(nffld).date;
      station.(htfld).data = (station.(nffld).data ./ Q0_factor) .* (60*60);
    end;

    if ( isfield(station,htfld) )
      %DEBUG:      disp(htfld);
      [ix1,ix2] = intersect_dates(station.(htfld).date, ...
                                  station.(UdotdelT).date);
      station.(dtfld).date = station.(htfld).date(ix1);
      station.(dtfld).data = station.(UdotdelT).data(ix2) + ...
          station.(htfld).data(ix1);
    else
      %DEBUG:
      disp(['No field ' htfld]);
    end;
  end;










  if ( isfield(station,'ndbc_hfbulk_net_heat_flux') && ...
       ~isfield(station,'ndbc_hfbulk_heat_flux_term') )
    [wz,az,pz,stz] = station_instrument_heights(station.station_name);
    Q0_factor = Q0factor(station.ndbc_sea_t,[],stz);
    station.ndbc_hfbulk_heat_flux_term.date = ...
        station.ndbc_hfbulk_net_heat_flux.date;
    station.ndbc_hfbulk_heat_flux_term.data = ...
        (station.ndbc_hfbulk_net_heat_flux.data ./ Q0_factor) .* (60*60);
  end;










  % HACK FOR SMKF1 ONLY!
  [station.lon,station.lat,station.depth] = get_station_coords(station.station_name);






  stokesu = rawstokesu;
  stokesv = rawstokesv;






  hycomfname = fullfile(datapath,[stnm '_hycom.mat']);
  if ( exist(hycomfname,'file') )
    disp(['Loading from presaved ' hycomfname]);
    load(hycomfname);
  else
    warning('Substituting all-zero surface currents for missing "%s"!', hycomfname);
    stn.global_hycom_u.date = ...
        ceil(station.(rawstokesu).date(1)) : 1 : ...
        floor(station.(rawstokesu).date(end));
    stn.global_hycom_u.data = repmat(0,size(stn.global_hycom_u.date));
    stn.global_hycom_v = stn.global_hycom_u;
  end;

  station.global_hycom_u = stn.global_hycom_u;
  station.global_hycom_v = stn.global_hycom_v;
  stn = []; clear stn;







      [PFX 'stokes_drift_u_7_day_lowpass'],...
      [PFX 'stokes_drift_v_7_day_lowpass'],...
      'global_hycom_u_7_day_lowpass','global_hycom_v_7_day_lowpass',...








  % station = verify_variable(station, [PFX 'stokes_drift_u_7_day_lowpass']);
  % station = verify_variable(station, [PFX 'stokes_drift_v_7_day_lowpass']);
  % station = verify_variable(station, 'global_hycom_u_7_day_lowpass');
  % station = verify_variable(station, 'global_hycom_v_7_day_lowpass');








  if ( ~ischar(bbox_or_stanm) )
    bbox = bbox_or_stanm;
    miny = (bbox(1) - minlon) / dlon;
    maxy = (bbox(2) - minlon) / dlon;
    minx = (maxlat - bbox(3)) / dlat;
    maxx = (maxlat - bbox(4)) / dlat;
    stn_x = minx + round((maxx+minx)/2);
    stn_y = miny + round((maxy+miny)/2);

  else
    stanm = bbox_or_stanm;

    [stn_lon, stn_lat, ig] = get_station_coords(stanm);

    stn_x = round((maxlat - stn_lat) / dlat);
    stn_y = round((stn_lon - minlon) / dlon);

    minx = stn_x - boxradius; maxx = stn_x + boxradius;
    miny = stn_y - boxradius; maxy = stn_y + boxradius;
  end;








  meansst = squeeze(nanmean(sst,1));
  [landrow,landcol] = find(~isfinite(meansst));
  badix = find(~isfinite(sst));
  [badwk,badrow,badcol] = ind2sub(size(sst),badix);
  badwk(ismember(badrow,landrow) & ismember(badcol,landcol)) = [];
  badwk = unique(badwk);
  dts(badwk) = [];
  sst(badwk,:,:) = [];

  alldts = dts(1):dts(end);
  station.avhrr_weekly_dates = dts;
  station.avhrr_weekly_sst.date = alldts;
  station.avhrr_weekly_sst.field = interp1(dts,sst,alldts,'pchip');

  sst = []; clear sst;











  badweeks = find(~isfinite(sst(:,1,1)));
  dts(badweeks) = [];
  sst(badweeks,:,:) = [];

  alldts = dts(1):dts(end);
  station.avhrr_weekly_dates = dts;
  station.avhrr_weekly_sst.date = alldts;
  station.avhrr_weekly_sst.field = interp1(dts,sst,alldts,'pchip');

  sst = []; clear sst;













    if ( ~isfield(stn, newfld) )
      stn.(newfld).date = [];
      stn.(newfld).data = [];
    end;





% Just combine all fields as loaded from each year-month MAT file, into one
% master struct WW3 for this station - we will merge with STN later
function ww3 = get_ww3_station_mash_all_fields(ww3,dat)
  if ( isempty(ww3) )
    ww3 = dat;
  else
    flds = fieldnames(dat);
    for fldix = 1:length(flds)
      fld = flds{fldix};
      if ( ~isfield(ww3, fld) )
        ww3.(fld).date = [];
        ww3.(fld).data = [];
      end;
      if ( isfield(dat.(fld),'data') )
        ndat = numel(dat.(fld).data);
        ww3.(fld).date(end+1:end+ndat,1) = dat.(fld).date;
        ww3.(fld).data(end+1:end+ndat,1) = dat.(fld).data;
      end;
    end;
  end;
return;


% Merge data from master WW3 struct with station struct STN
function stn = get_ww3_station_merge_fields(stn,ww3,dataset)

  flds = fieldnames(ww3);
  flds = flds(strmatch('ww3_',flds));
  newflds = strrep(flds,['ww3_' dataset '_' ],'ww3_');
  for fldix = 1:length(flds)
    fld = flds{fldix};
    newfld = newflds{fldix};
    if ( ~isfield(stn, newfld) )
      stn.(newfld).date = [];
      stn.(newfld).data = [];
    end;
    dts = ww3.(fld).date(1):(1/24):ww3.(fld).date(end);
    ndat = numel(dts);
    stn.(newfld).date(end+1:end+ndat,1) = dts;
    goodix = find(isfinite(ww3.(fld).data));
    stn.(newfld).data(end+1:end+ndat,1) = ...
        interp1(ww3.(fld).date(goodix),ww3.(fld).data(goodix),dts);
  end;

return;










function extract_ww3

  datapath = fullfile('..','data');

  dataset = 'wna';

  stncfg = {};
  cfgfname = fullfile(datapath,['ww3-' dataset '.cfg']);
  if ( exist(cfgfname,'file') )
    fid = fopen(cfgfname,'r');
    stncfg = textscan(fid,'%s%d%d', 'Delimiter',',', 'CommentStyle','#'); 
    fclose(fid);
  end;

  %ftp://polar.ncep.noaa.gov/pub/history/waves/wna.hs.200711.grb

  % tp - Peak Wave Period
  % dp - Peak Wave Direction
  % hs - Significant Wave Height
  vars = { 'tp',         'dp',         'hs', };
  flds = { 'peakwaveper','peakwavedir','sigwavehgt', };
  %DEBUG:  vars = { 'hs', };  flds = { 'sigwavehgt', };

  f = [];
  doFTP = false;
  if ( doFTP )
    fhost = 'polar.ncep.noaa.gov';
    fuser = 'anonymous';
    fpawd = 'lew.gramer@noaa.gov';
    fdir = '/pub/history/waves';

    f = ftp(fhost,fuser,fpawd);
    binary(f);
    cd(f, fdir);
  end;


  region = 'fknms';
  minlat = 24;
  maxlat = 27;
  minlon = 360-84;
  maxlon = 360-79;

  %for yr = 1999:2007
  for yr = 1999

    switch ( yr ),
      % case 1999, mos = 7:12;
     case 1999, mos = 7;
     case 2007, mos = 1:11;
     otherwise, mos = 1:12;
    end;

    for mo = mos(:)'

      tic,

      dat = [];
      clear dat;
      dat = [];

      for vix = 1:length(vars)

        var = vars{vix};

        fname = fullfile(datapath, ...
                         sprintf('%s.%s.%04d%02d.grb', dataset, var, yr, mo));
        %DEBUG:
        disp(fname);

        if ( ~isempty(f) )
          if ( ~exist(fname,'file') )
            mget(f, fname);
          end;
        end;

        x = read_grib(fname,-1, 'ScreenDiag',0);
        %DEBUG:            size(x),

        if ( ~exist('minlatix','var') )
          nlon = x(1).gds.Ni;
          nlat = x(1).gds.Nj;

          dlon = x(1).gds.Di;
          %dlat = x(1).gds.Dj;
          % Using dlat=gds.Di because some of these stupid hacker-
          % produced "GRIB" files apparently HAVE NO gds.Dj field!
          dlat = x(1).gds.Di;

          lon1 = x(1).gds.Lo1; lon2 = x(1).gds.Lo2;
          lat1 = x(1).gds.La1; lat2 = x(1).gds.La2;

          lons = lon1 :  dlon : lon2;
          lons(lons > 180) = lons(lons > 180) - 360;
          % Latitude indices are reversed (naturally!)
          lats = lat1 : -dlat : lat2;

          minlon(minlon < 0) = minlon(minlon < 0) + 360;
          maxlon(maxlon < 0) = maxlon(maxlon < 0) + 360;
          minlonix = round( (minlon - lon1) / dlon ) + 1;
          maxlonix = round( (maxlon - lon1) / dlon ) + 1;

          minlatix = round( (lat1 - maxlat) / dlat ) + 1;
          maxlatix = round( (lat1 - minlat) / dlat ) + 1;
        end;

        if ( ~isfield(dat,'n') )
          dat.n = length(x);
          dat.lon = lons(minlonix:maxlonix);
          dat.lat = lats(minlatix:maxlatix);
          dat.nlon = length(dat.lon);
          dat.nlat = length(dat.lat);
          for ix = 1:length(x)
            dat.date(ix) = datenum(x(ix).stime);
          end;
        end;

        for ix = 1:length(x)
          rawdat = reshape(x(ix).fltarray, [nlon nlat])';
          dat.(var)(ix,1:dat.nlat,1:dat.nlon) = rawdat(minlatix:maxlatix,minlonix:maxlonix);
          rawdat = []; clear rawdat;
        end;
        dat.(var)(dat.(var) > 360) = nan;

        x = [];
        clear x;

        % if ( ~isempty(f) )
        %     delete(fname);
        % end;

      end; %for vix

      %DEBUG:
      id=80; figure; contourf(squeeze(dat.hs(id,:,:))); maxigraph; caxis([0 3]); colorbar; hold on; set(gca,'yd','rev'); title(['sigwavehgt ' datestr(dat.date(id))]);

      % Subset our world list of stations to those inside our BBOX
      if ( ~exist('stns','var') || isempty(stns) )
        stns = get_all_station_metadata;
        XV = [ min(dat.lon) max(dat.lon) max(dat.lon) min(dat.lon) ];
        YV = [ max(dat.lat) max(dat.lat) min(dat.lat) min(dat.lat) ];
        goodix = find( inside(stns.lons, stns.lats, XV, YV) );
        clear XV YV;
        stns.codes  = stns.codes(goodix);
        stns.lons   = stns.lons(goodix);
        stns.lats   = stns.lats(goodix);
        stns.depths = stns.depths(goodix);
      end;

      for stix = 1:length(stns.lons)

        stnm = lower(stns.codes{stix});
        dat.(stnm).lon = stns.lons(stix);
        dat.(stnm).lat = stns.lats(stix);
        ididx = [];
        if ( ~isempty(stncfg) )
          ididx = find(strcmpi(stnm,stncfg{1}));
        end;
        if ( ~isempty(ididx) )
          lonix = double(stncfg{2}(ididx));
          latix = double(stncfg{3}(ididx));
        else
          [lonix,latix] = gridnbhd_km(dat.lon, dat.lat,...
                                      dat.(stnm).lon, dat.(stnm).lat, 0);
        end;
        dat.(stnm).lonix = lonix;
        dat.(stnm).latix = latix;
        %DEBUG:
        fprintf('%s,%d,%d\n',stnm,dat.(stnm).lonix,dat.(stnm).latix);

        for vix = 1:length(vars)
          var = vars{vix};
          fld = ['ww3_' dataset '_' flds{vix}];
          dat.(stnm).(fld).date = dat.date(:);
          dat.(stnm).(fld).data = squeeze(dat.(var)(:,dat.(stnm).latix,dat.(stnm).lonix));
        end; %for vix

        %DEBUG:
        plot(dat.(stnm).lonix,dat.(stnm).latix,'k*'); text(dat.(stnm).lonix,dat.(stnm).latix,stnm);

      end; %for stix

      %DEBUG: 
      figure; plot(dat.date,[dat.fwyf1.ww3_wna_sigwavehgt.data, dat.mlrf1.ww3_wna_sigwavehgt.data, dat.smkf1.ww3_wna_sigwavehgt.data, ]); maxigraph; datetick3; legend('FWYF1','MLRF1','SMKF1'); title('sigwavehgt');

      matfname = fullfile(datapath, ...
                          sprintf('ww3.%s.%04d%02d.mat', dataset, yr, mo));
      disp(['Saving to ' matfname]);
      save(matfname,'dat');  
      toc,

    end; %for mo

  end; %for yr

  if ( ~isempty(f) )
    close(f);
  end;

return;
















      dy=12; figure; contourf(squeeze(dat.hs(dy,:,:))); maxigraph; caxis([0 3]); colorbar; hold on; set(gca,'yd','rev'); title(sprintf('sigwavehgt day %d',dy));
      dy=14; figure; contourf(squeeze(dat.hs(dy,:,:))); maxigraph; caxis([0 3]); colorbar; hold on; set(gca,'yd','rev'); title(sprintf('sigwavehgt day %d',dy));
      dy=16; figure; contourf(squeeze(dat.hs(dy,:,:))); maxigraph; caxis([0 3]); colorbar; hold on; set(gca,'yd','rev'); title(sprintf('sigwavehgt day %d',dy));







function extract_ww3

  datapath = fullfile('..','data');

  dataset = 'wna';

  stncfg = {};
  cfgfname = fullfile(datapath,['ww3-' dataset '.cfg']);
  if ( exist(cfgfname,'file') )
    fid = fopen(cfgfname,'r');
    stncfg = textscan(fid,'%s,%d,%d','CommentStyle','%');
    fclose(fid);
  end;

  %ftp://polar.ncep.noaa.gov/pub/history/waves/wna.hs.200711.grb

  % tp - Peak Wave Period
  % dp - Peak Wave Direction
  % hs - Significant Wave Height
  vars = { 'tp',         'dp',         'hs', };
  flds = { 'peakwaveper','peakwavedir','sigwavehgt', };
  %DEBUG:  vars = { 'hs', };  flds = { 'sigwavehgt', };

  f = [];
  doFTP = false;
  if ( doFTP )
    fhost = 'polar.ncep.noaa.gov';
    fuser = 'anonymous';
    fpawd = 'lew.gramer@noaa.gov';
    fdir = '/pub/history/waves';

    f = ftp(fhost,fuser,fpawd);
    binary(f);
    cd(f, fdir);
  end;


  region = 'fknms';
  minlat = 24;
  maxlat = 27;
  minlon = 360-84;
  maxlon = 360-79;

  %for yr = 1999:2007
  for yr = 1999

    switch ( yr ),
      % case 1999, mos = 7:12;
     case 1999, mos = 7;
     case 2007, mos = 1:11;
     otherwise, mos = 1:12;
    end;

    for mo = mos(:)'

      tic,

      dat = [];
      clear dat;
      dat = [];

      for vix = 1:length(vars)

        var = vars{vix};

        fname = fullfile(datapath, ...
                         sprintf('%s.%s.%04d%02d.grb', dataset, var, yr, mo));
        %DEBUG:
        disp(fname);

        if ( ~isempty(f) )
          if ( ~exist(fname,'file') )
            mget(f, fname);
          end;
        end;

        x = read_grib(fname,-1, 'ScreenDiag',0);
        %DEBUG:            size(x),

        if ( ~exist('minlatix','var') )
          nlon = x(1).gds.Ni;
          nlat = x(1).gds.Nj;

          dlon = x(1).gds.Di;
          %dlat = x(1).gds.Dj;
          % Using dlat=gds.Di because some of these stupid hacker-
          % produced "GRIB" files apparently HAVE NO gds.Dj field!
          dlat = x(1).gds.Di;

          lon1 = x(1).gds.Lo1; lon2 = x(1).gds.Lo2;
          lat1 = x(1).gds.La1; lat2 = x(1).gds.La2;

          lons = lon1 :  dlon : lon2;
          lons(lons > 180) = lons(lons > 180) - 360;
          % Latitude indices are reversed (naturally!)
          lats = lat1 : -dlat : lat2;

          minlon(minlon < 0) = minlon(minlon < 0) + 360;
          maxlon(maxlon < 0) = maxlon(maxlon < 0) + 360;
          minlonix = round( (minlon - lon1) / dlon ) + 1;
          maxlonix = round( (maxlon - lon1) / dlon ) + 1;

          minlatix = round( (lat1 - maxlat) / dlat ) + 1;
          maxlatix = round( (lat1 - minlat) / dlat ) + 1;
        end;

        if ( ~isfield(dat,'n') )
          dat.n = length(x);
          dat.lon = lons(minlonix:maxlonix);
          dat.lat = lats(minlatix:maxlatix);
          dat.nlon = length(dat.lon);
          dat.nlat = length(dat.lat);
          for ix = 1:length(x)
            dat.date(ix) = datenum(x(ix).stime);
          end;
        end;

        for ix = 1:length(x)
          rawdat = reshape(x(ix).fltarray, [nlon nlat])';
          dat.(var)(ix,1:dat.nlat,1:dat.nlon) = rawdat(minlatix:maxlatix,minlonix:maxlonix);
          rawdat = []; clear rawdat;
        end;
        dat.(var)(dat.(var) > 360) = nan;

        x = [];
        clear x;

        % if ( ~isempty(f) )
        %     delete(fname);
        % end;

      end; %for vix

      %DEBUG:
      figure; contourf(squeeze(dat.hs(15,:,:))); maxigraph; caxis([0 4]); colorbar; hold on; set(gca,'yd','rev');

      % Subset our world list of stations to those inside our BBOX
      if ( ~exist('stns','var') || isempty(stns) )
        stns = get_all_station_metadata;
        XV = [ min(dat.lon) max(dat.lon) max(dat.lon) min(dat.lon) ];
        YV = [ max(dat.lat) max(dat.lat) min(dat.lat) min(dat.lat) ];
        goodix = find( inside(stns.lons, stns.lats, XV, YV) );
        clear XV YV;
        stns.codes  = stns.codes(goodix);
        stns.lons   = stns.lons(goodix);
        stns.lats   = stns.lats(goodix);
        stns.depths = stns.depths(goodix);
      end;

      for stix = 1:length(stns.lons)

        stnm = lower(stns.codes{stix});
        dat.(stnm).lon = stns.lons(stix);
        dat.(stnm).lat = stns.lats(stix);
        ididx = [];
        if ( ~isempty(stncfg) )
          ididx = find(strcmpi(stnm,stncfg{1}));
        end;
        if ( ~isempty(ididx) )
          lonix = double(stncfg{2}(ididx));
          latix = double(stncfg{3}(ididx));
        else
          [lonix,latix] = gridnbhd_km(dat.lon, dat.lat,...
                                      dat.(stnm).lon, dat.(stnm).lat, 0);
        end;
        dat.(stnm).lonix = lonix;
        dat.(stnm).latix = latix;
        %DEBUG:
        fprintf('%s,%d,%d\n',stnm,dat.(stnm).lonix,dat.(stnm).latix);

        for vix = 1:length(vars)
          var = vars{vix};
          fld = ['ww3_' dataset '_' flds{vix}];
          dat.(stnm).(fld).date = dat.date(:);
          dat.(stnm).(fld).data = squeeze(dat.(var)(:,dat.(stnm).latix,dat.(stnm).lonix));
        end; %for vix

        %DEBUG:
        plot(dat.(stnm).lonix,dat.(stnm).latix,'k*'); text(dat.(stnm).lonix,dat.(stnm).latix,stnm);

      end; %for stix

      %DEBUG: 
      figure; plot(dat.date,[dat.fwyf1.ww3_wna_sigwavehgt.data, dat.mlrf1.ww3_wna_sigwavehgt.data, dat.smkf1.ww3_wna_sigwavehgt.data, ]); maxigraph; datetick3; legend('FWYF1','MLRF1','SMKF1'); title('sigwavehgt');

      matfname = fullfile(datapath, ...
                          sprintf('ww3.%s.%04d%02d.mat', dataset, yr, mo));
      disp(['Saving to ' matfname]);
      save(matfname,'dat');  
      toc,

    end; %for mo

  end; %for yr

  if ( ~isempty(f) )
    close(f);
  end;

return;






















1;

minlat = 24;
maxlat = 27;
minlon = 360-84;
maxlon = 360-79;

%ftp://polar.ncep.noaa.gov/pub/history/waves/wna.hs.200711.grb

dataset = 'wna';

% tp - Peak Wave Period
% dp - Peak Wave Direction
% hs - Significant Wave Height
vars = { 'tp',         'dp',         'hs', };
flds = { 'peakwaveper','peakwavedir','sigwavehgt', };
% %%%% ??? DEBUG
% vars = { 'hs', };
% flds = { 'sigwavehgt', };

f = [];
doFTP = false;
if ( doFTP )
    fhost = 'polar.ncep.noaa.gov';
    fuser = 'anonymous';
    fpawd = 'lew.gramer@noaa.gov';
    fdir = '/pub/history/waves';

    f = ftp(fhost,fuser,fpawd);
    binary(f);
    cd(f, fdir);
end;

%for yr = 1999:2007
for yr = 1999

  switch ( yr ),
    % case 1999, mos = 7:12;
   case 1999, mos = 7;
   case 2007, mos = 1:11;
   otherwise, mos = 1:12;
  end;

  for mo = mos(:)'

    tic,

      for vix = 1:length(vars)

        var = vars{vix};

        fname = sprintf('%s.%s.%04d%02d.grb', dataset, var, yr, mo);
        %DEBUG:
        disp(fname);

        if ( ~isempty(f) )
          if ( ~exist(fname,'file') )
            mget(f, fname);
          end;
        end;

        x = read_grib(fname,-1, 'ScreenDiag',0);
        %DEBUG:            size(x),

        if ( ~exist('minlatix','var') )
          % Latitude indices are reversed (naturally!)
          minlatix = round( (x(1).gds.La1 - maxlat) / x(1).gds.Di ) + 2;
          maxlatix = round( (x(1).gds.La1 - minlat) / x(1).gds.Di ) + 2;
          lats = x(1).gds.La2:x(1).gds.Di:x(1).gds.La1;

          minlon(minlon < 0) = minlon(minlon < 0) + 360;
          % Using gds.Di here because SOME of these stupid meteo-hacker
          % produced "GRIB" files HAVE NO gds.Dj field in them... 
          minlonix = round( (minlon - x(1).gds.Lo1) / x(1).gds.Di ) + 1;
          maxlon(maxlon < 0) = maxlon(maxlon < 0) + 360;
          maxlonix = round( (maxlon - x(1).gds.Lo1) / x(1).gds.Di ) + 1;
          lons = x(1).gds.Lo1:x(1).gds.Di:x(1).gds.Lo2;
          lons(lons > 180) = lons(lons > 180) - 360;
        end;

        dat.n = length(x);
        dat.ny = x(1).gds.Ni;
        dat.nx = x(1).gds.Nj;
        dat.lon = lons(minlonix:maxlonix);
        dat.lat = lats(minlatix:maxlatix);
        for ix = 1:length(x)
          dat.date(ix) = datenum(x(ix).stime);
          dat.(var)(ix,1:dat.nx,1:dat.ny) = reshape(x(ix).fltarray, [dat.ny dat.nx])';
        end;
        dat.(var)(dat.(var) > 360) = nan;

        x = [];
        clear x;

        % if ( ~isempty(f) )
        %     delete(fname);
        % end;

      end; %for vix

      % Subset our world list of stations to those inside our BBOX
      if ( ~exist('stns','var') || isempty(stns) )
        stns = get_all_station_metadata;
        XV = [ min(dat.lon) max(dat.lon) max(dat.lon) min(dat.lon) ];
        YV = [ max(dat.lat) max(dat.lat) min(dat.lat) min(dat.lat) ];
        goodix = find( inside(stns.lons, stns.lats, XV, YV) );
        clear XV YV;
        stns.codes  = stns.codes(goodix);
        stns.lons   = stns.lons(goodix);
        stns.lats   = stns.lats(goodix);
        stns.depths = stns.depths(goodix);
      end;

      for ix = 1:length(stns.lons)

        stnm = lower(stns.codes{ix});
        disp(stnm);
        [lonix,latix] = gridnbhd_km(dat.lon,dat.lat,stns.lons(ix),stns.lats(ix),0);
        dat.(stnm).lonix = lonix;
        dat.(stnm).latix = latix;

        %dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'albedo', 'Albedo');
        for vix = 1:length(vars)
          var = vars{vix};
          fld = ['ww3_' dataset '_' flds{vix}];
          dat.(stnm).(fld).date = dat.date(:);
          % ??? OK - is this CORRECT this time???
          dat.(stnm).(fld).data = squeeze(dat.(var)(:,dat.(stnm).latix,dat.(stnm).lonix));
        end; %for vix

      end; %for stns

      matfname = sprintf('ww3.%s.%04d%02d.mat', dataset, yr, mo);
      disp(['Saving to ' matfname]);
      save(matfname,'dat');  
      dat = [];
      clear dat;
    toc,

  end; %for mo

end; %for yr

if ( ~isempty(f) )
    close(f);
end;












%%%%%%%%%%
%%%%%%%%%% PRIVATE FUNCTIONS
%%%%%%%%%%









function dotrial(stn)

  [ix1,ix2] = intersect_dates(stn.ndbc_sea_t.date,stn.ndbc_ncep_30a_net_heat_flux.date);

  seat.date = stn.ndbc_sea_t.date(ix1);
  seat.data = stn.ndbc_sea_t.data(ix1);
  flux.date = stn.ndbc_ncep_30a_net_heat_flux.date(ix2);
  flux.data = stn.ndbc_ncep_30a_net_heat_flux.data(ix2);

  badix = find(~isfinite(seat.data) | ~isfinite(flux.data));
  seat.date(badix) = [];
  seat.data(badix) = [];
  flux.date(badix) = [];
  flux.data(badix) = [];

  baseq0 = Q0factor(seat.data,[],1);

  ds = 1:0.1:3;
  for dix = 1:length(ds)
    q0 = baseq0 / ds(dix);
    term = real( (flux.data ./ q0) * (3600) );
    % R(dix) = corr2(seat.data,cumsum(term));
    cxy = mscohere(seat.data,cumsum(term));
    coh(dix) = nanmean(cxy);
  end;

  figure;
  plot(ds, coh);
  maxigraph;
  title(sprintf('%s - MS coher. T_s_e_a vs. \\Sigma Q_0', stn.station_name));

return;












%       'ndbc_air_sea_spechumid', ...





  %linspec = {'b.','g.','r.','c.','k.','y.','m.','bd','gd','rd','cd','kd','yd','md'};
  %linspec = {'k.-','k.:','k.-.','k.--','ko-','ko:','ko-.','ko--','ks-','ks:','ks-.','ks--','k^-','k^:','k^-.','k^--'};
  %linspec = {'k.-','k.:','k.--','ko-','ko:','ko--','kp-','kp:','kp--','k^-','k^:','k^--','k*-','k*:','k*--'};

  % linspec = {'k.-','k.:','ko-','ko:','kp-','kp:','k^-','k^:','k*-','k*:','kd-','kd:','ks-','ks:'};
  linspec = {'k.-','k.:','ks-','ks:','ko-','ko:','kp-','kp:','k^-','k^:','k*-','k*:','kd-','kd:'};








%   [stn,Q0_factor] = ...
%       station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
%                         'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lwhf',...
%                         'ndbc_ncep_26','ncep_dsrf','ncep_dlrf');
%   [stn,Q0_factor] = ...
%       station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
%                         'ndbc_barom','ndbc_sea_t','ncep_swhf','ncep_lwhf',...
%                         'ndbc_ncep_30','ncep_dsrf','ncep_dlrf','ncep_precip');
%   [stn,Q0_factor] = ...
%       station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
%                         'ndbc_barom','ndbc_sea_t','ncep_swhf','ncep_lwhf',...
%                         'ndbc_ncep_30a','ncep_dsrf','ncep_dlrf','ncep_precip',...
%                         'ndbc_wind1_dir','default','default','default','default',true);








  % [stn,Q0_factor] = ...
  %     station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,...
  %                       'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
  %                       'ndbc_bulk_26','ncep_dsrf','ndbc_bulk_rh_dlrf');





  % [stn,Q0_factor] = ...
  %     station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ncep_relhumid',...
  %                       'ndbc_barom','ndbc_sea_t','ncep_srf','ncep_lrf',...
  %                       'ncep_bulk');






%[m/s]->[kts]
        if (isfield(stn,'ncep_wind_speed'))
          stn.ncep_wind_speed.data = stn.ncep_wind_speed.data ./ 0.5144444444;
        end;
    if (isfield(stn,'ncep_wind1_speed')); stn.ncep_wind1_speed.data = stn.ncep_wind1_speed.data ./ 0.5144444444; end;







  if ( isempty(licor_valid_idx) && isempty(bic_valid_idx) )
    cldfld = 'ncep_cloud_cover';
    insfld = 'ncep_dsrf';
    uswfld = 'ncep_usrf';
    swfld = 'ncep_srf';
    warning('No in situ PAR! Using NCEP insolation...');
  elseif ( isempty(bic_valid_idx) )
    cldfld = 'licor_surf_cloud_cover';
    insfld = 'licor_surf_dsrf';
    uswfld = 'licor_surf_usrf';
    swfld = 'licor_surf_srf';
    if ( ~isfield(stn,insfld) )
      stn = station_par_to_insol(stn,'licor_surf_par','ndbc_wind1_speed', ...
                                 cldfld,insfld,uswfld,licor_valid_idx);
    end;
    if ( ~isfield(stn,swfld) )
      [ix1,ix2] = intersect_dates(stn.(insfld).date,stn.(uswfld).date);
      stn.(swfld).date = stn.(insfld).date(ix1);
      stn.(swfld).data = stn.(insfld).data(ix1) - stn.(uswfld).data(ix2);
    end;
  else
    cldfld = 'bic_surf_cloud_cover';
    insfld = 'bic_surf_dsrf';
    uswfld = 'bic_surf_usrf';
    swfld = 'bic_surf_srf';
    if ( ~isfield(stn,insfld) )
      stn = station_par_to_insol(stn,'bic_surf_par','ndbc_wind1_speed', ...
                                 cldfld,insfld,uswfld,bic_valid_idx);
    end;
    if ( ~isfield(stn,swfld) )
      [ix1,ix2] = intersect_dates(stn.(insfld).date,stn.(uswfld).date);
      stn.(swfld).date = stn.(insfld).date(ix1);
      stn.(swfld).data = stn.(insfld).data(ix1) - stn.(uswfld).data(ix2);
    end;
  end;







  stn.ncep_bulk_heat_flux_sum.date = stn.ncep_bulk_heat_flux_term.date;
  stn.ncep_bulk_heat_flux_sum.data = cumsum(stn.ncep_bulk_heat_flux_term.data);


  stn.ncep_heat_flux_sum.date = stn.ncep_heat_flux_term.date;
  stn.ncep_heat_flux_sum.data = cumsum(stn.ncep_heat_flux_term.data);


      stn.ndbc_bulk_rh_heat_flux_sum.date = stn.ndbc_bulk_rh_heat_flux_term.date;
      stn.ndbc_bulk_rh_heat_flux_sum.data = cumsum(stn.ndbc_bulk_rh_heat_flux_term.data);


      stn.ndbc_bulk_26_heat_flux_sum.date = stn.ndbc_bulk_26_heat_flux_term.date;
      stn.ndbc_bulk_26_heat_flux_sum.data = cumsum(stn.ndbc_bulk_26_heat_flux_term.data);


      stn.ndbc_bulk_30_heat_flux_sum.date = stn.ndbc_bulk_30_heat_flux_term.date;
      stn.ndbc_bulk_30_heat_flux_sum.data = cumsum(stn.ndbc_bulk_30_heat_flux_term.data);



    stn.ndbc_bulk_heat_flux_sum.date = stn.ndbc_bulk_heat_flux_term.date;
    stn.ndbc_bulk_heat_flux_sum.data = cumsum(stn.ndbc_bulk_heat_flux_term.data);

    stn.ncep_bulk_rh_heat_flux_sum.date = stn.ncep_bulk_rh_heat_flux_term.date;
    stn.ncep_bulk_rh_heat_flux_sum.data = cumsum(stn.ncep_bulk_rh_heat_flux_term.data);


  stn.sat_heat_flux_sum.date = stn.sat_heat_flux_term.date;
  stn.sat_heat_flux_sum.data = cumsum(stn.sat_heat_flux_term.data);



  % Diagnostic fields
  stn.ncep_dsrf_sum.date = stn.ncep_dsrf.date;
  stn.ncep_dsrf_sum.data = cumsum(stn.ncep_dsrf.date);
  if ( isfield(stn, 'licor_surf_dsrf') )
    stn.licor_surf_dsrf_sum.date = stn.licor_surf_dsrf.date;
    stn.licor_surf_dsrf_sum.data = cumsum(stn.licor_surf_dsrf.data);
  end;
  if ( isfield(stn, 'bic_surf_dsrf') )
    stn.bic_surf_dsrf_sum.date = stn.bic_surf_dsrf.date;
    stn.bic_surf_dsrf_sum.data = cumsum(stn.bic_surf_dsrf.data);
  end;

  stn.ncep_usrf_sum.date = stn.ncep_usrf.date;
  stn.ncep_usrf_sum.data = cumsum(stn.ncep_usrf.data);
  if ( isfield(stn, 'licor_surf_usrf') )
    stn.licor_surf_usrf_sum.date = stn.licor_surf_usrf.date;
    stn.licor_surf_usrf_sum.data = cumsum(stn.licor_surf_usrf.data);
  end;
  if ( isfield(stn, 'bic_surf_usrf') )
    stn.bic_surf_usrf_sum.date = stn.bic_surf_usrf.date;
    stn.bic_surf_usrf_sum.data = cumsum(stn.bic_surf_usrf.data);
  end;

  stn.ncep_dlrf_sum.date = stn.ncep_dlrf.date;
  stn.ncep_dlrf_sum.data = cumsum(stn.ncep_dlrf.data);
  if ( isfield(stn, 'ndbc_bulk_dlrf') )
    stn.ndbc_bulk_dlrf_sum.date = stn.ndbc_bulk_dlrf.date;
    stn.ndbc_bulk_dlrf_sum.data = cumsum(stn.ndbc_bulk_dlrf.data);
  end;

  stn.ncep_ulrf_sum.date = stn.ncep_ulrf.date;
  stn.ncep_ulrf_sum.data = cumsum(stn.ncep_ulrf.data);
  if ( isfield(stn, 'ndbc_bulk_ulrf') )
    stn.ndbc_bulk_ulrf_sum.date = stn.ndbc_bulk_ulrf.date;
    stn.ndbc_bulk_ulrf_sum.data = cumsum(stn.ndbc_bulk_ulrf.data);
  end;





















  q0f = Q0factor(stn.ndbc_sea_t.data,[],2);
  stn.ndbc_air_sea_spechumid_term.date = stn.ndbc_air_sea_spechumid.date;
  stn.ndbc_air_sea_spechumid_term.date = stn.ndbc_air_sea_spechumid.date;






  %DEBUG:    dat = query_ncep_narr_subset([], datenum(yr,mo,1), datenum(yr,mo,2));
  %%%% ??? DEBUG
  dat = query_ncep_narr_subset([], datenum(yr,mo,1), (datenum(yr,(mo+1),1)-1));
  disp(['Saving dat to ' fname]);
  save(fname, 'dat');




    pcolor(dat.lon,dat.lat,squeeze(dat.(fld)(ix,:,:)));








  minlon = -81;   maxlon = -80;
  minlat = 24.5;  maxlon = 25.5;

  [lonix








        if ( ~isempty(y) )
          % Unlike NARR or NAM, CFSR arranges latitude in rows (i.e., first)
          dat.(fld)(:,:,:) = permute(y,[1 3 2]);
        end;











       stn.sat_heat_flux_term.date, stn.sat_heat_flux_term.data, 'r-.', ...




figure; plot(1:52,[ cs.wkclim-cs.wkclim(1) ; cumsum([css.wkclim-css.wkclim(1) ; csn.wkclim-csn.wkclim(1) ; csh.wkclim-csh.wkclim(1) ; cst.wkclim-cst.wkclim(1)].*7) ]); maxigraph; xlim([0 53]); legend('\Delta T_s [^oC]','NOC Southampton','NCEP NARR','TOGA-COARE 3.0','Bulk (Smith 1988)'); titlename('SMKF1 Weekly Climatology: \Sigma_5_2_w_k Q_0 / \rho C_p h'); grid on;



figure; plot(1:52,[ cs.wkclim-cs.wkclim(1) ; cumsum([css.wkclim ; csn.wkclim ; csh.wkclim ; cst.wkclim].*7) ]); maxigraph; xlim([0 53]); legend('\Delta T_s [^oC]','NOC Southampton','NCEP NARR','TOGA-COARE 3.0','Bulk (Smith 1988)'); titlename('SMKF1 Weekly Climatology: \Sigma_5_2_w_k Q_0 / \rho C_p h'); grid on;




figure; plot(1:53,[ cs.wkclim-cs.wkclim(1) 0 ; repmat(0,[4 1]) cumsum([css.wkclim-css.wkclim(1) ; csn.wkclim-csn.wkclim(1) ; csh.wkclim-csh.wkclim(1) ; cst.wkclim-cst.wkclim(1)].*7) ]); maxigraph; xlim([0 54]); legend('\Delta T_s [^oC]','NOC Southampton','NCEP NARR','Bulk (Smith 1988)','TOGA-COARE 3.0'); titlename('SMKF1 Weekly Climatology: \Sigma_1_2_m_o Q_0 / \rho C_p h'); grid on;



figure; plot(1:366,[ cs.dyclim-cs.dyclim(1) 0 ; repmat(0,[4 1]) cumsum([css.dyclim-css.dyclim(1) ; csn.dyclim-csn.dyclim(1) ; csh.dyclim-csh.dyclim(1) ; cst.dyclim-cst.dyclim(1)]) ]); maxigraph; xlim([0 367]); legend('\Delta T_s [^oC]','NOC Southampton','NCEP NARR','Bulk (Smith 1988)','TOGA-COARE 3.0'); titlename('SMKF1 Climatology: \Sigma_3_6_5_d Q_0 / \rho C_p h'); grid on;



figure; plot(1:366,[ cs.dyclim-cs.dyclim(1) 0 ; repmat(0,[4 1]) cumsum([css.dyclim ; csn.dyclim ; csh.dyclim ; cst.dyclim].*24) ]); maxigraph; xlim([0 367]); legend('\Delta T_s [^oC]','NOC Southampton','NCEP NARR','Bulk (Smith 1988)','TOGA-COARE 3.0'); titlename('SMKF1 Climatology: \Sigma_3_6_5_d Q_0 / \rho C_p h'); grid on;




figure; plot(1:53,[ cs.wkclim-cs.wkclim(1) 0 ; repmat(0,[4 1]) cumsum([css.wkclim ; csn.wkclim ; csh.wkclim ; cst.wkclim].*24) ]); maxigraph; xlim([0 54]); legend('\Delta T_s [^oC]','NOC Southampton','NCEP NARR','Bulk (Smith 1988)','TOGA-COARE 3.0'); titlename('SMKF1 Weekly Climatology: \Sigma_5_2_w Q_0 / \rho C_p h'); grid on;




figure; plot(1:13,[ cs.moclim-cs.moclim(1) 0 ; repmat(0,[4 1]) cumsum([css.moclim ; csn.moclim ; csh.moclim ; cst.moclim].*24) ]); maxigraph; xlim([0 14]); legend('\Delta T_s [^oC]','NOC Southampton','NCEP NARR','Bulk (Smith 1988)','TOGA-COARE 3.0'); titlename('SMKF1 Monthly Climatology: \Sigma_1_2_m_o Q_0 / \rho C_p h'); grid on;




figure; plot(1:365,[ cs.dyclim-cs.dyclim(1) ; cumsum([css.dyclim-css.dyclim(1) ; csn.dyclim-csn.dyclim(1) ; csh.dyclim-csh.dyclim(1) ; cst.dyclim-cst.dyclim(1)].*24) ]); maxigraph; xlim([0 367]); legend('\Delta T_s [^oC]','NOC Southampton','NCEP NARR','Bulk (Smith 1988)','TOGA-COARE 3.0'); titlename('SMKF1 Climatology: \Sigma_3_6_5_d Q_0 / \rho C_p h'); grid on;






%     [stn,Q0_factor] = ...
%         station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ncep_relhumid',...
%                           'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
%                           'ncep_bulk_rh');

%     stn.ncep_bulk_rh_heat_flux_sum.date = stn.ncep_bulk_rh_heat_flux_term.date;
%     stn.ncep_bulk_rh_heat_flux_sum.data = cumsum(stn.ncep_bulk_rh_heat_flux_term.data);










  [ix1,ix2] = intersect_dates(stn.ndbc_bulk_rh_net_heat_flux.date,stn.ndbc_bulk_rh_latent_heat_flux.date);

  stn.ndbc_bulk_rh_net_heat_flux.data(ix1) = stn.ndbc_bulk_rh_net_heat_flux.data(ix1) - ...
      stn.ndbc_bulk_rh_latent_heat_flux.data(ix2);

  stn.ndbc_bulk_rh_latent_heat_flux.data = stn.ndbc_bulk_rh_latent_heat_flux.data .* 1e3;

  stn.ndbc_bulk_rh_net_heat_flux.data(ix1) = stn.ndbc_bulk_rh_net_heat_flux.data(ix1) + ...
      stn.ndbc_bulk_rh_latent_heat_flux.data(ix2);










  stnm = stn.station_name;
  Q0f = Q0factor(stn.ndbc_sea_t.data,[],2);


  stn.ndbc_bulk_rh_heat_flux_term.data = (stn.ndbc_bulk_rh_net_heat_flux.data ./ Q0_factor) .* (60*60);
  stn.ndbc_bulk_rh_latent_flux_term).data = (stn.ndbc_bulk_rh_latent_heat_flux.data ./ Q0_factor) .* (60*60);











  begyr = cs(1).yrs(1);
  endyr = cs(1).yrs(end);
  for cix = 2:length(cs)
    begyr = max(begyr,cs(cix).yrs(1));
    endyr = min(endyr,cs(cix).yrs(end));
  end;



  x = [];
  y = [];

  for cix = 1:length(cs)
    begyrix = find(cs(cix).yrs == begyr);
    endyrix = find(cs(cix).yrs == endyr);
    dat = [];
    for yix = begyrix:endyrix
      if ( cix == 1 )
        newyear = datenum(cs(1).yrs(yix),1,1) - 1;
        % x = [ x (cs(1).dtsvec(1,:) + newyear) ];
        x = [ x ([1:365] + newyear) ];
        dat = [ dat cs(cix).(meanfld)(yix,:)-cs(cix).(meanfld)(yix,1) ];
      else
        dat = [ dat cumsum(cs(cix).(meanfld)(yix,:)) ];
      end;
    end;
    % tmpdat = (cs(cix).(meanfld))';
    % dat = tmpdat(:)';

    y = [ y ; dat ];
  end;

  figure;
  plot(x,y);
  maxigraph;
  datetick3;
  legend(fldnms);
  titlename('Time series');










  x = datenum(begyr,1,1):(1/24):datenum(endyr,12,31);





  s = annocs(stnm);
  for fld = {'nocs_latent_heat_flux','nocs_lrf','nocs_sensible_heat_flux', ...
             'nocs_srf','nocs_net_heat_flux','nocs_heat_flux_term',...
             'nocs_heat_flux_sum',...
             'nocs_bulk_latent_heat_flux','nocs_bulk_lrf','nocs_bulk_sensible_heat_flux', ...
             'nocs_bulk_srf','nocs_bulk_net_heat_flux','nocs_bulk_heat_flux_term',...
             'nocs_bulk_heat_flux_sum',}
    stn.(fld) = s.(fld);
  end;
  s = []; clear;





  x = [];
  y = [];
  for cix = 1:length(cs)
    dts = [];
    dat = [];
    for yix = 1:length(cs(cix).yrs)
      newyear = datenum(cs(cix).yrs(yix),1,1) - 1;
      dts = [ dts(:) (cs(cix).dtsvec(:)+newyear) ];
    end;
    tmpdat = (cs(cix).datmtx)';
    dat = tmpdat(:)';

    x = [ x ; dts ];
    y = [ y ; dat ];
  end;











    stn.(fld).data = interp1(stn.(fld).date,stn.(fld).data,dts,'spline','extrap')';





   disp(['Interim save to ' matfname]);
   save(matfname,'stn');







[ ct.dyclim-ct.dyclim(1) ; cumsum(cn.dyclim) ; ...
                  cumsum(cr.dyclim) ; cumsum(c3.dyclim) ]



function stn = do30(stn)

%   [stn,Q0_factor] = ...
%       station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid',...
%                         'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
%                         'ndbc_bulk_30',...
%                         'ncep_dsrf','ndbc_bulk_rh_dlrf','ncep_precip');
  [stn,Q0_factor] = ...
      station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid',...
                        'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
                        'ndbc_bulk_30',...
                        'ncep_dsrf','ndbc_bulk_rh_dlrf');


  stn.ndbc_bulk_30_heat_flux_sum.date = stn.ndbc_bulk_30_heat_flux_term.date;
  stn.ndbc_bulk_30_heat_flux_sum.data = cumsum(stn.ndbc_bulk_30_heat_flux_term.data);

return;





function stn = annoc(stn)
%function stn = annoc(stn)
%
% Extract time series of monthly sea-surface heat fluxes from the National
% Oceanographic Centre at Southhampton (NOCS) climatology for station STN.
% Fields added to STN include .NOCS_CLOUD_COVER, .NOCS_LATENT_HEAT_FLUX,
% .NOCS_SENSIBLE_HEAT_FLUX, .NOCS_AIR_T, .NOCS_SEA_T, .NOCS_SPECHUMID,
% .NOCS_BAROM, .NOCS_LRF (net longwave), .NOCS_DRF (net shortwave),
% .NOCS_WIND_SPEED, .NOCS_NET_HEAT_FLUX, and .NOCS_HEAT_FLUX_TERM.
%
% Last Saved Time-stamp: <Fri 2010-03-26 17:12:40 Eastern Daylight Time lew.gramer>

  ALLVARS = { ...
      'at',	'air_t' ; ...
      'cldc',	'cloud_cover' ; ...
      'lhf',	'latent_heat_flux' ; ...
      'lw',	'lrf' ; ...
      'qair',	'spechumid' ; ...
      'shf',	'sensible_heat_flux' ; ...
      'slp',	'barom' ; ...
      'sst',	'sea_t' ; ...
      'sw',	'srf' ; ...
      'wspd',	'wind_speed' ; ...
            };

  for varix = 1:size(ALLVARS,1)
    var = ALLVARS{varix,1};
    fld = ['nocs_' ALLVARS{varix,2}];

    for yr = 1987:2006

      % ../data/nocs_v2_0_shf_2006.nc
      ncfname = fullfile(datapath, sprintf('nocs_v2_0_%s_%04d.nc', var, yr));
      nc = mDataset(ncfname);
      if ( ~exist('latix','var') )
        lons = nc{'lon'}(1:end);
        lats = nc{'lat'}(1:end);
      end;
      dat = nc{var}(1:end,1:end,1:end);
      close(nc); clear nc;

      if ( ~exist('latix','var') )
        [ig,latix] = min(abs(lats - stn.lat));
        [ig,lonix] = min(abs(lons - stn.lon));
      end;

      stn.(fld).date = datenum(yr,[1:12],1);
      stn.(fld).data = dat(:,latix,lonix);

      dat = []; clear dat;
    end;
  end;

return;











  yr = yr(yr ~= 2004 & yr ~= 2008);



  % Only compare means over MATCHING year-days of all good years
  for ix = 1:numel(yrdts)
    jd{ix} = unique(yrjds{ix});
  end;
  [ig,goodix,otherix] = intersect(jd{1},jd{2});
  for ix = 3:numel(yrdts)
    jd{1} = jd{1}(goodix);
    [goodix,otherix] = intersect(jd{1},jd{ix});
  end;
  goodix = find(ismember(yrjds{1},jd{1}));
  goodjds = yrjds{1}(goodix);
  yrmat(1,1:length(goodjds)) = yrdat{1}(goodix);
  for ix = 2:numel(yrdts)
    [goodix,otherix] = intersect(goodjds,yrjds{ix});
    yrmat(ix,:) = yrdat{ix}(otherix);
  end;









function [yrdts,yrdat,meandts,meandat] = anann(stn,fld)

  [yrs,mos,dys] = datevec(stn.(fld).date);
  newyrix = find(diff(yrs) > 0);

  % Use only complete years of data
  yrs = yrs(newyrix(1)+1:newyrix(end));
  mos = mos(newyrix(1)+1:newyrix(end));
  dys = dys(newyrix(1)+1:newyrix(end));

  uyrs = unique(yrs);

  for yrix = 1:length(uyrs)
    yr = uyrs(yrix);
    ix = find(yrs == yr);
    if ( length(ix) > 6000 && max(diff(stn.(fld).date(ix))) < 45 )
      yrdts{yrix} = stn.(fld).date(ix);
      yrdat{yrix} = stn.(fld).data(ix);
      meandts(yrix) = yr;
      meandat(yrix) = nanmean(yrdat{yrix});
    else
      meandat(yrix) = nan;
    end;
  end;

  badix = find(isnan(meandat));
  yrdts(badix) = [];
  yrdat(badix) = [];
  meandts(badix) = [];
  meandts(badix) = [];


  % Do analysis of variance (and show BOXPLOT) for all good years
  [P,ANOVATAB,STATS] = anova1(yrdat', cellstr(num2str(mndts')), 'on');
  titlename(['Boxplot: American Samoa MISST ' ttltag]);
  maximize_graph;
  print('-dpng', fullfile(figspath,['all-asam-boxplot-' ttltag '.png']));
  ylim(sstlims);
  P,

  figure;
  [CMPS,MNS,fh,NMS] = multcompare(STATS);
  maximize_graph(fh);
  xlim(sstlims);
  titlename([ttltag ' (Click station to test)']);
  print('-dpng', fullfile(figspath,['all-asam-multcompare-' ttltag '.png']));
  [NMS num2cell(MNS)],

  figure;
  hold on;
  envelope = [ nanmin(ssts(:)) nanmax(ssts(:)) ];
  envelope(3) = envelope(2) - envelope(1);
  plot(dts, [ nanmin(ssts) ; nanmean(ssts); nanmax(ssts) ; ...
              envelope(2) + envelope(3) + (nanmax(ssts) - nanmin(ssts)) ]);
  hold off;
  maximize_graph;
  datetick3;
  legend('Min','Mean','Max','Range');
  titlename(['American Samoa MISST inter-station ranges ' ttltag]);
  print('-dpng', fullfile(figspath,['all-asam-ranges-' ttltag '.png']));

return;









%%%% DEBUG - just for now
    stn.ncep_par.date = stn.ncep_dsrf.date;
    stn.ncep_parW.date = stn.ncep_dsrf.date;
    [stn.ncep_par.data, stn.ncep_parW.data] = ...
        insol_to_par(stn.ncep_dsrf.data);

    disp(['Saving to ' ncepmatfname '...']);
    station = stn;
    save(ncepmatfname,'station');
    station = []; clear station;
    disp('Done');
%%%% DEBUG - just for now








  else
    rhfld = 'ndbc_relhumid';
    if ( ~isfield(stn,rhfld) )
      stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t',rhfld);
    end;
    shfld = 'ndbc_spechumid';
    if ( ~isfield(stn,shfld) )
      stn = station_relhumid_to_spechumid(stn,'ndbc_air_t',rhfld,shfld);
    end;







  if ( ~isfield(stn,'ndbc_dew_t') || isempty(stn.ndbc_dew_t.date) )
    rhfld = 'ncep_relhumid';
    shfld = 'ncep_spechumid';
    warning('No in situ Relative Humidity! Using NCEP...');
  else
    rhfld = 'ndbc_relhumid';
    if ( ~isfield(stn,rhfld) )
      stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t',rhfld);
    end;
    shfld = 'ndbc_relhumid';
    if ( ~isfield(stn,shfld) )
      stn = station_relhumid_to_spechumid(stn,'ndbc_air_t',rhfld,shfld);
    end;
  end;







    % BASIC QUALITY CONTROL - we got some bad data along the line somewhere
    badix = find(stn.ncep_dlrf.data > 600);
    dts = stn.ncep_dlrf.date(badix);
    stn.ncep_dlrf.date(badix) = [];
    stn.ncep_dlrf.data(badix) = [];
    stn.ncep_dlrf.data(badix) = interp1(stn.ncep_dlrf.date,stn.ncep_dlrf.data,dts);
    stn.ncep_dlrf.date(badix) = dts;

    stn.ncep_lrf.date = stn.ncep_dlrf.date;
    stn.ncep_lrf.data = stn.ncep_dlrf.data - stn.ncep_ulrf.data;

    [ix1,ix2] = intersect_dates(stn.ncep_srf.date,stn.ncep_lrf.date);
    stn.ncep_net_heat_flux.date = stn.ncep_lrf.date(ix2);
    stn.ncep_net_heat_flux.data = stn.ncep_srf.data + stn.ncep_lrf.data + ...
        stn.ncep_latent_heat_flux.data + stn.ncep_sensible_heat_flux.data;








    % BASIC QUALITY CONTROL - we got some bad data along the line somewhere
    ncepflds = fieldnames(stn);
    ncepflds(~strmatch('ncep_',ncepflds)) = [];
    badix = find(stn.ncep_dlrf.data > 600);
    for ix = 1:length(ncepflds)
      fld = ncepflds{ix};
      stn.(fld).date(badix) = [];
      stn.(fld).data(badix) = [];
    end;






    % ndat = numel(ncep.(fld).data);
    % stn.(newfld).date(end+1:end+ndat,1) = ncep.(fld).date;
    % stn.(newfld).data(end+1:end+ndat,1) = ncep.(fld).data;





  if ( ~isfield(stn,'ncep_heat_flux_term') )
    stn.ncep_heat_flux_term.date = stn.ncep_net_heat_flux.date;
    stn.ncep_heat_flux_term.data = stn.ncep_net_heat_flux.data ./ Q0_factor;
  end;





    cldfld = 'licor_cloud_cover';


    cldfld = 'ncep_cloud_cover';


  stn = get_ncep_station(stn, 'nam');
  stn.ncep_nam_srf.date = stn.ncep_nam_dsrf.date;
  stn.ncep_nam_srf.data = stn.ncep_nam_dsrf.data + stn.ncep_nam_usrf.data;
  stn.ncep_nam_lrf.date = stn.ncep_nam_dlrf.date;
  stn.ncep_nam_lrf.data = stn.ncep_nam_dlrf.data + stn.ncep_nam_ulrf.data;
  stn = station_heat_flux(stn,'wind1_speed','air_t',rhfld,'ndbc_barom', ...
                          'sea_t','ncep_srf','ncep_lrf','ncep_nam');











  if ( ~isfield(stn,'ncep_dlrf') )
    stn = get_ncep_station(stn, 'narr');

    [ix1,ix2] = intersect_dates(stn.(insfld).date,stn.(uswfld).date);
    stn.ncep_srf.date = stn.(insfld).date(ix1);
    stn.ncep_srf.data = stn.(insfld).data(ix1) - stn.(uswfld).data(ix2);
    stn.ncep_lrf.date = stn.ncep_dlrf.date;
    stn.ncep_lrf.data = stn.ncep_dlrf.data - stn.ncep_ulrf.data;
  end;
  if ( ~isfield(stn,'ncep_net_heat_flux') )
    [ix1,ix2] = intersect_dates(stn.ncep_srf.date,stn.ncep_lrf.date);
    stn.ncep_net_heat_flux.date = stn.ncep_lrf.date(ix2);
    stn.ncep_net_heat_flux.data = stn.ncep_srf.data + stn.ncep_lrf.data + ...
        stn.ncep_latent_heat_flux.data + stn.ncep_sensible_heat_flux.data;
  end;
  if ( ~isfield(stn,'ncep_wind_speed') || ~isfield(stn,'ncep_wind_dir') )
    stn = redo_ncep
    stn = calc_ncep_wind(stn);

    [ix1,ix2] = intersect_dates(stn.ncep_wind_u.date, stn.ncep_wind_v.date);
    stn.ncep_wind_speed.date = stn.ncep_wind_u.date(ix1);
    stn.ncep_wind_speed.data = uv_to_spd(stn.ncep_wind_u.data(ix1), stn.ncep_wind_v.data(ix2));
    stn.ncep_wind_dir.date = stn.ncep_wind_u.date(ix1);
    stn.ncep_wind_dir.data = uv_to_dir(stn.ncep_wind_u.data(ix1), stn.ncep_wind_v.data(ix2));
  end;

  [stn,Q0_factor] = ...
      station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t',rhfld,'ndbc_barom', ...
                        'ndbc_sea_t','ncep_srf','ncep_lrf');

 if (0)
  [stn,Q0_factor] = ...
      station_heat_flux(stn,'ncep_wind_speed','ndbc_air_t','ncep_relhumid','ncep_barom', ...
                        'ndbc_sea_t','ncep_srf','ncep_lrf','ncep_calc');
 end;

  if ( ~isfield(stn,'ncep_heat_flux_term') )
    stn.ncep_heat_flux_term.date = stn.ncep_net_heat_flux.date;
    stn.ncep_heat_flux_term.data = stn.ncep_net_heat_flux.data ./ Q0_factor;
  end;


 if (0)
  stn = get_ncep_station(stn, 'nam');
  stn.ncep_nam_srf.date = stn.ncep_nam_dsrf.date;
  stn.ncep_nam_srf.data = stn.ncep_nam_dsrf.data + stn.ncep_nam_usrf.data;
  stn.ncep_nam_lrf.date = stn.ncep_nam_dlrf.date;
  stn.ncep_nam_lrf.data = stn.ncep_nam_dlrf.data + stn.ncep_nam_ulrf.data;
  stn = station_heat_flux(stn,'wind1_speed','air_t',rhfld,'barom', ...
                          'sea_t','ncep_srf','ncep_lrf','ncep_nam');

  dat = read_all_insolation_csvs(stnm);
  stn.sat_par.date = dat.date;
  stn.sat_par.data = dat.par;
  stn.sat_par_watt.date = dat.date;
  stn.sat_par_watt.data = dat.parW;
  stn.sat_insol_in.date = dat.date;
  stn.sat_insol_in.data = dat.sd;
  stn.sat_insol_out.date = dat.date;
  stn.sat_insol_out.data = dat.sd;
  stn.sat_net_insol.date = dat.date;
  stn.sat_net_insol.data = dat.sd - dat.su;
  stn.sat_longwave_in.date = dat.date;
  stn.sat_longwave_in.data = dat.ld;
  stn.sat_longwave_out.date = dat.date;
  stn.sat_longwave_out.data = dat.ld;
  stn.sat_net_longwave.date = dat.date;
  stn.sat_net_longwave.data = dat.ld - dat.lu;
  dat = []; clear dat;
  stn = station_heat_flux(stn,'wind1_speed','air_t',rhfld,'barom', ...
                          'sea_t','sat_net_insol','sat_net_longwave','sat');
 end;









%%%% ??? DEBUG
  flds = fieldnames(stn);
  flds = flds(strmatch('ncep_',flds));
  stn = rmfield(stn,flds);
%%%% ??? DEBUG


disp('Redo NCEP');







  stn = remove_diag_fields(stn);
  stn = remove_intrahour_winds(stn);




  flds = fieldnames(stn);
  flds = flds(strmatch('ncep_',flds));
  for fldix = 1:length(flds)
    fld = flds{fldix},
    dts = stn.(fld).date(1):(1/24):stn.(fld).date(end);
    goodix = find(isfinite(stn.(fld).data));
    stn.(fld).data = interp1(stn.(fld).date(goodix),stn.(fld).data(goodix),dts);
    stn.(fld).date = dts;
  end;




    stn.ncep_srf.date = stn.ncep_dsrf.date;
    stn.ncep_srf.data = stn.ncep_dsrf.data + stn.ncep_usrf.data;
    stn.ncep_lrf.date = stn.ncep_dlrf.date;
    stn.ncep_lrf.data = stn.ncep_dlrf.data + stn.ncep_ulrf.data;








    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'relhumid', 'Relative_humidity');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'dewp', 'Dew_point_temperature');

    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'dsrf', 'Downward_shortwave_radiation_flux');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'usrf', 'Upward_short_wave_radiation_flux_surface');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'dlrf', 'Downward_longwave_radiation_flux');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'ulrf', 'Upward_long_wave_radiation_flux_surface');

    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'ulrf', 'Latent_heat_flux');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'ulrf', 'Upward_long_wave_radiation_flux_surface');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'ulrf', 'Upward_long_wave_radiation_flux_surface');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'ulrf', 'Upward_long_wave_radiation_flux_surface');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'ulrf', 'Upward_long_wave_radiation_flux_surface');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'ulrf', 'Upward_long_wave_radiation_flux_surface');
    dat.(stnm) = _extract_ncep_field(dat.(stnm), 'ulrf', 'Upward_long_wave_radiation_flux_surface');






 dat.vars = { ...
      'Albedo', ...
      'Surface_friction_velocity', ...
      'Pressure', ...
      'Specific_humidity_height_above_ground', ...
      'Dew_point_temperature', ...
      'Relative_humidity', ...
      'u_wind_height_above_ground', ...
      'v_wind_height_above_ground', ...
      'Downward_longwave_radiation_flux', ...
      'Downward_shortwave_radiation_flux', ...
      'Evaporation', ...
      'Latent_heat_flux', ...
      'Sensible_heat_flux', ...
      'Total_cloud_cover', ...
      'Total_precipitation', ...
      'Upward_long_wave_radiation_flux_surface', ...
      'Upward_short_wave_radiation_flux_surface', ...
      };






  t.dsf(dix,1:nx,1:ny) = nc{'Downward_longwave_radiation_flux'}(1:end,1:end,1:end,1:end);
  t.dlf(dix,1:nx,1:ny) = nc{'Downward_longwave_radiation_flux'}(1:end,1:end,1:end,1:end);
  t.usf(dix,1:nx,1:ny) = nc{'Upward_short_wave_radiation_flux_surface'}(1:end,1:end,1:end,1:end);
  t.ulf(dix,1:nx,1:ny) = nc{'Upward_long_wave_radiation_flux_surface'}(1:end,1:end,1:end,1:end);




  url = sprintf( 'http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/%04d%02d/%04d%02d%02d/narr-a_221_%04d%02d%02d_%02d00_000.grb?var=Albedo&var=Surface_friction_velocity&var=Pressure&var=Specific_humidity_height_above_ground&var=Dew_point_temperature&var=Relative_humidity&var=u_wind_height_above_ground&var=v_wind_height_above_ground&var=Downward_longwave_radiation_flux&var=Downward_shortwave_radiation_flux&var=Evaporation&var=Latent_heat_flux&var=Sensible_heat_flux&var=Total_cloud_cover&var=Total_precipitation&var=Upward_long_wave_radiation_flux_surface&var=Upward_short_wave_radiation_flux_surface&latitude=25.01&longitude=-80.38&temporal=all&spatial=bb&north=27&south=24&west=-83&east=-79&east=-79&addLatLon=true', ...
                 yrs(dix), mns(dix), yrs(dix), mns(dix), dys(dix), ...
                 yrs(dix), mns(dix), dys(dix), hrs(dix) );








    url = sprintf( 'http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/%04d%02d/%04d%02d%02d/narr-a_221_%04d%02d%02d_%02d00_000.grb?var=Albedo&var=Surface_friction_velocity&var=Pressure&var=Specific_humidity_height_above_ground&var=Dew_point_temperature&var=Relative_humidity&var=u_wind_height_above_ground&var=v_wind_height_above_ground&var=Downward_longwave_radiation_flux&var=Downward_shortwave_radiation_flux&var=Evaporation&var=Latent_heat_flux&var=Sensible_heat_flux&var=Total_cloud_cover&var=Total_precipitation&var=Upward_long_wave_radiation_flux_surface&var=Upward_short_wave_radiation_flux_surface&latitude=25.01&longitude=-80.38&temporal=all&spatial=bb&north=27&south=24&west=-83&east=-79&east=-79&addLatLon=true&latitude=%g&longitude=%g&vertCoord=0&point=true', ...
                   yrs(dix), mns(dix), yrs(dix), mns(dix), dys(dix), ...
                   yrs(dix), mns(dix), dys(dix), hrs(dix), ...
                   stns.lats(stix), stns.lons(stix) );







%http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/198702/19870201/narr-a_221_19870201_1800_000.grb?var=Downward_shortwave_radiation_flux&var=Downward_longwave_radiation_flux&latitude=25.01&longitude=-80.38&time=1987-01-31T21:00:00Z&vertCoord=&accept=csv&point=true

% http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/198702/19870201/narr-a_221_19870201_2100_000.grb?var=Dew_point_temperature&var=Relative_humidity&latitude=25.01&longitude=-80.38&temporal=all&vertCoord=0&accept=csv&point=true

% http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/198702/19870201/narr-a_221_19870201_2100_000.grb?var=Dew_point_temperature&var=Relative_humidity&latitude=25.01&longitude=-80.38&temporal=all&vertCoord=1&accept=csv&point=true

% http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/198702/19870201/narr-a_221_19870201_2100_000.grb?var=Dew_point_temperature&var=Relative_humidity&latitude=25.01&longitude=-80.38&temporal=all&vertCoord=1&accept=csv&point=true








      'Albedo[0]', ...
      'Surface_friction_velocity[0]', ...

      'Pressure[0][0]', ...
      'Specific_humidity_height_above_ground[0][0]', ...

      'Dew_point_temperature[0][0]', ...
      'Relative_humidity[0][0]', ...

      'u_wind_height_above_ground[0][0]', ...
      'v_wind_height_above_ground[0][0]', ...

      'Downward_longwave_radiation_flux[0]', ...
      'Downward_shortwave_radiation_flux[0]', ...
      'Evaporation[0]', ...
      'Latent_heat_flux[0]', ...
      'Sensible_heat_flux[0]', ...
      'Total_cloud_cover[0]', ...
      'Total_precipitation[0]', ...
      'Upward_long_wave_radiation_flux_surface[0]', ...
      'Upward_short_wave_radiation_flux_surface[0]', ...







  ALLVARS = { ...
      'Dew_point_temperature[0][0]', ...
      'Downward_longwave_radiation_flux[0]', ...
      'Downward_shortwave_radiation_flux[0]', ...
      'Evaporation[0]', ...
      'Latent_heat_flux[0]', ...
      'Precipitation_rate[0]', ...
      'Pressure_nearest_grid_point[0]', ...
      'Relative_humidity[0][0]', ...
      'Sensible_heat_flux[0]', ...
      'Specific_humidity_height_above_ground[0][0]', ...
      'Total_cloud_cover[0]', ...
      'Total_precipitation[0]', ...
      'Upward_long_wave_radiation_flux_surface[0]', ...
      'Upward_short_wave_radiation_flux_surface[0]', ...
      'u_wind_height_above_ground[0][0]', ...
      'v_wind_height_above_ground[0][0]', ...
            };







  % NCEP_NARR_DAILIES variables, equivalent NCEP NARR (32km) variables
  ALLVARS = { ...
      { 'apcp',		'Total_precipitation'				}, ...
      { 'dlwrfsfc',	'Downward_longwave_radiation_flux'		}, ...
      { 'dswrfsfc',	'Downward_shortwave_radiation_flux'		}, ...
      { 'dpt2m',	'Dew_point_temperature[0]'			}, ...
      { 'evpsfc',	'Evaporation'					}, ...
      { 'fricvsfc',	'AccPrec'					}, ...
      { 'lhtfl',	'Latent_heat_flux'				}, ...
      { 'rh2m',		'Relative_humidity[0]'				}, ...
      { 'shtfl',	'Sensible_heat_flux'				}, ...
      { 'spfh2m',	'Specific_humidity_height_above_ground[0]'	}, ...
      { 'tcdc',		'Total_cloud_cover'				}, ...
      { 'ulwrfsfc',	'Upward_long_wave_radiation_flux_surface'	}, ...
      { 'uswrfsfc',	'Upward_short_wave_radiation_flux_surface'	}, ...
      { 'ugrd10m',	'u_wind_height_above_ground[0]'			}, ...
      { 'vgrd10m',	'v_wind_height_above_ground[0]'			}, ...
      { 'uflxsfc',	'u_wind_stress_surface?'			}, ...
      { 'vflxsfc',	'v_wind_stress_surface?'			}, ...
            };








    for vix = 1:length(flds.vars)
      dataquery = sprintf('%s?%s[%g:%g][%g:%g]', base_dataquery, flds.vars{vix}, ...
                          minlonix, maxlonix, minlatix, maxlatix);
      nc = mDataset(dataquery);
      flds.(lower(flds.varnames{vix}))(dix,:,:) = nc{flds.varnames{vix}}(:);
      close(nc);
      clear nc;
    end;






  vix = 1;
  flds.dataquery = sprintf('%s?%s[%g:%g][%g:%g]', flds.dataquery, flds.vars{vix}, ...
                           minlonix, maxlonix, minlatix, maxlatix);
  for vix = 2:length(vars)
    flds.dataquery = sprintf('%s,%s[%g:%g][%g:%g]', base_dataquery, flds.vars{vix}, ...
                             minlonix, maxlonix, minlatix, maxlatix);
  end;






  disp('');
  disp('Reviewing will start in 3 secs: Ctrl-C to stop...');
  pause(3);





station = verify_variable(station, 'wind1_speed_1_day_average');
station = filter_gaps(station, 'wind1_speed', 'wind1_speed_1_day_average',1,1);

station = verify_variable(station, 'air_t_1_day_average');
station = filter_gaps(station, 'air_t', 'air_t_1_day_average',1,1);

station = verify_variable(station, 'air_t_dewp_1_day_average');
station = filter_gaps(station, 'air_t_dewp', 'air_t_dewp_1_day_average',1,1);

station = verify_variable(station, 'barom_surf_1_day_average');
station = filter_gaps(station, 'barom_surf', 'barom_surf_1_day_average',1,1);

station = verify_variable(station, 'sea_t_1_day_average');
station = filter_gaps(station, 'sea_t', 'sea_t_1_day_average',1,1);

station = verify_variable(station, 'sea_t_5_hour_findiff');
station = filter_gaps(station, 'sea_t', 'sea_t_5_hour_findiff',1,(5/24));

station = verify_variable(station, 'licor_surf_par_1_day_average');
station = filter_gaps(station, 'licor_surf_par', 'licor_surf_par_1_day_average',1,1);

station = verify_variable(station, 'sat_par_1_day_average');
station = filter_gaps(station, 'sat_par', 'sat_par_1_day_average',1,1);

station = verify_variable(station, 'sat_insol_in_1_day_average');
station = filter_gaps(station, 'sat_insol_in', 'sat_insol_in_1_day_average',1,1);

station = verify_variable(station, 'sat_insol_out_1_day_average');
station = filter_gaps(station, 'sat_insol_out', 'sat_insol_out_1_day_average',1,1);

station = verify_variable(station, 'sat_net_insol_1_day_average');
station = filter_gaps(station, 'sat_net_insol', 'sat_net_insol_1_day_average',1,1);

station = verify_variable(station, 'sat_longwave_in_1_day_average');
station = filter_gaps(station, 'sat_longwave_in', 'sat_longwave_in_1_day_average',1,1);

station = verify_variable(station, 'sat_longwave_out_1_day_average');
station = filter_gaps(station, 'sat_longwave_out', 'sat_longwave_out_1_day_average',1,1);

station = verify_variable(station, 'sat_net_longwave_1_day_average');
station = filter_gaps(station, 'sat_net_longwave', 'sat_net_longwave_1_day_average',1,1);



station = verify_variable(station, 'relhumid_1_day_average');
station = filter_gaps(station, 'relhumid', 'relhumid_1_day_average',1,1);






% dat = read_insolation_csv(['../data/insolation_' stnm '.csv'],stnm);




%     % For now, just get the point nearest our site of interest
%     stnlat = 25.01; stnlon = -80.38;
%     ix = find( abs(stnlon-lons)<=(dlon/2) & abs(stnlat-lats)<=(dlat/2) );
%     result.lon = unique(lons(ix));
%     result.lat = unique(lats(ix));
%     result.date = dts(ix);
%     result.parW = parW(ix);
%     result.par = par(ix);
%     result.sd = sd(ix);
%     result.su = su(ix);
%     result.ld = ld(ix);
%     result.lu = lu(ix);
%     result.cld = cld(ix);

%     % For now, just take a median of the points around our site of interest
%     stnlat = 25.01; stnlon = -80.38;
%     ix = find( abs(stnlon-lons)<=(dlon*1.5) & abs(stnlat-lats)<=(dlat*1.5) );
%     result.lon = unique(lons(ix));
%     result.lat = unique(lats(ix));
%     result.date = unique(dts(ix));
%     for dtix = 1:length(result.date)
%       dt = result.date(dtix);
%       nowix = intersect(find( dts == dt ), ix);
%       result.parW = median(parW(nowix));
%       result.par = median(par(nowix));
%       result.sd = median(sd(nowix));
%       result.su = median(su(nowix));
%       result.ld = median(ld(nowix));
%       result.lu = median(lu(nowix));
%       result.cld = median(cld(nowix));
%     end;

%     % Form 3D meshgrids - one layer per date - of entire sample region
%     result.date = unique(dts);
%     result.lon = min(lons):dlon:max(lons);
%     result.lat = min(lats):dlat:max(lats);
%
%     [dategrd,longrd,latgrd]  = meshgrid(result.date, result.lon, result.lat);
%
%     result.parW = griddata3(dts, lons, lats, parW, dategrd, longrd, latgrd);
%     result.par = griddata3(dts, lons, lats, par, dategrd, longrd, latgrd);
%     result.sd = griddata3(dts, lons, lats, sd, dategrd, longrd, latgrd);
%     result.su = griddata3(dts, lons, lats, su, dategrd, longrd, latgrd);
%     result.ld = griddata3(dts, lons, lats, ld, dategrd, longrd, latgrd);
%     result.lu = griddata3(dts, lons, lats, lu, dategrd, longrd, latgrd);
%     result.cld = griddata3(dts, lons, lats, cld, dategrd, longrd, latgrd);










[ixis,ixsi] = intersect_dates(station.sat_insol_out_1_day_average.date, ...
                              station.sea_t_1_day_average.date);
[ixia,ixai] = intersect_dates(station.sat_insol_out_1_day_average.date, ...
                              station.air_t_1_day_average.date);
[ixsa,ixas] = intersect_dates(station.sea_t_1_day_average.date, ...
                              station.air_t_1_day_average.date);

ixi = intersect(ixis,ixia);
ixa = intersect(ixai,ixas);
ixs = intersect(ixsi,ixsa);

regresp = station.sat_insol_out_1_day_average.data(ixi)';
regdata = [ station.sea_t_1_day_average.data(ixs); ...
           station.air_t_1_day_average.data(ixa) ]';
stats = regstats(regresp, regdata, 'quadratic');







datetick2('x', 2, 'keepticks', 'keeplimits'); set_datetick_cursor;



% [ix1,ix2] = intersect_dates(station.licor_surf_par.date,dat.date);





    % For now, just get the point nearest our site of interest
    stnlat = 25.01; stnlon = -80.38;
    ix = find( abs(stnlon-lons)<=(dlon/2) & abs(stnlat-lats)<=(dlat/2) );
    result.lon = unique(lons(ix));
    result.lat = unique(lats(ix));
    result.date = dts(ix);
    result.parW = parW(ix);
    result.par = par(ix);
    result.sd = sd(ix);
    result.su = su(ix);
    result.ld = ld(ix);
    result.lu = lu(ix);
    result.cld = cld(ix);








    matsize = [numel(result.date) numel(result.lat) numel(result.lon)];
    result.parW = repmat(nan, matsize);
    result.par = repmat(nan, matsize);
    result.sd = repmat(nan, matsize);
    result.su = repmat(nan, matsize);
    result.ld = repmat(nan, matsize);
    result.lu = repmat(nan, matsize);
    result.cld = repmat(nan, matsize);




  % Allocate an array of NaNs cast/sized as single-floats. (The mediocre
  % MATLAB Map toolbox would force us to use doubles - TWICE as big...)
  [nlats, nlons, refvec] = sizem([minlat maxlat], [minlon maxlon], (1/dl));
  % grd = repmat(cast(nan,'single'), [nlats nlons]);
  grd = repmat(nan, [nlats nlons]);

  % Load data from file, but only 10,000 rows/points at a time
  fid = fopen(fname,'r');
  while ( isempty(ferror(fid)) && ~feof(fid) )
    rawd = fscanf(fid,'%g,%g,%g\n',[3 10000])';
    badix = find(rawd(:,1) == 0 | rawd(:,2) == 0);
    rawd(badix,:) = [];
    if ( ~isempty(rawd) )
      % d = cast(rawd, 'single');
      d = rawd;
      rawd = []; clear rawd;
      grd = imbedm(d(:,2), d(:,1), d(:,3), grd, refvec);
      d=[]; clear d;
    end;
  end;
  fclose(fid);







  lats = minlat:dlat:maxlat; nlats = length(lats);
  lons = minlon:dlon:maxlon; nlons = length(lons);

  [grd, refvec] = vec2mtx(lats, lons, (1/dlon));
  grd = repmat(cast(nan,'single'), [nlons nlats]);





1;

m = csvread('../data/isobaths-mlrf1.csv');

dlon = 8.3300e-004;
dlat = 8.3300e-004;

lons = min(m(:,1)):dlon:max(m(:,1));
lats = min(m(:,2)):dlat:max(m(:,2));

[LON,LAT] = meshgrid(lons, lats);

Z = griddata(m(:,1),m(:,2),m(:,3), LON,LAT);
Z = cast(Z, 'single');

figure;
contour(LON,LAT,Z,[-2:-2:-10 -15:-5:-30 -100:-50:-200]);






%y = cosd(bbox(3)) * abs(bbox(3)-bbox(4)) / abs(bbox(1)-bbox(2));





  fh = figure;
  set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% COMMENTED OUT FOR SPEED
%   ah(1) = subplot('position', [0.10 0.24 0.88 0.67]);
  ah(1) = subplot('position', [0.05 0.05 0.90 0.90]);
  hold on;

  % Assign colormap
  colormap(ah(1), cmap);

  % Surface-contour the property onto a map of Straits of Florida
%   pcolor(LONS, LATS, prp);
  pcolor(LONS, LATS, log(prp));
  % CAREFUL: Non-interp maps look awful - but this IS a smoothing operation!
  shading('interp');
  set_pcolor_cursor;

  % Show (rough) h=-30 and h=-220 isobaths of the bottom topography
  map_sofla([LONS(1,1) LONS(end,end) LATS(1,1) LATS(end,end)], [-30 -220]);
  % Plot locations of all SEAKEYS stations
  plot_seakeys;
  view(2);

  cbh = colorbar('EastOutside');
  set(cbh, 'YScale', 'log');

  set(ah(1), 'xlim', [LONS(1,1) LONS(end,end)]);
  set(ah(1), 'ylim', [LATS(1,1) LATS(end,end)]);
  if ( ~isempty(prplim) )
%     set(ah(1), 'clim', prplim);
    set(ah(1), 'clim', 'default');
  else
    set(ah(1), 'clim', 'default');
  end;


  % Outline contiguous gridpoints containing values of interest
  if ( ~isempty(prppeaks) )
    %peakixes = outline_peaks(fh, LONS, LATS, prp, prppeaks);
  end;

  ttlstr = sprintf('%s - %s', hr, ttl);
  th = title(strrep(ttlstr, '_', '\_'));
  set(fh, 'Name', ttlstr);
  figfname = fullfile(figspath, sprintf('%s.%s.png', hr, ylbl));
%   print('-dpng', figfname);













  % Surface-contour the property onto a map of Straits of Florida
  surf(LONS, LATS, prp);
  set(gca,'ZScale','log');
  % CAREFUL: Non-interp maps look awful - but this IS a smoothing operation!
  shading('interp');
  set_surf_cursor;






  % Plot actual vectors, color coded for magnitude
  fh = figure;
  set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
  hold on;
  plot3(lon, lat, spd, '.');
  quiver(lon, lat, u, v, 0.5);
  map_sofla([LONS(1,1) LONS(end,end) LATS(1,1) LATS(end,end)]);
  title(sprintf('%s - surface currents', hr));




        found = 0;
        for boxix = 1:length(bbox_or_stanm_list)
          pname = fullfile(datapath, [stanm{boxix} '-' fname]);
          if ( exist(pname, 'file') )
            found = found + 1;
          end;
        end;
        if ( found == length(bbox_or_stanm_list) )
          






        sst = read_avhrr_subset(sstbytes);


        % Calculate SST from 1-byte PNG (colormap) values
        sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;
        % Assume that BOTH 254 and 255 are cloud mask values(??)
        sst(sstbytes >= 254) = nan;
        % Just in case our assumptions about cloud mask are wrong!
        sst(minval > sst | sst > maxval) = nan;







    [ig,sat,ig,ig,ig,ig,ig,ig,ig,ig] = parseusfurl(pname);



          % % Make all SSTs relative to the station site, if possible
          % sst = sst - nanmean(ctrsst(:));





% % Reconstruct seasonal data from saved annual data...
% fname = sprintf('%s-ssts-%04d.mat', stanm, year);
% load(fname, 'ssts', 'dts');
% keepix = find(indts(1) <= dts & dts <= indts(end));
% ssts = ssts(keepix,:,:);
% dts = dts(keepix);






save(fname, 'sst');
load(fname, 'sst');



sstdims = [81 81];




% load(sprintf('%s-sst-%04d.mat', stanm, year), 'sst');





        % try, imwrite(sstbytes, pname);
        % catch, end;







          % Make all SSTs relative to the station site, if possible
          sst = sst - sst(ctr(1), ctr(2));










        if ( length(find(isnan(ctrsst))) < (0.50*numel(ctrsst)) )
          % Make all SSTs relative to the station site, if possible
          if ( ~isnan(sst(ctr(1), ctr(2))) )
            sst = sst - sst(ctr(1), ctr(2));
          else
            sst = sst - nanmean(ctrsst);
          end;

          ssts = { ssts{:} sst };
          urls = { urls{:} url };
        else
          skipped = skipped + 1;
        end;







        % Only take images with some cloud-free data surrounding our site!
        ctr = round(boxradius/2);
        ctrsst = sst((ctr-8+1):(ctr+8),(ctr-8+1):(ctr+8));
        if ( length(find(isnan(ctrsst))) < (0.66*numel(ctrsst)) )
          fprintf('Clear: "%s"\n', url);
          urls = { urls{:} url };
        else
          skipped = skipped + 1;
        end;
        clear sst;
        clear sstbytes;







        url = sprintf('%s/%s', baseurl, fname);
        %sstbytes = imread(url);

        fpath = fullfile(datapath, fname);

        if ( exist(fpath, 'file') )
          try
            sstbytes = imread(fpath);
          catch
            % Last ditch effort - delete corrupt file and load from URL anyway
            fprintf(2, '\n');
            warning('Corrupt cached file? Trying URL "%s"...', url);
            delete(fpath);
          end;
        end;
        if ( ~exist(fpath, 'file') )
          [fpath, fstatus] = urlwrite(url, fpath);
          if ( fstatus == 0 )
            fprintf(2, '\n');
            warning('SKIPPED! Failed to download "%s"...', url);
            try
              delete(fpath);
            catch
            end;
            skipped = skipped + 1;
            continue;
          end;
          sstbytes = imread(fpath);
        end;








          fprintf('Clear: "%s"\n', fname);

          fprintf('Cloudy! "%s"\n', fname);



  try
    sstbytes = imread(fpath);
  catch
    % Last ditch effort - delete corrupt file and just load from URL
    fprintf(2, '\n');
    warning('Corrupt cached file? Trying URL "%s"...', url);
    try
      delete(fpath);
      sstbytes = imread(url);
    catch
      fprintf(2, '\n');
      warning('SKIPPED! Failed to download "%s"...', url);
    end;
  end;








  try
    sstbytes = imread(fpath);
  catch
    % Last ditch effort - delete corrupt file and just load from URL
    fprintf(2, '\n');
    warning('Corrupt cached file? Trying URL "%s"...', url);
    try
      delete(fpath);
      sstbytes = imread(url);
    catch
      fprintf(2, '\n');
      warning('SKIPPED! Failed to download "%s"...', url);
    end;
  end;






  if ( ~exist('figspath', 'var') || isempty(figspath) )
    figspath = fullfile(pathroot, '../figs', '');
  end;
  datapath = '../data';





  % NOTE: Masked values NO LONGER 0 after we subtract the mean
  %goodsst = goodsst - tsmean;






  % SST at SEAKEYS station
  stn_sst = sst(stn_x, stn_y);




  % Return a field with zero spatial mean
  sst = sst - mean(sst(~isnan(sst)));






%   % Relate all SST to pixel at SEAKEYS station
%   sst = sst - stn_sst;

  %%%% DEBUG
  %%%% DEBUG
  %%%% DEBUG
  sst = sst - mean(sst(~isnan(sst)));







%
% Plot SOM results ("units") according to SOM Map Shape
%


figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
for ix = 1:nmaps
  subplot(mapdims(2), mapdims(1), ix);

  pcolor(LONS, LATS, reshape(results(ix, :), sstdims));
  shading('flat');
  % For pcolor and similar, let user see SST values when in datacursor mode
  h = datacursormode(gcf);
  set(h, 'UpdateFcn', @xyzc_select_cb);
%   set(gca, 'zlim', [minc maxc]);
%   set(gca, 'clim', [minc maxc]);

  % contourf(LONS, LATS, reshape(results(ix, :), sstdims));
  % set(gca, 'clim', [minc maxc]);

  title(sprintf('#%d (%d wks - %0.1f%%)', ...
                ix, wksmatched(ix), pctmatched(ix)));
end;
% suptitle from MATLAB online contrib
ax = suplabel(sprintf('SOM Units for mean(%s), 1993-2009 %s (%d weeks)', ...
                      meandim, period, size(sst,1)), 't');
% set(ax,'clim',[minc maxc]);
% colorbar('peer', ax);








% Calculate number of matched weeks for each "mode"
for ix = 1:size(results,1)
  wksmatched(ix) = length(find(bmus == ix));
  pctmatched(ix) = (wksmatched(ix) / length(bmus)) * 100;
end;


% Sort BMUs by number of weeks, from most-matched to least
[ign, mode_order] = sort(wksmatched, 2, 'descend');




%   % Set masked values BACK to 0 now that we have removed the mean
%   goodsst(~isfinite(sst)) = 0;





  % If a pixel is NaN for every week in our dataset, it's probably land
  % (Should save the actual land masks provided by PNG images instead!)
  landmask = false(size(sst, 2), 1);
  for ix = 1:size(sst, 2)
    if ( all(isnan(sst(:, ix))) )
      landmask(ix) = true;
    end;
  end;





  % imagesc(reshape(results(ix, :), sstdims), [minc maxc]);
  % imagesc(reshape(pcadat(ix, :), sstdims), [minc maxc]);




% for ix = 1:size(results,1)
%   fhs(ix) = figure;
%   imagesc(reshape(results(ix, :), sstdims));
%   set(fhs(ix), 'units', 'normalized', 'outerposition', [0 0 1 1]);
%   titlename(['Variable ' num2str(ix)]);
%   pause(0.5);
% end;

% while(1)
%   for fi = 1:length(fhs)
%     fh = fhs(fi);
%     figure(fh);
%     pause(0.5);
%   end;
%   pause(1);
% end;







% yrs = 1994:2008;
% wks = 24:38;


  yrs = repmat(1994:2008', [15 1]);
  wks = repmat(, [1 15]);
  yrs = [repmat(1993, [1 3]) yrs(:)'];
  wks = [50:52   wks];





  yrs = repmat(1994:2008', [52 1]);
  wks = repmat(1:52, [1 15]);
  yrs = [repmat(1993, [1 14]) yrs(:)'];
  wks = [39:52   wks];





  % mlrf1_x = 660; mlrf1_y = 1162;
  % smkf1_x = 705; smkf1_y = 1089;
  stn_x = smkf1_x; stn_y = smkf1_y;






  % Draw landmask in grey
  pcolor(LONS, LATS, (double(lnd).*max(prp(:))));




  % surf(LONS, LATS, prp);
  % contourf(LONS, LATS, prp);




  % Assign colormap, with gray added at the end
  colormap(ah(1), cmap);
  cmapclrs = colormap(ah(1));
  cmapclrs(end+1,1:3) = [0.5 0.5 0.5];
  colormap(ah(1), cmapclrs);






  th = title(sprintf('%s - %s', hr, ttl));
  set(th, 'Interpreter', 'none');





  chlix = find( sst_s.cloud_mask == 0 & ~bitget(chl_s.l2_flags,1) & ...
                ~bitget(chl_s.l2_flags,2) );
%                 ~bitget(chl_s.l2_flags,2) & ~bitget(chl_s.l2_flags,16) );




  sstix = find( sst_s.cloud_mask == 0 & ~bitget(sst_s.l2_flags,1) & ...
                ~bitget(sst_s.l2_flags,2) );
%                 ~bitget(sst_s.l2_flags,2) & ~bitget(sst_s.l2_flags,29) );
  sst(sstix) = (sst_s.seadas_sst(sstix) .* 0.005);
  sst(5 > sst | sst > 40) = nan;

  % Borrow cloud mask from SST data!
  chlix = find( sst_s.cloud_mask == 0 & ~bitget(chl_s.l2_flags,1) & ...
                ~bitget(chl_s.l2_flags,2) );
%                 ~bitget(chl_s.l2_flags,2) & ~bitget(chl_s.l2_flags,16) );
  % If we wanted to screen out obviously contaminated chl values








function s = query_dods(baseurl, queryargs, nrows, ncols)
%function s = query_dods(baseurl, queryargs, nrows, ncols)
%
% Submit a query to a DoDS (text access to binary data subset) server. Second
% arg 'queryargs' should have no file separators ("/" or "\") in it, so that
% it is suitable for conversion into a local filename where the ASCII query
% result will be cached for future requests. Returns a struct 's' containing
% one field for each variable succesfully queried by the 'queryargs' string.
% Each field will contain a numeric array - of double for most variables, but
% of 'uint32' if the variable name happens to be 'l2_flags' (HACK ALERT).
%
% Last Saved Time-stamp: <Thu 2009-02-12 13:19:31 Eastern Standard Time gramer>

  s = [];

  datapath = './data';
  queryfname = fullfile(datapath, [regexprep(queryargs, '\W+', '_') '.txt']);

  % DEBUG:  fprintf(1, 'querystr=%s, ncols=%d\n', querystr, ncols);
  % DEBUG:
  tic;
  if ( ~exist(queryfname, 'file') )
    querystr = sprintf('%s/%s', baseurl, queryargs);
    % DEBUG:
    disp(['Reloading ' querystr]);
    urlwrite(querystr, queryfname);
  end;
  % DEBUG:  disp('urlwrite'); toc;

  res = fileread(queryfname);
  % DEBUG:  disp('fileread'); toc;

  fmt = repmat(',%f', [1 ncols]); fmt = [ '%[^[][%d]' fmt ];
  flgfmt = repmat(',%d32', [1 ncols]); flgfmt = [ '%[^[][%d]' flgfmt ];

  nline = 0;
  while ( length(res) > 0 )
    nline = nline + 1;
    % DEBUG:
    if (mod(nline,100)==0); disp('start'); toc; end;
    [ln, res] = strtok(res, sprintf('\n'));
    % DEBUG:
    if (mod(nline,100)==0); disp('strtok'); toc; end;
    if ( ~isempty(ln) )
      clear cstrs;
      cstrs = textscan(ln, fmt);
      % DEBUG:
      if (mod(nline,100)==0); disp('textscan'); toc; end;
      if ( length(cstrs) >= 3 && ~isempty(cstrs{3}) )
        var = cstrs{1}{:};
        % row = cstrs{2} + 1;
        row = nrows - cstrs{2};

        % Handle bit flags specially
        % (Idiotic mazzafrazzin MATLAB)
        if ( strcmpi(var, 'l2_flags') )
          if ( ~isfield(s, var) )
            s.(var) = repmat(intmax('uint32'), [nrows ncols]);
          end;
          clear cstrs;
          cstrs = textscan(ln, flgfmt);
          s.(var)(row,1:length(cstrs)-2) = typecast([cstrs{3:end}],'uint32');

        else
          if ( ~isfield(s, var) )
            s.(var) = repmat(nan, [nrows ncols]);
          end;
          s.(var)(row,1:length(cstrs)-2) = [cstrs{3:end}];

        end;
      end;
      % DEBUG:
      if (mod(nline,100)==0); disp('typecast'); toc; end;
    end;
  end;

  % DEBUG:  disp('end while'); toc;

return;










  % Bounds for WERA's-eye view of Straits of Florida
  set(gca, 'xlim', [-80.50 -79.10]);
  set(gca, 'ylim', [+24.80 +25.90]);





  chlix = find(sst_s.cloud_mask == 0);
  % sstix = find(sst_s.cloud_mask == 0 & bitget(sst_s.l2_flags,29) == 0);





  % MLRF1: -80.38 +25.01; 1167:1169 330:332
  % lon1=-80.7; lon2=-80.1; lat1=24.7; lat2=25.3;



  %%%% DEBUG
  %%%% DEBUG
  %%%% DEBUG
  %latix1 = 200;  latix2 = 500;  lonix1 = 20;  lonix2 = 320;
  % LR corner of image
  % latix1 = 689;  latix2 = 989;  lonix1 = 1019;  lonix2 = 1319;




  % The subset (or intersection) of interest to us

  lon1 = boundingbox(1);
  lon2 = boundingbox(2);
  lat1 = boundingbox(3);
  lat2 = boundingbox(4);

  %%%% DEBUG
  %%%% DEBUG
  %%%% DEBUG
  % MLRF1: -80.38 +25.01; 1167:1169 330:332
  lon1=-80.7; lon2=-80.1; lat1=24.7; lat2=25.3;

  % 1:990 , 1:1320
  lonix1 = floor((lon1 - minlon) / dlon);
  lonix2 = ceil((lon2 - minlon) / dlon);
  latix1 = floor((lat1 - minlat) / dlat);
  latix2 = ceil((lat2 - minlat) / dlat);

  %%%% DEBUG
  %%%% DEBUG
  %%%% DEBUG
  % MLRF1: 330:332 1167:1169
  latix1 = 200;  latix2 = 500;  lonix1 = 1000;  lonix2 = 1300;

  % Shift user values to be at the appropriate gridpoints
  lon1 = minlon + (lonix1 * dlon);
  lon2 = minlon + (lonix2 * dlon);
  lat1 = minlat + (latix1 * dlat);
  lat2 = minlat + (latix2 * dlat);







  badflgs = [1 2 4 5 6 9 10 15 16];
  badflgs = [1 2 4 5 6 9 10];
  flgs = ones(size(sst_s.l2_flags));
  for flg = badflgs
    flgs(bitget(sst_s.l2_flags, flg) ~= 0) = 0;
  end
%   figure;
%   pcolor(LONS, LATS, flgs);
%   % map_sofla([LONS(1,1) LONS(end,end) LATS(1,1) LATS(end,end)]);
%   flgs,






  %%%% DEBUG
  %%%% DEBUG
  %%%% DEBUG
  % MLRF1: 330:332 1167:1169
  lonix1 = 1000;  lonix2 = 1300;  latix1 = 200;  latix2 = 500;







  % Subdirectories galore!
  %http://cyclops.marine.usf.edu/modis/level3/husf/florida/2008/222/1km/pass/final/MODIS.2008222.153656.florida.seadas_sst.png
  %http://www.imars.usf.edu/dods-bin/nph-dods/modis/level3/husf/florida/2009/003/1km/pass/final/MODIS.2009003.160747.florida.seadas_sst.hdf.html






  dlat = 0.00909090909090909;  dlon = 0.00909090909090909;
  minlat=22;  minlon=-91;  maxlat=31;  maxlon=-79;

  latix = round((lat - minlat) / dlat);
  lonix = round((lon - minlon) / dlon);

  % MLRF1: 330:332 1167:1169
  latix1 = 330; latix2 = 332;
  lonix1 =1167; lonix2 =1169;







  sst.seadas_sst = [];
  sst.l2_flags = [];
  sst.cloud_mask = [];

  chl.chlor_a = [];
  chl.l2_flags = [];
  chl.tsm_clark = [];




  while ( length(res) > 0 )
    [ln, res] = strtok(res, sprintf('\n'));
    if ( isempty(ln) )
      continue;
    end;
    [cstr, pos] = textscan(ln, '%[^[][%d],%*[^\n]')
    if ( ~isempty(cstr{2}) )
      var = cstr{1};
      row = cstr{2};
      ln = ln(pos+1:end);
      vals = textscan(ln, '%g', latix2 - latix1 + 1);
      sst.(var)(row,1:length(vals)) = [vals{:}];
    end;
  end;







  while ( length(res) > 0 )
    [ln, res] = strtok(res, sprintf('\n'));
  end;




  endl = [1 strfind(res, sprintf('\n'))];
  if ( endl(end) ~= length(res) )
    endl(end+1) = length(res);
  end;
  for ix = 2:length(endl)
    ln = res(endl(ix-1):endl(ix));
    fprintf(1, '"%s"\n', ln);
  end;







  position = get(0, 'ScreenSize');
  if ( position(4) >= 890 )
    % Lew's Desktop
    position = [1 61 1280 890];
  else
    % Lew's Laptop
    position = [1 1 1280 726];
  end;

% . . .

  % Use screen information to 'maximize' plot - stupid MATLAB
  set(gcf, 'Position', position);










  elseif ( numel(wdts) ~= numel(dts) )
    error('Size mismatch! Winds cannot be shown.');




  fprintf(2,  'Winds %d, OWP %d\n', numel(wdts), numel(dts));





  % Display local wind forcing beneath time series
  [wdts,wspd] = load_g2_data('FWYF1-WIND1-SPEED.csv', min(dts(:)), max(dts(:)));
  [wdts,wdir] = load_g2_data('FWYF1-WIND1-DIR.csv', min(dts(:)), max(dts(:)));
  wu = wspd .* sind(wind2currdir(wdir));
  wv = wspd .* cosd(wind2currdir(wdir));
  [taux, tauy] = wstress(wu, wv, 43);
  ah(2) = subplot('position', [0.10 0.05 0.88 0.16]);
  % feather(ah(2), wu, wv);
  featherdt(ah(2), wdts, taux, tauy);
  taux = taux .* 0.1; tauy = tauy .* 0.1;
%   xlim(ah(2), [1 length(wdts)]); ylim([-10 10]);
  xlim(ah(2), [min(wdts) max(wdts)]); ylim([-10 10]);
  set(ah(2), 'XTickLabel', []);
  xlabel('FWYF1 wind stress (\tau: N/m^2)');




  % Quiver-plot maximum current vector along each latitude
  [Vmax, VmaxIx] = max(SPD, [], 2);
  % Stupid f'in MATLAB; stupid f'in FORTRAN
  for ix = 1:length(VmaxIx)
    qlon(ix) = LONS(ix, VmaxIx(ix));
    qlat(ix) = LATS(ix, VmaxIx(ix));
    qu(ix) = U(ix, VmaxIx(ix));
    qv(ix) = V(ix, VmaxIx(ix));
  end;
  % Add a scale vector
  qlon(end+1) = -79.3;
  qlat(end+1) = 24.9;
  qu(end+1) = 100;
  qv(end+1) = 0;
  quiver(qlon, qlat, qu, qv, 0.25);
  text(qlon(end), qlat(end)+0.05, '1 m/s');



  % Bar-graph zonal maxima of mapped property, below each map
  ah(2) = subplot('position', [0.10 0.05 0.88 0.16]);
  bar('v6', LONS(1,:), max(prp));
  set(gca, 'xlim', [-80.50 -79.10]);
  if ( ~isempty(prplim) )
    set(gca, 'ylim', prplim);
  end;
  ylabel(sprintf('%s', ylbl));
  set(ah(2), 'XTickLabel', []);




  %%%% 
  %%%% DEBUG
  %%%% 
  if ( ~isempty(prppeaks) )
    figure;
    subplot('position', [0.10 0.24 0.88 0.67]);
    contour(LONS, LATS, abs(prp), [prppeaks prppeaks], 'r');
    set(gca, 'xlim', [-80.50 -79.10]);
    set(gca, 'ylim', [+24.80 +25.90]);
    set(gcf, 'Position', position);
  end;







  % Initial template size in # of gridpoints
  initialtemplatesize = 2;

  curlonix = 0;

  latix = 1;
  while ( latix <= size(owp,1) )

    lonix = curlonix + 1;
    while ( lonix <= size(owp,2) )

      tsz = initialtemplatesize;
      if ( abs(mean(mean(owp(latix:latix+tsz-1, lonix:lonix+tsz-1)))) >= cutoff )
        while (abs(mean(mean(owp(latix:latix+tsz-1, lonix:lonix+tsz-1)))) >= (cutoff/2) )
          tsz = tsz + 1;
        end;
        curlonix = lonix + tsz;
      end;

    end;
  end;









  endpts(1,1,:) = pts(1,2:end-1);
  endpts(1,2,:) = pts(2,2:end-1);
  endpts(2,1,:) = pts(1,2:end-1) + (Vmult(2:end-1).*nrm(1,:));
  endpts(2,2,:) = pts(2,2:end-1) + (Vmult(2:end-1).*nrm(2,:));



  figure;
  hold on;
  plot(lon, lat, 'rs');
  line(squeeze(endpts(:,1,:)), squeeze(endpts(:,2,:)));
  xlim([LON(1) LON(end)]); ylim([LAT(1) LAT(end)]);
  return;

  % Defaults: nautical miles, and degrees CCW from East. WTF, CSIRO??
  ix = 0;
  [rngkm, brgE] = sw_dist(endpts(:,2,ix), endpts(:,1,ix), 'km');
  rng = rngkm .* 1e3;
  brg(brgE <= 90) = 90 - brgE(brgE <= 90);
  brg(brg > 90) = 450 - brgE(brg > 90);

  % Defaults: nautical miles, and degrees CCW from East. WTF, CSIRO??
  [rngkm, brgE] = sw_dist(endpts(:,2,ix), endpts(:,1,ix), 'km');
  rng = rngkm .* 1e3;

  x = rng/hypot((lon2-lon1), (lat2-lat1));
  vs = maxV * sin(2*pi*(lambdakm/111));

  for ix = 1:numel(U)

    [rng, brg] = dist([lat LAT(ix)], [lon LON(ix)]);

    if ( 1e3 <= rng && rng <= radiusm )
      switch(vprof)
       case { 'solid' },
        % Solid-body rotation
        Vr = maxVr*(rng/radiusm);
       case { 'sin', 'sine' },
        % Maximum velocity at 1/2 radius - sine-bell profile
        Vr = maxVr*sin(pi*(rng/radiusm));
       otherwise,
        % Maximum velocity at 1/2 radius - sharp sech profile
        x = 10 * ( (rng/radiusm) - 0.5 );
        Vr = maxVr*sech(x);
      end;

      U(ix) = U(ix) + cosd(brg)*Vr;
      V(ix) = V(ix) - sind(brg)*Vr;

    end;

  end;








  lonW = fix(LON(1) - 0.01);
  lonE = -fix(-LON(1) - 0.01);
  latS = fix(LAT(1) - 0.01);
  latN = -fix(-LAT(1) - 0.01);





  for ix = 1:size(endpts,3)
    line(endpts(:,1,ix), endpts(:,2,ix));
  end;



  if ( numel(lon) == 2 )
    lonix = find( (abs(LON(1,:) - lon(1)) < resX), 1 );
    latix = find( (abs(LAT(:,1) - lat(1)) < resY), 1 );
    lon1 = LON(1, lonix);
    lat1 = LAT(latix, 1);
    lonix = find( (abs(LON(1,:) - lon(2)) < resX), 1 );
    latix = find( (abs(LAT(:,1) - lat(2)) < resY), 1 );
    lon2 = LON(1, lonix);
    lat2 = LAT(latix, 1);
  end;









minlon = min(lonl);
maxlon = max(lonl);
minlat = min(latl);
maxlat = max(latl);

minlonix = find(abs(lon - minlon) <= 0.01, 1);
maxlonix = find(abs(lon - maxlon) <= 0.01, 1);
minlatix = find(abs(lat - minlat) <= 0.01, 1);
maxlatix = find(abs(lat - maxlat) <= 0.01, 1);

lons = interp1(lonl, lonl, lon(minlonix:maxlonix));
lats = interp1(latl, latl, lat(minlatix:maxlatix));

% lons = lons([1 find(diff(lons) | diff(lats))]);
% lats = lats([1 find(diff(lons) | diff(lats))]);









lonsix = 0;
latsix = 0;
lons = [];
lats = [];
for lix = 1:length(lonl)
  newlonsix = find((abs(lon - lonl(lix)) <= 0.01), 1);
  newlatsix = find((abs(lat - latl(lix)) <= 0.01), 1);
  if ( newlonsix ~= lonsix || newlatsix ~= latsix )
    lonsix = newlonsix;
    latsix = newlatsix;
    lons(end+1) = lon(lonsix);
    lats(end+1) = lat(latsix);
  end;
end;





lonsix = 0;
latsix = 0;
lons = [];
lats = [];
for lix = 1:length(lonl)
  newlonsix = find((abs(lon - lonl(lix)) <= 0.01), 1);
  newlatsix = find((abs(lat - latl(lix)) <= 0.01), 1);
  if ( newlonsix ~= lonsix || newlatsix ~= latsix )
    lonsix = newlonsix;
    latsix = newlatsix;
    lons(end+1) = lon(lonsix);
    lats(end+1) = lat(latsix);
  end;
end;







lonix = round(1 + ( length(lon) * ( (max(lon) - lonl) / (max(lon) - min(lon)) ) ));
latix = round(1 + ( length(lat) * ( (max(lat) - latl) / (max(lat) - min(lat)) ) ));

lons = lon(lonix);
lats = lat(latix);





ds=0:0.02:1;
lonl = ((1-ds) * lonp(1)) + (ds * lonp(2));
latl = ((1-ds) * latp(1)) + (ds * latp(2));

R = fix(1000/length(lonl));
lons = interp(lonl, R);
lats = interp(latl, R);










  if ( numel(lon) == 2 )
    lon1 = LON(abs(LON - lon(1)) < resX);
    lon2 = LON(abs(LON - lon(2)) < resX);
    lat1 = LAT(abs(LAT - lat(1)) < resY);
    lat2 = LAT(abs(LAT - lat(2)) < resY);

    R = 1000/length(lons

    x = 
    lon
  end;










1;

lon = -80:0.02:-79;
lat = 24:0.02:25;
[LON, LAT] = meshgrid(lon, lat);

lonp=[-79.2 -79.3];
latp=[24.5 24.1];

ds=0:0.001:1;
lonl = ((1-ds) * lonp(1)) + (ds * lonp(2));
latl = ((1-ds) * latp(1)) + (ds * latp(2));


figure;
plot(lonl,  latl, 'b.');

lons = interp1(latl, lonl, LAT(:,1));
lats = interp1(lonl, latl, LON(1,:));

figure;
plot(lons,  lats, 'r.');




lonp=[3.2 5.3];
latp=[5.1 2.5];

ds=0:0.001:1;
lonl = ((1-ds) * lonp(1)) + (ds * lonp(2));
latl = ((1-ds) * latp(1)) + (ds * latp(2));

lwvgd = fit(lonp', latp', 'pchipinterp');
lons = feval(lwvgd, lon);
twvgd = fit(latp', lonp', 'pchipinterp');
lats = feval(twvgd, lat);





function [U, V] = meanflow(LON, LAT, u0, v0)

  lons = 1:10;
  lats = lons;

  [LON, LAT] = meshgrid(lons, lats);

  U = repmat(u0, size(LON));
  V = repmat(v0, size(LON));

return;




    if ( 1e3 <= rng && rng <= radiusm )
      switch(vprof)
       case { 'solid' },
        % Solid-body rotation
        U(ix) = U(ix) + cosd(brg)*maxVr*(rng/radiusm);
        V(ix) = V(ix) - sind(brg)*maxVr*(rng/radiusm);
       case { 'sin' },
        % Maximum velocity at 1/2 radius - sine-bell profile
        U(ix) = U(ix) + cosd(brg)*maxVr*sin((rng*pi)/radiusm);
        V(ix) = V(ix) - sind(brg)*maxVr*sin((rng*pi)/radiusm);
       otherwise,
        % Maximum velocity at 1/2 radius - double-sech profile
        U(ix) = U(ix) + cosd(brg)*maxVr*sech((4*rng*pi)/radiusm);
        V(ix) = V(ix) - sind(brg)*maxVr*sech((4*rng*pi)/radiusm);
      end;

      U(ix) = U(ix) + cosd(brg)*Vr;
      V(ix) = V(ix) - sind(brg)*Vr;

    end;










  % KISS
  rad = radiuskm / 111;

  % Grid points per radius
  pts = fix(rad/res);

  lons = linspace(lon-rad, lon+rad, (2*pts)-1);
  lats = linspace(lat-rad, lat+rad, (2*pts)-1);

  [LON, LAT] = meshgrid(lons, lats);







function [U, V, LON, LAT] = eddy(vort, solid, lon, lat, radiuskm)

  % KISS
  rad = radiuskm / 111;

  % Grid points per radius
  pts = fix(rad/res);

  lons = linspace(lon-rad, lon+rad, (2*pts)-1);
  lats = linspace(lat-rad, lat+rad, (2*pts)-1);

  [LON, LAT] = meshgrid(lons, lats);

  U = repmat(nan, size(LON));
  V = repmat(nan, size(LON));
  for lon = 1:length(lons)
    for lat = 1:length(lats)
      [rng, brg] = dist([0 lats(lat)], [0 lons(lon)]);
      if ( rng <= radiuskm )
        if ( solid )
          % Solid-body rotation
          U(lon, lat) = sind(brg)*vort*(rng/radiuskm);
          V(lon, lat) = cosd(brg)*vort*(rng/radiuskm);
        else
          % Maximum velocity near middle
          U(lon, lat) = sind(brg)*vort*sin((rng*pi)/radiuskm);
          V(lon, lat) = cosd(brg)*vort*sin((rng*pi)/radiuskm);
        end;
    end;
  end;


  circpts = fix((2*pi*rad)/res);
  ang = linspace(0, (2*pi), circpts);

  lons = linspace(-rad, +rad, (2*pts)-1);
  lats = linspace(-rad, +rad, (2*pts)-1);
  [LON, LAT] = meshgrid(lons, lats);


  U = (u .* sin(ang

  lons = lons + lon;
  lats = lats + lat;

return;












  % Solid-body rotation
  if ( solid )
    u = linspace(0, vort, pts);
    v = -u;
  % Maximum velocity near middle
  else
    u = vort .* sin(pi ./ pts);
    v = -u;
  end;





  circpts = fix((2*pi*rad)/res);
  ang = linspace(0, (2*pi), circpts);

  lons = linspace(-rad, +rad, (2*pts)-1);
  lats = linspace(-rad, +rad, (2*pts)-1);
  [LON, LAT] = meshgrid(lons, lats);


  U = (u .* sin(ang

  lons = lons + lon;
  lats = lats + lat;
  [LON, LAT] = meshgrid(lons, lats);






  if ( plotzetar )
    fh = figure;
    ah(1) = subplot('position', [0.10 0.24 0.88 0.67]);
    hold on;

    [Vmax, VmaxIx] = max(SPD, [], 2);
    % Stupid f'in MATLAB; stupid f'in FORTRAN
    for ix = 1:length(VmaxIx)
      qlon(ix) = LONS(ix, VmaxIx(ix));
      qlat(ix) = LATS(ix, VmaxIx(ix));
      qu(ix) = U(ix, VmaxIx(ix));
      qv(ix) = V(ix, VmaxIx(ix));
    end;
    % Add a scale vector
    qlon(end+1) = -79.3;
    qlat(end+1) = 24.9;
    qu = 50;
    qv = 0;
    qgh = quiver(qlon, qlat, qu, qv);
    text(qlon(end), qlat(end)+0.05, '50cm/s');

    surfc(lons, lats, zetar);
    map_sofla(boundingbox);
    view(2);
    colorbar('East');
    set(gca, 'xlim', [-80.50 -79.10]);
    set(gca, 'ylim', [+24.80 +25.90]);
    set(gca, 'clim', [-350 +350]);

    title(sprintf('%s - \\nabla\\timesU', hr));

    % Plot meridional maximum of mapped property, below each map
    ah(2) = subplot('position', [0.10 0.05 0.88 0.16]);
    bar('v6', lons, max(zetar));
    set(gca, 'xlim', [-80.50 -79.10]);
    set(gca, 'ylim', [-350 +350]);
    ylabel('\zeta_m_a_x');
    set(ah(2), 'XTickLabel', []);

    set(gcf, 'Position', position);
  end;


  if ( plotdivr )
    fh = figure;
    ah(1) = subplot('position', [0.10 0.24 0.88 0.67]);
    hold on;
    surfc(lons, lats, divr);
    map_sofla(boundingbox);
    view(2);
    colorbar('East');
    set(gca, 'xlim', [-80.50 -79.10]);
    set(gca, 'ylim', [+24.80 +25.90]);
    set(gca, 'clim', [-500 +1000]);
    title(sprintf('%s - \\nabla^.U', hr));

    % Plot meridional maximum of divergence below each map
    ah(2) = subplot('position', [0.10 0.05 0.88 0.16]);
    bar('v6', lons, max(divr));
    set(gca, 'xlim', [-80.50 -79.10]);
    set(gca, 'ylim', [-500 +1000]);
    ylabel('DIV_m_a_x');
    set(ah(2), 'XTickLabel', []);

    set(gcf, 'Position', position);
  end;

  if ( plotowp )
    fh = figure;
    ah(1) = subplot('position', [0.10 0.24 0.88 0.67]);
    hold on;
    surfc(lons, lats, owp);
    map_sofla(boundingbox);
    view(2);
    colorbar('East');
    set(gca, 'xlim', [-80.50 -79.10]);
    set(gca, 'ylim', [+24.80 +25.90]);
    set(gca, 'clim', [0e5 +5e5]);
    title(sprintf('%s - Okubo-Weiss', hr));

    % Plot meridional maximum of Okubo-Weiss below each map
    ah(2) = subplot('position', [0.10 0.05 0.88 0.16]);
    bar('v6', lons, max(owp));
    set(gca, 'xlim', [-80.50 -79.10]);
    set(gca, 'ylim', [-5e5 +10e5]);
    ylabel('OWP_m_a_x');
    set(ah(2), 'XTickLabel', []);

    set(gcf, 'Position', position);
  end;







  isobath = cs(:,(-81.0 < cs(1,:) & cs(1,:) < -79.5));




  isobath = repmat(nan, [2 size(sofla_topo_lats,1)]);
  for ltix = 1:size(sofla_topo_lats,1)
    zix = find(abs(sofla_topo(ltix,:) - isodepth) <= 10, 1);
    if ( ~isempty(zix) && sofla_topo_lons(1,zix) < -79.5 )
      isobath(1, ltix) = sofla_topo_lons(1,zix);
      isobath(2, ltix) = sofla_topo_lats(ltix,1);
    end;
  end;






    % [cs, hs] = contour(lg, lt, sofla_topo_low, [isodepth isodepth], ...
    %                   'LineColor', 'red');

    % sofla_topo = interp2(lg, lt, sofla_topo_low, lgi, lti, '*cubic');




  line(isobath(1,:), isobath(2,:), 'color', 'blue');







  sct = isobath(:,3:end) - isobath(:,1:end-2);
  nrm = -[sct(2,:) ; sct(1,:)];
  transects = repmat(nan, [2 2 length(isobath(1,1:end-2))]);
  transects(1,:,:) = [ (isobath(1,2:end-1) - (5.*nrm(1,:))) ; (isobath(2,2:end-1) - (5.*nrm(2,:))) ];
  transects(2,:,:) = [ (isobath(1,2:end-1) + (5.*nrm(1,:))) ; (isobath(2,2:end-1) + (5.*nrm(2,:))) ];






  [lg,lt] = meshgrid(lgs,lts);
  [cs, hs] = contour(lgs, lts, sofla_topo, [-300 -300], ...
                     'LineColor', [.6 .5 .4]);
  clabel(cs, hs);




  line(sofla_coast(:,1), sofla_coast(:,2), 'color', [.5 .4 .3]);





%     set(ah(1), 'position', [0.10 0.27 0.88 0.64]);



        % Sync time (X) axes, make subplots pretty together
        axes(ah(1)); datetick;
        axes(ah(2)); datetick;
        lims([1 2], 1) = xlim(ah(1)); lims([1 2], 2) = xlim(ah(2));
        minlim = min(lims(1,:)); maxlim = max(lims(2,:));
        xlim(ah(1), [minlim maxlim]); xlim(ah(2), [minlim maxlim]);






  if ( plotzetar )
    fh = figure;
    hold on;
    map_sofla(boundingbox);
    zetarex = zetar;
    % zetarex( -100 < zetarex & zetarex < 100 ) = nan;
    %zetarex = exp(zetar/100);
    %zetarex( (1/2.718) < zetarex & zetarex < 2.718 ) = nan;
    % zetarex = del2(zetarex);
    surfc(lons, lats, zetarex);
    view(2);
    colorbar;
    set(gca, 'xlim', [-80.50 -79.10]);
    set(gca, 'ylim', [+24.80 +25.90]);
    set(gca, 'clim', [-350 +350]);
    title(sprintf('%s - \\nabla\\timesU', hr));
    % title(sprintf('%s - \\nabla^2[\\nabla\\timesU]', hr));
    set(gcf, 'Position', [1 1 1280 726]);
  end;






    zetarex( LONS > -79.75 | (-100 < zetarex & zetarex < 100) ) = nan;






  persistent sofla_topo;


  if ( ~exist('sofla_topo', 'var') || isempty(sofla_topo) )
    sofla_topo = load('topo');
  end;

if rlong>360, rlong=rlong-360; llong=llong-360; end;
if llong<-360, rlong=rlong+360; llong=llong+360; end;

lts=(blat:tlat);
lgs=(llong:rlong);

if rlong<0,
  topo=topo(lts+90.5,lgs+360.5);
elseif llong<0 & rlong>=0,
  topo=topo(lts+90.5,[(360.5+llong:end) (1:rlong+0.5)]);
else
  topo=topo(lts+90.5,lgs+.5);
end;


  [dpths, dpthlon, dpthlat] = m_elev(boundingbox);
  contour(dpthlon, dpthlat, dpths, [30 50 100 150 300 600]);









  m_proj('Mercator');
  %m_proj('Transverse Mercator', 'lon', [boundingbox(1:2)], ...
  %       'lat', [boundingbox(3:4)]);
  %m_grid;
  %m_line(sofla_coast(:,1), sofla_coast(:,2));





  dd = 5;
  isod = 030; iso030 = find(isod-dd < dpths & dpths < isod+dd);
  isod = 050; iso050 = find(isod-dd < dpths & dpths < isod+dd);
  isod = 100; iso100 = find(isod-dd < dpths & dpths < isod+dd);
  isod = 200; iso200 = find(isod-dd < dpths & dpths < isod+dd);
  isod = 300; iso300 = find(isod-dd < dpths & dpths < isod+dd);
  isod = 600; iso600 = find(isod-dd < dpths & dpths < isod+dd);
  line(dpthlon(iso30), dpthlat(iso30));





    latix = find( abs(abs(lats - lat(ix)) - (dlat/2)) <= eps );
    lonix = find( abs(abs(lons - lon(ix)) - (dlon/2)) <= eps );




m_proj('Transverse Mercator', 'lon', [-80.50 -79.10], 'lat', [+24.80 +25.90]);



  %zeta = repmat(nan, size(U));
  %for lt = 2:length(lats)-1
  % for ln = 2:length(lons)-1
  %   zeta(lt, ln) = (V(lt,ln+1)-V(lt,ln-1)) - (U(lt+1,ln)-U(lt-1,ln));
  % end;
  %end;




    figure;
    hold on;
    quiver(lon(goodix), lat(goodix), u(goodix), v(goodix), 0.5); %0);
    quiver(lon(~goodix), lat(~goodix), u(~goodix), v(~goodix), 0.5, 'Color', 'red');




    LAT = unique(lat);
    infl = repmat(nan, size(LAT));
    grd = repmat(nan, size(LAT));
    LON = repmat(nan, size(LAT)); 
    for ix = 1:length(LAT)
      myline = find(lat == LAT(ix));
      % Don't bother near lower and upper boundaries
      if ( length(myline) >= 20 )
        [grd(ix) lonix] = max(diff(v(myline(2:end-1))));
        LON(ix) = lon(myline(lonix(1)));
        plot(LON(ix), LAT(ix), 'r*');
        [ign lonix] = min(diff(v(myline(2:end-1))));
        LONZ = lon(myline(lonix(1)));
%         plot(LONZ, LAT(ix), 'b*');
      end;
    end;











    LAT = unique(lat);
    infl = repmat(nan, size(LAT));
    grd = repmat(nan, size(LAT));
    LON = repmat(nan, size(LAT)); 
    for ix = 1:length(LAT)
      myline = find(lat == LAT(ix));
      % Don't bother near lower and upper boundaries
      if ( length(myline) >= 20 )
        revs = find(diff(v(myline)) < -1);
        if (isempty(revs))
          infl(ix) = length(myline);
        else
          infl(ix) = revs(1);
          LONZ = lon(myline(revs));
          plot(LONZ, repmat(LAT(ix), size(LONZ)), 'r*');
        end;
%         myline = myline(1:infl(ix));
%         [grd(ix) lonix] = max(diff(v(myline)));
%         LON(ix) = lon(myline(lonix(1)));
      end;
    end;
%     plot(LON, LAT, 'r*');
