1;

dts=datenum(2002,1,[1,16]);

if ( ~exist('ncs','var') || ~isfield(ncs,'tl2') )
  ncs=[]; clear ncs;
  ncs = get_ncore_2000;
end;

sfld='ndbc_sea_t';
bfld='b_ndbc_erai_erai_30a_avhrr_dt';
%% EXPERIMENT with NCORE TL temperatur
%sfld='tl_seatemp_hourly';
%bfld='b_tlseatemphourly_ndbc_erai_erai_30a_avhrr_dt';

if ( ~exist('stn','var') || ~isfield(stn,bfld) )
  error('No field "%s": Run a heat budget in STN first!',bfld);
  %% E.G., RUN WITH MINIMAL ASSUMPTIONS RE: ADVECTION OR LIGHT
  %stn = optimize_station_heat_budget('mlrf1','erai','avhrr_weekly','ndbc','tpxo_tide','erai',{'Q0_LOWPASS','','QE_LOWPASS','_72_h_lp','ufld','adcp_sfc_u','vfld','adcp_sfc_v'},struct('force_bic_insolation',true,'tidal_advection',false,'ocean_current_wind_multiplier',0,'force_adcp_currents',true,'add_alongshore_advection',true,'albedo',0));
end;
bffld=[bfld,'_flux'];

stn = station_reorient_vectors(stn,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');
stn = station_reorient_vectors(stn,'isobath_orientation','tpxo_tide_u','tpxo_tide_v');
stn = station_reorient_vectors(stn,'isobath_orientation','tmd_tide_u','tmd_tide_v');

[ig,t0ix]=min(abs(dts(1)-stn.(sfld).date));
t0=stn.(sfld).data(t0ix);
stn = station_heat_flux_term(stn,bffld,'shallow_term',sfld,stn.opts.default_salinity,stn.depth);
mean_depth = nanmean(stn.mean_tpxo_tide_i_depth.data);
stn = station_heat_flux_term(stn,bffld,'mean_term',sfld,stn.opts.default_salinity,mean_depth);
max_depth = 1.25e3*stn.ngdc_offshore_slope;
stn = station_heat_flux_term(stn,bffld,'deep_term',sfld,stn.opts.default_salinity,max_depth);
%qs=subset_ts(stn.b_ndbc_erai_erai_30a_avhrr_dt,@(x)(find(dts(1)<=x.date&x.date<=dts(2))));
qs=subset_ts(stn.shallow_term,@(x)(find(dts(1)<=x.date&x.date<=dts(2))));
qm=subset_ts(stn.mean_term,@(x)(find(dts(1)<=x.date&x.date<=dts(2))));
qd=subset_ts(stn.deep_term,@(x)(find(dts(1)<=x.date&x.date<=dts(2))));


fmg;
spt(3,1,1:2);
hold on;
plot_ts(ncs.cdp.cm_seatemp_3hlp,ncs.csh.cm_seatemp_3hlp,ncs.tl2.tl_seatemp_3hlp,ncs.tl3.tl_seatemp_3hlp,ncs.tl4.tl_seatemp_3hlp);
plot_ts(stn.(sfld),'k-','linew',3);
plot(qs.date,t0+cumsum(qs.data),'k:','linew',1);
plot(qm.date,t0+cumsum(qm.data),'k:','linew',2);
plot(qd.date,t0+cumsum(qd.data),'k-','linew',2,'Color',[.5,.5,.5]);
%%plot_ts(ts_op(ts_op(stn.tpxo_tide_xshore,10,'*'),21,'+'),'k:','linew',1.5);
%%plot_ts(ts_op(ts_op(ncs.csh.cm_xshore_3hlp,10,'/'),21,'+'),'m-','linew',1.5);
%xlim(dts); datetick3;
xlim(dts); datetick('x',6,'keeplimits');
set(gca,'xticklabel',[])
grid on;
lh=legend('NC reef crest 21m',...
          'NC reef crest 4m',...
          'NC outer shelf 7m',...
          'NC mid-shelf 4m',...
          'NC channel 5.6m',[upper(stn.station_name),' T_s ',num2str(stn.depth),'m'],...
          ['\partial_tT ',num2str(stn.depth,'%.1f'),'m'],...
          ['\partial_tT ',num2str(mean_depth,'%.1f'),'m'],...
          ['\partial_tT ',num2str(max_depth,'%.1f'),'m'],...
          'Location','SouthEast');
set(lh,'FontSize',8);
titlename(['Cold-snap Jan 2002: ',upper(stn.station_name),' heat budget and NCORE data']);
%ylim([0,28]);
ylim([16,25]);
ylabel('T_S [^oC]');

spt(3,1,3);
hold on;
plot_ts(stn.tpxo_tide_xshore,'k-','linew',1.5,'color',[.5,.5,.5]);
plot_ts(ts_op(ncs.csh.cm_xshore_3hlp,100,'/'),'r-','linew',1);
plot_ts(ts_op(ncs.cdp.cm_xshore_3hlp,100,'/'),'b-','linew',1);
plot_ts(ts_op(ncs.b.cm_xshore_3hlp,100,'/'),'g-','color',[0,.5,0],'linew',1);
%xlim(dts); datetick3;
xlim(dts); datetick('x',6,'keeplimits');
text(dts(1)+1,+0.25,'Offshore');
text(dts(1)+1,-0.25,'Onshore');
grid on;
lh=legend('TPXO tide model\bullet\nablah',...
          'NC crest 4m\bullet\nablah',...
          'NC crest 21m\bullet\nablah',...
          'NC channel 7m\bullet\nablah',...
          'Location','SouthEast');
set(lh,'FontSize',8);
%ylim([-0.4,+0.4]);
ylim([-0.28,+0.28]);
ylabel('u [m\bullets^-^1]');

%annotline([],21)
%ylim([19,25]);

%legend('C shallow','TL2 T_s','TL3 T_s','TL4 T_s',[upper(stn.station_name),' T_s'], 'TPXO XS','C sh XS','C dp XS', 'Location','NorthEast');

print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-ncore-cold-snap-Jan-2002.tif']));
