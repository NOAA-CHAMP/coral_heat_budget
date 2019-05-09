1;


%% First do some analysis of MLRF1
if ( ~exist('mlrf1','var') )
  %mlrf1 = optimize_station_heat_budget('mlrf1','erai','avhrr_weekly','ndbc','tpxo_tide','erai');
  mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1);
  % MLRF1 Time-series SOM/PCA
  [mlrf1,sm,pc,sc,framedts,fldnms,fhs] = som_vs_pca_ts(mlrf1,[],14,[],[],[],{'ndbc_sea_t'});
  mlrf1.ndbc_sea_t_14_d_pc1.date = framedts(:,1);
  mlrf1.ndbc_sea_t_14_d_pc1.data = sc(1,:);
  mlrf1.ndbc_sea_t_14_d_pc2.date = framedts(:,1);
  mlrf1.ndbc_sea_t_14_d_pc2.data = sc(2,:);
  mlrf1.ndbc_sea_t_14_d_pc3.date = framedts(:,1);
  mlrf1.ndbc_sea_t_14_d_pc3.data = sc(3,:);
  mlrf1.ndbc_sea_t_14_d_bmu.date = framedts(:,1);
  mlrf1.ndbc_sea_t_14_d_bmu.data = sm.bmus;
  mlrf1.ndbc_sea_t_14_d_qerr.date = framedts(:,1);
  mlrf1.ndbc_sea_t_14_d_qerr.data = sm.qerrs;

  mlrf1 = filter_gaps(mlrf1,'ndbc_sea_t','ndbc_sea_t_14_d_pc1',14);
  mlrf1 = filter_gaps(mlrf1,'ndbc_sea_t','ndbc_sea_t_14_d_pc2',14);
  mlrf1 = filter_gaps(mlrf1,'ndbc_sea_t','ndbc_sea_t_14_d_pc3',14);
  mlrf1 = filter_gaps(mlrf1,'ndbc_sea_t','ndbc_sea_t_14_d_bmu',14);
  mlrf1 = filter_gaps(mlrf1,'ndbc_sea_t','ndbc_sea_t_14_d_qerr',14);
end;


%% Then analyze SMKF1 in more detail

if ( ~exist('stn','var') )
  stn = optimize_station_heat_budget('smkf1','erai','avhrr_weekly','ndbc','tpxo_tide','erai');
end;

if ( ~isfield(stn,'MASS_CORAL_STRESS') )
    stn = read_station_ef_csv(stn,'MASS-CORAL-STRESS');
    stn.ef_mcs = stn.MASS_CORAL_STRESS.sri;
end;
if (0)
  fmg; 
  subplot_tight(7,1,3); 
  grpplot_ts(stn.ef_mcs,@get_year,@nansum,[],'s','MarkerSize',10,'MarkerFaceColor','k');
end;


if (0)

  % SMKF1 Time-series SOM/PCA
  [stn,sm,pc,sc,framedts,fldnms,fhs] = som_vs_pca_ts(stn,[],14,[],[],[],{'ndbc_sea_t'});
  %close(fhs(1)); % Close SOM modes figure
  if (0)
    fmg; for ix=1:3; spt(3,1,ix); plot(framedts(:,1),sc(ix,:),'.'); datetick3; xlabel(num2str(ix)); end; linkaxes;
    ylim([-1,1]*7000);
  end;
  stn.ndbc_sea_t_14_d_pc1.date = framedts(:,1);
  stn.ndbc_sea_t_14_d_pc1.data = sc(1,:);
  stn.ndbc_sea_t_14_d_pc2.date = framedts(:,1);
  % HOLEY HACKOLA BATMAN!
  stn.ndbc_sea_t_14_d_pc2.data = sc(3,:).*2;

  stn.ndbc_sea_t_14_d_bmu.date = framedts(:,1);
  stn.ndbc_sea_t_14_d_bmu.data = sm.bmus;
  stn.ndbc_sea_t_14_d_qerr.date = framedts(:,1);
  stn.ndbc_sea_t_14_d_qerr.data = sm.qerrs;


  % SMKF1 Extended SOM/PCA
  [stn,sm,pc,sc,framedts,fldnms,fhs] = som_vs_pca_ts(stn);
  %close(fhs);
  if (0)
    fmg; for ix=1:3; spt(3,1,ix); plot(framedts(:,1),sc(ix,:),'.'); datetick3; xlabel(num2str(ix)); end; linkaxes;
    ylim([-1,1]*5e4);
  end;
  stn.ndbc_sea_t_variability_14_d_pc1.date = framedts(:,1);
  stn.ndbc_sea_t_variability_14_d_pc1.data = sc(1,:)./10;
  stn.ndbc_sea_t_variability_14_d_pc2.date = framedts(:,1);
  stn.ndbc_sea_t_variability_14_d_pc2.data = sc(2,:)./10;
  stn.ndbc_sea_t_variability_14_d_pc3.date = framedts(:,1);
  stn.ndbc_sea_t_variability_14_d_pc3.data = sc(3,:)./10;

  stn.ndbc_sea_t_variability_14_d_bmu.date = framedts(:,1);
  stn.ndbc_sea_t_variability_14_d_bmu.data = sm.bmus;
  stn.ndbc_sea_t_variability_14_d_qerr.date = framedts(:,1);
  stn.ndbc_sea_t_variability_14_d_qerr.data = sm.qerrs;


  if (0)
    % SMKF1 Extended Heat-Budget SOM/PCA
    [stn,sm,pc,sc,framedts,fldnms,fhs] = som_vs_pca_ts(stn,[],1,[],[],[],{'ndbc_erai_erai_30a_avhrr_hc_dTdthc','fq_erai_avhrr_advected_heat','erai_ndbc_arf_term','b_erai_qbo_term','ndbc_erai_erai_30a_latent_flux_term','ndbc_sea_t_diff'});
  end;


  stn = filter_gaps(stn,'ndbc_sea_t','ef_mcs',14,[],14);

  stn = filter_gaps(stn,'ndbc_sea_t','ndbc_sea_t_14_d_pc1',14);
  stn = filter_gaps(stn,'ndbc_sea_t','ndbc_sea_t_14_d_pc2',14);
  stn = filter_gaps(stn,'ndbc_sea_t','ndbc_sea_t_14_d_bmu',14);
  stn = filter_gaps(stn,'ndbc_sea_t','ndbc_sea_t_14_d_qerr',14);

  stn = filter_gaps(stn,'ndbc_sea_t','ndbc_sea_t_variability_14_d_pc1',14);
  stn = filter_gaps(stn,'ndbc_sea_t','ndbc_sea_t_variability_14_d_pc2',14);
  stn = filter_gaps(stn,'ndbc_sea_t','ndbc_sea_t_variability_14_d_pc3',14);
  stn = filter_gaps(stn,'ndbc_sea_t','ndbc_sea_t_variability_14_d_bmu',14);
  stn = filter_gaps(stn,'ndbc_sea_t','ndbc_sea_t_variability_14_d_qerr',14);

end;


[hl,ax,fh] = multiplot_ts([upper(stn.station_name),' sea temperature anomaly analyses'],...
                          stn,'ndbc_sea_t','ndbc_sea_t_1_day_deviation_3_day_average',...
                          'ef_mcs','ndbc_sea_t_14_d_pc1','ndbc_sea_t_14_d_pc2',...
                          'ndbc_sea_t_variability_14_d_pc1','ndbc_sea_t_variability_14_d_pc2',...
                          'ndbc_sea_t_variability_14_d_bmu',...
                          mlrf1,'ndbc_sea_t_14_d_pc1','ndbc_sea_t_14_d_pc2');
hl = [hl{end:-1:1}];
set(hl,'LineStyle','none','Color','k');
%set(hl,'LineStyle','none','Color','k','MarkerSize',8);
ylabel(ax(1),'T_s');		ylim(ax(1),[16,34]);		xlim(ax(1),datenum([1995,2011],1,1));
ylabel(ax(2),'\Theta');		ylim(ax(2),[ 0,1.5]);		xlim(ax(2),datenum([1995,2011],1,1));
ylabel(ax(3),'S/RI');		ylim(ax(3),[ 0,85]);		xlim(ax(3),datenum([1995,2011],1,1));
ylabel(ax(4),'T_s PC1');	ylim(ax(4),[-1,+1]*6500);	xlim(ax(4),datenum([1995,2011],1,1));
ylabel(ax(5),'T_s PC2');	ylim(ax(5),[-1,+1]*6500);	xlim(ax(5),datenum([1995,2011],1,1));
ylabel(ax(6),'\Theta PC1');	ylim(ax(6),[-1,+1]*6500);	xlim(ax(6),datenum([1995,2011],1,1));
ylabel(ax(7),'\Theta PC2');	ylim(ax(7),[-1,+1]*6500);	xlim(ax(7),datenum([1995,2011],1,1));
ylabel(ax(8),'\Theta BMU');	ylim(ax(8),[ 0,10]);		xlim(ax(8),datenum([1995,2011],1,1));
ylabel(ax(9),'ML T_s PC1');	ylim(ax(9),[-1,+1]*6500);	xlim(ax(9),datenum([1995,2011],1,1));
ylabel(ax(10),'ML T_s PC2');	ylim(ax(10),[-1,+1]*6500);	xlim(ax(10),datenum([1995,2011],1,1));

disp('ANNOTATE figure! Then type "dbcont" to print...');
keyboard;
if ( ishandle(fh) )
  print(fh,'-dtiff','-r300',fullfile(get_thesis_path('../DISS'),[mfilename,'-',lower(stn.station_name),'.tif']));
else
  disp('Did you dismiss the Figure before printing??');
end;
