1;

if ( ~exist('x','var') || ~isfield(x,'ndbc_sea_t') )
  x=[]; clear x;
  x = get_station_from_station_name('smkf1');
  x = load_all_ndbc_data(x);
end;

flds = {'ndbc_sea_t',...
        'ndbc_sea_t_1_d_dev_3_d_avg','ndbc_sea_t_3_d_avg',...
        'ndbc_sea_t_3_d_avg_0_d_asof_diff_ndbc_sea_t_1_d_min',};

x = verify_variable(x,flds);

[ig,hl,ax,fh]=...
    multiplot_station(x,flds,...
                      [upper(x.station_name),': metrics derived from sea temperature'],...
                      [],{'T_s','\mu_3_d\sigma_1_dT_s','\mu_3_dT_s','anom^m^i^n_3_dT_s'},...
                      x.ndbc_sea_t.date([1,end]),[],false,{'k-','r-','g-','b-'});
datetick('x',10,'keeplimits','keepticks');
set(ax(end),'xtick',datenum([1990,1993,1996,1999,2002,2005],1,1))
set(hl{2},'Color',[0,0.5,0]);
set([hl{:}],'LineWidth',1.5);
set(ax,'FontSize',14);
xlh=get(ax,'YLabel'); set([xlh{:}],'FontSize',16);
print('-dtiff',fullfile(get_thesis_path('../DISS'),[lower(stn.station_name),'-',mfilename,'.tif']));

pause;

[ig,hl,ax,fh]=...
    multiplot_station(x,flds,...
                      [upper(x.station_name),': closeup of variability metrics - Summer 2001'],...
                      [],{'T_s','\mu_3_d\sigma_1_dT_s','\mu_3_dT_s','anom^m^i^n_3_dT_s'},...
                      datenum(2001,[3,11],[15,30]),[],false,{'k-','r-','g-','b-'});
datetick('x',3,'keeplimits');
set(hl{2},'Color',[0,0.5,0]);
set([hl{:}],'LineWidth',1.5);
set(ax,'FontSize',14);
xlh=get(ax,'YLabel'); set([xlh{:}],'FontSize',16);
print('-dtiff',fullfile(get_thesis_path('../DISS'),[lower(stn.station_name),'-',mfilename,'-Summer-2001.tif']));
