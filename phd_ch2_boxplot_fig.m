1;

tic,

figspath = get_thesis_path('../figs');

if ( ~exist('s','var') || ~isfield(s,'mlrf1') )
s=[]; clear s;
for cstnm={'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1'};
%for cstnm={'smkf1'};
    stnm=cstnm{:};
    s.(stnm)=get_station_from_station_name(stnm);
    s.(stnm)=load_all_ndbc_data(s.(stnm));
    if (isfield(s.(stnm),'ndbc_dew_t'));
        s.(stnm)=station_dewp_to_relhumid(s.(stnm),'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
        s.(stnm)=station_relhumid_to_spechumid(s.(stnm),'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
    end;

    s.(stnm)=station_optimal_isobath_orientation(s.(stnm));
    s.(stnm)=verify_variable(s.(stnm),{'ndbc_wind1_u','ndbc_wind1_v'});
    s.(stnm)=station_reorient_vectors(s.(stnm),'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');
end;
end;

%fh = figure; hold on;
fh = fmg;
ix=1;
subplot_tight(4,3,ix); boxplot_ts(s.mlrf1.ndbc_sea_t,[],'allcolors','k'); ylim([6,34]); ylabel(''); grid on; ix=ix+1;
subplot_tight(4,3,ix); boxplot_ts(s.lonf1.ndbc_sea_t,[],'allcolors','k'); ylim([6,34]); ylabel(''); grid on; ix=ix+1;
subplot_tight(4,3,ix); boxplot_ts(s.dryf1.ndbc_sea_t,[],'allcolors','k'); ylim([6,34]); ylabel(''); grid on; ix=ix+1;

subplot_tight(4,3,ix); boxplot_ts(s.mlrf1.ndbc_air_t,[],'allcolors','k'); ylim([6,34]); ylabel(''); grid on; ix=ix+1;
subplot_tight(4,3,ix); boxplot_ts(s.lonf1.ndbc_air_t,[],'allcolors','k'); ylim([6,34]); ylabel(''); grid on; ix=ix+1;
subplot_tight(4,3,ix); boxplot_ts(s.dryf1.ndbc_air_t,[],'allcolors','k'); ylim([6,34]); ylabel(''); grid on; ix=ix+1;

subplot_tight(4,3,ix); boxplot_ts(s.fwyf1.ndbc_wind1_speed,[],'allcolors','k'); ylim([0,40]); ylabel(''); grid on; ix=ix+1;
subplot_tight(4,3,ix); boxplot_ts(s.fwyf1.ndbc_wind1_xshore,[],'allcolors','k'); ylim([-40,40]); ylabel(''); grid on; ix=ix+1;
subplot_tight(4,3,ix); boxplot_ts(s.fwyf1.ndbc_wind1_lshore,[],'allcolors','k'); ylim([-40,40]); ylabel(''); grid on; ix=ix+1;

subplot_tight(4,3,ix); boxplot_ts(s.lonf1.ndbc_wind1_speed,[],'allcolors','k'); ylim([0,40]); ylabel(''); grid on; ix=ix+1;
subplot_tight(4,3,ix); boxplot_ts(s.sanf1.ndbc_wind1_xshore,[],'allcolors','k'); ylim([-40,40]); ylabel(''); grid on; ix=ix+1;
subplot_tight(4,3,ix); boxplot_ts(s.sanf1.ndbc_wind1_lshore,[],'allcolors','k'); ylim([-40,40]); ylabel(''); grid on; ix=ix+1;

%subplots_set('FontSize',12);
%subplots_set('FontSize',8);
%subplots_set('Xlim',[0.5,12.5]);
subplots_set('XTick',[2:2:12],'XTickLabel',{'2','4','6','8','10','12'})
print('-dtiff',fullfile(figspath,['phd-ch2-fig-5-boxplots.tif']));


clear figspath cstnm stnm
toc,
