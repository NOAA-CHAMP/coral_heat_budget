1;

disspath = get_thesis_path('../DISS');

pnl=0;
xl=[1e-1,2500];
yl=[9e-6,1.5e6];
for cst={'fwyf1','mlrf1','smkf1','lonf1','looe1'};
    x=get_station_from_station_name(cst{:});
    if(strcmpi(cst{:},'looe1'))
        x=get_looe1_microcat(x);
        [P,W,fh,lh]=plot_spec(x,'microcat_seatemp',[],[],xl,yl,[],true);
        set(lh,'LineWidth',1.5,'Color','k');
        text(2e-1,1e5,['(',char('a'+pnl),')'],'FontSize',24); pnl=pnl+1;
        print('-dtiff','-r300',fullfile(disspath,[mfilename,'-',cst{:},'_5m-plot_spec-seatemp.tif']));

        x=get_looe1_adcp(x);
        [P,W,fh,lh]=plot_spec(x,'adcp_seatemp',[],[],xl,yl,[],true);
    else
        x=load_all_ndbc_data(x);
        [P,W,fh,lh]=plot_spec(x,'ndbc_sea_t',[],[],xl,yl,[],true);
    end;
    set(lh,'LineWidth',1.5,'Color','k');
    text(2e-1,1e5,['(',char('a'+pnl),')'],'FontSize',24); pnl=pnl+1;
    print('-dtiff','-r300',fullfile(disspath,[mfilename,'-',cst{:},'-plot_spec-seatemp.tif']));
    x=[];
    clear x P W fh lh;
end;

clear cst xl disspath
