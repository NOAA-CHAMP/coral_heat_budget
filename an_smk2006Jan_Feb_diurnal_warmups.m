1;

xl = datenum(2006,[1,2],[20,7]);

if ( ~exist('mlrf1','var') )
    mlrf1 = load_all_ndbc_data([],'mlrf1');
end;
if ( ~exist('sanf1','var') )
    sanf1 = load_all_ndbc_data([],'sanf1');
end;


fmg; plot_ts(stn.ndbc_sea_t,stn.ndbc_air_t,mlrf1.ndbc_sea_t,mlrf1.ndbc_air_t); legend('SMK T_s','SMK T_a','MLR T_s','MLR T_a', 'Location','Best'); titlename('Temperatures');
xlim(xl); datetick3;

fmg; plot_ts(stn.ndbc_tide); annotline([],prctile(stn.ndbc_tide.data,1)); annotline([],prctile(stn.ndbc_tide.data,99)); titlename('Tide');
xlim(xl); datetick3;
fmg; fld='tpxo_tide_speed'; plot_ts(stn.(fld)); annotline([],prctile(stn.(fld).data(get_month(stn.(fld).date)==2),1)); annotline([],prctile(stn.(fld).data(get_month(stn.(fld).date)==2),99)); titlename([stn.station_name,'.',strrep(fld,'_','\_')]);
xlim(xl); datetick3;

fmg; fld='hourly_avhrr_weekly_sst_lshore'; plot_ts(stn.(fld)); annotline([],prctile(stn.(fld).data(get_month(stn.(fld).date)==2),1)); annotline([],prctile(stn.(fld).data(get_month(stn.(fld).date)==2),99)); titlename([stn.station_name,'.',strrep(fld,'_','\_')]);
xlim(xl); datetick3;
fmg; fld='hourly_avhrr_weekly_sst_xshore'; plot_ts(stn.(fld)); annotline([],prctile(stn.(fld).data(get_month(stn.(fld).date)==2),1)); annotline([],prctile(stn.(fld).data(get_month(stn.(fld).date)==2),99)); titlename([stn.station_name,'.',strrep(fld,'_','\_')]);
xlim(xl); datetick3;
fmg; fld='ww3_avhrr_advected_heat'; plot_ts(stn.(fld)); annotline([],prctile(stn.(fld).data(get_month(stn.(fld).date)==2),1)); annotline([],prctile(stn.(fld).data(get_month(stn.(fld).date)==2),99)); titlename([stn.station_name,'.',strrep(fld,'_','\_')]);
xlim(xl); datetick3;
ylim([-.5,.5]);


%fmg; plot_ts(stn.ndbc_sea_t,stn.ndbc_air_t); titlename([stn.station_name,'.ndbc\_sea\_t']);
fmg; plot_ts(stn.ndbc_sea_t,stn.ndbc_air_t,mlrf1.ndbc_sea_t,mlrf1.ndbc_air_t,sanf1.ndbc_sea_t,sanf1.ndbc_air_t); legend('SMK T_s','SMK T_a','MLR T_s','MLR T_a','SAN T_s','SAN T_a', 'Location','Best'); titlename('Temperatures: Jumps');
ax=gca;
jumpix = find(diff(stn.ndbc_sea_t.data) > 0.5)+1;
jumpix(find(diff(stn.ndbc_sea_t.date(jumpix))<16)+1) = [];
jumpix(get_year(stn.ndbc_sea_t.date(jumpix))<2000) = [];
for ix=jumpix(:)';
    if ( ~ishandle(ax) ); disp('Figure closed'); break; end;
    dt = stn.ndbc_sea_t.date(ix);
    vl = stn.ndbc_sea_t.data(ix);
    xlim(ax,[dt-10,dt+10]);
    ylim(ax,[vl-6,vl+3]);
    datetick3;
    pause;
end;
