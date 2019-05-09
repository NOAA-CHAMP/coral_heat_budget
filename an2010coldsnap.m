1;

%dts=datenum(1999,1,[1,20]);
dts=datenum(2010,1,[5,25]);
f=@(x)(find(dts(1)<=x.date&x.date<=dts(2)));

t=subset_ts(stn.ndbc_sea_t,f);
ta=subset_ts(stn.ndbc_air_t,f);
qt=subset_ts(stn.simple_ndbc_erai_erai_30a_net_flux_term,f);
bqt=subset_ts(stn.b_ndbc_erai_erai_30a_net_flux_term,f);
dt=subset_ts(stn.b_ndbc_erai_erai_30a_avhrr_dt,f);
hc=subset_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt,f);
fmg;
h=plot(t.date,t.data,...
    ta.date,ta.data,...
    qt.date,t.data(1)+cumsum(qt.data),...
    bqt.date,t.data(1)+cumsum(bqt.data),...
    dt.date,t.data(1)+cumsum(dt.data),...
    hc.date,t.data(1)+cumsum(hc.data));
if ( numel(h)>=6 )
    set(h(3:6),'LineW',3);
end;
datetick3;
ylim([12,28]);
legend('Ts','Ta','Q0','Q0+Qb','Q0+Qb+u^.\nablaT+K\nabla^2T','\partial_tT','Location','SouthEast');

fmg;
ax=plotyy(stn.ndbc_sea_t.date,stn.ndbc_sea_t.data,...
    stn.absorbed_erai_ndbc_srf.date,stn.absorbed_erai_ndbc_srf.data);
axes(ax(1)); hold on; axes(ax(2)); hold on;
%linkprop(ax,{'xlim','xticklabel'});
%%xlim(datenum(1999,1,[4,10])); datetick('x',2,'keeplimits');
%xlim([dts(1)+3,dts(2)-10]); datetick('x',2,'keeplimits');

axes(ax(1)); hold on;
xlim([dts(1)+3,dts(2)-10]); datetick('x',2,'keeplimits');
ylim([12,28]);
plot(stn.ndbc_air_t.date,stn.ndbc_air_t.data,'m');
legend('Ts','Ta','Location','NorthWest');

axes(ax(2)); hold on;
xlim([dts(1)+3,dts(2)-10]); datetick('x',2,'keeplimits');
ylim([0,700]);
plot(stn.erai_ndbc_srf.date,stn.erai_ndbc_srf.data,'r','linew',2);
plot(stn.simple_ndbc_erai_erai_30a_net_flux.date,stn.simple_ndbc_erai_erai_30a_net_flux.data,'g','linew',3);
%plot(stn.erai_ndbc_arf.date,stn.erai_ndbc_arf.data,'y','linew',2);
legend('\gammaQSW','QSW','Q0','Location','NorthEast');
