1;

more off;

if ( ~exist('stn','var') || ~isfield(stn,'station_name') || ~strcmpi(stn.station_name,'LOOE1') )
  stn=[]; clear stn
  stn = anlooe1([],false);
  stn = verify_variable(stn,{'ndbc_wind1_u_40_h_lp','ndbc_wind1_v_40_h_lp'});
  stn = verify_variable(stn,{'ndbc_wind1_u_3_d_lp','ndbc_wind1_v_3_d_lp'});
  stn = verify_variable(stn,{'adcp_sfc_u_40_h_lp','adcp_sfc_v_40_h_lp'});

  plot_spec(stn,'ndbc_wind1_u');
  plot_spec(stn,'ndbc_wind1_v');
  plot_spec(stn,'adcp_sfc_u');
  plot_spec(stn,'adcp_sfc_v');
  plot_spec(stn,'adcp_btm_u');
  plot_spec(stn,'adcp_btm_v');
end;

%{
[rgs,ixes]=find_date_ranges(stn.adcp_sfc_u.date,(18/24)); datestr(rgs)
ix=8670:28153;
ix=8670:30096;
%}

rawu=stn.adcp_sfc_u; rawu.date(~isfinite(rawu.data))=[]; rawu.data(~isfinite(rawu.data))=[];
rawv=stn.adcp_sfc_v; rawv.date(~isfinite(rawv.data))=[]; rawv.data(~isfinite(rawv.data))=[];
%find_date_ranges(rawu.date,(18/24));
ix=[7036:18700,18702:25850];
%find_date_ranges(rawu.date(ix),(1.1/24));

dts=rawu.date(ix(1)):(1/24):rawu.date(ix(end));
u=repmat(nan,size(dts)); [rawix,dtsix]=intersect_dates(rawu.date(ix),dts); u(dtsix)=rawu.data(ix(rawix));
v=repmat(nan,size(dts)); [rawix,dtsix]=intersect_dates(rawv.date(ix),dts); v(dtsix)=rawv.data(ix(rawix));

[T,O]=t_tide(u+(i.*v),'start time',dts(1),'latitude',stn.lat,'synthesis',0);

fmg; plot_ts(stn.adcp_sfc_u); plot(dts,real(O),'r:'); legend('ADCP','T\_Tide'); ylim([-2.5,2.5]); titlename('LOOE1 sfc. u');
fmg; plot_ts(stn.adcp_sfc_v); plot(dts,imag(O),'r:'); legend('ADCP','T\_Tide'); ylim([-2.5,2.5]); titlename('LOOE1 sfc. u');

fmg; plot(dts,[u-real(O);v-imag(O)]); datetick3; plot_ts(ts_op(stn.ndbc_wind1_u_3_d_lp,10,'/'),':',ts_op(stn.ndbc_wind1_v_3_d_lp,10,'/'),':'); plot(dts,[real(O);imag(O)],'o'); plot_ts(stn.adcp_u_40hlp,'x',stn.adcp_v_40hlp,'x'); legend('u-tide','v-tide','U_3_d/10','V_3_d/10','tide u','tide v','u_4_0_h','v_4_0_h'); axis([datenum(2007,7,20),datenum(2007,10,15),-2.5,+2.5]); datetick3('x',2,'keeplimits'); titlename('LOOE1 sfc ADCP vs. wind');


disp('CROSS-CORRELATIONS:');

disp('40h avg vs. 3d wind');
[uix,Uix]=intersect_dates(stn.adcp_u_40hlp.date,stn.ndbc_wind1_u_3_d_lp.date);
[vix,Vix]=intersect_dates(stn.adcp_v_40hlp.date,stn.ndbc_wind1_v_3_d_lp.date);
[R,P]=corrcoef(stn.adcp_u_40hlp.data(uix),stn.ndbc_wind1_u_3_d_lp.data(Uix)); disp({'u',R(1,2)^2, P(1,2),});
[R,P]=corrcoef(stn.adcp_v_40hlp.data(vix),stn.ndbc_wind1_v_3_d_lp.data(Vix)); disp({'v',R(1,2)^2, P(1,2),});

disp('40h avg vs. 40h wind');
[uix,Uix]=intersect_dates(stn.adcp_u_40hlp.date,stn.ndbc_wind1_u_40_h_lp.date);
[vix,Vix]=intersect_dates(stn.adcp_v_40hlp.date,stn.ndbc_wind1_v_40_h_lp.date);
[R,P]=corrcoef(stn.adcp_u_40hlp.data(uix),stn.ndbc_wind1_u_40_h_lp.data(Uix)); disp({'u',R(1,2)^2, P(1,2),});
[R,P]=corrcoef(stn.adcp_v_40hlp.data(vix),stn.ndbc_wind1_v_40_h_lp.data(Vix)); disp({'v',R(1,2)^2, P(1,2),});


disp('40h sfc vs. 3d wind');
[uix,Uix]=intersect_dates(stn.adcp_sfc_u_40_h_lp.date,stn.ndbc_wind1_u_3_d_lp.date);
[vix,Vix]=intersect_dates(stn.adcp_sfc_v_40_h_lp.date,stn.ndbc_wind1_v_3_d_lp.date);
[R,P]=corrcoef(stn.adcp_sfc_u_40_h_lp.data(uix),stn.ndbc_wind1_u_3_d_lp.data(Uix)); disp({'u',R(1,2)^2, P(1,2),});
[R,P]=corrcoef(stn.adcp_sfc_v_40_h_lp.data(vix),stn.ndbc_wind1_v_3_d_lp.data(Vix)); disp({'v',R(1,2)^2, P(1,2),});

disp('40h sfc vs. 40h wind');
[uix,Uix]=intersect_dates(stn.adcp_sfc_u_40_h_lp.date,stn.ndbc_wind1_u_40_h_lp.date);
[vix,Vix]=intersect_dates(stn.adcp_sfc_v_40_h_lp.date,stn.ndbc_wind1_v_40_h_lp.date);
[R,P]=corrcoef(stn.adcp_sfc_u_40_h_lp.data(uix),stn.ndbc_wind1_u_40_h_lp.data(Uix)); disp({'u',R(1,2)^2, P(1,2),});
[R,P]=corrcoef(stn.adcp_sfc_v_40_h_lp.data(vix),stn.ndbc_wind1_v_40_h_lp.data(Vix)); disp({'v',R(1,2)^2, P(1,2),});
fh=fmg; spt(1,2,1); scatter_fit_ts(stn.adcp_sfc_u_40_h_lp,stn.ndbc_wind1_u_40_h_lp,[],[],'Sfc u_40_h_l_p','U_40_h_l_p',fh); spt(1,2,2); scatter_fit_ts(stn.adcp_sfc_v_40_h_lp,stn.ndbc_wind1_v_40_h_lp,[],[],'Sfc v_40_h_l_p','V_40_h_l_p',fh);


disp('Raw sfc vs. raw wind');
Ufinix=find(isfinite(stn.adcp_sfc_u.data));
Vfinix=find(isfinite(stn.adcp_sfc_v.data));
[uix,Uix]=intersect_dates(stn.adcp_sfc_u.date(Ufinix),stn.ndbc_wind1_u.date);
[vix,Vix]=intersect_dates(stn.adcp_sfc_v.date(Vfinix),stn.ndbc_wind1_v.date);
[R,P]=corrcoef(stn.adcp_sfc_u.data(Ufinix(uix)),stn.ndbc_wind1_u.data(Uix)); disp({'u',R(1,2)^2, P(1,2),});
[R,P]=corrcoef(stn.adcp_sfc_v.data(Vfinix(vix)),stn.ndbc_wind1_v.data(Vix)); disp({'v',R(1,2)^2, P(1,2),});
fh=fmg; spt(1,2,1); scatter_fit_ts(stn.adcp_sfc_u,stn.ndbc_wind1_u,[],[],'Sfc u','U',fh); spt(1,2,2); scatter_fit_ts(stn.adcp_sfc_v,stn.ndbc_wind1_v,[],[],'Sfc v','V',fh);

R=corrcoef(stn.adcp_sfc_u.data(Ufinix(uix))+(i.*stn.adcp_sfc_v.data(Vfinix(vix))),stn.ndbc_wind1_u.data(Uix)+(i.*stn.ndbc_wind1_v.data(Vix))); disp({'complex',R(1,2)^2,});


more on;
