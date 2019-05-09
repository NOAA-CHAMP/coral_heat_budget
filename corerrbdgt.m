function [dts,hs,dwth,hl,dwq] = corerrbdgt(stn)
%function [dts,hs,dwth,hl,dwq] = corerrbdgt(stn)
%
% Calculate latent and sensible heat flux errors using COARE 3.0a method
%
% Last Saved Time-stamp: <Thu 2011-01-06 07:47:11  lew.gramer>

if ( ~isfield(stn,'ndbc_sea_spechumid') )
  stn.ndbc_sea_spechumid.date = stn.ndbc_sea_t.date;
  % 0.98 factor from Stommel - accounts for salinity
  stn.ndbc_sea_spechumid.data = 0.98 .* relhumid_to_spechumid(stn.ndbc_sea_t.data,100);
end;

%Input sensor heights
% zu=18;
% zt=18;
% zq=18;
[zu,zt] = station_instrument_heights(stn.station_name);
zq=zt;

%error specs  [intercept error   slope error]
%Xtrue=+-interror+(1+-sloperror)*Xmeas
ue=[.55 .03];%'true' wind speed relative to sea surface
% .03 slope error verified by duplicate wind data at MLRF1,FWYF1,SMKF1
te=[.09 0];%sea-air temperature difference
qe=[0.18 0];%sea-air spechumid difference

  [six,aix,qaix,qsix,Wix,hix,qlhix,qshix,swix,lwix] = ...
      intersect_all_dates([], stn.ndbc_sea_t.date,stn.ndbc_air_t.date,...
                          stn.ndbc_spechumid.date,stn.ndbc_sea_spechumid.date,...
                          stn.ndbc_wind1_speed.date,...
                          stn.tmd_tide_i_depth.date,...
                          stn.ndbc_ncep_30a_latent_heat_flux.date,...
                          stn.ndbc_ncep_30a_sensible_heat_flux.date,...
                          stn.ncep_dsrf.date,stn.ncep_dlrf.date);

  dts = stn.ndbc_sea_t.date(six);

  s = stn.ndbc_sea_t.data(six);
  a = stn.ndbc_air_t.data(aix);
  qss = stn.ndbc_sea_spechumid.data(qsix);
  qsa = stn.ndbc_spechumid.data(qaix);
  W = kts2mps(stn.ndbc_wind1_speed.data(Wix));
  h = stn.tmd_tide_i_depth.data(hix);
  qlh = real(stn.ndbc_ncep_30a_latent_heat_flux.data(qlhix));
  qsh = real(stn.ndbc_ncep_30a_sensible_heat_flux.data(qshix));
  dsrf = stn.ncep_dsrf.data(swix);
  dlrf = stn.ncep_dlrf.data(lwix);


if ( isempty(strfind(path, 'fairall')) )
  FAIRALLHOME = get_ecoforecasts_path('../fairall');
  addpath(FAIRALLHOME);
  clear FAIRALLHOME
  rehash
end;

disp('COR30A (no warm-layer)');
disp(length(W));
nby10 = floor(length(W)/10);

for ix=1:length(W)

ta=a(ix);
Ts=s(ix);
qa=qsa(ix).*1e3;  % Algorithm expects [g/kg]
qs=qss(ix).*1e3;  % Algorithm expects [g/kg]
u=W(ix);
qsw=dsrf(ix);
qlw=dlrf(ix);

bet=1.25;
k=0.4;
al=9.8/(ta+273.15);

ecu=.001/10;
ect=.001/10;
ecq=.001/10;

[wvper,wvhgt] = wind_to_wave(u);

% x=[u 0 Ts ta qs qa qsw qlw 0 600 1016 zu zt zq 25 1 0 5 1 ];%coare30, cool on, wave off
x=[u 0 Ts ta qs qa qsw qlw 0 600 1016 zu zt zq 25 1 1 wvper wvhgt];%coare30, cool on, wave on
dq=qs-qa;
dt=Ts-ta;

dth=dt+6.1e-4*(ta+273.16)*dq;

    zi=x(10);
    x(1)=u;
    x(3)=Ts;
    x(4)=ta;
    x(5)=qs;
    x(6)=qs-dq;

    eu=sqrt(ue(1)^2+(ue(2)*u)^2);
    et=sqrt(te(1)^2+(te(2)*dt)^2);
    eq=sqrt(qe(1)^2+(qe(2)*dq)^2);

    y=cor30a(x);
    if (mod(ix,nby10)==1); disp(ix); end;

   %y=[hsb hlb tau zo zot zoq L usr tsr qsr dter dqer tkt RF wbar Cd Ch Ce Cdn_10 Chn_10 Cen_10 ug];
   %    1   2   3   4  5   6  7  8   9  10   11   12  13  14  15  16 17 18    19      20    21  22

    hs(ix)=-y(1);
    hl(ix)=-y(2);
    tau=y(3);
    hb=hs(ix)+6.1e-4*300*hl(ix)/2.5;
    L=y(7);
    zetu=zu/L;
    zett=zt/L;
    zetq=zq/L;

    ztn=10/L;
    cunh=k/log(zu/y(4));
    ctnh=k/log(zt/y(5));
    cqnh=k/log(zq/y(6));

    dcu=(psiu_30(zetu)-psiu_30(.99*zetu))/.01;
    dct=(psit_30(zett)-psit_30(.99*zett))/.01;
    dcq=(psit_30(zetq)-psit_30(.99*zetq))/.01;;

    usr=y(8);
    tsr=y(9);
    qsr=y(10);
    ws=(al*usr*abs(tsr+6.1e-4*300*qsr)*zi)^.333;

    su=(1-cunh/k*psiu_30(zetu));
    st=(1-ctnh/k*psit_30(zett));
    sq=(1-cqnh/k*psit_30(zetq));
    uw=sqrt(u*u+(bet*ws)^2);
    aw=((bet*ws)/uw)^2;

    fs=cunh/k*dcu/su+aw/3;
    fu=cunh/k*dcu/su;
    ft=ctnh/k*dct/st;
    fq=cqnh/k*dcq/sq;

    %***** xfer coeff error contribution
    xs1=ecu/su/ctnh;
    xt1=ect/st/ctnh;
    xq1=ecq/sq/cqnh;
    xu1=ecu/su/cunh;
    %***** sensor error contribution

    xs2=eu/u*(u/uw)^2;
    xt2=et/dt;
    xq2=eq/dq;
    xu2=eu/u;
    
    %************** analytical error functions
    d=(1-ft)*(1-aw)+2*fs;%universal denominator
    dus=sqrt((1-ft)^2*(xs1^2+xs2^2)+fs^2*(xt1^2+xt2^2))/d;
    dzet=sqrt((1-aw)^2*(xt1^2+xt2^2)+4*(xs1^2+xs2^2))/d;

    dwth(ix)=sqrt(((1-aw)+3*fs)^2*(xt1^2+xt2^2)+(xs1^2+xs2^2)*(1-3*ft)^2)/d;
    dwq(ix)=sqrt(xq1^2+xq2^2+((fq*(1-aw)+fs)^2*(xt1^2+xt2^2)+((1-ft)-2*fq)^2*(xs1^2+xs2^2))/d.^2);
    dwu(ix)= sqrt((xu1*((1-ft)*(2-aw)+.67*aw))^2+(xu2*((1-ft)*(1+(u/uw)^2-aw)+2*fu*(1-(u/uw)^2)+.67*aw))^2+(xt1^2+xt2^2)*(fu*(2-aw)+.33*aw)^2)/d;

end;

% Above loop calculates relative (normalized) errors - change to absolute
dwth = dwth.*hs;
dwq = dwq.*hl;
