function trycor(stn)

dt = datenum(2003,07,15,16,00,00);

wix = find(abs(stn.ndbc_wind1_speed.date - dt) < (1/24), 1);
w = stn.ndbc_wind1_speed.data(wix);
aix = find(abs(stn.ndbc_air_t.date - dt) < (1/24), 1);
a = stn.ndbc_air_t.data(aix);
qix = find(abs(stn.ndbc_relhumid.date - dt) < (1/24), 1);
q = stn.ndbc_relhumid.data(qix);
pix = find(abs(stn.ndbc_barom.date - dt) < (1/24), 1);
p = stn.ndbc_barom.data(pix);
tix = find(abs(stn.ndbc_sea_t.date - dt) < (1/24), 1);
t = stn.ndbc_sea_t.data(tix);

dsrfix = find(abs(stn.ncep_dsrf.date - dt) < (1/24), 1);
dsrf = stn.ncep_dsrf.data(dsrfix);
dlrfix = find(abs(stn.ncep_dlrf.date - dt) < (1/24), 1);
dlrf = stn.ncep_dlrf.data(dlrfix);

prcp = 0;


[wz,az,pz,stz] = station_instrument_heights(stn.station_name);

p = p ./ exp( -pz ./ ( (a + 273.15) .* 29.263 ) );


MHOME = 'c:/Documents and Settings/gramer/My Documents/MATLAB';
rmpath([MHOME '/fairall']);
result = hfbulktc(w,wz,a,az,q,az,p,t);
rhf20 = result(:,1);
lhf20 = result(:,2) + result(:,3);
nhf20 = rhf20 + lhf20;
disp([20 rhf20 lhf20 nhf20]);


% Air specific humidity [kg/kg]
sa = relhumid_to_spechumid(a,q);
sa = sa .* 1e3;
% Saturated ("sea-surface") specific humidity [kg/kg]
ss = relhumid_to_spechumid(t,100);
ss = ss .* 1e3;

% (Atmospheric) Planetary Boundary Layer height == inversion height [m]
pblz = 600;

% Ocean current estimate (Shay et al. 1998)
ou = 0.015 .* w;

addpath([MHOME '/fairall']);
res = cor30a([w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,wz,az,az,stn.lat,1,0,0,0]);
result = res;
rhf26 = -result(:,1);
lhf26 = -result(:,2);
nhf26 = rhf26 + lhf26;
disp([26 rhf26 lhf26 nhf26]);

res = cor30a([w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,wz,az,az,stn.lat,1,1,2,2]);
result = res;
rhf30 = -result(:,1);
lhf30 = -result(:,2);
nhf30 = rhf30 + lhf30;
disp([30 rhf30 lhf30 nhf30]);
rmpath([MHOME '/fairall']);

return;
