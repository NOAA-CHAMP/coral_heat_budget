1;

load('../data/fwyf1_ms.mat');

stn = station; station = []; clear station;

stn = get_gom_hycom(stn);

[ud,dx,dy] = calc_udotdelt(stn.gom_hycom_u.data,stn.gom_hycom_v.data,stn.gom_hycom_seatemp_field);

stn.ud.date = stn.gom_hycom_u.date;
stn.ud.data = ud;

stn.uds.date = stn.gom_hycom_u.date(1):(1/24):stn.gom_hycom_u.date(end);
stn.uds.data = spline(stn.ud.date,stn.ud.data,stn.uds.date);

[ix1,ix2] = intersect_dates(stn.netqf.date,stn.uds.date);
stn.dt.date = stn.netqf.date(ix1);
stn.dt.data = stn.netqf.data(ix1) + stn.uds.data(ix2)';

plot_fluxes(stn,2003);
