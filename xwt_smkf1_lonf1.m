1;

smkf1 = load_station_data('smkf1');
lonf1 = load_station_data('lonf1');

[ix1,ix2] = intersect_dates(smkf1.air_t.date,lonf1.air_t.date);
xwt(smkf1.air_t.data(ix1),lonf1.air_t.data(ix2),'BlackandWhite');
titlename('SMKF1 vs. LONF1 T_a_i_r: 2002-2009');
print('-dpng','../figs/smkf1_air_t_xwt_lonf1_air_t.png');

[ix1,ix2] = intersect_dates(smkf1.air_t_dewp.date,lonf1.air_t_dewp.date);
xwt(smkf1.air_t_dewp.data(ix1),lonf1.air_t_dewp.data(ix2),'BlackandWhite');
titlename('SMKF1 vs. LONF1 T_d_e_w: 2004-2009');
print('-dpng','../figs/smkf1_air_t_dewp_xwt_lonf1_air_t_dewp.png');

[ix1,ix2] = intersect_dates(smkf1.barom.date,lonf1.barom.date);
xwt(smkf1.barom.data(ix1),lonf1.barom.data(ix2),'BlackandWhite');
titlename('SMKF1 vs. LONF1 P_a: 2002-2009');
print('-dpng','../figs/smkf1_barom_xwt_lonf1_barom.png');

smkf1 = verify_variable(smkf1,'wind1_u');
smkf1 = verify_variable(smkf1,'wind1_v');
lonf1 = verify_variable(lonf1,'wind1_u');
lonf1 = verify_variable(lonf1,'wind1_v');

[ix1,ix2] = intersect_dates(smkf1.wind1_u.date,lonf1.wind1_u.date);
xwt(smkf1.wind1_u.data(ix1),lonf1.wind1_u.data(ix2),'BlackandWhite');
titlename('SMKF1 vs. LONF1 W_U: 2002-2009');
print('-dpng','../figs/smkf1_wind1_u_xwt_lonf1_wind1_u.png');

[ix1,ix2] = intersect_dates(smkf1.wind1_v.date,lonf1.wind1_v.date);
xwt(smkf1.wind1_v.data(ix1),lonf1.wind1_v.data(ix2),'BlackandWhite');
titlename('SMKF1 vs. LONF1 W_V: 2002-2009');
print('-dpng','../figs/smkf1_wind1_v_xwt_lonf1_wind1_v.png');
